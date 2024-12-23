from abc import ABC, abstractmethod
import collections
import datetime
import itertools
import logging
import time
from typing import Iterable, List, Optional, Tuple

import numpy as np

from fastai.basic_data import DeviceDataLoader
import torch
from torch import Tensor
import torch.nn as nn
from torch.optim.lbfgs import LBFGS
from torch.optim.adamw import AdamW
from torch.utils.data import TensorDataset

from apu import TrainingLogger, config
from apu.datasets import synthetic
# noinspection PyUnresolvedReferences
from apu.datasets.types import BaseFFModule, Labels, TensorGroup
# noinspection PyUnresolvedReferences
from apu.loss import PN, log_loss, loss_0_1, ramp_loss, sigmoid_loss
from apu.nnpu import NNPU
from apu.purr import PURR
from apu.two_step import APNU, Step1Method, WUU
from apu.types import LearnerParams
import apu.utils
from apu.utils import ClassifierBlock, NUM_WORKERS, TORCH_DEVICE, shuffle_tensors

from generate_results import calculate_results
from puc_learner import PUcLearner

SplitTensorInfo = collections.namedtuple("SplitTensorInfo", ["x", "y", "sigma", "tk"])

ALL_STEP1_METHODS = False


class _BaseLearner(nn.Module, ABC):
    def __init__(self, name: str):
        super(_BaseLearner, self).__init__()
        self._logger = self._train_start = None
        self._name = name

    def _configure_fit_vars(self, modules: nn.ModuleDict):
        r""" Set initial values/construct all variables used in a fit method """
        # Fields that apply regardless of loss method
        tb_dir = apu.utils.BASE_DIR / "tb"
        TrainingLogger.create_tensorboard(tb_dir)
        names, sizes = [], []
        for _, module in modules.items():
            _name, _size = module.logger_field_info()
            names.extend(_name)
            sizes.extend(_size)
        # Always log the time in number of seconds
        names.append("Time")
        sizes.append(10)
        self._logger = TrainingLogger(names, sizes, logger_name=apu.utils.LOGGER_NAME,
                                      tb_grp_name=self._name)

    def _log_epoch(self, ep: int, modules: nn.ModuleDict) -> None:
        r"""
        Log the results of the epoch
        :param ep: Epoch number
        :param modules: Modules to log
        """
        flds = []
        for _, module in modules.items():
            flds.extend(module.epoch_log_fields())

        flds.append(time.time() - self._train_start)
        self._logger.log(ep, flds)

    def train_start_time(self) -> str:
        r""" Returns the training start time as a string """
        assert self._train_start is not None, "Training never started"
        return datetime.datetime.fromtimestamp(self._train_start).strftime("%Y-%m-%d-%H-%M-%S")

    @abstractmethod
    def fit(self, tg: TensorGroup):
        r""" Fit all models """

    def _fit(self, modules: nn.ModuleDict, train_dl: DeviceDataLoader, valid_dl: DeviceDataLoader):
        r""" Fits \p modules' learners to the training and validation \p DataLoader objects """
        self._configure_fit_vars(modules)

        for mod_name, module in modules.items():
            lr = config.get_learner_val(mod_name, LearnerParams.Attribute.LEARNING_RATE)
            wd = config.get_learner_val(mod_name, LearnerParams.Attribute.WEIGHT_DECAY)
            is_lin_ff = config.DATASET.is_synthetic() and module.module.num_hidden_layers == 0
            if is_lin_ff:
                module.optim = LBFGS(module.parameters(), lr=lr)
            else:
                module.optim = AdamW(module.parameters(), lr=lr, weight_decay=wd, amsgrad=True)
            logging.debug(f"{mod_name} Optimizer: {module.optim.__class__.__name__}")

        for ep in range(1, config.NUM_EPOCH + 1):
            # noinspection PyUnresolvedReferences
            for _, module in modules.items():
                module.epoch_start()

            for batch in train_dl:
                for _, module in modules.items():
                    module.process_batch(batch)

            for _, module in modules.items():
                module.calc_valid_loss(valid_dl)
            self._log_epoch(ep, modules)
        self._restore_best_model(modules)
        self.eval()

    @staticmethod
    def _restore_best_model(modules: nn.ModuleDict):
        r""" Restores the best trained model from disk """
        for _, module in modules.items():
            module.restore_best()


class CalibratedLearner(_BaseLearner):
    r""" Simple module used to represent a calibrated module """
    def __init__(self, prior: float, gamma: float, base_module: nn.Module):
        super().__init__(name="Calibrated_nnPU")

        train_loss = log_loss
        val_loss = log_loss

        nnpu = NNPU(prior=prior, gamma=gamma, train_loss=train_loss, valid_loss=val_loss,
                    only_u_train=True, abs_nn=True)

        self._mod = ClassifierBlock(net=base_module, estimator=nnpu)
        self.block = self._mod.module

        self._blocks = nn.ModuleDict()
        self._blocks[self._name] = self._mod

        self._sigmoid = nn.Sigmoid()
        self._is_cal = torch.full((), False, dtype=torch.bool, device=TORCH_DEVICE)

        self._prior = prior
        self.cal_err = None

        self.to(device=TORCH_DEVICE)

    def forward(self, x: Tensor) -> Tensor:
        r""" Passes through the calibrated module and if relevant adds the sigmoid function """
        y = self._mod.forward(x)
        if self.is_calibrated():
            y = self._sigmoid.forward(y)
        return y

    def set_calibrated(self) -> None:
        r""" Mark the module as calibrated"""
        self._is_cal = torch.full(self._is_cal.shape, True, dtype=torch.bool, device=TORCH_DEVICE)

    def is_calibrated(self) -> bool:
        r""" Returns \p True if the module is calibrated """
        return bool(self._is_cal.item())

    def freeze(self) -> None:
        r""" Disable training on the block """
        self.eval()
        for module in self.modules():
            module.eval()
            for param in module.parameters():
                param.requires_grad = False
        # self.set_calibrated()

    def fit(self, tg: TensorGroup) -> None:
        r""" Fit all models """
        self._train_start = time.time()

        size_p, p_in_u = tg.p_x.shape[0], self._prior * tg.u_tr_x.shape[0]
        self.cal_err = (size_p + p_in_u) / size_p

        train_dl, valid_dl = create_pu_dataloaders(p_x=tg.p_x, u_x=tg.u_tr_x,
                                                   bs=config.SIGMA_BATCH_SIZE)
        self._fit(self._blocks, train_dl=train_dl, valid_dl=valid_dl)

        # Mark calibrated and frozen
        self.set_calibrated()
        self.freeze()

    def calc_cal_weights(self, x: Tensor) -> Tensor:
        r"""
        Calculate the calibrated weights.
        :return: Tuple of the sigma probability
        """
        assert self.is_calibrated(), "Sigma learner not calibrated"
        return self.forward(x)


class SklearnCalibrated(_BaseLearner):
    r""" Simple module used to represent a calibrated module """
    def __init__(self, prior: float, base, base_name: str):
        super().__init__(name=f"Calibrated_Sk-{base_name}")

        self._is_cal = False
        self._base = base

        self._prior = prior
        self.cal_err = None

    def forward(self, x: Tensor) -> Tensor:
        r""" Passes through the calibrated module and if relevant adds the sigmoid function """
        x = x.cpu().numpy()
        try:
            y = self._base.predict_proba(x)
            # In case of two dimensional probability predictions take the second dimension
            if y.shape[1] == 2:
                y = y[:, 1]
        except AttributeError:
            y = self._base.predict(x)

        y = torch.from_numpy(y).squeeze().to(TORCH_DEVICE)
        # y = (y * self.cal_err).clamp(0., 1.)
        return y

    def set_calibrated(self) -> None:
        r""" Mark the module as calibrated"""
        self._is_cal = True

    def is_calibrated(self) -> bool:
        r""" Returns \p True if the module is calibrated """
        return self._is_cal

    def fit(self, tg: TensorGroup) -> None:
        r""" Fit all models """
        self._train_start = time.time()

        size_p, p_in_u = tg.p_x.shape[0], self._prior * tg.u_tr_x.shape[0]
        self.cal_err = (size_p + p_in_u) / size_p

        msg = "Fitting Sklearn probabilistic classifier"
        logging.debug(f"Starting: {msg}")
        x = torch.cat((tg.p_x.cpu(), tg.u_tr_x.cpu()), dim=0).float().numpy()
        y = torch.cat((torch.full([tg.p_x.shape[0]], fill_value=Labels.POS),
                       torch.full([tg.u_tr_x.shape[0]], fill_value=Labels.NEG)),
                      dim=0).int().squeeze().numpy()

        self._base.fit(x, y)

        # Mark calibrated and unfrozen
        self.set_calibrated()
        logging.debug(f"COMPLETED: {msg}")

    def freeze(self):
        r""" No effect of freeze.  Added for common API """
        pass

    @staticmethod
    def _restore_best_model(modules: nn.ModuleDict):
        assert False, "Restore not supported"


# noinspection PyPep8Naming
class APU_Learner(_BaseLearner):
    r"""
    Single learner that trains multiple modules using the wUU, aPNU, PURR, nnPU, and PN risk
    estimators.
    """
    def __init__(self, base_module: nn.Module, sigma: Optional[CalibratedLearner] = None,
                 rho_vals: Optional[Iterable] = np.linspace(start=0.0, stop=1.0, num=11),
                 bias_priors: bool = False):
        super().__init__("aPU")

        # train_loss = log_loss
        # train_loss = ramp_loss
        train_loss = sigmoid_loss
        valid_loss = sigmoid_loss
        # valid_loss = loss_0_1

        self._pu_blocks = nn.ModuleDict()

        tr_priors, te_priors = [config.TRAIN_PRIOR], [config.TEST_PRIOR]
        if bias_priors:
            biases = [0.8, 1.2]
            tr_priors.extend([bias * config.TRAIN_PRIOR for bias in biases])
            te_priors.extend([bias * config.TEST_PRIOR for bias in biases])

        self.s1_soft_acc = self.s1_hard_acc = self.s1_topk_acc = None
        self._sigma = sigma
        if not self._sigma.is_calibrated():
            raise ValueError("Sigma method is not calibrated")
        self._sigma.freeze()

        s1_methods = self.get_step1_methods()
        wuu_itrs = itertools.product(s1_methods, tr_priors, te_priors)
        for s1_method, tr_prior, te_prior in wuu_itrs:
            wuu = WUU(train_prior=tr_prior, test_prior=te_prior,
                      u_tr_label=Labels.Training.U_TRAIN.value,
                      u_te_label=Labels.Training.U_TEST.value,
                      train_loss=train_loss, valid_loss=valid_loss, abs_nn=True,
                      s1_method=s1_method)

            wuu.gamma = config.get_learner_val(wuu.name(), LearnerParams.Attribute.GAMMA)
            self._pu_blocks[wuu.name()] = apu.utils.ClassifierBlock(base_module, wuu)

        purr_itrs = itertools.product(tr_priors, te_priors)
        for tr_prior, te_prior in purr_itrs:
            l_purr = PURR(train_prior=tr_prior, test_prior=te_prior,
                          train_loss=train_loss, valid_loss=valid_loss, abs_nn=True)

            l_purr.gamma = config.get_learner_val(l_purr.name(), LearnerParams.Attribute.GAMMA)
            self._pu_blocks[l_purr.name()] = apu.utils.ClassifierBlock(base_module, l_purr)

        apnu_itrs = itertools.product(s1_methods, rho_vals, tr_priors, te_priors)
        for s1_method, rho, tr_prior, te_prior in apnu_itrs:
            nn_pnu = APNU(train_prior=tr_prior, test_prior=te_prior, rho=rho,
                          train_loss=train_loss, valid_loss=valid_loss,
                          abs_nn=True, s1_method=s1_method)

            nn_pnu.gamma = config.get_learner_val(nn_pnu.name(), LearnerParams.Attribute.GAMMA)
            classifier = apu.utils.ClassifierBlock(base_module, nn_pnu)
            self._pu_blocks[nn_pnu.name()] = classifier

        # Construct the nnPU learners
        for i in range(2):
            if i == 0:
                num_unlabeled_pos = (config.TRAIN_PRIOR * config.N_U_TRAIN
                                     + config.TEST_PRIOR * config.N_U_TEST)
                tot_u_size = config.N_U_TRAIN + config.N_U_TEST
                prior = num_unlabeled_pos / tot_u_size

                only_u_test = False
            elif i == 1:
                prior = config.TEST_PRIOR
                only_u_test = True
            else:
                raise ValueError("Unknown configuration")
            nnpu = NNPU(prior=prior, train_loss=train_loss, valid_loss=valid_loss,
                        only_u_test=only_u_test)
            nnpu.gamma = config.get_learner_val(nnpu.name(), LearnerParams.Attribute.GAMMA)
            self._pu_blocks[nnpu.name()] = ClassifierBlock(base_module, nnpu)

        pn = PN(train_loss=train_loss, valid_loss=valid_loss, is_train=False,
                prior=config.TEST_PRIOR)
        self._pn_test = nn.ModuleDict({"te_pn": ClassifierBlock(base_module, pn)})

        self.to(device=TORCH_DEVICE)

    def fit(self, tg: TensorGroup):
        r""" Fit all models """
        self._train_start = time.time()

        # If weights object exists, define the weights
        if self._sigma is not None:
            tg.set_sigmas(sigma_module=self._sigma, config=config, device=TORCH_DEVICE)
            self._calc_step1_acc(tg=tg)

        train, valid = create_apu_dataloaders(tg=tg, bs=config.BATCH_SIZE)
        self._fit(self._pu_blocks, train_dl=train, valid_dl=valid)

        tot_tensor_size = config.N_P + config.N_U_TRAIN + config.N_U_TEST
        # Scale batch size to match only TEST unlabeled set
        te_bs = int(config.BATCH_SIZE * tg.u_te_x.shape[0] / tot_tensor_size)
        train, valid = create_pn_dataloaders(tg.u_te_x, tg.u_te_y, bs=te_bs)
        self._fit(self._pn_test, train_dl=train, valid_dl=valid)

    def _calc_step1_acc(self, tg: TensorGroup) -> None:
        r""" Helper method that calculates the accuracy of the step 1 methods """
        apu.utils.export_step1_results(tg=tg)

        n_u_tr = tg.u_tr_y.numel()

        def _convert_to_lbl(_w: Tensor) -> Tensor:
            p_mask = _w == Labels.POS
            _w[~p_mask] = Labels.NEG
            return _w

        def _calc_acc(y_hat: Tensor) -> float:
            # Accuracy of step 1
            n_correct = int((y_hat == tg.u_tr_y).sum())
            return n_correct / n_u_tr

        err_vec = (2 * tg.u_tr_sigma - 1 - tg.u_tr_y.float()) / 2
        self.s1_soft_acc = 1. - float(err_vec.abs().mean())
        logging.debug(f"Step #1 - Soft Accuracy: {100 * self.s1_soft_acc:.1f}")

        sigma_lbl = tg.u_tr_sigma.round().clamp(0, 1).type(tg.u_tr_y.dtype)
        sigma_lbl = _convert_to_lbl(sigma_lbl)
        self.s1_hard_acc = _calc_acc(sigma_lbl)
        logging.debug(f"Step #1 - Hard Accuracy: {100 * self.s1_hard_acc:.1f}")

        tk_lbl = _convert_to_lbl(tg.u_tr_tk.type(tg.u_tr_y.dtype))
        self.s1_topk_acc = _calc_acc(tk_lbl)
        logging.debug(f"Step #1 - Top-K Accuracy: {100 * self.s1_topk_acc:.1f}")

    def forward(self, x: Tensor) -> dict:
        # noinspection PyUnresolvedReferences
        return {key: block.forward(x) for key, block in self.blocks()}

    def blocks(self) -> Iterable:
        r""" Iterates through all the blocks """
        modules, n = nn.ModuleDict(), 0
        for _mod_dict in (self._pu_blocks, self._pn_test):
            # noinspection PyTypeChecker
            modules.update(_mod_dict.items())
            assert len(modules) == n + len(_mod_dict), "Duplicate key detected"
            n += len(_mod_dict)
        return modules.items()

    @staticmethod
    def get_step1_methods() -> List[Step1Method]:
        r""" Method to determine which step1 methods to try """
        if ALL_STEP1_METHODS:
            return [method for method in Step1Method]

        if not config.DATASET.is_spam():
            return [Step1Method.SOFT]
        else:
            return [Step1Method.TOP_K]


def _split_tensor(x: Tensor, lbl: int, sigma: Tensor, tk: Tensor) -> \
        Tuple[SplitTensorInfo, SplitTensorInfo]:
    r""" Splits a tensor into train and validation """
    frac_train = 1 - config.VALIDATION_SPLIT_RATIO

    tr_size, n_ele = int(frac_train * x.shape[0]), x.shape[0]
    # Size of the validation set
    val_size = n_ele - tr_size

    # Randomly select the x and sigma elements
    idx = torch.randperm(n_ele)
    tr_idx, val_idx = idx[:tr_size], idx[tr_size:]

    x_tr, x_val = x[tr_idx], x[val_idx]
    sigma_tr, sigma_val = sigma[tr_idx], sigma[val_idx]
    tk_tr, tk_val = tk[tr_idx], tk[val_idx]

    def _build_y(size: int) -> Tensor:
        r""" Standardizes constructing \p y Tensor"""
        return torch.full([size], lbl, dtype=torch.int32, device="cpu")

    # y tensors are constant throughout
    spl_tr = SplitTensorInfo(x_tr, _build_y(tr_size), sigma_tr, tk_tr)
    spl_val = SplitTensorInfo(x_val, _build_y(val_size), sigma_val, tk_val)
    return spl_tr, spl_val


def create_pu_dataloaders(p_x: Tensor, u_x: Tensor, bs: int) \
        -> Tuple[DeviceDataLoader, DeviceDataLoader]:
    r"""
    Simple method that splits the positive and unlabeled sets into stratified training and
    validation \p DataLoader objects

    :param p_x: Feature vectors for the positive (labeled) examples
    :param u_x: Feature vectors for the unlabeled examples
    :param bs: \p DataLoader's batch size
    :return: Training and validation \p DataLoader objects respectively
    """
    tr_x, tr_y, val_x, val_y = [], [], [], []
    for x, lbl in ((p_x, Labels.Training.POS), (u_x, Labels.Training.U_TRAIN)):
        num_ex = x.shape[0]
        tr_size = int((1 - config.VALIDATION_SPLIT_RATIO) * num_ex)
        x = shuffle_tensors(x)

        tr_x.append(x[:tr_size])
        tr_y.append(torch.full([tr_size], lbl.value, dtype=torch.int))

        val_x.append(x[tr_size:])
        val_y.append(torch.full([num_ex - tr_size], lbl.value, dtype=torch.int))

    def _cat_tensors(lst_tensors: List[Tensor]) -> Tensor:
        return torch.cat(lst_tensors, dim=0).cpu()

    tr_x, tr_y = shuffle_tensors(_cat_tensors(tr_x), _cat_tensors(tr_y))
    val_x, val_y = shuffle_tensors(_cat_tensors(val_x), _cat_tensors(val_y))

    # Create the training and validation dataloaders respectively
    dls = []
    for x, y, shuffle in ((tr_x, tr_y, True), (val_x, val_y, False)):
        tk = w = -torch.ones_like(y).float().cpu()
        dls.append(DeviceDataLoader.create(TensorDataset(x.cpu(), y.cpu(), w.cpu(), tk.cpu()),
                                           shuffle=shuffle,
                                           drop_last=shuffle and not config.DATASET.is_synthetic(),
                                           bs=bs, num_workers=NUM_WORKERS, device=TORCH_DEVICE))
    # noinspection PyTypeChecker
    return tuple(dls)


def create_apu_dataloaders(tg: TensorGroup, bs: int) -> Tuple[DeviceDataLoader, DeviceDataLoader]:
    r"""
    Creates the training and validation dataloaders
    :param tg: Stores the raw tensor information
    :param bs: \p DataLoader's batch size
    :return: Training, validation and optionally calibration \p DataLoader objects respectively
    """
    all_tr, all_val = [], []
    # Split the tensors into train/validation
    flds = (("p", Labels.Training.POS), ("u_tr", Labels.Training.U_TRAIN),
            ("u_te", Labels.Training.U_TEST))
    for ds_name, lbl in flds:
        x = tg.__getattribute__(f"{ds_name}_x")
        sigma = tg.__getattribute__(f"{ds_name}_sigma")
        if sigma is None:
            sigma = -torch.ones(x.shape[:1], dtype=torch.float, device=TORCH_DEVICE)
        tk = tg.__getattribute__(f"{ds_name}_tk")
        if tk is None:
            tk = -torch.ones(x.shape[:1], dtype=torch.float, device=TORCH_DEVICE)

        spl_tr, spl_val = _split_tensor(x, lbl.value, sigma, tk)
        all_tr.append(spl_tr)
        all_val.append(spl_val)

    # Construct the individual dataloaders
    dls = []
    flds = ((all_tr, True), (all_val, False))
    for spl_info, shuffle in flds:
        if not spl_info:
            dls.append(None)
            continue

        x = torch.cat([info.x for info in spl_info], dim=0).cpu()
        y = torch.cat([info.y for info in spl_info], dim=0).cpu()
        sigma = torch.cat([info.sigma for info in spl_info], dim=0).cpu()
        tk = torch.cat([info.tk for info in spl_info], dim=0).cpu()

        dl = DeviceDataLoader.create(dataset=TensorDataset(x, y, sigma, tk), shuffle=shuffle,
                                     drop_last=shuffle, bs=bs, num_workers=NUM_WORKERS,
                                     device=TORCH_DEVICE)
        dls.append(dl)
    # train, validation, and calibration dataloaders respectively
    # noinspection PyTypeChecker
    return tuple(dls)


# noinspection DuplicatedCode
def create_pn_dataloaders(x: Tensor, y: Tensor, bs: int) \
        -> Tuple[DeviceDataLoader, DeviceDataLoader]:
    r"""
    Creates training and validation PN \p DataLoader objects from the specified tensors
    :param x: Feature values tensor
    :param y: Labels tensor
    :param bs: \p DataLoader's batch size
    :return: Tuple of the training and validation \p DataLoader objects respectively
    """
    num_ele = x.shape[0]
    assert num_ele == y.shape[0], "Mismatch in number of elements"

    tr_size = int(round(num_ele * (1. - config.VALIDATION_SPLIT_RATIO)))
    x, y = shuffle_tensors(x, y)
    tk = sigma = -torch.ones(y.shape).float().cpu()

    tensor_ds = TensorDataset(x[:tr_size], y[:tr_size], sigma[:tr_size], tk[:tr_size])
    train_dl = DeviceDataLoader.create(dataset=tensor_ds, shuffle=True, drop_last=True, bs=bs,
                                       num_workers=NUM_WORKERS, device=TORCH_DEVICE)

    tensor_ds = TensorDataset(x[tr_size:], y[tr_size:], sigma[tr_size:], tk[tr_size:])
    valid_dl = DeviceDataLoader.create(dataset=tensor_ds,
                                       shuffle=False, drop_last=False, bs=bs,
                                       num_workers=NUM_WORKERS, device=TORCH_DEVICE)

    return train_dl, valid_dl


def execute_test(tg: TensorGroup, module: nn.Module, bias_priors: bool = False,
                 use_sklearn_cal: bool = False):
    if use_sklearn_cal:
        raise ValueError("Install LightGBM to use scikit-learn mode, then comment out this line")
        # import lightgbm as lgb
        # base = lgb.LGBMClassifier()
        # cal_learner = SklearnCalibrated(prior=config.TRAIN_PRIOR, base=base, base_name="LGBM")
    else:
        if config.DATASET.is_synthetic():
            cal_module = synthetic.Module()
        else:
            cal_module = BaseFFModule(x=tg.p_x, num_hidden_layers=config.NUM_SIGMA_LAYERS)
        cal_learner = CalibratedLearner(gamma=config.GAMMA, prior=config.TRAIN_PRIOR,
                                        base_module=cal_module)
    cal_learner.fit(tg)
    if config.DATASET.is_synthetic():
        apu.utils.log_decision_boundary(cal_learner.block, name="Sigma")

    rho_vals = [0.5]
    apu_learners = APU_Learner(base_module=module, sigma=cal_learner, rho_vals=rho_vals,
                               bias_priors=bias_priors)
    apu_learners.fit(tg)

    priors = [config.TRAIN_PRIOR]
    if bias_priors:
        scalars = [0.8, 1.2]
        priors.extend(scale * config.TRAIN_PRIOR for scale in scalars)

    pucs = []
    for prior in priors:
        puc_learner = PUcLearner(prior=prior)
        puc_learner.fit(tg)
        pucs.append(puc_learner)

    calculate_results(tg, apu_learners, pucs)
