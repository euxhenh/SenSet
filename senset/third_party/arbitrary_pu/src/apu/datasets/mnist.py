# -*- utf-8 -*-
r"""
    mnist.py
    ~~~~~~~~~~~~~~~

    Downloads and returns MNIST and MNIST-variant datasets (e.g.,~FashionMNIST and KMNIST)

    :copyright: (c) 2020 by Zayd Hammoudeh.
    :license: MIT, see LICENSE file for more details.
"""

__all__ = ["load_data"]

from pathlib import Path

import torchvision

from .. import _config as config
from .utils import shared_tensor_dataset_importer
from .types import APU_Dataset, TensorGroup

MNIST_NORMALIZE_FACTOR = 255


def load_data(dest: Path) -> TensorGroup:
    r"""
    Loads the MNIST, FashionMNIST, and KMNIST datasets
    :return: \p TensorGroup of extracted data
    """
    if config.DATASET == APU_Dataset.MNIST:
        torchvision.datasets.MNIST(dest, download=True)
        name = "MNIST"
    elif config.DATASET == APU_Dataset.FASHION_MNIST:
        torchvision.datasets.FashionMNIST(dest, download=True)
        name = "FashionMNIST"
    elif config.DATASET == APU_Dataset.KMNIST:
        torchvision.datasets.KMNIST(dest, download=True)
        name = "KMNIST"
    else:
        raise ValueError("Unknown MNIST(-variant) dataset")

    dest /= f"{name}/processed"
    dest.mkdir(parents=True, exist_ok=True)
    return shared_tensor_dataset_importer(dest=dest, normalize_factor=MNIST_NORMALIZE_FACTOR)
