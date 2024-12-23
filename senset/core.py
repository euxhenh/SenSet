from typing import Dict, List, Tuple

import anndata
import numpy as np
import scipy.sparse as sp
import torch
from scipy.stats import ranksums
from sklearn.decomposition import PCA

from .third_party.arbitrary_pu.src.my_types import TensorGroup
from .third_party.arbitrary_pu.src.puc_learner import PUcLearner
from .utils import _compute_log2fc, pqvals


def PUC_ct(
    adata: anndata.AnnData,
    cell_type: str,
    age_thresh: Tuple[float, float] = (30, 50),
    n_components: int = 10,
    prior: float = 0.9,
    min_n_cells: int = 50,
    known_markers: List[str] | None = None,
    age_key: str = 'age_or_mean_of_age_range',
):
    """Runs PUc results for cell_type.

    Parameters
    ----------
    adata: Annotated Data object
    cell_type: string
        The cell type to consider in adata.
    age_thresh: tuple[int, int]
        Age thresholds used to split the data into 3 age groups. These are
        upper limits (exclusive).
    n_components: int
        Number of components for PCA embeddings.
    prior: float
        The mixing coefficient. Corresponds to the fraction of
        healthy cells expected in the older age groups.
    min_n_cells: int
        Cell types with fewer than `min_n_cells` cells for any age group
        are filtered out.
    known_markers: array-like
        If not None, will use only these genes.
    """
    argidx = np.argwhere(adata.obs['cell_type'] == cell_type).flatten()
    if known_markers is not None:
        ct_adata = adata[argidx, np.in1d(adata.var['feature_name'], known_markers)]
    else:
        ct_adata = adata[argidx]
    genes = ct_adata.var['feature_name'].to_numpy()
    age = ct_adata.obs[age_key]

    # y consists of 0, 1, 2 depending on age group
    y = np.where(age < age_thresh[0], 0, 1)  # 0 for <age_thresh[0] else 1
    y[age >= age_thresh[1]] = 2  # 2 for >=age_thresh[1]

    _, counts = np.unique(y, return_counts=True)
    # skip if few cells in any age group
    if np.any(counts < min_n_cells):
        print(f"Skipping low count {cell_type}")
        return

    X = ct_adata.X.toarray() if sp.issparse(ct_adata.X) else ct_adata.X
    pca = PCA(n_components)
    x_emb = pca.fit_transform(X)
    p_x = torch.from_numpy(x_emb[y == 0])  # train positive
    u_tr_x = torch.from_numpy(x_emb[y == 1])  # train mixed unshifted
    u_te_x = torch.from_numpy(x_emb[y == 2])  # test mixed shifted
    data = TensorGroup(p_x=p_x, u_tr_x=u_tr_x, u_te_x=u_te_x)

    puc = PUcLearner(prior)
    puc.fit(data)

    # get predictions for each class
    pred_p_x = puc.predict(p_x)
    pred_u_tr_x = puc.predict(u_tr_x)
    pred_u_te_x = puc.predict(u_te_x)

    result = {
        'cell_type': cell_type,
        'argidx': argidx,
        'y': y,
        'p_x_pred': pred_p_x,
        'u_tr_x_pred': pred_u_tr_x,
        'u_te_x_pred': pred_u_te_x,
        'puc': puc,
        'prior': prior,
        'n_components': n_components,
        'age_thresh': age_thresh,
        'counts': counts,
        'genes': genes,
    }

    return result


def DE_test_ct(
    adata: anndata.AnnData,
    cell_type: str,
    result: Dict,
    known_markers: List[str] | None = None,
    min_n_snc_cells: int = 50,
    is_logged: bool = True,  # hlca is logged
):
    """Run DE tests via ranksums between senescent cells and healthy.

    Parameters
    ----------
    adata: Annotated Data object
    cell_type: string
        The cell type to consider in adata.
    result: Dict
        Dict containing PU learning results.
    known_markers: array-like
        If not None, will use only these genes.
    min_n_snc_cells: int
        Minimum number of senescent cells to use. If less, will skip.
    is_logged: bool
        If True, will assume data is logged. Needed for log2fc computation.
    """
    if known_markers is not None:
        ct_adata = adata[result['argidx'], np.in1d(adata.var['feature_name'], known_markers)]
    else:
        ct_adata = adata[result['argidx']]
    genes = ct_adata.var['feature_name'].to_numpy()

    y = result['y']  # age groups
    # pred == 1 if cell is in positive class (healthy) and -1 otherwise
    pred = result['u_te_x_pred']  # test distribution >=50

    # Skip cell type if there are fewer than min_n_snc_cells snc cells
    # or if the fraction of snc cells in train data is large (>30%) as
    # this may suggest a subtype rather than snc.
    if (pred != 1).sum() < min_n_snc_cells or (result['p_x_pred'] == 1).mean() < 0.7:
        print(f"Skipping {cell_type} with {(pred != 1).sum()} points")
        return

    old_adata = ct_adata[y == 2]  # take age group >=50
    X = old_adata.X.toarray() if sp.issparse(old_adata.X) else old_adata.X
    x_snc = X[pred != 1]
    x_healthy = X[pred == 1]

    statistic, pvals = ranksums(x_snc, x_healthy)
    pvals, qvals = pqvals(pvals)
    log2fc = _compute_log2fc(np.asarray(x_snc.mean(0)).ravel(),
                             np.asarray(x_healthy.mean(0)).ravel(),
                             is_logged=is_logged)
    test_result = {
        'cell_type': cell_type,
        'statistic': statistic,
        'pvals': pvals,
        'qvals': qvals,
        'log2fc': log2fc,
        'genes': genes,
    }
    return test_result


def get_top_markers(
    test_results,
    K: float = 0.0,
    Q: float = 0.05,
    use_statistic: bool = True,
    skip_cell_types: list | None = None,
):
    """Select SenSet genes based on the DE analysis.

    Parameters
    ----------
    test_results: Dict
        A dict mapping a cell type to a `test_result` from the function
        `DE_test_ct`.
    K: float
        The threshold to use for 'statistic' or 'log2fc' depending on the
        value of `use_statistic`. E.g., if this is 1 and
        `use_statistic=True`, will take all genes with abs(statistic) >= 1.
    Q: float
        The significance alpha for FDR.
    use_statistic: bool
        If True, will use the statistic from ranksum test to select genes.
        Otherwise, will use log2fc. If statistic > 0, will assume that the
        gene is enriched in senescent cells. Check order in ranksums.
    skip_cell_types: list, None
        If a list, then will skip the cell types in this list.
    """
    de_results = {}

    for cell_type, test_result in test_results.items():
        if test_result is None:
            continue
        if skip_cell_types is not None and cell_type in skip_cell_types:
            continue

        score = test_result['statistic' if use_statistic else 'log2fc']

        mask = (np.abs(score) >= K) & (test_result['qvals'] < Q)
        qvals = test_result['qvals'][mask]
        sig_statistic = score[mask]
        sig_genes = test_result['genes'][mask]

        mask_up = sig_statistic > 0
        mask_down = sig_statistic < 0

        de_results[cell_type] = {
            'all': sig_genes,
            'all_statistic': sig_statistic,
            'all_qvals': qvals,
            'up': sig_genes[mask_up],
            'up_statistic': sig_statistic[mask_up],
            'up_qvals': qvals[mask_up],
            'down': sig_genes[mask_down],
            'down_statistic': sig_statistic[mask_down],
            'down_qvals': qvals[mask_down],
        }

    return de_results
