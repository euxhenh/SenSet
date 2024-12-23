from typing import Tuple

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def _correct(pvals, method='fdr_bh'):
    """Simple wrapper for multipletests."""
    return multipletests(
        pvals=pvals,
        alpha=0.05,
        method=method,
        is_sorted=False,
        returnsorted=False,
    )


def _compute_log2fc(mean1, mean2, base='e', is_logged=False):
    """Computes log2 fold change and converts base if data is already logged."""
    if is_logged:
        log2fc = mean1 - mean2
        # Convert base
        if base is not None and base != 2:
            base = np.e if base == 'e' else float(base)
            log2fc *= np.log2(base)
    else:
        log2fc = np.log2(mean1 + 1) - np.log2(mean2 + 1)
    return log2fc


def pqvals(pvals):
    pvals = pvals.copy()
    not_nan_mask = ~np.isnan(pvals)
    qvals = np.full_like(pvals, 1.0)
    qvals[not_nan_mask] = _correct(pvals[not_nan_mask])[1]
    pvals[~not_nan_mask] = 1.0
    return pvals, qvals


def construct_PU_results_table(
    results: dict, age_thresh: Tuple[float, float] = (30, 50)
) -> pd.DataFrame:
    """Combines PU results from all cell types into a DataFrame.
    """
    _table = []
    a1, a2 = age_thresh

    for cell_type, result in results.items():
        if result is None:
            continue
        mpl = (result['p_x_pred'] == 1).mean()
        mpu = (result['u_tr_x_pred'] == 1).mean()
        mnu = (result['u_te_x_pred'] == 1).mean()
        _table.append({
            'Cell Type': cell_type,
            f'N. Cells <{a1}': len(result['p_x_pred']),
            f'N. Cells {a1}-{a2}': len(result['u_tr_x_pred']),
            f'N. Cells >={a2}': len(result['u_te_x_pred']),
            f'Healthy Prop. <{a1}': mpl,
            f'Healthy Prop. {a1}-{a2}': mpu,
            f'Healthy Prop. >={a2}': mnu,
        })

    table_df = pd.DataFrame.from_dict(_table)

    table_df[f'N. Snc >={a2}'] = np.round(
        (1 - table_df[f'Healthy Prop. >={a2}']) * table_df[f'N. Cells >={a2}']
    ).astype(int)
    table_df[f'Snc. Prop. >={a2}'] = 1 - table_df[f'Healthy Prop. >={a2}']
    table_df[f'Snc. Prop. >={a2}'] = table_df[f'Snc. Prop. >={a2}'].map(
        lambda x: '{:.2g}'.format(x)).astype(float)  # 2 decimal points

    return table_df
