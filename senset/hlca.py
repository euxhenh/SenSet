import os
from functools import cached_property
from typing import List, Tuple

import anndata
import numpy as np
from numpy.typing import NDArray


class HLCA:
    """Human Lung Cell Atlas"""

    def __init__(
        self,
        hlca_path: str = 'data/HLCA.h5ad',
        backed: str | None = None,
        age_key: str = 'age_or_mean_of_age_range',
        remove_nan_age: bool = True,
        remove_smokers: bool = True,
    ):
        self.hlca_path = os.path.expanduser(hlca_path)
        self.age_key = age_key
        self._hlca = anndata.read_h5ad(self.hlca_path, backed)

        mask = np.zeros(self._hlca.shape[0]).astype(bool)

        if remove_nan_age:
            mask |= self._hlca.obs[age_key].isna()
            print(f"Removing {mask.sum()} cells with age='NaN'")
        if remove_smokers:
            prev = mask.sum()
            mask |= self._hlca.obs['smoking_status'] != 'never'
            print(f"Removing {mask.sum() - prev} cells from smokers")

        if mask.sum() > 0:
            self._hlca = self._hlca[~mask]

    @property
    def adata(self) -> anndata.AnnData:
        return self._hlca

    @cached_property
    def age(self) -> NDArray[np.float_]:
        """Return an array of ages per cell"""
        return self.adata.obs[self.age_key].to_numpy()

    @cached_property
    def donor_age(self) -> List[Tuple[str, float]]:
        """Return a list of (subject, age) pairs"""
        donors, indices = np.unique(self.adata.obs['donor_id'], return_index=True)
        ages = self.age[indices]
        return list(zip(donors, ages))

    @cached_property
    def unique_age(self) -> NDArray[np.float_]:
        """Return each subjects age"""
        return np.asarray([v[1] for v in self.donor_age])

    @cached_property
    def cell_type(self) -> NDArray[np.str_]:
        """Return each cell's type"""
        return self.adata.obs['cell_type'].to_numpy().astype(str)

    @cached_property
    def unique_cell_type(self) -> NDArray[np.str_]:
        """Return all unique cell types"""
        return np.unique(self.adata.obs['cell_type']).astype(str)

    @cached_property
    def genes(self) -> NDArray[np.str_]:
        """Return gene names"""
        return self.adata.var['feature_name'].to_numpy().astype(str)
