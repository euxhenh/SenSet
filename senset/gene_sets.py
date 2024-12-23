import json
from functools import cached_property, reduce
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy.typing import NDArray


class SncGeneSets:
    """Read senescence gene sets."""

    if TYPE_CHECKING:
        SenSet: NDArray[np.str_]

    def __init__(
        self, snc_markers_path: str = 'data/senescence_list.xlsx',
        remove_downregulated: bool = False,
        **marker_paths,
    ):
        self.snc_markers_path = snc_markers_path
        self.remove_downregulated = remove_downregulated
        self._go = pd.read_excel(snc_markers_path, sheet_name=0, index_col=0)
        self._fridman = pd.read_excel(snc_markers_path, sheet_name=1, index_col=0)
        self._senmayo = pd.read_excel(snc_markers_path, sheet_name=2, index_col=0)
        self._cellage = pd.read_excel(snc_markers_path, sheet_name=3, index_col=0)
        self._gene_sets = ['GO', 'Fridman', 'SenMayo', 'CellAge']

        self.keys = marker_paths.keys()
        for key, val in marker_paths.items():
            if val.endswith('.json'):
                with open(val, "r") as f:
                    markers = np.asarray(json.load(f)).astype(str)
            elif val.endswith('.csv'):
                markers = pd.read_csv(val, index_col=0, header=None)[1].to_numpy().ravel()
            elif val.endswith('.xlsx'):
                markers = pd.read_excel(val, header=None).to_numpy().ravel()
            elif val.endswith('.txt'):
                markers = open(val, "r").readlines()
                markers = np.asarray([m.strip() for m in markers])
            else:
                raise ValueError
            self._gene_sets.append(key)
            setattr(self, key, markers)

    def __repr__(self):
        s = (
            "Gene Set module with keys:\n"
            f"\tGO: {self.GO.size}\n"
            f"\tFridman: {self.Fridman.size}\n"
            f"\tSenMayo: {self.SenMayo.size}\n"
            f"\tCellAge: {self.CellAge.size}\n"
            f"\tunion: {self.union.size}"
        )

        if len(self.keys) > 0:
            s += "\n"
        for key in self.keys:
            s += f"\t{key}: {getattr(self, key).size}\n"
        return s

    @property
    def gene_sets(self):
        return self._gene_sets

    @staticmethod
    def _prep_symbols(symbols):
        symbols = np.char.upper(symbols.to_numpy().astype(str))
        symbols.sort()
        return symbols

    @cached_property
    def GO(self) -> NDArray[np.str_]:
        """Returns all genes in the GO term 0090398"""
        symbol = self._go['symbol']
        symbol = symbol[~symbol.isna()]
        return self._prep_symbols(symbol)

    @cached_property
    def Fridman(self) -> NDArray[np.str_]:
        """Returns all upregulated marker genes in Fridman et al."""
        if self.remove_downregulated:
            df = self._fridman[self._fridman['Regulation'] == 'UP']
        else:
            df = self._fridman
        return self._prep_symbols(df['symbol'])

    @cached_property
    def SenMayo(self) -> NDArray[np.str_]:
        """Returns all genes in the SenMayo list."""
        return self._prep_symbols(self._senmayo['symbol'])

    @cached_property
    def CellAge(self) -> NDArray[np.str_]:
        """Returns all upregulated genes in CellAge."""
        if self.remove_downregulated:
            df = self._cellage[self._cellage['Senescence Effect'] == 'Induces']
        else:
            df = self._cellage
        return self._prep_symbols(df['Symbol'])

    @cached_property
    def union(self) -> NDArray[np.str_]:
        """Return the union of all lists."""
        return reduce(np.union1d, [self.GO, self.Fridman, self.SenMayo, self.CellAge])
