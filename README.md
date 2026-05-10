# SenSet

This repository contains the code used in the preprint "SenSet defines cell-type specific senescence signatures in the aged human lung" published in the [EMBO Journal](https://link.springer.com/article/10.1038/s44318-026-00762-8).

Senset is a 106-gene senescence signature for human lung tissue, derived by applying positive-unlabeled (PU) learning under covariate shift to the Human Lung Cell Atlas (HLCA). (see [notebooks/PUc_Learning.ipynb](notebooks/PUc_Learning.ipynb))

The full list of genes can be downloaded from [resources/SenSet.txt](resources/SenSet.txt).

### Derivation of SenSet

<p align="center">
  <img src="data/Figure1.png" width="600">
  <br>
  <em>Given a large single-cell lung cohort from young and old donors (A), we designate cells from young donors as labeled positives (non-senescent), whereas all cells from old donors are initially unlabeled (B). A covariate-shift-aware positive-unlabeled classifier (PUc) assigns senescence probabilities to the unlabeled cells and flags those with high scores as SnCs (C). Predicted SnCs are aggregated within each type to derive refined, cell-type-specific senescence-marker profiles (D).</em>
</p>

### Citation

```
@ARTICLE{Hasanaj2026-sq,
  title     = "{SenSet} defines cell-type specific senescence signatures in the
               aged human lung",
  author    = "Hasanaj, Euxhen and Beaulieu, Delphine and Wang, Cankun and Hu,
               Qianjiang and Rosas, Lorena and Bueno, Marta and Sembrat, John C
               and Pineda, Ricardo H and Melo-Narvaez, Maria Camila and
               Cardenes, Nayra and Yanwu, Zhao and Yingze, Zhang and Lafyatis,
               Robert and Morris, Alison and Mora, Ana and Rojas, Mauricio and
               Li, Dongmei and Rahman, Irfan and Pryhuber, Gloria S and Lehmann,
               Mareike and Alder, Jonathan and Gurkar, Aditi and Finkel, Toren
               and Ma, Qin and Lugo-Martinez, Jose and Póczos, Barnabás and
               Bar-Joseph, Ziv and Eickelberg, Oliver and Königshoff, Melanie",
  journal   = "EMBO J.",
  publisher = "Springer Science and Business Media LLC",
  pages     = "1--50",
  month     =  apr,
  year      =  2026,
  language  = "en"
}
```
