# CANDiT

Code and derived data for Sinha S, Alcantara J, Perry K, et al., "CANDiT: a
machine learning framework for differentiation therapy in colorectal cancer",
Cell Reports Medicine 6(11):102421, 2025.
https://doi.org/10.1016/j.xcrm.2025.102421

## What the study does

Differentiation therapy works in some haematologic malignancies but has not
translated to solid tumours, largely because tumour heterogeneity obscures the
stem-cell compartment. CANDiT approaches this by building a transcriptomic
network anchored on CDX2, an intestinal lineage transcription factor lost in
poorly differentiated colorectal cancers, and using it to nominate nodes whose
modulation reinstates lineage commitment.

The analysis prioritises PRKAB1, a regulatory subunit of AMPK, and the
prediction was then tested in CRC cell lines, mouse xenografts, and a
prospective cohort of patient-derived organoids. A 50-gene response signature
derived from those platforms is evaluated against clinical outcome data.

## Notebooks

- `univariate analysis_Fig3F_5E.ipynb` - univariate and multivariate regression
  against clinical covariates, coefficient plots, survival analysis
- `multivariate analysis_Fig3E.ipynb` - multivariate models on the AM cohort
- `prodiff_ROC_sig.ipynb` - ROC and AUC for the response signature
- `corr_plot_Fig3.ipynb` - correlation structure between signature components
- `PRODIFF_Fig4_A_G.ipynb` - panels for Figure 4
- `PRODIFF_Heatmap.ipynb` - expression heatmaps
- `CANDiT/Prodiff-Paper.ipynb` - the network construction and node ranking

## Data

Derived tables used by the notebooks are committed here: cell line and xenograft
differential expression (`cell_xeno_deg.txt`, `SW480.txt`, `HCT114.txt`),
organoid results (`organoid.txt`, `ROC_PDO.txt`), in vivo results
(`invivo.txt`), IC50 values (`IC50_heatmap.txt`), signature definitions (the two
xlsx files) and network output (`node-2(1).txt`, `node-3(1).txt`,
`pgsig-res-1.txt`, `pgsig-res-2.txt`).

Primary expression data is not included. Some notebook cells query a
lab-internal expression database and will not run outside that environment;
those calls are listed in MIGRATION.md.

## Setup

    pip install -r requirements.txt
    pip install git+https://github.com/sinha7290/bioutils.git

Shared helper functions (thresholding, regression tables, survival, plotting)
live in [bioutils](https://github.com/sinha7290/bioutils). This repository
previously carried copies of lab-internal scripts that only ran against one
filesystem; they have been removed. See MIGRATION.md for the mapping and for the
calls that still require the internal database.

## License

MIT
