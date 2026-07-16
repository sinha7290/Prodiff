# Prodiff (CANDiT)

**Companion code for our *Cell Reports Medicine* 2025 paper introducing CANDiT, a machine learning framework for AI-guided differentiation therapy in colorectal cancer.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Publication](https://img.shields.io/badge/Published-Cell%20Reports%20Medicine%202025-blue)](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(25)00494-X)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.xcrm.2025.102421-informational)](https://doi.org/10.1016/j.xcrm.2025.102421)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter)](https://jupyter.org/)

## Overview

This repository contains the code and derived data underlying our study introducing **CANDiT (Cancer-Associated Nodes for Differentiation Targeting)**, a machine learning framework that identifies transcriptomic vulnerabilities for **differentiation therapy in colorectal cancer (CRC)**.

Differentiation therapy has transformed the management of some hematologic malignancies (most notably ATRA in acute promyelocytic leukemia) but has never translated to solid tumors, largely because of the intra and inter tumoral heterogeneity that obscures cancer stem cells (CSCs). CANDiT addresses this gap by:

* Centering on **CDX2**, a master intestinal lineage transcription factor and tumor suppressor that is lost in high risk, poorly differentiated CRCs
* Building a transcriptomic network anchored on CDX2 to nominate therapeutic strategies that reinstate its expression
* Prioritizing **PRKAB1**, a stress polarity sensor and regulatory subunit of AMPK, as a top druggable node
* Validating a clinical grade PRKAB1 agonist across three complementary models: CRC cell lines, xenografts in mice, and a prospective cohort of patient derived organoids (PDOs)

Across all three platforms, PRKAB1 agonism reactivates lineage commitment, dismantles Wnt / YAP driven stemness programs, and **selectively eliminates CDX2 low cancer stem cells** while sparing normal cells. In xenografts, this translated to a 68 percent reduction in tumor volume and elimination of mortality (hazard ratio 0.09). A **50 gene response signature** derived from integrated modeling across all platforms predicts approximately a 50 percent reduction in recurrence and mortality risk.

## Repository contents

### Analysis notebooks

| Notebook | Purpose |
| :--- | :--- |
| `PRODIFF_Fig4_A_G.ipynb` | Panels A through G of Figure 4 (in vivo and PDO efficacy) |
| `PRODIFF_Heatmap.ipynb` | Heatmap visualizations of the 50 gene response signature |
| `corr_plot_Fig3.ipynb` | Correlation analyses supporting Figure 3 |
| `multivariate analysis_Fig3E.ipynb` | Multivariate model for Figure 3E |
| `univariate analysis_Fig3F_5E.ipynb` | Univariate response analyses for Figures 3F and 5E |
| `prodiff_ROC_sig.ipynb` | ROC and AUC analyses of the response signature |

### CANDiT framework

* `CANDiT/` — the CANDiT machine learning framework implementation

### Data files

| File | Content |
| :--- | :--- |
| `HCT114.txt`, `SW480.txt` | Transcriptomic profiles from CRC cell lines |
| `organoid.txt` | Patient derived organoid (PDO) response data |
| `invivo.txt` | Xenograft in vivo response data |
| `IC50_heatmap.txt` | Dose response IC50 measurements |
| `ROC_PDO.txt` | ROC input for PDO response classification |
| `cell_xeno_deg.txt` | Differentially expressed genes across cell line and xenograft models |
| `pgsig-res-1.txt`, `pgsig-res-2.txt` | Signature scoring results |
| `node-2(1).txt`, `node-3(1).txt` | Network node outputs |
| `AM_MV3.txt` | Multivariate model inputs |
| `Gene signatures used to predict Rx synergy.xlsx` | Curated signatures for drug response prediction |
| `Signatures of PDO and PDX response in CRC.xlsx` | Curated PDO and PDX response signatures |

### Utilities

| File | Purpose |
| :--- | :--- |
| `bone.py` | Boolean Network Explorer (BoNE) helper functions |
| `StepMiner.py` | StepMiner tool for identifying step patterns in gene expression |
| `HegemonUtil.py` | Python utility for the Hegemon Boolean analysis framework |
| `hegemonutils.pl` | Perl companion utilities for Hegemon |
| `explore.conf` | Configuration file for the Hegemon exploration workflow |

## Reproducing the analysis

1. Clone the repository and install dependencies (Python 3.9+ with the scientific stack, plus `lifelines` for survival analyses, `scikit-learn` for ROC/AUC, and Perl for the Hegemon utility scripts).
2. Explore the CANDiT framework in `CANDiT/` to build or apply the network prioritization model.
3. Run the figure notebooks in numerical order (Fig3 → Fig4 → Fig5) to reproduce panels from the paper.
4. Use `prodiff_ROC_sig.ipynb` to evaluate signature performance in your own cohorts.

## Publication

> Sinha S, Alcantara J, Perry K, Castillo V, Ondersma AK, Banerjee S, McLaren E, Espinoza CR, Taheri S, Vidales E, Tindle C, Adel A, Amirfakhri S, Sawires JR, Yang J, Bouvet M, Ghosh P.
> **CANDiT: A machine learning framework for differentiation therapy in colorectal cancer.**
> *Cell Reports Medicine*, 2025.
> doi: [10.1016/j.xcrm.2025.102421](https://doi.org/10.1016/j.xcrm.2025.102421)
> PMCID: PMC12711679

If you use CANDiT, the 50 gene response signature, or any derived analyses in your own work, please cite the manuscript.

## Press coverage

* **Newswise**: *Precision reprogramming: How AI tricks cancer's toughest cells*
* **EurekAlert!**: *Precision reprogramming: How AI tricks cancer's toughest cells*
* **Mirage News**: *AI reprograms cancer cells with precision*

## Contact

**Saptarshi Sinha, Ph.D.**
Assistant Project Scientist, Department of Cellular and Molecular Medicine
Interim Director, Center for Precision Computational Systems Network (PreCSN)
University of California San Diego
Email: sasinha@health.ucsd.edu

## License

Released under the MIT License. See `LICENSE`.
