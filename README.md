
## Introduction

These scripts were used to analyze data and make figures for the paper "Probing dermal immunity to mycobacteria through a controlled human infection model".

## Required packages

General data wrangling and plotting: 

`tidyverse, readxl, here, knitr, rstatix, scales, cowplot`
`ggpubr, ggbeeswarm, ggh4x, grDevices, gridExtra, svglite`

* CyTOF: `CATALYST, flowCore, diffcyt, FlowSOM`
* Single-cell: `Seurat, scDblFinder, SingleR, celldex, Azimuth, CellChat`
* Manuscript plots: `harmony, ComplexHeatmap`

## Input data

#### CyTOF

* Folder of FCS files: `'BCG Skin Biopsy 1 Live CD45+ cells FCS'/`
* Metadata: `2022_human_BCG_challenge_metadata_for_flowsom.xlsx`
* Panel: `2022_human_BCG_challenge_CyTOF_panel.xlsx`
* Manual annotations: `BCG_challenge_cluster_annotations_24.xlsx`

#### Single-cell

* Day 3 cellranger output: `filtered_feature_bc_matrix.h5`
* Day 15 cellranger output: `sample_filtered_feature_bc_matrix.h5`

#### Manuscript plots

* Processed CyTOF data: `BCG_Skin_Biopsy_CD45_Subsets.csv`
* Processed single-cell data: 
  * `8_annot_sub_final_d3.rds`
  * `8_annot_sub_final_d15.rds`
  * `9_final_annot_d3.rds`
  * `9_final_annot_d15.rds`
* Micriobiology data: `Copy of Combined CFU MVT RS data.xlsx`
* DEGs from Andrew Fiore-Gartland: `deg_heatmap_values.csv`

## Output

These scripts create an `output/` folder and subfolders within this repo if it doesn't exist.

* `plots/`: PDFs and PNGs of manuscript plots
* `processsed_data/`: RDS objects and CSVs of intermediate data

```
output/
├── plots/
└── processsed_data/
```

## Contact

Contact Emma Bishop with questions (emmab5@uw.edu) or file an issue.

