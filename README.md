# scriabin
Single-cell resolution interaction analysis through binning

NB 23 June 2023: I am currently the only person maintaining this package and am in medical school rotations. I appreciate your patience as there may be times I respond slowly to issues that are opened. 

## Goal
Scriabin aims to provide a comprehensive view of cell-cell communication (CCC) at the single-cell level without requiring subsampling or aggregation. 

## Summary
Scriabin is a computational framework for analysis of cell-cell communication at single-cell resolution. Scriabin consists of 3 main workflows depending on dataset size and analytical goals: 1. the cell-cell interaction matrix workflow, optimal for smaller datasets, analyzes communication methods used for each cell-cell pair in the dataset; 2. the summarized interaction graph workflow, designed for large comparative analyses, identifies cell-cell pairs with different total communicative potential between samples; and, 3) the interaction program discovery workflow, suitable for any dataset size, finds modules of co-expressed ligand-receptor pairs. 

## Installation
Scriabin is implemented in R. To install: 
```
devtools::install_github("BlishLab/scriabin", ref = "main")
```

## Vignettes
Vignettes are available in the vignettes directory of this repo.
