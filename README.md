# Neutrophil Retention Pipeline

This R script calculates **neutrophil retention time** and **interaction dynamics** with pre-neoplastic cells (PNCs) using 3D time-lapse imaging data.

## 🚀 Overview

This pipeline:
1. Calculates 3D distances between neutrophils and PNCs.
2. Assigns each neutrophil to its closest PNC (within a defined threshold).
3. Computes interaction durations (retention time) per cell pair.

Future versions will integrate:
- 3D segmentation
- Distance transformation-based analysis

## ⚙Key Assumptions

- **PNCs are spatially fixed** over time (`Time == 0` is used).
- The **distance threshold** for interaction is manully setted by measurement
- Neutrophils are matched with their **nearest PNC within threshold**.
- Interaction durations are based on **consecutive frames** of proximity.
- Time is recorded as **frame indices**; optional conversion to seconds available.

## Dependencies

R packages used:
- `tidyverse`
- `dplyr`
- `data.table`

Install them via:

```r
install.packages(c("tidyverse", "dplyr", "data.table"))
