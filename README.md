# MW_effects_repli_and_exten

This project replicates the findings of **Dube, Lester, and Reich (2010)** regarding the effects of the minimum wage on earnings and employment and extends the analysis using modern causal inference methods for staggered treatment adoption.

## Replication & Extension Guide

To replicate the results of this study, please execute the scripts in the following order:

### Phase 1: Original Replication
These scripts reproduce the primary tables and specifications from the original Dube et al. (2010) paper using the traditional fixed effects framework.

1.  **`replicate_dube2010.R`**: This script handles the core data processing and foundational models required for the replication.

---

### Phase 2: Modern Extensions
After replicating the original results, these scripts apply state-of-the-art methods to address potential biases in staggered adoption designs (e.g., negative weights in TWFE).

2.  **Goodman-Bacon Decomposition**: 
    * Execute the scripts related to the Goodman-Bacon decomposition (**`gb.R`**) to visualize the variation contributing to the TWFE estimate and identify potential "forbidden comparisons."
3.  **Callaway & Sant’Anna (CS) Estimator**: 
    * Execute the CS-related scripts (**`cs2021.R`**) to obtain group-time average treatment effects ($ATT(g,t)$). This step accounts for heterogeneous treatment effects over time, which the original 2010 specifications may not fully capture.

---

## Requirements
* **Language**: R
* **Key Libraries**: `did` (for Callaway & Sant’Anna), `bacondecomp`, `fixest`, `tidyverse`.

## Data

To run the R scripts, you will need the following two datasets:

* `QCEWindustry_minwage_contig.dta`
* `QCEWindustry_minwage_all.dta`

Both files can be downloaded from the [Harvard Dataverse](https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/L4DUZ7/8ODET4&version=2.0).


