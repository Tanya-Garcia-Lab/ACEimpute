# Mission Imputable: Correcting for Berkson Error When Imputing a Censored Covariate

This repository contains the data and scripts needed to reproduce results from the manuscript by Grosser, Lotspeich, and Garcia (2023+). This manuscript introduces the novel method "actice correction for error in imputation" (ACE imputation), which adjusts for the errors incurred when imputing censored covariates.

# Installation 

The `ACEimpute` package, which corrects for imputation error, can be found in this repo [here](ACEimpute/). The `imputeCensoRd` package, which implements the conditional mean imputation approaches from the paper, can be found in its own repo [here](https://github.com/Tanya-Garcia-Lab/Imputing-Censored-Covariates/tree/main/imputeCensoRd).

```{r}
# Run once: install.packages("devtools")
devtools::install_github(repo = "Tanya-Garcia-Lab/ACEimpute/ACEimpute")
library(ACEimpute)
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates/imputeCensoRd")
library(imputeCensoRd)
```

## Tables 

The scripts in this repo are coded to run 5 replication of each simulation setting for demonstration purposes. In the manuscript listed above, all simulation settings were run with 1,000 simulations each.

**Tables 1 and 4.** Simulation results comparing i) restricted maximum likelihood estimation (REML) with the full data, ii) conditional mean imputation (with a correctly specified imputaiton model) plus REML, and iii) conditional mean imputation (with a correctly specified imputaiton model) plus ACE imputation to correct for imputation error.

  - [Script (run simulation)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/0_run_all.R)
  - [Script (make table)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/5_make_tables.R)
  - [Data (simulation replicates)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/sim_data)
  - [Data (simulation results)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/sim_results)
  - [Tables (manuscript Table 1 and Table 4)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/tables)
  
 **Tables 2 and 5.** Simulation results comparing i) restricted maximum likelihood estimation (REML) with the full data, ii) conditional mean imputation (with a misspecified imputaiton model) plus REML, and iii) conditional mean imputation (with a misspecified imputaiton model) plus ACE imputation to correct for imputation error. 

  - [Script (run simulation)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/0_run_all.R)
  - [Script (make table)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/5_make_tables.R)
  - [Data (simulation replicates)](Manuscript_Simulations/Misspecified_Imputation_Model/sim_data)
  - [Data (simulation results)](Manuscript_Simulations/Misspecified_Imputation_Model/sim_results)
  - [Tables (manuscript Table 2 and Table 5)](Manuscript_Simulations/Misspecified_Imputation_Model/tables)

## Figure

**Figure 2.** Plot power curves comparing sample sizes for a clinical trial based on estimates from i) complete case analysis, ii) conditional mean imputaiton plus REML, and iii) ACE imputation

  - [Script (plot power curves)](Manuscript_Simulations/Power_Curves/compare_power_curves.R)
  - [Figure](Manuscript_Simulations/Power_Curves/compare_power_curves.png)
