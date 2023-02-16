# Mission Imputable: Correcting for Berkson Error When Imputing a Censored Covariate

This repository contains the data and scripts needed to reproduce results from the manuscript by Grosser, Lotspeich, and Garcia (2023+). This manuscript introduces the novel method "actice correction for error in imputation" (ACE imputation), which adjusts for the errors incurred when imputing censored covariates.

The `ACEimpute` package, which corrects for imputation error, can be found in this repo [here](ACEimpute/)

The `imputeCensRd` package, which implements the conditional mean imputation approaches from the paper, can be found in its own repo [here](https://github.com/sarahlotspeich/imputeCensRd).

Each of the "Script X" files is coded to run 5 replication of each setting for demonstration purposes. In the manuscript listed above, all simulation settings were run with 1,000 simulations each.

## Tables 

**Tables 1 and 4.** Simulation results comparing i) restricted maximum likelihood estimation (REML) with the full data, ii) conditional mean imputation (with a correctly specified imputaiton model) plus REML, and iii) conditional mean imputation (with a correctly specified imputaiton model) plus ACE imputation to correct for imputation error.

  - [Script 1 (generate simulation data)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/1_data_generation.R)
  - [Script 2 (full data analysis)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/2_full_data_reml_analysis.R)
  - [Script 3 (conditional mean imputation + REML analysis)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/3_cmi_mi_reml_analysis.R)
  - [Script 4 (ACE imputation analysis)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/4_ace_analysis.R)
  - [Script 5 (make table)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/5_make_tables.R)
  - [Data (Simulation Results)](Manuscript_Simulations/Correctly_Specified_Imputation_Model/sim_data)
  
  
 **Tables 2 and 5.** Simulation results comparing i) restricted maximum likelihood estimation (REML) with the full data, ii) conditional mean imputation (with a misspecified imputaiton model) plus REML, and iii) conditional mean imputation (with a misspecified imputaiton model) plus ACE imputation to correct for imputation error.

  - [Script 1 (generate simulation data)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/1_data_generation.R)
  - [Script 2 (full data analysis)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/2_full_data_reml_analysis.R)
  - [Script 3 (conditional mean imputation + REML analysis)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/3_cmi_mi_reml_analysis.R)
  - [Script 4 (ACE imputation analysis)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/4_ace_analysis.R)
  - [Script 5 (make table)](Manuscript_Simulations/Misspecified_Imputation_Model/R_scripts/5_make_tables.R)
  - [Data (Simulation Results)](Manuscript_Simulations/Misspecified_Imputation_Model/sim_data)
