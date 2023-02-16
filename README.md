# Mission Imputable: Correcting for Berkson Error When Imputing a Censored Covariate

This repository contains the data and scripts needed to reproduce results from the manuscript by Grosser, Lotspeich, and Garcia (2023+)

The `ACEimpute` package, which corrects for imputation error, can be found in this repo [here] (ACEimpute)

The `imputeCensRd` package, which implements the conditional mean imputation approaches from the paper, can be found in its own repo [here](https://github.com/sarahlotspeich/imputeCensRd).

Each of the "Script (Run Simulations)" files is coded to run 5 replication of each setting for demonstration purposes. In the manuscript listed above, all simulation settings were run with 1,000 simulations each.

## Tables 

**Table 1.** Simulation results comparing i) restricted maximum likelihood estimation (REML) with the full data, ii) conditional mean imputation (with a correctly specified imputaiton model) plus REML, and iii) conditional mean imputation (with a correctly specified imputaiton model) plus ACE imputation to correct for imputation error.

<!-- ![](Tables/Table1.png) -->

  - [Script (Run Simulations)](Manuscripts_Simulations/Correctly_Specified_Imputation_Model/R_scripts)
  - [Script (Make Figure)](Manuscripts_Simulations/Correctly_Specified_Imputation_Model/R_scripts/make_tables)
  - [Data (Simulation Results)](Manuscripts_Simulations/Correctly_Specified_Imputation_Model/sim_data)
