# clear workspace
rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model")

# # Run once: 
# install.packages("devtools")
# devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates/imputeCensoRd")
# devtools::install_github(repo = "Tanya-Garcia-Lab/ACEimpute/ACEimpute")

# Load packages
library(tidyverse)
library(lme4)
library(imputeCensoRd)
library(ACEimpute)
library(geex)

source("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/1_data_generation.R", echo=TRUE)

source("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/2_full_data_reml_analysis.R", echo=TRUE)

source("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/3_cmi_mi_reml_analysis.R", echo=TRUE)

source("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/4_ace_analysis.R", echo=TRUE)

source("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model/R_scripts/5_make_tables.R", echo=TRUE)