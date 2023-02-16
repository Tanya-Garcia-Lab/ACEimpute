# clear workspace
rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model")

# # Run once: 
# install.packages("devtools")
# devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates/imputeCensoRd")
# devtools::install_github(repo = "Tanya-Garcia-Lab/ACEimpute/ACEimpute")

# install.packages("tidyverse")
# install.packages("lme4")
# install.packages("geex")
# install.packages("xtable")
# install.packages("latex2exp")

# Load packages
library(tidyverse)
library(lme4)
library(imputeCensoRd)
library(ACEimpute)
library(geex)

source("R_scripts/1_data_generation.R", echo=TRUE)

source("R_scripts/2_full_data_reml_analysis.R", echo=TRUE)

source("R_scripts/3_cmi_mi_reml_analysis.R", echo=TRUE)

source("R_scripts/4_ace_analysis.R", echo=TRUE)

source("R_scripts/5_make_tables.R", echo=TRUE)