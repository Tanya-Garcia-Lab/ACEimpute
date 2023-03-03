# clear workspace
rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Misspecified_Imputation_Model/")

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

# generate data for simulation
source("R_scripts/1_data_generation.R", echo=TRUE)

# analyze full (uncensored data) using restriced maximum likelihood estimation (REML)
source("R_scripts/2_full_data_reml_analysis.R", echo=TRUE)

# impute censored values using conditional mean imputaiton (CMI), then analyze imputed data with REML
source("R_scripts/3_cmi_mi_reml_analysis.R", echo=TRUE)

# impute censored values using CMI, then analyze imputed data using "active correction for error in imputation"
source("R_scripts/4_ace_analysis.R", echo=TRUE)

# create xtables of simulation results
# source("R_scripts/5_make_tables.R", echo=TRUE)