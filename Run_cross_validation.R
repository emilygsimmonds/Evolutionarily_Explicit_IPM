#### Script to run cross validation

# COLLAPSE ALL FOR EASIER VIEWING

#### Packages ####

library(dplyr)
library(tidyr)
library(MASS)
library(lattice)
library(mvtnorm)
library(scales)

#### Scripts/functions ####

source('EE_IPM.R')
source('Standard_IPM.R')
source('Extract_param_values.R')
source('Demographic_functions.R')
source('Moments_function.R')

#### Data ####
# Total data
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled

# Observed data
cross_val <- read.csv("cross_validation_data.csv", header=T)

#### Initialisation ####

# list test dataset years
list_K_fold <- list(1961:1965, 1966:1970, 1971:1975, 1976:1980,
                    1981:1985, 1986:1990, 1991:1995, 1996:2000,
                    2001:2005, 2006:2010)

#### Import parameter values ####
# SURVIVAL

# resident pop size only
# in two parts due to large size
load('survival_resident_only_K1.RData')
load('survival_resident_only_K2.RData')
survival_RO <- c(model_outputs_full_K1, model_outputs_full_K2)

# tot pop
load('survival_tot_pop_K1.RData')
load('survival_tot_pop_K2.RData')
survival_TP <- c(model_outputs_tot_pop_K1, model_outputs_tot_pop_K2)

# Extract parameter values
# order = 
# int, slope.z, slope.sync, slope.sync2, slope.spring p, slope.spring t, slope.winter t, slope.pop s, 
# slope.b
s.params_RO <- mapply(extract_surv, model = survival_RO, 
                      MoreArgs = list(type="res_only"), SIMPLIFY = F)
s.params_TP <- mapply(extract_surv, model = survival_TP, 
                      MoreArgs = list(type="tot_pop"), SIMPLIFY = F)

# DEVELOPMENT

# resident pop size only
load('Development_resident_only_K.RData')

# total pop size
load('Development_tot_pop_K.RData')

# Extract parameter values
d.params_RO <- mapply(extract_dev, model = Dev_final_K, SIMPLIFY = F)
d.params_TP <- mapply(extract_dev, model = Dev_final_K_tot_pop, SIMPLIFY = F)

# RECRUITMENT

# resident pop size only
load('Recruitment_resident_only_K.RData')

# total pop size
load('Recruitment_tot_pop_K.RData')

# Extract parameter values
r.params_RO <- mapply(extract_recru, model = Recru_final_K, SIMPLIFY = F)

r.params_TP <- mapply(extract_recru, model = Recru_final_K_tot_pop, SIMPLIFY = F)

# INHERITANCE

# STANDARD IPM (parent-offspring)
load('Inheritance_STANDARD_K.RData')

# resident pop size only - QUANT GEN IPM
load('Inheritance_QG_resident_only_K.RData') 

# total pop size
load('Inheritance_QG_tot_pop_K.RData')

# Extract parameter values

ipo.params_TP <- mapply(extract_inher, model = Inher_PO_tot_pop, 
                        MoreArgs = list(method = "PO"), SIMPLIFY = F)

iqg.params_RO <- mapply(extract_inher, model = Inher_QG_full, 
                        MoreArgs = list(method = "QG"), SIMPLIFY = F)

iqg.params_TP <- mapply(extract_inher, model = Inher_QG_full_tot_pop, 
                        MoreArgs = list(method = "QG"), SIMPLIFY = F)

# CATERPILLAR TIMING

load("HF_output_K.RData")

# extract parameter values

HF.params <- mapply(extract_HF, model = HF_output_K, SIMPLIFY = F)

#### Run models - K-fold CROSS VALIDATION ####

m.N <- mean(descale$R_pop_size) # Assign mean of population size to scale and unscale
sd.N <- sd(descale$R_pop_size) # Assign sd of population size

# Run K-fold cross validation for resident only population size

### ••• SHOULD REMOVE POPULATION SIZE CORRECTION FOR IMMIGRANTS L259 in 'EE_IPM.R'

cross_val_RESULT_RO <- mapply(run_CV_IPM, test_years = list_K_fold, s.params = s.params_RO, 
                              d.params = d.params_RO, r.params = r.params_RO, h.params = iqg.params_RO,
                              hf.params = HF.params, MoreArgs = list(n = 100, n.gens = 6, g.range = c(-2,2), 
                                                                     m.N = m.N, sd.N = sd.N, CV_data = cross_val,
                                                                     type = "res_only", method = "CV",
                                                                     descale = descale), SIMPLIFY = F)
# save result
save(cross_val_RESULT_RO, file = "cross_val_RESULT_RO.RData")

m.N_TP <- mean(descale$pop_size) # Re-assign mean of TOTAL population size to scale and unscale
sd.N_TP <- sd(descale$pop_size) 

# set 1960 cross validation data to be the total pop size
cross_val$R_only_PS[1] <- cross_val$Pop_size[1]   

### ••• PUT BACK POPULATION SIZE CORRECTION FOR IMMIGRANTS L259 in 'EE_IPM.R'

# Run K-fold cross validation for total population size
cross_val_RESULT_TP <- mapply(run_CV_IPM, test_years = list_K_fold, s.params = s.params_TP, 
                              d.params = d.params_TP, r.params = r.params_TP, h.params = iqg.params_TP,
                              hf.params = HF.params, MoreArgs = list(n = 100, n.gens = 6, g.range = c(-2,2), 
                                                                     m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                                                                     type = "tot_pop", method = "CV",
                                                                     descale = descale), SIMPLIFY = F)
# save result
save(cross_val_RESULT_TP, file = "cross_val_RESULT_TP.RData")

# Run K-fold cross validation for standard IPM
cross_val_RESULT_STANDARD <- mapply(run_CV_IPM_PO, test_years = list_K_fold, s.params = s.params_TP, 
                              d.params = d.params_TP, r.params = r.params_TP, h.params = ipo.params_TP,
                              hf.params = HF.params, MoreArgs = list(n = 10000, n.gens = 6, g.range = c(-4,4), 
                                                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                                                     CV_data = cross_val,
                                                                     descale = descale,
                                                                     method="CV"), SIMPLIFY = F)
save(cross_val_RESULT_STANDARD, file = "cross_val_RESULT_STANDARD.RData")

#### Run models - 50-yr CROSS VALIDATION ####

m.N_TP <- mean(descale$pop_size) # Re-assign mean of TOTAL population size to scale and unscale
sd.N_TP <- sd(descale$pop_size) 

# Run for 50 years - EE IPM
cross_val_RESULT_50 <- run_CV_IPM(test_years = 1961:2009, s.params = s.params_final, 
                                  d.params = d.params_final, r.params = r.params_final, 
                                  h.params = iqg.params_final,
                                  hf.params = half_fall_params, n = 100, n.gens = 50, 
                                  g.range = c(-2,2), 
                                  m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                                  type = "tot_pop", method = "CV",
                                  descale = descale, decoupled=FALSE)
# save result
save(cross_val_RESULT_50, file = "cross_val_RESULT_50.RData")

# Run for 50 years - Standard IPM
cross_val_RESULT_STANDARD_50 <- run_CV_IPM_PO(test_years = 1961:2009, s.params = s.params_final, 
                                     d.params = d.params_final, r.params = r.params_final, h.params = ipo.params_final,
                                     hf.params = half_fall_params, n = 10000, n.gens = 50, g.range = c(-4,4), 
                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                     CV_data = cross_val,
                                     descale = descale, 
                                     method = "CV")
# save result
save(cross_val_RESULT_STANDARD_50, file = "cross_val_RESULT_STANDARD_50.RData")
