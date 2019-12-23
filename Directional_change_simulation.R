#### Script to run directional change scenario

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

#### Import parameter values ####
# SURVIVAL

load('All_data_survival.RData')

# Extract parameter values
# order = 
# int, slope.z, slope.sync, slope.sync2, slope.spring p, slope.spring t, slope.winter t, slope.pop s, 
# slope.b
s.params_final <- extract_surv(model=final_survival, type ='tot_pop')

# DEVELOPMENT

load('All_data_development.RData')

# Extract parameter values
d.params_final <- extract_dev(model=final_development)

# RECRUITMENT

load('All_data_recruitment.RData')

# Extract parameter values
r.params_final <- extract_recru(final_recruitment)

# INHERITANCE

# STANDARD IPM (parent-offspring)
load('All_data_inheritance_STANDARD.RData') 


# loading all exhausts memory
load('All_data_inheritance_QG.RData')

# Extract parameter values
iqg.params_final <- extract_inher(final_inheritance_QG, method='QG')

ipo.params_final <- extract_inher(final_inheritance_PO, method='PO')

# CATERPILLAR TIMING

half_fall_params <- extract_HF(run_half_fall(test_years = 2020, datafile=datafile))

#### Initialisation ####

m.N_TP <- mean(descale$pop_size) # Re-assign mean of TOTAL population size to scale and unscale
sd.N_TP <- sd(descale$pop_size) 

# set up a vector of spring temperature
summary(lm(Spring_temp ~ Year, data = descale))
means <- seq(mean(descale$Spring_temp)+sd(descale$Spring_temp), 
             mean(descale$Spring_temp)+sd(descale$Spring_temp)+(0.04*50),
             length.out=51)
set.seed(10)
cross_val$s_temp <- rnorm(51, mean=means, sd=sd(descale$Spring_temp))

# Run Quant Gen IPM for 50 years
cross_val_RESULT_DIRECTIONAL <- run_CV_IPM(test_years = 1961:2010, s.params = s.params_final, 
                                  d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                  hf.params = half_fall_params, n = 100, n.gens = 50, g.range = c(-2,2), 
                                  m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                                  type = "tot_pop", method = "Predict",
                                  descale = descale)
# save result
save(cross_val_RESULT_DIRECTIONAL, file = "cross_val_RESULT_DIRECTIONAL.RData")

# Run Standard IPM for 50 years
cross_val_RESULT_STANDARD_DIRECTIONAL <- run_CV_IPM_PO(test_years = 1961:2010, s.params = s.params_final, 
                                     d.params = d.params_final, r.params = r.params_final, h.params = ipo.params_final,
                                     hf.params = half_fall_params, n = 10000, n.gens = 50, g.range = c(-4,4), 
                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                     CV_data = cross_val,
                                     descale = descale, method = "Predict")
# save result
save(cross_val_RESULT_STANDARD_DIRECTIONAL, file = "cross_val_RESULT_STANDARD_DIRECTIONAL.RData")
