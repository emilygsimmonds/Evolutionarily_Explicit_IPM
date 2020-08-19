#### Script to run projection scenarios

# COLLAPSE ALL FOR EASIER VIEWING

#### Packages ####

library(dplyr)
library(tidyr)
library(MASS)
library(lattice)
library(mvtnorm)
library(scales)
library(furrr)

#### Scripts/functions ####

source('EE_IPM.R')
source('Extract_param_values.R')
source('Demographic_functions.R')
source('Moments_function.R')

#### Data ####
# Total data
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled

# Climate projections

# These have same columns as cross validation data
load("low_predictions.RData") # change to your own
load("mid_predictions.RData") # change to your own
load("high_predictions.RData") # change to your own


#### Import parameter values ####
# SURVIVAL

load('All_data_survival.RData')

# Extract parameter values
# order = 
# int, slope.z, slope.sync, slope.sync2, slope.spring p, slope.spring t, slope.winter t, slope.pop s, 
# slope.b
s.params_final <- extract_surv(model=final_survival, type ='tot_pop')

# DEVELOPMENT

load('Final_development_model_maximum_temperature.RData') # using maximum temperature

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

load('caterpillar_timing_model.RData')

half_fall_params <- as.vector(coef(caterpillar_timing_model))

#### Initialisation ####

m.N_TP <- mean(descale$pop_size) # Re-assign mean of TOTAL population size to scale and unscale
sd.N_TP <- sd(descale$pop_size) 


#### Parallel running of simulations using furrr ####

plan("multiprocess")

low_results <- future_map(low_predictions, run_CV_IPM,
                          n = 100, n.gens = 91, g.range = c(-2,2), 
                          m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                          descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                          d.params = d.params_final, r.params = r.params_final, 
                          h.params = iqg.params_final,
                          hf.params = half_fall_params_max,
                          decoupled = TRUE, .progress=TRUE)

#save results
save(low_results, file="low_emissions_results.RData")

plan("multiprocess")
mid_results <- future_map(mid_predictions, run_CV_IPM,
                          n = 100, n.gens = 91, g.range = c(-2,2), 
                          m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                          descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                          d.params = d.params_final, r.params = r.params_final, 
                          h.params = iqg.params_final,
                          hf.params = half_fall_params_max,
                          decoupled = TRUE, .progress=TRUE)

#save results
save(mid_results, file="medium_emissions_results.RData")

plan("multiprocess")
high_results <- future_map(high_predictions, run_CV_IPM,
                           n = 100, n.gens = 91, g.range = c(-2,2), 
                           m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                           descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                           d.params = d.params_final, r.params = r.params_final, 
                           h.params = iqg.params_final,
                           hf.params = half_fall_params_max,
                           decoupled = TRUE, .progress=TRUE)

#save results
save(high_results, file="high_emissions_results.RData")

