---
title: "README"
output: html_document
---

## Code for manuscripts: 

## 1. Testing the effect of quantitative genetic inheritance in structured models on projections of population dynamics (https://doi.org/10.1111/oik.06985)

## 2. Phenological asynchrony: a ticking time-bomb for seemingly stable populations? ()

## DOI of data: [10.5281/zenodo.3601072](https://zenodo.org/record/3601072#.XhW8zRdKjUI)

This repository contains:

### Data files:

Full description of all files and contents in [Data_description.md](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Data_description.md).

[bio_data.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data.csv) : biological data used in the analyses, updated to include the separate cues for great tits and caterpillars. 

[descale.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/descale.csv) : biological data used in the analyses. Rows are individual breeding attempts. Continuous variables are not scaled. All variables are the same as in bio_data.csv, updated to include the separate cues for great tits and caterpillars.

[bio_data_inheritance.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data_inheritance.csv) : biological data used in the inheritance analyses. 

[bio_data_t.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data_t.csv) : biological data used in the development analyses.

[bio_data_t1.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data_t1.csv) : biological data used in the development analyses. All variables for year t+1.

[column_names_indiv_time_varying.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/column_names_indiv_time_varying.csv) : Vector of column names required for capture-mark-recapture analysis. Single column, row names indicate column names.

[cross_validation_data.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/cross_validation_data.csv) : data used to direct cross validation. Contains observed values for environmental and population drivers. No variables are scaled - all in raw units. 

[ped_data.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/ped_data.csv) : pedigree data for the population, based on social information i.e. which parents birds were caught or observed at the nest. 

[year_variables.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/year_variables.csv) : datafile containing year varying environmental drivers required for capture-mark-recapture analysis. Each row is a year and each column contains a different variable. Continuous variables are scaled. 

Environmental data and projections note: ©Crown Copyright 2009. The UK Climate Projections (UKCP09) have been made available by the Department for Environment, Food and Rural Affairs (Defra) and the Department of Climate Change (DECC) under licence from the Met
Office, UK Climate Impacts Programme, British Atmospheric Data Centre, Newcastle University, University of East Anglia,
Environment Agency, Tyndall Centre and Proudman Oceanographic Laboratory. These organisations give no warranties,
express or implied, as to the accuracy of the UKCP09 and do not accept any liability for loss or damage, which may arise
from reliance upon the UKCP09 and any use of the UKCP09 is undertaken entirely at the users risk.

### Scripts relating to statistical model selection for individual functions:

[Survival_analysis.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Survival_analysis.R) : Script to run model selection for the survival function using a capture mark recapture analysis in R package 'rmarked'

[Survival_functions.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Survival_functions.R) : Script containing two wrapper functions to format survival data into capture histories and create design data

[Development_analysis.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Development_analysis.R) : Script to run model selection for the development function using a linear mixed effect model in R package 'lme4'

[Recruitment_analysis.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Recruitment_analysis.R) : Script to run model selection for the recruitment function using a linear mixed effect model in R package 'lme4'

[Inheritance_analysis_standard.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Inheritance_analysis_standard.R) : Script to run model selection for the inheritance function in the standard IPM using a linear mixed effect model in R package 'lme4'

[Inheritance_analysis_quantgen.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Inheritance_analysis_quantgen.R) : Script to run model selection for the inheritance function in the quantitative genetic IPM using an animal model using R package 'MCMCglmm'

[Statistical_model_functions.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Statistical_model_functions.R) : Script containing wrapper functions to run each demographic function and the caterpillar timing function for all of the cross validation test/training dataset combinations

[Run_demographic_functions_cross_validation.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Run_demographic_functions_cross_validation.R) : Script to run each demographic function and the caterpillar timing function using Statistical_model_functions.R, produces parameter values for each function for the cross validation datasets

[Moment_function.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing function to extract four moments from model simulation outputs

[Extract_param_values.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Extract_param_values.R) : Script containing function to extract parameter values for each function from model object

### Scripts relating to the running of the integral projection models (standard and quantitative genetic):

[Demographic_functions.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Demographic_functions.R) : Script containing set up for all demographic functional forms

[EE_IPM.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/EE_IPM.R) : Script containing the function to run the evolutionarily explicity integral projection model 

[Standard_IPM.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Standard_IPM.R) : Script containing the function to run the standard integral projection model 

[Run_cross_validation.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Run_cross_validation.R) : Script to run model cross validation for both the standard and the evolutionarily explicit integral projection models

[Directional_change_simulation.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Directional_change_simulation.R) : Script to run the directional environmental change simulation for the evolutionarily explicit integral projection model

### Scripts to create Figures from the manuscript:

[Reformat_data_for_plotting.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Reformat_data_for_plotting.R) : Script containing a function to reformat model simulation outputs ready for plotting

[Figure_1_Oikos.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Figure_1_Oikos.R) : Script to plot Figure 1

[Figure_2_Oikos.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Figure_2_Oikos.R) : Script to plot Figure 2

[Figure_3_Oikos.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Figure_3_Oikos.R) : Script to plot Figure 3

[Figure_4_Oikos.R](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Figure_4_Oikos.R) : Script to plot Figure 4

### Model outputs:

Many of these are too large to store on GitHub but can be found on our Zenodo repository

### Specific scripts for 'Phenological asynchrony: a ticking time-bomb for seemingly stable populations?'

#### Model outputs

[Final_development_model_maximum_temperature.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/Final_development_model_maximum_temperature.RData) : the final model object for the development function, including the maximum temperature cue for great tits

[caterpillar_timing_model.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/caterpillar_timing_model.RData) : the final model object for the caterpillar timing function, including the caterpillar specific maximum temperature cue

#### Model scripts

[Predictive_model_code.R]() : script to run the IPM to create projections of future population size and structure 


#### Results

[Change_in_weather.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/Change_in_weather.csv) : Csv file containing a summary of weather changes over the projection period (2010-2100)

[Perturbations_results.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/Perturbations_results.csv) : Csv file of the results of the perturbation analysis

[high_emissions_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/high_emissions_results.RData) : Results of simulations using the high emissions scenario

[medium_emissions_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/medium_emissions_results.RData) : Results of simulations using the medium emissions scenario

[low_emissions_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/low_emissions_results.RData) : Results of simulations using the low emissions scenario

[only_ST_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/only_ST_results.RData) : Results of simulations with spring temperature held at 1961-2010 levels

[only_SP_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/only_SP_results.RData) : Results of simulations with spring precipitation held at 1961-2010 levels

[only_WP_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/only_WP_results.RData) : Results of simulations with winter precipitation held at 1961-2010 levels

[only_WT_results.RData](https://github.com/emilygsimmonds/Evolutionarily_Explicit_IPM/blob/master/only_WT_results.RData) : Results of simulations with winter temperature held at 1961-2010 levels

