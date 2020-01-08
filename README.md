---
title: "README"
output: html_document
---

## Code for manuscript: Testing the effect of quantitative genetic inheritance in structured models on projections of population dynamics

## DOI: https://doi.org/10.1111/oik.06985

This repository contains:

### Data files:

Full description of all files and contents in [Data_description.md](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/Data_description.md).

[bio_data.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data.csv) : biological data used in the analyses. 

[descale.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/descale.csv) : biological data used in the analyses. Rows are individual breeding attempts. Continuous variables are not scaled. All variables are the same as in bio_data.csv.

[bio_data_inheritance.csv](https://github.com/emilygsimmonds/Evolutionarily_Explicity_IPM/blob/master/bio_data_inheritance.csv) : biological data used in the inheritance analyses. 

[bio_data_t.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the development analyses.

[bio_data_t1.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the development analyses. All variables for year t+1.

[column_names_indiv_time_varying.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Vector of column names required for capture-mark-recapture analysis. Single column, row names indicate column names.

[cross_validation_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : data used to direct cross validation. Contains observed values for environmental and population drivers. No variables are scaled - all in raw units. 

[ped_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : pedigree data for the population, based on social information i.e. which parents birds were caught or observed at the nest. 

[year_variables.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : datafile containing year varying environmental drivers required for capture-mark-recapture analysis. Each row is a year and each column contains a different variable. Continuous variables are scaled. 

### Scripts relating to statistical model selection for individual functions:

[Survival_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the survival function using a capture mark recapture analysis in R package 'rmarked'

[Survival_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing two wrapper functions to format survival data into capture histories and create design data

[Development_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the development function using a linear mixed effect model in R package 'lme4'

[Recruitment_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the recruitment function using a linear mixed effect model in R package 'lme4'

[Inheritance_analysis_standard.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the inheritance function in the standard IPM using a linear mixed effect model in R package 'lme4'

[Inheritance_analysis_quantgen.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the inheritance function in the quantitative genetic IPM using an animal model using R package 'MCMCglmm'

[Statistical_model_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing wrapper functions to run each demographic function and the caterpillar timing function for all of the cross validation test/training dataset combinations

[Run_demographic_functions_cross_validation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run each demographic function and the caterpillar timing function using Statistical_model_functions.R, produces parameter values for each function for the cross validation datasets

[Moment_function.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing function to extract four moments from model simulation outputs

### Scripts relating to the running of the integral projection models (standard and quantitative genetic):

[Demographic_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing set up for all demographic functional forms

[EE_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing the function to run the evolutionarily explicity integral projection model 

[Standard_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing the function to run the standard integral projection model 

[Run_cross_validation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model cross validation for both the standard and the evolutionarily explicit integral projection models

[Directional_change_simulation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run the directional environmental change simulation for the evolutionarily explicit integral projection model

### Scripts to create Figures from the manuscript:

[Reformat_data_for_plotting.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing a function to reformat model simulation outputs ready for plotting

[Figure_1.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 1

[Figure_2.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 2

[Figure_3.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 3
