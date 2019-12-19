---
title: "README"
output: html_document
---

## Code for manuscript: Testing the effect of quantitative genetic inheritance in structured models on projections of population dynamics

## DOI: 

This repository contains:

### Data files:

[bio_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the analyses. Rows are years of study. Columns are; Year (year of study), hatch_date (annual population mean hatch date in YYYY-MM-DD format), MORE COLUMN NAMES - AND DETAILS OF EACH VARIABLE (INC. TABLE FROM PAPER).

### Scripts relating to statistical model selection for individual functions:

[Survival_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the survival function using a capture mark recapture analysis in R package 'rmarked'

[Survival_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing two wrapper functions to format survival data into capture histories and create design data

[Development_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the development function using a linear mixed effect model in R package 'lme4'

[Recruitment_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the recruitment function using a linear mixed effect model in R package 'lme4'

[Caterpillar_timing_analysis.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the caterpillar timing function using a linear model

[Inheritance_analysis_standard_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the inheritance function in the standard IPM using a linear mixed effect model in R package 'lme4'

[Inheritance_analysis_quantgen_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model selection for the inheritance function in the quantitative genetic IPM using an animal model using R package 'MCMCglmm'

[Statistical_model_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing wrapper functions to run each demographic function and the caterpillar timing function for all of the cross validation test/training dataset combinations

[Run_demographic_functions_cross_validation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run each demographic function and the caterpillar timing function using Statistical_model_functions.R, produces parameter values for each function for the cross validation datasets

### Scripts relating to the running of the integral projection models (standard and quantitative genetic):

[Demographic_functions.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing set up for all demographic functional forms

[EE_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing the function to run the evolutionarily explicity integral projection model 

[Standard_IPM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing the function to run the standard integral projection model 

[Cross_validation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run model cross validation for both the standard and the evolutionarily explicit integral projection models

[Directional_change_simulation.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to run the directional environmental change simulation for the evolutionarily explicit integral projection model

### Scripts to create Figures from the manuscript:

[Moment_function.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing function to extract four moments from model simulation outputs

[Reformat_data_for_plotting.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script containing a function to reformat model simulation outputs ready for plotting

[Figure_1.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 1

[Figure_2.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 2

[Figure_3.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot Figure 3

[Figure_SOM.R](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Script to plot supporting information figures and tables
