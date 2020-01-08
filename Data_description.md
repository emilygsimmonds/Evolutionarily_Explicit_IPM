---
title: "DATA DESCRIPTION"
output: html_document
---

## Description of datasrets held on this repository

## DOI: 

### Data files:

[bio_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the analyses. Rows are individual breeding attempts. All continuous variables are scaled (mean subtracted and divided by the standard deviation).

Columns are: 

* F_ID (the unique British Trust for Ornithology ring number for each female)
* Year (year of study)
* Nest_box (the ID of the nest box used for the breeding attempt) 
* Survival (binary variable indicating survival to next year (1) or not (0))
* Num_recruited (number of offspring produced from this breeding attempt that were caught as breeding adults in the population)
* April_hatch (date that first egg in the nest associated with this breeding attempt hatched. In days since April 1st - SCALED) 
* Synchrony (April_hatch - date of caterpillar timing, also in days since April 1st - SCALED)
* Clutch_size (Maximum number of eggs counted in the nest of this breeding attempt - SCALED)
* Age (binary variable of first year bird (1) or older (2))
* Immigrant (binary variable indicating if a bird was recorded as born in the woodland (1) or not (0))
* beech (factor indicating amount of beech mast in winter from none (0), little (1), to high (2))
* pop_size (number of unique breeding females counted in the breeding season - SCALED)
* Spring_temp (mean of mean daily temperatures (ºC) from 1st march to 9th April (Met Office, 2009, Hollis and McCarthy, 2017) - SCALED)
* winter_temp (mean of mean daily temperatures (ºC) from 1st December to 28/29th February Met Office, 2009; Hollis and McCarthy, 2017) - SCALED)
* winter_precip_t (sum of precipitation (mm) from 1st April to 31st May (Radcliffe Meterological Station, 2016) - SCALED)
* spring_precip_t (sum of precipitation (mm) from 1st December to 28/29th February (Radcliffe Meterological Station, 2016) - SCALED)
* Round (Section of the woodland based on habitat type (factor))
* Half_fall (date of peak caterpillar abundance in days since 1st April	Median date on which first instar winter moth larvae descend to ground to pupate - SCALED)
* R_pop_size (number of unique breeding females that were recorded as born in the woodland and counted in the breeding season - SCALED)

[descale.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the analyses. Rows are individual breeding attempts. Continuous variables are not scaled. All variables are the same as in bio_data.csv.

[bio_data_inheritance.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the inheritance analyses. Rows are individual breeding attempts. All continuous variables are scaled (mean subtracted and divided by the standard deviation). Restricted to birds that appear in the pedigree. Columns are the same as bio_data.csv with three additions.

Columns are: 

* MOTHER (the unique British Trust for Ornithology ring number for the mother of the focal female)
* animal (the unique British Trust for Ornithology ring number for the focal female)
* ID (the unique British Trust for Ornithology ring number for the focal female)

[bio_data_t.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the development analyses. Rows are individual breeding attempts. All continuous variables are scaled (mean subtracted and divided by the standard deviation). Restricted to birds that appear in at least two consecutive years in the data.

[bio_data_t1.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : biological data used in the development analyses. Rows are individual breeding attempts. All continuous variables are scaled (mean subtracted and divided by the standard deviation). Restricted to birds that appear in at least two consecutive years in the data. This datafile contains the data for the year following that in bio_data_t.csv.

Columns are: 

* F_ID (the unique British Trust for Ornithology ring number for each female)
* April_hatch_t1 (date that first egg in the nest associated with this breeding attempt hatched. In days since April 1st - SCALED - for time t+1)
* pop_size_t1 (number of unique breeding females counted in the breeding season - SCALED - for time t+1)
* Spring_temp_t1 (mean of mean daily temperatures (ºC) from 1st march to 9th April (Met Office, 2009, Hollis and McCarthy, 2017) - SCALED - for time t+1)
* winter_temp_t1 (mean of mean daily temperatures (ºC) from 1st December to 28/29th February Met Office, 2009; Hollis and McCarthy, 2017) - SCALED - for time t+1)
* winter_precip_t1 (sum of precipitation (mm) from 1st April to 31st May (Radcliffe Meterological Station, 2016) - SCALED - for time t+1)
* spring_precip_t1 (sum of precipitation (mm) from 1st December to 28/29th February (Radcliffe Meterological Station, 2016) - SCALED - for time t+1) 
* Round_t1 (Section of the woodland based on habitat type (factor) - for time t+1)
* R_pop_size_t1 (number of unique breeding females that were recorded as born in the woodland and counted in the breeding season - SCALED - for time t+1)

[column_names_indiv_time_varying.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : Vector of column names required for capture-mark-recapture analysis. Single column, row names indicate column names.

[cross_validation_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : data used to direct cross validation. Contains observed values for environmental and population drivers. No variables are scaled - all in raw units. 

Columns are: 

* Year (year of study)
* R_Mean_hatch (Annual mean hatch date. In days since April 1st, for resident birds only) 
* R_Min_hatch (Annual minimum hatch date. In days since April 1st, for resident birds only) 
* R_Max_hatch (Annual maximum hatch date. In days since April 1st, for resident birds only) 
* R_Var_hatch (Annual variance in hatch date. In days since April 1st, for resident birds only) 
* R_Clutch_size (Mean of the clutch sizes - for resident birds only)
* R_Var_clutch_size (Variance of the clutch sizes - for resident birds only)
* Beech (factor indicating amount of beech mast in winter from none (0), little (1), to high (2))
* Pop_size (number of unique breeding females counted in the breeding season)
* s_temp (mean of mean daily temperatures (ºC) from 1st march to 9th April (Met Office, 2009, Hollis and McCarthy, 2017))
* w_temp (mean of mean daily temperatures (ºC) from 1st December to 28/29th February Met Office, 2009; Hollis and McCarthy, 2017))
* w_precip (sum of precipitation (mm) from 1st April to 31st May (Radcliffe Meterological Station, 2016))
* s_precip (sum of precipitation (mm) from 1st December to 28/29th February (Radcliffe Meterological Station, 2016))
* half_fall (date of peak caterpillar abundance in days since 1st April	Median date on which first instar winter moth larvae descend to ground to pupate)
* R_pop_size (number of unique breeding females that were recorded as born in the woodland and counted in the breeding season)

[ped_data.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : pedigree data for the population, based on social information i.e. which parents birds were caught or observed at the nest. 

Columns are: 

* ID (the unique British Trust for Ornithology ring number for each female)
* MOTHER (the unique British Trust for Ornithology ring number for the father of each female)
* FATHER (the unique British Trust for Ornithology ring number for the mother of each female)

[year_variables.csv](https://github.com/emilygsimmonds/Cue_Identification/blob/master/bio_data.csv) : datafile containing year varying environmental drivers required for capture-mark-recapture analysis. Each row is a year and each column contains a different variable. Continuous variables are scaled. 

Columns are: 

* Year (year of study)
* Beech (factor indicating amount of beech mast in winter from none (0), little (1), to high (2))
* Pop_size (number of unique breeding females counted in the breeding season - SCALED)
* Spring_temp (mean of mean daily temperatures (ºC) from 1st march to 9th April (Met Office, 2009, Hollis and McCarthy, 2017) - SCALED)
* Winter_temp (mean of mean daily temperatures (ºC) from 1st December to 28/29th February Met Office, 2009; Hollis and McCarthy, 2017) - SCALED)
* Winter_precip (sum of precipitation (mm) from 1st April to 31st May (Radcliffe Meterological Station, 2016) - SCALED)
* Spring_precip (sum of precipitation (mm) from 1st December to 28/29th February (Radcliffe Meterological Station, 2016) - SCALED)
* R_pop_size (number of unique breeding females that were recorded as born in the woodland and counted in the breeding season - SCALED)


















