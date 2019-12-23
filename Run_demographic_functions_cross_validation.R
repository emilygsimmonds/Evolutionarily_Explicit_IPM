#### Script to run statistical models for each demographic function 
# and the caterpillar timing function using Statistical_model_functions.R, 
# produces parameter values for each cross validation dataset

# COLLAPSE ALL FOR EASIER VIEWING

#### Packages ####

library(BaSTA)
library(lme4)
library(marked)
library(MCMCglmm)
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

#### Data ####
# Total data
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled

# Survival datasets
year_variables <- read.csv("year_variables.csv", header=T)
column_names <- read.csv("column_names_indiv_time_varying.csv", header=T)

# Development datasets 
datafile_t <- read.csv("bio_data_t.csv", header=T)
datafile_t1 <- read.csv("bio_data_t1.csv", header=T)
datafile_t <- cbind(datafile_t, datafile_t1)
datafile_t$Year <- datafile_t$Year + 1

# Inheritance datasets
datafile_I <- read.csv("bio_data_inheritance.csv", header=T)
# ensure correct structure of data
datafile_I$Age[which(datafile_I$Age > 1)] <- 2
datafile_I$Age <- as.factor(datafile_I$Age)
datafile_I$beech <- as.factor(datafile_I$beech)
datafile_I$Immigrant <- as.factor(datafile_I$Immigrant)
datafile_I$animal <- as.factor(datafile_I$animal)
ped <- read.csv("ped_data.csv", header=T)

# Observed data
cross_val <- read.csv("cross_validation_data.csv", header=T)

#### SURVIVAL ####

basta_input <- CensusToCaptHist(datafile$F_ID, datafile$Year, dformat = "%Y", timeInt = "Y")
# creates a data frame with row = ID and column = year, entries = capture history

# to use in 'marked' need to collapse capture history to single column and add ; at end
# want to do this for sections of 5 years
# need to save as text files are re-import

for(i in 1:10){
capture_history <- create_marked(test_years = list_K_fold[[i]], datafile = datafile,
basta_input = basta_input, full_years = 1960:2009)
filename <- paste("capture_history", i, ".txt", sep = "")
write.table(capture_history, file = filename)
}

# list the filenames to re-import below
list_filenames <- as.list(paste("capture_history", seq(1,10,1), ".txt", sep = ""))

library(RMark) # use RMark to import then detach to use marked for analyses
list_capture_history <- mapply(import.chdata, filename = list_filenames,
                               MoreArgs = list(field.types=c("f", "f", "f", rep("n", 225), "f"), 
                                               header=T), 
                               SIMPLIFY = F)

detach("package:RMark", unload=TRUE)

# Next process the data
list_processed <- mapply(process.data, data=list_capture_history, 
                       MoreArgs = list(model="CJS", nocc = 45, groups = c("Immigrant", "BYH", "Round"), 
                                       begin.time = 1960),
                      SIMPLIFY = F)
save(list_processed, file = "list_processed.RData")

# re-load saved list of processed data
load("list_processed.RData")

# create design data and run model
design.Phi = list(static = c("Immigrant", "Round"), time.varying = c("cs", "sn", "ah", "sns", "ahs"))
design.p = list(time.varying=c("ah"))
design.parameters <- list(Phi = design.Phi, p = design.p)
list_ddl <- mapply(create_design, processed_data = list_processed[1:5], test_years = list_K_fold[1:5],
                   MoreArgs = list(full_years = 1960:2009, year_variables = year_variables,
                                   parameters = design.parameters),
                   SIMPLIFY = F)

# need to complete in two different R sessions to preserve memory
list_ddl2 <- mapply(create_design, processed_data = list_processed[6:10], test_years = list_K_fold[6:10],
                    MoreArgs = list(full_years = 1960:2009, year_variables = year_variables,
                                    parameters = design.parameters),
                    SIMPLIFY = F)


# now have all of the data, just need to run the models

#### RUN Survival models K-fold ####
# for total population size and resident only population size
Phi.tot_pop = list(formula=~Round+Immigrant+ah+sns+sn+Pop_size+Beech+Spring_precip+Spring_temp+Winter_temp) 
Phi.res_only = list(formula=~Round+Immigrant+ah+sns+sn+R_pop_size+Beech+Spring_precip+Winter_temp) 
pars.tot_pop= list(Phi=Phi.tot_pop, p=p.full)
pars.res_only= list(Phi=Phi.res_only, p=p.full)

# First part of K-fold cross validation (data is too big for memory in one part)
# RESIDENT ONLY
model_outputs_full_K1 <- mapply(run_survival, data_processed = list_processed[1:5],
                                ddl = list_ddl,
                                MoreArgs = list(pars = pars.res_only), SIMPLIFY = F)

save(model_outputs_full_K2, file = "survival_resident_only_K1.RData")

# Second part of K-fold cross validation (data is too big for memory in one part)
model_outputs_full_K2 <- mapply(run_survival, data_processed = list_processed[6:10],
                                ddl = list_ddl2,
                                MoreArgs = list(pars = pars.res_only), SIMPLIFY = F)

save(model_outputs_full_K2, file = "survival_resident_only_K2.RData")

# First part of K-fold cross validation (data is too big for memory in one part)
# TOTAL POPULATION
model_outputs_tot_pop_K1 <- mapply(run_survival, data_processed = list_processed[1:5],
                                   ddl = list_ddl,
                                   MoreArgs = list(pars = pars.tot_pop), SIMPLIFY = F)

save(model_outputs_tot_pop_K1, file = "survival_tot_pop_K1.RData")

# Part 2
model_outputs_tot_pop_K2 <- mapply(run_survival, data_processed = list_processed[6:10],
                                   ddl = list_ddl2,
                                   MoreArgs = list(pars = pars.tot_pop), SIMPLIFY = F)

save(model_outputs_tot_pop_K2, file = "survival_tot_pop_K2.RData.RData")

#### DEVELOPMENT ####

# Run development model for resident only K-fold cross validation
Dev_final_K <- mapply(run_development, test_years = list_K_fold,
                      MoreArgs = list(datafile_t = datafile_t,
                                      datafile_t1 = datafile_t1, type = "res_only"),
                      SIMPLIFY = F)

save(Dev_final_K, file = "Development_resident_only_K.RData")

# Run development model for total pop size K-fold cross validation
Dev_final_K_tot_pop <- mapply(run_development, test_years = list_K_fold,
                              MoreArgs = list(datafile_t = datafile_t,
                                              datafile_t1 = datafile_t1, type = "tot_pop"),
                              SIMPLIFY = F)

save(Dev_final_K_tot_pop, file = "Development_tot_pop_K.RData")

#### RECRUITMENT ####

# Run recruitment models for resident only K-fold cross validation
Recru_final_K <- mapply(run_recruitment, test_years = list_K_fold,
                        MoreArgs = list(datafile, type = "res_only"), SIMPLIFY = F)

save(Recru_final_K, file="Recruitment_resident_only_K.RData")

# Run recruitment models for total pop size K-fold cross validation
Recru_final_K_tot_pop <- mapply(run_recruitment, test_years = list_K_fold,
                                MoreArgs = list(datafile, type = "tot_pop"), SIMPLIFY = F)

save(Recru_final_K, file="Recruitment_tot_pop_K.RData")

#### INHERITANCE QUANT GEN ####

# choose random effect structure - VA (G1), PE (G2), MO (G3), residual
# need to partition pheno variance as prior assumption
p.var <- var(datafile_I$April_hatch, na.rm=T) # variance of phenotype
# 0.16 split into Va
# V is the variance, n is the amount of confidence or weight of prior
prior.rand <- list(G=list(G1=list(V=matrix(p.var*0.16), n=1), 
                          G2=list(V=matrix(p.var*(0.84/3)), n=1),
                          G3=list(V=matrix(p.var*(0.84/3)), n=1)),
                   R=list(V=matrix(p.var*(0.84/3)), n=1))
datafile_I <- datafile_I[-which(is.na(datafile_I$MOTHER)==T),] # restrict to those in ped

# Run inheritance model for resident only K-fold cross validation
Inher_QG_full <- mapply(run_inheritance_QG, test_years = list_K_fold,
                        MoreArgs = list(datafile_I = datafile_I,
                                        ped = ped,
                                        type = "res_only", prior.rand=prior.rand), SIMPLIFY = F)

save(Inher_QG_full2, file = "Inheritance_QG_resident_only_K.RData")

# Run inheritance model for total pop size K-fold cross validation
Inher_QG_full_tot_pop <- mapply(run_inheritance_QG, test_years = list_K_fold,
                                MoreArgs = list(datafile_I = datafile_I,
                                                ped = ped,
                                                type = "tot_pop", prior.rand=prior.rand), SIMPLIFY = F)

save(Inher_QG_full_tot_pop, file = "Inheritance_QG_tot_pop_K.RData")

#### INHERITANCE STANDARD ####

# need to extract hatch date of mother and include
# create loop to pull mother hatch dates

mother_hatch <- rep(NA, length(datafile_I[,1]))
for(j in 1:length(datafile_I[,1])){
  temp <- filter(datafile_I, F_ID == datafile_I$F_ID[j])
  marker <- which(as.character(datafile_I$F_ID) 
                  == as.character(temp$MOTHER[1]))
  ifelse(length(marker) == 0, mother_hatch[j] <- NA,
         mother_hatch[j] <- mean(datafile_I$April_hatch[marker]))
}

datafile_I$Mother_hatch <- mother_hatch
datafile_I2 <- datafile_I[-which(is.na(datafile_I$Mother_hatch)),]
datafile_I2 <- filter(datafile_I2, Immigrant == 0)

# Run inheritance model for K-fold cross validation
# For standard IPM, only use total population size
Inher_PO_tot_pop <- mapply(run_inheritance_PO, test_years = list_K_fold,
                           MoreArgs = list(datafile_I = datafile_I2,
                                           type = "tot_pop"))

save(Inher_PO_tot_pop, file = "Inheritance_STANDARD_K.RData")


#### CATERPILLAR TIMING ####

# interpolated using everything 
# extrapolate using fewer covariates (i.e. no beech - held at 0)

# Run model for catperillar timing for all training datasets for K-fold cross validation
HF_output_K <- mapply(run_half_fall, test_years = list_K_fold,
                      MoreArgs = list(datafile), SIMPLIFY = F)

save(HF_output_K, file = "HF_output_K.RData")
