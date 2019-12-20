# Script to run statistical analyses to capture fundamental processes

# Beginning: 14.6.18
# Author: EGS

# Contains:
# - recapture analysis
# - survival analysis
# - recruitment analysis
# - development analysis
# - inheritance analysis
#-----------------------------------------------------------------------------

#### Packages and data ####
#datafile <- read.csv("updated_model_input_11.10.18.csv", header=T) # now scaled correctly
datafile <- read.csv("updated_model_input_MAX.csv", header=T) 
head(datafile)
summary(datafile)
#descale <- read.csv("updated_descale_11.10.18.csv", header=T)
descale <- read.csv("updated_descale_MAX.csv", header=T)
# Development datasets 
#datafile_t <- read.csv("datafile_t_11.10.18.csv", header=T)
#datafile_t1 <- read.csv("datafile_t1_11.10.18.csv", header=T)
datafile_t <- read.csv("datafile_t_MAX.csv", header=T)
datafile_t1 <- read.csv("datafile_t1_MAX.csv", header=T)
datafile_t <- cbind(datafile_t, datafile_t1)
datafile_t$Year <- datafile_t$Year + 1
year_variables <- read.csv("year_variables_Jan19.csv", header=T)
# Inheritance dataset
#datafile_I <- read.csv("datafile_I_11.10.18.csv", header=T)
datafile_I <- read.csv("datafile_I_MAX.csv", header=T)
# need to have an animal and ID and a BYEAR column.
datafile_I$Age[which(datafile_I$Age > 1)] <- 2
datafile_I$Age <- as.factor(datafile_I$Age)
datafile_I$beech <- as.factor(datafile_I$beech)
datafile_I$Immigrant <- as.factor(datafile_I$Immigrant)
datafile_I$animal <- as.factor(datafile_I$animal)
ped <- read.csv("ped_inheritance_11.10.16.csv", header=T)
cross_val <- read.csv("cross_validation_data_1961to2010_19.10.csv", header=T)

# K-fold cross validation years
list_K_fold <- list(1961:1965, 1966:1970, 1971:1975, 1976:1980,
                    1981:1985, 1986:1990, 1991:1995, 1996:2000,
                    2001:2005, 2006:2010)

#list_K_fold <- list(1961:1970,1971:1980,1981:1990,1991:2000,2001:2010)

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

source("Functions_IPM1.R")

#half_fall_params <- extract_HF(run_half_fall(test_years = 2020, datafile=datafile))
half_fall_params_max <- as.vector(coef(lm(Half_fall ~ Spring_temp_cat + spring_precip_t + winter_temp + 
                            winter_precip_t + as.factor(beech), data=datafile)))

#### SURVIVAL ####
# to do this need to perform analysis in 'marked'
# here we want to restrict to create the capture history using 'BaSTA'
basta_input <- CensusToCaptHist(datafile$F_ID, datafile$Year, dformat = "%Y", timeInt = "Y")
# creates a data frame with row = ID and column = year, entries = capture history

# to use in 'marked' need to collapse capture history to single column and add ; at end
# want to do this for sections of 5 years
# need to save as text files are re-import

#for(i in 1:10){
#capture_history <- create_marked(test_years = list_K_fold[[i]], datafile = datafile,
#basta_input = basta_input, full_years = 1960:2009)
#filename <- paste("capture_history", i, ".txt", sep = "")
#write.table(capture_history, file = filename)
#}

list_filenames <- as.list(paste("capture_history", seq(1,10,1), ".txt", sep = ""))


library(RMark)
list_capture_history <- mapply(import.chdata, filename = list_filenames,
                               MoreArgs = list(field.types=c("f", "f", "f", rep("n", 225), "f"), 
                                               header=T), 
                               SIMPLIFY = F)

detach("package:RMark", unload=TRUE)

# now have a list of capture histories
# need to do the next step for all variables so that can do it for null and full models

# can process the data
#list_processed <- mapply(process.data, data=list_capture_history, 
 #                        MoreArgs = list(model="CJS", nocc = 45, groups = c("Immigrant", "BYH", "Round"), 
  #                                       begin.time = 1960),
   #                      SIMPLIFY = F)
#save(list_processed, file = "list_processed.RData")
load("list_processed.RData")

# create design data and run
design.Phi = list(static = c("Immigrant", "BYH", "Round"), time.varying = c("cs", "sn", "ah", "sns", "ahs"))
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

# from nulls and full to final

# NULLS
p.full = list(formula=~time)

# FULL
Phi.tot_pop = list(formula=~Round+Immigrant+ah+sns+sn+Pop_size+Beech+Spring_precip+Spring_temp+Winter_temp) 
Phi.res_only = list(formula=~Round+Immigrant+ah+sns+sn+R_pop_size+Beech+Spring_precip+Winter_temp) 
pars.tot_pop= list(Phi=Phi.tot_pop, p=p.full)
pars.res_only= list(Phi=Phi.res_only, p=p.full)


# Now need to run all different parameterisations on full dataset
# THEN - run full model and nulls for the K fold cross validation

# actual survival analysis
model_outputs_full_K1 <- mapply(run_survival, data_processed = list_processed[1:5],
                               ddl = list_ddl,
                               MoreArgs = list(pars = pars.res_only), SIMPLIFY = F)

# 7 and full are the final pars
model_outputs_full_K2 <- mapply(run_survival, data_processed = list_processed[6:10],
                        ddl = list_ddl2,
                        MoreArgs = list(pars = pars.res_only), SIMPLIFY = F)
#summary(model_outputs_full_K[[5]]$results$beta)
#save(model_outputs_full_K2, file = "model_outputs_full_K2.RData")

# actual survival analysis
model_outputs_tot_pop_K1 <- mapply(run_survival, data_processed = list_processed[1:5],
                                ddl = list_ddl,
                                MoreArgs = list(pars = pars.tot_pop), SIMPLIFY = F)

# 7 and full are the final pars
model_outputs_tot_pop_K2 <- mapply(run_survival, data_processed = list_processed[6:10],
                                ddl = list_ddl2,
                                MoreArgs = list(pars = pars.tot_pop), SIMPLIFY = F)

#summary(model_outputs_tot_pop_K[[5]]$results$beta)
save(model_outputs_tot_pop_K1, file = "model_outputs_tot_pop_K1.RData")
#summary(model_outputs_tot_pop_K2[[5]]$results$beta)
save(model_outputs_tot_pop_K2, file = "model_outputs_tot_pop_K2.RData")

### run full model
capture_history <- create_marked(test_years = NULL, datafile = datafile, basta_input = basta_input, full_years = 1960:2009, 
                                 type="Full")
write.table(capture_history, file = "capture_history_full.txt")
library(RMark)
capture_history <- import.chdata(filename = "capture_history_full.txt", 
                                      field.types=c("f", "f", rep("n", 250), "f"), 
                                               header=T)

detach("package:RMark", unload=TRUE)

data_processed <- process.data(capture_history, model="CJS", accumulate = F, begin.time = 1960, 
                               groups = c("Immigrant", "Round"), nocc = 50)

# Create design data with static and time-varying covariates
design.Phi = list(static = c("Immigrant", "Round"), time.varying = c("cs", "sn", "ah", "sns", "ahs"))
design.p = list(time.varying=c("ah"))
design.parameters <- list(Phi = design.Phi, p = design.p)
ddl <- make.design.data(data_processed, parameters = design.parameters)
year_variables$time <- 1960:2009
ddl$Phi <- merge_design.covariates(ddl$Phi, year_variables, bytime=TRUE)
ddl$Phi$Beech <- as.factor(ddl$Phi$Beech)
year_variables2 <- year_variables
year_variables2$time[1:50] <- 1961:2010
ddl$p <- merge_design.covariates(ddl$p, year_variables2, bytime=TRUE)

Phi.1 = list(formula=~1)
p.full = list(formula=~time)
# FULL
Phi.full = list(formula=~Immigrant+Round+ah+ahs+sns+sn+cs+Pop_size+Beech+Spring_precip+Winter_precip+Spring_temp+Winter_temp)

null_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.1, p=p.full), accumulate = F)
full_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.full, p=p.full), accumulate = F)

AIC_store <- data.frame(AIC=rep(NA,5), delta_AIC=rep(NA,5), model_name=rep(NA,5))
AIC_store$AIC <- as.numeric(AIC_store$AIC)

AIC_store[1,1] <- null_model$results$AIC
AIC_store[1,2] <- AIC_store$AIC[1] - AIC_store$AIC[1]
AIC_store[1,3] <- "null"
AIC_store[2,1] <- full_model$results$AIC
AIC_store[2,2] <- AIC_store$AIC[2] - AIC_store$AIC[1]
AIC_store[2,3] <- "full"
AIC_store

removal_matrix <- matrix(NA, ncol=4, nrow=100)
removal_matrix[1:length(coef(full_model)$Estimate),1] <- row.names(coef(full_model))
removal_matrix[1:length(coef(full_model)$Estimate),2] <- coef(full_model)$Estimate
removal_matrix[1:length(coef(full_model)$Estimate),3] <- 2*coef(full_model)$se
removal_matrix[1:length(coef(full_model)$Estimate),4] <- (abs(coef(full_model)$Estimate) - 2*coef(full_model)$se)/abs(coef(full_model)$Estimate)

x2 <- removal_matrix[which(removal_matrix[,4] < 0), ] 

# REMOVE: clutch size
Phi.2 = list(formula=~Immigrant+Round+ah+ahs+sns+sn+Pop_size+Beech+Spring_precip+Winter_precip+Spring_temp+Winter_temp)

red1_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.2, p=p.full), accumulate = F)

AIC_store[3,1] <- red1_model$results$AIC
AIC_store[3,2] <- AIC_store$AIC[3] - AIC_store$AIC[2]
AIC_store[3,3] <- "- clutch size"
AIC_store

removal_matrix <- matrix(NA, ncol=4, nrow=100)
removal_matrix[1:length(coef(red1_model)$Estimate),1] <- row.names(coef(red1_model))
removal_matrix[1:length(coef(red1_model)$Estimate),2] <- coef(red1_model)$Estimate
removal_matrix[1:length(coef(red1_model)$Estimate),3] <- 2*coef(red1_model)$se
removal_matrix[1:length(coef(red1_model)$Estimate),4] <- (abs(coef(red1_model)$Estimate) - 2*coef(red1_model)$se)/abs(coef(red1_model)$Estimate)

x2 <- removal_matrix[which(removal_matrix[,4] < 0), ] 

# REMOVE: winter precip
Phi.3 = list(formula=~Immigrant+Round+ah+ahs+sns+sn+Pop_size+Beech+Spring_precip+Spring_temp+Winter_temp)

red2_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.3, p=p.full), accumulate = F)

AIC_store[4,1] <- red2_model$results$AIC
AIC_store[4,2] <- AIC_store$AIC[4] - AIC_store$AIC[3]
AIC_store[4,3] <- "- winter precip"
AIC_store

coef(red2_model)

# REMOVE: quadratic hatch
Phi.4 = list(formula=~Immigrant+Round+ah+sns+sn+Pop_size+Beech+Spring_precip+Spring_temp+Winter_temp)

red3_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.4, p=p.full), accumulate = F)

AIC_store[5,1] <- red3_model$results$AIC
AIC_store[5,2] <- AIC_store$AIC[5] - AIC_store$AIC[4]
AIC_store[5,3] <- "- quadratic hatch"
AIC_store

coef(red3_model)

removal_matrix <- matrix(NA, ncol=4, nrow=100)
removal_matrix[1:length(coef(red3_model)$Estimate),1] <- row.names(coef(red3_model))
removal_matrix[1:length(coef(red3_model)$Estimate),2] <- coef(red3_model)$Estimate
removal_matrix[1:length(coef(red3_model)$Estimate),3] <- 2*coef(red3_model)$se
removal_matrix[1:length(coef(red3_model)$Estimate),4] <- round((abs(coef(red3_model)$Estimate) - 2*coef(red3_model)$se)/abs(coef(red3_model)$Estimate),2)

x2 <- removal_matrix[which(removal_matrix[,4] < 0), ] 

# REMOVE: Spring temp
Phi.5 = list(formula=~Immigrant+Round+ah+sns+sn+Pop_size+Beech+Spring_precip+Winter_temp)

red4_model <- crm(data_processed, ddl, hessian = T, model="CJS", model.parameters=list(Phi=Phi.5, p=p.full), accumulate = F)

AIC_store[6,1] <- red4_model$results$AIC
AIC_store[6,2] <- AIC_store$AIC[6] - AIC_store$AIC[5]
AIC_store[6,3] <- "- spring temp"
AIC_store

coef(red4_model)
coef(final_survival)

AIC_store[2:6,2] <- round(AIC_store[1:5,1] - AIC_store[2:6,1],2)
write.csv(AIC_store, "AIC_output_survival_APR19.csv")

final_survival <- run_survival(data_processed = data_processed, ddl = ddl, pars = pars.tot_pop)
final_survival$results$beta$Phi
save(final_survival, file="Final_survival_Feb19.RData")
write.csv(coef(final_survival), "Coefs_survival_APR19.csv")

#### RECRUITMENT ####

# Here just want to automate the K-fold cross validation
# so use different data subsets in same models - save whole model and extract params later

# FINAL
Recru_final_K <- mapply(run_recruitment, test_years = list_K_fold,
                     MoreArgs = list(datafile, type = "res_only"), SIMPLIFY = F)
Recru_final_K_tot_pop <- mapply(run_recruitment, test_years = list_K_fold,
                        MoreArgs = list(datafile, type = "tot_pop"), SIMPLIFY = F)

save(Recru_final_K, file="Recru_final_K_Feb19.RData")

save(Recru_final_K_tot_pop, file="Recru_final_K_TP_Feb19.RData")

### run full model

final_recruitment <- run_recruitment(test_years = NULL, datafile, type = "tot_pop")
save(final_recruitment, file="Final_recruitment_Feb19.RData")

#### DEVELOPMENT ####

# Here just want to automate the K-fold cross validation
# so use different data subsets in same models - just save parameter values

Dev_final_K <- mapply(run_development, test_years = list_K_fold,
                   MoreArgs = list(datafile_t = datafile_t,
                                   datafile_t1 = datafile_t1, type = "res_only"),
                   SIMPLIFY = F)

Dev_final_K_tot_pop <- mapply(run_development, test_years = list_K_fold,
                      MoreArgs = list(datafile_t = datafile_t,
                                      datafile_t1 = datafile_t1, type = "tot_pop"),
                      SIMPLIFY = F)

save(Dev_final_K, file = "Dev_final_K_11.10.RData")

save(Dev_final_K_tot_pop, file = "Dev_final_K_TP_Feb19.RData")

### run full model

final_development <- run_development(test_years = NULL, datafile_t = datafile_t, datafile_t1 = datafile_t1, type = "tot_pop")
save(final_development, file="Final_development_Feb19.RData")

#### INHERITANCE QUANT ####
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

Inher_QG_full <- mapply(run_inheritance_QG, test_years = list_K_fold,
                      MoreArgs = list(datafile_I = datafile_I,
                                      ped = ped,
                                      type = "res_only", prior.rand=prior.rand), SIMPLIFY = F)

Inher_QG_full_tot_pop <- mapply(run_inheritance_QG, test_years = list_K_fold,
                        MoreArgs = list(datafile_I = datafile_I,
                                        ped = ped,
                                        type = "tot_pop", prior.rand=prior.rand), SIMPLIFY = F)

save(Inher_QG_full2, file = "Inher_QG_full.RData")

save(Inher_QG_full_tot_pop, file = "Inher_QG_K_tot_pop_Feb19.RData")

### run full model

final_inheritance_QG <- run_inheritance_QG(test_years = NULL, datafile_I = datafile_I, ped = ped,
                                               type = "tot_pop", prior.rand=prior.rand)
save(final_inheritance_QG, file = "Final_inheritance_QG_Feb19.RData")

# testing inclusion of WT - crosses 0 so remove
colnames(ped) <- c("animal", "FATHER", "MOTHER")
model <- MCMCglmm(April_hatch ~ Mother_hatch + Spring_temp +
                    winter_precip_t + spring_precip_t + Round + pop_size, random=~animal+ID+MOTHER, pedigree=ped, 
                  data=datafile_I2, prior=prior.rand2, verbose=F, nitt=1500, burnin=100, thin=10)
summary(model)
#save(model, "Inher_all_data.RData")

#### INHERITANCE NORM ####

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

Inher_PO_tot_pop <- mapply(run_inheritance_PO, test_years = list_K_fold,
                      MoreArgs = list(datafile_I = datafile_I2,
                                      type = "tot_pop"))
save(Inher_PO_tot_pop, file = "Inher_PO_tot_pop.RData")

Inher_PO_res_only <- mapply(run_inheritance_PO, test_years = list_K_fold,
                           MoreArgs = list(datafile_I = datafile_I2,
                                           type = "res_only"))
save(Inher_PO_res_only, file = "Inher_PO_res_only.RData")

### run full model

final_inheritance_PO <- run_inheritance_PO(test_years = NULL, datafile_I = datafile_I2,
                                           type = "tot_pop")
save(final_inheritance_PO, file = "Final_inheritance_PO_Feb19.RData")

# testing PO
model <- lmer(April_hatch ~ Mother_hatch + pop_size + Spring_temp +
                winter_temp + spring_precip_t + winter_precip_t + Round +
                (1|Year), data=datafile_I2) # AIC = 8317.371 new = 8409.131
summary(model)

# - WT
model2 <- lmer(April_hatch ~ Mother_hatch + pop_size + Spring_temp +
                spring_precip_t + winter_precip_t + Round +
                (1|Year), data=datafile_I2) # AIC = 8311.812 new = 8403.311
summary(model2)

# - SP
model3 <- lmer(April_hatch ~ Mother_hatch + pop_size  + Spring_temp +
                 spring_precip_t + Round +
                (1|Year), data=datafile_I2) # AIC = 8306.746 new = 8397.863
summary(model3)

# - WP
model4 <- lmer(April_hatch ~ Mother_hatch + pop_size  + Spring_temp +
                 Round +
                 (1|Year), data=datafile_I2) # AIC = 8393.28
summary(model4)

null <- lmer(April_hatch ~ 
               (1|Year), data=datafile_I2) # AIC = 8393.28


# write output AIC for PO
AIC <- data.frame(AIC = c(AIC(null), AIC(model), AIC(model2), AIC(model3), AIC(model4)),
                  delta = rep(NA,5),
                  removal = c("Null", 'Full', "WT", "WP", "WP"))
AIC$delta <- round(c(0,AIC[2,1]-AIC[1,1], AIC[3,1]-AIC[2,1], AIC[4,1]-AIC[3,1], AIC[5,1]-AIC[4,1]),2)
#write.csv(AIC, 'Inher_po_AIC_Mar19.csv', row.names=F)
#write.csv(summary(model4)$coef, "model_output_inher_po_MAX_Mar19.csv", row.names=T)
#save(model4,file="MAX_T_inher_PO.RData")

#### HALF FALL ####
# interpolated using everything 
# extrapolate using fewer (i.e. no beech)
HF_output_K <- mapply(run_half_fall, test_years = list_K_fold,
                      MoreArgs = list(datafile), SIMPLIFY = F)
save(HF_output_K, file = "HF_output_K.RData")

#### Import parameter values ####
# SURVIVAL
# import
# res only
load('model_outputs_full_K1_11.10.18.RData')
load('model_outputs_full_K2_11.10.18.RData')
survival_RO <- c(model_outputs_full_K1, model_outputs_full_K2)

# tot pop
load('model_outputs_tot_pop_K1.RData')
load('model_outputs_tot_pop_K2.RData')
survival_TP <- c(model_outputs_tot_pop_K1, model_outputs_tot_pop_K2)

load('final_survival_Feb19.RData')


# Extract parameter values
# order = 
# int, slope.z, slope.sync, slope.sync2, slope.spring p, slope.spring t, slope.winter t, slope.pop s, 
# slope.b
s.params_RO <- mapply(extract_surv, model = survival_RO, 
                   MoreArgs = list(type="res_only"), SIMPLIFY = F)
s.params_TP <- mapply(extract_surv, model = survival_TP, 
                        MoreArgs = list(type="tot_pop"), SIMPLIFY = F)

s.params_final <- extract_surv(model=final_survival, type ='tot_pop')


# DEVELOPMENT
# import
# res only
load('Dev_final_K_11.10.RData')

# tot pop
load('Dev_final_K_TP_Feb19.RData')

#load('final_development_Feb19.RData')
load('final_development_MAX.RData')

# Extract parameter values
d.params_RO <- mapply(extract_dev, model = Dev_final_K, SIMPLIFY = F)
d.params_TP <- mapply(extract_dev, model = Dev_final_K_tot_pop, SIMPLIFY = F)

d.params_final <- extract_dev(model=final_development)

# RECRUITMENT
# import
# res only
load('Recru_final_K_Feb19.RData')

# tot pop
load('Recru_final_K_TP_Feb19.RData')

load('final_recruitment_Feb19.RData')

# Extract parameter values
r.params_RO <- mapply(extract_recru, model = Recru_final_K, SIMPLIFY = F)

r.params_TP <- mapply(extract_recru, model = Recru_final_K_tot_pop, SIMPLIFY = F)

r.params_final <- extract_recru(final_recruitment)

# INHERITANCE
# import
load('final_inheritance_PO_Feb19.RData') # keep it to resident only as only 14 immigrants


# loading all exhausts memory
load('final_inheritance_QG_Feb19.RData') # keep to resident only as only 20 immigrants with mother

load('Inher_QG_K_tot_pop_Feb19.RData')
load('Inher_QG_full_19.10.18.RData')

load('Inher_PO_tot_pop.RData')

# Extract parameter values

ipo.params_TP <- mapply(extract_inher, model = Inher_PO_tot_pop, 
                        MoreArgs = list(method = "PO"), SIMPLIFY = F)

iqg.params_RO <- mapply(extract_inher, model = Inher_QG_full, 
                          MoreArgs = list(method = "QG"), SIMPLIFY = F)
iqg.params_TP <- mapply(extract_inher, model = Inher_QG_full_tot_pop, 
                          MoreArgs = list(method = "QG"), SIMPLIFY = F)

iqg.params_final <- extract_inher(final_inheritance_QG, method='QG')
#iqg.params_final[8] <- 0.0001 # test with VA = 0
ipo.params_final <- extract_inher(final_inheritance_PO, method='PO')

# HALF FALL
load("HF_output_K.RData")

# extract

HF.params <- mapply(extract_HF, model = HF_output_K, SIMPLIFY = F)

#### Functions set up in Function script ####

#### Run models - CROSS VALIDATION ####

m.N <- mean(descale$R_pop_size) # mean of pop size so I can scale and unscale
sd.N <- sd(descale$R_pop_size) # sd of population size

cross_val_RESULT_RO <- mapply(run_CV_IPM, test_years = list_K_fold, s.params = s.params_RO, 
                              d.params = d.params_RO, r.params = r.params_RO, h.params = iqg.params_RO,
                              hf.params = HF.params, MoreArgs = list(n = 100, n.gens = 6, g.range = c(-2,2), 
                                                                     m.N = m.N, sd.N = sd.N, CV_data = cross_val,
                                                                     type = "res_only", method = "CV",
                                                                     descale = descale), SIMPLIFY = F)
save(cross_val_RESULT_RO, file = "cross_val_RESULT_RO_Nov19.RData")

m.N_TP <- mean(descale$pop_size) # mean of pop size so I can scale and unscale
sd.N_TP <- sd(descale$pop_size)
# set 1960 cross val to total pop size
cross_val$R_only_PS[1] <- cross_val$Pop_size[1]   

source('Functions_IPM1.R')
cross_val_RESULT_TP <- mapply(run_CV_IPM, test_years = list_K_fold, s.params = s.params_TP, 
                                d.params = d.params_TP, r.params = r.params_TP, h.params = iqg.params_TP,
                                hf.params = HF.params, MoreArgs = list(n = 100, n.gens = 6, g.range = c(-2,2), 
                                                                       m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                                                                       type = "tot_pop", method = "CV",
                                                                       descale = descale), SIMPLIFY = F)
#save(cross_val_RESULT_TP, file = "cross_val_RESULT_TP_Nov19_N_double.RData")
save(cross_val_RESULT_TP, file = "cross_val_RESULT_TP_Nov19.RData")

cross_val_RESULT_po <- mapply(run_CV_IPM_PO, test_years = list_K_fold, s.params = s.params_TP, 
                              d.params = d.params_TP, r.params = r.params_TP, h.params = ipo.params_TP,
                              hf.params = HF.params, MoreArgs = list(n = 10000, n.gens = 6, g.range = c(-4,4), 
                                                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                                                     CV_data = cross_val,
                                                                     descale = descale,
                                                                     method="CV"), SIMPLIFY = F)
#save(cross_val_RESULT_po, file = "cross_val_RESULT_po_Nov19_N_double.RData")
save(cross_val_RESULT_po, file = "cross_val_RESULT_po_Nov19.RData")

cross_val_RESULT_TP <- run_CV_IPM(test_years = 1961:2009, s.params = s.params_final, 
                              d.params = d.params_final, r.params = r.params_final, 
                              h.params = iqg.params_final,
                              hf.params = half_fall_params, n = 100, n.gens = 50, 
                              g.range = c(-2,2), 
                              m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                              type = "tot_pop", method = "CV",
                              descale = descale, decoupled=FALSE)
save(cross_val_RESULT_TP, file = "cross_val_RESULT_TP_50_Nov19.RData")

cross_val_RESULT_po <- run_CV_IPM_PO(test_years = 1961:2009, s.params = s.params_final, 
                              d.params = d.params_final, r.params = r.params_final, h.params = ipo.params_final,
                              hf.params = half_fall_params, n = 10000, n.gens = 50, g.range = c(-4,4), 
                                                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                                                     CV_data = cross_val,
                                                                     descale = descale, 
                              method = "CV")
save(cross_val_RESULT_po, file = "cross_val_RESULT_po_50_Nov19.RData")

### DIRECTIONAL CHANGE

# set up a vector of spring temperature
summary(lm(Spring_temp ~ Year, data = descale))
means <- seq(mean(descale$Spring_temp)+sd(descale$Spring_temp), 
             mean(descale$Spring_temp)+sd(descale$Spring_temp)+(0.04*50),
             length.out=51)
set.seed(10)
cross_val$s_temp <- rnorm(51, mean=means, sd=sd(descale$Spring_temp))

cross_val_RESULT_TP <- run_CV_IPM(test_years = 1961:2010, s.params = s.params_final, 
                                  d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                  hf.params = half_fall_params, n = 100, n.gens = 50, g.range = c(-2,2), 
                                  m.N = m.N_TP, sd.N = sd.N_TP, CV_data = cross_val,
                                  type = "tot_pop", method = "Predict",
                                  descale = descale)
save(cross_val_RESULT_TP, file = "cross_val_RESULT_TP_50_DIRECTIONAL_Nov.RData")

cross_val_RESULT_po <- run_CV_IPM_PO(test_years = 1961:2010, s.params = s.params_final, 
                                     d.params = d.params_final, r.params = r.params_final, h.params = ipo.params_final,
                                     hf.params = half_fall_params, n = 10000, n.gens = 50, g.range = c(-4,4), 
                                     m.N = m.N_TP, sd.N = sd.N_TP, 
                                     CV_data = cross_val,
                                     descale = descale, method = "Predict")
save(cross_val_RESULT_po, file = "cross_val_RESULT_po_50_DIRECTIONAL_Nov.RData")



#### PERTURBATION ANALYSIS ####
# here will perturb one environmental variable at a time
# run model for 50 generations and then calculate pop size
# will perturb from -0.5 to +0.5 (all centered on 0)
# also summarise how this perturb relates to true values

# first create perturbation vector
changer <- seq(-4,4,length.out=51) # 0.05 change < 1 standard deviation. So do 0.2
# change of 2 standard deviations
context <- data.frame(Variable = c("ST", "SP", "WT", "WP", "HF"),
                      Range_U = rep(NA, 5),
                      Range_L = rep(NA, 5),
                      Range_UC = rep(NA, 5),
                      Range_UL = rep(NA, 5))
context[1,2:3] <- range(descale$Spring_temp)
context[1,4:5] <- range((changer*sd(descale$Spring_temp)) + mean(descale$Spring_temp))
context[2,2:3] <- range(descale$spring_precip_t)
context[2,4:5] <- range(((changer)*sd(descale$spring_precip_t)) + mean(descale$spring_precip_t))
context[3,2:3] <- range(descale$winter_temp)
context[3,4:5] <- range((changer*sd(descale$winter_temp)) + mean(descale$winter_temp))
context[4,2:3] <- range(descale$winter_precip_t)
context[4,4:5] <- range((changer*sd(descale$winter_precip_t)) + mean(descale$winter_precip_t))
context[5,2:3] <- range(descale$Half_fall)
context[5,4:5] <- range((changer*sd(descale$Half_fall))+ mean(descale$Half_fall))

# make data order = ST, WP, WT, SP, HF
# HF must be NA when not focal change
perturbation <-data.frame(ST = c(changer, rep(mean(datafile$Spring_temp_max), 51*5)),
                          WP = c(rep(mean(datafile$winter_precip_t), 51), changer, 
                                 rep(mean(datafile$winter_precip_t), 51*4)),
                          WT = c(rep(mean(datafile$winter_temp), 51*2), changer, 
                                 rep(mean(datafile$winter_temp), 51*3)),
                          SP = c(rep(mean(datafile$spring_precip_t), 51*3), changer, 
                                 rep(mean(datafile$spring_precip_t), 51*2)),
                          #HF = c(rep(NA, 51*4), -0.0178+changer))
                          #HF = c(rep(NA, 51*4), -0.113+changer)) # for max
                          HF = c(rep(NA, 51*4), changer, rep(NA, 51)),
                          ST.cat = c(rep(mean(datafile$Spring_temp_cat), 51*5), changer))
                          
#perturbation <-data.frame(ST = (as.numeric(mean(low_predictions[[1]]$s_temp))-mean(descale$Spring_temp))/sd(descale$Spring_temp),
 #                         WP = (as.numeric(mean(low_predictions[[1]]$w_precip))-mean(descale$winter_precip_t))/sd(descale$winter_precip_t),
  #                        WT = (as.numeric(mean(low_predictions[[1]]$w_temp))-mean(descale$winter_temp))/sd(descale$winter_temp),
   #                       SP = (as.numeric(mean(low_predictions[[1]]$s_precip))-mean(descale$spring_precip_t))/sd(descale$spring_precip_t),
    #                      HF = NA)

# run the model
m.N <- mean(descale$pop_size) # mean of pop size so I can scale and unscale
sd.N <- sd(descale$pop_size) # sd of population size
descale$Pop_size <- descale$pop_size

perturb_results <- apply(perturbation, 1, FUN = run_CV_IPM, test_years = NULL, s.params = s.params_final,
      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
      hf.params = half_fall_params_max, n = 100, n.gens = 50, g.range = c(-2,2), 
      method = "Perturb", m.N = m.N, sd.N = sd.N, type = "tot_pop", descale = descale, decoupled = TRUE)

write.csv(perturb_results, "perturb_results_ExtremeHF_extra_SEPT.csv", row.names=F)

#### IPM2 PREDICTION - COUPLED ####
# make sure to load all parameter values first
# then load predictions

load("low_predictions_FINAL.RData")
load("mid_predictions_Oct18.RData")
load("high_predictions_Oct18.RData")

# try running for low predictions
m.N <- mean(descale$pop_size) # mean of pop size so I can scale and unscale
sd.N <- sd(descale$pop_size) # sd of population size
low_results <- mapply(run_CV_IPM, CV_data = low_predictions[[1]], 
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-3,3.6), 
                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                      hf.params = half_fall_params), SIMPLIFY = F)

save(low_results, file="low_results_VA0.RData")

mid_results <- mapply(run_CV_IPM, CV_data = mid_predictions2, 
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-3,3.6), 
                                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                      hf.params = half_fall_params), SIMPLIFY = F)

save(mid_results, file="mid_results_VA0.RData")

high_results <- mapply(run_CV_IPM, CV_data = high_predictions2, 
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-3,3.6), 
                                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                      hf.params = half_fall_params), SIMPLIFY = F)

save(high_results, file="high_results_VA0.RData")

#### IPM2 PREDICTION - DE-COUPLED ####
# make sure to load all parameter values first
# then load predictions

datafile <- read.csv("updated_model_input_MAX.csv", header=T) 
descale <- read.csv("updated_descale_MAX.csv", header=T)
datafile_t <- read.csv("datafile_t_MAX.csv", header=T)
datafile_t1 <- read.csv("datafile_t1_MAX.csv", header=T)
datafile_t <- cbind(datafile_t, datafile_t1)
datafile_t$Year <- datafile_t$Year + 1
datafile_I <- read.csv("datafile_I_MAX.csv", header=T)
# need to have an animal and ID and a BYEAR column.
datafile_I$Age[which(datafile_I$Age > 1)] <- 2
datafile_I$Age <- as.factor(datafile_I$Age)
datafile_I$beech <- as.factor(datafile_I$beech)
datafile_I$Immigrant <- as.factor(datafile_I$Immigrant)
datafile_I$animal <- as.factor(datafile_I$animal)
ped <- read.csv("ped_inheritance_11.10.16.csv", header=T)
cross_val <- read.csv("cross_validation_data_1961to2010_19.10.csv", header=T)
half_fall_params_max <- as.vector(coef(lm(Half_fall ~ Spring_temp_cat + spring_precip_t + winter_temp + 
                                            winter_precip_t + as.factor(beech), data=datafile)))

load('final_development_MAX.RData')

# Extract parameter values
d.params_final <- extract_dev(model=final_development)

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")


# try running for low predictions
m.N <- mean(descale$pop_size) # mean of pop size so I can scale and unscale
sd.N <- sd(descale$pop_size) # sd of population size
low_results <- mapply(run_CV_IPM, CV_data = low_predictions,
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                      hf.params = half_fall_params_max,
                                      decoupled = TRUE), SIMPLIFY = F)

save(low_results, file="low_results_decoupled_VA0.RData")
#save(low_results, file="low_results_decoupled_SEPT.RData")

mid_results <- mapply(run_CV_IPM, CV_data = mid_predictions,
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                      hf.params = half_fall_params_max,
                                      decoupled = TRUE), SIMPLIFY = F)

save(mid_results, file="mid_results_decoupled_VA0.RData")
#save(mid_results, file="mid_results_decoupled_SEPT.RData")

high_results <- mapply(run_CV_IPM, CV_data = high_predictions,
                       MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                       m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                       descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                       d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                       hf.params = half_fall_params_max,
                                       decoupled = TRUE ), SIMPLIFY = F)

save(high_results, file="high_results_decoupled_VA0.RData")
#save(high_results, file="high_results_decoupled_SEPT.RData")

#### IPM2 PREDICTION - NO CHANGE ####
# make sure to load all parameter values first
# then load predictions

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")

# try running for low predictions
m.N <- mean(descale$pop_size) # mean of pop size so I can scale and unscale
sd.N <- sd(descale$pop_size) # sd of population size
no_ST_results <- mapply(run_CV_IPM, CV_data = high_predictions,
                      MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                      m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                      descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                      d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                      hf.params = half_fall_params_max,
                                      decoupled = TRUE, weather_var = "ST"), SIMPLIFY = F)
save(no_ST_results, file="only_ST_results_SEPT.RData")

no_SP_results <- mapply(run_CV_IPM, CV_data = high_predictions,
                        MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                        m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                        descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                        d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                        hf.params = half_fall_params_max,
                                        decoupled = TRUE, weather_var = "SP"), SIMPLIFY = F)
save(no_SP_results, file="only_SP_results_SEPT.RData")

no_WT_results <- mapply(run_CV_IPM, CV_data = high_predictions,
                        MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                        m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                        descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                        d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                        hf.params = half_fall_params_max,
                                        decoupled = TRUE, weather_var = "WT"), SIMPLIFY = F)
save(no_WT_results, file="only_WT_results_SEPT.RData")

no_WP_results <- mapply(run_CV_IPM, CV_data = high_predictions,
                        MoreArgs = list(n = 100, n.gens = 91, g.range = c(-2,2), 
                                        m.N = m.N, sd.N = sd.N, type = "tot_pop", method = "Predict",
                                        descale = descale, test_years = 2010:2100, s.params = s.params_final, 
                                        d.params = d.params_final, r.params = r.params_final, h.params = iqg.params_final,
                                        hf.params = half_fall_params_max,
                                        decoupled = TRUE, weather_var = "WP"), SIMPLIFY = F)
save(no_WP_results, file="only_WP_results_SEPT.RData")



