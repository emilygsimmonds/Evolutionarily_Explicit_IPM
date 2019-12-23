#### Script containing wrapper functions to run each 
# demographic function and the caterpillar timing function 
# for all of the cross validation test/training dataset combinations

# COLLAPSE ALL FOR EASIER VIEWING

#### run surv model ####

#Inputs: processed data, design data, parameters for Phi and p

run_survival <- function(data_processed, ddl, pars){
  model <- crm(data_processed, ddl, 
               hessian = T, model="CJS", 
               model.parameters=list(Phi=pars[[1]], p=pars[[2]]), 
               accumulate = F) # don't want to accumulate histories
  return(model)
}

#### run recru model ####

#Inputs: scaled data, which population size variable to use (total or resident only)
# which years are test years

run_recruitment <- function(datafile, type = c("tot_pop", "res_only"), test_years = NULL){
  if(is.null(test_years)==TRUE){temp <- datafile}else{marker <- c(which(datafile$Year == test_years[1]), which(datafile$Year == test_years[2]),
                                                                  which(datafile$Year == test_years[3]), which(datafile$Year == test_years[4]),
                                                                  which(datafile$Year == test_years[5]))
  temp <- datafile[-marker,]}
  if(type == "tot_pop"){   
    model <- glmer(Num_recruited ~ April_hatch + Synchrony + I(Synchrony^2) + Clutch_size + I(Clutch_size^2) +
                     as.factor(Round) + pop_size + winter_precip_t + spring_precip_t + as.factor(beech) +
                     + (1|Year) + (1|Nest_box),
                   data=temp, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))}
  if(type == "res_only"){  
    model <- glmer(Num_recruited ~ April_hatch + Synchrony + I(Synchrony^2) + Clutch_size + I(Clutch_size^2) +
                     as.factor(Round) + R_pop_size + winter_precip_t + spring_precip_t + as.factor(beech) +
                     + (1|Year) + (1|Nest_box),
                   data=temp, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))}
  return(model)
}

#### run development model ####

#Inputs: scaled data for t and t+1, which years are test years, 
# which population size variable to use (total or resident only)

run_development <- function(datafile_t, datafile_t1, test_years = NULL, type = c("tot_pop", "res_only")){
  datafile <- cbind(datafile_t, datafile_t1)
  if(is.null(test_years)==TRUE){datafile <- datafile}else{marker <- c(which(datafile$Year == test_years[1]), which(datafile$Year == test_years[2]),
                                                                      which(datafile$Year == test_years[3]), which(datafile$Year == test_years[4]),
                                                                      which(datafile$Year == test_years[5]))
  datafile$Year <- datafile$Year + 1
  datafile <- datafile[-marker,]}
  if(type == "tot_pop"){model <- lmer(April_hatch_t1 ~ April_hatch + Synchrony + I(Synchrony^2) 
                                      + Spring_temp_t1 + as.factor(Age) + pop_size_t1 + 
                                        as.factor(Round_t1) + (1|Year), data=datafile)}
  if(type == "res_only"){model <- lmer(April_hatch_t1 ~ April_hatch + Synchrony + I(Synchrony^2) 
                                       + Spring_temp_t1 + as.factor(Age) + R_pop_size_t1 + 
                                         as.factor(Round_t1) + (1|Year), data=datafile)}
  return(model)
}

#### run inheritance model ####

#Inputs: scaled data for inheritance, which years are test years, 
# which population size variable to use (total or resident only),
# random priors

### quant genetics
run_inheritance_QG <- function(datafile_I, ped, test_years = NULL, type = c("tot_pop", "res_only"),
                               prior.rand){
  if(is.null(test_years) == TRUE){datafile_I <- datafile_I}else{marker <- c(which(datafile_I$Year == test_years[1]), which(datafile_I$Year == test_years[2]),
                                                                            which(datafile_I$Year == test_years[3]), which(datafile_I$Year == test_years[4]),
                                                                            which(datafile_I$Year == test_years[5]))
  datafile_I <- datafile_I[-marker,]}
  if(type == "tot_pop"){model <- MCMCglmm(April_hatch ~ Spring_temp + winter_temp + winter_precip_t + 
                                            spring_precip_t + 
                                            Round + pop_size, random=~animal+ID+MOTHER, pedigree=ped, 
                                          data=datafile_I, prior=prior.rand, verbose=F, 
                                          nitt=50000, burnin=5000, thin=50)
  }
  if(type == "res_only"){model <- MCMCglmm(April_hatch ~ Spring_temp + winter_temp + winter_precip_t + 
                                             spring_precip_t + 
                                             Round + R_pop_size, random=~animal+ID+MOTHER, pedigree=ped, 
                                           data=datafile_I, prior=prior.rand, verbose=F, nitt=50000, burnin=5000, thin=50)
  }
  return(model)
}

### standard IPM

# PO = parent-offspring as it uses a parent-offspring regression

#Inputs: scaled data for inheritance, which years are test years, 
# which population size variable to use (total or resident only),

run_inheritance_PO <- function(datafile_I, test_years = NULL, type = c("tot_pop", "res_only")){
  if(is.null(test_years)==TRUE){datafile_I <- datafile_I}else{marker <- c(which(datafile_I$Year == test_years[1]), which(datafile_I$Year == test_years[2]),
                                                                          which(datafile_I$Year == test_years[3]), which(datafile_I$Year == test_years[4]),
                                                                          which(datafile_I$Year == test_years[5]))
  datafile_I <- datafile_I[-marker,]}
  if(type == "tot_pop"){model <- lmer(April_hatch ~ Mother_hatch + pop_size + Spring_temp +
                                        spring_precip_t + Round +
                                        (1|Year), data=datafile_I)}
  if(type == "res_only"){model <- lmer(April_hatch ~ Mother_hatch + R_pop_size + Spring_temp +
                                         spring_precip_t + Round +
                                         (1|Year), data=datafile_I)}
  return(model)
}

#### Run caterpillar timing model ####

# Measure of caterpillar timing is the half-fall of fifth instar larvae

# Inputs: scaled data, indication of the test years

run_half_fall <- function(datafile, test_years){
  marker <- c(which(datafile$Year == test_years[1]), which(datafile$Year == test_years[2]),
              which(datafile$Year == test_years[3]), which(datafile$Year == test_years[4]),
              which(datafile$Year == test_years[5]))
  if(length(marker > 0)){datafile <- datafile[-marker,]}
  model <- lm(Half_fall ~ Spring_temp + spring_precip_t + winter_temp + 
                winter_precip_t + as.factor(beech), data=datafile)
  return(model)
}

