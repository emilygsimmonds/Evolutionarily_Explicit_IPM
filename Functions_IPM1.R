# Script with functions
# to run statistical analyses to capture fundamental processes

# Beginning: 14.6.18
# Author: EGS
#-----------------------------------------------------------------------------

#### create marked data ####

create_marked <- function(datafile, basta_input, full_years, test_years, type=c("K_fold","Full")){
  # split capture history by removing test years
  if(type=="K_fold"){marker <- seq(which(colnames(basta_input) == test_years[1]),
                which(colnames(basta_input) == test_years[5]),1)
  temp <- basta_input[,-marker]}else{temp <- basta_input}
  # remove individuals will all 0 occurances
  marker2 <- which(apply(temp[,2:length(temp)], 1, sum) == 0)
  if(length(marker2) > 0){temp <- temp[-marker2,]}
  
  temp$capture_history <- apply(temp[,2:length(temp[1,])] , 1 , paste , collapse = "" )
  end <- rep("1;", length(temp[,1]))
  
  # add non time varying covariates to file
  IDs <- row.names(temp)
  temp$Immigrant <- rep(NA, length(temp[,1]))
  temp$Round <- rep(NA, length(temp[,1]))
  
  # for each individual add their Immigrant status, birth year hatch date, and section of woodland
  for(l in IDs){
    temp_data <- subset(datafile, datafile$F_ID == l)
    marker <- which(row.names(temp)==l)
    ifelse(temp_data$Immigrant[1] == 1, temp$Immigrant[marker] <- "Immigrant", temp$Immigrant[marker] <- "Resident")
    temp$Round[marker] <- as.character(temp_data$Round[which.max(temp_data$Round)])
  }

  # now we want the individual covariates
  if(type=="K_fold"){marker <- seq(which(full_years == min(test_years)-1), which(full_years == max(test_years)-1), 1)
  years <- full_years[-marker]}else{years <- full_years}
  output_store <- matrix(NA, ncol=length(years)*5, nrow=length(temp[,1]))
  
  # this loop should create output that has one row for each individual and columns for 
  # every value for each variable for 50 years
  
  for(l in 1:length(years)){
    
    year_data <- subset(datafile, datafile$Year == years[l]+1) # first subset to the year you want
    
    for(m in row.names(temp)){ # now for each ID
      
      marker <- which(as.character(year_data$F_ID) == m) # find when the first ID appears in year data
      marker2 <- which(row.names(temp) == m) # generates marker in each dataset that references the correct ID
      
      # Now want for each variable to fill in the first year for each individual so var 1 = 1 var 2 = 37 etc
      
      ifelse(length(marker) > 0, output_store[marker2,l] <- year_data[marker,6], output_store[marker2,l] <- mean(year_data[,6], na.rm=T))
      ifelse(length(marker) > 0, output_store[marker2,l+length(years)] <- year_data[marker,7], output_store[marker2,l+(length(years))] <- mean(year_data[,7], na.rm=T))
      ifelse(length(marker) > 0, output_store[marker2,l+(length(years)*2)] <- year_data[marker,8], output_store[marker2,l+(length(years)*2)] <- mean(year_data[,8], na.rm=T))
      ifelse(length(marker) > 0, output_store[marker2,l+(length(years)*3)] <- (year_data[marker,6]^2), output_store[marker2,l+(length(years)*3)] <- (mean(year_data[,6], na.rm=T)^2))
      ifelse(length(marker) > 0, output_store[marker2,l+(length(years)*4)] <- (year_data[marker,7]^2), output_store[marker2,l+(length(years)*4)] <- (mean(year_data[,7], na.rm=T)^2))                   
    }
    
  }
  
  col_name_data <- data.frame(April_hatch = "ah", Synchrony = "sn", Clutch_size = "cs",
                              April_hatch2 = "ahs", Synchrony2 = "sns", Year = seq(1960, 2004, 1)) 
  if(type=="Full"){col_name_data <- data.frame(April_hatch = "ah", Synchrony = "sn", Clutch_size = "cs",
                                               April_hatch2 = "ahs", Synchrony2 = "sns", Year = seq(1960, 2009, 1)) }
  # need to name them in specific way
  marker1 <- c(1,6)
  marker2 <- c(2,6)
  marker3 <- c(3,6)
  marker4 <- c(4,6)
  marker5 <- c(5,6)
  

  mark_input <- cbind(as.character(temp$capture_history), 
                      temp$Immigrant, 
                      temp$Round, output_store, end)
  column_names <- c(apply(col_name_data[,marker1] , 1 , paste , collapse = "" ), apply(col_name_data[,marker2] , 1 , paste , collapse = "" ), 
                    apply(col_name_data[,marker3] , 1 , paste , collapse = "" ), apply(col_name_data[,marker4] , 1 , paste , collapse = "" ), 
                    apply(col_name_data[,marker5] , 1 , paste , collapse = "" ))
  colnames(mark_input) <- c("ch", "Immigrant", "Round", column_names, "end")
  mark_input <- as.data.frame(mark_input)
  mark_input$ch <- as.character(mark_input$ch)
  return(mark_input)
}

#### create design data with year covariates ####

create_design <- function(processed_data, test_years, full_years, year_variables, parameters){
  ddl <- make.design.data(processed_data, parameters = parameters)
  marker <- seq(which(full_years == min(test_years)-1), which(full_years == max(test_years)-1), 1)
  years <- full_years[-marker]
  # Now have dataframe organised by year and individual. Here is where we want to add our year varying variables.
  year_variables_temp1 <- year_variables[-marker,]
  year_variables_temp1$time[1:45] <- 1960:2004
  ddl$Phi <- merge_design.covariates(ddl$Phi, year_variables_temp1, bytime=TRUE)
  ddl$Phi$Beech <- as.factor(ddl$Phi$Beech)
  year_variables_temp2 <- year_variables_temp1
  year_variables_temp2$time[1:45] <- 1961:2005
  ddl$p <- merge_design.covariates(ddl$p, year_variables_temp2, bytime=TRUE)
  return(ddl)
}

#### logistic transform ####
logistic <- function(u){
  x <- exp(u)
  return(x/(1+x))
}

#### run surv model ####

run_survival <- function(data_processed, ddl, pars){
  model <- crm(data_processed, ddl, 
                  hessian = T, model="CJS", 
                  model.parameters=list(Phi=pars[[1]], p=pars[[2]]), 
                  accumulate = F) # don't want to accumulate histories
  return(model)
}

#### run recru model ####

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
# quant genetics
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

# parent-offspring

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

#### Run half fall model ####

run_half_fall <- function(datafile, test_years){
  marker <- c(which(datafile$Year == test_years[1]), which(datafile$Year == test_years[2]),
              which(datafile$Year == test_years[3]), which(datafile$Year == test_years[4]),
              which(datafile$Year == test_years[5]))
  if(length(marker > 0)){datafile <- datafile[-marker,]}
  model <- lm(Half_fall ~ Spring_temp + spring_precip_t + winter_temp + 
                winter_precip_t + as.factor(beech), data=datafile)
  return(model)
}

#### EXTRACT PARAM VALUES ####

extract_surv <- function(model, type = c("tot_pop", "res_only")){
  if(type == "tot_pop"){s.params <- c(model$results$beta$Phi[1], # intercept
                                      model$results$beta$Phi[11], # hatch date 
                                      model$results$beta$Phi[14], # pop slope
                                      model$results$beta$Phi[13], # synchrony 
                                      model$results$beta$Phi[12], # synchrony^2 
                                      model$results$beta$Phi[17], # spring p 
                                      model$results$beta$Phi[19], # winter t
                                      model$results$beta$Phi[18], # spring t
                                      mean(model$results$beta$Phi[2:9]), # round
                                      model$results$beta$Phi[10], # resident
                                      model$results$beta$Phi[15:16] # beech
                                   )}
  if(type == "res_only"){s.params <- c(model$results$beta$Phi[1], # intercept
                                   model$results$beta$Phi[11], # hatch date 
                                   model$results$beta$Phi[14], # pop slope
                                   model$results$beta$Phi[13], # synchrony 
                                   model$results$beta$Phi[12], # synchrony^2 
                                   model$results$beta$Phi[17], # spring p 
                                   model$results$beta$Phi[18], # winter t
                                   mean(model$results$beta$Phi[2:9]), # round
                                   model$results$beta$Phi[10], # resident
                                   model$results$beta$Phi[15:16] # beech
                                   )}
  return(as.numeric(s.params))
}

extract_dev <- function(model){
               d.params <- c(summary(model)$coefficients[1], # intercept
                                   summary(model)$coefficients[2], # hatch date
                                   summary(model)$coefficients[7], # pop size
                                   summary(model)$coefficients[3], # synchrony 
                                   summary(model)$coefficients[4], # synchrony^2 
                                   summary(model)$coefficients[5], # spring t
                                   summary(model)$coefficients[6], # age
                                   mean(summary(model)$coefficients[8:15]), # round
                                   sd(resid(model)) # sigma
                                   )
  return(d.params)
}

extract_recru <- function(model){
                 r.params <- c(summary(model)$coefficients[1], # intercept
                                   summary(model)$coefficients[2], # hatch date 
                                   summary(model)$coefficients[15], # pop size
                                   summary(model)$coefficients[3], # synchrony 
                                   summary(model)$coefficients[4], # synchrony^2
                                   summary(model)$coefficients[16], # winter precip
                                   summary(model)$coefficients[17], # spring precip
                                   summary(model)$coefficients[5], # clutch size
                                   summary(model)$coefficients[6], # clutch size^2
                                   mean(summary(model)$coefficients[7:14]), # round
                                   summary(model)$coefficients[18:19] # beech
                                   )
  return(r.params)
}

extract_inher <- function(model, method = c("PO", "QG")){
  if(method == "PO"){i.params <- c(summary(model)$coefficients[1], # intercept
                                   summary(model)$coefficients[3], # pop size
                                   summary(model)$coefficients[2], # mother hatch
                                   summary(model)$coefficients[4], # spring t 
                                   summary(model)$coefficients[5], # spring precip
                                   mean(summary(model)$coefficients[6:13]), # round
                                   sd(resid(model)) # sigma
                                   )
  }
  if(method == "QG"){results <- posterior.mode(mcmc(cbind(model$Sol, model$VCV)))
                                  i.params <- c(results[1], # intercept
                                                results[14], # pop slope
                                                results[2], # Spring temp
                                                results[5], # Spring p
                                                results[4], # Winter p
                                                results[3], # winter t
                                                mean(results[6:13]), # round
                                                results[15], # animal
                                                results[18] # units - sigma 
                                                )
  }
  return(as.numeric(i.params))
}

extract_HF <- function(model){
  HF.params <- c(summary(model)$coefficients[1], # int
                 summary(model)$coefficients[2], # Spring temp
                 summary(model)$coefficients[4], # Winter temp
                 summary(model)$coefficients[3], # Spring precip
                 summary(model)$coefficients[5], # Winter precip
                 summary(model)$coefficients[6], # beech
                 summary(model)$coefficients[7]) 
 return(HF.params)
}

#### FUNDAMENTAL FUNCTIONS ####

# Survival - minimal form so most versatile
# covariates all included as 'extra'
S.fun <- function(z, int, slope.z, slope.PS, extra, Ns){ # extra = all other slope and sum
  u <- exp(int+(z*slope.z)+(Ns*slope.PS)+extra)
  # need to add on resident value at end
  u/(1+u)
}

# Development
D.fun <- function(z, zz, mu.int, mu.zs, mu.ps, extra, sigma.d, Ns) {
  mu.z <- mu.int+(mu.zs*z)+(mu.ps*Ns)+extra
  sigma.d2 <- sigma.d^2
  temp1 <- sqrt(2*pi)*sigma.d
  temp2 <- ((zz-mu.z)^2)/(2*sigma.d2)
  return(exp(-temp2)/temp1)
}

# Recruitment
R.fun <- function(z, int, slope.z, slope.ps, extra, Ns){
  u <- exp(int+(z*slope.z)+(Ns*slope.ps)+extra)
  u
}

# Inheritance
H.fun <- function(z, zz, mu.int, mu.N, sigma.h, extra, Ns) {
  mu.z <- mu.int+(mu.N*Ns)+extra
  sigma.h2 <- sqrt(sigma.h)
  temp1 <- sqrt(2*pi)*sigma.h2
  temp2 <- ((zz-mu.z)^2)/(2*sigma.h)
  return(exp(-temp2)/temp1)
}
H.fun.PO <- function(z, zz, mu.int, mu.N, mu.mZ, sigma.h, extra, Ns) {
  mu.z <- mu.int+(mu.N*Ns)+extra
  sigma.h2 <- sigma.h^2
  temp1 <- sqrt(2*pi)*sigma.h
  temp2 <- ((zz-mu.z)^2)/(2*sigma.h2)
  return(exp(-temp2)/temp1)
}


# Half fall 
HF.fun <- function(int, slope.ST, slope.WT, slope.SP, slope.WP, ST, WT, WP, SP, B) {
  # need to subtract the difference in intercept for beech not beech value
  u <- int+(ST*slope.ST)+(WT*slope.WT)+(WP*slope.WP)+(SP*slope.SP)+B
  u
}

# calculates moments of the results
moments.fun <- function(x,z){
  mu <- sum(z*x)/sum(x) # calculates mean
  va <- sum(z*z*x)/sum(x)-mu^2 # calculates variance
  sk <- (sum(z*z*z*x)/sum(x)-3*mu*va-mu^3)/(sqrt(va)^3) # calculates skew
  m4 <- sum(z*z*z*z*x)/sum(x)
  m3 <- sum(z*z*z*x)/sum(x)
  m2 <- sum(z*z*x)/sum(x)
  ku <- (m4-4*mu*m3+6*(mu^2)*m2-3*(mu^4))/(va^2)-3
  res <- c(mu,va,sk,ku,m2,m3,m4)
  names(res) <- c('Mean','Variance','Skew','Kurtosis','Moment 2','Moment 3','Moment 4')
  return(res)
}

#### RUN IPM ####

run_CV_IPM <- function(n, n.gens, g.range,
                       m.N, sd.N, CV_data, type = c("tot_pop", "res_only"),
                       method = c("CV", "Perturb", "Predict"),
                       test_years, descale, s.params, d.params, r.params, h.params,
                       hf.params,
                       decoupled = FALSE, 
                       weather_var = "NO",
                       IPM1 = FALSE){
  # num classes = n
  # num iterations = n.gens
  # min and max G and E values = g.range
  # mean and sd pop size for scaling = m.N, sd.N
  # cross validation data = CV_data
  # type = type of model
  # years we are predicting = test_years
  # data for scaling = descale
  # parameters for functions = s,d,r,h.params 

  # results vectors 
  res.A <- array(NA,c(n*n,n.gens)) # adults
  res.J <- array(NA,c(n*n,n.gens)) # juveniles
  res <- array(NA,c(n*n,n.gens)) # overall results
  res.fit.A <- array(NA,c(n*n,n.gens)) # storage vector to save fitness at each stage
  res.fit.J <- array(NA,c(n*n,n.gens)) # storage vector to save fitness at each stage
  res.d <- array(NA,c(n*n,n.gens)) # D%S
  res.d.A <- array(NA,c(n*n,n.gens)) # D%S
  res.d.J <- array(NA,c(n*n,n.gens)) # D%S
  res.r <- array(NA,c(n*n,n.gens)) # storage vector for recruitment
  res.h <- array(NA,c(n*n,n.gens)) # R%H
  s.store <- array(NA,c(n*n,n.gens)) # to save synchrony distribution
  ps <- rep(NA, n.gens) # population size
  ps.A <- rep(NA, n.gens) # population size
  ps.J <- rep(NA, n.gens) # population size
  G.store <- rep(NA, n.gens) # save genotype values
  E.store <- rep(NA, n.gens) # save environment values
  Z.store <- rep(NA, n.gens) # save phenotype values
  surv_store <- array(NA,c(n*n,n.gens)) # save phenotype values
  recru_store <- array(NA,c(n*n,n.gens)) # save phenotype values
  h_store <- array(NA,c(n*n,n.gens))# save phenotype values
  d_store <- array(NA,c(n*n,n.gens)) # save phenotype values
  
  # Initial distributions
  g <- seq(g.range[1],g.range[2],length.out = n) # set up range of Gs
  e <- seq(g.range[1],g.range[2],length.out = n) # set up range of Es
  G <- rep(g, each=n) 
  E <- rep(e, times=n)
  z.cols <- cbind(G,E) # binds together all combinations of E and G to create Zs
  G.count <- rep(1:n, each=n) # used later to reference individual Gs or Es
  E.count <- rep(1:n, times=n)

  # scale the environmental variables for cross val
  if(method == "CV"||method == "Predict"){
  if(method == "CV"){test_years <- c(test_years[1]-1, test_years)}
  marker_yr <- c(which(CV_data$Year == min(test_years)),which(CV_data$Year == max(test_years)))
  enviro <- CV_data[(marker_yr[1]:marker_yr[2]),]
  hatch <- (enviro$R_Mean_hatch-mean(descale$April_hatch))/sd(descale$April_hatch)
  ST.t <- (enviro$s_temp-mean(descale$Spring_temp))/sd(descale$Spring_temp)
  WP.t <- (enviro$w_precip-mean(descale$winter_precip_t))/sd(descale$winter_precip_t)
  WT.t <- (enviro$w_temp-mean(descale$winter_temp))/sd(descale$winter_temp)
  SP.t <- (enviro$s_precip-mean(descale$spring_precip_t))/sd(descale$spring_precip_t)
  B.t <- enviro$Beech
  if(decoupled == TRUE){
    ST.t_max <- (enviro$s_temp_max-mean(descale$Spring_temp_max))/sd(descale$Spring_temp_max)
    ST.t_cat <- (enviro$s_temp_cat-mean(descale$Spring_temp_cat))/sd(descale$Spring_temp_cat)
    if(weather_var=="ST"){ST.t_max <- rnorm(100,mean(datafile$Spring_temp_max),sd(datafile$Spring_temp_max))
    ST.t_cat <- rnorm(100,mean(datafile$Spring_temp_cat),sd(datafile$Spring_temp_cat))
    ST.t <- rnorm(100,mean(datafile$Spring_temp),sd(datafile$Spring_temp))}
    if(weather_var=="SP"){SP.t <- rnorm(100,mean(datafile$spring_precip_t),sd(datafile$spring_precip_t))}
    if(weather_var=="WT"){WT.t <- rnorm(100,mean(datafile$winter_temp),sd(datafile$winter_temp))}
    if(weather_var=="WP"){WP.t <- rnorm(100,mean(datafile$winter_precip_t),sd(datafile$winter_precip_t))}
  }}
  
  if(method == "Perturb"){
  ST.t <- rep(CV_data[1], n.gens)
  if(decoupled == TRUE){
    ST.t_max <- rep(CV_data[1], n.gens)
    ST.t_cat <- rep(CV_data[6], n.gens)}
  WP.t <- rep(CV_data[2], n.gens)
  WT.t <- rep(CV_data[3], n.gens)
  SP.t <- rep(CV_data[4], n.gens)
  HF.t <- rep(CV_data[5], n.gens)
  B.t <- rep(2, n.gens)
  enviro <- data.frame(Pop_size = descale$pop_size[which(descale$Year==2010)[1]])
  hatch <- (descale$April_hatch[which(descale$Year==2010)[1]]-mean(descale$April_hatch))/sd(descale$April_hatch)}
  
  cs <- mean(datafile$Clutch_size)#0 # always using mean clutch size across data which = 0 
  B.ref.S <- c(0, s.params[10:11])
  if(type == "tot_pop"){B.ref.S <- c(0, s.params[11:12])}
  B.ref.R <- c(0, r.params[11:12])
  B.ref.HF <- c(0,0,0)
  
  # then create initial distribution
  N <- enviro$R_only_PS[1]
  if(type == "tot_pop"){N <- enviro$Pop_size[1]}
  if(method == "Predict"){N <- enviro$Pop_size[1]
  HF.t <- rep(NA, n.gens)}
  ps[1] <- N
  N.A <- round(N*0.45)
  N.J <- round(N*0.55)
  mu.A <- c(hatch[1]/2,hatch[1]/2) # from previous year
  mu.J <- c(hatch[1]/2,hatch[1]/2) # from previous year
  sigma <- matrix(c(h.params[length(h.params)-1],0,0,(h.params[length(h.params)])),2,2) 
  res.A[,1] <- dmvnorm(z.cols,mu.A,sigma)
  res.A[,1] <- res.A[,1]/sum(res.A[,1])
  res.J[,1] <- dmvnorm(z.cols,mu.J,sigma)
  res.J[,1] <- res.J[,1]/sum(res.J[,1])
  res[,1] <- ((res.A[,1]*N.A)+(res.J[,1]*N.J))/N # keep as proportions of whole population
  z <- seq(min(G+E), max(G+E), length.out = n) # meshpoints, created by me
  # predict half fall before based on enviro values
  if(method == "Perturb"){
    B <- B.t[1]
    ifelse(is.na(HF.t[1]) == T, 
              HF.t[1:50] <- c(HF.fun(hf.params[1], hf.params[2], hf.params[3], 
                                                 hf.params[4], hf.params[5], ST.t[2], WT.t[1], 
                                                 WP.t[1], SP.t[2], B.ref.HF[B+1]),rep(0,49)), 
           HF.t[1] <- HF.t[1])
    if(decoupled == TRUE){ifelse(is.na(HF.t[1]) == T, 
                                 HF.t[1:50] <- c(HF.fun(hf.params[1], hf.params[2], hf.params[3], 
                                                          hf.params[4], hf.params[5], ST.t_cat[2], WT.t[1], 
                                                          WP.t[1], SP.t[2], B.ref.HF[B+1]),rep(0,49)), 
                                 HF.t[1] <- HF.t[1])}
    }
  if(method == "Predict"){B <- B.t[1]
    HF.t[1] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], 
                                            hf.params[4], hf.params[5], ST.t[2], WT.t[1], 
                                            WP.t[1], SP.t[2], B.ref.HF[B+1])
    if(decoupled == TRUE){HF.t[1] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], 
                                            hf.params[4], hf.params[5], ST.t_cat[2], WT.t[1], 
                                            WP.t[1], SP.t[2], B.ref.HF[B+1])}}
  if(method == "CV"){HF.t <- (enviro$half_fall-mean(descale$Half_fall))/sd(descale$Half_fall)}
  
  s.store[,1] <- (G+E)-HF.t[1] # starting synchrony
  
    
    for (i in 2:n.gens){
      
      # cross validation
      # predicting test_years[2] initially from test_years[1]
      B <- B.t[i-1]
      ST1 <- ST.t[i-1] # lagged effect
      ST2 <- ST.t[i]
      WT <- WT.t[i-1]
      WP <- WP.t[i-1]
      SP1 <- SP.t[i-1] # lagged effect
      SP2 <- SP.t[i]
      if(decoupled==TRUE){
      ST2_HF <- ST.t_cat[i]
      ST2_max <- ST.t_max[i]}
      
      Ns <- (N-m.N)/sd.N
      s <- s.store[,i-1]
      
      # have 'extra' for all effects beyond density dependence
      if(type == "res_only"){
        # z, ps, s, s^2, sp, st, wt, Round, resident, B
        extra.S <- (s*s.params[4])+((s^2)*s.params[5])+(SP1*s.params[6])+
        (WT*s.params[7])+s.params[8]+s.params[9]+B.ref.S[B+1]}else{
        extra.S <- (s*s.params[4])+((s^2)*s.params[5])+(SP1*s.params[6])+
        (WT*s.params[7])+(ST1*s.params[8])+s.params[9]+s.params[10]+B.ref.S[B+1]  
        }
      
        # mu.S, mu.S2, mu.ST, Ns, A, Round
        extra.DJ <- (s*d.params[4])+((s^2)*d.params[5])+(ST2*d.params[6])+0+d.params[8]
        extra.DA <- (s*d.params[4])+((s^2)*d.params[5])+(ST2*d.params[6])+d.params[7]+d.params[8]
        if(decoupled == TRUE){
        extra.DJ <- (s*d.params[4])+((s^2)*d.params[5])+(ST2_max*d.params[6])+0+d.params[8]
        extra.DA <- (s*d.params[4])+((s^2)*d.params[5])+(ST2_max*d.params[6])+d.params[7]+d.params[8]  
        }
        # s, s2, sp, cs, cs2, Round, B
        extra.R <- (s*r.params[4])+((s^2)*r.params[5])+(WP*r.params[6])+(SP1*r.params[7])+(cs*r.params[8])+
          ((cs^2)*r.params[9])+r.params[10]+B.ref.R[B+1]
        # mu.ST, mu.SP, mu.WT, mu.WP, Round
        extra.H <- (ST2*h.params[3])+(SP2*h.params[4])+(WP*h.params[5])+(WT*h.params[6])+h.params[7]
      
      ### S%D ###
      # Want to keep survival first - fit phenotype to the current environment and pop size
      # estimate survival for all combinations but then multiply fitness by proportion of those in the population
      surv <- S.fun(G+E, s.params[1], s.params[2], s.params[3], extra=extra.S, Ns) 
      surv_store[,i] <- surv
      res.fit.J[,i-1] <- surv*res.J[,i-1] # multiply chance of survival by % of pop in each z
      res.fit.A[,i-1] <- surv*res.A[,i-1] # multiply chance of survival by % of pop in each z

      for(j in G.count[1]: length(unique(G.count))){
        marker <- which(G.count == unique(G.count)[j])
        D.J <- t(outer(G[marker]+e,G[marker]+e, D.fun, d.params[1], d.params[2], d.params[3], extra = extra.DJ, sigma.d = d.params[length(d.params)], Ns)) # define for all z values
        D.J <- D.J/matrix(as.vector(apply(D.J,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        res.d.J[marker, i-1] <- D.J%*%res.fit.J[marker,i-1]
      } # does not change G value

      for(j in G.count[1]: length(unique(G.count))){
        marker <- which(G.count == unique(G.count)[j])
        D <- t(outer(G[marker]+e,G[marker]+e, D.fun, d.params[1], d.params[2], d.params[3], extra = extra.DA, sigma.d = d.params[length(d.params)], Ns))
        D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        res.d.A[marker, i-1] <- D%*%res.fit.A[marker,i-1]
      } # does not change G value
      d_store[,i] <- D
      
      res.d[,i-1] <- ((res.d.J[,i-1]*N.J)+(res.d.A[,i-1]*N.A))/N # combined to create adult population at t+1, both on scale of whole pop so ok.
      
      #CHECKS
      #moments.fun(res.A[,i-1],G+E)[1:2]
      #moments.fun(res.fit.A[,i-1],G+E)[1:2]
      
      #moments.fun(res.A[,i-1],E)[1:2]
      #moments.fun(res.fit.A[,i-1],E)[1:2]
      
      ### R%I ###
      # now recruitment
      recru <- R.fun(G+E, r.params[1], r.params[2], r.params[3], extra=extra.R, Ns)
      recru_store[,i] <- recru
      #plot(recru ~ G)
      res.r[,i-1] <- ((recru*(res[,i-1]*N))/N) # actual recruitment for population scaled to be proportion of pop in each z
      # have proportion of new recruits, want to assign Es, 
      # Do this the same way as development
      
      ## BIT TO KEEP VA SAME
      # set mean hatch date for G and E - beginning of new distribution of recruits
      mu.R <- c(moments.fun(res.r[,i-1],G)[1],moments.fun(res.r[,i-1],E)[1:2])
      R.ps <- sum(res.r[,i-1])*N
      sigma.R <- matrix(c(h.params[length(h.params)-1],0,0,mu.R[3]),2,2) 
      res.r[,i-1] <- dmvnorm(z.cols,mu.R[1:2],sigma.R)
      res.r[,i-1] <- (res.r[,i-1]/sum(res.r[,i-1]))*(R.ps/N) # scaled by num recruits/total N
      
      #CHECKS
      #moments.fun(res[,i-1],G+E)[1:2]
      #moments.fun(res.r[,i-1],G+E)[1:2]
      
      for(k in G.count[1]: length(unique(G.count))){
        marker <- which(G.count == unique(G.count)[k])
        H <- t(outer(G[marker]+e,G[marker]+e, H.fun, h.params[1], h.params[2], sigma.h = h.params[length(h.params)], extra = extra.H, Ns)) 
        H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        res.h[marker, i-1] <- H%*%res.r[marker,i-1]
      }
      
      h_store[,i] <- H
      
      # Now combine everything
      # (S%*%D) + (R%*%H)
      res[,i] <- res.d[,i-1] + res.h[,i-1] # now have complete population at t+1
      res.A[,i] <- res.d[,i-1]
      res.J[,i] <- res.h[,i-1] # all as proportions of whole population
      
      ps[i] <- sum(res[,i])*N # population size
      res[,i] <- res[,i]/sum(res[,i]) # save next distribution scaled to 1
      N.A <- ps.A[i] <- sum(res.A[,i])*N
      res.A[,i] <- res.A[,i]/sum(res.A[,i])
      N.J <- ps.J[i] <- sum(res.J[,i])*N
      res.J[,i] <- res.J[,i]/sum(res.J[,i]) # scale these back to proportions of new pop size
      N <- (ps[i]/0.58) # save new pop size for next t
      
      # set up HF and synchrony for t+1 
      if(method == "Predict"){
      if(decoupled == TRUE){HF.t[i] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], hf.params[4], hf.params[5], ST2_HF, WT, WP, SP2, B.ref.HF[B+1])}else{
        HF.t[i] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], hf.params[4], hf.params[5], ST2, WT, WP, SP2, B.ref.HF[B+1])
      }}
      s <- (G+E) - HF.t[i] # checked and it is ok to calculate using scaled values
      if(method == "Perturb"){ s <- (s - HF.t[i])}
      s.store[,i] <- s 
    }
    
    # need to save out moments from each run
    moment.res.s <- moment.res.g <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
    for (m in 1:n.gens){
      moment.res.g[,m] <- moments.fun(res[,m],G)[1:4]
      moment.res.z[,m] <- moments.fun(res[,m],G+E)[1:4]
      moment.res.e[,m] <- moments.fun(res[,m],E)[1:4]
    }
    
    if(method == "CV"){return(list(data.frame(PS = ps,
                      Hatch_mean = moment.res.z[1,],
                      Hatch_var = moment.res.z[2,], 
                      G_mean = moment.res.g[1,],
                      E_mean = moment.res.e[1,],
                      HF = HF.t), 
                      S = surv_store,
                      R = recru_store,
                      H = h_store,
                      D = d_store))}
    if(method == "Perturb"){return(ps[50])}
    if(method == "Predict"){return(list(data.frame(PS = ps,
                                         Hatch_mean = moment.res.z[1,],
                                         Hatch_var = moment.res.z[2,],
                                         G_mean = moment.res.g[1,],
                                         G_var = moment.res.g[2,],
                                         E_mean = moment.res.e[1,],
                                         E_var = moment.res.e[2,],
                                         HF = HF.t),
                                         S = surv_store,
                                         R = recru_store,
                                         H = h_store,
                                         D = d_store))}

}

run_CV_IPM_PO <- function(n, n.gens, g.range, m.N, sd.N, CV_data, 
                          test_years, descale, 
                          s.params, d.params, r.params, h.params,
                          hf.params,
                          method=c("CV", "Predict")){
  # num classes = n
  # num iterations = n.gens
  # min and max G and E values = g.range
  # mean and sd pop size for scaling = m.N, sd.N
  # cross validation data = CV_data
  # years we are predicting = test_years
  # data for scaling = descale
  # parameters for functions = s,d,r,h.params 
  # type = either cross validation or prediction
  
  # results vectors 
  res.A <- array(NA,c(n,n.gens)) # adults
  res.J <- array(NA,c(n,n.gens)) # juveniles
  res <- array(NA,c(n,n.gens)) # overall results
  res.fit.A <- array(NA,c(n,n.gens)) # storage vector to save fitness at each stage
  res.fit.J <- array(NA,c(n,n.gens)) # storage vector to save fitness at each stage
  res.d <- array(NA,c(n,n.gens)) # D%S
  res.d.A <- array(NA,c(n,n.gens)) # D%S
  res.d.J <- array(NA,c(n,n.gens)) # D%S
  res.r <- array(NA,c(n,n.gens)) # storage vector for recruitment
  res.h <- array(NA,c(n,n.gens)) # R%H
  s.store <- array(NA,c(n,n.gens)) # to save synchrony distribution
  ps <- rep(NA, n.gens) # population size
  ps.A <- rep(NA, n.gens) # population size
  ps.J <- rep(NA, n.gens) # population size
  Z.store <- rep(NA, n.gens) # save phenotype values
  surv_store <- array(NA,c(n,n.gens)) # save phenotype values
  recru_store <- array(NA,c(n,n.gens)) # save phenotype values
  h_store <- array(NA,c(n,n.gens))# save phenotype values
  d_store <- array(NA,c(n,n.gens)) # save phenotype values
  
  # Initial distributions
  z <- seq(g.range[1],g.range[2],length.out = n) # set up range of Gs
  
  # scale the environmental variables for cross val
  B.ref.S <- c(0, s.params[11:12])
  B.ref.R <- c(0, r.params[11:12])
  B.ref.HF <- c(0,0,0)
  if(method == "CV"){test_years <- c(test_years[1]-1, test_years)}
  marker_yr <- c(which(CV_data$Year == min(test_years)),which(CV_data$Year == max(test_years)))
  enviro <- CV_data[(marker_yr[1]:marker_yr[2]),]
  hatch <- (enviro$R_Mean_hatch-mean(descale$April_hatch))/sd(descale$April_hatch)
  ST.t <- (enviro$s_temp-mean(descale$Spring_temp))/sd(descale$Spring_temp)
  WP.t <- (enviro$w_precip-mean(descale$winter_precip_t))/sd(descale$winter_precip_t)
  WT.t <- (enviro$w_temp-mean(descale$winter_temp))/sd(descale$winter_temp)
  SP.t <- (enviro$s_precip-mean(descale$spring_precip_t))/sd(descale$spring_precip_t)
  B.t <- enviro$Beech
  HF.t <- (enviro$half_fall-mean(descale$Half_fall))/sd(descale$Half_fall)
  if(method == "Predict"){HF.t <- rep(NA, n.gens)
  B <- B.t[1]
  HF.t[1] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], 
                    hf.params[4], hf.params[5], ST.t[2], WT.t[1], 
                    WP.t[1], SP.t[2], B.ref.HF[B+1])}
  cs <- mean(datafile$Clutch_size)#0 # always using mean clutch size across data which = 0 

  
  # then create initial distribution
  N <- enviro$Pop_size[1]

  ps[1] <- N
  N.A <- round(N*0.45)
  N.J <- round(N*0.55)
  mu.A <- c(hatch[1]) # from previous year
  mu.J <- c(hatch[1]) # from previous year
  sigma <- 0.61 
  z <- seq(min(g.range), max(g.range), length.out = n) # meshpoints, created by me
  res.A[,1] <- dnorm(z,mu.A,sigma)
  res.A[,1] <- res.A[,1]/sum(res.A[,1])
  res.J[,1] <- dnorm(z,mu.J,sigma)
  res.J[,1] <- res.J[,1]/sum(res.J[,1])
  res[,1] <- ((res.A[,1]*N.A)+(res.J[,1]*N.J))/N # keep as proportions of whole population
  
  s.store[,1] <- z-HF.t[1] # starting synchrony
  HF.store <- rep(NA, n.gens)
  HF.store[1] <- HF.t[1]
  
  
  for (i in 2:n.gens){
    
    # cross validation
    # predicting test_years[2] initially from test_years[1]
    B <- B.t[i-1]
    ST1 <- ST.t[i-1] # lagged effect
    ST2 <- ST.t[i]
    WT <- WT.t[i-1]
    WP <- WP.t[i-1]
    SP1 <- SP.t[i-1] # lagged effect
    SP2 <- SP.t[i]
    
    Ns <- (N-m.N)/sd.N
    s <- s.store[,i-1]
    
    # have 'extra' for all effects beyond density dependence
    extra.S <- (s*s.params[4])+((s^2)*s.params[5])+(SP1*s.params[6])+
            (WT*s.params[7])+(ST1*s.params[8])+s.params[9]+s.params[10]+B.ref.S[B+1]  
    
    # mu.S, mu.S2, mu.ST, Ns, A, Round
    extra.DJ <- (s*d.params[4])+((s^2)*d.params[5])+(ST2*d.params[6])+0+d.params[8]
    extra.DA <- (s*d.params[4])+((s^2)*d.params[5])+(ST2*d.params[6])+d.params[7]+d.params[8]

    # s, s2, sp, cs, cs2, Round, B
    extra.R <- (s*r.params[4])+((s^2)*r.params[5])+(WP*r.params[6])+(SP1*r.params[7])+(cs*r.params[8])+
      ((cs^2)*r.params[9])+r.params[10]+B.ref.R[B+1]
    # mother z, ST, SP, WT, WP, Round
    extra.H <- (ST2*h.params[4])+(SP2*h.params[5])+h.params[6]
    
    ### S%D ###
    # Want to keep survival first - fit phenotype to the current environment and pop size
    # estimate survival for all combinations but then multiply fitness by proportion of those in the population
    surv <- S.fun(z, s.params[1], s.params[2], s.params[3], extra=extra.S, Ns) 
    res.fit.J[,i-1] <- surv*res.J[,i-1] # multiply chance of survival by % of pop in each z
    res.fit.A[,i-1] <- surv*res.A[,i-1] # multiply chance of survival by % of pop in each z
    # do across all z here
    surv_store[,i] <- surv
    D.J <- (t(outer(z,z,D.fun, d.params[1], d.params[2], d.params[3], 
                    extra = extra.DJ, sigma.d = d.params[length(d.params)], Ns)))
    # Here run on all values not split by G
    D.J <- D.J/matrix(as.vector(apply(D.J,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    res.d.J[, i-1] <- D.J%*%res.fit.J[,i-1]
    D <- t(outer(z,z, D.fun, d.params[1], d.params[2], d.params[3],  
                 extra = extra.DA, sigma.d = d.params[length(d.params)], Ns))
    D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    res.d.A[, i-1] <- D%*%res.fit.A[,i-1]
    res.d[,i-1] <- ((res.d.J[,i-1]*N.J)+(res.d.A[,i-1]*N.A))/N # combined to create adult population at t+1, both on scale of whole pop so ok.
    d_store[,i] <- D[,1]
    
    ### R%I ###
    # now recruitment
    recru <- R.fun(z, r.params[1], r.params[2], r.params[3], extra=extra.R, Ns)
    res.r[,i-1] <- ((recru*(res[,i-1]*N))/N) # actual recruitment for population scaled to be proportion of pop in each z
    recru_store[,i] <- recru
    # have proportion of new recruits, want to assign Es, G assignment is clonal
    # Do this the same way as development
    H <- t(outer(z, z, H.fun.PO, mu.int=h.params[1], mu.N=h.params[2], mu.mZ=h.params[3], sigma.h = h.params[length(h.params)],
                 extra = extra.H, Ns)) 
    H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    res.h[, i-1] <- H%*%res.r[,i-1]
    h_store[,i] <- H[,1]
    
    # Now combine everything
    # (S%*%D) + (R%*%H)
    res[,i] <- res.d[,i-1] + res.h[,i-1] # now have complete population at t+1
    res.A[,i] <- res.d[,i-1]
    res.J[,i] <- res.h[,i-1] # all as proportions of whole population
    
    ps[i] <- sum(res[,i])*N # population size
    res[,i] <- res[,i]/sum(res[,i]) # save next distribution scaled to 1
    N.A <- ps.A[i] <- sum(res.A[,i])*N
    res.A[,i] <- res.A[,i]/sum(res.A[,i])
    N.J <- ps.J[i] <- sum(res.J[,i])*N
    res.J[,i] <- res.J[,i]/sum(res.J[,i]) # scale these back to proportions of new pop size
    N <- (ps[i]/0.58) # save new pop size for next t
    
    # set up HF and synchrony for t+1 
    HF <- HF.t[i]
    
    if(method=="Predict"){HF <- HF.fun(hf.params[1], hf.params[2], 
                                       hf.params[3], hf.params[4], hf.params[5], 
                                       ST2, WT, WP, SP2, B.ref.HF[B+1])}
    s <- z - HF # checked and it is ok to calculate using scaled values
    s.store[,i] <- s
  }
  
  # need to save out moments from each run
  moment.res.s <- moment.res.g <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
  for (m in 1:n.gens){
    moment.res.z[,m] <- moments.fun(res[,m],z)[1:4]
  }
  
  return(list(data.frame(PS = ps,
                    Hatch_mean = moment.res.z[1,],
                    Hatch_var = moment.res.z[2,],
                    HF = HF.t),
                    S = surv_store,
                    R = recru_store,
                    H = h_store,
                    D = d_store
                    ))
  
}

#### plot cross val ####

for_plotting_CV <- function(data_list, test_years,
                            descale = descale){
  datafile <- data_list[[1]][2:6,]
  datafile$Year <- test_years[[1]]
  for(i in 2:length(data_list)){
    data_list[[i]]$Year <- c(NA,test_years[[i]])
    datafile <- rbind(datafile, data_list[[i]][2:6,])
  }
  datafile$Hatch_mean <- (datafile$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)
  datafile$Hatch_var <- datafile$Hatch_var*var(descale$April_hatch)
  return(datafile)
}

#### plots IPM2 ####

plotter <- function(PS.store, y_label, x_label, main_label, x.axis=TRUE, y.axis=TRUE, adj,
                    ylim=c(0,450)){
  plot(1:91, PS.store[[1]]$PS, type="l", xlab="", ylab="", 
       main="",ylim=ylim, axes=F, col='grey', las=1, cex.lab=1.5)
  if(y.axis == TRUE){axis(2, las=1, cex.axis=1.5)}
  if(x.axis==TRUE){axis(1, at=c(1,91), labels=c("",2100), cex.axis=0.75)}
  title(main=main_label, adj=adj, cex.main=2)
  title(xlab=x_label, line=2.5, cex.lab=1.5)
  title(ylab=y_label, line=3.5, cex.lab=1.5)
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){
    test <- which(PS.store[[i]]$PS<1)
    if(length(test)>0){lines(1:test[1], PS.store[[i]]$PS[1:test[1]], type='l', col=alpha("grey",0.2))
      extinctions[i,] <- c(test[1], PS.store[[i]]$PS[test[1]])}
      else{lines(1:91, PS.store[[i]]$PS, type='l', col=alpha("grey",0.2))}}
  points(extinctions[,1], extinctions[,2], pch=4, cex=1, col=alpha("black",1))
  abline(h=131, col="black", lty=1)
  abline(h=219, col="black", lty=2)
  abline(h=40, col="black", lty=2)
  
}

plotter2 <- function(results, y_label, x_label, main_label, datafile, x.axis=TRUE, y.axis=FALSE){
  plot(1:91,seq(1,75,length.out = 91), type='n', xlab="", ylab="", ylim=c(5,70), axes=F, cex.lab=1.5)
  if(y.axis==TRUE){axis(2, las=1, at=c(5,30),labels=rep("",2), tck=0.05)
  axis(2, las=1, at=c(35,70),labels=rep("",2), tck=0.05)
  axis(2, las=1, at=seq(10,30,10), labels=seq(10,30,10), cex.axis=1.5)
  axis(2, las=1, at=seq(40,70,10), labels=seq(10,40,10), cex.axis=1.5)}
  if(x.axis==TRUE){axis(1, at=c(1,91), labels=c("",2100), cex.axis=0.75)}
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){
    ds_E.store <- (results[[i]]$E_mean*sd(datafile$April_hatch))+(mean(datafile$April_hatch)/2)
    test <- which(results[[i]]$PS < 1)
    if(length(test)>0){lines(1:test[1], ds_E.store[1:test[1]]+30, type='l', col=alpha("purple3",0.05))
      extinctions[i,] <- c(test[1], ds_E.store[test[1]]+30)}
      else{lines(1:91, ds_E.store+30, type='l', col=alpha("purple3",0.05))}
  }
  #points(extinctions[,1], extinctions[,2], pch=4, cex=0.75, col=alpha("black",1))
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){ds_G.store <- (results[[i]]$G_mean*sd(datafile$April_hatch))+(mean(datafile$April_hatch)/2)
    test <- which(results[[i]]$PS < 1)
    if(length(test)>0){lines(1:test[1], ds_G.store[1:test[1]], type='l', col=alpha("darkorange",0.05))
    extinctions[i,] <- c(test[1], ds_G.store[test[1]])}
    else{lines(1:91, ds_G.store, type='l', col=alpha("darkorange",0.05))}
  }
  #points(extinctions[,1], extinctions[,2], pch=4, cex=0.75, col=alpha("black",1))
  abline(h=20, col="grey70", lty=2)
  abline(h=(20+30), col="grey70", lty=2)
  title(xlab=x_label, line=2.5, cex.lab=1.5)
  mtext(y_label, side=2, line=3.5, cex=1.5)
  title(main=main_label, cex.main=2)
}

plotter3 <- function(results, y_label, x_label, main_label, datafile, x.axis=TRUE, y.axis=FALSE, adj=0, adj2=0.5,
                     ylim=ylim){
  ds_Z.store <- (results[[1]]$Hatch_mean*sd(datafile$April_hatch))+mean(datafile$April_hatch)
  plot(1:91,seq(5,55,length.out = 91), type='n', xlab="", ylab="", ylim=c(10,60), axes=F)
  if(y.axis==TRUE){axis(2, las=1, at=c(10,60),labels=rep("",2), tck=0.05)
    axis(2, las=1, at=seq(20,50,10), labels=seq(20,50,10), cex.axis=1.5)}
  if(x.axis==TRUE){axis(1, at=c(1,91), labels=c("",2100), cex.axis=0.75)}
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){
    ds_Z.store <- (results[[i]]$Hatch_mean*sd(datafile$April_hatch))+mean(datafile$April_hatch)
    test <- which(results[[i]]$PS < 1)
    if(length(test)>0){lines(1:test[1], ds_Z.store[1:test[1]], type='l', col=alpha("darkgreen",0.05))
      extinctions[i,] <- c(test[1], ds_Z.store[test[1]])}
    else{lines(1:91, ds_Z.store[1:91], type='l', col=alpha("darkgreen",0.01))}
  }
  #points(extinctions[,1], extinctions[,2], pch=4, cex=0.75, col=alpha("black",1))
  abline(h=40, col="black", lty=2)
  mtext(x_label, line=2.5, side=1, cex=1.5, adj=adj2)
  title(ylab=y_label, line=3.5, cex.lab=1.5)
  title(main=main_label, adj=adj, cex.main=2)
}

plotter4 <- function(SYNC, y_label, x_label, main_label, datafile, x.axis=TRUE, y.axis=TRUE, adj=0){
  plot(1:91, SYNC[[1]]$SYNC, type='n', xlab="", ylab="", ylim=c(-30,25), axes=F, cex.lab=1.5)
  if(x.axis==TRUE){axis(1, at=c(1,91), labels=c("",2100), cex.axis=0.75)}
  if(y.axis==TRUE){axis(2, las=1, at=seq(-30,25,5), cex.axis=1.5)}
  extinctions <- cbind(rep(NA,1000), rep(NA,1000))
  for(i in 1:1000){
    test <- which(SYNC[[i]]$PS<1)
    if(length(test)>0){lines(1:test[1], SYNC[[i]]$SYNC[1:test[1]], type='l', col=alpha("green",0.05))
      extinctions[i,] <- c(test[1], SYNC[[i]]$SYNC[test[1]])}
    else{lines(1:91, SYNC[[i]]$SYNC, type='l', col=alpha("green",0.05))}
  }
  extinctions <- extinctions[!is.na(extinctions[,1]),]
  points(extinctions[,1], extinctions[,2], pch=4, cex=1, col=alpha("black",1))
  abline(h=-12, col="black", lty=2)
  title(ylab=y_label, line=3.5, cex.lab=1.5)
  title(xlab=x_label, line=2.5, cex.lab=1.5)
  title(main=main_label, adj=adj, cex.main=2)
}

output_predictions_ps <- function(PS.store){
  PGR <- rep(NA,1000)
  SDs <- rep(NA,1000)
 for(j in 1:1000){PGR[j] <- mean((PS.store[[j]]$PS[2:91]-PS.store[[j]]$PS[1:90])/PS.store[[j]]$PS[1:90])
 SDs[j] <- sd((PS.store[[j]]$PS[2:91]-PS.store[[j]]$PS[1:90])/PS.store[[j]]$PS[1:90])}
  PGR_O <- sort(PGR)
  SD_O <- sort(SDs)
  return(cbind(PGR_O, SD_O))
}

output_predictions_z <- function(PS.store, descale){
  slope <- rep(NA,1000)
  SDs <- rep(NA,1000)
  for(j in 1:1000){PS.store[[j]]$Hatch_mean <- (PS.store[[j]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)
    slope[j] <- mean(PS.store[[j]]$Hatch_mean[2:91]-PS.store[[j]]$Hatch_mean[1:90])
    SDs[j] <- sd(PS.store[[j]]$Hatch_mean[2:91]-PS.store[[j]]$Hatch_mean[1:90])}
  slope_O <- sort(slope)
  SD_O <- sort(SDs)
  return(cbind(slope_O, SD_O))
}

unscale <- function(x, descale, variable){
  newx <- (x*sd(descale[,which(colnames(descale)==variable)]))+
    mean(descale[,which(colnames(descale)==variable)])
  return(newx)
}

Mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}
