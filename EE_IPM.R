#### Script containing the function to run the evolutionarily explicit/
# quant gen IPM

### INPUTS: n = number of classes for the phenotype
#           n.gens = number of generations
#           g.range = range of genetic values (approx. half total phenotype)
#           m.N = mean Population size (unscaled)
#           sd.N = standard deviation of Population size (unscaled)
#           g.range = range of genetic values (approx. half total phenotype)
#           CV_data = cross validation data
#           type = which variable to be used for population size
#           method = whether the model is simulating for cross validation, perturbation or prediction
#           test_years = which years are test years for cross validation
#           descale = unscaled biological data
#           s.params etc = parameters for each function (survival, development, recruitment, inheritance)
#           hf.params = caterpillar timing parameters
#           decoupled = TRUE or FALSE, whether caterpillars and great tits share a cue or not
#           weather_var = which weather variable is held at 1960-2010 levels or "NO"


run_EE_IPM <- function(n, n.gens, g.range,
                       m.N, sd.N, CV_data, type = c("tot_pop", "res_only"),
                       method = c("CV", "Perturb", "Predict"),
                       test_years, descale, s.params, d.params, r.params, h.params,
                       hf.params,
                       decoupled = FALSE, 
                       weather_var = "NO"){
  
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
  
  # Initial distributions
  g <- seq(g.range[1],g.range[2],length.out = n) # set up range of Gs
  e <- seq(g.range[1],g.range[2],length.out = n) # set up range of Es
  G <- rep(g, each=n) 
  E <- rep(e, times=n)
  z.cols <- cbind(G,E) # binds together all combinations of E and G to create Zs
  G.count <- rep(1:n, each=n) # used later to reference individual Gs or Es
  E.count <- rep(1:n, times=n)
  
  # scale the environmental variables for cross validation input
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
    
    # If decoupled cues pull spring temperature for each species
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
  
  # If perturbing, the environmental drivers remain the same across the simulation
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
  
  cs <- mean(datafile$Clutch_size) # always using mean clutch size across data which = 0 
  B.ref.S <- c(0, s.params[10:11])
  if(type == "tot_pop"){B.ref.S <- c(0, s.params[11:12])} # different survival parameters for total population
  B.ref.R <- c(0, r.params[11:12])
  B.ref.HF <- c(0,0,0)
  
  # then create initial distribution
  N <- enviro$R_only_PS[1]
  if(type == "tot_pop"){N <- enviro$Pop_size[1]}
  if(method == "Predict"){N <- enviro$Pop_size[1]
  HF.t <- rep(NA, n.gens)}
  ps[1] <- N # population size
  N.A <- round(N*0.45) # assign a proportion of pop to adults or juveniles (taken from 1960-2010 data)
  N.J <- round(N*0.55)
  mu.A <- c(hatch[1]/2,hatch[1]/2) # set up mean hatch date, take from previous year of data
  mu.J <- c(hatch[1]/2,hatch[1]/2) 
  # set up variance of hatch date based on phenotypic variance and
  # additive genetic variance
  sigma <- matrix(c(h.params[length(h.params)-1],0,0,(h.params[length(h.params)])),2,2) 
  res.A[,1] <- dmvnorm(z.cols,mu.A,sigma) # create initial adult distribution
  res.A[,1] <- res.A[,1]/sum(res.A[,1]) # and scale so it sums to 1
  res.J[,1] <- dmvnorm(z.cols,mu.J,sigma)
  res.J[,1] <- res.J[,1]/sum(res.J[,1])
  res[,1] <- ((res.A[,1]*N.A)+(res.J[,1]*N.J))/N # keep as proportions of whole population
  z <- seq(min(G+E), max(G+E), length.out = n) # meshpoints
  
  # predict caterpillar timing based on environmental values
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
    
    # Assign environmental driver values for t+1
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
    
    Ns <- (N-m.N)/sd.N # scale population size to be used in demographic functions
    s <- s.store[,i-1] # store synchrony
    
    # Calculate the effect of environmental drivers for different functions.
    # Have 'extra' for all effects beyond density dependence
    # allows for flexibility in the functions
    
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
    # Estimate survival for all Z values and then multiply fitness by proportion of those in the population
    surv <- S.fun(G+E, s.params[1], s.params[2], s.params[3], extra=extra.S, Ns) 
    res.fit.J[,i-1] <- surv*res.J[,i-1] # multiply chance of survival by % of pop in each z
    res.fit.A[,i-1] <- surv*res.A[,i-1] # multiply chance of survival by % of pop in each z
    
    # Estimate development of phenotype for all Z values. 
    # Estimates chance of change to all Z values from all other Z values based on environmental conditions
    # First, split to do this WITHIN each G value
    for(j in G.count[1]: length(unique(G.count))){
      marker <- which(G.count == unique(G.count)[j]) # find the focal G value
      # Estimate development
      D.J <- t(outer(G[marker]+e,G[marker]+e, D.fun, d.params[1], d.params[2], d.params[3], extra = extra.DJ, sigma.d = d.params[length(d.params)], Ns)) # define for all z values
      # scale to sum to the original proportion of the population with that G value
      D.J <- D.J/matrix(as.vector(apply(D.J,2,sum)),nrow=n,ncol=n,byrow=TRUE) 
      res.d.J[marker, i-1] <- D.J%*%res.fit.J[marker,i-1]
    } # does not change G value
    
    # Repeat above for juveniles
    for(j in G.count[1]: length(unique(G.count))){
      marker <- which(G.count == unique(G.count)[j])
      D <- t(outer(G[marker]+e,G[marker]+e, D.fun, d.params[1], d.params[2], d.params[3], extra = extra.DA, sigma.d = d.params[length(d.params)], Ns))
      D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      res.d.A[marker, i-1] <- D%*%res.fit.A[marker,i-1]
    } # does not change G value
    
    # combine to create the adult population at t+1, both on scale of whole pop so works.
    res.d[,i-1] <- ((res.d.J[,i-1]*N.J)+(res.d.A[,i-1]*N.A))/N 
    
    ### R%H ###
    # Estimate recruitment for each Z value
    recru <- R.fun(G+E, r.params[1], r.params[2], r.params[3], extra=extra.R, Ns)
    # actual recruitment for population scaled to be proportion of pop in each z
    res.r[,i-1] <- ((recru*(res[,i-1]*N))/N) 
    
    ## Assign G values to new recruits
    # set mean hatch date for G and E - beginning of new distribution of recruits
    # the mean comes from the selected parents
    mu.R <- c(moments.fun(res.r[,i-1],G)[1],moments.fun(res.r[,i-1],E)[1:2])
    R.ps <- sum(res.r[,i-1])*N # calculate population size of recruits
    # Sigma is set by additive genetic variance and environmental variance of new recruits
    # so variance of E does not change
    sigma.R <- matrix(c(h.params[length(h.params)-1],0,0,mu.R[3]),2,2) 
    res.r[,i-1] <- dmvnorm(z.cols,mu.R[1:2],sigma.R)
    res.r[,i-1] <- (res.r[,i-1]/sum(res.r[,i-1]))*(R.ps/N) # scaled by num recruits/total N
    
    # Next, redistribute E values within in each G value
    for(k in G.count[1]: length(unique(G.count))){
      marker <- which(G.count == unique(G.count)[k]) # split by each G value
      # Estimate probability of changing from one Z value to another
      H <- t(outer(G[marker]+e,G[marker]+e, H.fun, h.params[1], h.params[2], sigma.h = h.params[length(h.params)], extra = extra.H, Ns)) 
      # Scale to remain the same proportion og the population
      H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      res.h[marker, i-1] <- H%*%res.r[marker,i-1]
    }
    
    # Now combine everything
    # (S%*%D) + (R%*%H)
    res[,i] <- res.d[,i-1] + res.h[,i-1] # now have complete population at t+1
    res.A[,i] <- res.d[,i-1]
    res.J[,i] <- res.h[,i-1] # all as proportions of whole population
    
    ps[i] <- sum(res[,i])*N # population size
    res[,i] <- res[,i]/sum(res[,i]) # save next distribution scaled to 1
    N.A <- ps.A[i] <- sum(res.A[,i])*N # new adult population
    res.A[,i] <- res.A[,i]/sum(res.A[,i]) # scaled to 1
    N.J <- ps.J[i] <- sum(res.J[,i])*N # new juvenile population
    res.J[,i] <- res.J[,i]/sum(res.J[,i]) # scale these back to proportions of new pop size
    N <- (ps[i]/0.58) # save new pop size for next t
    
    # set up caterpillar timing and synchrony for t+1 
    if(method == "Predict"){
      if(decoupled == TRUE){HF.t[i] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], hf.params[4], hf.params[5], ST2_HF, WT, WP, SP2, B.ref.HF[B+1])}else{
        HF.t[i] <- HF.fun(hf.params[1], hf.params[2], hf.params[3], hf.params[4], hf.params[5], ST2, WT, WP, SP2, B.ref.HF[B+1])
      }}
    s <- (G+E) - HF.t[i] 
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
                                            HF = HF.t)))}
  if(method == "Perturb"){return(ps[50])}
  if(method == "Predict"){return(list(data.frame(PS = ps,
                                                 Hatch_mean = moment.res.z[1,],
                                                 Hatch_var = moment.res.z[2,],
                                                 G_mean = moment.res.g[1,],
                                                 G_var = moment.res.g[2,],
                                                 E_mean = moment.res.e[1,],
                                                 E_var = moment.res.e[2,],
                                                 HF = HF.t)))}
  
}

