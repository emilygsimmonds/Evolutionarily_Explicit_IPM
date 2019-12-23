#### Script containing the function to run the standard IPM

### INPUTS: n = number of classes for the phenotype
#           n.gens = number of generations
#           g.range = range of genetic values (approx. half total phenotype)
#           m.N = mean Population size (unscaled)
#           sd.N = standard deviation of Population size (unscaled)
#           g.range = range of genetic values (approx. half total phenotype)
#           CV_data = cross validation data
#           method = whether the model is simulating for cross validation (CV) or prediction
#           test_years = which years are test years for cross validation
#           descale = unscaled biological data
#           s.params etc = parameters for each function (survival, development, recruitment, inheritance)
#           hf.params = caterpillar timing parameters

run_CV_IPM_PO <- function(n, n.gens, g.range, m.N, sd.N, CV_data, 
                          test_years, descale, 
                          s.params, d.params, r.params, h.params,
                          hf.params,
                          method=c("CV", "Predict")){
  
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
  
  # Initial distributions
  z <- seq(g.range[1],g.range[2],length.out = n) # set up range of Zs
  
  # scale the environmental variables for cross validation
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
  
  cs <- mean(datafile$Clutch_size)# always using mean clutch size across data which = 0 
  
  
  # then create initial distribution
  N <- enviro$Pop_size[1] # set population size
  ps[1] <- N # population size
  N.A <- round(N*0.45) # assign a proportion of pop to adults or juveniles (taken from 1960-2010 data)
  N.J <- round(N*0.55)
  mu.A <- c(hatch[1]/2,hatch[1]/2) # set up mean hatch date, take from previous year of data
  mu.J <- c(hatch[1]/2,hatch[1]/2) 
  # set up standard deviation based on the standard deviation from observed data
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
    
    # Assign environmental driver values for t+1
    B <- B.t[i-1]
    ST1 <- ST.t[i-1] # lagged effect
    ST2 <- ST.t[i]
    WT <- WT.t[i-1]
    WP <- WP.t[i-1]
    SP1 <- SP.t[i-1] # lagged effect
    SP2 <- SP.t[i]
    
    Ns <- (N-m.N)/sd.N # scale population size to be used in demographic functions
    s <- s.store[,i-1] # store synchrony
    
    # Calculate the effect of environmental drivers for different functions.
    # Have 'extra' for all effects beyond density dependence
    # allows for flexibility in the functions
    
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
    # Estimate survival for all Z values and then multiply fitness by proportion of those in the population
    surv <- S.fun(z, s.params[1], s.params[2], s.params[3], extra=extra.S, Ns) 
    res.fit.J[,i-1] <- surv*res.J[,i-1] # multiply chance of survival by % of pop in each z
    res.fit.A[,i-1] <- surv*res.A[,i-1] # multiply chance of survival by % of pop in each z
    
    # Estimate development of phenotype for all Z values. 
    # Estimates chance of change to all Z values from all other Z values based on environmental conditions
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
    # Estimate recruitment for each Z value
    recru <- R.fun(z, r.params[1], r.params[2], r.params[3], extra=extra.R, Ns)
    res.r[,i-1] <- ((recru*(res[,i-1]*N))/N) # actual recruitment for population scaled to be proportion of pop in each z
    recru_store[,i] <- recru

    # Estimate inheritance of phenotype for all Z values. 
    # Estimates chance of change to all Z values from all other Z values based on environmental conditions
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
    
    # set up caterpillar timing and synchrony for t+1 
    HF <- HF.t[i]
    
    if(method=="Predict"){HF <- HF.fun(hf.params[1], hf.params[2], 
                                       hf.params[3], hf.params[4], hf.params[5], 
                                       ST2, WT, WP, SP2, B.ref.HF[B+1])}
    s <- z - HF 
    s.store[,i] <- s
  }
  
  # need to save out moments from each run
  moment.res.s <- moment.res.g <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
  for (m in 1:n.gens){
    moment.res.z[,m] <- moments.fun(res[,m],z)[1:4]
  }
  
  return(data.frame(PS = ps,
                         Hatch_mean = moment.res.z[1,],
                         Hatch_var = moment.res.z[2,],
                         HF = HF.t))
  
}