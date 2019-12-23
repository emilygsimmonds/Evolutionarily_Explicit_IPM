#### Script to run each demographic function 
#and the caterpillar timing function using Statistical_model_functions.R, 
# produces parameter values for each function for the cross validation datasets

## These are needed for K-fold cross validation

# COLLAPSE ALL FOR EASIER VIEWING

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
