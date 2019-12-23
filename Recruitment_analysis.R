#### Script to run model selection for the recruitment function

# COLLAPSE ALL FOR EASIER VIEWING

library(dplyr)
library(lme4)

#### Import data scaled data ####

datafile <- read.csv("bio_data.csv", header=T)

#### Run recruitment models ####

## Begin with the null model
null_model <- glmer(Num_recruited ~ (1|Year) + (1|Nest_box), family="poisson", data=datafile)
summary(null_model)

# For the full model the optimizer is customised to aid model estimation
#control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

## Then the full model
full_model <- glmer(Num_recruited ~ April_hatch +
                      I(April_hatch^2) +
                      Spring_temp +
                      Synchrony +
                      I(Synchrony^2) +
                      Clutch_size +
                      I(Clutch_size^2) + 
                      as.factor(Age) +
                      as.factor(Immigrant) + 
                      as.factor(beech) +
                      R_pop_size +
                      winter_temp +
                      winter_precip_t +
                      spring_precip_t +
                      as.factor(Round) +
                      (1|Year) + (1|Nest_box), family="poisson", data=datafile, 
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# Look at results
summary(full_model)

### MODEL SELECTION

# Create data frame to save AIC results in
output_AIC <- data.frame(AIC = rep(NA, 5), delta_AIC = rep(NA, 5), variable_removed = rep(NA, 5))
output_AIC[1,1] <- AIC(full_model)
output_AIC[1,3] <- "full model"

# Want to remove any variables with estimate < 2*SE
# so looking for negative values of x
temp <- as.data.frame(summary(full_model)$coef)
temp$x <- abs(temp[,1]) - (2*(temp[,2]))
temp[which(temp$x < 0),]

# Remove until no variables remain with estimate < 2*SE
# For our study the final model was:
red_1 <- update(full_model, .~. - I(April_hatch^2))
red_2 <- update(red_1, .~. - as.factor(Immigrant))
red_3 <- update(red_2, .~. - Spring_temp)
red_4 <- update(red_3, .~. - as.factor(Age))
red_5 <- update(red_4, .~. - winter_temp)
summary(red_5)

output_AIC[2,1] <- AIC(red_5)
output_AIC[2,2] <- output_AIC[2,1] - output_AIC[1,1]
output_AIC[2,3] <- "final model"
output_AIC

save(red_5, file="All_data_recruitment.RData")
