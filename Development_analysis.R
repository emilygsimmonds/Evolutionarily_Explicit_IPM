#### Script to run model selection for the development function

# COLLAPSE ALL FOR EASIER VIEWING

library(dplyr)
library(lme4)

#### Import data scaled data and create development data ####

datafile <- read.csv("bio_data.csv", header=T)

# For this analysis need to restrict to only individuals
# observed in the year after the focal year

# Loops checks which individuals have a t+1 and removes those that do not
u_FID <- unique(datafile$F_ID)

# Re-assign the data as the original is needed again later
datafilet1 <- datafile

for(i in 1:length(u_FID)){
  temp_check <- subset(datafile, datafile$F_ID == as.character(u_FID[i]))
  temp_check <- temp_check[order(temp_check$Year),]
  # want to check for each entry for each individual
  for(j in 1:length(temp_check[,1])){
    # if there is no row above then we want to remove them from the datafile
    marker1 <- which(datafile$F_ID == as.character(u_FID[i]) & datafile$Year == temp_check$Year[1])
    temp_check <- temp_check[-1,]
    if(length(temp_check[,1]) < 1){datafile <- datafile[-marker1,]}
  }
}

# This leaves a datafile of individuals at t who will also have an entry at t+1
write.csv(datafile, "bio_data_t.csv", row.names=F)

## Now need to extract the trait values and environmental conditions from t+1
# for each individual

# Create a new datafile to store results in
new_datafile <- data.frame(F_ID = datafile$F_ID, April_hatch_t1 = rep(NA, length(datafile$F_ID)), 
                           Spring_temp_t1 = rep(NA, length(datafile$F_ID)), 
                           Spring_precip_t1 = rep(NA, length(datafile$F_ID)), 
                           Round_t1 = rep(NA, length(datafile$F_ID)), 
                           age_t1 = rep(NA, length(datafile$F_ID)), 
                           pop_size_t1 = rep(NA, length(datafile$F_ID)), 
                           R_pop_size_t1 = rep(NA, length(datafile$F_ID)))

# Loop extracts the t+1 values for each variable in the new_datafile
for(l in 1:length(datafile[,1])){
  # subset full datafile to just the F_ID needed
  # find which row is the focal year (t), then take the next one (t+1)
  temp1 <- subset(datafile, datafile$F_ID == as.character(datafile$F_ID[l]))
  temp1 <- temp1[order(temp1$Year),]
  marker <- which(temp1$Year == datafile$Year[l]) # mark where in the datafile the focal individual is
  # extract the covariate values
  new_datafile$April_hatch_t1[l] <- temp1$April_hatch[marker+1]
  new_datafile$Spring_temp_t1[l] <- temp1$Spring_temp[marker+1]
  new_datafile$Spring_precip_t1[l] <- temp1$spring_precip_t[marker+1]
  new_datafile$Round_t1[l] <- as.character(temp1$Round[marker+1])
  new_datafile$age_t1[l] <- temp1$Age[marker+1]
  new_datafile$pop_size_t1[l] <- temp1$pop_size[marker+1]
  new_datafile$R_pop_size_t1[l] <- temp1$R_pop_size[marker+1]
}

# save the new datafile
write.csv(new_datafile, "bio_data_t1.csv", row.names=F)

#### Import formatted data ####

datafile_t <- read.csv("bio_data_t.csv", header=T)
datafile_t1 <- read.csv("bio_data_t1.csv", header=T)
datafile_t <- cbind(datafile_t, datafile_t1)

#### Run development models ####

# Begin with a null model with only a random effect of year
null_model <- lmer(April_hatch_t1 ~ 1 + (1|Year), data=datafile_t)
# Look at the rests
summary(null_model)

# Then try the full model
full_model <- lmer(April_hatch_t1 ~ April_hatch + Synchrony + I(Synchrony^2) + 
                     Spring_temp_t1 + Spring_precip_t1 + 
                     as.factor(Age) + as.factor(Immigrant) + 
                     pop_size_t1 + as.factor(Round_t1) +
                     winter_temp + winter_precip_t + 
                     as.factor(beech) + (1|Year), data=datafile_t)
# Look at results
summary(full_model)

# Create dataframe for storing AIC results
AIC_output <- data.frame(AIC=rep(NA,5), delta_AIC=rep(NA,5), model_name=rep(NA,5))
AIC_output$AIC[1] <- AIC(null_model)
AIC_output$model_name[1] <- "null model"
AIC_output$AIC[2] <- AIC(full_model)
AIC_output$model_name[2] <- "full model"
AIC_output$delta_AIC[2] <- AIC_output[2,1]-AIC_output[1,1]
AIC_output

### MODEL SELECTION

# Want to remove any variables with estimate < 2*SE
# so looking for negative values of x
temp <- as.data.frame(summary(full_model)$coef)
temp$x <- abs(temp[,1]) - (2*(temp[,2]))
temp[which(temp$x < 0),]

## Continue until no variables have x < 0

# For our study the final model was:
red1_model <- update(full_model, . ~ . - as.factor(beech))
red2_model <- update(red1_model, . ~ . -winter_temp)
red3_model <- update(red2_model, . ~ . -winter_precip_t)
red4_model <- update(red3_model, . ~ . -as.factor(Immigrant))
red5_model <- update(red4_model, . ~ . -Spring_precip_t1)
AIC_output$AIC[3] <- AIC(red5_model)
AIC_output$model_name[3] <- "final model"
AIC_output$delta_AIC[3] <- AIC_output$AIC[3] - AIC_output$AIC[2]
AIC_output

save(red5_model, file="All_data_development.RData")