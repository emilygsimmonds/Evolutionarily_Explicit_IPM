#### Script to run model selection for the inheritance function
# for the standard IPM

# COLLAPSE ALL FOR EASIER VIEWING

library(dplyr)

#### Import data inheritance data and add mother hatch date ####

# This dataset is restricted to only females that have an 
# identifed mother in the pedigree
# it also includes exact age in years of the individual
# not as a binary variable

datafile <- read.csv("bio_data_inheritance.csv", header=T)

# Create vector to store results, same length as datafile
mother_hatch <- rep(NA, length(datafile[,1]))
# Loop to extract mother ID and find the mother's hatch date 
for(j in 1:length(datafile[,1])){
  temp <- filter(datafile, F_ID == datafile$F_ID[j])
  marker <- which(as.character(datafile$F_ID) 
                  == as.character(temp$MOTHER[1]))
  ifelse(length(marker) == 0, mother_hatch[j] <- NA,
         mother_hatch[j] <- mean(datafile$April_hatch[marker]))
}

# Add to datafile
datafile$Mother_hatch <- mother_hatch
# Remove any NAs
datafile2 <- datafile[-which(is.na(datafile$Mother_hatch)),]
# Remove any immigrants
datafile2 <- filter(datafile2, Immigrant == 0)

#### Run inheritance models ####

## Begin with the null model
null <- lmer(April_hatch ~ 
               (1|Year), data=datafile2)

## Then the full model
full_model <- lmer(April_hatch ~ Mother_hatch + pop_size + Spring_temp +
                winter_temp + spring_precip_t + winter_precip_t + Round +
                (1|Year), data=datafile2)

# Look at the results
summary(full_model)

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
final_model <- lmer(April_hatch ~ Mother_hatch + pop_size  + Spring_temp +
                 Round +
                 (1|Year), data=datafile2)

output_AIC[2,1] <- AIC(final_model)
output_AIC[2,1] <- output_AIC[2,1]-output_AIC[1,1]
output_AIC[2,3] <- "final model"

save(final_model, file="All_data_inheritance_STANDARD.RData")