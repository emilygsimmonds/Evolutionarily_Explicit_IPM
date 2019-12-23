#### Script to run model selection for the inheritance function
# for the quantitative genetic IPM

# COLLAPSE ALL FOR EASIER VIEWING

library(dplyr)
library(MCMCglmm) 

#### Import data inheritance data and load pedigree ####

# This dataset is restricted to only females that have an 
# identifed mother in the pedigree
# it also includes exact age in years of the individual
# not as a binary variable

datafile <- read.csv("bio_data_inheritance.csv", header=T)

ped <- read.csv("ped_data.csv", header=T)

# Remove any individuals missing from the pedigree
F_IDs <- unique(datafile$F_ID)
missing <- rep(NA, length(F_IDs))
F_IDs <- toupper(F_IDs)
ped$ID<- toupper(ped$ID)
ped$MOTHER <- toupper(ped$MOTHER)
datafile$F_ID <- toupper(datafile$F_ID)

for(j in 1:length(F_IDs)){
  
  marker <- which(ped$ID == as.character(F_IDs[j]))
  if(length(marker) > 0){missing[j] <- "F"}
  marker2 <- which(datafile$F_ID == as.character(F_IDs[j]))
  if(length(marker) == 0){datafile <- datafile[-marker2,]}
}

# need to have an animal and ID and a BYEAR column.
str(datafile) # check structure incase need to change some variables
# Change age back to a binary variable
datafile$Age[which(datafile$Age > 1)] <- 2
# Ensure all factors have correct structure
datafile$Age <- as.factor(datafile$Age)
datafile$beech <- as.factor(datafile$beech)
datafile$Immigrant <- as.factor(datafile$Immigrant)
datafile$animal <- as.factor(datafile$animal)
str(datafile) # Check again

#### Run inheritance models ####

# Need to set up priors based on previous assessment of heritability from van der Jeugd and McCleery 2002 (0.16)
# so prior for additive genetic should be 0.16*phenotypic variance

# priors, setting up one for each random effect:
# as including Va (G1), permanent environment (G2), mother ID (G3) and residual (R)
p.var <- var(datafile$April_hatch, na.rm=T) # variance of phenotype # 0.16 split into Va
prior.rand <- list(G=list(G1=list(V=matrix(p.var*0.16), n=1), 
                            G2=list(V=matrix(p.var*(0.84/3)), n=1),
                            G3=list(V=matrix(p.var*(0.84/3)), n=1)),
                     R=list(V=matrix(p.var*(0.84/3)), n=1))  


## Begin with a null model
null_model <- MCMCglmm(April_hatch ~ 1, random=~animal+ID+MOTHER, pedigree=ped, 
                            data=datafile, prior=prior.rand, verbose=F, 
                            nitt=50000, burnin=5000, thin=50) 
# using high iteration number and burn in to try and avoid autocorrelation
summary(null_model) # look at result

### THERE WILL BE WARNINGS BUT IT IS OK

## Then run the full model
full_model <- MCMCglmm(April_hatch ~ Spring_temp + 
                         winter_temp + winter_precip_t + 
                         spring_precip_t + 
                         Round + pop_size, 
                         random=~animal+ID+MOTHER, pedigree=ped, 
                  data=datafile, prior=prior.rand, verbose=F, 
                  nitt=50000, burnin=5000, thin=50)

### MODEL SELECTION

# Create output matrix to store DIC
DIC_output <- matrix(NA, ncol=3, nrow=10)
colnames(DIC_output) <- c("DIC", "DIC diff", "variable removed")
DIC_output[1,1] <- summary(full_model)$DIC
DIC_output[1,3] <- "Full model"

# Check which variables have a p-value >0.05
table.coef <- summary(full_model)$solutions
non_sig <- table.coef[which(table.coef[,5] > 0.05),]

# Keep removing variables until none have p-value >0.05
# In our study the final model was:

final_model <- MCMCglmm(April_hatch ~ Spring_temp + winter_temp + winter_precip_t + 
                          spring_precip_t + 
                          Round + pop_size, random=~animal+ID+MOTHER, pedigree=ped, 
                        data=datafile_I, prior=prior.rand, verbose=F, 
                        nitt=50000, burnin=5000, thin=50)

DIC_output[1,1] <- summary(final_model)$DIC
DIC_output[1,3] <- "Final model"

save(final_model, file="All_data_inheritance_QG.RData")
