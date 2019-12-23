#### Script to run model selection for the survival function

# COLLAPSE ALL FOR EASIER VIEWING

# The recapture rate for our population is approx 81 % 
# so parameterisation of the survival function uses a
# capture-mark-recapture model

#### Import data scaled data and create capture history ####

datafile <- read.csv("bio_data.csv", header=T)

## Load BaSTA package to create a capture history from survival data
library(BaSTA)
## Use ID and Year columns to indicate survival and year as a time step
CaptHist <- CensusToCaptHist(datafile$F_ID, datafile$Year, 
                                dformat = "%Y", timeInt = "Y")

# Determine how long the capture history is and collapse to a single column
mark <- length(CaptHist[,1])
CaptHist$capture_history <- apply(CaptHist[,2:mark] , 1 , 
                                  paste , collapse = "" )

# Set up all factors, which are individually varying but not 
# time varying. We have only two: Immigrant status and Section of the woodland (Round)
CaptHist$Immigrant <- rep(NA, length(CaptHist[,1]))
CaptHist$Round <- rep(NA, length(CaptHist[,1]))

# Set up a vector of the unique female IDs
IDs <- row.names(CaptHist)

# Loop to assign values for the individually varying covariates
for(l in IDs){
  # Reduce input data to one individual only
  temp_data <- subset(datafile, datafile$F_ID == l)
  # Assign a marker for this individual in the capture history dataset
  marker <- which(row.names(CaptHist)==l)
  # If the input data shows an individual is an Immigrant (1) give that
  # value, otherwise, assign value 'Resident'
  ifelse(temp_data$Immigrant[1] == 1, CaptHist$Immigrant[marker] <- "Immigrant", 
         CaptHist$Immigrant[marker] <- "Resident")
  # Take the section of the woodland from the input data for this
  # individual and assign to capture history
  CaptHist$Round[marker] <- as.character(temp_data$Round[which.max(temp_data$Round)])
}

# Next we add covariates which vary by individual and by time.
# Here we have 5: hatch date (April_hatch), 
# synchrony, clutch size,
# hatch date quadratic, synchrony quadratic

# Look at the column names in the input data to find correct columns 
# for the covariates
colnames(datafile) # 6,7,8, (6^2), (7^2) 

# Create a vector of all unique years in the data
years <- unique(datafile$Year)

# Create a matrix with ncol = number of years* number of covariates
# and nrow = number of individuals
output_store <- matrix(NA, ncol=50*5, nrow=length(CaptHist[,1]))

# Loop to create output that has one row for each individual
# and columns for every value for each variable for 50 years
for(l in 1:50){
  
  # first subset to only the year required
  year_data <- subset(datafile, datafile$Year == years[l]) 
  
  # Nested loop for each ID
  for(m in row.names(CaptHist)){ 
    
    # find when the first ID appears in year data
    marker <- which(as.character(year_data$F_ID) == m) 
    
    # generates marker in each dataset that references the correct ID
    marker2 <- which(row.names(CaptHist) == m) 
    
    # Now want for each variable to fill in the first year for each individual:
    # so var 1 = 1 var 2 = 36+1... etc...
    
    ifelse(length(marker) > 0, 
           output_store[marker2,l] <- year_data[marker,6], 
           output_store[marker2,l] <- mean(year_data[,6], na.rm=T))
    ifelse(length(marker) > 0, 
           output_store[marker2,l+50] <- year_data[marker,7], 
           output_store[marker2,l+50] <- mean(year_data[,7], na.rm=T))
    ifelse(length(marker) > 0,
           output_store[marker2,l+100] <- year_data[marker,8], 
           output_store[marker2,l+100] <- mean(year_data[,8], na.rm=T))
    ifelse(length(marker) > 0, 
           output_store[marker2,l+150] <- (year_data[marker,6]^2), 
           output_store[marker2,l+150] <- (mean(year_data[,6], na.rm=T)^2))
    ifelse(length(marker) > 0, 
           output_store[marker2,l+200] <- (year_data[marker,7]^2), 
           output_store[marker2,l+200] <- (mean(year_data[,7], na.rm=T)^2))                   
  }
  
}

# Create a vector of column names for the individual, time varying covariates
# need to name them in specific way: hatch date in 1960 = ah1960 ... etc...
col_name_data <- data.frame(April_hatch = "ah", 
                            Synchrony = "sn", 
                            Clutch_size = "cs", 
                            April_hatch2 = "ahs", 
                            Synchrony2 = "sns", 
                            Year = seq(1960, 2009, 1)) 

# Set up markers for which column names to use where
marker1 <- c(1,6)
marker2 <- c(2,6)
marker3 <- c(3,6)
marker4 <- c(4,6)

column_names <- c(apply(col_name_data[,marker1] , 1 , paste , collapse = "" ), 
                  apply(col_name_data[,marker2] , 1 , paste , collapse = "" ), 
                  apply(col_name_data[,marker3] , 1 , paste , collapse = "" ), 
                  apply(col_name_data[,marker4] , 1 , paste , collapse = "" ), 
                  apply(col_name_data[,5:6] , 1 , paste , collapse = "" ))

# Save the column names
write.csv(column_names, "column_names_indiv_time_varying.csv", row.names = F)

# Give the column names to the output matrix
colnames(output_store) <- column_names

# Combine all components of the dataset (capture history, individual covariates, and
# individual time varying covariates, the end column)

# Column must end in 1;
end <- rep("1;", mark)

# This will be the input to the marked and RMark packages
marked_input <- cbind(as.character(CaptHist$capture_history), 
                      CaptHist$Immigrant, 
                      CaptHist$Round, output_store, end)

# Rename the columns, this is needed for importing
colnames(marked_input) <- c("ch", "Immigrant", "Round", column_names, "end")
marked_input <- as.data.frame(marked_input)
marked_input$ch <- as.character(marked_input$ch)

# Save it 
write.table(marked_input, "capture_history.txt", row.names=F)

#### Import data and begin capture-mark-recapture dataset ####

# Import dataframe of variables that vary by year but not by individual
# These are: temperature, rainfall, beech mast index, and population size
year_variables <- read.csv("year_variables.csv", header=T)

# Load the RMark package to import out capture history
# this was made with the steps above
library(RMark)

# Importing three factors ("f") and 250 numeric variables ("n")
capture_history <- import.chdata(filename = "capture_history.txt", 
                                 field.types=c("f", "f", rep("n", 250), "f"), 
                                 header=T)

# detach package as different functions needed below
detach("package:RMark", unload=TRUE)

#### Run survival model ####

# Load marked package to run the model
library(marked)

# Create processed dataframe model types is CJS,
# factors indicated using groups argument,
# begin.time = first year of study,
# nocc = number of occasions i.e. number of years in study
data_processed <- process.data(capture_history, model="CJS", 
                               begin.time = 1960,
                               groups = c("Immigrant", "Round"), 
                               nocc = 50)

# Create design data with static and time-varying covariates
# Static are the factors (individual covariates)
# Time varying individual covariates must be named as their column names
# cs = clutch size, sn = synchrony, sns = quadratic synchrony,
# ah = hatch date, ahs = quadratic hatch date

# Phi = survival rate, p = recapture rate
design.Phi = list(static = c("Immigrant", "Round"), time.varying = c("cs", "sn", "ah", "sns", "ahs"))
design.p = list(time.varying=c("ah"))

design.parameters <- list(Phi = design.Phi, p = design.p)

# This line names the design data, called ddl
ddl <- make.design.data(data_processed, parameters = design.parameters)

# Now we want to add our year varying variables to Phi and p
# Phi and p are slightly out of sink, p only starts in year 2 of our data
# So we have to set up a time variable for our year variables to reflect this
year_variables$time <- 1960:2009
# Then we merge the year data with our design data
ddl$Phi <- merge_design.covariates(ddl$Phi, year_variables, bytime=TRUE)
ddl$Phi$Beech <- as.factor(ddl$Phi$Beech)
year_variables2 <- year_variables
year_variables2$time[1:50] <- 1961:2010
ddl$p <- merge_design.covariates(ddl$p, year_variables2, bytime=TRUE)

# Got design data so now try a model

## Begin with a null model with no variables in Phi or p
Phi.1 = list(formula=~1) 
p.1 = list(formula=~1)

null_model <- crm(data_processed, ddl, 
                  hessian = T, model="CJS", 
                  model.parameters=list(Phi=Phi.1, p=p.1), 
                  accumulate = F) # don't want to accumulate histories
# Take a look at the results
summary(null_model$results)
null_model$results$beta$Phi

## Then try the full model - this will take much longer
Phi.full = list(formula=~Immigrant+Round+
                  ah+ahs+sns+sn+cs+
                  Pop_size+Beech+Spring_precip+Winter_precip+Spring_temp+Winter_temp)
p.full = list(formula=~time)

full_model <- crm(data_processed, ddl, 
                  hessian = T, model="CJS", 
                  model.parameters=list(Phi=Phi.full, p=p.full), 
                  accumulate = F) # don't want to accumulate histories
# Take a look at the results
summary(full_model$results)
full_model$results$beta$Phi

# Store the AIC result
AIC_store[1,1] <- as.numeric(full_model$results$AIC)
AIC_store[1,3] <- "full model"

### MODEL SELECTION

# Extract estimates that are less than 2 times their standard error
# Remove from the next model the one with the smallest estimate to standard error ratio
row.names(coef(full_model))[which((abs(coef(full_model)$Estimate) - 2*coef(full_model)$se) < 0)]

# Re-run the new model until you have no variables with estimates < 2*SE

# For our study we get this final model
Phi.5 = list(formula=~Immigrant+Round+
               ah+sns+sn+
               Pop_size+Beech+Spring_precip+Spring_temp+Winter_temp) 

red5_model <-  crm(data_processed, ddl, hessian = T, model="CJS", 
                   model.parameters=list(Phi=Phi.5, p=p.full), accumulate = F)
coef(red5_model)

AIC_store[2,1] <- as.numeric(red5_model$results$AIC)
AIC_store[2,2] <- AIC_store[2,1]-AIC_store[1,1]
AIC_store[2,3] <- "final model"

save(red5_model, file="All_data_survival.RData")
