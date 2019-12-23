#### Script containing wrapper functions for the survival analysis
## These are needed for K-fold cross validation

# COLLAPSE ALL FOR EASIER VIEWING

library(BaSTA)

#### create marked data ####

# Inputs: datafile = full biological data (scaled),
#         basta_input = capture history
#         full_years = number of years in full dataset
#         test_years = years on which to test cross validation
#         type = what kind of analysis? K_fold or Full data

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
  
  # for each individual add their Immigrant status and section of woodland
  for(l in IDs){
    temp_data <- subset(datafile, datafile$F_ID == l)
    marker <- which(row.names(temp)==l)
    ifelse(temp_data$Immigrant[1] == 1, temp$Immigrant[marker] <- "Immigrant", 
           temp$Immigrant[marker] <- "Resident")
    temp$Round[marker] <- as.character(temp_data$Round[which.max(temp_data$Round)])
  }
  
  # now the individual covariates
  if(type=="K_fold"){marker <- seq(which(full_years == min(test_years)-1), 
                                   which(full_years == max(test_years)-1), 1)
  years <- full_years[-marker]}else{years <- full_years}
  output_store <- matrix(NA, ncol=length(years)*5, nrow=length(temp[,1]))
  
  # this loop should create an output that has one row for each individual and columns for 
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
                              April_hatch2 = "ahs", Synchrony2 = "sns", 
                              Year = seq(1960, 2004, 1)) 
  if(type=="Full"){col_name_data <- data.frame(April_hatch = "ah", Synchrony = "sn", 
                                               Clutch_size = "cs",
                                               April_hatch2 = "ahs", Synchrony2 = "sns",
                                               Year = seq(1960, 2009, 1)) }
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

# Inputs: processed_data = output from data processing, 
#         full_years = number of years in full dataset,
#         test_years = years on which to test cross validation,
#         year_variables = dataframe of all year varying variables,
#         parameters = design parameters for Phi and p

create_design <- function(processed_data, test_years, full_years, year_variables, parameters){
  ddl <- make.design.data(processed_data, parameters = parameters)
  marker <- seq(which(full_years == min(test_years)-1), which(full_years == max(test_years)-1), 1)
  years <- full_years[-marker]
  # Now have dataframe organised by year and individual. Here is where to add the year varying variables.
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
