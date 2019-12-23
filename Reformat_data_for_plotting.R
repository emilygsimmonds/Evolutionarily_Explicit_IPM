#### Function to reformat data for plotting AND unscaling data

# Inputs: data_list = all output datafiles from cross validation
#         test_years = years predicted by model
#         descale = unscaled data

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

#### Function to descale data for plotting ####

# Inputs: x = scaled data (vector), descale = unscaled data frame, variable = variable name
unscale <- function(x, descale, variable){
  newx <- (x*sd(descale[,which(colnames(descale)==variable)]))+
    mean(descale[,which(colnames(descale)==variable)])
  return(newx)
}