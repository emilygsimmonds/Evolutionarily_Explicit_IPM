#### FIGURE 4 ####

#### Packages and function scripts ####
source('Reformat_data_for_plotting.R')

#### IMPORT AND FORMAT DATA ####
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled
cross_val <- read.csv("cross_validation_data.csv", header=T)

load("cross_val_RESULT_DIRECTIONAL.RData")
cross_val_RESULT_TP[[1]]$SYNC <- ((cross_val_RESULT_TP[[1]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
  ((cross_val_RESULT_TP[[1]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))

#### PLOT FIGURE ####

layout(matrix(c(1,2), 1,2, byrow=T), heights=c(4))

X <- 1:50

plot(PS ~ X, data = cross_val_RESULT_TP[[1]], type = 'l',
     xlab = "Years", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "black", lwd=2, cex.lab=1.5)
axis(1)
axis(2, las=1, cex.axis=1.25)

plot(SYNC ~ X, data = cross_val_RESULT_TP[[1]], type = 'l',
     xlab = "Years", ylab = "Synchrony", axes=F,
     ylim=c(-35,5), col = "green", lwd=2, cex.lab=1.5)
axis(1)
axis(2, las=1, cex.axis=1.25, at = seq(-35,5,5))
abline(h=-13, col="grey", lty = 2)
