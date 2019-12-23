#### FIGURE 3 ####

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

plot(1:50,seq(1,75,length.out = 50), type='n', xlab="", ylab="", ylim=c(5,60), axes=F, cex.lab=1.5)
axis(2, las=1, at=c(5,30),labels=rep("",2), tck=0.05)
axis(2, las=1, at=c(35,60),labels=rep("",2), tck=0.05)
axis(2, las=1, at=seq(10,30,10), labels=seq(10,30,10), cex.axis=1.5)
axis(2, las=1, at=seq(40,60,10), labels=seq(10,30,10), cex.axis=1.5)
axis(1, at=seq(1,50,5), labels=seq(1,50,5), cex.axis=0.75)
ds_E.store <- (cross_val_RESULT_TP2$E_mean*sd(descale$April_hatch))+(mean(descale$April_hatch)/2)
lines(1:50, ds_E.store[1:50]+30, type='l', col=alpha("purple3",1))
ds_G.store <- (cross_val_RESULT_TP2$G_mean*sd(descale$April_hatch))+(mean(descale$April_hatch)/2)
lines(1:50, ds_G.store[1:50], type='l', col=alpha("darkorange",1))
abline(h=20, col="grey70", lty=2)
abline(h=(20+30), col="grey70", lty=2)
title(xlab="Years", line=2.5, cex.lab=1.5)
mtext("G and E", side=2, line=3, cex=1.5)