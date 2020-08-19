# !diagnostics off

# FIGURES - IPM CODE
# code to run all figures
######################################################################


#### Packages and Data ####
datafile <- read.csv("updated_model_input_11.10.18.csv", header=T) # now scaled correctly
head(datafile)
summary(datafile)
descale <- read.csv("updated_descale_11.10.18.csv", header=T)
# Development datasets 
datafile_t <- read.csv("datafile_t_11.10.18.csv", header=T)
datafile_t1 <- read.csv("datafile_t1_11.10.18.csv", header=T)
datafile_t <- cbind(datafile_t, datafile_t1)
datafile_t$Year <- datafile_t$Year + 1
year_variables <- read.csv("year_variables_Jan19.csv", header=T)
# Inheritance dataset
datafile_I <- read.csv("datafile_I_11.10.18.csv", header=T)
# need to have an animal and ID and a BYEAR column.
datafile_I$Age[which(datafile_I$Age > 1)] <- 2
datafile_I$Age <- as.factor(datafile_I$Age)
datafile_I$beech <- as.factor(datafile_I$beech)
datafile_I$Immigrant <- as.factor(datafile_I$Immigrant)
datafile_I$animal <- as.factor(datafile_I$animal)
ped <- read.csv("ped_inheritance_11.10.16.csv", header=T)
cross_val <- read.csv("cross_validation_data_1961to2010_19.10.csv", header=T)

# K-fold cross validation years
list_K_fold <- list(1961:1965, 1966:1970, 1971:1975, 1976:1980,
                    1981:1985, 1986:1990, 1991:1995, 1996:2000,
                    2001:2005, 2006:2010)

library(lme4)
library(MCMCglmm)
library(dplyr)
library(mvtnorm)
library(scales)
library(ggplot2)
library(viridis)
library(data.table)
library(egg)
library(grid)
library(gridExtra)
library(cowplot)

source("Functions_IPM1.R")


# Code to draw all figures needed for chapter 4

# all linked to correct directories and datafiles
# includes supplementary


### FIGURE 1 IPM-1 ####
# cross validation

cross_val <- read.csv("cross_validation_data_1961to2010_19.10.csv", header=T)
load("cross_val_RESULT_RO_Nov19.RData")
cross_val_RESULT_RO1 <- list(cross_val_RESULT_RO[[1]][[1]])
for(i in 2:10){cross_val_RESULT_RO1 <- 
  c(cross_val_RESULT_RO1, list(cross_val_RESULT_RO[[i]][[1]]))}
cross_val_RESULT_RO2 <- for_plotting_CV(cross_val_RESULT_RO1, list_K_fold, descale)

load("cross_val_RESULT_TP_Nov19_N_double.RData")
cross_val_RESULT_TP1 <- list(cross_val_RESULT_TP[[1]][[1]])
for(i in 2:10){cross_val_RESULT_TP1 <- 
  c(cross_val_RESULT_TP1, list(cross_val_RESULT_TP[[i]][[1]]))}
cross_val_RESULT_TP2 <- for_plotting_CV(cross_val_RESULT_TP1, list_K_fold, descale)

load("cross_val_RESULT_po_Nov19_N_double.RData")
cross_val_RESULT_PO1 <- list(cross_val_RESULT_po[[1]][[1]])
for(i in 2:10){cross_val_RESULT_PO1 <- 
  c(cross_val_RESULT_PO1, list(cross_val_RESULT_po[[i]][[1]]))}
cross_val_RESULT_PO2 <- for_plotting_CV(cross_val_RESULT_PO1, list_K_fold, descale)

png(filename="Figure_1_IPM1.png", units="cm", height = 20, width = 16,
    res=400)
layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))

# two population size first
par(mar=c(2,4.5,2,1))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, cex.lab=1.5,
     main="Quantitative genetic IPM")
title(main='a.', adj=0, cex.main=1)
#axis(1)
axis(2, las=1, cex.axis=1.25)
points(PS ~ Year, data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main="Standard IPM")
title(main='c.', adj=0, cex.main=1)
#axis(1)
#axis(2, las=1)
points(PS ~ Year, data = cross_val_RESULT_PO2,
       type = 'l', col = 'black', lty = 4)

# then hatch dates
par(mar=c(4,4.5,2,1))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "Mean hatch date", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2, cex.lab=1.5)
title(main='b.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
axis(2, las=1, cex.axis=1.25, at = seq(30,60,10))
points(Hatch_mean ~ Year, data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(4,1,2,4.5))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "Mean hatch date", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2, cex.lab=1.5)
title(main='d.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
#axis(2, las=1, cex.axis = 1.25)
points(Hatch_mean ~ Year, data = cross_val_RESULT_PO2,
       type = 'l', col = 'black', lty = 4)

plot.new()
legend("top", c("Model projections", "Observations"), 
       bty="n", col = c("black",  "grey"), lty=c(1,1),
       cex=1.5, ncol=2, lwd=c(2,2))

dev.off()



# calculate mean squared error

MSE <- data.frame(#ps.RO = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_RO2$PS)),2),
  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP2$PS)),2),
  ps.po = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_PO2$PS)),2),
  #z.RO = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_RO2$Hatch_mean)),2),
  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP2$Hatch_mean)),2),
  z.po = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_PO2$Hatch_mean)),2)
)

colnames(MSE) <- rep(c( "Quantitative genetic IPM", "Standard IPM"), 2)
MSE
#write.csv(MSE, "MSE_results_Nov19.csv", row.names=F)

### FIGURE 2 IPM-1 ####
# cross validation - 50 years

load("cross_val_RESULT_TP_50_Nov19.RData")
cross_val_RESULT_TP2 <- cross_val_RESULT_TP[[1]]
cross_val_RESULT_TP2$Year <- 1961:2010
cross_val_RESULT_TP2$Hatch_mean <- (cross_val_RESULT_TP2$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)

load("cross_val_RESULT_po_50_Nov19.RData")
cross_val_RESULT_PO2 <- cross_val_RESULT_po[[1]]
cross_val_RESULT_PO2$Year <- 1961:2010
cross_val_RESULT_PO2$Hatch_mean <- (cross_val_RESULT_PO2$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)


png(filename="~/Dropbox/Working PhD chapter drafts/IPM/Finals to submit/Oikos/Revisions/finals to resubmit/Next revisions - 11.19/To Submit Dec 2019/Figure2.png",
    units="cm", height=20, width=16, res=400)
layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))

par(mar=c(2,4.5,2,1))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, cex.lab=1.5,
     main = "Quantitative genetic IPM")
title(main='a.', adj=0, cex.main=1)
#axis(1)
axis(2, las=1, cex.axis=1.25)
points(PS[2:51] ~ Year[1:50], data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2,
     main = "Standard IPM")
title(main='c.', adj=0, cex.main=1)
#axis(1)
#axis(2, las=1)
points(PS[2:51] ~ Year[1:50], data = cross_val_RESULT_PO2,
       type = 'l', col = 'black', lty = 4)

par(mar=c(4,4.5,2,1))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "Mean hatch date", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2, cex.lab=1.5)
title(main='b.', adj=0, cex.main=1)
axis(1)
axis(2, las=1, cex.axis=1.25, at = seq(30,60,10))
points(Hatch_mean[2:51] ~ Year[1:50], data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(4,1,2,4.5))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2, cex.lab=1.5)
title(main='d.', adj=0, cex.main=1)
axis(1)
#axis(2, las=1)
points(Hatch_mean[2:51] ~ Year[1:50], data = cross_val_RESULT_PO2,
       type = 'l', col = 'black', lty = 4)

plot.new()
legend("center", c("Model projections", "Observations"), 
       bty="n", col = c("black", "grey"), lty=c(1,1),
       cex=2, ncol=2, lwd=c(2,2))

dev.off()



# calculate mean squared error

MSE <- data.frame(#ps.RO = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_RO2$PS)),2),
  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP2$PS)),2),
  ps.po = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_PO2$PS)),2),
  #z.RO = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_RO2$Hatch_mean)),2),
  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP2$Hatch_mean)),2),
  z.po = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_PO2$Hatch_mean)),2)
)

colnames(MSE) <- rep(c( "Quantitative genetic IPM", "Standard IPM"), 2)
MSE
#write.csv(MSE, "MSE_results_50_Nov19.csv", row.names=F)

# check sign
pop1 <- cross_val_RESULT_TP2$PS[1:49]
pop2 <- cross_val_RESULT_TP2$PS[2:50]

change_pred <- pop2 - pop1
change_obs <- cross_val$R_only_PS[2:50] - cross_val$R_only_PS[1:49]

plot(change_pred, change_obs, pch=16, ylab="Observed population change", xlab = "Projected population change")
abline(h=0, v=0)
polygon(x=c(-100,0,0,-100), y=c(0,0,150,150), col=2)
polygon(x=c(0,100,100,0), y=c(0,0,-150,-150), col=2)
points(change_pred, change_obs, pch=16)

### FIGURE 3 IPM-1 ####
# Plot of G and E
load("cross_val_RESULT_TP_50_DIRECTIONAL.RData")
cross_val_RESULT_TP2 <- cross_val_RESULT_TP



png(filename="Figure_GE.png", units="cm", height = 20, width = 16,
    res=400)
#layout(matrix(c(1,2), 2,1, byrow=T), heights=c(4,2))

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

dev.off()


# calculate slopes
summary(lm(ds_E.store ~ seq(1,50,1)))
summary(lm(ds_G.store ~ seq(1,50,1)))

### FIGURE 4 IPM-1 ####
# Directional change
load("cross_val_RESULT_TP_50_DIRECTIONAL.RData")
cross_val_RESULT_TP$SYNC <- ((cross_val_RESULT_TP$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
  ((cross_val_RESULT_TP$HF*sd(descale$Half_fall))+mean(descale$Half_fall))

png(filename="Figure_directional.png", units="cm", height = 16, width = 16,
    res=400)
layout(matrix(c(1,2), 1,2, byrow=T), heights=c(4))

X <- 1:50

plot(PS ~ X, data = cross_val_RESULT_TP, type = 'l',
     xlab = "Years", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "black", lwd=2, cex.lab=1.5)
axis(1)
axis(2, las=1, cex.axis=1.25)

plot(SYNC ~ X, data = cross_val_RESULT_TP, type = 'l',
     xlab = "Years", ylab = "Synchrony", axes=F,
     ylim=c(-35,5), col = "green", lwd=2, cex.lab=1.5)
axis(1)
axis(2, las=1, cex.axis=1.25, at = seq(-35,5,5))
abline(h=-13, col="grey", lty = 2)

dev.off()


#### SOM FIGURE 1 IPM-1 ####

# cross validation

cross_val <- read.csv("cross_validation_data_1961to2010_19.10.csv", header=T)
load("cross_val_RESULT_RO_Nov19.RData")
cross_val_RESULT_RO1 <- list(cross_val_RESULT_RO[[1]][[1]])
for(i in 2:10){cross_val_RESULT_RO1 <- 
  c(cross_val_RESULT_RO1, list(cross_val_RESULT_RO[[i]][[1]]))}
cross_val_RESULT_RO3 <- for_plotting_CV(cross_val_RESULT_RO1, list_K_fold, descale)

load("cross_val_RESULT_TP_Nov19.RData")
cross_val_RESULT_TP1 <- list(cross_val_RESULT_TP[[1]][[1]])
for(i in 2:10){cross_val_RESULT_TP1 <- 
  c(cross_val_RESULT_TP1, list(cross_val_RESULT_TP[[i]][[1]]))}
cross_val_RESULT_TP3 <- for_plotting_CV(cross_val_RESULT_TP1, list_K_fold, descale)

png(filename="IPM1_SOM_figure.png", units="cm", height = 20, width = 16,
    res=400)
layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))

# two pop sizes first with no X axis
par(mar=c(2,4.5,2,1))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main = "Resident only",
     cex.lab=1.5)
title(main='a.', adj=0, cex.main=1)
#axis(1)
axis(2, las=1, cex.axis=1.25)
points(PS ~ Year, data = cross_val_RESULT_RO3,
       type = 'l', col = 'black', lty = 3)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main = "Total population")
title(main='c.', adj=0, cex.main=1)
#axis(1)
#axis(2, las=1)
points(PS ~ Year, data = cross_val_RESULT_TP3,
       type = 'l', col = 'black', lty = 1)


# two mean hatch with X axes
par(mar=c(4,4.5,1,1))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "Mean hatch date", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2,
     cex.lab=1.5)
title(main='b.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
axis(2, las=1, at = seq(30,60,10), cex.axis=1.25)
points(Hatch_mean ~ Year, data = cross_val_RESULT_RO3,
       type = 'l', col = 'black', lty = 3)
par(mar=c(4,1,1,4.5))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2,
     cex.lab=1.5)
title(main='d.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
#axis(2, las=1)
points(Hatch_mean ~ Year, data = cross_val_RESULT_TP3,
       type = 'l', col = 'black', lty = 1)


plot.new()
legend("top", c("Model projections", "Observations"), 
       bty="n", col = c("black", "grey"), lty=c(1,1),
       cex=1.5, ncol=2, lwd=c(2,2))

dev.off()



# calculate mean squared error
load("cross_val_RESULT_po_Sept19.RData")
cross_val_RESULT_PO2 <- for_plotting_CV(cross_val_RESULT_po, list_K_fold, descale)

MSE <- data.frame(ps.RO = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_RO3$PS)),2),
                  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP3$PS)),2),
                  #ps.po = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_PO2$PS)),2),
                  z.RO = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_RO3$Hatch_mean)),2),
                  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP3$Hatch_mean)),2)
                  #z.po = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_PO2$Hatch_mean)),2)
)

colnames(MSE) <- rep(c("Resident only", "Total population"), 2)
MSE
write.csv(MSE, "MSE_results_SOM_Nov19.csv", row.names=F)

#### SOM FIGURE 2 IPM-1 ####

load("cross_val_RESULT_TP_Nov19.RData")
cross_val_RESULT_TP3 <- for_plotting_CV(cross_val_RESULT_TP, list_K_fold, descale)
load("cross_val_RESULT_TP_Nov19_N_double.RData")
cross_val_RESULT_TP2 <- for_plotting_CV(cross_val_RESULT_TP, list_K_fold, descale)

png(filename="IPM1_SOM_figure_correction.png", units="cm", height = 20, width = 16,
    res=400)
layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))

# two pop sizes first with no X axis
par(mar=c(2,4.5,2,1))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main = "Total population")
title(main='a.', adj=0, cex.main=1)
#axis(1)
axis(2, las=1, cex.axis=1.25)
points(PS ~ Year, data = cross_val_RESULT_TP3,
       type = 'l', col = 'black', lty = 1)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main = "Total population \nwith correction")
title(main='c.', adj=0, cex.main=1)
#axis(1)
#axis(2, las=1)
points(PS ~ Year, data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)


# two mean hatch with X axes
par(mar=c(4,4.5,1,1))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "Mean hatch date", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2,
     cex.lab=1.5)
title(main='b.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
axis(2, las=1, at = seq(30,60,10), cex.axis=1.25)
points(Hatch_mean ~ Year, data = cross_val_RESULT_TP3,
       type = 'l', col = 'black', lty = 1)
par(mar=c(4,1,1,4.5))
plot(R_Mean_hatch ~ Year, data = cross_val, type = 'l',
     xlab = "Time (years)", ylab = "", axes=F,
     ylim=c(30,60), col = "grey70", lwd=2,
     cex.lab=1.5)
title(main='d.', adj=0, cex.main=1)
axis(1, cex.axis=1.25)
#axis(2, las=1)
points(Hatch_mean ~ Year, data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)


plot.new()
legend("top", c("Model projections", "Observations"), 
       bty="n", col = c("black", "grey"), lty=c(1,1),
       cex=1.5, ncol=2, lwd=c(2,2))

dev.off()



# calculate mean squared error

MSE <- data.frame(#ps.RO = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_RO2$PS)),2),
  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP3$PS)),2),
  ps.TPC = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP2$PS)),2),
  #z.RO = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_RO2$Hatch_mean)),2),
  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP3$Hatch_mean)),2),
  z.TPC = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP2$Hatch_mean)),2)
)

colnames(MSE) <- rep(c("Total population", "Corrected"), 2)
MSE
write.csv(MSE, "MSE_results_SOM_C_Nov19.csv", row.names=F)

#### SOM FIGURE 3 IPM-1 ####

residents <- group_by(datafile, Immigrant, Year) %>% summarise(n=n())
proportion <- residents$n[1:50]/(residents$n[1:50]+residents$n[51:100])

png(filename="IPM1_SOM_figure_prop_resident.png", units="cm", height = 16, width = 16,
    res=400)

plot(x=1961:2010, proportion, type='l', axes=F, 
     ylab = "Proportion of residents", xlab="Time (years)",
     ylim=c(0,1))
axis(1)
axis(2,las=1)

dev.off()


#### SOM FIGURE 4 IPM-1 ####

load("cross_val_RESULT_TP_Sept19_N_double.RData")
cross_val_RESULT_TP2 <- for_plotting_CV(cross_val_RESULT_TP, list_K_fold, descale)
load("cross_val_RESULT_po_Feb19_N_double.RData")
cross_val_RESULT_PO2 <- for_plotting_CV(cross_val_RESULT_po, list_K_fold, descale)

png(filename="Figure_Variance.png", units="cm", height = 20, width = 16,
    res=400)
layout(matrix(c(1,2), 2,1, byrow=T), heights=c(4,2))

plot(cross_val_RESULT_PO2$Year, sqrt(cross_val_RESULT_PO2$Hatch_var),
     type = 'l', lty = 4, axes = F, xlab = "Year", ylab = "Standard deviation in hatch date",
     ylim = c(0,10))
axis(1)
axis(2, las=1)
lines(cross_val_RESULT_TP2$Year, sqrt(cross_val_RESULT_TP2$Hatch_var))

plot.new()
legend("center", c("Standard IPM", "Quantitative genetic IPM"), 
       bty="n", col = c("black",  "black"), lty=c(4,1),
       cex=1, ncol=2, lwd=c(2,2))

dev.off()

### FUNCTION PLOTS SOM IPM-1 ####
### survival ####
load('final_survival_Feb19.RData')
surv_output <- coef(final_survival)
SEs <- data.frame(est = surv_output$Estimate[1:19],
                  ucl = surv_output$Estimate[1:19]+(2*surv_output$se[1:19]),
                  lcl = surv_output$Estimate[1:19]-(2*surv_output$se[1:19]))

# save out tslope values
slopes <- as.numeric(final_survival$results$beta$Phi[11:19])
slopes <- slopes[-5] # remove those that are factors
slopes <- slopes[-5] 
xs <- c(mean(datafile$April_hatch), mean(datafile$Synchrony),mean(datafile$Synchrony)^2,mean(datafile$pop_size), mean(datafile$spring_precip_t), mean(datafile$Spring_temp), mean(datafile$winter_temp))
sum_of_slope <- sum(slopes*xs) 
int <- as.numeric(final_survival$results$beta$Phi[1]) 
B.0 <- as.numeric(final_survival$results$beta$Phi[15]) 
B.1 <- as.numeric(final_survival$results$beta$Phi[16]) 
Resident <- as.numeric(final_survival$results$beta$Phi[10]) 
Rounds <- as.numeric(final_survival$results$beta$Phi[2:9]) 
# matrix of slope values to be used for predicting
# number of rows = number of variables *30
# number of columns = continuous variables
test_mat_xs <- matrix(NA, nrow=(30*10), ncol=7)
for(i in 1:300){test_mat_xs[i,] <- xs}


# For each variable want to calculate number of duplicate values
# Then want to predict values and plot
png(filename="SOM1_Figure1.png", units="cm", height = 20, width = 18,
    res=400)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=T), heights=rep(7.5,3), widths = rep(6,3))

par(mar=c(5, 4, 1, 1))

# Hatch
datafile <- datafile %>% group_by(Survival, April_hatch) %>% mutate(coder = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[1:30,1] <- seq(min(datafile$April_hatch), max(datafile$April_hatch), length.out = 30)
predictions_H <- rep(NA, 30)
predictions_Hu <- rep(NA, 30)
predictions_Hl <- rep(NA, 30)
for(i in 1:30){predictions_H[i] <- logistic(u=(int+sum(slopes*test_mat_xs[i,])+0))
predictions_Hu[i] <- logistic(u=(int+sum(c(SEs$ucl[11],slopes[2:7])*test_mat_xs[i,])+0))
predictions_Hl[i] <- logistic(int+sum(c(SEs$lcl[11],slopes[2:7])*test_mat_xs[i,])+0)
}
# for CIs need to order the CI predictions
predictions_CI <- cbind(c(test_mat_xs[1:30,1], seq(max(datafile$April_hatch), min(datafile$April_hatch), length.out = 30)), c(predictions_Hl, sort(predictions_Hu)))

plot(Survival ~ April_hatch, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Hatch date", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-3,3.5,0.5), labels=c(-3.0,"",-2.0,"",-1.0,"",0.0,"",1.0,"",2.0,"",3.0,""))
axis(2, las=1, cex.axis=1)
title(ylab="Survival", line=2.5, cex.lab=1.5)
points(datafile$April_hatch, datafile$Survival, cex=datafile$coder/50, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[1:30,1], predictions_H, col="red", lwd=2)

# Synchrony^2
datafile <- datafile %>% group_by(Survival, Synchrony) %>% mutate(coder2 = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[31:60,2] <- (seq(min(datafile$Synchrony), max(datafile$Synchrony), length.out = 30))^2
test_mat_xs[31:60,3] <- (seq(min(datafile$Synchrony), max(datafile$Synchrony), length.out = 30))
predictions_SN <- rep(NA, 30)
predictions_SNu <- rep(NA, 30)
predictions_SNl <- rep(NA, 30)
for(i in 1:30){predictions_SN[i] <- logistic(int+sum(slopes*test_mat_xs[i+30,]))
predictions_SNu[i] <- logistic(int+sum(c(slopes[1], SEs$ucl[12:13], slopes[4:7])*test_mat_xs[i+30,]))
predictions_SNl[i] <- logistic(int+sum(c(slopes[1], SEs$lcl[12:13], slopes[4:7])*test_mat_xs[i+30,]))
}

predictions_CI <- cbind(c(test_mat_xs[31:60,3], rev(test_mat_xs[31:60,3])), c(predictions_SNl, rev(predictions_SNu)))

plot(Survival ~ Synchrony, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Synchrony", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, at=seq(-3.3, 4.3, 1), cex.axis=1)
axis(2, las=1, cex.axis=1)
points(datafile$Synchrony, datafile$Survival, cex=datafile$coder2/50, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[31:60,3], predictions_SN, col="red", lwd=2)

# Beech
datafile <- datafile %>% group_by(Survival, beech) %>% mutate(coder3 = length(Survival))
datafile <- ungroup(datafile)

predictions_B <- rep(NA, 30)
predictions_Bu <- rep(NA, 30)
predictions_Bl <- rep(NA, 30)
factors <- c(rep(0,10), rep(B.1, 10), rep(B.0, 10))
factorsu <- c(rep(0,10), rep(SEs$ucl[16], 10), rep(SEs$ucl[15], 10))
factorsl <- c(rep(0,10), rep(SEs$lcl[16], 10), rep(SEs$lcl[15], 10))
factors2 <- c(rep(2,10), rep(1, 10), rep(0, 10))
for(i in 1:30){predictions_B[i] <- logistic(int+sum(slopes*test_mat_xs[i+60,])+factors[i])
predictions_Bu[i] <- logistic(int+sum(slopes*test_mat_xs[i+60,])+factorsu[i])
predictions_Bl[i] <- logistic(int+sum(slopes*test_mat_xs[i+60,])+factorsl[i])
}

plot(Survival ~ beech, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Beech", ylab="", 
     xlim=c(-0.5,2.5), axes=F, cex.lab=1.5)
axis(2, las=1, cex.axis=1)
axis(side=1, at=c(0,1,2), labels=c(0,1,2), cex.axis=1)
points(datafile$beech, datafile$Survival, cex=datafile$coder3/500, pch=16, col="grey70")
points(factors2, predictions_B, col="red", lwd=2, pch=16)
lines(x=c(2,2), y=c(predictions_Bl[1], predictions_Bu[1]), col=rgb(1, 0, 0,0.5), lwd=2)
lines(x=c(1,1), y=c(predictions_Bl[16], predictions_Bu[16]), col=rgb(1, 0, 0,0.5), lwd=2)
lines(x=c(0,0), y=c(predictions_Bl[25], predictions_Bu[25]), col=rgb(1, 0, 0,0.5), lwd=2)

# ST
datafile <- datafile %>% group_by(Survival, Spring_temp) %>% mutate(coder4 = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[91:120,6] <- seq(min(datafile$Spring_temp), max(datafile$Spring_temp), length.out=30)
predictions_ST <- rep(NA, 30)
predictions_STu <- rep(NA, 30)
predictions_STl <- rep(NA, 30)
for(i in 1:30){predictions_ST[i] <- logistic(int+sum(slopes*test_mat_xs[i+90,]))
predictions_STu[i] <- logistic(int+sum(c(slopes[1:5], SEs$ucl[18], slopes[7])*test_mat_xs[i+90,]))
predictions_STl[i] <- logistic(int+sum(c(slopes[1:5], SEs$lcl[18], slopes[7])*test_mat_xs[i+90,]))
}

predictions_CI <- cbind(c(test_mat_xs[91:120,6], rev(test_mat_xs[91:120,6])), c(predictions_STl, rev(predictions_STu)))

plot(Survival ~ Spring_temp, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Spring temperature", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-2,2,0.5))
axis(2, las=1, cex.axis=1)
title(ylab="Survival", line=2.5, cex.lab=1.5)
points(datafile$Spring_temp, datafile$Survival, cex=datafile$coder4/50, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[91:120,6], predictions_ST, col="red", lwd=2)

# SP
datafile <- datafile %>% group_by(Survival, spring_precip_t) %>% mutate(coder5 = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[121:150,5] <- seq(min(datafile$spring_precip_t), max(datafile$spring_precip_t), length.out=30)
predictions_SP <- rep(NA, 30)
predictions_SPu <- rep(NA, 30)
predictions_SPl <- rep(NA, 30)
for(i in 1:30){predictions_SP[i] <- logistic(int+sum(slopes*test_mat_xs[i+120,]))
predictions_SPu[i] <- logistic(int+sum(c(slopes[1:4], SEs$ucl[17], slopes[6:7])*test_mat_xs[i+120,]))
predictions_SPl[i] <- logistic(int+sum(c(slopes[1:4], SEs$lcl[17], slopes[6:7])*test_mat_xs[i+120,]))
}

predictions_CI <- cbind(c(test_mat_xs[121:150,5], rev(test_mat_xs[121:150,5])), c(predictions_SPl, rev(predictions_SPu)))

plot(Survival ~ spring_precip_t, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Spring precipitation", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-2,2,0.5))
axis(2, las=1, cex.axis=1)
points(datafile$spring_precip_t, datafile$Survival, cex=datafile$coder5/50, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[121:150,5], predictions_SP, col="red", lwd=2)

# WT
datafile <- datafile %>% group_by(Survival, winter_temp) %>% mutate(coder6 = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[151:180,7] <- seq(min(datafile$winter_temp), max(datafile$winter_temp), length.out=30)
predictions_WT <- rep(NA, 30)
predictions_WTu <- rep(NA, 30)
predictions_WTl <- rep(NA, 30)
for(i in 1:30){predictions_WT[i] <- logistic(int+sum(slopes*test_mat_xs[i+150,]))
predictions_WTu[i] <- logistic(int+sum(c(slopes[1:6], SEs$ucl[19])*test_mat_xs[i+150,]))
predictions_WTl[i] <- logistic(int+sum(c(slopes[1:6], SEs$lcl[19])*test_mat_xs[i+150,]))
}

predictions_CI <- cbind(c(test_mat_xs[151:180,7], rev(test_mat_xs[151:180,7])), c(predictions_WTl, rev(predictions_WTu)))

plot(Survival ~ winter_temp, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Winter temperature", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-4,1,0.5))
axis(2, las=1, cex.axis=1)
points(datafile$winter_temp, datafile$Survival, cex=datafile$coder6/50, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[151:180,7], predictions_WT, col="red", lwd=2)

# Immigrant
datafile <- datafile %>% group_by(Survival, Immigrant) %>% mutate(coder7 = length(Survival))
datafile <- ungroup(datafile)

predictions_I <- rep(NA, 30)
predictions_Iu <- rep(NA, 30)
predictions_Il <- rep(NA, 30)
factors <- c(rep(Resident,15), rep(0, 15))
factorsu <- c(rep(SEs$ucl[10],15), rep(0, 15))
factorsl <- c(rep(SEs$lcl[10],15), rep(0, 15))
factors2 <- c(rep(0,15), rep(2, 15))
for(i in 1:30){predictions_I[i] <- logistic(int+sum(slopes*test_mat_xs[i+180,])+factors[i])
predictions_Iu[i] <- logistic(int+sum(slopes*test_mat_xs[i+180,])+factorsu[i])
predictions_Il[i] <- logistic(int+sum(slopes*test_mat_xs[i+180,])+factorsl[i])
}

datafile$Immigrant2 <- datafile$Immigrant*2
plot(Survival ~ Immigrant2, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Immigrant", ylab="", xlim=c(-1,3),
     axes=F, cex.lab=1.5)
axis(2, las=1, cex.axis=1)
axis(side=1, at=c(0,2), labels=c("Resident", "Immigrant"), cex.axis=1)
title(ylab="Survival", line=2.5, cex.lab=1.5)
points(datafile$Immigrant2, datafile$Survival, cex=datafile$coder7/700, pch=16, col="grey70")
lines(x=c(0,0), y=c(predictions_Il[1], predictions_Iu[1]), col=rgb(1, 0, 0,0.5))
points(factors2, predictions_I, col="red", lwd=2, pch=16)

# Section
datafile <- datafile %>% group_by(Survival, Round) %>% mutate(coder8 = length(Survival))
datafile <- ungroup(datafile)

predictions_SE <- rep(NA, 9)
predictions_SEu <- rep(NA, 9)
predictions_SEl <- rep(NA, 9)
factors <- c(0,Rounds)
factorsu <- c(0,SEs$ucl[2:9])
factorsl <- c(0,SEs$lcl[2:9])
factors2 <- c("B","C", "CP", "EX", "MP", "O", "P", "SW", "W")
for(i in 1:9){predictions_SE[i] <- logistic(int+sum(slopes*test_mat_xs[i+210,])+factors[i])
predictions_SEu[i] <- logistic(int+sum(slopes*test_mat_xs[i+210,])+factorsu[i])
predictions_SEl[i] <- logistic(int+sum(slopes*test_mat_xs[i+210,])+factorsl[i])
}

plot(Survival ~ as.numeric(Round), data=datafile, type="n", ylim=c(-0.05,1.05), xlim=c(0,18), xlab="Section of woodland", ylab="", 
     axes=F, cex.lab=1.5)
axis(2, las=1, cex.axis=1)
axis(side=1, at=seq(1,18,2), labels=factors2, cex.axis=1, las=2)
points((as.numeric(datafile$Round)*2-1), datafile$Survival, cex=datafile$coder8/500, pch=16, col="grey70")
for(j in seq(1,18,2))lines(x=c(j,j), y=c(predictions_SEl[(j+1)/2], predictions_SEu[(j+1)/2]), col=rgb(1, 0, 0,0.5))
points(seq(1,18,2), predictions_SE, col="red", lwd=2, pch=16)

# Pop size
datafile <- datafile %>% group_by(Survival, pop_size) %>% mutate(coder9 = length(Survival))
datafile <- ungroup(datafile)

test_mat_xs[241:270,4] <- seq(min(datafile$pop_size), max(datafile$pop_size), length.out=30)
predictions_PS <- rep(NA, 30)
predictions_PSu <- rep(NA, 30)
predictions_PSl <- rep(NA, 30)
for(i in 1:30){predictions_PS[i] <- logistic(int+sum(slopes*test_mat_xs[i+240,]))
predictions_PSu[i] <- logistic(int+sum(c(slopes[1:3], SEs$ucl[14], slopes[5:7])*test_mat_xs[i+240,]))
predictions_PSl[i] <- logistic(int+sum(c(slopes[1:3], SEs$lcl[14], slopes[5:7])*test_mat_xs[i+240,]))
}

predictions_CI <- cbind(c(test_mat_xs[241:270,4], rev(test_mat_xs[241:270,4])), c(predictions_PSl, rev(predictions_PSu)))

plot(Survival ~ pop_size, data=datafile, type="n", ylim=c(-0.05,1.05), xlab="Population size", ylab="", 
     axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=round(seq(-1.8, 2, length.out = 6), 1))
axis(2, las=1, cex.axis=1)
points(datafile$pop_size, datafile$Survival, cex=datafile$coder9/100, pch=16, col="grey70")
polygon(x=predictions_CI[,1], y=predictions_CI[,2], col=rgb(1, 0, 0,0.25), border=F)
lines(test_mat_xs[241:270,4], predictions_PS, col="red", lwd=2)

dev.off()

### recruitment ####
load('final_recruitment_Feb19.RData')

# issues with plotting Round so change to number
Rx <- c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")
changer <- seq(1,9,1)
datafile$Round2 <- rep(NA, length(datafile$Round))
for(i in 1:length(Rx)){
  marker.x <- which(datafile$Round == Rx[i])
  datafile$Round2[marker.x] <- changer[i]
}

# save out slope values and mean values of each variable
slopes <- data.frame(fit = summary(final_recruitment)$coef[,1], 
                     ucl = summary(final_recruitment)$coef[,1] + 2*(summary(final_recruitment)$coef[,2]), 
                     lcl = summary(final_recruitment)$coef[,1] - 2*(summary(final_recruitment)$coef[,2]))

xs <- c(mean(datafile$April_hatch),
        mean(datafile$Synchrony),
        mean(datafile$Synchrony)^2,
        mean(datafile$Clutch_size),
        mean(datafile$Clutch_size)^2,
        mean(datafile$pop_size),
        mean(datafile$winter_precip_t),
        mean(datafile$spring_precip_t))

predict_recruit <- function(slopes, xs, factors){
  sum_of_slope <- sum(slopes[2:length(slopes)]*xs)
  u <- exp(slopes[1] + sum_of_slope + factors)
  return(u)
}


# like survival, use group.by instead of subsetting
png(filename="SOM1_Figure3.png", units="cm", height = 20, width = 18,
    res=400)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=T), heights=rep(7.5,3), widths = rep(6,3))
#layout.show(9)

par(mar=c(4,4,1,1))
# Hatch
datafile <- datafile %>% group_by(Num_recruited, April_hatch) %>% mutate(coder = length(Num_recruited))
datafile <- ungroup(datafile)

new_data <- data.frame(Num_recruited = rep(1:5,6),
                       April_hatch = seq(min(datafile$April_hatch),max(datafile$April_hatch),length=30), 
                       spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                       winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                       Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                       pop_size = mean(datafile$pop_size, na.rm=T),
                       beech = factor(2, levels=c(0,1,2)), Round = factor("B", levels=Rx), Year = 2000, Nest_box = factor("B14", levels=unique(datafile$Nest_box)))
predictions <- predict(final_recruitment, new_data, type="response", re.form = NA)

# can try and predict CIs manually

predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1,1],slopes[2,2],slopes[3:6,1],slopes[15:17,1]), xs=c(new_data$April_hatch[i], xs[2:8]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1,1],slopes[2,3],slopes[3:6,1],slopes[15:17,1]), xs=c(new_data$April_hatch[i], xs[2:8]), factors = slopes[19,1])}

plot(Num_recruited ~ April_hatch, data=datafile, ylab = "", 
     ylim=c(-0.5, 4), xlab="Hatch date", type='n', axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-3,3.5,0.5))
axis(2, las=1, cex.axis=1)
title(ylab="Number of recruits", line=2, cex.lab=1.5)
points(datafile$April_hatch, datafile$Num_recruited, cex=datafile$coder/100, pch=16, col="grey")
polygon(x=c(new_data$April_hatch, rev(new_data$April_hatch)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data$April_hatch, predictions, col="red", lwd=2)

# Synchrony
datafile <- datafile %>% group_by(Num_recruited, Synchrony) %>% mutate(coder2 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data2 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = seq(min(datafile$Synchrony), max(datafile$Synchrony), length=30), 
                        Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = 2, Round = "B", Year = 2000, Nest_box = "B14")
predictions2 <- predict(final_recruitment, new_data2, type="response", re.form = NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:3,1],slopes[4,2],slopes[5:6,1],slopes[15:17,1]), xs=c(xs[1],new_data2$Synchrony[i],new_data2$Synchrony[i]^2,xs[4:8]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:3,1],slopes[4,3],slopes[5:6,1],slopes[15:17,1]), xs=c(xs[1],new_data2$Synchrony[i],new_data2$Synchrony[i]^2,xs[4:8]), factors = slopes[19,1])}

plot(Num_recruited ~ Synchrony, data=datafile, ylab = "",ylim=c(-0.5, 4),  
     xlab="Synchrony", type='n', axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-3,4,0.5))
axis(2, las=1, cex.axis=1)
points(datafile$Synchrony, datafile$Num_recruited, cex=datafile$coder2/100, pch=16, col="grey")
polygon(x=c(new_data2$Synchrony, rev(new_data2$Synchrony)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data2$Synchrony, predictions2, col="red", lwd=2)

# Beech
datafile <- datafile %>% group_by(Num_recruited, beech) %>% mutate(coder3 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data3 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = seq(0,2,1), Round = "B", Year = 2000, Nest_box = "B14")
predictions3 <- predict(final_recruitment, new_data3, type="response", re.form = NA)
factors_u <- c(0,slopes[18,2], slopes[19,2])
factors_l <- c(0,slopes[18,3], slopes[19,3])
predictions_ucl <- rep(NA,3) 
for(i in 1:3){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,1]), xs=xs, factors = factors_u[i])}
predictions_lcl <- rep(NA,3)
for(i in 1:3){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,1]), xs=xs, factors = factors_l[i])}

plot(Num_recruited ~ beech, data=datafile, ylab = "", xlab="Beech", type='n', 
     ylim=c(-0.5, 4), xlim=c(-0.5,2.5), axes=F, cex.lab=1.5)
axis(side=1, labels=seq(0,2,1), at=seq(0,2,1), cex.axis=1)
axis(2, las=1, cex.axis=1)
points(datafile$beech, datafile$Num_recruited, cex=datafile$coder3/500, pch=16, col="grey")
for(j in 1:3){lines(x=rep(new_data3$beech[j],2),y=c(predictions_lcl[j], predictions_ucl[j]), col=rgb(1, 0, 0,0.5), lwd=2)}
points(new_data3$beech, predictions3, col="red", lwd=2, pch=16)

# Spring precip
datafile <- datafile %>% group_by(Num_recruited, spring_precip_t) %>% mutate(coder4 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data4 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = seq(min(datafile$spring_precip_t), max(datafile$spring_precip_t),length=30),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = 2, Round = "B", Year = 2000, Nest_box = "B14")
predictions4 <- predict(final_recruitment, new_data4, type="response", re.form = NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15,1],slopes[16,2], slopes[17,1]), xs=c(xs[1:6],new_data4$spring_precip_t[i], xs[8]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15,1],slopes[16,3], slopes[17,1]), xs=c(xs[1:6],new_data4$spring_precip_t[i], xs[8]), factors = slopes[19,1])}


plot(Num_recruited ~ spring_precip_t, data=datafile, ylab = "", ylim=c(-0.5, 4), 
     xlab="Spring precipitation", type='n', axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-2,2.5,0.5))
axis(2, las=1, cex.axis=1)
title(ylab="Number of recruits", line=2, cex.lab=1.5)
points(datafile$spring_precip_t, datafile$Num_recruited, cex=datafile$coder4/50, pch=16, col="grey")
polygon(x=c(new_data4$spring_precip_t, rev(new_data4$spring_precip_t)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data4$spring_precip_t, predictions4, col="red", lwd=2)

# Winter precip
datafile <- datafile %>% group_by(Num_recruited, winter_precip_t) %>% mutate(coder4 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data5 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t),
                        winter_precip_t = seq(min(datafile$winter_precip_t), max(datafile$winter_precip_t), length=30),
                        Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = 2, Round = "B", Year = 2000, Nest_box = "B14")
predictions5 <- predict(final_recruitment, new_data5, type="response", re.form = NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:16,1],slopes[17,2]), xs=c(xs[1:7],new_data5$winter_precip_t[i]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:16,1],slopes[17,3]), xs=c(xs[1:7],new_data5$winter_precip_t[i]), factors = slopes[19,1])}


plot(Num_recruited ~ winter_precip_t, data=datafile, ylab = "", ylim=c(-0.5, 4), 
     xlab="Winter precipitation", type='n', axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-2,2.5,0.5))
axis(2, las=1, cex.axis=1)
title(ylab="Number of recruits", line=2, cex.lab=1.5)
points(datafile$winter_precip_t, datafile$Num_recruited, cex=datafile$coder5/50, pch=16, col="grey")
polygon(x=c(new_data5$winter_precip_t, rev(new_data5$winter_precip_t)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data5$winter_precip_t, predictions5, col="red", lwd=2)

# Clutch size
datafile <- datafile %>% group_by(Num_recruited, Clutch_size) %>% mutate(coder6 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data6 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = 0, Clutch_size = seq(min(datafile$Clutch_size), max(datafile$Clutch_size), length=30), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = 2, Round = "B", Year = 2000, Nest_box = "B14")
predictions6 <- predict(final_recruitment, new_data6, type="response", re.form = NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:5,1],slopes[6,2],slopes[15:17,1]), xs=c(xs[1:3],new_data6$Clutch_size[i],new_data6$Clutch_size[i]^2, xs[6:8]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:5,1],slopes[6,3],slopes[15:17,1]), xs=c(xs[1:3],new_data6$Clutch_size[i],new_data6$Clutch_size[i]^2, xs[6:8]), factors = slopes[19,1])}

plot(Num_recruited ~ Clutch_size, data=datafile, ylab = "", ylim=c(-0.5, 4), 
     xlab="Clutch size", type='n', axes=F, cex.lab=1.5)
axis(1, cex.axis=1, at=seq(-4,7.5,1))
axis(2, las=1, cex.axis=1)
points(datafile$Clutch_size, datafile$Num_recruited, cex=datafile$coder6/500, pch=16, col="grey")
polygon(x=c(new_data6$Clutch_size, rev(new_data6$Clutch_size)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data6$Clutch_size, predictions6, col="red", lwd=2)

# Section of woodland
datafile <- datafile %>% group_by(Num_recruited, Round) %>% mutate(coder7 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data7 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = mean(datafile$pop_size, na.rm=T),
                        beech = 2, Round = c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W"), Year = 2000, Nest_box = "B14")
predictions7 <- predict(final_recruitment, new_data7, type="response", re.form = NA)
factors_u <- c(0,slopes[7:14,2])
factors_l <- c(0,slopes[7:14,3])
predictions_ucl <- rep(NA,9) 
for(i in 1:9){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,1]), xs=xs, factors = factors_u[i] + slopes[19,1])}
predictions_lcl <- rep(NA,9)
for(i in 1:9){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,1]), xs=xs, factors = factors_l[i] + slopes[19,1])}

plot(Num_recruited ~ Round2, data=datafile, ylab = "", xlab="Section of woodland",ylim=c(-0.5, 4), 
     type='n', axes=F, cex.lab=1.5)
axis(side=1, at=changer, labels=Rx, cex.axis=1, las=2)
axis(2, las=1, cex.axis=1)
points(datafile$Round2, datafile$Num_recruited, cex=datafile$coder7/200, pch=16, col="grey")
for(j in 1:9){lines(x=rep(new_data7$Round[j],2),y=c(predictions_lcl[j], predictions_ucl[j]), col=rgb(1, 0, 0,0.5), lwd=2)}
points(changer, predictions7, col="red", lwd=2, pch=16)

# Population size
datafile <- datafile %>% group_by(Num_recruited, pop_size) %>% mutate(coder8 = length(Num_recruited))
datafile <- ungroup(datafile)

new_data8 <- data.frame(April_hatch = mean(datafile$April_hatch), 
                        spring_precip_t = mean(datafile$spring_precip_t, na.rm=T),
                        winter_precip_t = mean(datafile$winter_precip_t, na.rm=T),
                        Synchrony = 0, Clutch_size = mean(datafile$Clutch_size, na.rm=T), 
                        pop_size = seq(min(datafile$pop_size),max(datafile$pop_size), length=30),
                        beech = 2, Round = "B", Year = 2000, Nest_box = "B14")
predictions8 <- predict(final_recruitment, new_data8, type="response", re.form = NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,2]), xs=c(xs[1:5], new_data8$pop_size[i], xs[7:8]), factors = slopes[19,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_recruit(slopes=c(slopes[1:6,1],slopes[15:17,3]), xs=c(xs[1:5], new_data8$pop_size[i], xs[7:8]), factors = slopes[19,1])}

plot(Num_recruited ~ pop_size, data=datafile, ylab = "", ylim=c(-0.5, 4), 
     xlab="Population size", type='n', axes=F, cex.lab=1.5)
axis(side=1, cex.axis=1, at=round(seq(-1.8,2,length.out = 6),1))
axis(2, las=1, cex.axis=1)
title(ylab="Number of recruits", line=2, cex.lab=1.5)
points(datafile$pop_size, datafile$Num_recruited, cex=datafile$coder8/200, pch=16, col="grey")
polygon(x=c(new_data8$pop_size, rev(new_data8$pop_size)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data8$pop_size, predictions8, col="red", lwd=2)

dev.off()

### development ####
load('final_development_Feb19.RData')

slopes <- data.frame(fit = summary(final_development)$coef[,1], 
                     ucl = summary(final_development)$coef[,1] + 2*(summary(final_development)$coef[,2]), 
                     lcl = summary(final_development)$coef[,1] - 2*(summary(final_development)$coef[,2]))

xs <- c(mean(datafile_t$April_hatch),
        mean(datafile_t$Synchrony),
        mean(datafile_t$Synchrony)^2,
        mean(datafile_t$Spring_temp_t1),
        mean(datafile_t$pop_size_t1))

datafile_t <- datafile_t[,-1]

predict_devel <- function(slopes, xs, factors){
  sum_of_slope <- sum(slopes[2:length(slopes)]*xs)
  u <- slopes[1] + sum_of_slope + factors
  return(u)
}

png(filename="SOM1_Figure2.png", units="cm", height = 20, width = 18,
    res=400)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=T), heights=rep(7.5,3), widths = rep(6,3))
#layout.show(9)

par(mar=c(4,4,1,1))

# Hatch
datafile_t <- datafile_t %>% group_by(April_hatch_t1, April_hatch) %>% mutate(coder = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data1 <- data.frame(April_hatch = seq(min(datafile_t$April_hatch), max(datafile_t$April_hatch), length=30), 
                        Spring_temp_t1 = rep(mean(datafile_t$Spring_temp_t1), 30), 
                        Synchrony = rep(mean(datafile_t$Synchrony), 30),
                        Round_t1 = rep("B"), Age = rep(2, 30), 
                        pop_size_t1 = rep(mean(datafile_t$pop_size_t1), 30))
prediction1 <- predict(final_development, new_data1, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,2],slopes[3:5,1],slopes[7,1]), xs=c(new_data1$April_hatch[i], xs[2:5]), factors = slopes[6,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,3],slopes[3:5,1],slopes[7,1]), xs=c(new_data1$April_hatch[i], xs[2:5]), factors = slopes[6,1])}

plot(April_hatch_t1 ~ April_hatch, data=datafile_t, xlab="Scaled hatch date (t)", ylab="", 
     type='n', axes=F, cex.lab =1.5)
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date (t+1)", cex.lab=1.5, line=2.5)
points(datafile_t$April_hatch, datafile_t$April_hatch_t1, cex=datafile_t$coder/5, pch=16, col="grey")
polygon(x=c(new_data1$April_hatch, rev(new_data1$April_hatch)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data1$April_hatch, prediction1, col="red", lwd=2)

# Synchrony
datafile_t <- datafile_t %>% group_by(April_hatch_t1, Synchrony) %>% mutate(coder2 = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data2 <- data.frame(April_hatch = rep(mean(datafile_t$April_hatch), 30), 
                        Synchrony = seq(min(datafile_t$Synchrony), max(datafile_t$Synchrony), length=30),
                        Spring_temp_t1 = rep(mean(datafile_t$Spring_temp_t1), 30),
                        Round_t1 = rep("B"), Age = rep(2, 30), pop_size_t1 = rep(mean(datafile_t$pop_size_t1), 30))
prediction2 <- predict(final_development, new_data2, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:2,1],slopes[3,2],slopes[4:5,1],slopes[7,1]), xs=c(xs[1], new_data2$Synchrony[i], new_data2$Synchrony[i]^2,xs[4:5]), factors = slopes[6,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:2,1],slopes[3,3],slopes[4:5,1],slopes[7,1]), xs=c(xs[1], new_data2$Synchrony[i], new_data2$Synchrony[i]^2,xs[4:5]), factors = slopes[6,1])}

plot(April_hatch_t1 ~ Synchrony, data=datafile_t, xlab="Synchrony (t)", ylab="", 
     type='n', axes=F, cex.lab =1.5)
axis(1)
axis(2, las = 1)
points(datafile_t$Synchrony, datafile_t$April_hatch_t1, cex=datafile_t$coder2/5, pch=16, col="grey")
polygon(x=c(new_data2$Synchrony, rev(new_data2$Synchrony)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data2$Synchrony, prediction2, col="red", lwd=2)

# Spring temperature
datafile_t <- datafile_t %>% group_by(April_hatch_t1, Spring_temp_t1) %>% mutate(coder3 = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data3 <- data.frame(April_hatch = rep(mean(datafile_t$April_hatch), 30), Synchrony = rep(mean(datafile_t$Synchrony), 30),
                        Spring_temp_t1 = seq(min(datafile_t$Spring_temp_t1), max(datafile_t$Spring_temp_t1), length=30),
                        Round_t1 = rep("B"), Age = rep(2, 30),pop_size_t1 = rep(mean(datafile_t$pop_size_t1), 30))
prediction3 <- predict(final_development, new_data3, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:4,1],slopes[5,2],slopes[7,1]), xs=c(xs[1:3], new_data3$Spring_temp_t1[i], xs[5]), factors = slopes[6,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:4,1],slopes[5,3],slopes[7,1]), xs=c(xs[1:3], new_data3$Spring_temp_t1[i], xs[5]), factors = slopes[6,1])}

plot(April_hatch_t1 ~ Spring_temp_t1, data=datafile_t, xlab="Spring temperature (t+1)", ylab="", 
     type='n', axes=F, cex.lab=1.5)
axis(1)
axis(2, las=1)
points(datafile_t$Spring_temp_t1, datafile_t$April_hatch_t1, cex=datafile_t$coder3/10, pch=16, col="grey")
polygon(x=c(new_data3$Spring_temp_t1, rev(new_data3$Spring_temp_t1)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data3$Spring_temp_t1, prediction3, col="red", lwd=2)

# Section of woodland
datafile_t <- datafile_t %>% group_by(April_hatch_t1, Round_t1) %>% mutate(coder4 = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data4 <- data.frame(April_hatch = rep(mean(datafile_t$April_hatch), length=9), 
                        Synchrony = rep(mean(datafile_t$Synchrony), 9),
                        Spring_temp_t1 = rep(mean(datafile_t$Spring_temp_t1), 9),
                        Round_t1 = c("EX", "B", "MP", "SW", "C", "CP", "O", "P", "W"), Age = rep(2, 9), 
                        pop_size_t1 = rep(mean(datafile_t$pop_size_t1), 9))
prediction4 <- predict(final_development, new_data4, type="response", re.form=NA)
factors_u <- c(0, slopes[8:15,2])
factors_l <- c(0, slopes[8:15,3])
predictions_ucl <- rep(NA,9) 
for(i in 1:9){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,1]), xs=c(xs[1:5]), factors = factors_u[i] + slopes[6,1])}
predictions_lcl <- rep(NA,9)
for(i in 1:9){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,1]), xs=c(xs[1:5]), factors = factors_l[i] + slopes[6,1])}

plot(April_hatch_t1 ~ as.numeric(Round_t1), data=datafile_t, xlab="Section of woodland (t+1)", 
     ylab="", type='n', axes=F, cex.lab=1.5)
axis(side=1, at=seq(1,9,1), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W"), las=2)
axis(2, las=1)
title(ylab="Scaled hatch date (t+1)", line=2.5, cex.lab=1.5)
points(as.numeric(datafile_t$Round_t1), datafile_t$April_hatch_t1, cex=datafile_t$coder4/50, col="grey", pch=16)
for(j in 1:9){lines(x=rep(j, 2), y=c(predictions_lcl[j], predictions_ucl[j]), col = rgb(1, 0, 0,0.25), lwd=2)}
points(new_data4$Round_t1, prediction4, col="red", pch=16)

# Population size
datafile_t <- datafile_t %>% group_by(April_hatch_t1, pop_size) %>% mutate(coder5 = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data5 <- data.frame(April_hatch = rep(mean(datafile_t$April_hatch), length=30), 
                        Synchrony = rep(mean(datafile_t$Synchrony), 30),
                        Spring_temp_t1 = rep(mean(datafile_t$Spring_temp_t1), 30),
                        Round_t1 = rep("B"), Age = rep(2, 30),
                        pop_size_t1 = seq(min(datafile_t$pop_size_t1), max(datafile_t$pop_size_t1), length=30))
prediction5 <- predict(final_development, new_data5, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,2]), xs=c(xs[1:4], new_data5$pop_size_t1[i]), factors = slopes[6,1])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,3]), xs=c(xs[1:4], new_data5$pop_size_t1[i]), factors = slopes[6,1])}

plot(April_hatch_t1 ~ pop_size_t1, data=datafile_t, xlab="Population size (t+1)", ylab="", 
     type="n", axes=F, cex.lab=1.5)
axis(1, at=seq(-1.8,2.2,0.2))
axis(2, las=1)
points(datafile_t$pop_size_t1, datafile_t$April_hatch_t1, cex=datafile_t$coder5/10, pch=16, col="grey")
polygon(x=c(new_data5$pop_size_t1, rev(new_data5$pop_size_t1)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, col = rgb(1, 0, 0,0.25))
lines(new_data5$pop_size_t1, prediction5, col="red", lwd=2)

# Age
datafile_t <- datafile_t %>% group_by(April_hatch_t1, Age) %>% mutate(coder6 = length(April_hatch_t1))
datafile_t <- ungroup(datafile_t)

new_data6 <- data.frame(April_hatch = rep(mean(datafile_t$April_hatch), length=30), 
                        Synchrony = rep(mean(datafile_t$Synchrony), 30),
                        Spring_temp_t1 = rep(mean(datafile_t$Spring_temp_t1), 30),
                        Round_t1 = "B", Age = c(rep(1,15), rep(2,15)), 
                        pop_size_t1 = rep(mean(datafile_t$pop_size_t1), 30))
prediction6 <- predict(final_development, new_data6, type="response", re.form=NA)
factors_u <- c(rep(1,15), rep(slopes[6,2], 15))
factors_l <- c(rep(1,15), rep(slopes[6,3], 15))
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,1]), xs=c(xs[1:5]), factors = factors_u[i])}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:5,1],slopes[7,1]), xs=c(xs[1:5]), factors = factors_l[i])}

plot(April_hatch_t1 ~ as.numeric(Age), data=datafile_t, xlab="Age", 
     ylab="", type='n', axes=F, xlim=c(0.5,2.5), cex.axis=1.5, cex.lab=1.5)
axis(side=1, at=seq(1,2,1), labels=seq(1,2,1))
axis(2, las=1)
points(as.numeric(datafile_t$Age), datafile_t$April_hatch_t1, cex=datafile_t$coder6/50, col="grey", pch=16)
#lines(x=rep(1, 2), y=c(predictions_lcl[1], predictions_ucl[1]), col = rgb(1, 0, 0,0.25), lwd=2)
lines(x=rep(2, 2), y=c(predictions_lcl[16], predictions_ucl[16]), col = rgb(1, 0, 0,0.5), lwd=2)
points(new_data6$Age, prediction6, col="red", cex=0.75, pch=16)

dev.off()

### inheritance ####
load('final_inheritance_QG_Feb19.RData')

inher_output <- summary(final_inheritance_QG)$solutions

load('final_inheritance_PO_Feb19.RData')

PO_output <- read.csv('model_output_inher_po_tot_pop.csv', header=T)

SEs <- data.frame(est = inher_output[1:14,2],
                  ucl = inher_output[1:14,4],
                  lcl = inher_output[1:14,3])
SEsPO <- data.frame(est = PO_output[1:13,2],
                    ucl = PO_output[1:13,2]+(2*PO_output[1:13,3]),
                    lcl = PO_output[1:13,2]-(2*PO_output[1:13,3]))

# GLMM basis too so also predict like development

xs <- c(mean(datafile$Spring_temp),
        mean(datafile$winter_temp),
        mean(datafile$winter_precip_t),
        mean(datafile$spring_precip_t),
        mean(datafile$pop_size))

predict_inher <- function(slopes, xs, factors){
  sum_of_slope <- sum(slopes[2:length(slopes)]*xs)
  u <- slopes[1] + sum_of_slope + factors
  return(u)
}

png(filename="SOM1_Figure4.png", units="cm", height = 20, width = 18,
    res=400)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=T), heights=rep(7.5,3), widths = rep(6,3))
#layout.show(9)

par(mar=c(5,4,1,1))

# Spring temperature
datafile_I <- datafile_I %>% group_by(April_hatch, Spring_temp) %>% mutate(coder = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data1 <- data.frame(April_hatch = 0,
                        Spring_temp = seq(min(datafile_I$Spring_temp), max(datafile_I$Spring_temp), length=30),
                        winter_temp = mean(datafile_I$winter_temp, na.rm=T),
                        winter_precip_t = mean(datafile_I$winter_precip_t),
                        spring_precip_t = mean(datafile_I$spring_precip_t),
                        pop_size = mean(datafile_I$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal = factor(datafile_I$animal[1000], levels = unique(datafile_I$ID)), 
                        ID = factor(datafile_I$ID[1000], levels = unique(datafile_I$ID)), 
                        MOTHER = factor(datafile_I$MOTHER[1], levels = unique(datafile_I$MOTHER)))
predictions1 <- predict.MCMCglmm(final_inheritance_QG, new_data1, marginal=final_inheritance_QG$Random$formula,
                                 type="response", interval="confidence", level=0.95, it=NULL, 
                                 posterior="all", verbose=FALSE)

plot(April_hatch ~ Spring_temp, data=datafile_I, ylab="", xlab="Spring temperature", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I$Spring_temp, datafile_I$April_hatch, cex=datafile_I$coder/20, pch=16, col="grey")
polygon(x=c(new_data1$Spring_temp, rev(new_data1$Spring_temp)), y=c(predictions1[,2], rev(predictions1[,3])), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data1$Spring_temp, predictions1[,1], col="red", lwd=2)

# Spring precipitation
datafile_I <- datafile_I %>% group_by(April_hatch, spring_precip_t) %>% mutate(coder2 = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data2 <- data.frame(April_hatch = seq(min(datafile_I$April_hatch), max(datafile_I$April_hatch), length=30),
                        Spring_temp = mean(datafile_I$Spring_temp),
                        winter_temp = mean(datafile_I$winter_temp, na.rm=T),
                        winter_precip_t = mean(datafile_I$winter_precip_t),
                        spring_precip_t = seq(min(datafile_I$spring_precip_t), max(datafile_I$spring_precip_t), length=30),
                        pop_size = mean(datafile_I$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal=datafile_I$animal[1000], Age = factor(2, levels=c(1,2)), ID = datafile_I$ID[1000], 
                        MOTHER = datafile_I$MOTHER[1])
predictions2 <- predict(final_inheritance_QG, new_data2, marginal=final_inheritance_QG$Random$formula,
                        type="response", interval="confidence", level=0.95, it=NULL, 
                        posterior="all", verbose=FALSE)

plot(April_hatch ~ spring_precip_t, data=datafile_I, ylab="", xlab="Spring precipitation", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1, at=seq(-2,2.5,0.5))
axis(2, las=1)
points(datafile_I$spring_precip_t, datafile_I$April_hatch, cex=datafile_I$coder2/20, pch=16, col="grey")
polygon(x=c(new_data2$spring_precip_t, rev(new_data2$spring_precip_t)), y=c(predictions2[,2], rev(predictions2[,3])), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data2$spring_precip_t, predictions2[,1], col="red", lwd=2)

# Winter temperature
datafile_I <- datafile_I %>% group_by(April_hatch, winter_temp) %>% mutate(coder3 = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data3 <- data.frame(April_hatch = seq(min(datafile_I$April_hatch), max(datafile_I$April_hatch), length=30),
                        Spring_temp = mean(datafile_I$Spring_temp),
                        winter_temp = seq(min(datafile_I$winter_temp), max(datafile_I$winter_temp), length=30),
                        winter_precip_t = mean(datafile_I$winter_precip_t),
                        spring_precip_t = mean(datafile_I$spring_precip_t, na.rm=T),
                        pop_size = mean(datafile_I$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal=datafile_I$animal[1000], Age = factor(2, levels=c(1,2)), ID = datafile_I$ID[1000], 
                        MOTHER = datafile_I$MOTHER[1])
predictions3 <- predict(final_inheritance_QG, new_data3, marginal=final_inheritance_QG$Random$formula,
                        type="response", interval="confidence", level=0.95, it=NULL, 
                        posterior="all", verbose=FALSE)

plot(April_hatch ~ winter_temp, data=datafile_I, ylab="", xlab="Winter temperature", 
     pch=16, col="grey", type='n', axes=F, cex.lab=1.5)
axis(1, at=seq(-4,1.5,0.5))
axis(2, las=1)
points(datafile_I$winter_temp, datafile_I$April_hatch, cex=datafile_I$coder3/20, pch=16, col="grey")
polygon(x=c(new_data3$winter_temp, rev(new_data3$winter_temp)), y=c(predictions3[,2], rev(predictions3[,3])), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data3$winter_temp, predictions3[,1], col="red", lwd=2)

# Winter precipitation
datafile_I <- datafile_I %>% group_by(April_hatch, winter_precip_t) %>% mutate(coder4 = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data4 <- data.frame(April_hatch = seq(min(datafile_I$April_hatch), max(datafile_I$April_hatch), length=30),
                        Spring_temp = mean(datafile_I$Spring_temp),
                        winter_temp = mean(datafile_I$winter_temp, na.rm=T),
                        winter_precip_t = seq(min(datafile_I$winter_precip_t), max(datafile_I$winter_precip_t), length=30),
                        spring_precip_t = mean(datafile_I$spring_precip_t, na.rm=T),
                        pop_size = mean(datafile_I$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal=datafile_I$animal[1000], Age = factor(2, levels=c(1,2)), ID = datafile_I$ID[1000], 
                        MOTHER = datafile_I$MOTHER[1])
predictions4 <- predict(final_inheritance_QG, new_data4, marginal=final_inheritance_QG$Random$formula,
                        type="response", interval="confidence", level=0.95, it=NULL, 
                        posterior="all", verbose=FALSE)

plot(April_hatch ~ winter_precip_t, data=datafile_I, ylab="", xlab="Winter precipitation", 
     pch=16, axes=F, cex.lab=1.5, type='n')
axis(1, seq(-2,2.25,0.5))
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I$winter_precip_t, datafile_I$April_hatch, cex=datafile_I$coder4/20, pch=16, col="grey")
polygon(x=c(new_data4$winter_precip_t, rev(new_data4$winter_precip_t)), y=c(predictions4[,2], rev(predictions4[,3])), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data4$winter_precip_t, predictions4[,1], col="red", lwd=2)

# Section of woodland
datafile_I <- datafile_I %>% group_by(April_hatch, Round) %>% mutate(coder5 = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data5 <- data.frame(April_hatch = seq(min(datafile_I$April_hatch), max(datafile_I$April_hatch), length=9),
                        Spring_temp = mean(datafile_I$Spring_temp),
                        winter_temp = mean(datafile_I$winter_temp, na.rm=T),
                        winter_precip_t = mean(datafile_I$winter_precip_t),
                        spring_precip_t = mean(datafile_I$spring_precip_t, na.rm=T),
                        pop_size = mean(datafile_I$pop_size, na.rm=T), 
                        Round = factor(seq(1:9),levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal=datafile_I$animal[1000], Age = factor(2, levels=c(1,2)), ID = datafile_I$ID[1000], 
                        MOTHER = datafile_I$MOTHER[1])
predictions5 <- predict(final_inheritance_QG, new_data5, marginal=final_inheritance_QG$Random$formula,
                        type="response", interval="confidence", level=0.95, it=NULL, 
                        posterior="all", verbose=FALSE)

plot(April_hatch ~ as.numeric(Round), data=datafile_I, ylab="", xlab="Section of woodland", pch=16, 
     type='n',axes=F, cex.lab=1.5, xaxt='n')
axis(side=1, at=seq(1,9,1), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W"), las=2)
axis(2, las=1)
points(as.numeric(datafile_I$Round), datafile_I$April_hatch, cex=datafile_I$coder5/50, pch=16, col="grey")
for(j in 1:9){lines(x=rep(new_data5$Round[j], 2), y=c(predictions5[j,2], predictions5[j,3]), col=rgb(1, 0, 0,0.5), lwd=2)}
points(new_data5$Round, predictions5[,1], col="red", pch=16)

# Population size
datafile_I <- datafile_I %>% group_by(April_hatch, pop_size) %>% mutate(coder6 = length(April_hatch))
datafile_I <- ungroup(datafile_I)

new_data6 <- data.frame(April_hatch = seq(min(datafile_I$April_hatch), max(datafile_I$April_hatch), length=30),
                        Spring_temp = mean(datafile_I$Spring_temp),
                        winter_temp = mean(datafile_I$winter_temp, na.rm=T),
                        winter_precip_t = mean(datafile_I$winter_precip_t),
                        spring_precip_t = mean(datafile_I$spring_precip_t, na.rm=T),
                        pop_size = seq(min(datafile_I$pop_size), max(datafile_I$pop_size), length=30), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")), 
                        animal=datafile_I$animal[1000], Age = factor(2, levels=c(1,2)), ID = datafile_I$ID[1000], 
                        MOTHER = datafile_I$MOTHER[1])
predictions6 <- predict(final_inheritance_QG, new_data6, marginal=final_inheritance_QG$Random$formula,
                        type="response", interval="confidence", level=0.95, it=NULL, 
                        posterior="all", verbose=FALSE)

plot(April_hatch ~ pop_size, data=datafile_I, ylab="", xlab="Population size", 
     pch=16, type='n', cex.lab=1.5, axes=F)
axis(1, at=seq(-1.8,2.2,0.2))
axis(2, las=1)
points(datafile_I$pop_size, datafile_I$April_hatch, cex=datafile_I$coder6/20, pch=16, col="grey")
polygon(x=c(new_data6$pop_size, rev(new_data6$pop_size)), y=c(predictions6[,2], rev(predictions6[,3])), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data6$pop_size, predictions6[,1], col="red", lwd=2)


dev.off()

### PO ####
mother_hatch <- rep(NA, length(datafile_I$April_hatch))
for(j in 1:length(datafile_I$April_hatch)){
  temp <- filter(datafile_I, F_ID == datafile_I$F_ID[j])
  marker <- which(as.character(datafile_I$F_ID) 
                  == as.character(temp$MOTHER[1]))
  ifelse(length(marker) == 0, mother_hatch[j] <- NA,
         mother_hatch[j] <- mean(datafile_I$April_hatch[marker]))
}

datafile_I$Mother_hatch <- mother_hatch
datafile_I2 <- datafile_I[-which(is.na(datafile_I$Mother_hatch)),]
datafile_I2 <- filter(datafile_I2, Immigrant == 0)

slopes <- data.frame(est = PO_output[1:13,2],
                     ucl = PO_output[1:13,2]+(2*PO_output[1:13,3]),
                     lcl = PO_output[1:13,2]-(2*PO_output[1:13,3]))

xs <- c(mean(datafile_I2$Mother_hatch, na.rm=T),
        mean(datafile_I2$pop_size),
        mean(datafile_I2$Spring_temp),
        mean(datafile_I2$spring_precip_t))

png(filename="SOM1_Figure5.png", units="cm", height = 20, width = 18,
    res=400)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=T), heights=rep(7.5,3), widths = rep(6,3))
#layout.show(9)

par(mar=c(5,4,1,1))

# Mother hatch
datafile_I2 <- datafile_I2 %>% group_by(April_hatch, Mother_hatch) %>% mutate(coder = length(April_hatch))
datafile_I2 <- ungroup(datafile_I2)

new_data1 <- data.frame(April_hatch = seq(min(datafile_I2$April_hatch, na.rm=T), max(datafile_I2$April_hatch, na.rm=T), length=30),
                        Mother_hatch = seq(min(datafile_I2$Mother_hatch, na.rm=T), max(datafile_I2$Mother_hatch, na.rm=T), length=30),
                        Spring_temp = mean(datafile_I2$Spring_temp, na.rm=T),
                        spring_precip_t = mean(datafile_I2$spring_precip_t),
                        pop_size = mean(datafile_I2$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")))
prediction1 <- predict(final_inheritance_PO, new_data1, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,2],slopes[3:5,1]), 
                                                   xs=c(new_data1$Mother_hatch[i], xs[2:4]), factors=0)}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,3],slopes[3:5,1]), 
                                                   xs=c(new_data1$Mother_hatch[i], xs[2:4]), factors=0)}

plot(April_hatch ~ Mother_hatch, data=datafile_I2, ylab="", xlab="Mother hatch", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I2$Mother_hatch, datafile_I2$April_hatch, cex=datafile_I2$coder/20, pch=16, col="grey")
polygon(x=c(new_data1$Mother_hatch, rev(new_data1$Mother_hatch)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data1$Mother_hatch, prediction1, col="red", lwd=2)


# Spring temperature
datafile_I2 <- datafile_I2 %>% group_by(April_hatch, Spring_temp) %>% mutate(coder = length(April_hatch))
datafile_I2 <- ungroup(datafile_I2)

new_data2 <- data.frame(April_hatch = seq(min(datafile_I2$April_hatch), max(datafile_I2$April_hatch), length=30),
                        Mother_hatch = mean(datafile_I2$Mother_hatch, na.rm=T),
                        Spring_temp = seq(min(datafile_I2$Spring_temp, na.rm=T), max(datafile_I2$Spring_temp), length=30),
                        spring_precip_t = mean(datafile_I2$spring_precip_t),
                        pop_size = mean(datafile_I2$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")))
prediction2 <- predict(final_inheritance_PO, new_data2, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:2,1],slopes[3,2],slopes[4:5,1]), 
                                                   xs=c(xs[1:2],new_data2$Spring_temp[i], xs[4]), factors=0)}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:2,1],slopes[3,3],slopes[4:5,1]), 
                                                   xs=c(xs[1:2],new_data2$Spring_temp[i], xs[4]), factors=0)}

plot(April_hatch ~ Spring_temp, data=datafile_I2, ylab="", xlab="Spring temperature", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I2$Spring_temp, datafile_I2$April_hatch, cex=datafile_I2$coder/20, pch=16, col="grey")
polygon(x=c(new_data2$Spring_temp, rev(new_data2$Spring_temp)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data2$Spring_temp, prediction2, col="red", lwd=2)

# Spring precipitation
datafile_I2 <- datafile_I2 %>% group_by(April_hatch, spring_precip_t) %>% mutate(coder = length(April_hatch))
datafile_I2 <- ungroup(datafile_I2)

new_data3 <- data.frame(April_hatch = seq(min(datafile_I2$April_hatch), max(datafile_I2$April_hatch), length=30),
                        Mother_hatch = mean(datafile_I2$Mother_hatch, na.rm=T),
                        Spring_temp = mean(datafile_I2$Spring_temp),
                        spring_precip_t = seq(min(datafile_I2$spring_precip_t), max(datafile_I2$spring_precip_t), length=30),
                        pop_size = mean(datafile_I2$pop_size, na.rm=T), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")))
prediction3 <- predict(final_inheritance_PO, new_data3, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1:3,1],slopes[4,2],slopes[5,1]), 
                                                   xs=c(xs[1:3],new_data3$spring_precip_t[i]), factors=0)}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1:3,1],slopes[4,3],slopes[5,1]), 
                                                   xs=c(xs[1:3],new_data3$spring_precip_t[i]), factors=0)}

plot(April_hatch ~ spring_precip_t, data=datafile_I2, ylab="", xlab="Spring precipitation", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I2$spring_precip_t, datafile_I2$April_hatch, cex=datafile_I2$coder/20, pch=16, col="grey")
polygon(x=c(new_data3$spring_precip_t, rev(new_data3$spring_precip_t)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data3$spring_precip_t, prediction3, col="red", lwd=2)

# Section of woodland
datafile_I2 <- datafile_I2 %>% group_by(April_hatch, Round) %>% mutate(coder5 = length(April_hatch))
datafile_I2 <- ungroup(datafile_I2)

new_data5 <- data.frame(April_hatch = seq(min(datafile_I2$April_hatch), max(datafile_I2$April_hatch), length=9),
                        Mother_hatch = mean(datafile_I2$Mother_hatch, na.rm=T),
                        Spring_temp = mean(datafile_I2$Spring_temp),
                        spring_precip_t = mean(datafile_I2$spring_precip_t, na.rm=T),
                        pop_size = mean(datafile_I2$pop_size, na.rm=T), 
                        Round = factor(seq(1:9),levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")))
prediction4 <- predict(final_inheritance_PO, new_data5, type="response", re.form=NA)
factors_u <- c(0, slopes[6:13,2])
factors_l <- c(0, slopes[6:13,3])
predictions_ucl <- rep(NA,9) 
for(i in 1:9){predictions_ucl[i] <- predict_devel(slopes=slopes[,1], 
                                                  xs=xs, factors=factors_u[i])}
predictions_lcl <- rep(NA,9)
for(i in 1:9){predictions_lcl[i] <- predict_devel(slopes=slopes[,1], 
                                                  xs=xs, factors=factors_l[i])}

plot(April_hatch ~ as.numeric(Round), data=datafile_I2, ylab="", xlab="Section of woodland", pch=16, 
     type='n',axes=F, cex.lab=1.5, xaxt='n')
axis(side=1, at=seq(1,9,1), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W"), las=2)
axis(2, las=1)
points(as.numeric(datafile_I2$Round), datafile_I2$April_hatch, cex=datafile_I2$coder5/50, pch=16, col="grey")
for(j in 1:9){lines(x=rep(new_data5$Round[j], 2), y=c(predictions_ucl[j], predictions_lcl[j]), col=rgb(1, 0, 0,0.5), lwd=2)}
points(new_data5$Round, prediction4, col="red", pch=16)

# Population size
datafile_I2 <- datafile_I2 %>% group_by(April_hatch, pop_size) %>% mutate(coder6 = length(April_hatch))
datafile_I2 <- ungroup(datafile_I2)

new_data6 <- data.frame(April_hatch = seq(min(datafile_I2$April_hatch), max(datafile_I2$April_hatch), length=30),
                        Mother_hatch = mean(datafile_I2$Mother_hatch, na.rm=T),
                        Spring_temp = mean(datafile_I2$Spring_temp),
                        spring_precip_t = mean(datafile_I2$spring_precip_t),
                        pop_size = seq(min(datafile_I2$pop_size, na.rm=T), max(datafile_I2$pop_size,na.rm=T), length=30), 
                        Round = factor(1,levels=c(1:9), labels=c("B", "C", "CP", "EX", "MP", "O", "P", "SW", "W")))
prediction5 <- predict(final_inheritance_PO, new_data6, type="response", re.form=NA)
predictions_ucl <- rep(NA,30) 
for(i in 1:30){predictions_ucl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,2],slopes[3:5,1]), 
                                                   xs=c(xs[1],new_data6$pop_size[i], xs[3:4]), factors=0)}
predictions_lcl <- rep(NA,30)
for(i in 1:30){predictions_lcl[i] <- predict_devel(slopes=c(slopes[1,1],slopes[2,3],slopes[3:5,1]), 
                                                   xs=c(xs[1],new_data6$pop_size[i], xs[3:4]), factors=0)}

plot(April_hatch ~ pop_size, data=datafile_I2, ylab="", xlab="Population size", 
     pch=16, col="grey", axes=F, cex.lab=1.5, type='n')
axis(1)
axis(2, las=1)
title(ylab="Scaled hatch date", line=2.5, cex.lab=1.5)
points(datafile_I2$pop_size, datafile_I2$April_hatch, cex=datafile_I2$coder/20, pch=16, col="grey")
polygon(x=c(new_data6$pop_size, rev(new_data6$pop_size)), y=c(predictions_ucl, rev(predictions_lcl)), border=F, 
        col = rgb(1, 0, 0,0.25))
lines(new_data6$pop_size, prediction5, col="red", lwd=2)

dev.off()

######################################################################


### FIGURE SET UP IPM-2 ####
# Plot of the predictions to the end of 21st Century
# load data
load("low_results_decoupled_SEPT.RData")
load("mid_results_decoupled_SEPT.RData")
load("high_results_decoupled_SEPT.RData")

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")

# now need to try plots
extinction_checks <- as.data.frame(matrix(NA,nrow=3000, ncol=5))
colnames(extinction_checks) <- c("max_sync", "ex_prob", "max_temp", "actual_sync", "year")

# SYNCHRONY
for(i in 1:1000){
  low_results[[i]]$SYNC <- ((low_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((low_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  low_results[[i]]$Hatch_mean-low_results[[i]]$HF
  mid_results[[i]]$SYNC <- ((mid_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((mid_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  high_results[[i]]$SYNC <- ((high_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((high_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  extinction_checks[i,1] <- max(low_results[[i]]$SYNC)
  extinction_checks[i+1000,1] <- max(mid_results[[i]]$SYNC)
  extinction_checks[i+2000,1] <- max(high_results[[i]]$SYNC)
  extinction_checks[i,3] <- max(low_predictions[[i]]$s_temp_max)
  extinction_checks[i+1000,3] <- max(mid_predictions[[i]]$s_temp_max)
  extinction_checks[i+2000,3] <- max(high_predictions[[i]]$s_temp_max)
  if(length(which(low_results[[i]]$PS<1))>0){extinction_checks[i,2] <- 1
  extinction_checks[i,4] <- low_results[[i]]$SYNC[which(low_results[[i]]$PS<1)[1]]
  extinction_checks[i,5] <- which(low_results[[i]]$PS<1)[1]}
  if(length(which(mid_results[[i]]$PS<1))>0){extinction_checks[i+1000,2] <- 1
  extinction_checks[i+1000,4] <- mid_results[[i]]$SYNC[which(mid_results[[i]]$PS<1)[1]]
  extinction_checks[i+1000,5] <- which(mid_results[[i]]$PS<1)[1]}
  if(length(which(high_results[[i]]$PS<1))>0){extinction_checks[i+2000,2] <- 1
  extinction_checks[i+2000,4] <- high_results[[i]]$SYNC[which(high_results[[i]]$PS<1)[1]]
  extinction_checks[i+2000,5] <- which(high_results[[i]]$PS<1)[1]}
}

### FIGURE 1 IPM-2 ####

## LOW EMISSIONS SYNCRHONY

low_plotting <- rbindlist(low_results, idcol=T)
low_plotting$Year <- rep(2010:2100,1000)

low_synchrony <- ggplot(low_plotting,aes(Year, SYNC))+
  theme_minimal() +
  coord_cartesian(ylim= c(-28,28)) + 
  ylab("Synchrony") +
  ggtitle("Low emissions") +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100) +
  scale_fill_viridis_c(begin=0.1, trans="sqrt")+
  geom_line(data=filter(low_plotting, .id < 11), 
            aes(x=Year, y=SYNC, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_point(data=extinction_checks[1:1000,], 
             #aes(year+2009, actual_sync), 
             #shape=4, colour="orangered", size=0.5)+
  theme(legend.position='none')+  
  expand_limits(y=c(-28,28), x=2101)

## HIGH EMISSIONS SYNCRHONY

high_plotting <- rbindlist(high_results, idcol=T)
high_plotting$Year <- rep(2010:2100,1000)

high_synchrony <- ggplot(high_plotting,aes(Year, SYNC))+
  theme_minimal() +
  coord_cartesian(ylim= c(-28,28)) + 
  ylab("") +
  ggtitle("High emissions") +
  theme(plot.title = element_text(hjust = 1)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100, show.legend = F) +
  scale_fill_viridis_c(begin=0.1, trans="sqrt")+
  geom_line(data=filter(high_plotting, .id < 11), 
            aes(x=Year, y=SYNC, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  #geom_point(data=extinction_checks[2001:3000,], 
             #aes(year+2009, actual_sync, colour="Extinction event"), 
             #shape=4, size=0.5,
             #show.legend = TRUE)+
  #scale_shape_manual(values=4) +
  #scale_colour_manual(values="orangered") +
  theme(legend.position='bottom')+  
  expand_limits(y=c(-28,28), x=2101)

#leg2 <- get_legend(high_synchrony)

high_synchrony <- high_synchrony+theme(legend.position = 'none')


### Now ARRANGE!

Figure1 <- ggarrange(low_synchrony, 
          high_synchrony,
          widths = c(1.2,1.2), nrow = 1)

ggsave(filename = "FIGURE1_IPM2.png", plot = Figure1, device = NULL, path = NULL,
       scale = 1, width = 14, height = 8, units = c("cm"),
       dpi = 300)
x

### FIGURE 2 IPM-2 ####
# Trait change

# GGPLOT CODE

high_plotting$Hatch_mean <- unscale(high_plotting$Hatch_mean, 
                                    descale, "April_hatch")
high_plotting$G_mean <- unscale(high_plotting$G_mean, 
                                    descale, "April_hatch")/2
high_plotting$E_mean <- unscale(high_plotting$E_mean, 
                                    descale, "April_hatch")/2

Figure2 <- ggplot(high_plotting,aes(Year, Hatch_mean))+
  theme_minimal() +
  ylab("Mean hatch date") +
  xlab("Year") +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100, show.legend = FALSE) +
  scale_fill_viridis_c(begin=0.1, option="E", trans="sqrt")+
  geom_line(data=filter(high_plotting, .id < 11), 
            aes(x=Year, y=Hatch_mean, group=.id, colour="a"), 
            size = 0.1, show.legend = TRUE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_line(data=filter(high_plotting), 
            aes(x=Year, y=E_mean, group=.id, colour="b"), 
            size = 0.1, alpha=0.2, show.legend = TRUE)+
  geom_line(data=filter(high_plotting), 
            aes(x=Year, y=G_mean, group=.id, colour="c"), 
            size = 0.1, alpha=0.2, show.legend = TRUE)+
  scale_colour_manual(values=c("black", "grey", "orange"),
                      labels=c("Phenotype", "Plasticity", "Evolution")) +
  theme(legend.position='bottom', legend.title = element_blank(),
        legend.key.width = unit(1,"cm"))+ 
  guides(colour=guide_legend(override.aes = list(size = 2)))+
  expand_limits(y=c(10,55), x=2101)

ggsave(filename = "FIGURE2_IPM2.png", plot = Figure2, device = NULL, path = NULL,
       scale = 1, width = 10, height = 10, units = c("cm"),
       dpi = 300)
x

### FIGURE 3 IPM-2 ####

## EXINCTION PROBABILITY
extinction_checks[which(is.na(extinction_checks[,2])),2] <- 0
(c(sum(extinction_checks[1:1000,2]),
   sum(extinction_checks[1001:2000,2]),
   sum(extinction_checks[2001:3000,2]))/1000)*100

# 2% (low), 10.3% (mid), 29.8% (high)

# NO VA 2.7, 11.3, 32.9

model <- glm(ex_prob ~ max_sync, 
             data = extinction_checks, 
             family='binomial')

# Create a grid of X and Y values to predict
# can then fill whole area
preds = expand.grid(
  max_sync = seq(0, 30, length.out = 100),
  ex_prob = seq(0, 1, length.out = 100)
)
preds$pred_p <- predict(model, preds, type = "response")

pred_line <- predict(model, data.frame(max_sync = seq(0,30,length.out=100)),
                     type="response")
pred_line <- as.data.frame(cbind(pred_line, as.numeric(seq(0,30,length.out=100))))
colnames(pred_line) <- c("pred", "x")

# choose colours
colours <- heat.colors(10)
colfunc<-colorRampPalette(c(colours[10], colours[2]))

Figure3 <- 
  ggplot(preds, aes(x = max_sync, y = ex_prob)) +
  xlim(c(0,30))+
  ylab("Extinction probability")+
  xlab("Maximum annual synchrony")+
  theme_minimal() +
  geom_tile(aes(fill = pred_p)) +
  scale_fill_gradientn(colors=colfunc(5), 
                       name="Extinction probability") +
  expand_limits(x=c(0,30))+
  geom_jitter(data=extinction_checks, aes(x=max_sync,
                                           y=ex_prob), 
               colour="grey50", width=0.01, height=0.01,
               alpha=0.5, pch=16)+
  geom_line(data=pred_line, aes(x=x,y=pred))+
  theme(legend.position = "bottom")

ggsave(filename = "FIGURE3_IPM2.png", plot = Figure3, 
       device = NULL, path = NULL,
        scale = 1, width = 10, height = 10, units = c("cm"),
        dpi = 300)



### FIGURE 4 IPM-2 ####
# DECOUPLED

extinction_checks$PS <- rep(0, 3000)

## LOW EMISSIONS POP SIZE
  
low_popsize <- ggplot(low_plotting,aes(Year, PS))+
  theme_minimal() +
  coord_cartesian(ylim=c(0,450)) + 
  ylab("Population size") +
  ggtitle("Low emissions") +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100) +
  scale_fill_viridis(begin=0.1, option="B", trans="sqrt")+
  geom_line(data=filter(low_plotting, .id < 11), 
            aes(x=Year, y=PS, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_rug(data=extinction_checks[1:1000,], 
           mapping=aes(x=year+2009), color="orangered",
           na.rm=T, show.legend=T, position ="jitter",
           sides="t")+
  #geom_point(data=extinction_checks[1:1000,], 
   #          aes(year+2009, ex_prob-1), 
   #          shape=4, colour="orangered", size=0.5)+
  theme(legend.position='none')+  
  expand_limits(y=c(-1,500), x=2101)+
  geom_segment(aes(x=2010,xend=2020,y=131,yend=131), color="white", size=0.75)+
  geom_segment(aes(x=2010,xend=2020,y=40,yend=40), color="white", 
               size=0.75, linetype=3)+
  geom_segment(aes(x=2010,xend=2020,y=219,yend=219), color="white", 
               size=0.75, linetype=3)

## HIGH EMISSIONS POP SIZE

high_popsize <- ggplot(high_plotting,aes(Year, PS))+
  theme_minimal() +
  coord_cartesian(ylim=c(0,450)) + 
  ylab("") +
  ggtitle("High emissions") +
  theme(plot.title = element_text(hjust = 1)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100, show.legend = F) +
  scale_fill_viridis_c(begin=0.1, option="B", trans="sqrt")+
  geom_line(data=filter(high_plotting, .id < 11), 
            aes(x=Year, y=PS, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  geom_rug(data=extinction_checks[2001:3000,], 
           mapping=aes(x=year+2009, colour="Extinction events"),
           na.rm=T, show.legend=T, position ="jitter",
           sides="t")+
  scale_colour_manual(values="orangered") +
  theme(legend.position='bottom', legend.title = element_blank()) +  
  expand_limits(y=c(-1,500), x=2101)+
  geom_segment(aes(x=2090,xend=2101,y=131,yend=131), color="white", size=0.75)+
  geom_segment(aes(x=2090,xend=2100,y=40,yend=40), color="white", 
               size=0.75, linetype=3)+
  geom_segment(aes(x=2090,xend=2100,y=219,yend=219), color="white", 
               size=0.75, linetype=3)

legp <- get_legend(high_popsize)

high_popsize <- high_popsize+theme(legend.position = 'none')

### SAVE

Figure4 <- ggarrange(low_popsize, 
                     high_popsize, widths= c(1.2,1.2), nrow = 1)

ggsave(filename = "FIGURE4_IPM2.png", plot = Figure4, device = NULL, path = NULL,
       scale = 1, width = 14, height = 8, units = c("cm"),
       dpi = 300)
ggsave(filename = "FIGURE4_IPM2.1.png", plot = legp, device = NULL, path = NULL,
       scale = 1, width = 14, height = 8, units = c("cm"),
       dpi = 300)

### SOM FIGURE 1 IPM-2  perturb####
# perturbation
perturb_results<-read.csv("perturb_results_ExtremeHF_extra_SEPT.csv", header=T)

# finding boundaries
x <- perturb_results[154:204,1]-perturb_results[(154+25),1]
y <- perturb_results[205:255,1]-perturb_results[(205+25),1]

z <- abs(x)-abs(y)
which(z < 0)

which(perturb_results[205:255,1] < 1)
changer[14]

png(filename="SOM1_IPM2.png", units="cm", height = 16, width = 16,
    res=400)

layout(matrix(c(1,2), 2, 1, byrow=T), heights=c(3,1.5))
changer <- seq(-4, 4, length.out=51) 

par(mar=c(4, 5, 2, 2))
plot(perturb_results[1:51,1] ~ changer, type="l", col="red", 
     xlab="", ylab="", axes=F, ylim=c(0,600))
title(ylab="Population size", cex.lab=1.5, line=3.5)
title(xlab="Change in driver", cex.lab=1.5, line=2.5)
axis(1, at=seq(-4,4,length.out = 9), labels=seq(-4,4,length.out = 9), cex.axis=1.5)
axis(2, las=1, cex.axis=1.5)
#abline(v=changer[5], col="black", lty=1)
#abline(v=changer[39], col="black", lty=1)
points(changer, perturb_results[52:102,1], col="blue", lty=3, type='l')
points(changer, perturb_results[103:153,1], col="blue", lty=1, type='l')
points(changer, perturb_results[154:204,1], col="red", lty=3, type='l', lwd=1.5)
points(changer, perturb_results[205:255,1], col="green", type='l')

plot.new()
par(mar=c(1,1,1,1))
legend("bottom", c("Spring temperature", "Spring precipitation", "Winter temperature", "Winter precipitation", "Synchrony"),
       lty=c(1,3,1,3,1), bty='n', col=c("red", "red", "blue", "blue", "green"), cex=0.7, ncol=3, y.intersp = 5)

dev.off()

(changer[35]*sd(descale$Synchrony))+mean(descale$Synchrony)

### FIGURE 2 IPM
### SOM FIGURE 2 IPM-2 SYNC~ST ####
# load data
load("low_results_decoupled_SEPT.RData")
load("mid_results_decoupled_SEPT.RData")
load("high_results_decoupled_SEPT.RData")

# now need to try plots
extinction_checks <- matrix(NA,nrow=3000, ncol=8)

# SYNCHRONY
for(i in 1:1000){
  low_results[[i]]$SYNC <- ((low_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((low_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  low_results[[i]]$Hatch_mean-low_results[[i]]$HF
  mid_results[[i]]$SYNC <- ((mid_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((mid_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  high_results[[i]]$SYNC <- ((high_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((high_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  extinction_checks[i,1] <- max(low_results[[i]]$SYNC)
  extinction_checks[i+1000,1] <- max(mid_results[[i]]$SYNC)
  extinction_checks[i+2000,1] <- max(high_results[[i]]$SYNC)
  extinction_checks[i,4] <- sd(low_results[[i]]$SYNC)
  extinction_checks[i+1000,4] <- sd(mid_results[[i]]$SYNC)
  extinction_checks[i+2000,4] <- sd(high_results[[i]]$SYNC)
  extinction_checks[i,5] <- max(low_predictions[[i]]$s_temp_max)
  extinction_checks[i+1000,5] <- max(mid_predictions[[i]]$s_temp_max)
  extinction_checks[i+2000,5] <- max(high_predictions[[i]]$s_temp_max)
  extinction_checks[i,6] <- min(low_predictions[[i]]$s_precip)
  extinction_checks[i+1000,6] <- min(mid_predictions[[i]]$s_precip)
  extinction_checks[i+2000,6] <- min(high_predictions[[i]]$s_precip)
  extinction_checks[i,7] <- min(low_predictions[[i]]$w_precip)
  extinction_checks[i+1000,7] <- min(mid_predictions[[i]]$w_precip)
  extinction_checks[i+2000,7] <- min(high_predictions[[i]]$w_precip)
  extinction_checks[i,8] <- min(low_predictions[[i]]$w_temp)
  extinction_checks[i+1000,8] <- min(mid_predictions[[i]]$w_temp)
  extinction_checks[i+2000,8] <- min(high_predictions[[i]]$w_temp)
  if(length(which(low_results[[i]]$PS<1))>0){extinction_checks[i,2] <- 1}
  if(length(which(mid_results[[i]]$PS<1))>0){extinction_checks[i+1000,2] <- 1}
  if(length(which(high_results[[i]]$PS<1))>0){extinction_checks[i+2000,2] <- 1}
}

extinction_checks[which(is.na(extinction_checks[,2])),2] <- 0
(c(sum(extinction_checks[1:1000,2]),
   sum(extinction_checks[1001:2000,2]),
   sum(extinction_checks[2001:3000,2]))/1000)*100
extinction_checks <- as.data.frame(extinction_checks)

modelST <- lm(max_sync ~ max_temp, data = extinction_checks)
summary(modelST)

### Keep this code incase

plot(max_sync ~ max_temp, data = extinction_checks,
     pch=16, xlab="Maximum spring temperature", 
     ylab="Maximum positive synchrony", 
     axes=F)
axis(1)
axis(2, at=seq(8,28,2))
abline(modelST, col=2)

### But plot this for now
model <- glm(ex_prob ~ max_sync, 
             data = extinction_checks, 
             family='binomial')

# Create a grid of X and Y values to predict
# can then fill whole area
preds = expand.grid(
  max_sync = seq(-28, 28, length.out = 100)
)
preds$pred_l <- predict(model, preds, type = "response")

png(filename="SOM2_IPM2.png", units="cm", height = 10, width = 10,
    res=400)

plot(jitter(extinction_checks$ex_prob, 0.15) ~ extinction_checks$max_sync, 
     axes=F, ylab="Extinction probability",
     xlab="Maximum positive synchrony", pch=16, col=alpha("black",0.2))
points(preds$max_sync, preds$pred_l, type='l', col=2)
axis(1, at = seq(8,28,2))
axis(2, las=1)


dev.off()

### SOM FIGURE 3 IPM-2 weather held ####

# load data
load("only_ST_results_SEPT.RData")
load("only_SP_results_SEPT.RData")
load("only_WP_results_SEPT.RData")
load("only_WT_results_SEPT.RData")

# SYNCHRONY

# now need to try plots
extinction_checks <- matrix(NA, nrow=4000, ncol=1)

for(i in 1:1000){
  no_ST_results[[i]]$SYNC <- ((no_ST_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((no_ST_results[[i]]$HF*sd(descale$April_hatch))+mean(descale$April_hatch))
  no_ST_results[[i]]$Hatch_mean-no_ST_results[[i]]$HF
  no_SP_results[[i]]$SYNC <- ((no_SP_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((no_SP_results[[i]]$HF*sd(descale$April_hatch))+mean(descale$April_hatch))
  no_SP_results[[i]]$Hatch_mean-no_SP_results[[i]]$HF
  no_WT_results[[i]]$SYNC <- ((no_WT_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((no_WT_results[[i]]$HF*sd(descale$April_hatch))+mean(descale$April_hatch))
  no_WT_results[[i]]$Hatch_mean-no_WT_results[[i]]$HF
  no_WP_results[[i]]$SYNC <- ((no_WP_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((no_WP_results[[i]]$HF*sd(descale$April_hatch))+mean(descale$April_hatch))
  no_WP_results[[i]]$Hatch_mean-no_WP_results[[i]]$HF
  if(length(which(no_ST_results[[i]]$PS<1))>0){extinction_checks[i,1] <- 1}
  if(length(which(no_SP_results[[i]]$PS<1))>0){extinction_checks[i+1000,1] <- 1}
  if(length(which(no_WT_results[[i]]$PS<1))>0){extinction_checks[i+2000,1] <- 1}
  if(length(which(no_WP_results[[i]]$PS<1))>0){extinction_checks[i+3000,1] <- 1}
}

sum(extinction_checks[1:1000,1], na.rm=T) # 0
sum(extinction_checks[1001:2000,1], na.rm=T) # 508
sum(extinction_checks[2001:3000,1], na.rm=T) # 329
sum(extinction_checks[3001:4000,1], na.rm=T) # 361

# DECOUPLED
# line for sync is mean of data
png(filename="SOM3_IPM2.png", units="cm", height = 16, width = 16,
    res=400)

layout(matrix(c(1,2,3,4,5,6,7,8),2,4, byrow=T), heights=c(2.8,2.9), widths=c(3.7,2.8, 2.8,3))

par(mar=c(2,5,2,0))
# Population size
plotter(no_ST_results, "Population size", "", "held ST", x.axis=F, y.axis=T, adj=0.5,
        ylim=c(0,600))
par(mar=c(2,1,2,1))
plotter(no_SP_results, "", "", "held SP", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))
par(mar=c(2,1,2,1))
plotter(no_WT_results, "", "", "held WT", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))
par(mar=c(2,0,2,2.5))
plotter(no_WP_results, "", "", "held WP", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))


par(mar=c(4,5,0,0))
# Hatch
plotter3(no_ST_results, "Hatch date", "", "", descale, x.axis=T, y.axis=T)
par(mar=c(4,1,0,1))
plotter3(no_SP_results, "", "Time", "", descale, x.axis=T, y.axis=F, adj2=1)
par(mar=c(4,1,0,1))
plotter3(no_WT_results, "", "(years)", "", descale, x.axis=T, y.axis=F, adj2=0)
par(mar=c(4,0,0,2.5))
plotter3(no_WP_results, "", "", "", descale, x.axis=T, y.axis=F)

dev.off()


### SOM FIGURE 4-6 IPM-2 no evo ####
# load data
load("low_results_decoupled_VA0.RData")
load("mid_results_decoupled_VA0.RData")
load("high_results_decoupled_VA0.RData")

# now need to try plots
extinction_checks <- as.data.frame(matrix(NA,nrow=3000, ncol=5))
colnames(extinction_checks) <- c("max_sync", "ex_prob", "max_temp", "actual_sync", "year")

# SYNCHRONY
for(i in 1:1000){
  low_results[[i]]$SYNC <- ((low_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((low_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  low_results[[i]]$Hatch_mean-low_results[[i]]$HF
  mid_results[[i]]$SYNC <- ((mid_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((mid_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  high_results[[i]]$SYNC <- ((high_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((high_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  extinction_checks[i,1] <- max(low_results[[i]]$SYNC)
  extinction_checks[i+1000,1] <- max(mid_results[[i]]$SYNC)
  extinction_checks[i+2000,1] <- max(high_results[[i]]$SYNC)
  extinction_checks[i,3] <- max(low_predictions[[i]]$s_temp_max)
  extinction_checks[i+1000,3] <- max(mid_predictions[[i]]$s_temp_max)
  extinction_checks[i+2000,3] <- max(high_predictions[[i]]$s_temp_max)
  if(length(which(low_results[[i]]$PS<1))>0){extinction_checks[i,2] <- 1
  extinction_checks[i,4] <- low_results[[i]]$SYNC[which(low_results[[i]]$PS<1)[1]]
  extinction_checks[i,5] <- which(low_results[[i]]$PS<1)[1]}
  if(length(which(mid_results[[i]]$PS<1))>0){extinction_checks[i+1000,2] <- 1
  extinction_checks[i+1000,4] <- mid_results[[i]]$SYNC[which(mid_results[[i]]$PS<1)[1]]
  extinction_checks[i+1000,5] <- which(mid_results[[i]]$PS<1)[1]}
  if(length(which(high_results[[i]]$PS<1))>0){extinction_checks[i+2000,2] <- 1
  extinction_checks[i+2000,4] <- high_results[[i]]$SYNC[which(high_results[[i]]$PS<1)[1]]
  extinction_checks[i+2000,5] <- which(high_results[[i]]$PS<1)[1]}
}

# DECOUPLED
# line for sync is mean of data
# split the pop dynamics and synchrony
png(filename="SOM4_IPM2.png", units="cm", height = 10, width = 16,
    res=400)

#layout(matrix(c(1,2,3),1,3, byrow=T), widths=c(3.7,2.8,3))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(3.7,3))

par(mar=c(4,5,2,0))
# Population size
plotter(low_results, "Population size", "Time (years)", "Low", x.axis=T, y.axis=T,  adj=0.5)
#par(mar=c(4,1,2,1))
#plotter(mid_results, "", "Time (years)", "Mid", x.axis=T, y.axis=F,  adj=0.5)
par(mar=c(4,0,2,2.5))
plotter(high_results, "", "Time (years)", "High", x.axis=T,  y.axis=F, adj=0.5)

dev.off()

png(filename="SOM5_IPM2.png", units="cm", height = 10, width = 16,
    res=400)

#layout(matrix(c(1,2,3),1,3, byrow=T), widths=c(3.7,2.8,3))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(3.7,3))

par(mar=c(4,5,2,0))
# Synchrony
plotter4(low_results, "Synchrony", "Time (Years)", "Low", descale, x.axis=T, y.axis=T, adj=0.5)
#par(mar=c(4,1,2,1))
#plotter4(mid_results, "", "Time (Years)", "Mid", descale, x.axis=T, y.axis=F, adj=0.5)
par(mar=c(4,0,2,2.5))
plotter4(high_results, "", "Time (Years)", "High", descale, x.axis=T, y.axis=F, adj=0.5)

dev.off()

### Trait change
png(filename="SOM6_IPM2.png", units="cm", height = 10, width = 15,
    res=400)

#layout(matrix(c(1,2,3,4,5,6),2,3, byrow=T), heights=c(2.8,2.9), widths=c(3.7,2.8,3))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(3,3))

par(mar=c(4,5,2,2))
# phenotypic change
plotter3(high_results, "Hatch date (days)", "Time (years)", "Z", descale, x.axis=T, y.axis=T, adj=0.5)
#par(mar=c(2,1,2,1))
#plotter3(mid_results, "", "", "Mid", descale, x.axis=F, y.axis=F, adj=0.5)
#par(mar=c(2,0,2,2.5))
#plotter3(high_results, "", "", "High", descale, x.axis=F, y.axis=F, adj=0.5)

par(mar=c(4,2,2,5))
# geno evo change
plotter2(high_results, "", "Time (years)", "G and E", descale, x.axis=T, y.axis=T)
#par(mar=c(4,1,0,1))
#plotter2(mid_results, "", "Time (years)", "", descale, x.axis=T, y.axis=F)
#par(mar=c(4,0,0,2.5))
#plotter2(high_results, "", "", "", descale, x.axis=T, y.axis=F)


dev.off()

### SOM FIGURE 7 IPM-2 change points ####

### Change point checking
library(changepoint)
library(chngpt)
library(strucchange)

## EVOLUTION

# Try looping through
changes <- rep(NA, 1000)
changes2 <- rep(NA, 1000)

for(i in 1:1000){
  test <- chngpt.test(formula.null=y~1, formula.chngpt=~x, 
                      dat= data.frame(y = high_results[[i]]$G_mean[6:91],
                                      x = seq(2015,2100,1)), 
                      type="segmented", family='gaussian', chngpts.cnt = 80,
                      mc.n=10)
  changes[i] <- test$chngpt
  if(test$p.value > 0.05){changes[i] <- NA}
  
  test2 <- breakpoints(high_results[[i]]$G_mean[6:91]~
                         seq(2015,2100,1), breaks=1)
  changes2[i] <- test2$breakpoints
}

changes_total <- data.frame(change = c(changes,
                                       changes2+2015),
                            index = c(rep(1,1000),
                                      rep(0,1000)))
plot(changes_total$change, col=changes_total$index+1)
summary(changes_total)
length(which(is.na(changes)==T))


## TESTING RATES
rate_store <- rep(NA, 1000)

for(i in 1:1000){
  rates <- high_results[[i]]$G_mean[6:91]-high_results[[i]]$G_mean[5:90]
  model <- lm(rates~c(1:86))
  if(coef(model)[2]<0 & summary(model)$coefficients[8]<0.05){rate_store[i] <- T}else
  {rate_store[i] <- F}
}

length(which(rate_store==T)) # 771
length(which(rate_store[which(extinction_checks[2001:3000,2]==1)]==T)) # 288


## PLASTICITY

# Little change in plasticity. No change in rate of plasticity change. 

# Try looping through and check variance
changes_E <- matrix(NA, nrow = 1000, ncol=2)

for(i in 1:1000){
  if(is.na(changes[i])==T){changes_E[i,] <- c(NA, NA)}else{
    temp1 <- high_results[[i]]$E_mean[6:(changes[i]-2010)]
    temp2 <- high_results[[i]]$E_mean[(changes[i]-2009):91]
    changes_E[i,] <- c(var(temp1), var(temp2))}
}

plot(c(changes_E[,1], changes_E[,2]), c(rep(1,1000), rep(2,1000)))

summary(changes_E[,1])
summary(changes_E[,2])

# make into dataframe for anova
change_E_dataframe <- data.frame(changes = c(changes_E[,1], changes_E[,2]),
                                 group = c(rep(1,1000), rep(2,1000)),
                                 pair = rep(seq(1,1000,1),2))
summary(lm(changes ~ as.factor(group)+as.factor(pair), 
           data=change_E_dataframe))


png(filename="SOM7_IPM2.png", units="cm", height = 10, width = 16,
    res=400)

layout(matrix(c(1,2),1,2, byrow=T), widths=c(3.7,3.7))

par(mar=c(4,4,2,0))

## Evolutionary change

hist(changes, main="",
     ylab="Frequency", las=1,
     xlab="Year", xlim=c(2020,2100))
abline(v=summary(changes)[3], lwd=2, col=2)
abline(v=summary(changes)[2], lwd=1, col=2, lty=2)
abline(v=summary(changes)[5], lwd=1, col=2, lty=2)
length(which(changes==2071))
title(main='a. Evolution', adj=0, cex=0.5)
# between 2069 to 2073 

## Plasticity change
par(mar=c(4,5,2,1))

plot(changes ~ jitter(group, 0.03), data=change_E_dataframe,
     pch=16, col=group, axes=F,
     ylab="Variance in mean", 
     xlab="Time relative to change point",
     ylim=c(0,0.8),
     xlim=c(0.8,2.2))
axis(1, at=c(1,2), labels=c("Prior", "Post"))
axis(2, las=1, at=seq(0, 0.8, 0.2))
lines(x=c(0.9,1.4),y=c(summary(changes_E[,1])[3], 
                       summary(changes_E[,1])[3]))
lines(x=c(0.9,1.4),y=c(summary(changes_E[,1])[2], 
                       summary(changes_E[,1])[2]), lty=2)
lines(x=c(0.9,1.4),y=c(summary(changes_E[,1])[5], 
                       summary(changes_E[,1])[5]), lty=2)

lines(x=c(1.6,2.1),y=c(summary(changes_E[,2])[3], 
                       summary(changes_E[,2])[3]), 
      col=2)
lines(x=c(1.6,2.1),y=c(summary(changes_E[,2])[2], 
                       summary(changes_E[,2])[2]), 
      col=2, lty=2)
lines(x=c(1.6,2.1),y=c(summary(changes_E[,2])[5], 
                       summary(changes_E[,2])[5]), 
      col=2, lty=2)
title(main='b. Plasticity', adj=0, cex=0.5)

dev.off()

#### Plot climate predictions IPM-2 ####

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")

png(filename="Climate_predictions.png", units="cm", height = 16, width = 16,
    res=400)
layout(matrix(c(1,2,3,4), 2,2, byrow=T), heights=c(5,5))

plot(low_predictions[[1]]$s_temp_max~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="Spring temp", ylim=c(5,25))
axis(1)
axis(2, las=1)
for(i in 1:1000){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_temp_max, type='l', col=alpha("red",0.01))}
abline(h=mean(descale$Spring_temp))

plot(low_predictions[[1]]$s_precip~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="Spring precip", ylim=c(0,400))
axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$s_precip, type='l', lty=2, col=alpha("red",0.01))}
abline(h=mean(descale$spring_precip_t))

plot(low_predictions[[1]]$w_temp~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="Winter temp", ylim=c(0,15))
axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$w_temp, type='l', col=alpha("blue",0.01))}
abline(h=mean(descale$winter_temp))

plot(low_predictions[[1]]$w_precip~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="Winter precip", ylim=c(0,400))
axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$w_precip, type='l', lty=2, col=alpha("blue",0.01))}
abline(h=mean(descale$winter_precip_t))

dev.off()

#### Summaries IPM-2 ####

# load data
load("low_results_decoupled_SEPT.RData")
load("mid_results_decoupled_SEPT.RData")
load("high_results_decoupled_SEPT.RData")

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")


output_predictions_L <- output_predictions_ps(low_results)
output_predictions_L[c(50,950),]
output_predictions_M <- output_predictions_ps(mid_results)
output_predictions_M[c(50,950),]
output_predictions_H <- output_predictions_ps(high_results)
output_predictions_H[c(50,950),]

output_predictions_zL <- output_predictions_z(low_results, descale)
output_predictions_zL[c(50,950),]
output_predictions_zM <- output_predictions_z(mid_results, descale)
output_predictions_zM[c(50,950),]
output_predictions_zH <- output_predictions_z(high_results, descale)
output_predictions_zH[c(50,950),]

# summary

# amount of change
G_changeL <- rep(NA,3000)
E_changeL <- rep(NA,3000)
Z_change <- rep(NA, 1000)
x <- 1:91
for(i in 1:1000){
  G_changeL[i] <- coef(lm(((low_results[[i]]$G_mean*sd(descale$April_hatch))+
                             mean(descale$April_hatch))~x))[2]
  E_changeL[i] <- coef(lm(((low_results[[i]]$E_mean*sd(descale$April_hatch))+
                             mean(descale$April_hatch))~x))[2]
  G_changeL[i+1000] <- coef(lm(((mid_results[[i]]$G_mean*sd(descale$April_hatch))+
                                  mean(descale$April_hatch))~x))[2]
  E_changeL[i+1000] <- coef(lm(((mid_results[[i]]$E_mean*sd(descale$April_hatch))+
                                  mean(descale$April_hatch))~x))[2]
  G_changeL[i+2000] <- coef(lm(((high_results[[i]]$G_mean*sd(descale$April_hatch))+
                                  mean(descale$April_hatch))~x))[2]
  E_changeL[i+2000] <- coef(lm(((high_results[[i]]$E_mean*sd(descale$April_hatch))+
                                  mean(descale$April_hatch))~x))[2]
  Z_change[i] <- coef(lm(((high_results[[i]]$Hatch_mean*sd(descale$April_hatch))+
                            mean(descale$April_hatch))~x))[2]
}

summary(G_changeL) # -0.05 days per year to 0.002 days per year (-4.55 days max, 0.182 min, across 91 years)
summary(E_changeL) # -0.24 days per year to +0.09 days per year (-21.84 days max, +8.19 min)
length(which(E_changeL[1:1000] < 0)) # 816
length(which(E_changeL[1001:2000] < 0)) # 930
length(which(E_changeL[2001:3000] < 0)) # 972
length(which(G_changeL[1:1000] < 0)) # 902
length(which(G_changeL[1001:2000] < 0)) # 965
length(which(G_changeL[2001:3000] < 0)) # 985
sort(G_changeL[1:1000])[c(225,975)] # -0.0067, -0.0008
sort(G_changeL[1001:2000])[c(225,975)] # -0.0108, -0.0004
sort(G_changeL[2001:3000])[c(225,975)] # -0.0199, -0.0003

sort(E_changeL[1:2000])[c(225,975)] # -0.0867, 0.0457
sort(E_changeL[1001:2000])[c(225,975)] # -0.1177, 0.0236
sort(E_changeL[2001:3000])[c(225,975)] # -0.1498, 0.0052

summary(Z_change)

#### Standard deviation of population size

sd_pop <- rep(NA, 3000)

for(i in 1:1000){
  sd_pop[i] <- sd(low_results[[i]]$PS)
  sd_pop[i+1000] <- sd(mid_results[[i]]$PS)
  sd_pop[i+2000] <- sd(high_results[[i]]$PS)
}

sort(sd_pop)[c(75,2925)] # 25 to 81
length(which(sort(sd_pop)>46))

#### SOM TREND TABLE IPM-2 ####

load("low_results_decoupled_SEPT.RData")
load("mid_results_decoupled_SEPT.RData")
load("high_results_decoupled_SEPT.RData")

load("low_predictions_DECOUPLED.RData")
load("mid_predictions_DECOUPLED.RData")
load("high_predictions_DECOUPLED.RData")

for(i in 1:1000){
  low_results[[i]]$SYNC <- ((low_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((low_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  low_results[[i]]$Hatch_mean-low_results[[i]]$HF
  mid_results[[i]]$SYNC <- ((mid_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((mid_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))
  high_results[[i]]$SYNC <- ((high_results[[i]]$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)) - 
    ((high_results[[i]]$HF*sd(descale$Half_fall))+mean(descale$Half_fall))}

trends <- data.frame(Hatch = rep(NA,3000),
                     G = rep(NA,3000),
                     E = rep(NA,3000),
                     PS = rep(NA,3000),
                     Sync = rep(NA,3000))

# loop to calculate trends

for(i in 1:1000){
  trends$Hatch[i] <- coef(lm(unscale(low_results[[i]]$Hatch_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$G[i] <- coef(lm(unscale(low_results[[i]]$G_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$E[i] <- coef(lm(unscale(low_results[[i]]$E_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$PS[i] <- coef(lm(unscale(low_results[[i]]$PS, descale, "April_hatch")~c(1:91)))[2]
  trends$Hatch[i+1000] <- coef(lm(unscale(mid_results[[i]]$Hatch_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$G[i+1000] <- coef(lm(unscale(mid_results[[i]]$G_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$E[i+1000] <- coef(lm(unscale(mid_results[[i]]$E_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$PS[i+1000] <- coef(lm(unscale(mid_results[[i]]$PS, descale, "April_hatch")~c(1:91)))[2]
  trends$Hatch[i+2000] <- coef(lm(unscale(high_results[[i]]$Hatch_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$G[i+2000] <- coef(lm(unscale(high_results[[i]]$G_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$E[i+2000] <- coef(lm(unscale(high_results[[i]]$E_mean, descale, "April_hatch")~c(1:91)))[2]
  trends$PS[i+2000] <- coef(lm(unscale(high_results[[i]]$PS, descale, "April_hatch")~c(1:91)))[2]
}

# summarising results, how many increase, decrease or < 1% change?

#Hatch
length(which(trends$Hatch[1:1000] < 0)) # 828/1000 decrease
length(which(trends$Hatch[1001:2000] < 0)) # 916/1000 decrease
length(which(trends$Hatch[2001:3000] < 0)) # 983/1000 decrease

#E
length(which(trends$E[1:1000] < 0)) # 106/1000 decrease
length(which(trends$E[1001:2000] < 0)) # 296/1000 decrease
length(which(trends$E[2001:3000] < 0)) # 517/1000 decrease

#G
length(which(trends$G[1:1000] < 0)) # 1000 decrease
length(which(trends$G[1001:2000] < 0)) # 1000 decrease
length(which(trends$G[2001:3000] < 0)) # 1000 decrease

#PS
length(which(trends$PS[1:1000] < 0)) # 42 decrease
length(which(trends$PS[1001:2000] < 0)) # 55 decrease
length(which(trends$PS[2001:3000] < 0)) # 179 decrease

quantile(sort(trends$PS[which(trends$PS[1:1000] > 0)]), c(0.025,0.975)) # 0.45 to 3.7
quantile(sort(trends$PS[which(trends$PS[1001:2000] > 0)+1000]), c(0.025,0.975)) # 0.339 to 4.249
quantile(sort(trends$PS[which(trends$PS[2001:3000] > 0)+2000]), c(0.025,0.975)) # 0.1995 to 4.200


### standard deviation

standdev <- data.frame(Hatch = rep(NA,3000),
                       G = rep(NA,3000),
                       E = rep(NA,3000),
                       PS = rep(NA,3000),
                       Sync = rep(NA,3000))

# loop to calculate trends

for(i in 1:1000){
  standdev$Hatch[i] <- sd(unscale(low_results[[i]]$Hatch_mean, descale, "April_hatch"))
  standdev$G[i] <- sd(unscale(low_results[[i]]$G_mean, descale, "April_hatch"))
  standdev$E[i] <- sd(unscale(low_results[[i]]$E_mean, descale, "April_hatch"))
  standdev$PS[i] <- sd(low_results[[i]]$PS)
  standdev$SYNC[i] <- sd(low_results[[i]]$SYNC)
  standdev$Hatch[i+1000] <- sd(unscale(mid_results[[i]]$Hatch_mean, descale, "April_hatch"))
  standdev$G[i+1000] <- sd(unscale(mid_results[[i]]$G_mean, descale, "April_hatch"))
  standdev$E[i+1000] <- sd(unscale(mid_results[[i]]$E_mean, descale, "April_hatch"))
  standdev$PS[i+1000] <- sd(mid_results[[i]]$PS)
  standdev$SYNC[i+1000] <- sd(mid_results[[i]]$SYNC)
  standdev$Hatch[i+2000] <- sd(unscale(high_results[[i]]$Hatch_mean, descale, "April_hatch"))
  standdev$G[i+2000] <- sd(unscale(high_results[[i]]$G_mean, descale, "April_hatch"))
  standdev$E[i+2000] <- sd(unscale(high_results[[i]]$E_mean, descale, "April_hatch"))
  standdev$PS[i+2000] <- sd(high_results[[i]]$PS)
  standdev$SYNC[i+2000] <- sd(high_results[[i]]$SYNC)
}

quantile(sort(standdev$PS), c(0.025,0.975)) # 66.676 to 127.266
quantile(sort(standdev$Hatch), c(0.025,0.975)) # 2.14 to 5.97
quantile(sort(standdev$SYNC), c(0.025,0.975)) # 4.6 to 9.19

## data
summary(tapply(descale$R_pop_size, descale$Year, mean)) # 45.88
sd(tapply(descale$R_pop_size, descale$Year, mean)) # 45.88

### Summary table
# var, scenario, trend lower, trend upper, sd lower, sd upper
output <- data.frame(Variable = rep(c("Z", "G", "E", "PS"), 3),
                     Scenario = rep(c("low", "mid", "high"), each = 4),
                     trend_lower = c(quantile(sort(trends$Hatch[1:1000]), 0.025),
                                     quantile(sort(trends$G[1:1000]), 0.025), 
                                     quantile(sort(trends$E[1:1000]), 0.025),
                                     quantile(sort(trends$PS[1:1000]), 0.025),
                                     quantile(sort(trends$Hatch[1001:2000]), 0.025),
                                     quantile(sort(trends$G[1001:2000]), 0.025), 
                                     quantile(sort(trends$E[1001:2000]), 0.025),
                                     quantile(sort(trends$PS[1001:2000]), 0.025),
                                     quantile(sort(trends$Hatch[2001:3000]), 0.025),
                                     quantile(sort(trends$G[2001:3000]), 0.025), 
                                     quantile(sort(trends$E[2001:3000]), 0.025),
                                     quantile(sort(trends$PS[2001:3000]), 0.025)),
                     trend_upper = c(quantile(sort(trends$Hatch[1:1000]), 0.975),
                                     quantile(sort(trends$G[1:1000]), 0.975), 
                                     quantile(sort(trends$E[1:1000]), 0.975),
                                     quantile(sort(trends$PS[1:1000]), 0.975),
                                     quantile(sort(trends$Hatch[1001:2000]), 0.975),
                                     quantile(sort(trends$G[1001:2000]), 0.975), 
                                     quantile(sort(trends$E[1001:2000]), 0.975),
                                     quantile(sort(trends$PS[1001:2000]), 0.975),
                                     quantile(sort(trends$Hatch[2001:3000]), 0.975),
                                     quantile(sort(trends$G[2001:3000]), 0.975), 
                                     quantile(sort(trends$E[2001:3000]), 0.975),
                                     quantile(sort(trends$PS[2001:3000]), 0.975)),
                     SD_lower = c(quantile(sort(standdev$Hatch[1:1000]), 0.025),
                                  quantile(sort(standdev$G[1:1000]), 0.025), 
                                  quantile(sort(standdev$E[1:1000]), 0.025),
                                  quantile(sort(standdev$PS[1:1000]), 0.025),
                                  quantile(sort(standdev$Hatch[1001:2000]), 0.025),
                                  quantile(sort(standdev$G[1001:2000]), 0.025), 
                                  quantile(sort(standdev$E[1001:2000]), 0.025),
                                  quantile(sort(standdev$PS[1001:2000]), 0.025),
                                  quantile(sort(standdev$Hatch[2001:3000]), 0.025),
                                  quantile(sort(standdev$G[2001:3000]), 0.025), 
                                  quantile(sort(standdev$E[2001:3000]), 0.025),
                                  quantile(sort(standdev$PS[2001:3000]), 0.025)),
                     SD_upper = c(quantile(sort(standdev$Hatch[1:1000]), 0.975),
                                  quantile(sort(standdev$G[1:1000]), 0.975), 
                                  quantile(sort(standdev$E[1:1000]), 0.975),
                                  quantile(sort(standdev$PS[1:1000]), 0.975),
                                  quantile(sort(standdev$Hatch[1001:2000]), 0.975),
                                  quantile(sort(standdev$G[1001:2000]), 0.975), 
                                  quantile(sort(standdev$E[1001:2000]), 0.975),
                                  quantile(sort(standdev$PS[1001:2000]), 0.975),
                                  quantile(sort(standdev$Hatch[2001:3000]), 0.975),
                                  quantile(sort(standdev$G[2001:3000]), 0.975), 
                                  quantile(sort(standdev$E[2001:3000]), 0.975),
                                  quantile(sort(standdev$PS[2001:3000]), 0.975)))

write.csv(output, "Trends_IPM2,csv", row.names=F)

