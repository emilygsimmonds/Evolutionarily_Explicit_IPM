#### FIGURE 2 ####

#### Packages and function scripts ####
source('Reformat_data_for_plotting.R')

#### IMPORT AND FORMAT DATA ####
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled
cross_val <- read.csv("cross_validation_data.csv", header=T)

load("cross_val_RESULT_50.RData")
cross_val_RESULT_TP2 <- cross_val_RESULT_TP[[1]]
cross_val_RESULT_TP2$Year <- 1961:2010
cross_val_RESULT_TP2$Hatch_mean <- (cross_val_RESULT_TP2$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)

load("cross_val_RESULT_STANDARD_50.RData")
cross_val_RESULT_PO2 <- cross_val_RESULT_po[[1]]
cross_val_RESULT_PO2$Year <- 1961:2010
cross_val_RESULT_PO2$Hatch_mean <- (cross_val_RESULT_PO2$Hatch_mean*sd(descale$April_hatch))+mean(descale$April_hatch)

#### PLOT FIGURE ####

layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))

par(mar=c(2,4.5,2,1))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, cex.lab=1.5,
     main = "Quantitative genetic IPM")
title(main='a.', adj=0, cex.main=1)
axis(2, las=1, cex.axis=1.25)
points(PS[2:51] ~ Year[1:50], data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2,
     main = "Standard IPM")
title(main='c.', adj=0, cex.main=1)
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

#### CALCULATE MAE ####

MAE <- data.frame(
  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP2$PS)),2),
  ps.po = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_PO2$PS)),2),
  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP2$Hatch_mean)),2),
  z.po = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_PO2$Hatch_mean)),2)
)

colnames(MAE) <- rep(c( "Quantitative genetic IPM", "Standard IPM"), 2)
MAE