#### FIGURE 1 ####

#### Packages and function scripts ####
source('Reformat_data_for_plotting.R')

#### IMPORT AND FORMAT DATA ####
datafile <- read.csv("bio_data.csv", header=T) # scaled
descale <- read.csv("descale.csv", header=T) # unscaled
cross_val <- read.csv("cross_validation_data.csv", header=T)

cross_val <- read.csv("cross_validation_data.csv", header=T)
load("cross_val_RESULT_RO.RData")
cross_val_RESULT_RO1 <- list(cross_val_RESULT_RO[[1]][[1]])
for(i in 2:10){cross_val_RESULT_RO1 <- 
  c(cross_val_RESULT_RO1, list(cross_val_RESULT_RO[[i]][[1]]))}
cross_val_RESULT_RO2 <- for_plotting_CV(cross_val_RESULT_RO1, list_K_fold, descale)

load("cross_val_RESULT_TP.RData")
cross_val_RESULT_TP1 <- list(cross_val_RESULT_TP[[1]][[1]])
for(i in 2:10){cross_val_RESULT_TP1 <- 
  c(cross_val_RESULT_TP1, list(cross_val_RESULT_TP[[i]][[1]]))}
cross_val_RESULT_TP2 <- for_plotting_CV(cross_val_RESULT_TP1, list_K_fold, descale)

load("cross_val_RESULT_STANDARD.RData")
cross_val_RESULT_PO1 <- list(cross_val_RESULT_po[[1]][[1]])
for(i in 2:10){cross_val_RESULT_PO1 <- 
  c(cross_val_RESULT_PO1, list(cross_val_RESULT_po[[i]][[1]]))}
cross_val_RESULT_PO2 <- for_plotting_CV(cross_val_RESULT_PO1, list_K_fold, descale)

#### PLOT FIGURE ####

layout(matrix(c(1,2,3,4,5,5), 3,2, byrow=T), heights=c(4,4,2))
par(mar=c(2,4.5,2,1))
# Population size
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "Population size", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, cex.lab=1.5,
     main="Quantitative genetic IPM")
title(main='a.', adj=0, cex.main=1)
axis(2, las=1, cex.axis=1.25)

points(PS ~ Year, data = cross_val_RESULT_TP2,
       type = 'l', col = 'black', lty = 1)
par(mar=c(2,1,2,4.5))
plot(R_only_PS ~ Year, data = cross_val[2:51,], type = 'l',
     xlab = "", ylab = "", axes=F,
     ylim=c(0,350), col = "grey70", lwd=2, main="Standard IPM")
title(main='c.', adj=0, cex.main=1)
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
points(Hatch_mean ~ Year, data = cross_val_RESULT_PO2,
       type = 'l', col = 'black', lty = 4)

# Add legend
plot.new()
legend("top", c("Model projections", "Observations"), 
       bty="n", col = c("black",  "grey"), lty=c(1,1),
       cex=1.5, ncol=2, lwd=c(2,2))

#### CALCULATE MAE ####

MAE <- data.frame(
  ps.TP = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_TP2$PS)),2),
  ps.po = round(mean(abs(cross_val$R_only_PS[2:51]-cross_val_RESULT_PO2$PS)),2),
  z.TP = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_TP2$Hatch_mean)),2),
  z.po = round(mean(abs(cross_val$R_Mean_hatch[2:51]-cross_val_RESULT_PO2$Hatch_mean)),2)
)

colnames(MAE) <- rep(c( "Quantitative genetic IPM", "Standard IPM"), 2)
MAE