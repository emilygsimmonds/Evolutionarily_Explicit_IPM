# code to run all figures from Ecology Letters paper
######################################################################

### SET UP ####

library(tidyverse)
library(data.table)
library(viridis)
library(patchwork)

source('Reformat_data_for_plotting.R')
source('Plotter_function.R')

######################################################################

# Plot of the predictions to the end of 21st Century
# load data
load("low_emissions_results.RData")
load("medium_emissions_results.RData")
load("high_emissions_results.RData")

load("low_emission_predictions.RData")
load("mid_emission_predictions.RData")
load("high_emission_predictions.RData")

descale <- read.csv("descale.csv", header=T)

# create data frame to store output
extinction_checks <- as.data.frame(matrix(NA, nrow=3000, ncol=5))
colnames(extinction_checks) <- c("max_sync", "ex_prob", "max_temp", "actual_sync", "year")

# SYNCHRONY: calculate synchrony values for each predicted year on original scale
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

# create plotting datasets
low_plotting <- rbindlist(low_results, idcol=T)
low_plotting$Year <- rep(2010:2100,1000)
low_plotting$SYNC <- low_plotting$SYNC + 13

mid_plotting <- rbindlist(mid_results, idcol=T)
mid_plotting$Year <- rep(2010:2100,1000)
mid_plotting$SYNC <- mid_plotting$SYNC + 13

high_plotting <- rbindlist(high_results, idcol=T)
high_plotting$Year <- rep(2010:2100,1000)
high_plotting$SYNC <- high_plotting$SYNC + 13

######################################################################

### FIGURE 1 ####

## LOW EMISSIONS SYNCRHONY

low_synchrony <- ggplot(low_plotting,aes(Year, SYNC))+
  theme_minimal() +
  coord_cartesian(ylim= c(-15,41)) + 
  ylab("Asynchrony (days)") +
  xlab(NULL)+
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
  theme(legend.position='none')+  
  expand_limits(y=c(-15,41), x=2101)

## MID EMISSIONS SYNCRHONY

mid_synchrony <- ggplot(mid_plotting,aes(Year, SYNC))+
  theme_minimal() +
  coord_cartesian(ylim= c(-15,41)) + 
  ylab(NULL) +
  ggtitle("Mid emissions") +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100) +
  scale_fill_viridis_c(begin=0.1, trans="sqrt")+
  geom_line(data=filter(mid_plotting, .id < 11), 
            aes(x=Year, y=SYNC, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position='none', axis.text.y=element_blank())+  
  expand_limits(y=c(-15,41), x=2101)

## HIGH EMISSIONS SYNCRHONY

high_synchrony <- ggplot(high_plotting,aes(Year, SYNC))+
  theme_minimal() +
  coord_cartesian(ylim= c(-15,41)) + 
  ylab(NULL) +
  xlab(NULL)+
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
  theme(legend.position='bottom')+  
  expand_limits(y=c(-15,41), x=2101)

### Now combine

Figure1 <- low_synchrony + mid_synchrony + high_synchrony

######################################################################

### FIGURE 2 ####
# Trait change

mid_plotting$Hatch_mean <- unscale(mid_plotting$Hatch_mean, 
                                    descale, "April_hatch")
mid_plotting$G_mean <- (unscale(mid_plotting$G_mean, 
                                    descale, "April_hatch")/2)
mid_plotting$G_mean <- mid_plotting$G_mean-mean(mid_plotting$G_mean)
mid_plotting$E_mean <- (unscale(mid_plotting$E_mean, 
                                    descale, "April_hatch")/2)
mid_plotting$E_mean <- mid_plotting$E_mean - mean(mid_plotting$E_mean)

Figure2 <- ggplot(mid_plotting,aes(Year, Hatch_mean))+
  theme_minimal() +
  ylab("Mean hatch date") +
  xlab("Year") +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100, show.legend = FALSE) +
  scale_fill_viridis_c(begin=0.1, option="E", trans="sqrt")+
  geom_line(data=filter(mid_plotting, .id < 11), 
            aes(x=Year, y=Hatch_mean, group=.id, colour="a"), 
            size = 0.1, show.legend = TRUE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_line(data=filter(mid_plotting), 
            aes(x=Year, y=E_mean, group=.id, colour="b"), 
            size = 0.1, alpha=0.1, show.legend = TRUE)+
  geom_line(data=filter(mid_plotting), 
            aes(x=Year, y=G_mean, group=.id, colour="c"), 
            size = 0.1, alpha=0.1, show.legend = TRUE)+
  scale_colour_manual(values=c("black", "grey", "orange"),
                      labels=c("Phenotype", "Plasticity", "Evolution")) +
  theme(legend.position='bottom', legend.title = element_blank(),
        legend.key.width = unit(1,"cm"))+ 
  guides(colour=guide_legend(override.aes = list(size = 2)))+
  expand_limits(y=c(-10,55), x=2101)

######################################################################

### FIGURE 3 ####

## EXINCTION PROBABILITY
extinction_checks[which(is.na(extinction_checks[,2])),2] <- 0
(c(sum(extinction_checks[1:1000,2]),
   sum(extinction_checks[1001:2000,2]),
   sum(extinction_checks[2001:3000,2]))/1000)*100

extinction_checks$max_sync <- extinction_checks$max_sync + 13

model <- glm(ex_prob ~ max_sync, 
             data = extinction_checks, 
             family='binomial')

summary(model)

# Create a grid of X and Y values to predict
# can then fill whole area
preds = expand.grid(
  max_sync = seq(0, 50, length.out = 100),
  ex_prob = seq(0, 1, length.out = 100)
)
preds$pred_p <- predict(model, preds, type = "response")

pred_line <- predict(model, data.frame(max_sync = seq(0,50,length.out=100)),
                     type="response")
pred_line <- as.data.frame(cbind(pred_line, as.numeric(seq(0,50,length.out=100))))
colnames(pred_line) <- c("pred", "x")

# choose colours
colours <- heat.colors(10)
colfunc<-colorRampPalette(c(colours[10], colours[2]))

Extinction_plot<- 
  ggplot(preds, aes(x = max_sync, y = ex_prob)) +
  xlim(c(0,50))+
  ylab("Extinction probability")+
  xlab("Maximum annual asynchrony (days)")+
  theme_minimal() +
  geom_tile(aes(fill = pred_p)) +
  scale_fill_gradientn(colors=colfunc(5), 
                       name="Extinction probability") +
  expand_limits(x=c(0,50))+
  geom_jitter(data=extinction_checks, aes(x=max_sync,
                                           y=ex_prob), 
               colour="grey50", width=0.01, height=0.01,
               alpha=0.5, pch=16)+
  geom_line(data=pred_line, aes(x=x,y=pred))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        plot.margin = margin(t = 1, r = 1, b = 0.5, l = 1, unit = "cm"))

# perturbation
perturb_results<-read.csv("Perturbations_results.csv", header=T)

# finding boundaries
x <- perturb_results[154:204,1]-perturb_results[(154+25),1]
y <- perturb_results[205:255,1]-perturb_results[(205+25),1]

z <- abs(x)-abs(y)
which(z < 0)

which(perturb_results[205:255,1] < 1)
changer <- seq(-4,4,length.out=51) # 0.05 change < 1 standard deviation. So do 0.2
mean(descale$Half_fall)-((changer[20]*sd(descale$Half_fall))+mean(descale$Half_fall))
mean(descale$Half_fall)-((changer[36]*sd(descale$Half_fall))+mean(descale$Half_fall))

perturb_results_new0 <- data.frame(changer = changer,
                                  ST = perturb_results[1:51,1],
                                  WP = perturb_results[52:102,1],
                                  WT = perturb_results[103:153,1],
                                  SP = perturb_results[154:204,1],
                                  S = perturb_results[205:255,1])

perturb_results_new <- pivot_longer(perturb_results_new0, -changer, 
                                    names_to="variable", 
                                    values_to="pop_size_50") %>% 
  mutate(
         variable = factor(variable, levels = c("ST", "SP", "WT", "WP", "S"),
                            labels = c(
                                       "ST" = "Spring\nTemperature",
                                       "SP" = "Spring\nPrecipitation", 
                                       "WT" = "Winter\nTemperature",
                                       "WP" = "Winter\nPrecipitation", 
                                       "S" = "Synchrony")))

colours <- c("Spring\nTemperature" = "red", 
             "Spring\nPrecipitation" = "red", 
             "Winter\nTemperature" = "blue", 
             "Winter\nPrecipitation" = "blue", 
             "Synchrony" = "green")
linetypes <- c("Spring\nTemperature" = 1, 
               "Spring\nPrecipitation" = 3, 
               "Winter\nTemperature" = 1, 
               "Winter\nPrecipitation" = 3, 
               "Synchrony" = 1)

Perturb <- ggplot(perturb_results_new, aes(x = changer, y = pop_size_50)) +
  ylab("Population size")+
  xlab("Change in scaled driver")+
  geom_line(aes(changer, pop_size_50, color=variable, linetype=variable)) +
  scale_colour_manual(values = colours, 
                      guide = guide_legend(ncol=3)) +
  scale_linetype_manual(values = linetypes)+
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=7),
        plot.margin = margin(t = 1, r = 1, b = 0.5, l = 0.5, unit = "cm"))
  

Figure3 <- Extinction_plot + Perturb


### FIGURE 4 ####

extinction_checks$PS <- rep(0, 3000)
minPS <- 40
maxPS <- 219
meanPS <- 113

## LOW EMISSIONS POP SIZE
  
low_popsize <- ggplot(low_plotting,aes(Year, PS))+
  theme_minimal() +
  coord_cartesian(ylim=c(0,450)) + 
  ylab("Population size") +
  xlab(NULL)+
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
  theme(legend.position='none')+  
  expand_limits(y=c(-1,500), x=2101)+
  geom_segment(aes(x=2010,xend=2020,y=meanPS,yend=meanPS), color="white", size=0.75)+
  geom_segment(aes(x=2010,xend=2020,y=minPS,yend=minPS), color="white", 
               size=0.75, linetype=3)+
  geom_segment(aes(x=2010,xend=2020,y=maxPS,yend=maxPS), color="white", 
               size=0.75, linetype=3)

## MID EMISSIONS POP SIZE

mid_popsize <- ggplot(mid_plotting,aes(Year, PS))+
  theme_minimal() +
  coord_cartesian(ylim=c(0,450)) + 
  ylab(NULL) +
  ggtitle("Mid emissions") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", 
                  contour = FALSE, 
                  n=100, show.legend = F) +
  scale_fill_viridis_c(begin=0.1, option="B", trans="sqrt")+
  geom_line(data=filter(mid_plotting, .id < 11), 
            aes(x=Year, y=PS, group=.id), 
            colour="white", size = 0.1)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_rug(data=extinction_checks[1001:2000,], 
           mapping=aes(x=year+2009, colour="Extinction events"),
           na.rm=T, show.legend=T, position ="jitter",
           sides="t")+
  scale_colour_manual(values="orangered") +
  theme(legend.position='bottom', legend.title = element_blank(), axis.text.y = element_blank()) +  
  expand_limits(y=c(-1,500), x=2101)


## HIGH EMISSIONS POP SIZE

high_popsize <- ggplot(high_plotting,aes(Year, PS))+
  theme_minimal() +
  coord_cartesian(ylim=c(0,450)) + 
  ylab("") +
  xlab(NULL)+
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
  geom_segment(aes(x=2090,xend=2101,y=meanPS,yend=meanPS), color="white", size=0.75)+
  geom_segment(aes(x=2090,xend=2100,y=minPS,yend=minPS), color="white", 
               size=0.75, linetype=3)+
  geom_segment(aes(x=2090,xend=2100,y=maxPS,yend=maxPS), color="white", 
               size=0.75, linetype=3)


# combine

Figure4 <- low_popsize +
  mid_popsize + high_popsize

######################################################################

### Figure S1 ####

set.seed(10)
marker <- sample(1:1000,1)

markers <- 1961:2010
test_data2 <- descale[which(descale$Year==markers[1])[1],2:21]
for(i in 2:51){
  test_data2 <- rbind(test_data2, descale[which(descale$Year==markers[i])[1],2:21])
}
summary(test_data2)
test2 <- test_data2[-51,]


par(mar=c(2,5,2,1))
## SPRING TEMP GT *3
plot(low_predictions[[1]]$s_temp_max~low_predictions[[1]]$Year, axes=F, type='n',
     ylab=expression(paste("Spring temperature - \nPredator (",degree,"C)")), ylim=c(5,25), main="Low", xlab="")
#axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$s_temp_max, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_max))

par(mar=c(2,2,2,2))
plot(mid_predictions[[1]]$s_temp_max~mid_predictions[[1]]$Year, axes=F, type='n',
     ylim=c(5,25), main="Mid", ylab="", xlab="")
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$s_temp_max, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_max))

par(mar=c(2,2,2,2))
plot(high_predictions[[1]]$s_temp_max~high_predictions[[1]]$Year, axes=F, type='n',
     ylim=c(5,25), main = "High", ylab="", xlab="")
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_temp_max, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_max))

par(mar=c(2,1,2,4))
plot(high_predictions[[1]]$s_temp_max~high_predictions[[1]]$Year, axes=F, type='n',
     ylim=c(5,25), main = "Single", ylab="", xlab="")
#axis(1)
#axis(2, las=1)
for(i in marker){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_temp_max, type='l', col=alpha("red",1))}
abline(h=mean(test2$Spring_temp_max))

## SPRING TEMP CAT
par(mar=c(2,5,2,1))
plot(low_predictions[[1]]$s_temp_cat~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab=expression(paste("Spring temperature - \nPrey (",degree,"C)")), ylim=c(5,25))
#axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$s_temp_cat, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_cat))

par(mar=c(2,2,2,2))
plot(mid_predictions[[1]]$s_temp_cat~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(5,25))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$s_temp_cat, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_cat))

par(mar=c(2,2,2,2))
plot(high_predictions[[1]]$s_temp_cat~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(5,25))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_temp_cat, type='l', col=alpha("red",0.01))}
abline(h=mean(test2$Spring_temp_cat))

par(mar=c(2,1,2,4))
plot(high_predictions[[1]]$s_temp_cat~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(5,25))
#axis(1)
#axis(2, las=1)
for(i in marker){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_temp_cat, type='l', col=alpha("red",1))}
abline(h=mean(test2$Spring_temp_cat))

## SPRING PRECIP

par(mar=c(2,5,2,1))
plot(low_predictions[[1]]$s_precip~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="Spring precipitation (mm)", ylim=c(0,400))
#axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$s_precip, type='l', lty=2, col=alpha("red",0.01))}
abline(h=mean(test2$spring_precip_t))

par(mar=c(2,2,2,2))
plot(mid_predictions[[1]]$s_precip~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,400))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$s_precip, type='l', lty=2, col=alpha("red",0.01))}
abline(h=mean(test2$spring_precip_t))

par(mar=c(2,2,2,2))
plot(high_predictions[[1]]$s_precip~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,400))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_precip, type='l', lty=2, col=alpha("red",0.01))}
abline(h=mean(test2$spring_precip_t))

par(mar=c(2,1,2,4))
plot(high_predictions[[1]]$s_precip~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,400))
#axis(1)
#axis(2, las=1)
for(i in marker){lines(high_predictions[[i]]$Year, high_predictions[[i]]$s_precip, type='l', lty=2, col=alpha("red",1))}
abline(h=mean(test2$spring_precip_t))

## WINTER TEMP
par(mar=c(2,5,2,1))
plot(low_predictions[[1]]$w_temp~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab=expression(paste("Winter temperature (",degree,"C)")), ylim=c(0,15))
#axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$w_temp, type='l', col=alpha("blue",0.01))}
abline(h=mean(test2$winter_temp))

par(mar=c(2,2,2,2))
plot(mid_predictions[[1]]$w_temp~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,15))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$w_temp, type='l', col=alpha("blue",0.01))}
abline(h=mean(test2$winter_temp))

par(mar=c(2,2,2,2))
plot(high_predictions[[1]]$w_temp~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,15))
#axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(high_predictions[[i]]$Year, high_predictions[[i]]$w_temp, type='l', col=alpha("blue",0.01))}
abline(h=mean(test2$winter_temp))

par(mar=c(2,1,2,4))
plot(high_predictions[[1]]$w_temp~high_predictions[[1]]$Year, axes=F, type='n',
     xlab="", ylab="", ylim=c(0,15))
#axis(1)
#axis(2, las=1)
for(i in marker){lines(high_predictions[[i]]$Year, high_predictions[[i]]$w_temp, type='l', col=alpha("blue",1))}
abline(h=mean(test2$winter_temp))

## WINTER PRECIP
par(mar=c(4,5,2,1))
plot(low_predictions[[1]]$w_precip~low_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="Winter precipitation (mm)", ylim=c(0,400))
axis(1)
axis(2, las=1)
for(i in 1:1000){lines(low_predictions[[i]]$Year, low_predictions[[i]]$w_precip, type='l', lty=2, col=alpha("blue",0.01))}
abline(h=mean(test2$winter_precip_t))

par(mar=c(4,2,2,2))
plot(mid_predictions[[1]]$w_precip~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="", ylim=c(0,400))
axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$w_precip, type='l', lty=2, col=alpha("blue",0.01))}
abline(h=mean(test2$winter_precip_t))

par(mar=c(4,2,2,2))
plot(mid_predictions[[1]]$w_precip~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="", ylim=c(0,400))
axis(1)
#axis(2, las=1)
for(i in 1:1000){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$w_precip, type='l', lty=2, col=alpha("blue",0.01))}
abline(h=mean(test2$winter_precip_t))

par(mar=c(4,1,2,4))
plot(mid_predictions[[1]]$w_precip~mid_predictions[[1]]$Year, axes=F, type='n',
     xlab="Year", ylab="", ylim=c(0,400))
axis(1)
#axis(2, las=1)
for(i in marker){lines(mid_predictions[[i]]$Year, mid_predictions[[i]]$w_precip, type='l', lty=2, col=alpha("blue",1))}
abline(h=mean(test2$winter_precip_t))

######################################################################

### Figure S2 ####

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

## TESTING RATES
rate_store <- rep(NA, 1000)

for(i in 1:1000){
  rates <- high_results[[i]]$G_mean[6:91]-high_results[[i]]$G_mean[5:90]
  model <- lm(rates~c(1:86))
  if(coef(model)[2]<0 & summary(model)$coefficients[8]<0.05){rate_store[i] <- T}else
  {rate_store[i] <- F}
}


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

# make into dataframe for anova
change_E_dataframe <- data.frame(changes = c(changes_E[,1], changes_E[,2]),
                                 group = c(rep(1,1000), rep(2,1000)),
                                 pair = rep(seq(1,1000,1),2))

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

plot(changes ~ jitter(group, 0.03), data=change_E_dataframe,
     pch=16, col=group, axes=F,
     ylab="Variance in mean", 
     xlab="Time relative to change point",
     ylim=c(0,1),
     xlim=c(0.8,2.2))
axis(1, at=c(1,2), labels=c("Prior", "Post"))
axis(2, las=1, at=seq(0, 1, 0.2))
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

######################################################################

### Figure S3 ####

# load data
load("only_ST_results.RData")
load("only_SP_results.RData")
load("only_WP_results.RData")
load("only_WT_results.RData")

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


par(mar=c(2,5,2,0))
# Population size
plotter(no_ST_results, "Population size", "", "a. held ST", x.axis=F, y.axis=T, adj=0.5,
        ylim=c(0,600))
par(mar=c(2,1,2,1))
plotter(no_SP_results, "", "", "b. held SP", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))
par(mar=c(2,1,2,1))
plotter(no_WT_results, "", "", "c. held WT", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))
par(mar=c(2,0,2,2.5))
plotter(no_WP_results, "", "", "d. held WP", x.axis=F,y.axis=F, adj=0.5, ylim=c(0,600))


par(mar=c(4,5,0,0))
# Hatch
plotter3(no_ST_results, "Hatch date", "", "", descale, x.axis=T, y.axis=T)
par(mar=c(4,1,0,1))
plotter3(no_SP_results, "", "Time", "", descale, x.axis=T, y.axis=F, adj2=1)
par(mar=c(4,1,0,1))
plotter3(no_WT_results, "", "(years)", "", descale, x.axis=T, y.axis=F, adj2=0)
par(mar=c(4,0,0,2.5))
plotter3(no_WP_results, "", "", "", descale, x.axis=T, y.axis=F)


######################################################################

### SOM TREND TABLE ####

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
summary(trends$G) # 106/1000 decrease
summary(trends$E[1001:2000]) # 296/1000 decrease
summary(trends$E[2001:3000]) # 517/1000 decrease

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



