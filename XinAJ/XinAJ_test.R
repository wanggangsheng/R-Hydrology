## Xin'AnJiang Hydrologic Model
## Example code to conduct hydrologic modeling
## Gangsheng Wang (wang.gangsheng@gmail.com)
## March 15, 2021

##------------------------------------------------------------------------------
## check data only:
# install.packages('devtools')
# devtools::install_github('Sibada/XAJ')

# library(XAJ)
# data(SB)
# write.table(SB, file="XinAJ_input.csv",sep="\t",row.names = F)
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)  ## ggscatter
library(tidyr)  ## gather
#Clean up Environment
rm(list = ls())
setwd("/Users/wgs/ownCloud/Rcode/R-Hydrology/XinAJ")
source("XinAJ.R")

## Precipitation + Evaporation (pan evaporation or potential evapotranspiration)
# input_data <- read.table("XinAJ_PE.csv", sep="\t",head=T)
## Precipitation (mm) + Evaporation (mm) + Observed streamflow (m3/s)
input_data <- read.table("XinAJ_PEQ.csv", sep="\t",head=T)
head(input_data)
nt <- nrow(input_data)

XAJpar <- XinAJpar()
str(XAJpar)

param <- XAJpar$param_value
npar <- length(param)
param
##-----------------------------------------------------------------------------
XAJpar
## Sketch for the saturation-excess runoff generation mechanism
WM <- sum(XAJpar$param_value[3:5])  ## watershed-average soil moisture storage (SMS) capacity 
B <- XAJpar$param_value[7]          ## exponent of the SMS capacity curve
print(paste0("WM= ",WM, "; B= ",B)) 
WMM <- WM * (1. + B)  ## watershed max SMS capacity

xpos <- 0.2
W <- seq(0,WMM,by=1)
fsat <- SatFrac(B,WMM,W)
plot(x=fsat,y=W,type='l',xlab="Saturation Fraction",ylab="Max Soil Moisture Storage (mm)")
abline(h=WMM,col="blue",lty=2,lwd=3)
text(x=xpos,y=WMM*0.95,col="blue","Wmm: watershed Max SMS capacity")
abline(h=WM,col="red",lty=2,lwd=2)
text(x=xpos,y=WM*1.05,col="blue","Wm: : watershed Mean SMS capacity")
abline(h=WM/2,col="black",lty=1,lwd=2)
text(x=xpos,y=WM/2*1.1,col="blue","PE+a: Current watershed Max SMS")
abline(h=WM/6,col="black",lty=1,lwd=2)
text(x=xpos,y=WM/6*1.15,col="blue","a: Initial watershed Max SMS")
axis(side=2, at=seq(0, WMM, by=20))

##-----------------------------------------------------------------------------
## Randomly Generate Parameter Values from their ranges 
for (i in 1:npar) {
  param[i] <- runif(1, min = XAJpar$param_lower[i], max = XAJpar$param_upper[i])
}
param[10] <- 0.70 - param[11]  ## KI = 0.70 - KG
names(param) <- XAJpar$param_names
##-----------------------------------------------------------------------------
param


dt <- 24 ## days?
area <- 1000 ## watershed area, km^2
UH <- IUH(param[14], param[15],dt)
plot(UH, type="l",xlab="Time", ylab="Unit Hydrograph")
points(UH)

XAJ_test <- XinAJrun(input_data$PREC,input_data$EVAP,param,UH, area, dt)
str(XAJ_test)

plot(XAJ_test$Qsim,type="l", xlab="Time", ylab="Qsim")
head(XAJ_test)
mean_sim <- colMeans(XAJ_test)
ETsim <- mean_sim[2:4]
Rsim <- mean_sim[10:12]
Qsim <- mean_sim[14:16]

pie(ETsim,main="ET")
pie(Rsim, main="Rsim")
pie(Qsim, main="Qsim")

## Save the Qsim as Qobs for model evaluation
# input <- data.frame(input_data,XAJ_test$Qsim*runif(1,0.8,1.2))
# write.table(round(input,digits=3), file="XinAJ_input2.csv",sep="\t",row.names = F)

R2_XAJ <- R2(input_data$Qobs,XAJ_test$Qsim)
R2_XAJ

##------------------------------------------------------------------------------
## scatter plot
size_axis <- 20
size_title <- size_axis+2
size_regression <- size_axis/3

df <- data.frame(seq(1,nt),input_data$Qobs,XAJ_test$Qsim)
legend <- c("Qobs","Qsim")
colnames(df) <- c("Time",legend)

xmax1 <- ceiling(max(df$Qobs)/2)*2
xmax2 <- ceiling(max(df$Qsim)/2)*2
ymax <- max(xmax1,xmax2)
ylim1 <- c(0, ymax)

title1 <- bquote("Xin'AnJiang Model: "*R^2*" = "*.(round(R2_XAJ,digits=3)))
xlab1 <- bquote(bold("Observed Streamflow ")*bold(" ("*m^3*" "*s^bold("-1")*")"))
ylab1 <- bquote(bold("Simulated Streamflow ")*bold(" ("*m^3*" "*s^bold("-1")*")"))

ggscatter(df, x="Qobs",y="Qsim",color="blue",shape=21,size=3, add="reg.line") +
  geom_abline(intercept = 0, slope = 1) + 
  stat_cor(aes(label = paste('r = ', ..r..)),
           label.x = 3, label.y =ymax*1, size=size_regression) +
  stat_regline_equation(label.x = 3, label.y =ymax*0.95, size=size_regression) +
  coord_cartesian(xlim=ylim1, ylim=ylim1) +
  scale_x_continuous(breaks=seq(ylim1[1],ylim1[2],10)) +
  scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],10)) +
  ggtitle(title1) +
  xlab(xlab1) +
  ylab(ylab1) +
  theme(
    panel.border = element_rect(size=2, colour = "black",fill=NA),
    plot.title = element_text(hjust =0.5,size=size_title),
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

##------------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Monte Carlo simulations
## Randomly Generate multiple sets of Parameter Values from their ranges 
nrun <- 100
param_sample <- matrix(nrow = nrun, ncol = npar)
R2_XAJ <- rep(0.,nrun)
for (i in 1:npar) {
  param_sample[,i] <- runif(nrun, min = XAJpar$param_lower[i], max = XAJpar$param_upper[i])
}
colnames(param_sample) <- XAJpar$param_names
param_sample[,10] <- 0.70 - param_sample[,11]  ## KI = 0.70 - KG
head(param_sample)

for(i in 1:nrun) {
  print(paste0("No. of model run = ",i," of ",nrun))
  UH <- IUH(param_sample[i,14], param_sample[i,15],dt)
  XAJ_test <- XinAJrun(input_data$PREC,input_data$EVAP,param_sample[i,],UH, area, dt)
  R2_XAJ[i] <- R2(input_data$Qobs,XAJ_test$Qsim)
}

##---------------------------------------------------------
## empirical cumulative probability 
cumulative_probability <- ecdf(R2_XAJ)

plot(cumulative_probability, 
     main="", xlim=c(-2,1),
     xlab=bquote(R^bold("2")),ylab="Cumulative probability")
axis(side=1, at=seq(0, 1, by=0.2))
##---------------------------------------------------------

## visualization of the sensitivity using dotty plots:
par_names <- colnames(parameters)
par_R2 <- data.frame(param_sample,R2_XAJ)

par_R2_gather <- gather(par_R2,Parameter,Value,-R2_XAJ)

summary(par_R2_gather)
ylab1 <-bquote("Coefficient of Determination: "*R^bold("2"))

ggplot(data = par_R2_gather, aes(x = Value, y = R2_XAJ)) + 
  geom_point(aes(color = Parameter)) + 
  ylab(ylab1) +
  scale_y_continuous(breaks=seq(0,1, by=0.2), limits=c(0,1)) +
  theme_bw() +
  facet_wrap(~Parameter,scales = "free",  ncol=3) 

##------------------------------------------------------------------------------

max(R2_XAJ)
iR2max <- which.max(R2_XAJ)
iR2max

## Re-run the model with the best parameters------------------------------------
param <- param_sample[iR2max,]
UH <- IUH(param_sample[iR2max,14], param_sample[iR2max,15],dt)
XAJ_test <- XinAJrun(input_data$PREC,input_data$EVAP,param,UH, area, dt)
R2_XAJ <- R2(input_data$Qobs,XAJ_test$Qsim)
R2_XAJ

##------------------------------------------------------------------------------
## calculate mean values and ratios
head(XAJ_test)
mean_sim <- colMeans(XAJ_test)
ETmean <- mean_sim[1]
Rmean <- mean_sim[9]

ETsim <- mean_sim[2:4]
Rsim <- mean_sim[10:12]
Qsim <- mean_sim[14:16]

Pmean <- mean(input_data$PREC)
PETmean <- mean(input_data$EVAP)
ET_P <- ETmean/Pmean ## ratio of ET to Precipitation
R_P <- Rmean/Pmean
PET_P <- PETmean/Pmean
ET_PET <- ETmean/PETmean

df1 <- data.frame(ET_P = ET_P, R_P = R_P, PET_P = PET_P, ET_PET = ET_PET)
rownames(df1) <- "Ratio"
df1

##------------------------------------------------------------------------------
## Time-Series Plot

size_axis <- 20
size_title <- size_axis+2

title1 <- bquote("Xin'AnJiang Model: "*R^2*" = "*.(round(R2_XAJ,digits=3)))

xlab1 <- 'Time'
ylab1 <- bquote('Streamflow ('*m^3*' '*s^"-1"*')')

df <- data.frame(seq(1,nt),input_data$Qobs,XAJ_test$Qsim)
legend <- c("Qobs","Qsim")
colnames(df) <- c("Time",legend)

color1 <- c("blue","red")

sp1 <- ggplot(data=df, aes(x=Time)) + 
  ggtitle(title1) +
  xlab(xlab1) + ylab(ylab1) + 
  geom_line(aes(y = Qobs, color=legend[1]), size=1) +
  geom_line(aes(y = Qsim, color=legend[2]), size=1) + 
  
  scale_color_manual(name=NULL,values = color1,limits=legend) +
  guides(col = guide_legend(nrow=1)) +
  theme(
    plot.title=element_text(hjust=0.5, size=size_title,face="bold"),
    plot.subtitle=element_text(hjust=-0.03, size=size_title, face="bold"),
    strip.text.x=element_blank(),
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold"),
    legend.text=element_text(size=size_axis),
    legend.title=element_text(size=size_axis),
    legend.position=c(0.5,0.9)
  )

print(sp1)
fn1 <- paste0("XAJ_eval.pdf")
ggsave(filename=fn1, plot=sp1,width=50,height=20,units="cm", dpi=600)

##------------------------------------------------------------------------------
## scatter plot
size_axis <- 20
size_title <- size_axis+2
size_regression <- size_axis/3

xmax1 <- ceiling(max(df$Qobs)/2)*2
xmax2 <- ceiling(max(df$Qsim)/2)*2
ymax <- max(xmax1,xmax2)
ylim1 <- c(0, ymax)

title1 <- bquote("Xin'AnJiang Model: "*R^2*" = "*.(round(R2_XAJ,digits=3)))
xlab1 <- bquote(bold("Observed Streamflow ")*bold(" ("*m^3*" "*s^bold("-1")*")"))
ylab1 <- bquote(bold("Simulated Streamflow ")*bold(" ("*m^3*" "*s^bold("-1")*")"))

sp2 <- ggscatter(df, x="Qobs",y="Qsim",color="blue",shape=21,size=3) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_cartesian(xlim=ylim1, ylim=ylim1) +
  scale_x_continuous(breaks=seq(ylim1[1],ylim1[2],10)) +
  scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],10)) +
  ggtitle(title1) +
  xlab(xlab1) +
  ylab(ylab1) +
  theme(
    panel.border = element_rect(size=2, colour = "black",fill=NA),
    plot.title = element_text(hjust =0.5,size=size_title),
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

print(sp2)
fn2 <- paste0("XAJ_eval_scatter.pdf")
ggsave(filename=fn2, plot=sp2,width=30,height=30,units="cm", dpi=600)

##------------------------------------------------------------------------------

xlab1 <- 'Time'
ylab1 <- 'Q'

df <- data.frame(seq(1,length(XAJ_test$Qsim)),XAJ_test$Qsim,XAJ_test$QSsim,XAJ_test$QIsim,XAJ_test$QGsim)
legend <- c("Qsim","QSsim","QIsim","QGsim")
colnames(df) <- c("Time",legend)

color1 <- c("red","cornflowerblue","blue","darkgreen")

sp3 <- ggplot(data=df, aes(x=Time)) + 
  ggtitle(title1) +
  xlab(xlab1) + ylab(ylab1) + 
  geom_line(aes(y = Qsim, color=legend[1]), size=1) +
  geom_line(aes(y = QSsim, color=legend[2]), size=1) + 
  geom_line(aes(y = QIsim, color=legend[3]), size=1) +
  geom_line(aes(y = QGsim, color=legend[4]), size=1) +
  
  
  scale_color_manual(name=NULL,values = color1,limits=legend) +
  guides(col = guide_legend(nrow=1)) +
  theme(
    plot.title=element_text(hjust=0.5, size=size_title,face="bold"),
    plot.subtitle=element_text(hjust=-0.03, size=size_title, face="bold"),
    strip.text.x=element_blank(),
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold"),
    legend.text=element_text(size=size_axis),
    legend.title=element_text(size=size_axis),
    legend.position=c(0.5,0.9)
  )

print(sp3)
fn3 <- paste0("XAJ_Qsim.pdf")
ggsave(filename=fn3, plot=sp3,width=50,height=20,units="cm", dpi=600)

##------------------------------------------------------------------------------