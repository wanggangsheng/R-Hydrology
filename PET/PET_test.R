## Calculate potential/reference evapotranspiration (ETo)
## (1) Penman-Monteith Method
## (2) Hargreaves Method
## (3) Priestley Taylor Method

## Author: Gangsheng Wang (wang.gangsheng@gmail.com)
## March 1, 2021

library(ggplot2)
library(ggpubr) #ggscatter, stat_cor
library(scales) #as.POSIXct, date_format
#Clean up Environment
rm(list = ls())

## Please change the following to your own Working Directory
setwd("/Users/wgs/ownCloud/Rcode/R-Hydrology/PET")
## All functions are included in "PET.R"
source("PET.R") 

path <- getwd()

Latitude <- 46.74 ## degrees
Altitude <- 758   ## [m]
WindHeight <- 10  ## [m]
albedo <- 0.23

file0 <- "PET_input.txt"
file0 <- paste0(path,"/", file0)

## read input data from file0
df0 <- read.table(file0, sep="\t",header=TRUE)
names(df0)
head(df0)

df <- df0
df$Tmean <- 0.5*(df$Tmax + df$Tmin)
df$esTmax <- SatVP(df$Tmax)
df$esTmin <- SatVP(df$Tmin)
df$esTmean <- 0.5*(df$esTmax + df$esTmin)
df$ea <-   VP(df$Tmax, df$Tmin, df$RHmean, df$RHmean)
df$VPD <- VPD(df$Tmax, df$Tmin, df$RHmean, df$RHmean)

## inverse relative distance from Earth to Sun
df$DR <- InvDistance(df$DOY)
## solar declination
df$SolDec <- SolarDec(df$DOY)
## sunset hour angle
df$SunsetAngle <- SunsetHourAngle(Latitude,df$DOY)
## aerodynamic resistance
df$ra <- AeroRes(df$WindSpeed,WindHeight)
## slope of saturation vapor pressure curve
df$Delta <- DeltaEs(df$Tmax, df$Tmin) 
## latent heat of vaporization
df$Lambda <- Lambda(df$Tmax, df$Tmin)
## psychrometrric constant
df$Gamma <- Gamma(Altitude, df$Tmax, df$Tmin)
## potential extraterrestrail radiation
df$PotRad <- PotRad(Latitude,df$DOY)

## shortsave net radiation
df$NetRad <- NetRad(albedo, Latitude, Altitude,
                    df$DOY, df$SolarRad, df$Tmax, df$Tmin, df$RHmean, df$RHmean)

## reference ET radiation term
df$ETRad <- ETRadTerm(albedo, Latitude, Altitude,
                      df$DOY, df$SolarRad, df$Tmax, df$Tmin, 
                      df$RHmean, df$RHmean,df$WindSpeed,WindHeight)

## reference ET aerodynamic or wind term
df$ETAero <- ETAeroTerm(Altitude,
                        df$Tmax, df$Tmin, df$RHmean, df$RHmean,
                        df$WindSpeed,WindHeight)
##-----------------------------------------------------------------------------
## Penman-Monteith Reference ET
df$ET_PenmanMonteith <- ETPenmanMonteith(albedo, Latitude, Altitude, 
                           df$DOY, df$SolarRad, df$Tmax, df$Tmin, 
                           df$RHmean, df$RHmean,df$WindSpeed,WindHeight)

##-----------------------------------------------------------------------------
## Hargreaves Reference ET
df$ET_Hargreaves <- ETHargreaves(Latitude,df$DOY, df$Tmax, df$Tmin)

##-----------------------------------------------------------------------------
## Priestley-Taylor Reference ET
## Alfa ranging from 1.26 to 1.70
alfa <- 1.26
df$ET_PriestleyTaylor <- ETPriestleyTaylor(alfa,albedo, Latitude, Altitude, 
                                           df$DOY, df$SolarRad, df$Tmax, df$Tmin, 
                                           df$RHmean, df$RHmean)

##-----------------------------------------------------------------------------
size_axis <- 20
size_title <- size_axis+2
size_regression <- size_axis/3

## ET_PenmanMonteith vs. ET_PriestleyTaylor
cor1 <- cor(df$ET_PenmanMonteith,df$ET_PriestleyTaylor)
R2v <- R2(df$ET_PenmanMonteith,df$ET_PriestleyTaylor)
cor1
cor1^2
R2v

xmax1 <- ceiling(max(df$ET_PenmanMonteith)/2)*2  ## round to the ceiling even number
xmax2 <- ceiling(max(df$ET_PriestleyTaylor)/2)*2
ymax <- max(xmax1,xmax2)
ylim1 <- c(0, ymax)

xlab1 <- bquote(bold("Penman-Monteith ETo ")*bold(" (mm "*d^bold("-1")*")"))
ylab1 <- bquote(bold("Priestley-Taylor ETo ")*bold(" (mm "*d^bold("-1")*")"))

## scatter plot
sp1 <- ggscatter(df, x="ET_PenmanMonteith",y="ET_PriestleyTaylor",
          add="reg.line",  ## add regression line
          conf.int = TRUE, ## add confidence interval
          add.params = list(color = "red", fill = "darkgray")) +
  coord_cartesian(xlim=ylim1, ylim=ylim1) +
  scale_x_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  xlab(xlab1) +
  ylab(ylab1) +
  stat_cor(label.x = 3, label.y =ymax*1, size=size_regression) +
  stat_regline_equation(label.x = 3, label.y =ymax*0.95, size=size_regression) +
  theme(
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

print(sp1)
fn1 <- paste0(path,"/ETo_PM-PT.pdf")
ggsave(filename=fn1, plot=sp1,width=30,height=30,units="cm", dpi=600)
##------------------------------------------------------------------------------

## ET_PenmanMonteith vs. ET_Hargreaves
cor1 <- cor(df$ET_PenmanMonteith,df$ET_Hargreaves)
R2v <- R2(df$ET_PenmanMonteith,df$ET_Hargreaves)
cor1
cor1^2
R2v

xmax1 <- ceiling(max(df$ET_PenmanMonteith)/2)*2  ## round to the ceiling even number
xmax2 <- ceiling(max(df$ET_Hargreaves)/2)*2
ymax <- max(xmax1,xmax2)
ylim1 <- c(0, ymax)

ylab1 <- bquote(bold("Hargreaves ETo ")*bold(" (mm "*d^bold("-1")*")"))
sp2 <- ggscatter(df, x="ET_PenmanMonteith",y="ET_Hargreaves",
          add="reg.line",  ## add regression line
          conf.int = TRUE, ## add confidence interval
          add.params = list(color = "red", fill = "darkgray")) +
  coord_cartesian(xlim=ylim1, ylim=ylim1) +
  scale_x_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  xlab(xlab1) +
  ylab(ylab1) +
  stat_cor(label.x = 3, label.y =ymax*1, size=size_regression) +
  stat_regline_equation(label.x = 3, label.y =ymax*0.95, size=size_regression) +
  theme(
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

print(sp2)
fn2 <- paste0(path,"/ETo_PM-HG.pdf")
ggsave(filename=fn2, plot=sp2,width=30,height=30,units="cm", dpi=600)
##------------------------------------------------------------------------------
## Time-Series Plot
xlab1 <- 'Date (mm/dd/yy)'
ylab1 <- bquote(bold("ETo ")*bold(" (mm "*d^bold("-1")*")"))

legend <- c("Penman-Monteith","Priestley-Taylor","Hargreaves")
color1 <- c("cornflowerblue","darkgreen", "red")
df$Date1 <- as.Date(df$Date,'%m/%d/%y')
df$Date1 <- as.POSIXct(df$Date1)


  
  title1 <- bquote(bold("Reference Evapotranspiration (ETo)"))
  
  sp3 <- ggplot(data=df,aes(x = Date1)) + 
    ggtitle(title1) +
    xlab(xlab1) + ylab(ylab1) + 
    geom_line(aes(y = ET_PenmanMonteith, color=legend[1]), size=1) + 
    geom_line(aes(y = ET_PriestleyTaylor, color=legend[2]), size=1) +
    geom_line(aes(y = ET_Hargreaves, color=legend[3]), size=1) +

    scale_x_datetime(date_break="2 month",labels=date_format("%m/%d/%y")) +
    scale_color_manual(name="Model",values = color1) +
    guides(col = guide_legend(nrow=3)) +
    theme(
      plot.title=element_text(hjust=0.5, size=size_title,face="bold"),
      plot.subtitle=element_text(hjust=-0.03, size=size_title, face="bold"),
      strip.text.x=element_blank(),
      axis.text=element_text(size=size_axis),
      axis.title=element_text(size=size_title,face="bold"),
      legend.text=element_text(size=size_axis),
      legend.title=element_text(size=size_axis),
      legend.position=c(0.2,0.8)
    )
  
  print(sp3)
  fn3 <- paste0(path,"/ETo_TimeSeries.pdf")
  ggsave(filename=fn3, plot=sp3,width=50,height=20,units="cm", dpi=600)

##==============================================================================
## Validate the above Penman-Monteith algorithm by the results from packages("Evapotranspiration")
## please install the package if you have NOT done it
# install.packages("Evapotranspiration")
library(Evapotranspiration)
data("processeddata")
data("constants")

constants$lat
constants$Elev
constants$z

## analyze the structure of the data
str(processeddata)

results <- ET.PenmanMonteith(processeddata, constants, ts="daily", solar="sunshine hours",
                             wind="yes", crop = "short", message="no", 
                             AdditionalStats="yes", save.csv="no")
str(results)

#create data frame with 0 rows and 4 columns
df1 <- data.frame(matrix(ncol=4,nrow=length(processeddata$Date.daily)))

#provide column names
colnames(df1) <- c("Date","DOY", "ETo1", "ETo2")

df1$Date <- processeddata$Date.daily
df1$DOY <- processeddata$J

df1$ETo1 <- results$ET.Daily

albedo = 0.23

## Estimate Rs (Solar radiation) when measured data are NOT available
## Rs data are NOT provided in this package ("Evapotranspiration")
Rs <- SolarRadiation(constants$as, constants$bs, 
                     constants$lat, processeddata$J, processeddata$n)

df1$ETo2 <- ETPenmanMonteith(albedo, constants$lat, constants$Elev, 
                             processeddata$J, Rs, 
                             processeddata$Tmax, processeddata$Tmin, 
                             processeddata$RHmax, processeddata$RHmin,
                             processeddata$uz,constants$z)


## my-Package vs. Evapotranspiration-Package
xmax1 <- ceiling(max(df1$ETo1)/2)*2
xmax2 <- ceiling(max(df1$ETo2)/2)*2
ymax <- max(xmax1,xmax2)
ylim1 <- c(0, ymax)


xlab1 <- bquote(bold("Evapotranspiration-Package ETo ")*bold(" (mm "*d^bold("-1")*")"))
ylab1 <- bquote(bold("my-Package ETo ")*bold(" (mm "*d^bold("-1")*")"))

sp4 <- ggscatter(df1, x="ETo1",y="ETo2",
          add="reg.line",  ## add regression line
          conf.int = TRUE, ## add confidence interval
          add.params = list(color = "red", fill = "darkgray")) +
  coord_cartesian(xlim=ylim1, ylim=ylim1) +
  scale_x_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  xlab(xlab1) +
  ylab(ylab1) +
  stat_cor(label.x = 3, label.y =ymax*1, size=size_regression) +
  stat_regline_equation(label.x = 3, label.y =ymax*0.95,size=size_regression) +
  theme(
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

print(sp4)
fn4 <- paste0(path,"/ETo_Validate.pdf")
ggsave(filename=fn4, plot=sp4,width=30,height=30,units="cm", dpi=600)


##------------------------------------------------------------------------------
## plot saturation vapor pressure vs. temperature
temp1 <- seq(from=0, to=60, by=1)
es1 <- SatVP(temp1)

##-----------------------------------------
## example to use 'for loop' 
## Not recommend: too tedious
temp1 <- vector(mode="numeric", length=61)
es1 <- vector(mode="numeric", length=61)
for(i in 1:61) {
  temp1[i] <- i - 1 
  es1[i] <- SatVP(temp1[i])
}
##-----------------------------------------

df_es <- data.frame(temp1, es1)

xlab1 <- bquote(bold("Temperature (°C)"))
ylab1 <- bquote(bold("Saturation vapor pressure (kPa)"))

xlim1 <-c (0, 60)
## scatter plot
sp1 <- ggscatter(df_es, x="temp1",y="es1", shape=1, size=3, color="blue") +
  geom_line(size=1, color="red") + 
  scale_x_continuous(breaks=seq(xlim1[1],xlim1[2],10)) +
  # scale_y_continuous(breaks=seq(ylim1[1],ylim1[2],2)) +
  xlab(xlab1) +
  ylab(ylab1) +
  theme(
    axis.text=element_text(size=size_axis),
    axis.title=element_text(size=size_title,face="bold")
  )

print(sp1)
fn1 <- paste0(path,"/SatVP-Temp.pdf")
ggsave(filename=fn1, plot=sp1,width=30,height=30,units="cm", dpi=600)
##------------------------------------------------------------------------------
##==============================================================================

