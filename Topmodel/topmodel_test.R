## Reference Code: https://paramo.cc.ic.ac.uk/topmodel_tutorial
## There are some problems in the original code
## Code corrected and added by Gangsheng Wang (wang.gangsheng@gmail.com)

## install and load the required packages:
# install.packages("topmodel")
# install.packages("Hmisc")
# install.packages("raster")
# install.packages("rasteVis")
# install.packages("tidyr")

library(topmodel)
library(Hmisc)
library(raster)
library(rasterVis) ## gplot
library(tidyr)  ## gather
#Clean up Environment
rm(list = ls())

## please change to your own directory
setwd("/Users/wgs/ownCloud/Rcode/R-Hydrology/Topmodel")
##------------------------------------------------------------------------------
############ PART 0: topographical analysis ##############

# Topmodel requires to types of topographical information: the topographic index distribution, and a channel flow delay function
# As topmodel is a semidistributed model, it does not use the maps directly, but instead uses histograms.
# Many GIS packages include routines to calculate these maps, but you can do it directly in R as well,
# starting from a digital elevation model (DEM).

# The DEM has to be imported as a matrix, which can then be processed by topidx().
# Take for instance this minimalistic DEM, saved in a test file called "DEM.txt"
# Values outside the catchment are given the value -9999 (this can be any other value that, obviously, does not occur in the DEM values):

## Topmodel parameters:
##------------------------------------------------------------------------------
# qs0	Initial subsurface flow per unit area [m]
# lnTe	log of the areal average of T0 [m2/h]
# m	Model parameter controlling the rate of decline of transmissivity in the soil profile, see Beven, 1984
# Sr0	Initial root zone storage deficit [m]
# Srmax	Maximum root zone storage deficit [m]
# td	Unsaturated zone time delay per unit storage deficit [h/m]
# vch	channel flow outside the catchment [m/h] (currently not used)
# vr	channel flow inside catchment [m/h]
# k0	Surface hydraulic conductivity [m/h]
# CD	capillary drive, see Morel-Seytoux and Khanji (1974)
# dt	The timestep [h]
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
dem2 <- read.table("dem2.csv", sep="\t",head=T)
dem2 <- as.matrix(dem2)

# Remove the values outside the catchment:

dem2[dem2==-9999] <- NA

# You may want to plot the DEM to see whether everything looks OK:
image(dem2)
image(dem2, useRaster=TRUE)
##------------------------------------------------------------------------------
# Then calculate the topographic index, the resolution should be in [m].
# Here we use the DEM from Huagrahuma as an example:

data(huagrahuma.dem)
DEM <- sinkfill(huagrahuma.dem, res=25, degree=0.1)

# write.table(huagrahuma.dem, file="dem3.csv",sep="\t",row.names = F)

topindex0 <- topidx(DEM, resolution=25)
str(topindex0)

area <- as.matrix(topindex0$area)
image(area)

topindex <- as.matrix(topindex0$atb)
image(topindex)

topindex2 <- raster(topindex)
theme_set(theme_bw())
gplot(topindex2) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal()


# The values need to be split into a set of classes, 
# since topmodel() is a semidistributed model that 
# lumps hydrologically similar areas into the same hydrological response units.
# Here we define 16 hydrological response units:

topidx <- make.classes(topindex,16)
##------------------------------------------------------------------------------
## Flow is routed through a delay function which represents the time spent in the channel system. 
## The parameter delay is used for this. 
## Delay is a matrix with 2 columns. 
## The first column gives the cumulative relative area. 
## The second column gives the average distance towards the outlet (m).

# the delay function is a bit more tricky because this requires cumulative fractions, but you generate it as follows:

n <- 5 # number of classes; a higher number will result in a smoother histogram
delay <- flowlength(huagrahuma.dem)*25 # TODO: add the outlet coordinates; otherwise the flowlength will be calculated to the edge of the map.
delay <- make.classes(delay, n)
delay <- delay[n:1,]
delay[,2] <- c(0, cumsum(delay[1:(n-1),2]))
##------------------------------------------------------------------------------

############ PART 1: running the rainfall-runoff model ##############

## Load the example dataset from the Huagrahuma catchment
## and attach it to the search path

data(huagrahuma)
attach(huagrahuma)

## Initial exploration of the data:

str(huagrahuma)
topidx
parameters
#rain

plot(rain, type="h")

## run the model and visualise the outcome:
multiplier <- 1
rain <- rain*multiplier
ETp <- ETp*multiplier
Qobs <- Qobs*multiplier

Qsim <- topmodel(parameters, topidx, delay, rain, ETp)
str(Qsim)
n <- nrow(Qsim)
if(is.null(n)) {
  Qsim <- Qsim
} else {
    Qsim <- Qsim[,1]
  }


plot(Qsim, type="l", col="red")
points(Qobs)

## Evaluate the model with a performance metric

NSeff(Qobs, Qsim)

##------------------------------------------------------------------------------
############ PART 2: Sensitivity analysis ##############

## Now sample all parameters at random. We take a sample size of 100

n <- 100

qs0 <- runif(n, min = 0.0001, max = 0.00025)*multiplier
lnTe <- runif(n, min = -2, max = 3)
m <- runif(n, min = 0, max = 0.1)
Sr0 <- runif(n, min = 0, max = 0.2)
Srmax <- runif(n, min = 0, max = 0.1)
td <- runif(n, min = 0, max = 3)
vch <- runif(n, min = 100, max = 2500)
vr <- runif(n, min = 100, max = 2500)
k0 <- runif(n, min = 0, max = 10)
CD <- runif(n, min = 0, max = 5)
dt <- 0.25

parameters <- cbind(qs0,lnTe,m,Sr0,Srmax,td,vch,vr,k0,CD,dt)

## run the model and evaluate with the Nash – Sutcliffe efficiency metric:
## Note: the function accepts a table of parameter sets
## (one parameter set per row)

R2 <- topmodel(parameters, topidx, delay, rain, ETp, Qobs = Qobs)
max(R2)

## empirical cumulative probability 
cumulative_probability <- ecdf(R2)

plot(cumulative_probability, 
     main="", 
     xlab=bquote(R^bold("2")),ylab="Cumulative probability")
axis(side=1, at=seq(0, 1, by=0.1))

## visualisation of the sensitivity using dotty plots:
npar <- ncol(parameters)
par_names <- colnames(parameters)
## remove the last parameter 'dt' since it's not changed
par_R2 <- data.frame(parameters[,1:(npar-1)],R2)

par_R2_gather <- gather(par_R2,Parameter,Value,-R2)

summary(par_R2_gather)
ylab1 <-bquote("Coefficient of Determination: "*R^bold("2"))

ggplot(data = par_R2_gather, aes(x = Value, y = R2)) + 
  geom_point(aes(color = Parameter)) + 
  ylab(ylab1) +
  scale_y_continuous(breaks=seq(0,1, by=0.2), limits=c(0,1)) +
  theme_bw() +
  facet_wrap(~Parameter,scales = "free",  ncol=2) 

##------------------------------------------------------------------------------

############ PART 3: GLUE uncertainty analysis ##############

## choose a behavioural threshold and remove the “bad” parameter sets:

parameters <- parameters[R2 > 0.3,]
R2 <- R2[R2 > 0.3]

## generate predictions for the behavioural parameter sets:

Qsim <- topmodel(parameters,topidx,delay,rain,ETp)

## (have a look at the predictions for the first time step:)

hist(Qsim[1,])

## construct weights based on the performance measure

weights <- R2 - 0.3
weights <- weights / sum(weights)

## make prediction boundaries by weighted quantile calculation
## (we need the Hmisc package for that)

limits <- apply(Qsim, 1, "wtd.quantile", weights = weights,
                probs = c(0.05,0.95), normwt=T)

plot(limits[2,], type="l", xlab="Time", ylab="Qobs")
points(limits[1,], type="l")
points(Qobs, col="red")

## how many measurements fall outside of the prediction limits?

outside <- (Qobs > limits[2,]) | (Qobs < limits[1,])
summary(outside)

## width of the prediction boundaries

mean(limits[2,] - limits[1,]) / mean(Qobs, na.rm=T)

##------------------------------------------------------------------------------
## Time-Series Plot

size_axis <- 20
size_title <- size_axis+2

xlab1 <- 'Time'
ylab1 <- 'Qobs'

df <- data.frame(seq(1,length(Qobs)),Qobs,limits[1,],limits[2,])
legend <- c("Qobs","Qsim_low","Qsim_high")
colnames(df) <- c("Time","Qobs","Qsim_low","Qsim_high")

color1 <- c("cornflowerblue","darkgreen","red")


title1 <- bquote(bold("Topmodel Simulation"))

sp3 <- ggplot(data=df, aes(x=Time)) + 
  ggtitle(title1) +
  xlab(xlab1) + ylab(ylab1) + 
  geom_point(aes(y = Qobs, color=legend[1]), size=2, shape = 1) +
  geom_line(aes(y = Qsim_low, color=legend[2]), size=1) + 
  geom_line(aes(y = Qsim_high, color=legend[3]), size=1) +
  
  
  scale_color_manual(name="Data",values = color1,limits=legend) +
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
fn3 <- paste0("Topmodel.pdf")
ggsave(filename=fn3, plot=sp3,width=50,height=20,units="cm", dpi=600)

##--------------------------------------------------------------------------