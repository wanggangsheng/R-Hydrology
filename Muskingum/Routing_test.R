## River flow routing
## Muskingum Method

## Author: Gangsheng Wang (wang.gangsheng@gmail.com)
## March 1, 2021

library(ggplot2)
library(ggpubr) #ggscatter, stat_cor
library(scales) #as.POSIXct, date_format
library(Metrics) ## rmse
#Clean up Environment
rm(list = ls())

## Please change the following to your own Working Directory
setwd("/Users/wgs/ownCloud/Rcode/R-Hydrology/Muskingum")
##-----------------------------------------------------------------------------
Muskingum_Routing <- function(dt, K, X, I1, I2, O1){
  C0 <- dt/K + 2*(1-X)
  C1 <- (dt/K + 2*X)/C0
  C2 <- (dt/K - 2*X)/C0
  C3 <- 1 - C1 -C2
  O2 <- C1*I1 + C2*I2 + C3*O1
  Muskingum_Routing <- list(C0=C0, coef=c(C1,C2,C3),Outflow2=O2)
} 
##-----------------------------------------------------------------------------
R2 <- function(obs,sim){
  obs_mean <- mean(obs)
  R2 <-  1-sum((sim-obs)^2)/sum((obs-obs_mean)^2)
}
  
##-----------------------------------------------------------------------------
## dt: time interval
## K: storage constant
## X: weighting factor

## O2 = C1*I1 + C2*I2 + C3*O1
##-----------------------------------------------------------------------------
## derive K and X from C1/C2/C3 
## Let Y = dt/K
## C1 + C2 = 2Y/C0
## C1 - C2 = 4X/C0
## a = Y/(2X) = (C1 + C2)/(C1 - C2)

## X = C1/[a + 1 - C1*(a - 1)]
## Y = 2aX
## K = dt/Y
##-----------------------------------------------------------------------------
path <- getwd()

file0 <- "Routing_input.txt"
file0 <- paste0(path,"/", file0)

## read input data from file0
df0 <- read.table(file0, sep="\t",header=TRUE)
names(df0)
colnames(df0)
head(df0)


dt <- 1
Xlimit <- c(0, 0.5)

nt <- nrow(df0)
ncol(df0)
df <- df0[2:nt,]
colnames(df) <- c("Time","Inflow2", "Outflow2")
df$Inflow1 <- df0$Inflow[1:(nt-1)]
df$Outflow1 <- df0$Outflow[1:(nt-1)]
## order the columns
df <- df[,c("Time","Inflow1","Inflow2","Outflow1","Outflow2")]
## add random noise 
## df$Outflow2 <- df$Outflow2 + runif(nt-1, min=-20, max=20)  
##-----------------------------------------------------------------------------
## Method 1: use linear model to derive 3 coefficient C1,C2,C3
## Attention: Method 1 does NOT consider the constraint: C1+C2+C3 = 1
model <- lm(Outflow2 ~ 0 + Inflow1 + Inflow2 + Outflow1, df)
summary(model)
coef <- as.numeric(model$coefficients)
## you will see sum(coef) is NOT = 1 exactly
coef
sum(coef)
##-----------------------------------------------------------------------------
## Method 2: use linear model to derive 2 coefficient C1,C2
## Consider the constraint: C1+C2+C3 = 1 
## O2 = C1*I1 + C2*I2 + C3*O1
## O2 = C1*I1 + C2*I2 + (1-C1-C2)*O1
## O2-O1 = C1*(I1-O1) + C2*(I2-O1)

df$O2O1 <- df$Outflow2 - df$Outflow1
df$I1O1 <- df$Inflow1 - df$Outflow1
df$I2O1 <- df$Inflow2 - df$Outflow1
model <- lm(O2O1 ~ 0 + I1O1 + I2O1, df)
summary(model)
model$coefficients
coef2 <- as.numeric(model$coefficients)
coef <- c(coef2,1-sum(coef2))
coef
sum(coef)
##-----------------------------------------------------------------------------
a <- (coef[1] + coef[2])/(coef[1] - coef[2])
X <- coef[1]/(a+1-coef[1]*(a-1))
Y <- 2*a*X
K <- dt/Y

print(paste0("C1+C2+C3 = ",sum(coef)))
print("C1, C2, C3 = ")
print(coef)
print(paste0("X = ",X,"; K = ",K))

##-----------------------------------------------------------------------------
## validate the Muskingum Method
results <- Muskingum_Routing(dt, K, X, df$Inflow1, df$Inflow2, df$Outflow1)
df$Outflow2_sim <- results$Outflow2

results$coef - coef

cor_Model <- cor(df$Outflow2,df$Outflow2_sim)
R2_Model <- R2(df$Outflow2,df$Outflow2_sim)
RMSE_Model <- rmse(df$Outflow2,df$Outflow2_sim)
RMSE2 <- sqrt(sum((df$Outflow2 - df$Outflow2_sim)^2)/nrow(df))
print(paste0("Pearson Correlation Coefficient           (r) = ",cor_Model))
print(paste0("Squared Pearson Correlation Coefficient (r^2) = ",cor_Model^2))
print(paste0("Coefficient of Determination            (R^2) = ",R2_Model))
print(paste0("Root Mean Square Error (RMSE) = ",RMSE_Model, "; ",RMSE2))
##-----------------------------------------------------------------------------
## Time-Series Plot
xlab1 <- 'Time (h)'
ylab1 <- bquote(bold("Discharge ") * bold(" ("*m^bold(3) * s ^ bold("-1") * ")"))

legend <- c("Inflow", "Observed Outflow", "Simulated Outflow")
color1 <- c("cornflowerblue", "darkgreen", "red")

size_axis <- 20
size_title <- size_axis+2
size_regression <- size_axis/3


title1 <- bquote(bold("Muskingum Flow Routing"))

sp3 <- ggplot(data = df, aes(x = df$Time)) +
  ggtitle(title1) +
  xlab(xlab1) + ylab(ylab1) +
  geom_line(aes(y = Inflow2, color = legend[1]), size = 1) +
  geom_line(aes(y = Outflow2, color = legend[2]), size = 1) +
  geom_line(aes(y = Outflow2_sim, color = legend[3]), size = 1) +
  
  geom_point(aes(y = Inflow2, color = legend[1]), size=4, shape=21) + 
  geom_point(aes(y = Outflow2, color = legend[2]), size=6, shape=21) + 
  geom_point(aes(y = Outflow2_sim, color = legend[3]), size=4, shape=24) + 
  
  # scale_x_datetime(date_break = "2 month", labels = date_format("%m/%d/%y")) +
  scale_color_manual(name = "", values = color1) +
  guides(col = guide_legend(nrow = 3)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = size_title, face = "bold"),
    plot.subtitle = element_text(
      hjust = -0.03,
      size = size_title,
      face = "bold"
    ),
    strip.text.x = element_blank(),
    axis.text = element_text(size = size_axis),
    axis.title = element_text(size = size_title, face = "bold"),
    legend.text = element_text(size = size_axis),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.8)
  )

print(sp3)
fn3 <- paste0(path, "/Muskingum_Routing.pdf")
ggsave(
  filename = fn3,
  plot = sp3,
  width = 40,
  height = 20,
  units = "cm",
  dpi = 600
)


##-----------------------------------------------------------------------------


##==============================================================================

##==============================================================================

