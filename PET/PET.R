## Potential/Reference Evapotranspiration
## Penman-Monteith Equation 
## Author: Gangsheng Wang (wang.gangsheng@gmail.com)
## March 1, 2021

## Reference:
## (1) Monteith, J. L. 1965. Evaporation and environment. Symposia of the
## Society for Experimental Biology. 19, pp. 205-234. 
## (2) Zotarellli L, Dukes MD, Romero CC, Migliaccio KW, and Morgan KT. 
## Step by step calculation of the Penman-Monteith evapotranspiration (FAO-56 Method)
## University of Florida IFAS Extension

##------------------------------------------------------------------------------
SatVP <- function(T) {
  ## T [°C]: air temperature 
  ## SatVP [kPa]: saturation vapor pressure 
  SatVP <- 0.6108 * exp(17.27 * T / (T + 237.3))
  return(SatVP)
}

##------------------------------------------------------------------------------

VP <- function(Tmax, Tmin, RHmax, RHmin){
  ## VP [kPa]: actual vapor pressure
  ## esTmax and esTmin [kPa]: saturation vapor pressure at daily maximum and minimum temperature, respectively
  ## RHmax and RHmin [%]: maximum and minimum relative humidity, respectively
  esTmax <- SatVP(Tmax)
  esTmin <- SatVP(Tmin)
  VP <- (esTmin * RHmax / 100 + esTmax * RHmin / 100) / 2
}

##------------------------------------------------------------------------------
VPD <- function(Tmax, Tmin, RHmax, RHmin){
  ## VPD [kPa]: saturation vapor pressure deficient 
  ## esTmax and esTmin [kPa]: saturation vapor pressure at daily maximum and minimum temperature, respectively
  ## ea [kPa]: actural vapor pressure
  esTmax <- SatVP(Tmax)
  esTmin <- SatVP(Tmin)
  ea <- VP(Tmax, Tmin, RHmax, RHmin)
  VPD <- (esTmax + esTmin) / 2 - ea
}

##------------------------------------------------------------------------------
InvDistance <- function(DOY){
  ## InvDistance: inverse relative distance from Earth to Sun
  ## DOY: julian day
  InvDistance <- 1 + 0.033 * cos(2 * pi * DOY / 365)
}

##------------------------------------------------------------------------------
SolarDec <- function(DOY){
  ## SolarDec: solar declination
  ## DOY: julian day
  SolarDec <- 0.409 * sin(2 * pi * (DOY - 80) / 365)
}

##------------------------------------------------------------------------------
SunsetHourAngle <- function(Lat, DOY){
  ## SunsetHourangle [rad]: sunset hour angle
  ## Lat [m]: Latitude
  ## DOY: julian day
  
  ## convert from degrees to radians
  Lat_Rad <- Lat * pi / 180
  
  ## solar declination
  SolDec <- SolarDec(DOY)
  SunsetHourAngle <- acos(-tan(Lat_Rad) * tan(SolDec))
}

##------------------------------------------------------------------------------
PotRad <- function(Lat, DOY){
  ## PotRad [MJ m-2 d-1]: potential extraterrestrial radiation
  ## Lat: latitude in degrees
  ## DOY: Julian day, i.e., number of the day in the year between 1 and 365 or 366 
  
  Solar_Constant <- 118.08  ## [MJ m-2 d-1]
  ## Lat_Rad: convert latitude in degress to radians
  Lat_Rad <- Lat * pi / 180
  ## dr: inverse relative distance from Earth to Sun
  dr <- InvDistance(DOY) ##1 + 0.033 * cos(2 * pi * DOY / 365)
  ## SolDec: solar declination
  SolDec <- SolarDec(DOY) ##0.409 * sin(2 * pi * (DOY - 80) / 365)
  ## sun set hour angle
  SunsetAngle <- SunsetHourAngle(Lat, DOY)  ## acos(-tan(Lat_Rad) * tan(SolDec))

  Term <- SunsetAngle * sin(Lat_Rad) * sin(SolDec) + cos(Lat_Rad) * cos(SolDec) * sin(SunsetAngle)
  PotRad <- Solar_Constant * dr * Term / pi
}


##------------------------------------------------------------------------------
NetRad <- function(albedo, Lat, Elev, DOY, SolarRad, Tmax, Tmin, RHmax, RHmin){
  ## Calculate shortwave net radiation
  ## NetRad [MJ m-2 d-1]: net radiation
  ## ER [MJ m-2 d-1]: potential radiation
  ## Rs [MJ m-2 d-1]: incoming solar radiation
  ## ea [kPa]: actual vapor pressure
  ## albedo: albedo or canopy reflection coefficient, 
  ## albedo = 0.23 for the hypothetical grass reference crop 
  
  ER <- PotRad(Lat, DOY)
  Rs <- SolarRad
  ea <- VP(Tmax, Tmin, RHmax, RHmin)
  

  ## 
  Rns <- (1 - albedo) * Rs
  
  ## Calculate cloud factor
  Rs_ClearSky <- (0.75 + 2E-5*Elev)*ER
  F_Cloud <- 1.35 * (Rs /Rs_ClearSky) - 0.35
  
  ## Calculate humidity factor
  F_Hum <- (0.34 - 0.14 * sqrt(ea))
  
  ## Calculate Isothermal Longwave net radiation
  ## Stefan-Boltzmann constant
  Boltzmann <- 4.903E-09 ## [MJ K-4 m-2 d-1]
  ## longwave radiation
  LWR <-  Boltzmann* ((Tmax + 273) ^ 4 + (Tmin + 273) ^ 4) / 2
  ## Net outgoing long wave radiation
  Rnl <- LWR * F_Cloud * F_Hum
  
  ## Calculate Rn
  NetRad <- Rns - Rnl
}

##------------------------------------------------------------------------------
SolarRadiation <- function(as, bs, Lat, DOY, SunshineHour){
  ## when measured solar radiation data are NOT available, 
  ## Estimate Rs (incoming Solar radiation) by Sunshine hour data 
  
  ## Reference:
  ## Black JN, Bonython CW, Prescott JA. 
  ## Solar radiation and the duration of sunshine. 
  ## Quarterly Journal of the Royal Meteorological Society. 1954;80:231-235.

  
  ## SunshineHour: actual duration of sunshine received
  ## as, bs: regression coefficients, typical values: as = 0,23, bs = 0.5 
  ## Rs / Ra = as + bs*SunshineHour/N
  ## Rs: SolarRadiation, i.e., actual incoming solar radiation  
  ## Ra: potential radiation received if the atmosphere is perfectly transparent
  
  ## N: maximum possible duration of sunshine
  
  SunsetAngle <- SunsetHourAngle(Lat,DOY)
  Ra <- PotRad(Lat,DOY)
  N <- 24/pi * SunsetAngle
  SolarRadiation <- (as + bs * (SunshineHour/N)) * Ra
}
##------------------------------------------------------------------------------
WindSpeed2m <- function(WindSpeed,WindHeight){
  ## WindSpeed2m: Wind speed at 2 m above the ground level.
  ## WindSpeed: wind speed measured at a height (WindHeight)
  if(WindHeight == 2) {
    WindSpeed2m <- WindSpeed} else {
      WindSpeed2m <- WindSpeed * (4.87 / (log(67.8 * WindHeight - 5.42)))
    }
}

##------------------------------------------------------------------------------
AeroRes <- function(WindSpeed, WindHeight){
  ## AeroRes [d/m]: Aerodynamic resistance
  Uz <- WindSpeed
  z <- WindHeight
  U2 <- WindSpeed2m(Uz,z)
  U2 <- U2 * 86400 ##Convert to m/day
  
  d <- 0.08        ## zero plane displacement height [m]
  zom <- 0.01476   ## roughness length governing momentum transfer [m]
  zoh <- 0.001476  ## roughness length governing transfer of heat and vapor [m]
  zm <- 2          ## height of wind measurement [m]
  zh <- 2          ## height of humidity measurements [m]
  VK <- 0.41       ## von Karman's constant [-]
  Term1 <- log((zm - d) / zom)
  Term2 <- log((zh - d) / zoh)
  AeroRes <- Term1 * Term2 / (VK * VK * U2)
}
##------------------------------------------------------------------------------

DeltaEs <- function(Tmax, Tmin){
## DeltaEs: slope of saturation vapor pressure curve
  Tmean <- 0.5*(Tmax + Tmin)
  esTmean <- SatVP(Tmean)
  DeltaEs <- 4098.17 * esTmean / (Tmean + 237.3) ^ 2
}
##------------------------------------------------------------------------------
 
Lambda <- function(Tmax, Tmin){
  ## Lambda [MJ kg-1]: latent heat of vaporization
  Tmean <- 0.5*(Tmax + Tmin)
  Lambda <- 2.501 - 0.002361 * Tmean
}
# End Function
##------------------------------------------------------------------------------
Pressure <- function(Elev){
  ## Atmospheric pressure [kPa]
  ## Elev [m]: elevation / altitude
  Pressure <- 101.3*((293 - 0.0065*Elev)/293)^5.26  
}
##------------------------------------------------------------------------------
Gamma <- function(Elev,Tmax, Tmin){
  ## Gamma [kPa C-1]: psychrometric constant
  
  ## specific heat at constant pressure [MJ kg-1 C-1]
  Cp <- 0.001013
  ## ratio molecular weight of water vapor/dry air
  epsilon <- 0.622
  
  P <- Pressure(Elev)
  Lambda <- Lambda(Tmax, Tmin)
  Gamma <- Cp * P / (epsilon * Lambda)
}


# End Function
##------------------------------------------------------------------------------
ETRadTerm <- function(albedo, Lat, Elev,
                      DOY, SolarRad, Tmax, Tmin, 
                      RHmax, RHmin, WindSpeed, WindHeight){
  ## ETRadTerm [mm d-1]: radiation term
  ## surface resistance
  rs <- 70/86400 ## [d m-1]
  ## aerodynamic resistance
  ra <- AeroRes(WindSpeed, WindHeight)
  
  Delta <- DeltaEs(Tmax, Tmin)
  Rn <- NetRad(albedo, Lat, Elev, DOY, SolarRad, Tmax, Tmin, RHmax, RHmin)
  Lambda1 <- Lambda(Tmax, Tmin)
  Gamma1 <- Gamma(Elev,Tmax, Tmin)
  
  ETRadTerm <- Delta * Rn / (Delta + Gamma1 * (1 + rs / ra))
  ETRadTerm <- ETRadTerm / Lambda1
}
##  End Function
##------------------------------------------------------------------------------
ETAeroTerm <- function(Elev, Tmax, Tmin, RHmax, RHmin, 
                       WindSpeed, WindHeight){
  ## ETAeroTerm [mm d-1]: aerodynamic or wind term
  ## specific heat at constant pressure [MJ kg-1 C-1]
  Cp <- 0.001013
  ## surface resistance
  rs <- 70/86400 ## [d m-1]
  ## aerodynamic resistance
  ra <- AeroRes(WindSpeed, WindHeight)
  
  Delta <- DeltaEs(Tmax, Tmin)
  Lambda1 <- Lambda(Tmax, Tmin)
  Gamma1 <- Gamma(Elev,Tmax, Tmin)
  VPD1 <- VPD(Tmax, Tmin, RHmax, RHmin)
  
  P <- Pressure(Elev)
  
  ## virtual temperature [K]
  Tkv <- 1.01 * (0.5*(Tmax+Tmin) + 273)
  AirDensity <- 3.486 * P / Tkv
  VolHeatCap <- Cp * AirDensity
  ETAeroTerm <- (VolHeatCap * VPD1 / ra) / (Delta + Gamma1 * (1 + rs / ra))
  ETAeroTerm <- ETAeroTerm / Lambda1
}
# End Function
##------------------------------------------------------------------------------
ETPenmanMonteith <- function(albedo, Lat, Elev, 
                             DOY, SolarRad, Tmax, Tmin, 
                             RHmax, RHmin, WindSpeed, WindHeight){
  ## Penman-Monteith Reference ET
  ETRad <- ETRadTerm(albedo, Lat, Elev, 
                     DOY, SolarRad, Tmax, Tmin, 
                     RHmax, RHmin, WindSpeed, WindHeight)
  
  ETAero <- ETAeroTerm(Elev, Tmax, Tmin, 
                       RHmax, RHmin, WindSpeed, WindHeight)
  ETPenmanMonteith <- ETRad + ETAero
}


##------------------------------------------------------------------------------
ETHargreaves <- function(Lat, DOY, Tmax, Tmin){
  ## Hargreaves Method
  
  ## Hargreaves GH, Samani ZA. Reference crop evapotranspiration from temperature. 
  ## Applied engineering in agriculture. 1985;1:96-99.
  Tmean <- 0.5 * (Tmax + Tmin)
  Lambda1 <- Lambda(Tmax, Tmin)
  Pot_Rad <- PotRad(Lat, DOY)
  ETHargreaves <- 0.0023 * (Tmean + 17.8) * (Tmax - Tmin) ^ 0.5 * Pot_Rad / Lambda1
}
# End Function

##------------------------------------------------------------------------------
ETPriestleyTaylor <- function(Alfa, albedo, Lat, Elev, DOY, SolarRad, Tmax, Tmin, RHmax, RHmin){
  ## Priestley Taylor Method
  ## Alfa ranging from 1.26 to 1.70
  
  ## Akumaga, U. and Alderman, P.D. (2019), 
  ## Comparison of Penman–Monteith and Priestley‐Taylor Evapotranspiration Methods for Crop Modeling in Oklahoma. 
  ## Agronomy Journal, 111: 1171-1180. https://doi.org/10.2134/agronj2018.10.0694
  Delta1 <- DeltaEs(Tmax, Tmin)
  Lambda1 <- Lambda(Tmax, Tmin)
  Gamma1 <- Gamma(Elev,Tmax, Tmin)

  Rn <- NetRad(albedo, Lat, Elev,DOY, SolarRad, Tmax, Tmin, RHmax, RHmin)

  ETPriestleyTaylor <- Alfa * Delta1 * Rn / (Delta1 + Gamma1)
  ETPriestleyTaylor <- ETPriestleyTaylor / Lambda1
}
##------------------------------------------------------------------------------


R2 <- function(obs,sim){
  obs_mean <- mean(obs)
  R2 <-  1-sum((sim-obs)^2)/sum((obs-obs_mean)^2)
##-----------------------------------------------------------------------------
}