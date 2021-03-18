## Xin'AnJiang Hydrologic Model
## Modified from C++ code in 'Sibada/XAJ'
## (Errors exist in 'Sibada/XAJ': runoff routing by unit hydrograph)
## Gangsheng Wang (wang.gangsheng@gmail.com)
## March 15, 2021

XinAJrun <- function(PREC, EVAP, parameters, UH, Area, dt) {
  ## Area [km^2]: watershed area
  # /* **************************************************************************
  #   * Model parameters.
  # ************************************************************************* */
  nt <- length(PREC) ## Number of steps
  
  i <- 0
  j <- 0
  
  KC <- parameters[1]          ## 1.  Ratio of potential evap to pan evap
  IM <- parameters[2]          ## 2.  Fraction of impermeable area
  WUM <- parameters[3]         ## 3.  Soil moisture capacity of upper layer
  WLM <- parameters[4]         ## 4.  Soil moisture capacity of lower layer
  WDM <- parameters[5]         ## 5.  Soil moisture capacity of deep layer
  C <- parameters[6]           ## 6.  Coefficient of deep evap
  B <- parameters[7]           ## 7.  Exponent of the soil moisture storage capacity curve
  
  SM <- parameters[8]          ## 8.  Areal mean free water capacity of the surface soil layer
  EX <- parameters[9]          ## 9.  Exponent of the free water capacity curve
  KI <- parameters[10]          ## 10. outflow coefficients of the free water storage to groundwater
  KG <- parameters[11]         ## 11. outflow coefficients of the free water storage to interflow
  
  CI <- parameters[12]         ## 12. recession constant of groundwater storage
  CG <- parameters[13]         ## 13. recession constant of the lower interflow storage
  
  WM <- WUM + WLM + WDM        ## Mean water sotrage of the basin
  WMM <- WM * (1. + B) / (1. - IM)   ## Maximum water storage in the basin
  
  # /** *************************************************************************
  #   * Output and process variables.
  # ************************************************************************* */
  ##int nvars = 16
  ##NumericMatrix out(nt, nvars)
  ## Simulated Variabes
  E_s <- rep(0, nt)    ## Total ET (mm)
  EU_s <- rep(0, nt)   ## ET of upper soil layer
  EL_s <- rep(0, nt)   ## ET of lower soil layer
  ED_s <- rep(0, nt)   ## ET of deep soil layer
  W_s <- rep(0, nt)    ## Total soil water/moisture (mm)
  WU_s <- rep(0, nt)   ## soil water moisture of upper soil layer
  WL_s <- rep(0, nt)   ## soil water moisture of lower soil layer
  WD_s <- rep(0, nt)   ## soil water moisture of deep soil layer
  R_s <- rep(0, nt)     ## Total runoff (mm) of each time step
  RS_s <- rep(0, nt)    ## Surface runoff (mm) of each time step
  RI_s <- rep(0, nt)    ## Interflow (mm) of each time step
  RG_s <- rep(0, nt)    ## Underground runoff (mm) of each time step
  Q_s <- rep(0, nt)     ## Total runoff (m^3/s) at the outlet of the basin
  QS_s <- rep(0, nt)    ## Surface runoff (m^3/s) at the outlet of the basin
  QI_s <- rep(0, nt)    ## Interflow runoff (m^3/s) at the outlet of the basin
  QG_s <- rep(0, nt)    ## Underground runoff (m^3/s) at the outlet of the basin
  
  U <- Area / 3.6 / dt   ## Convert runoff from mm to m^3/s
  
  ## Model state variables for each step
  PE <- 0.    ## net prec when > 0  insu evap when < 0
  Ep <- 0.    ## Potential evapotranspiration, KC * EVAP[i]
  P <- 0.
  R <- 0.     ## Total runoff yield (mm)
  RB <- 0.    ## Runoff yield of the impermeable area
  RG <- 0.    ## Underground runoff yield
  RI <- 0.    ## Interflow runoff yield
  RS <- 0.    ## Surface runoff yield
  A <- 0.     ##
  E <- 0.     ## Total evap
  EU <- 0.     ## Evap at upper soil layer
  EL <- 0.     ## Evap at lower soil layer
  ED <- 0.     ## Evap at deep soil layer
  
  FR <- 0.
  AU <- 0.
  WU <- WUM * 0.8       ## Soil moisture (mm) of upper layer
  WL <- WLM * 0.8      ## Soil moisture (mm) of lower layer
  WD <- WDM * 0.8      ## Soil moisture (mm) of deep layer
  W <- WU + WL + WD    ## Soil moisture (mm) of all layer
  
  INUL <- 0.    ## Infiltration from upper layer to lower layer
  INLD <- 0.    ## Infiltration from lower layer to deep layer
  S <- SM * 0.2
  
  MS <- SM * (1. + EX)
  tmpA <- 0.   ## temporary variable in calculation of free water
  
  # /** ************************************************************************
  #   * Runoff yield
  # ************************************************************************* */
  for (i in 1:nt) {

    # /**   Evaporation   */
    RB <- PREC[i] * IM       ## RB: precipitation of the impermeable area
    P <- PREC[i] * (1. - IM)
    Ep <- KC * EVAP[i]
    
    if ((WU + P) >= Ep) {
      EU <- Ep
      EL <- 0
      ED <- 0
    } else if ((WU + P) < Ep) {
      EU <- WU + P
      ED <- 0.
      if (WL >= (C * WLM)) {
        EL <- (Ep - EU) * WL / WLM
      } else if (WL >= C * (Ep - EU)) {
        EL <- C * (Ep - EU)
      } else if (WL < C * (Ep - EU)) {
        EL <- WL
        ED <- C * (Ep - EU) - EL
      }
      if (ED > WD) {
        ED <- WD
      }
    } ## else if ((WU + P) < Ep)
    
    E <- EU + EL + ED
    PE <- P - E
    
    # /**  Infiltration and runoff yeild */
    if (PE <= 0) {
      R <- 0.00
      W <- W + PE
    } else {
      # /**  Use soil water capacity curve to cal change of soil water change */
      A <- WMM * (1 - '^'((1.0 - W / WM), 1. / (1. + B)))
      
      ## Depth of soil moisture + net prec < maximum soil water storage
      if ((A + PE) < WMM) {
        R <- PE - (WM - WM * '^'((1 - (PE + A) / WMM), (1 + B)) - W)
      } else {
        R <- PE - (WM - W)
      }
    }
    
    ## Infiltration
    if (WU + P - EU - R <= WUM) {
      WU = WU + (P - EU - R)
      WL = WL - EL
      WD = WD - ED
      INUL <- 0.
    } else {
      WU <- WUM
      INUL <- WU + P - EU - R - WUM
      if (WL - EL + INUL <= WLM) {
        WL <- WL - EL + INUL
        WD <- WD - ED
        INLD <- 0.
      } else {
        WL <- WLM
        INLD <- WL + INUL - EL - WLM
        if (WD - ED + INLD <= WDM) {
          WD <- WD - ED + INLD
        } else {
          WD <- WDM
        }
      }
    }
    
    W <- WU + WL + WD
    R_s[i] <- R
    
    # /**  Runoff source division */
    if (PE > 0.) {
      FR <- R / PE
      AU <- MS * (1 - '^'((1 - S / SM), 1 / (1 + EX)))
      if (PE + AU < MS) {
        tmpA <- SM * '^'((1 - (PE + AU) / MS), 1 + EX)
      } else {
        tmpA <- 0.
      }
      RS <- FR * (PE + S - SM + tmpA)
      S <- SM - tmpA
      RI <- FR * KI * S
      RG <- FR * KG * S
      S <- S * (1 - KI - KG)
      RS = RS + RB
    } else {
      FR <- 1 - '^'((1 - W / WM), B / (1 + B))
      RI <- FR * KI * S
      RG <- FR * KG * S
      S <- S * (1 - KI - KG)
      RS <- RB
    }
    ## Rcout << i << '\t' << FR<< '\t' << RS << '\t' << RI<<'\t' << RG << '\n'   ## !
    R <- RS + RI + RG
    
    # /** Save process variables */
    R_s[i] <- R
    RS_s[i] <- RS    ## Surface runoff (mm) of each time step
    RI_s[i] <- RI    ## Interflow (mm) of each time step
    RG_s[i] <- RG    ## Underground runoff (mm) of each time step
    
    E_s[i] <- E
    EU_s[i] <- EU
    EL_s[i] <- EL
    ED_s[i] <- ED
    W_s[i] <- W
    WU_s[i] <- WU
    WL_s[i] <- WL
    WD_s[i] <- WD
  } ## for (i in 1: nt)
  
  # /***************************************************************************
  #   * Routing
  # ************************************************************************* */
  
  ## Routing of surface runoff
  nUH <- length(UH)
  ## nconv
  for (i in 1:nt) {
    nconv <- nUH
    if (nUH > i) {
      nconv <- i
    }
    for (j in 1:nconv)
    {
      QS_s[i] = QS_s[i] + RS_s[i - j + 1] * UH[j] * U
    }
  }
  
  ## Routing of interflow and underground runoff
  for (i in 2:nt)	{
    QI_s[i] <- CI * QI_s[i - 1] + (1.0 - CI) * RI_s[i] * U
    QG_s[i] <- CG * QG_s[i - 1] + (1.0 - CG) * RG_s[i] * U
  }
  
  Q_s <- QS_s + QI_s + QG_s
  
  
  ## Rcout << "Modeling complete.\n"
  XinAJrun <-data.frame(
      Esim = E_s,
      EUsim = EU_s,
      ELsim = EL_s,
      EDsim = ED_s,
      Wsim = W_s,
      WUsim = WU_s,
      WLsim = WL_s,
      WDsim = WD_s,
      Rsim = R_s,
      RSsim = RS_s,
      RIsim = RI_s,
      RGsim = RG_s,
      Qsim = Q_s,
      QSsim = QS_s,
      QIsim = QI_s,
      QGsim = QG_s
    )
}

##==============================================================================

# Create IUH
IUH <- function(N, NK, len) {
  UH <- pgamma(seq(0, 100, length.out = len + 1),
               N, scale = NK)
  UH <- diff(UH)
  UH <- UH/sum(UH)
  UH
}
##==============================================================================

XinAJpar <- function() {
  ## KI + KG = 0.70
  param_names <- c('KC', 'IM', 'WUM', 'WLM', 'WDM', 'C', 'B', 'SM',  'EX', 'KI', 'KG', 'CI', 'CG', 'N', 'NK')
  param_value <- c(0.90, 0.03,  20.0,  60.0,  40.0, 0.15, 0.3, 30.0, 1.50, 0.40, 0.30, 0.80, 0.96, 4.0, 3.0)
  param_lower <- c(0.80, 0.00,   5.0,  10.0,  10.0, 0.05, 0.1, 10.0, 0.50, 0.01, 0.01, 0.50, 0.95, 0.1, 1.0)
  param_upper <- c(1.20, 0.05,  20.0,  90.0,  60.0, 0.20, 0.6, 60.0, 2.00, 0.70, 0.70, 0.90, 0.998,5.0, 6.0)

  XinAJpar <- data.frame(param_names = param_names, param_value = param_value, 
                         param_lower = param_lower, param_upper = param_upper)
  ## The lumped XAJ model has 13 parameters, including:
  #'
  #'          KC,   Ratio of potential evap to pan evap
  #'
  #'          IM,   Fraction of impermeable area
  #'
  #'          WUM,  Soil moisture capacity of upper layer
  #'
  #'          WLM,  Soil moisture capacity of lower layer
  #'
  #'          WDM,  Soil moisture capacity of deep layer
  #'
  #'          C,    Coefficient of deep evap
  #'
  #'          B,    Exponent of the soil moisture storage capacity curve
  #'
  #'          SM,   Areal mean free water capacity of the surface soil layer
  #'
  #'          EX,   Exponent of the free water capacity curve
  #'
  #'          KI,   outflow coefficients of the free water storage to interflow
  #'
  #'          KG,   outflow coefficients of the free water storage to groundwater
  #'
  #'          CI,   recession constant of the lower interflow storage
  #'
  #'          CG,   recession constant of groundwater storage
  #'
  #'          If use the instantaneous unit hydrograph (IUH) of Nash for routing
  #'          of surface runoff, should provided two other parameters:
  #'
  #'          N,    number of reservoirs in the instantaneous unit hydrograph
  #'
  #'          NK,   common storage coefficient in the instantaneous unit hydrograph
}


##-----------------------------------------------------------------------------
R2 <- function(obs,sim){
  obs_mean <- mean(obs)
  R2 <-  1-sum((sim-obs)^2)/sum((obs-obs_mean)^2)
}
##-----------------------------------------------------------------------------
SatFrac <- function(B, WMM, W){
  SatFrac <- 1-(1-W/WMM)**B
}
  
##-----------------------------------------------------------------------------
