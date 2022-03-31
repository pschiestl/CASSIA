PRELES_GPP <- function(photoparameters) {
  # Preles-parameters
  preles_parameters = c(
    413.0, ## 1 soildepth
    0.450, ## 2 ThetaFC
    0.118, ## 3 ThetaPWP
    3, ## 4 tauDrainage
    ## GPP_MODEL_PARAMETERS
    0.7457, ## 5 betaGPP
    10.93, ## 6 tauGPP
    -3.063, ## 7 S0GPP
    17.72, ## 8 SmaxGPP
    -0.1027, ## 9 kappaGPP
    0.03673, ## 10 gammaGPP
    0.7779, ## 11 soilthresGPP
    0.500, ## 12 b.CO2, cmCO2
    -0.364, ## 13 x.CO2, ckappaCO2
    ## EVAPOTRANSPIRATION_PARAMETERS
    0.2715, ## 14 betaET
    0.8351, ## 15 kappaET
    0.07348, ## 16 chiET
    0.9996, ## 17 soilthresET
    0.4428, ## 18 nu ET
    ## SNOW_RAIN_PARAMETERS
    1.2, ## 19 Meltcoef
    0.33, ## 20 I_0
    4.970496, ## 21 CWmax, i.e. max canopy water
    0, ## 22 SnowThreshold,
    0, ## 23 T_0,
    160, ## 24 SWinit, ## START INITIALISATION PARAMETERS
    0, ## 25 CWinit, ## Canopy water
    0, ## 26 SOGinit, ## Snow on Ground
    20, ## 27 Sinit ##CWmax
    -999, ## t0 fPheno_start_date_Tsum_accumulation; conif -999, for birch 57
    -999, ## tcrit, fPheno_start_date_Tsum_Tthreshold, 1.5 birch
    -999 ##tsumcrit, fPheno_budburst_Tsum, 134 birch
  )

  growth_photo_coef = 1
  if (photoparameters == 1) {
    PF = Rprebasso::PRELES(PAR = PAR, TAir = Temp, VPD = VPD, Precip = Rain, CO2=CO2, fAPAR=rep(0.85,365))$GPP
  }
  if (photoparameters == 2) {
    PF = Rprebasso::PRELES(PAR = PAR, TAir = Temp, VPD = VPD, Precip = Rain, CO2=CO2, fAPAR=rep(0.85,365), p = c(rep(NA,4), 0.8957, rep(NA, 4), 0.03473, rep(NA, 20)))$GPP
    growth_photo_coef = 1.24
  }
  if (photoparameters == 3) {
    PF = Rprebasso::PRELES(PAR = PAR, TAir = Temp, VPD = VPD, Precip = Rain, CO2=CO2, fAPAR=rep(0.85,365), p = c(rep(NA,4), 0.9357, rep(NA, 4), 0.03273, rep(NA, 20)))$GPP
    growth_photo_coef = 1.35
  }
  if (photoparameters == 4) {
    PF = Rprebasso::PRELES(PAR = PAR, TAir = Temp, VPD = VPD, Precip = Rain, CO2=CO2, fAPAR=rep(0.85,365), p = c(rep(NA,4), 0.9357, rep(NA, 4), 0.03273, rep(NA, 20)))$GPP
  }

  out <- list(PF, growth_photo_coef)
  names(out) = c("PF", "growth_photo_coef")

  return(out)
}


