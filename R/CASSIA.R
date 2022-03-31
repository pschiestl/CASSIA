CASSIA <- function(
  #####
  ## Weather Inputs - input in a dataframe with date, temperature, Photosynthesis, soil temperature a and b horizon, soil moisture and precipitation
  #####
  # TODO: is this a sensible input?

  weather,

  #####
  ## Site
  #####

  site,

  #####
  ## Parameters
  #####
  ratios = ratios_p,
  parameters = parameters_p,
  common = common_p,
  sperling = sperling_p,
  repo = repo_p,

  #####
  ## Default values of the set up
  #####

  storage.reset = TRUE,			# storage.reset<-TRUE=Same initial storage each year, storage.reset<-False, The storage on the last day of year X is  postponded to the first day of the year X+1
  storage.grows = FALSE,			# TRUE if the critical storage level increases with tree size.

  LN.estim = TRUE,				# LN depends on the GPP during previous july-august
  mN.varies = TRUE,				# needle mass (in maintenance respiration) is 2/3 of the total during period 1.10 - 31.5.

  LD.estim = TRUE,				# LD depends on the GPP during March-August
  sD.estim.T.count = FALSE,			# sD depends on the number of days when g in growing window - analogue to needles

  LH.estim = TRUE,

  trees_grow = FALSE,				# can be false if mature trees are modelled and not for a very long period
  growth_decreases = FALSE,			# the height and diameter growth (alfa_S and alfaD) decrease during the simulation
  needle_mass_grows = FALSE,		# Is needle mass dynamic i.e. the modelled growth is also respiring etc and following for some years? If true, note that root mass is related to needle mass

  mychorrhiza = TRUE, 			# If allocation to mychorrhiza is taken into account
  root_as_Ding = TRUE,

  sperling_model = FALSE,       # Dynamic sugar model using Sperling's enzyme dynamics
  xylogenesis = FALSE,

  PRELES_GPP = FALSE,
  environment_effect_xylogenesis = FALSE,

  photoparameters = 3,
  temp_rise = FALSE,
  drought = FALSE,
  Rm_acclimation = TRUE,

  s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
  s.H0 = 1					# and for shoot grwoth
  ) {

  #####
  ## Input tests!
  #####
  # weather
  if (ncol(weather) != 7) {stop("Incomplete weather data: variable number incorrect")}
  if (sum(names(weather) == c("date", "T", "P", "TSA", "TSB", "MB", "Rain")) != 7) {stop("Incomplete weather data - incorrect variables, or named incorrectly")}

  # Check that the sites are within the sites allowed
  if ((site %in% c("Hyde", "Lettosuo", "Flakaliden_c")) == F) {stop("Unknown site: Please pick between Hyde and Lettosuo")}
  # TODO: Flakaliden_c is a possible site, but J doesn't have any more information on it - need to check with Pauliina

  if (sperling_model == TRUE) {if (mychorrhiza == T) {
    mychorrhiza = FALSE
    warning("Mycorrhiza has been changed to mycorrhiza = false as mycorrhiza is included explicitly in the Sperling submodel")
    }
  }

  if (xylogenesis == TRUE) {
    # TODO: this is not set in stone, but fits the initial setup of Lettosuo
    warning("As xylogenesis is TRUE, LN.estim and trees_grow set to FALSE and TRUE respectively")
    LN.estim = FALSE   # LN depends on the GPP during previous july-august
    trees_grow = TRUE  # can be false if mature trees are modelled and not for a very long period
    mycorrhiza = FALSE
  }

  #if (nrow(sperling) != nrow(sperling_p)) {stop("Sperling input is the wrong size!")}
  #if (nrow(parameters) != nrow(parameters_p)) {stop("Parameters input is the wrong size!")}
  #if (nrow(ratios) != nrow(ratios_p)) {stop("Ratios input is the wrong size!")}
  #if (length(common) != length(common_p)) {stop("Common input is the wrong size!")}
  #if (length(repo) != length(repo_p)) {stop("Common input is the wrong size!")}
  #if (rownames(sperling) != names(sperling_p)) {stop("Sperling has the wrong row names")}
  #if (rownames(parameters) != rownames(parameters_p)) {stop("Parameters have the wrong row names")}
  #if (rownames(ratios) != rownames(ratios_p)) {stop("Ratios have the wrong row names")}
  #if (names(common) != names(common_p)) {stop("Ratios have the wrong row names")}
  #if (names(repo) != names(repo_p)) {stop("Ratios have the wrong row names")}

  #####
  ## Model conditions derived from model inputs
  #####
  # years from weather data
  years <- seq(as.numeric(substring(weather$date[end(weather$date)], 1, 4))[2],
               as.numeric(substring(weather$date[end(weather$date)], 1, 4))[1])
  if (sum(is.na(years)) != 0) {stop("NA in years vector: Check the date format in weather input")}

  #####
  ## Creating the output vectors
  #####

  total.days <- 366 * sum(years %in% leap_years) + 365 * (length(years) - sum(years %in% leap_years))

  if (sperling_model == T) {
    export_yearly <- data.frame(matrix(ncol=26, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=34, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts",
                              "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                              "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot",
                             "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                             "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
  } else if (xylogenesis == T) {
    export_yearly <- data.frame(matrix(ncol=20, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=27, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "ring_width", "needle_mass", "sum.needle.cohorts",
                              "ew_width", "ring_density", "ew.cells_tot", "lw.cells_tot")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "daily.consumption", "ring_width", "GD.tot", "n.E.tot", "n.W.tot", "n.M.tot")
  } else {
    export_yearly <- data.frame(matrix(ncol=16, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=24, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot")
  }

  n.days.export <- 0
  n.year <- 1

  #####
  ## Non yearly dependent coefficients
  #####
  # height growth coefficient

  height_growth_coefficient <- diameter_growth_coefficient <- NULL
  height_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("height_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))
  diameter_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("diameter_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))

  if (growth_decreases == TRUE) {
    height_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("height_growth_coefficient_max"),c(site)], ratios[c("height_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
    diameter_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("diameter_growth_coefficient_max"),c(site)], ratios[c("diameter_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
  }

  if (xylogenesis == TRUE) {
    cell.l<-(parameters[c("cell.l.ew"),c(site)]+parameters[c("cell.l.lw"),c(site)])/2
    cell.d<-(parameters[c("cell.d.ew"),c(site)]+parameters[c("cell.d.lw"),c(site)])/2

    cell.volume.ew<-(parameters[c("cell.d.ew"),c(site)])^2*cell.l
    cell.volume.lw<-(parameters[c("cell.d.lw"),c(site)])^2*cell.l

    wall.volume.ew<-cell.volume.ew-(parameters[c("cell.d.ew"),c(site)]-2*parameters[c("wall.thickness.ew"),c(site)])^2*(cell.l-2*parameters[c("wall.thickness.ew"),c(site)])	# m3
    wall.volume.lw<-cell.volume.lw-(parameters[c("cell.d.lw"),c(site)]-2*parameters[c("wall.thickness.lw"),c(site)])^2*(cell.l-2*parameters[c("wall.thickness.lw"),c(site)])	# m3

    cell.wall.volume.growth.per.day.ew <- wall.volume.ew / parameters[c("tau.We"),c(site)]
    cell.wall.volume.growth.per.day.lw <- wall.volume.lw / parameters[c("tau.Wl"),c(site)]
  } else {
    cell.volume.ew<-(parameters[c("cell.d.ew"), c(site)])^2*parameters[c("cell.l.ew"),c(site)]
    cell.volume.lw<-(parameters[c("cell.d.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)])*parameters[c("cell.l.lw"),c(site)]		# the tangential width of the cell is the same for early and late wood

    wall.volume.ew<-cell.volume.ew-(parameters[c("cell.d.ew"), c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])^2*(parameters[c("cell.l.ew"),c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])	# m3
    wall.volume.lw<-cell.volume.lw-(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.l.lw"),c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])	# m3
  }

  if (environment_effect_xylogenesis == FALSE) {
    cell.density.ew<-parameters[c("cell.wall.density.ew"),c(site)]*wall.volume.ew/cell.volume.ew		# to calculate the wood density
    cell.density.lw<-parameters[c("cell.wall.density.lw"),c(site)]*wall.volume.lw/cell.volume.lw

    daily.rate.ew<-wall.volume.ew/parameters[c("tau.We"),c(site)]		#m3 cell wall day-1
    daily.rate.lw<-wall.volume.lw/parameters[c("tau.Wl"),c(site)]		#m3 cell wall day-1

    Carbon.daily.rate.ew<-daily.rate.ew*parameters[c("cell.wall.density.ew"),c(site)]	# kg C day-1	carbon to early wood wall formation
    Carbon.daily.rate.lw<-daily.rate.lw*parameters[c("cell.wall.density.lw"),c(site)]	# kg C day-1	carbon to late wood wall formation
  }

  if (site == "Hyde") {
    stem.no <- cbind(1997 : 2020, rep(1010, length.out = length(1997 : 2020)))	# Photosynthesis calculated with SPP-model for tree class 15-20 cm. Then compared with eddy GPP, determined with which stem no. the portion of eddy GPP is same as SPP estimate.
  } else if (site == "Lettosuo") {
    stem.no <- cbind(1997 : 2020, rep(500, length.out = length(1997 : 2020)))	# Photosynthesis calculated with PRELES or SPP-model.
  }

  # Carbon to height growth
  CH<-parameters[c("density_tree"),c(site)]*parameters[c("carbon_share"),c(site)]
  M.suc<-12*common[[c("M.C")]]+22*common[c("M.H")]+11*common[[c("M.O")]]

  #####
  ## Year loop
  #####

  LAI <- needle_mass <- NULL
  count <- 1

  for (year in if (sperling_model == F) {years} else {rep(years, each = 2)}) {

    n.days <- if (year %in% leap_years) 366 else 365

    if (n.year == 1) {
      repol = repola(parameters[c("D0"), c(site)], parameters[c("h0"), c(site)], n.year, needle_mas = NULL, ste = site, params = parameters, reps = repo)
      needle_mass[1] <- repol[[c("needle_mass")]]
    } else {
      repol = repola(parameters[c("D0"), c(site)], parameters[c("h0"), c(site)], n.year, needle_mas = needle_mass, ste = site, params = parameters, reps = repo)
      if (needle_mass_grows==FALSE) {
        needle_mass[n.year]=needle_mass[n.year-1] 		# constant needle mass, needlemass determined in sitespecific parameters (e.g. parameters_hyde.r)
      } else {
        needle_mass <- repol[[c("needle_mass")]] # calculated using biomass equations, reached height and diameter
      }
    }

    # Needle age classes: oldest fall off, the others age with one year
    # The youngest needle class is the (needle_length[year-]/average_length) * 1/3*(estimated needle mass by the biomass equations)
    needle_cohorts <- matrix(ncol = parameters[c("n_age"), c(site)], nrow = length(years))				# Needle age classes (assumed three classes as in central Finland)

    if (xylogenesis == TRUE) {
      if (n.year > 1) {
        needle_cohorts[n.year, 1] <- needle.tot[365]
        for (i in 2 : n_age) {
          needle_cohorts[n.year, i] <- needle_cohorts[n.year-1, i-1]
        }
      }
    } else {
      if (n.year > 1) {
        needle_cohorts[n.year, 1] <- cum.Daily.N.tot[365] / parameters[c("n_lenght"), c(site)] * needle_mass[n.year] / 3
        needle_cohorts[n.year, 2] <- needle_cohorts[n.year-1, 1]
        needle_cohorts[n.year, 3] <- needle_cohorts[n.year-1, 2]
      }
    }

    # Not used for anything and something is wrong - don´t use the LAI estimates unless corrected!
    LAI[n.year] <- stem.no[which(stem.no[,1] == year), 2] * sum(needle_cohorts[n.year,]) * parameters[c("SLA"), c(site)] / 10000

    ## Weather inputs for the year are in vector form

    Temp <- PF <- Tsa <- Tsb <- M.soil <- Rain <- NULL
    Temp <- weather[substring(weather$date, 1, 4) == year, c("T")]		    # Temperature, C
    PF <- weather[substring(weather$date, 1, 4) == year, c("P")]            # g C m-2 day-1
    Tsa <- weather[substring(weather$date, 1, 4) == year, c("TSA")]			    # Soil temperature in A-horizon, C
    Tsb <- weather[substring(weather$date, 1, 4) == year, c("TSB")]			    # Soil temperature in B-horizon, C
    M.soil <- weather[substring(weather$date, 1, 4) == year, c("MB")]		# Soil moisture (m3 /m3) in B-horizon
    Rain <- weather[substring(weather$date, 1, 4) == year, c("Rain")]		    # mm day-1

    # CO2, VPD and PAR preles

    growth_photo_coef <- 1
    if (PRELES_GPP == TRUE) {
      growth_photo_coef = PRELES_GPP(photoparameters, Temp, PF, Tsa, Tsb, M.soil, Rain)
    }

    # Initalising the basic values for these variables
    B0<-pi/4*parameters[c("D0"), c(site)]^2		# basal area in the beginning
    h00<-parameters[c("h0"), c(site)]
    D00 <- parameters[c("D0"), c(site)]
    if (xylogenesis == TRUE) {
      parameters[c("LH0"),c(site)] <- h_increment/(0.5*sHc)
      parameters[c("LN0"),c(site)] <- n_lenght/(0.5*sNc)
      parameters[c("LR0"),c(site)] <- 2 * m.R.tot / sRc
    }

    ### Photosynthesis
    # Photosynthesis of one tree per day (g C / m2 --> kg C / tree)
    P <- tot.P <- NULL

    P <- PF / stem.no[which(stem.no[,1] == year), 2] * 10000 / 1000

    # Total photosynthesis (kg C)
    tot.P <- cumsum(P)

    #################  Growth potential #######################
    ## g
    g <- NULL
    g <- ((Temp > 0) * (1 / (1 + exp(-common[[c("a")]] * (Temp - common[[c("b")]])))))
    g[is.na(g)] <- 0

    ## Length growth
    g.sH <- sH <- height.pot.growth <- GH <- HH <- fH <- height.pot <- NULL
    for (yy in 1 : n.days) if (yy < s.H0) g.sH[yy] <- 0 else (g.sH[yy] <- g[yy])
    sH <- parameters[c("sH0"), c(site)] + cumsum(g.sH)										# The phase of the annual cycle of length growth, sH:

    fH <- (sH > 0) * (sH < parameters[c("sHc"), c(site)]) * (sin(2 * pi / parameters[c("sHc"), c(site)] * (sH - parameters[c("sHc"), c(site)] / 4)) + 1) / 2		# A function driven by the phase of the annual cycle [0,1]

    LH <- parameters[c("LH0"), c(site)] * height_growth_coefficient[which(height_growth_coefficient[,1] == year), 2]
    if (LH.estim == TRUE) LH <- LH * GPP_previous_sum[which(GPP_previous_sum[,1] == year), 2] / mean(GPP_previous_sum[,2])
    LH = LH * growth_photo_coef

    # Daily potential length growth (mm/day)  and length (mm)
    GH <- g.sH * fH * LH
    HH <- parameters[c("HH0"),c(site)] + cumsum(GH)

    # Daily use of carbon to length growth (m2 * kg C m-3 day-1 * = kg C day-1)
    height.pot.growth <- B0 * CH * GH / 1000 * ratios[c("form_factor"),c(site)]

    ## Needle growth

    sN <- fN <- GN <- HN <- needle.pot.growth <- needle.pot <- NULL

    # Similarly as shoots (fN same as with secondary growth)

    sN <- parameters[c("sN0"), c(site)] + cumsum(g.sH)								# The phase of the annual cycle of length growth, sN:

    fN <- (sN > 0) * (sN < parameters[c("sNc"), c(site)]^2) * (parameters[c("sNc"), c(site)] * sN^(1 / 2) - sN) / (parameters[c("sNc"), c(site)]^2 / 4)	# A function driven by the phase of the annual cycle [0,1] (annual pattern of growth)
    fN[is.na(fN)] <- 0

    if (LN.estim==TRUE) LN <- parameters[c("LN0"), c(site)] * GPP_previous_sum[which(GPP_previous_sum[,1] == year), 2] / mean(GPP_previous_sum[,2])
    LN = LN * growth_photo_coef

    # Daily potential length growth (mm/day)  and length (mm)  of needles
    GN <- g * fN * LN
    HN <- parameters[c("HN0"), c(site)] + cumsum(GN)

    # Carbon to needle growth per day (if there's no carbon limitation) (kg C / day)
    needle.pot.growth <- repol[[c("m.N")]] * GN * (max(HH) / parameters[c("h_increment"), c(site)])		# Simulated shoot length/average shoot length

    ## Root growth
    sR <- fR <- GR <- root.pot.growth <- NULL

    if (site == "Hyde") {
      LR <- LR0 <- 2 * repol[[c("m.R.tot")]] / parameters[c("sRc"), c(site)] # root growth rate
    } else {
      LR <- parameters[c("LR0"), c(site)] / parameters[c("root.lifetime"),c(site)] * growth_photo_coef
    }

    if (root_as_Ding == TRUE) {
      fib_coef = 0.25    # 0.25, 2.3  # Determines the proportion of fibrous roots (1 leads to 37 % of fibrous roots, 0.25 to 13 % of fibrous roots and 2.3 to 63 % of fibrous roots)
      sR <- cumsum((1 : n.days) >= 150)   # Growth begins earliest in the end of May
      fR <- ifelse(sR == 0, 0, 1/(1+exp(-0.038014*(sR-56.06243))))      # function determines by stage of development [0,1] (parameters from excel file)
      fR[(1 : n.days) >= 320] = 0
      gR_fib = -0.84 + 0.13 * Tsa -0.44 + 2.11 * M.soil     # growth of fibrous roots from Ding et al. 2019 (model 5)
      gR_fib[gR_fib<0] = 0
      gR_pio = -0.84 + 0.13 * Tsb + 0.32 -0.16 + 0.78 * M.soil    # growth of pioneer roots from Ding et al. 2019 (model 5)
      gR_pio[gR_pio<0] = 0
      gR <- fib_coef * gR_fib + gR_pio    # if fib_coef = 1 this leads in year 2018 to 37 % fibrous roots of all roots
      LR <- 0.0049 * 1 / (0.37 * fib_coef + 0.63)  # if fib_coef = 1 this leads to (roughly and on average) same total root growth as original. If fib_coef is changed, L is changed accordingly

      # Carbon to root growth per day (if there's no carbon limitation) (kg C / day)
      GR <- fR * LR * gR
    } else { # TODO: check conditions
      if (xylogenesis == FALSE) {
        gR <- (Tsb > common[[c("TR0")]]) * (1 / (1 + exp(-common[[c("a")]] * ((Tsb - common[[c("TR0")]]) - common[[c("b")]])))) * (1 - 1 / exp(M.soil * 10)) # Temp and M driving the phase of the annual cycle of root growth
      } else {
        soil_moisture_effect <- 1 # TODO: Find the original value
        gR <- (Tsb > TR0) * (1 / (1 + exp(-a * ((Tsb - TR0) - b)))) * soil_moisture_effect # Temp and M driving the phase of the annual cycle of root growth
      }
      sR <- parameters[c("sR0"), c(site)] + cumsum(gR)										# The phase of the annual cycle of root growth
      fR <- (sR < parameters[c("sRc"), c(site)]) * (sR > 0) * (sin(2 * pi / parameters[c("sRc"), c(site)] * (sR - parameters[c("sRc"), c(site)] / 4)) + 1) / 2	# A function driven by the phase of the annual cycle (annual pattern of growth) [0,1]


      # Carbon to root growth per day (if there's no carbon limitation) (kg C / day)
      GR <- fR * LR * gR

    }

    root.pot.growth <- GR								# Carbon to root growth per day (if there?s no carbon limitation) (kg C / day)
    root.pot.growth[is.na(root.pot.growth)] <- 0

    ## Diameter growth
    g.sD.T <- g.sD.GPP <- NULL
    sDA <- sD <- fD <- GD <- tot.cells <- NULL

    for (yy in 1:n.days) if(yy<s.D0) g.sD.T[yy] <- 0 else (g.sD.T[yy] <- g[yy])

    sD <- parameters[c("sD0.Trad"), c(site)] + cumsum(g.sD.T)						# The phase of the annual cycle of diameter growth, sD:

    if (sD.estim.T.count == TRUE){
      sDA <- if (xylogenesis == TRUE) sD0.T.count + cumsum(g.sD.T) else sD
      sD <- cumsum(sDA > 0)
      if (xylogenesis == TRUE) {parameters[c("sDc"), c(site)] <- sDc.T.count}
    }

    fD <- (sD > 0) * (sD < parameters[c("sDc"), c(site)]^2) * (parameters[c("sDc"), c(site)] * sD^(1/2) - sD)/(parameters[c("sDc"), c(site)]^2 / 4) 	# A function driven by the phase of the annual cycle (annual pattern of growth) [0,1]
    fD[is.na(fD)] <- 0

    LD <- parameters[c("LD0"), c(site)] * diameter_growth_coefficient[which(diameter_growth_coefficient[,1] == year), 2]
    if (LD.estim == TRUE) {
      S.GPP <- S.GPP_ref <- dS.GPP <- dS.GPP_ref <- NULL
      LD <- rep(parameters[c("LD0"), c(site)] * diameter_growth_coefficient[which(height_growth_coefficient[,1] == year), 2], length.out = n.days)
      for (yy in 1 : n.days) {
        S.GPP[yy] <- if (yy == 1) 0 else S.GPP[yy-1] + dS.GPP[yy-1]
        dS.GPP[yy] <- (PF[yy] - S.GPP[yy]) / parameters[c("tau.GPP"), c(site)]
        S.GPP_ref[yy] <- if (yy == 1) 0 else S.GPP_ref[yy-1] + dS.GPP_ref[yy-1]
        dS.GPP_ref[yy] <- (GPP_ref[yy] - S.GPP_ref[yy]) / parameters[c("tau.GPP"), c(site)]
      }
      # Daily LD depends on the GPP of five previous days:
      for (yy in (s.D0 + 1) : n.days) LD[yy] <- parameters[c("LD0"), c(site)] * diameter_growth_coefficient[which(height_growth_coefficient[,1] == year), 2] * S.GPP[yy] / S.GPP_ref[yy]
    }

    # The number of forming cell rows in the tree
    n.rows <- if (xylogenesis == FALSE) {
      # In the orginal hyde version cell.l was one value,
      # in the new parameters I have seperated it into lw and ew so the value of this formula is still the same
      ratios[1,c(site)] * parameters[c("h0"), c(site)] / parameters[c("cell.l.ew"),c(site)] * pi * parameters[c("D0"), c(site)] / parameters[c("cell.d.ew"), c(site)]
    } else {
      # The number of forming cell rows in the tree
      n.rows <- parameters[c("form_factor"),c(site)] * parameters[c("h0"), c(site)] / parameters[c("cell.l.ew"),c(site)] * pi * parameters[c("D0"), c(site)] / parameters[c("cell.d.ew"), c(site)]
    }

    GD <- g.sD.T * fD * LD						# Daily potential number of new cells per day (in one radial cell row), used to be called division
    GD[is.na(GD)] <- 0
    tot.cells.pot <- cumsum(GD)

    storage_reduction = rep(1, n.days)

    if (xylogenesis == TRUE) {
      GD <- g.sD.T * fD * LD	* storage_reduction					# Daily potential number of new cells per day (in one radial cell row), used to be called division
      GD[is.na(GD)] <- 0
      tot.cells <- cumsum(GD)


      if (environment_effect_xylogenesis == TRUE) {        # Constant duration of enlargement and wall formation - end size depends on water availability, final wall thickness depends on temperature
        temperature_effect_on_wall_formation = (g + 1)/1.5
        cell_volume = cell_wall_volume = carbon_to_enlargement = carbon_release_enlargement = carbon_to_wall_formation = phase = matrix(NA, ncol = floor(max(tot.cells)), nrow = n.days)
        ew_lw = cell.wall.density = NULL    # ew = 1, lw = 2
        for (i in 1:floor(max(tot.cells))) {
          division_day = min(which(tot.cells >= i))
          cell_volume[division_day,i] = cell_wall_volume[division_day,i] = 0
          ew_lw[i] = ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), 1, 2)
          tau.E <- ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), round(parameters[c("tau.Ee"), c(site)], 0), round(parameters[c("tau.El"), c(site)], 0))
          tau.W <- ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), round(parameters[c("tau.We"), c(site)], 0), round(parameters[c("tau.Wl"), c(site)], 0))
          cell.wall.density[i] <- ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), cell.wall.density.ew, cell.wall.density.lw)
          cell.volume.growth.per.day = ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), cell.volume.growth.per.day.ew, cell.volume.growth.per.day.lw)
          cell.wall.volume.growth.per.day = ifelse(sD[division_day] < (parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]), cell.wall.volume.growth.per.day.ew, cell.wall.volume.growth.per.day.lw)
          for (d in 1:tau.E) {    # time of enlargement
            cell_volume[(division_day + d), i] = cell_volume[(division_day + d - 1), i] + cell.volume.growth.per.day * soil_moisture_effect[division_day + d] * storage_reduction[division_day + d]
            cell_wall_volume[(division_day + d), i] = 0
            carbon_to_enlargement[(division_day + d), i] = cell.volume.growth.per.day * common[[c("osmotic.sugar.conc")]] * M.suc / (common[[c("gas.const")]] * (Temp[division_day + d] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc  * storage_reduction[division_day + d]
          }
          for (d in (division_day + tau.E + 1) : n.days) {    # time after enlargement
            cell_volume[d,i] = cell_volume[(division_day + tau.E),i]
          }
          for (d in 1:tau.W) {    # time of wall formation
            cell_wall_volume[(division_day + tau.E + d), i] = cell_wall_volume[(division_day + tau.E + d - 1),i] + cell.wall.volume.growth.per.day * temperature_effect_on_wall_formation[division_day + tau.E + d] * storage_reduction[division_day + tau.E + d]
            carbon_to_wall_formation[(division_day + tau.E + d), i] = cell.wall.volume.growth.per.day * temperature_effect_on_wall_formation[division_day + tau.E + d] * storage_reduction[division_day + tau.E + d] * cell.wall.density[i]
          }
          for (d in (division_day + tau.E + tau.W + 1) : n.days) {    # time after wall formation
            cell_wall_volume[d, i] = cell_wall_volume[(division_day + tau.E + tau.W),i]
          }
          carbon_release_enlargement[division_day + tau.E + tau.W, i] = sum(carbon_to_enlargement[,i], na.rm = TRUE)
          phase[1:division_day, i] = 0
          phase[((division_day + 1) : (division_day + tau.E)), i] = 1
          phase[((division_day + tau.E + 1) : (division_day + tau.E + tau.W)), i] = 2
          phase[((division_day + tau.E + tau.W + 1) : n.days), i] = 3
        }
        cell_volume[is.na(cell_volume)] = cell_wall_volume[is.na(cell_wall_volume)] = carbon_to_enlargement[is.na(carbon_to_enlargement)] =
          carbon_release_enlargement[is.na(carbon_release_enlargement)] = carbon_to_wall_formation[is.na(carbon_to_wall_formation)] = 0
        cell_d = sqrt(cell_volume/cell.l)     # matrix ndays * ncells
        cell_d_final = cell_d[365,]           # for each cell (tot.cells.pot values)
        ring_width = rowSums(cell_d, na.rm = TRUE) * 1000 # mm, 365 values
        ew_width = rowSums(cell_d[, which(ew_lw == 1)], na.rm = TRUE) * 1000 # mm, 365 values
        lw_width = rowSums(cell_d[, which(ew_lw == 2)], na.rm = TRUE) * 1000 # mm, 365 values
        ew_cells = sum(ew_lw==1)
        lw_cells = sum(ew_lw==2)
        cell_wall_thickness = cell_wall_volume[365,]/(2 * cell_d[365,]*(2 * cell.l + cell_d[365,])) # for each cell (tot.cells.pot values)
        cell_density = cell_wall_volume[365,] * cell.wall.density / cell_volume[365,]        # for each cell (tot.cells.pot values)
        ring_density = sum(cell_wall_volume[365,] * cell.wall.density) / sum(cell_volume[365,])      # single value

        carbon_to_wall_cellrow = rowSums(carbon_to_wall_formation) # kg C day-1
        carbon_to_enlargement_cellrow = rowSums(carbon_to_enlargement)  # kg C day-1
        enlargement_release_cellrow = rowSums(carbon_release_enlargement)  # kg C day-1

        n.E = rowSums(phase == 1)
        n.W = rowSums(phase == 2)
        n.M = rowSums(phase == 3)

        # Carbon to enlargement per day (kg C /day) (in the whole tree)
        en.growth <- n.rows * carbon_to_enlargement_cellrow
        en.release <- n.rows * enlargement_release_cellrow

        # The use of carbon to wall growth kg C per day
        wall.growth <- n.rows * carbon_to_wall_cellrow

        # Total use of carbon to wall growth so far kg C
        wall.tot <- cumsum(wall.growth)

      }
      if (environment_effect_xylogenesis == FALSE) {
        ## "Uggla" divides cells to early and late wood
        tau.E <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.Ee"), c(site)]) + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.El"), c(site)])       # Duration of cell enlargement of earlywood cells and late wood cells
        tau.W <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.We"), c(site)]) + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.Wl"), c(site)])       # Duration of cell wall formation of earlywood cells and late wood cells

        # Number of cells in enlargement, wall formation and mature phases on each day
        n.E <- n.W <- n.M <- NULL
        for (i in 1 : n.days) {
          n.E[i] <- if (i < 2) 0 else (n.E[i-1] + GD[i] - n.E[i-1] / tau.E[i])
          n.W[i] <- if (i < 2) 0 else (n.W[i-1] + n.E[i-1] / tau.E[i] - n.W[i-1] / tau.W[i-1])
          n.M[i] <- if (i < 2) 0 else (n.M[i-1] + n.W[i-1] / tau.W[i-1])
        }
        # n.E<-round(n.E, digits=0)
        # n.W<-round(n.W, digits=0)
        # n.M<-round(cumsum(GD) - n.W - n.E, digits=0)

        # carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
        CE.ew <- CE.lw <- NULL
        for (i in 1 : n.days){
          CE.ew[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.ew"), c(site)] / 2)^2 * (cell.l) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
          CE.lw[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.ew"), c(site)] / 2)^2 * (cell.l) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
        }
        CE.ew <- unlist(CE.ew)
        CE.lw <- unlist(CE.lw)

        carbon.enlargement <- carbon.wall <- en.growth <- en.release <- wall.growth <- d.growth <- NULL
        # The carbon needed to enlarge the cells  (kg C /day) (in one radial cell row)
        carbon.enlargement <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * CE.ew * n.E + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * CE.lw * n.E
        # Carbon to enlargement per day (kg C /day) (in the whole tree)
        en.growth <- n.rows * carbon.enlargement

        for(i in 1:round(parameters[c("tau.Ee"), c(site)])) en.release[i] <- 0	# the carbon used in enlargement is released after some days.
        for(i in (round(parameters[c("tau.Ee"), c(site)]) + 1) : n.days) en.release[i] <- en.growth[i-round(parameters[c("tau.Ee"), c(site)])]

        # Carbon to wall formation
        # Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!
        CW <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * Carbon.daily.rate.ew + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * Carbon.daily.rate.lw
        CW <- rep(2.9 * 10^-11, length.out = n.days)

        # The use of carbon to wall growth kg C per day
        wall.growth <- n.rows * CW * n.W

        en.growth[is.na(en.growth)] <- 0
        en.release[is.na(en.release)] <- 0
        wall.growth[is.na(wall.growth)] <- 0

        # Total use of carbon to wall growth so far kg C
        wall.tot <- cumsum(wall.growth)

        ew_cells <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * tot.cells
        ew_cells[sD > parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]] = max(ew_cells)
        lw_cells <- tot.cells - ew_cells

        # Calculation of the daily annual ring width, starting from 16.12.2013.
        # Note! This does not include the size of enlarging cells, only wall forming and mature
        diameter_cells = n.M + n.W
        diameter_ew_cells = ifelse(diameter_cells <= ew_cells, diameter_cells, max(ew_cells))
        diameter_lw_cells = diameter_cells - diameter_ew_cells
        ew_width <- diameter_ew_cells * parameters[c("cell.d.ew"), c(site)] * 1000 #+ ifelse(sD < sDc^2 / Uggla, n.E * parameters[c("cell.d.ew"), c(site)] / 2)
        lw_width <- diameter_lw_cells * parameters[c("cell.d.ew"), c(site)] * 1000 #+ ifelse(sD >= sDc^2 / Uggla, n.E * parameters[c("cell.d.ew"), c(site)] / 2)
        ring_width = ew_width + lw_width
        cell_wall_thickness = c(rep(parameters[c("wall.thickness.ew"),c(site)], round(ew_cells[365],0)), rep(parameters[c("wall.thickness.lw"),c(site)], round(lw_cells[365], 0)))
        cell_d_final = c(rep(parameters[c("cell.d.ew"), c(site)], round(ew_cells[365],0)), rep(parameters[c("cell.d.ew"), c(site)], round(lw_cells[365], 0)))
        cell_density = c(rep(cell.density.ew, round(ew_cells[365],0)), rep(cell.density.lw, round(lw_cells[365], 0)))
        ring_density = cumsum(CW * n.W)[365] / (ew_cells[365] * parameters[c("cell.d.ew"), c(site)]^2 * parameters[c("cell.l.ew"),c(site)] + lw_cells[365] * parameters[c("cell.d.ew"), c(site)]^2 * parameters[c("cell.l.lw"),c(site)])

        print(ring_density)
        print(ew.cells_tot)
        print(lw.cells_tot)

        ew_cells = max(ew_cells)
        lw_cells = max (lw_cells)
      }

      list_xylogenesis = list("GD" = GD, "tot.cells" = tot.cells, "n.E" = n.E, "n.W" = n.W, "n.M" = n.M, "ew_cells" = ew_cells, "lw_cells" = lw_cells,
                              "ring_width" = ring_width, "ew_width" = ew_width, "lw_width" = lw_width, "en.growth" = en.growth,
                              "en.release" = en.release, "wall.growth" = wall.growth, "wall.tot" = wall.tot, "cell_d_final" = cell_d_final,
                              "cell_wall_thickness" = cell_wall_thickness, "cell_density" = cell_density, "ring_density" = ring_density)


      en.pot.growth = list_xylogenesis$en.growth
      en.pot.release = list_xylogenesis$en.release
      wall.pot.growth = list_xylogenesis$wall.growth
    } else {
      ## "Uggla" divides cells to early and late wood
      tau.W <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.We"), c(site)]) + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.Wl"), c(site)])
      tau.E <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.Ee"), c(site)]) + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * (parameters[c("tau.El"), c(site)])

      # Number of cells in enlargement, wall formation and mature phases on each day
      n.E.pot <- n.W.pot <- n.M.pot <- NULL
      for (i in 1 : n.days) {
        n.E.pot[i] <- if (i < 2) 0 else (n.E.pot[i-1] + GD[i] - n.E.pot[i-1] / tau.E[i])
        n.W.pot[i] <- if (i < 2) 0 else (n.W.pot[i-1] + n.E.pot[i-1] / tau.E[i] - n.W.pot[i-1] / tau.W[i-1])
        n.M.pot[i] <- if (i < 2) 0 else (n.M.pot[i-1] + n.W.pot[i-1] / tau.W[i-1])
      }

      # carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
      CE.ew <- CE.lw <- NULL
      for (i in 1 : n.days){
        CE.ew[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.ew"), c(site)] / 2)^2 * (parameters[c("cell.l.ew"), c(site)]) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
        CE.lw[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.lw"), c(site)] / 2)^2 * (parameters[c("cell.l.lw"), c(site)]) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
      }
      CE.ew <- unlist(CE.ew)
      CE.lw <- unlist(CE.lw)

      carbon.enlargement.pot <- carbon.wall.pot <- en.pot.growth <- en.pot.release <- wall.pot.growth <- d.pot.growth <- NULL
      # The carbon needed to enlarge the cells  (kg C /day) (in one radial cell row)
      carbon.enlargement.pot <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * CE.ew * n.E.pot + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * CE.lw * n.E.pot
      # Carbon to enlargement and wall forming per day (kg C /day) (in the whole tree)
      en.pot.growth <- n.rows * carbon.enlargement.pot

      for(i in 1:round(parameters[c("tau.Ee"), c(site)])) en.pot.release[i] <- 0	# the carbon used in enlargement is released after some days.
      for(i in (round(parameters[c("tau.Ee"), c(site)]) + 1) : n.days) en.pot.release[i] <- en.pot.growth[i-round(parameters[c("tau.Ee"), c(site)])]

      # Carbon to wall formation
      # Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!
      CW <- (sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * Carbon.daily.rate.ew + (sD >= parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * Carbon.daily.rate.lw

      # The use of carbon to wall growth kg C per day
      wall.pot.growth <- n.rows * CW * n.W.pot

      en.pot.growth[is.na(en.pot.growth)] <- 0
      en.pot.release[is.na(en.pot.release)] <- 0
      wall.pot.growth[is.na(wall.pot.growth)] <- 0

      # Total use of carbon to wall growth so far kg C
      wall.pot <- cumsum(wall.pot.growth)

      # Calculation of the daily annual ring width, starting from 16.12.2013.
      # Note! This does not include the size of enlarging cells, only wall forming and mature
      cells_pot <- ew.cells_pot <- lw.cells_pot <- ew.width_pot <- lw.width_pot <- pot.mm <- NULL
      cells_pot <- (n.W.pot + n.M.pot)
      ew.cells_pot <- ((sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * cells_pot)
      for (i in 1 : n.days) if (sD[i] > parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) lw.cells_pot[i] <- cells_pot[i] - ew.cells_pot[i] - max(ew.cells_pot) else lw.cells_pot[i] <- 0
      ew.width_pot <- ew.cells_pot * parameters[c("cell.d.ew"), c(site)] * 1000
      lw.width_pot <- lw.cells_pot * parameters[c("cell.d.lw"), c(site)] * 1000
      for (i in 1:n.days) if (sD[i] < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) pot.mm[i] <- ew.width_pot[i] else pot.mm[i] <- lw.width_pot[i] + max(ew.width_pot)
    }

    ## New bud growth

    sB <- fB <- GB <- CB <- bud.pot.growth <- bud.pot <- NULL

    # Factor driving the phase of the annual cycle of buds (no. of days after needle growth onset).
    day.no <- 1 : n.days

    # xylogenesis version below
    sB <- cumsum(day.no>parameters[c("sB0"), c(site)])								# The phase of the annual cycle of bud growth

    fB <- (sB < parameters[c("sBc"), c(site)]) * (sin(2 * pi / parameters[c("sBc"), c(site)] * (sB - parameters[c("sBc"), c(site)] / 4)) + 1) / 2		# A function driven by the phase of the annual cycle (annual pattern of growth) [0,1]

    bud.pot.growth <- g * fB * parameters[c("LB"), c(site)]
    bud.pot <- cumsum(bud.pot.growth)

    ############### Maintenance respiration (without carbon limitation) ##########################
    RmN.a <- RmS.a <- RmR.a <- Rm.a <- RmN <- RmS <- RmR <- NULL

    m.S.tot <- ratios[1,c(site)] * (B0) * (parameters[c("h0"), c(site)]) * parameters[c("density_tree"),c(site)] * parameters[c("carbon_share"),c(site)]		# woody carbon mass
    Ra.share <- -0.0007 * Tsa^2 + 0.0424 * Tsa + 0.0273					# share of autotrophic soil respiration
    Ra.share[Ra.share<0] = 0

    Rm_accl = 1
    if (temp_rise == TRUE & Rm_acclimation == TRUE) {
      Rm_accl = 0.85
    }

    RmS.a <- parameters[c("Rm.S"),c(site)] * (exp(log(parameters[c("Q10.S"),c(site)]) / 10 * (Temp)) - 0.7) * m.S.tot * ratios[c("sapwood.share"),c(site)] * Rm_accl # Maintenance respiration of wood
    RmR.a <- parameters[c("Rm.R"),c(site)] * (exp(log(parameters[c("Q10.R"),c(site)]) / 10 * (Tsa)) - exp(-log(parameters[c("Q10.R"),c(site)]) / 2)) * repol[[c("m.R.tot")]] * Ra.share * Rm_accl	# Maintenance respiration of roots
    RmN.a <- parameters[c("Rm.N"),c(site)] * (exp(log(parameters[c("Q10.N"),c(site)]) / 10 * (Temp)) - 0.7) * repol[[c("m.N.tot")]] * Rm_accl						# Maintenance respiration of needles

    if (mN.varies == TRUE){
      m.N.tot2 = NULL
      for(yy in 1 : n.days) {
        m.N.tot2[yy] = if(yy > 150 & yy < 285) repol[[c("m.N.tot")]] else repol[[c("m.N.tot")]] * 2/3
      }
      RmN.a <- parameters[c("Rm.N"),c(site)] * (exp(log(parameters[c("Q10.N"),c(site)]) / 10 * (Temp)) - 0.7) * m.N.tot2 * Rm_accl
    }

    RmN <- (RmN.a > 0) * RmN.a
    RmN[is.na(RmN)] <- 0
    RmS <- (RmS.a>0) * RmS.a
    RmS[is.na(RmS)] <- 0
    RmR <- (RmR.a>0) * RmR.a
    RmR[is.na(RmR)] <- 0

    Rm.a <- RmN + RmS + RmR

    ############ Potential growth and carbon consumption #############
    pot.growth.sofar = cumsum(root.pot.growth + height.pot.growth + wall.pot.growth + needle.pot.growth + bud.pot.growth)

    pot.consumption = (1 + common[[c("Rg.N")]]) * bud.pot.growth + (1 + common[[c("Rg.R")]]) * root.pot.growth +
      (1 + common[[c("Rg.S")]]) * height.pot.growth + (1 + common[[c("Rg.S")]]) * wall.pot.growth +
      (1 + common[[c("Rg.N")]]) * needle.pot.growth + Rm.a

  if (sperling_model == FALSE) {
    # Initial storage values
    starch <- sugar <- storage <- to_sugar <- to_starch <- storage_term <- storage_term_Rm <- NULL
    if (storage.reset == TRUE) {
      sperling[c("starch0"), c(site)] = sperling[c("starch00"), c(site)]
      sperling[c("sugar0"), c(site)] = sperling[c("sugar00"), c(site)]
    }

    starch[1] <- sperling[c("starch0"), c(site)]
    sugar[1] <- sperling[c("sugar0"), c(site)]
    if (n.year == 1) {
      optimal.level <- W.crit <- starch[1] + sugar[1]
      sugar.level <- sugar[1]
      optimal.level.myco <- optimal.level + 0.3
    }
    storage[1] <- starch[1] + sugar[1]
    to_sugar[1] <- to_starch[1] <- 0
    storage_term[1] <- 1
    storage_term_Rm[1] <- 1

    # Modifying storage parameters if storage is assumed to grow as trees grow
    if (storage.grows == TRUE) {
      sperling[c("Wala"),c(site)] <- parameters[c("D0"), c(site)] / D00 * sperling[c("Wala"),c(site)]
      W.crit <- parameters[c("D0"), c(site)] / D00 * optimal.level
      sugar.level <- parameters[c("D0"), c(site)] / D00 * sugar[1]
    }
    a.k <- 1 / (1 - 1/exp(sperling[c("alfa"),c(site)] * (W.crit - sperling[c("Wala"),c(site)])))

    if(mychorrhiza == FALSE){
      parameters[c("growth.myco"),c(site)] = 0
    }

    root.tot.growth <- height.tot.growth <- needle.tot.growth <- wall.tot.growth <- bud.tot.growth <- GD.tot <- NULL
    root.tot <- height.tot <- needle.tot <- wall.tot <- Daily.H.tot <- Daily.N.tot <- cum.Daily.H.tot <- cum.Daily.N.tot <- wall.tot <- tot.cells.tot <- Rm.tot <- RmR.tot <- NULL

    for(i in 2 : n.days) {
      storage_term[i] <- max(0 , min(1 , a.k * (1 - 1 / exp(sperling[c("alfa"),c(site)] * (storage[i-1] - sperling[c("Wala"),c(site)])))))
      storage_term_Rm[i] <- if (storage[i-1] < 0.1) 0 else 1
      sugar[i] <- sugar[i-1] + P[i] - en.pot.growth[i] + en.pot.release[i] - storage_term_Rm[i] * Rm.a[i]-
        (1 + common[[c("Rg.S")]]) * storage_term[i] * height.pot.growth[i]-
        (1 + common[[c("Rg.S")]]) * storage_term[i] * wall.pot.growth[i]-
        (1 + common[[c("Rg.N")]]) * storage_term[i] * needle.pot.growth[i]-
        (1 + common[[c("Rg.R")]]) * storage_term[i] * root.pot.growth[i]-
        (1 + common[[c("Rg.N")]]) * storage_term[i] * bud.pot.growth[i]-
        (sH[i] > parameters[c("sHc"), c(site)]) * (storage[i-1] > optimal.level.myco) * P[i] * parameters[c("growth.myco"),c(site)]
      to_sugar[i] <- if (sugar[i] < sugar.level) (min(starch[i-1], (sugar.level - sugar[i]) / sperling[c("tau.t"), c(site)])) else 0
      to_starch[i] <- if (sugar[i] > sugar.level) (sugar[i] - sugar.level) / sperling[c("tau.s"), c(site)] else 0
      starch[i] <- starch[i-1] + to_starch[i] - to_sugar[i]
      sugar[i] <- sugar[i] + to_sugar[i] - to_starch[i]
      storage[i] <- starch[i] + sugar[i]
    }

  } else if (sperling_model == TRUE) {

    ########### Sperling Model

    #####
    ## Set up equilibrium points and initial values
    #####

    storage_term <- storage_term_Rm <- storage_term_needles <- storage_term_phloem <- storage_term_roots <- storage_term_xylem.sh <- storage_term_xylem.st <- NULL
    root.tot.growth <- height.tot.growth <- needle.tot.growth <- wall.tot.growth <- bud.tot.growth <- GD.tot <- NULL
    root.tot <- height.tot <- needle.tot <- wall.tot <- Daily.H.tot <- Daily.N.tot <- cum.Daily.H.tot <- cum.Daily.N.tot <- wall.tot <- tot.cells.tot <- Rm.tot <- RmR.tot <- NULL
    Ad.needles <- Ad.phloem <- Ad.roots <- Ad.xylem.sh <- Ad.xylem.st <- As.needles <- As.phloem <- As.roots <- As.xylem.sh <- As.xylem.st <- NULL
    Ks.needles <- Ks.phloem <- Ks.roots <- Ks.xylem.sh <- Ks.xylem.st <- Kd.needles <- Kd.phloem <- Kd.roots <- Kd.xylem.sh <- Kd.xylem.st <-  NULL
    sugar.needles <- sugar.phloem <- sugar.roots <- sugar.xylem.sh <- sugar.xylem.st <- starch.needles <- starch.phloem <- starch.roots <- starch.xylem.sh <- starch.xylem.st <- NULL
    to_sugar.needles <- to_sugar.phloem <- to_sugar.roots <- to_sugar.xylem.sh <- to_sugar.xylem.st <- to.mycorrhiza <-  NULL
    sugar.to_phloem <- sugar.to_roots <- starch.to_phloem <- starch.to_roots <- sugar.to_xylem.sh <- sugar.to_xylem.st <- sugar.to_myco <- NULL
    DF_np <- DF_pr <- DF_rm <- DF_pxsh <- DF_pxst <- NULL

    As0 <- Te0 <- NULL

    starch0 <- sperling[c("starch.needles0"),c(site)] + sperling[c("starch.phloem0"),c(site)] + sperling[c("starch.roots0"),c(site)] + sperling[c("starch.xylem.sh0"),c(site)] + sperling[c("starch.xylem.st00"),c(site)]
    sugar0 <- sperling[c("sugar.needles0"),c(site)] + sperling[c("sugar.phloem0"),c(site)] + sperling[c("sugar.roots0"),c(site)] + sperling[c("sugar.xylem.sh0"),c(site)] + sperling[c("sugar.xylem.st00"),c(site)]
    W.crit.needles <- sperling[c("starch.needles0"),c(site)] + sperling[c("sugar.needles0"),c(site)]
    W.crit.phloem <- sperling[c("starch.phloem0"),c(site)] + sperling[c("sugar.phloem0"),c(site)]
    W.crit.roots <- sperling[c("starch.roots0"),c(site)] + sperling[c("sugar.roots0"),c(site)]
    W.crit.xylem.sh <- sperling[c("starch.xylem.sh0"),c(site)] + sperling[c("sugar.xylem.sh0"),c(site)]
    W.crit.xylem.st <- sperling[c("starch.xylem.st0"),c(site)] + sperling[c("sugar.xylem.st0"),c(site)]
    a.k.needles <- 1 / (1 - 1/exp(sperling[c("alfa.needles"),c(site)] * (W.crit.needles - sperling[c("Wala.needles"),c(site)])))
    a.k.phloem <- 1 / (1 - 1/exp(sperling[c("alfa.phloem"),c(site)] * (W.crit.phloem - sperling[c("Wala.phloem"),c(site)])))
    a.k.roots <- 1 / (1 - 1/exp(sperling[c("alfa.roots"),c(site)]* (W.crit.roots - sperling[c("Wala.roots"),c(site)])))
    a.k.xylem.st <- 1 / (1 - 1/exp(sperling[c("alfa.xylem.st"),c(site)]* (W.crit.xylem.st - sperling[c("Wala.xylem.st"),c(site)])))
    a.k.xylem.sh <- 1 / (1 - 1/exp(sperling[c("alfa.xylem.sh"),c(site)]* (W.crit.xylem.sh - sperling[c("Wala.xylem.sh"),c(site)])))

    if (storage.grows == TRUE) {
      sperling[c("Wala.needles"),c(site)] <- HN0 / D00 * sperling[c("Wala.needles"),c(site)]
      sperling[c("Wala.phloem"),c(site)] <- D0 / D00 * sperling[c("Wala.phloem"),c(site)]
      sperling[c("Wala.roots"),c(site)] <- LR0 / D00 * sperling[c("Wala.roots"),c(site)]
      sperling[c("Wala.xylem.st"),c(site)] <- LR0 / D00 * sperling[c("Wala.xylem.st"),c(site)]
      sperling[c("Wala.xylem.sh"),c(site)] <- LR0 / D00 * sperling[c("Wala.xylem.sh"),c(site)]
      W.crit.needles <- HN0 / D00 * W.crit.needles
      W.crit.phloem <- D0 / D00 * W.crit.phloem
      W.crit.roots <- LR0 / D00 * W.crit.roots
      W.crit.xylem.st <- LR0 / D00 * W.crit.xylem.st
      W.crit.xylem.sh <- LR0 / D00 * W.crit.xylem.sh
    }

    storage_term_needles <- storage_term_phloem <- storage_term_roots <- storage_term_xylem.sh <- storage_term_xylem.st <- rep(1, length = n.days)

    Bd=log(sperling[c("Q10d"),c(site)])/10 		# Compute energies of activation
    Bs=log(sperling[c("Q10s"),c(site)])/10
    Te0=mean(Temp[274:281], na.rm = T)+3 	# Compute initial Te by the mean temperature for the first week of # October plus 3C (for the exponential nature of the curves)
    # The indexes changed here as CASSIA is from Jan not Oct, also changed to represent days not hours
    As0.needles=sperling[c("Ad0.needles"),c(site)]*exp(Te0*(Bd-Bs)) 	# Initial As by Ad0 and initial Te
    As0.phloem=sperling[c("Ad0.phloem"),c(site)]*exp(Te0*(Bd-Bs))
    As0.roots=sperling[c("Ad0.roots"),c(site)]*exp(Te0*(Bd-Bs))
    As0.xylem.sh=sperling[c("Ad0.xylem.sh"),c(site)]*exp(Te0*(Bd-Bs))
    As0.xylem.st=sperling[c("Ad0.xylem.st"),c(site)]*exp(Te0*(Bd-Bs))
    As.needles[1] <- As0.needles
    As.phloem[1] <- As0.phloem
    As.roots[1] <- As0.roots
    As.xylem.sh[1] <- As0.xylem.sh
    As.xylem.st[1] <- As0.xylem.st
    Ad.needles[1] <- sperling[c("Ad0.needles"),c(site)]
    Ad.phloem[1] <- sperling[c("Ad0.phloem"),c(site)]
    Ad.roots[1] <- sperling[c("Ad0.roots"),c(site)]
    Ad.xylem.sh[1] <- sperling[c("Ad0.xylem.sh"),c(site)]
    Ad.xylem.st[1] <- sperling[c("Ad0.xylem.st"),c(site)]
    Kd.needles=sperling[c("Ad0.needles"),c(site)]*exp(Bd*Temp)
    Kd.phloem=sperling[c("Ad0.phloem"),c(site)]*exp(Bd*Temp)
    Kd.roots=sperling[c("Ad0.roots"),c(site)]*exp(Bd*Temp)
    Kd.xylem.sh=sperling[c("Ad0.xylem.sh"),c(site)]*exp(Bd*Temp)
    Kd.xylem.st=sperling[c("Ad0.xylem.st"),c(site)]*exp(Bd*Temp)
    Ks.needles=As0.needles*exp(Bs*Temp)
    Ks.phloem=As0.phloem*exp(Bs*Temp)
    Ks.roots=As0.roots*exp(Bs*Temp)
    Ks.xylem.sh=As0.xylem.sh*exp(Bs*Temp)
    Ks.xylem.st=As0.xylem.st*exp(Bs*Temp)

    # initial values from 17.3
    starch.needles[1] <- sperling[c("starch.needles0"),c(site)] # kg C / raw material
    starch.phloem[1] <- sperling[c("starch.phloem0"),c(site)]
    starch.roots[1] <- sperling[c("starch.roots0"),c(site)]
    starch.xylem.sh[1] <- sperling[c("starch.xylem.sh0"),c(site)]
    starch.xylem.st[1] <- sperling[c("starch.xylem.st0"),c(site)]
    sugar.needles[1] <- sperling[c("sugar.needles0"),c(site)]
    sugar.phloem[1] <- sperling[c("sugar.phloem0"),c(site)]
    sugar.roots[1] <- sperling[c("sugar.roots0"),c(site)]
    sugar.xylem.sh[1] <- sperling[c("sugar.xylem.sh0"),c(site)]
    sugar.xylem.st[1] <- sperling[c("sugar.xylem.st0"),c(site)]

    DF_np[1] = (sugar.needles[1]+starch.needles[1] - (1/3)*(sugar.phloem[1]+starch.phloem[1]))
    DF_pr[1] = (sugar.phloem[1]+starch.phloem[1] - 11*(sugar.roots[1]+starch.roots[1]))
    DF_rm[1] = (sugar.roots[1]+starch.roots[1] - sperling[c("myco.thresh"),c(site)])
    DF_pxsh[1] = max(sugar.phloem[1]+starch.phloem[1], 0) - 11*max(sugar.xylem.sh[1]+starch.xylem.sh[1], 0)
    DF_pxst[1] = max(sugar.phloem[1]+starch.phloem[1], 0) - 11*max(sugar.xylem.st[1]+starch.xylem.st[1], 0)

    for (i in 2 : n.days) { # Sperling model has been added here, SC = sugar and ST = starch to match the variable names already in CASSIA
      storage_term_needles[i] <- max(0 , min(1, a.k.needles * (1 - 1 / exp(sperling[c("alfa.needles"),c(site)] * (starch.needles[i-1] + sugar.needles[i-1] - sperling[c("Wala.needles"),c(site)])))))
      storage_term_phloem[i] <- max(0 , min(1, a.k.phloem * (1 - 1 / exp(sperling[c("alfa.phloem"),c(site)] * (starch.phloem[i-1] + sugar.phloem[i-1] - sperling[c("Wala.phloem"),c(site)])))))
      storage_term_roots[i] <- max(0 , min(1, a.k.roots * (1 - 1 / exp(sperling[c("alfa.roots"),c(site)] * (starch.roots[i-1] + sugar.roots[i-1] - sperling[c("Wala.roots"),c(site)])))))
      storage_term_xylem.sh[i] <- max(0 , min(1, a.k.xylem.sh * (1 - 1 / exp(sperling[c("alfa.xylem.sh"),c(site)] * (starch.xylem.sh[i-1] + sugar.xylem.sh[i-1] - sperling[c("Wala.xylem.sh"),c(site)])))))
      storage_term_xylem.st[i] <- max(0 , min(1, a.k.xylem.st * (1 - 1 / exp(sperling[c("alfa.xylem.st"),c(site)] * (starch.xylem.st[i-1] + sugar.xylem.st[i-1] - sperling[c("Wala.xylem.st"),c(site)])))))

      Ks.needles[i]=As.needles[i-1]*exp(Bs*Temp[i]) # Compute activity (K) mg g-1 DW day -1 by the previous frequency (A)
      Ks.phloem[i]=As.phloem[i-1]*exp(Bs*Temp[i]) # Compute activity (K) mg g-1 DW day -1 by the previous frequency (A)
      Ks.roots[i]=As.roots[i-1]*exp(Bs*Temp[i]) # Compute activity (K) mg g-1 DW day -1 by the previous frequency (A)
      Ks.xylem.sh[i]=As.xylem.sh[i-1]*exp(Bs*Temp[i]) # Compute activity (K) mg g-1 DW day -1 by the previous frequency (A)
      Ks.xylem.st[i]=As.xylem.st[i-1]*exp(Bs*Temp[i]) # Compute activity (K) mg g-1 DW day -1 by the previous frequency (A)
      Kd.needles[i]=Ad.needles[i-1]*exp(Bd*Temp[i])
      Kd.phloem[i]=Ad.phloem[i-1]*exp(Bd*Temp[i])
      Kd.roots[i]=Ad.roots[i-1]*exp(Bd*Temp[i])
      Kd.xylem.sh[i]=Ad.xylem.sh[i-1]*exp(Bd*Temp[i])
      Kd.xylem.st[i]=Ad.xylem.st[i-1]*exp(Bd*Temp[i])

      # if there is a surplus goes to next organ down - concentration driven model
      # The differences are normalised by a multiplier which represents the average difference in magnitude between the two stores
      # otherwise all of the sugar would just immediately go to the roots
      # NOTE: the forces and therefore sugar transfered are worked out for all organs based on the amount of sugar there in the beginning
      # this could lead to a slight error, but should be corrected by the starch latter just have to imagine that all of the sugar
      # goes to the allocated organs simultaneously

      DF_np[i] = 3*max(sugar.needles[i-1]+starch.needles[i-1], 0) - max(sugar.phloem[i-1]+starch.phloem[i-1], 0)
      DF_pr[i] = max(sugar.phloem[i-1]+starch.phloem[i-1], 0) - 11*max(sugar.roots[i-1]+starch.roots[i-1], 0)
      DF_pxsh[i] = max(sugar.phloem[i-1]+starch.phloem[i-1], 0) - 8*max(sugar.xylem.sh[i-1]+starch.xylem.sh[i-1], 0)
      DF_pxst[i] = max(sugar.phloem[i-1]+starch.phloem[i-1], 0) - 2*max(sugar.xylem.st[i-1]+starch.xylem.st[i-1], 0)
      # This one works from a threshold as mycrorhiza is not considered as an organ in the model
      DF_rm[i] = min(max(sugar.roots[i-1]+starch.roots[i-1] - sperling[c("myco.thresh"),c(site)],0), sugar.roots[i-1])

      # Rm.a matainence resperation seperated into organs

      sugar.needles[i] <- sugar.needles[i-1] + P[i] -
        RmN[i] * storage_term_needles[i] - # maintenance respiration
        (1 + common[[c("Rg.N")]]) * storage_term_needles[i] * (needle.pot.growth[i] + bud.pot.growth[i]) - # growth
        en.pot.growth[i] + en.pot.release[i] - # growth use and release and to the rest of the organs
        sperling[c("k_np"),c(site)] * DF_np[i] + # transfer between organs
        (Kd.needles[i] - Ks.needles[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * needle_mass[n.year] # links to the needle growth process

      ### coefficient is from mass ratio in starch and sugar 2015 xls
      sugar.phloem[i] <- sugar.phloem[i-1] -
        0.082179938 * RmS[i] * storage_term_phloem[i] - # maintenance respiration
        0.082179938 * (1 + common[[c("Rg.S")]]) * storage_term_phloem[i] * (wall.pot.growth[i] + height.pot.growth[i]) + # growth
        sperling[c("k_np"),c(site)] * DF_np[i] - # transfer between organs
        sperling[c("k_pr"),c(site)] * DF_pr[i] - # transfer between organs
        sperling[c("k_pxsh"),c(site)] * DF_pxsh[i] -
        sperling[c("k_pxst"),c(site)] * DF_pxst[i] +
        (Kd.phloem[i] - Ks.phloem[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 7.4

      sugar.roots[i] <- sugar.roots[i-1] +
        sperling[c("k_pr"),c(site)] * DF_pr[i] - # transfer between organs
        DF_rm[i] + # transfer between organs, no multiplier as this is for mycorhiza and the model just takes the extra sugar
        (Kd.roots[i] - Ks.roots[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8 -
        (1 + common[[c("Rg.R")]]) * storage_term_roots[i] * root.pot.growth[i] - # growth
        RmR[i] * storage_term_roots[i] # maintenance respiration

      ### coefficient is from mass ratio in starch and sugar 2015 xls
      sugar.xylem.sh[i] <- sugar.xylem.sh[i-1] -
        0.096020683 * RmS[i] * storage_term_xylem.sh[i] - # maintenance respiration
        0.096020683 * (1 + common[[c("Rg.S")]]) * storage_term_xylem.sh[i] * (wall.pot.growth[i] + height.pot.growth[i]) + # growth
        sperling[c("k_pxsh"),c(site)] * DF_pxsh[i] +
        (Kd.xylem.sh[i] - Ks.xylem.sh[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8

      ### coefficient is from mass ratio in starch and sugar 2015 xls
      sugar.xylem.st[i] <- sugar.xylem.st[i-1] -
        0.821799379 * RmS[i] * storage_term_xylem.st[i] - # maintenance respiration
        0.821799379 * (1 + common[[c("Rg.S")]]) * storage_term_xylem.st[i] * (wall.pot.growth[i] + height.pot.growth[i]) + # growth
        sperling[c("k_pxst"),c(site)] * DF_pxst[i] +
        (Kd.xylem.st[i] - Ks.xylem.st[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8

      to.mycorrhiza[i] <- DF_rm[i]
      # As the tree now has the bucket model I have changed photosynthesis derived allocation to just excess in roots as it should get here from the tree
      # I am assuming that that this is why photosynthesis was driving it before, it was just a measure of the excess photosynthates
      # (sH[i] > sHc) * ((sugar.roots[i-1] + starch.roots[i-1]) > optimal.level.myco) * P[i] * growth.myco # belowground allocation

      # carbon sugar is used here as the model was in terms of sugar and is being transformed to kg C
      starch.needles[i] <- starch.needles[i-1] + (- Kd.needles[i] + Ks.needles[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * needle_mass[n.year] # Subtract starch degradation and add synthase to ST
      starch.phloem[i] <- starch.phloem[i-1] + (- Kd.phloem[i] + Ks.phloem[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 7.4 # Subtract starch degradation and add synthase to ST
      starch.roots[i] <- starch.roots[i-1] + (- Kd.roots[i] + Ks.roots[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8 # Subtract starch degradation and add synthase to ST
      # TOOD: are the densities right here?
      starch.xylem.sh[i] <- starch.xylem.sh[i-1] + (- Kd.xylem.sh[i] + Ks.xylem.sh[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8 # Subtract starch degradation and add synthase to ST
      starch.xylem.st[i] <- starch.xylem.st[i-1] + (- Kd.xylem.st[i] + Ks.xylem.st[i]) * sperling[c("carbon.sugar"),c(site)] * 0.001 * 2.8 # Subtract starch degradation and add synthase to ST
      # Sugar already done above

      # If sugar is below a certain value starch is released so the sugar doesn't go negative before starch
      # This is a proxy for a starch metabolism system, which seems to be present under stress in literature
      # but I can't find a mechanism for scots pine
      # values are below the lowest recorded value
      to_sugar.needles[i] <- if (sugar.needles[i] < 0.05) (min(starch.needles[i], max((0.05 - sugar.needles[i]) / sperling[c("tau.t"),c(site)], 0))) else 0
      to_sugar.phloem[i] <- if (sugar.phloem[i] < 0.12) (min(starch.phloem[i], max((0.12 - sugar.phloem[i]) / sperling[c("tau.t"),c(site)], 0))) else 0
      to_sugar.roots[i] <- if (sugar.roots[i] < 0.005) (min(starch.roots[i], max((0.005 - sugar.roots[i]) / sperling[c("tau.t"),c(site)], 0))) else 0
      to_sugar.xylem.sh[i] <- if (sugar.xylem.sh[i] < 0.009) (min(starch.xylem.sh[i], max((0.009 - sugar.xylem.sh[i]) / sperling[c("tau.t"),c(site)], 0))) else 0
      to_sugar.xylem.st[i] <- if (sugar.xylem.st[i] < 0.001) (min(starch.xylem.st[i], max((0.001 - sugar.xylem.st[i]) / sperling[c("tau.t"),c(site)], 0))) else 0

      # storage update, both bucket and emergency as the level shouldn't be in both
      starch.needles[i] <- starch.needles[i] - to_sugar.needles[i]
      starch.phloem[i] <- starch.phloem[i] - to_sugar.phloem[i]
      starch.roots[i] <- starch.roots[i] - to_sugar.roots[i]
      starch.xylem.sh[i] <- starch.xylem.sh[i] - to_sugar.xylem.sh[i]
      starch.xylem.st[i] <- starch.xylem.st[i] - to_sugar.xylem.st[i]
      sugar.needles[i] <- sugar.needles[i]  + to_sugar.needles[i]
      sugar.phloem[i] <- sugar.phloem[i]  + to_sugar.phloem[i]
      sugar.roots[i] <- sugar.roots[i]  + to_sugar.roots[i]
      sugar.xylem.sh[i] <- sugar.xylem.sh[i]  + to_sugar.xylem.sh[i]
      sugar.xylem.st[i] <- sugar.xylem.st[i]  + to_sugar.xylem.st[i]

      As.needles[i]=(1-sperling[c("lamda.needles"),c(site)])*As.needles[i-1]
      As.phloem[i]=(1-sperling[c("lamda.phloem"),c(site)])*As.phloem[i-1]
      As.roots[i]=(1-sperling[c("lamda.roots"),c(site)])*As.roots[i-1]
      As.xylem.sh[i]=(1-sperling[c("lamda.xylem.sh"),c(site)])*As.xylem.sh[i-1]
      As.xylem.st[i]=(1-sperling[c("lamda.xylem.st"),c(site)])*As.xylem.st[i-1]
      Ad.needles[i]=(1-sperling[c("lamda.needles"),c(site)])*Ad.needles[i-1]
      Ad.phloem[i]=(1-sperling[c("lamda.phloem"),c(site)])*Ad.phloem[i-1]
      Ad.roots[i]=(1-sperling[c("lamda.roots"),c(site)])*Ad.roots[i-1]
      Ad.xylem.sh[i]=(1-sperling[c("lamda.xylem.sh"),c(site)])*Ad.xylem.sh[i-1]
      Ad.xylem.st[i]=(1-sperling[c("lamda.xylem.st"),c(site)])*Ad.xylem.st[i-1]
      # Induce starch synthase if SC is high or degradation if it is low
      # These numbers are from september 2018
      # xylem not changed as no data to support it
      if  (sugar.needles[i]>0.1351581) {As.needles[i]=As.needles[i]+sperling[c("delta.needles"),c(site)]}
      else if (sugar.needles[i]<0.1351581 && starch.needles[i]> 0) {Ad.needles[i]=Ad.needles[i]+sperling[c("delta.needles"),c(site)]}
      if  (sugar.phloem[i]>0.5269149) {As.phloem[i]=As.phloem[i]+sperling[c("delta.phloem"),c(site)]}
      else if (sugar.phloem[i]<0.5269149 && starch.phloem[i]> 0) {Ad.phloem[i]=Ad.phloem[i]+sperling[c("delta.phloem"),c(site)]}
      if  (sugar.roots[i]>0.04911007) {As.roots[i]=As.roots[i]+sperling[c("delta.roots"),c(site)]}
      else if (sugar.roots[i]<0.04911007 && starch.roots[i]> 0) {Ad.roots[i]=Ad.roots[i]+sperling[c("delta.roots"),c(site)]}
      if  (sugar.xylem.sh[i]>0.0199) {As.xylem.sh[i]=As.xylem.sh[i]+sperling[c("delta.xylem.sh"),c(site)]}
      else if (sugar.xylem.sh[i]<0.0199 && starch.xylem.sh[i]> 0) {Ad.xylem.sh[i]=Ad.xylem.sh[i]+sperling[c("delta.xylem.sh"),c(site)]}
      if  (sugar.xylem.st[i]>0.0199) {As.xylem.st[i]=As.xylem.st[i]+sperling[c("delta.xylem.st"),c(site)]}
      else if (sugar.xylem.st[i]<0.0199 && starch.xylem.st[i]> 0) {Ad.xylem.st[i]=Ad.xylem.st[i]+sperling[c("delta.xylem.st"),c(site)]}

    }

    if (count%%2 != 0) {
      if (length(which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])) == 0) {
        warning(paste("Never cold enough for the model to trigger bud burst, sugar never lower than", sperling[c("SCb"),c(site)], "kg C"))
      } else {
        sB0 <- which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])[1]
      }
    }

  }

    ########### Total growth and carbon consumption
    #  Occurred growth kg C day-1 (potential growth * storage effect)
    # TODO: should this be in vector form or for loop?
    root.tot.growth <- if (sperling_model == FALSE) storage_term * root.pot.growth else {storage_term_roots * root.pot.growth}
    height.tot.growth <- if (sperling_model == FALSE) storage_term * height.pot.growth else height.pot.growth * (0.082179938 * storage_term_phloem + 0.821799379 * storage_term_xylem.st + 0.096020683 * storage_term_xylem.sh)
    needle.tot.growth <- if (sperling_model == FALSE) storage_term * needle.pot.growth else storage_term_needles * needle.pot.growth
    wall.tot.growth <- if (sperling_model == FALSE & xylogenesis == FALSE) {storage_term * wall.pot.growth}
    # TODO: This should be done with xylogenesis using starch storage, make xylogenesis a seperate function
    else if (sperling_model == FALSE & xylogenesis == TRUE) {list_xylogenesis$wall.growth}
    else {wall.pot.growth * (0.082179938 * storage_term_phloem + 0.821799379 * storage_term_xylem.st + 0.096020683 * storage_term_xylem.sh)}
    bud.tot.growth <- if (sperling_model == FALSE) storage_term * bud.pot.growth else storage_term_needles * bud.pot.growth
    # These shouldn't be for intervdal
    average_storage <- (storage_term_needles+storage_term_phloem+storage_term_roots + storage_term_xylem.sh + storage_term_xylem.sh)/5
    GD.tot <-  if (sperling_model == FALSE & xylogenesis == FALSE) {storage_term * GD}
    # TODO: This should be done with xylogenesis using starch storage, make xylogenesis a seperate function
    else if (sperling_model == FALSE & xylogenesis == TRUE) {list_xylogenesis$GD}
    else if (sperling_model == TRUE & xylogenesis == FALSE) {GD * average_storage} # Average from all of the organs
    Rm.tot <- if (sperling_model == FALSE) storage_term_Rm * Rm.a else storage_term_roots * RmR + storage_term_needles * RmN + 0.082179938 * storage_term_phloem * RmS + 0.821799379 * storage_term_xylem.st * RmS + 0.096020683 * storage_term_xylem.sh * RmS
    RmR.tot <- if (sperling_model == FALSE) storage_term_Rm * Rm.a else storage_term_roots * RmR

    # TODO: add back in when the model_xylogenesis function works
    # if (xylogenesis == TRUE) {tot_xylogenesis = model_xylogenesis(storage_term)}

    # Total occurred growth so far kg C

    root.tot <- cumsum(root.tot.growth)
    height.tot <- cumsum(height.tot.growth)
    needle.tot <- cumsum(needle.tot.growth)
    wall.tot <- cumsum(wall.tot.growth)
    bud.tot <- cumsum(bud.tot.growth)

    # Dimensional growth
    Daily.H.tot <- height.tot.growth / (B0 * CH / 1000 * ratios[1,c(site)])		# mm day-1
    cum.Daily.H.tot <- parameters[c("HH0"), c(site)] + cumsum(Daily.H.tot)					# mm
    Daily.N.tot <- needle.tot.growth / (repol[[c("m.N")]] * max(cum.Daily.H.tot) / parameters[c("h_increment"), c(site)])		# mm day-1
    cum.Daily.N.tot <- parameters[c("HN0"), c(site)] + cumsum(Daily.N.tot)							# mm

    if (xylogenesis == FALSE) {

      # Secondary cells and ring width
      GD.tot[is.na(GD.tot)] <- 0
      tot.cells.tot <- cumsum(GD.tot)							# cells in new ring

      n.E.tot <- n.W.tot <- n.M.tot <- NULL 		# Number of cells in phases on each day
      for(i in 1 : n.days) n.E.tot[i] <- if(i < 2) 0 else (n.E.tot[i-1] + GD.tot[i] - n.E.tot[i-1] / tau.E[i])
      for(i in 1 : n.days) n.W.tot[i] <- if(i < 2) 0 else (n.W.tot[i-1] + n.E.tot[i-1] / tau.E[i] - n.W.tot[i-1] / tau.W[i-1])
      for(i in 1 : n.days) n.M.tot[i] <- if(i < 2) 0 else (n.M.tot[i-1] + n.W.tot[i-1] / tau.W[i-1])

      CW_new <- daily_carbon_to_walls <- NULL
      for(i in 1 : n.days) if (n.W.tot[i] < common[[c("ypsilon")]]) daily_carbon_to_walls[i] <- 0 else daily_carbon_to_walls[i] <- wall.tot.growth[i] / (n.W.tot[i] * n.rows)
      for(i in 1 : n.days) CW_new[i] <- if (daily_carbon_to_walls[i] < CW[i]) daily_carbon_to_walls[i] else CW_new[i] <- CW[i]

      cells_tot <- ew.cells_tot <- lw.cells_tot <- ew.width_tot <- lw.width_tot <- tot.mm <- NULL
      cells_tot <- (n.W.tot + n.M.tot)
      ew.cells_tot <- ((sD < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) * cells_tot)
      for (i in 1 : n.days) if (sD[i] > parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) lw.cells_tot[i] <- cells_tot[i] - ew.cells_tot[i] - max(ew.cells_tot) else lw.cells_tot[i] <- 0
      ew.width_tot <- ew.cells_tot * parameters[c("cell.d.ew"), c(site)] * 1000
      lw.width_tot <- lw.cells_tot * parameters[c("cell.d.lw"), c(site)] * 1000
      for (i in 1 : n.days) if (sD[i] < parameters[c("sDc"), c(site)]^2 / parameters[c("Uggla"), c(site)]) tot.mm[i] <- ew.width_tot[i] else tot.mm[i] <- lw.width_tot[i] + max(ew.width_tot)

      # wood_density_model<-max(ew.width_tot)/tot.mm[n.days]*cell.density.ew+lw.width_tot[n.days]/tot.mm[n.days]*cell.density.lw	# not used
      wood_density.CW <- wall.tot[365] / (ratios[1,c(site)] * parameters[c("h0"), c(site)] * pi * ((parameters[c("D0"), c(site)] / 2 + tot.mm[n.days] / 1000)^2 - (parameters[c("D0"), c(site)] / 2)^2))
      #print(paste(year, "Proportion of earlywood:", round(max(ew.width_tot) / max(tot.mm), digits=2)," wood density:", round(wood_density.CW)," ring width:", round(tot.mm[n.days], digits=2)))

    } else if (xylogenesis == TRUE) {

      tot.cells.tot <- cumsum(GD.tot)							# cells in new ring
      n.E.tot = list_xylogenesis$n.E
      n.W.tot = list_xylogenesis$n.W
      n.M.tot = list_xylogenesis$n.M
      ew.cells_tot = list_xylogenesis$ew_cells
      lw.cells_tot = list_xylogenesis$lw_cells
      ring_width = list_xylogenesis$ring_width
      ew_width = list_xylogenesis$ew_width
      lw_width = list_xylogenesis$lw_width
      ring_density = list_xylogenesis$ring_density

    }

    # Carbon to mychorrhiza
    if (sperling_model == FALSE) {
      to.mycorrhiza = (Tsb > 10) * (storage > optimal.level.myco) * P * parameters[c("growth.myco"),c(site)]
      mycorrhiza.tot = cumsum(to.mycorrhiza)
    } else {
      mycorrhiza.tot = cumsum(to.mycorrhiza)
    }

    # Carbon used for respiration
    Rg.tot <- tot.Rm <- tot.Rg <- NULL
    # TODO: Already asked Pauliina: should the bud growth respiration be assumed to be the same as the needle respiration
    Rg.tot <- common[[c("Rg.N")]] * needle.tot.growth + common[[c("Rg.N")]] * bud.tot.growth + common[[c("Rg.R")]] * root.tot.growth + common[[c("Rg.S")]] * height.tot.growth + common[[c("Rg.S")]] * wall.tot.growth
    Rg.root <- common[[c("Rg.R")]] * root.tot.growth
    # NOTE! These are the same for sperling and non-sperling as this is the sum of possible Rm and Rg not what actually happened
    tot.Rm <- cumsum(Rm.tot)
    tot.Rg <- cumsum(Rg.tot)

    if (sperling_model == FALSE) {
      parameters[c("starch0"), c(site)] <- starch[n.days]		# initial value of starch for the following year
      parameters[c("sugar0"), c(site)] <- sugar[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
    } else {
      parameters[c("starch.needles0"), c(site)] <- starch.needles[n.days]		# initial value of starch for the following year
      parameters[c("sugar.needles0"), c(site)] <- sugar.needles[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
      parameters[c("starch.phloem0"), c(site)] <- starch.phloem[n.days]		# initial value of starch for the following year
      parameters[c("sugar.phloem0"), c(site)] <- sugar.phloem[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
      parameters[c("starch.xylem.sh0"), c(site)] <- starch.xylem.sh[n.days]		# initial value of starch for the following year
      parameters[c("sugar.xylem.sh0"), c(site)] <- sugar.xylem.sh[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
      parameters[c("starch.xylem.st0"), c(site)] <- starch.xylem.st[n.days]		# initial value of starch for the following year
      parameters[c("sugar.xylem.st0"), c(site)] <- sugar.xylem.st[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
      parameters[c("starch.roots0"), c(site)] <- starch.roots[n.days]		# initial value of starch for the following year
      parameters[c("sugar.roots0"), c(site)] <- sugar.roots[n.days]		# the initial value for labile sugars for the following year (~optimal.level)
    }

    # New tree dimensions if trees grow
    if (trees_grow == TRUE) {
      parameters[c("D0"), c(site)] <- if (xylogenesis == FALSE) {parameters[c("D0"), c(site)] + tot.mm[n.days] * 2 / 1000}	else {parameters[c("D0"), c(site)] + ring_width[n.days] * 2 / 1000}		# The diameter in the begining of next year
      B0 <- pi / 4 * parameters[c("D0"), c(site)]^2							# Basal area
      parameters[c("h0"), c(site)] <- parameters[c("h0"), c(site)] + (cum.Daily.H.tot[n.days] / 1000)	# The top shoot probably grows faster than the average of other shoots
    }

    #### Total daily growth, total growth so far, total carbon consumption and total respiration (kg C day-1 or kg C)
    daily.tot.growth <- daily.tot <- daily.consumption <- tot.consumption <- NULL
    daily.tot.growth <- height.tot.growth + wall.tot.growth + needle.tot.growth + root.tot.growth + bud.tot.growth  	# Total daily growth
    if (sperling_model == F) {
      daily.consumption <- storage_term_Rm * Rm.a +
        (1 + common[[c("Rg.S")]]) * storage_term * height.pot.growth +
        (1 + common[[c("Rg.S")]]) * storage_term * wall.pot.growth +
        (1 + common[[c("Rg.N")]]) * storage_term * needle.pot.growth +
        (1 + common[[c("Rg.R")]]) * storage_term * root.pot.growth +
        (1 + common[[c("Rg.N")]]) * storage_term * bud.pot.growth +
        to.mycorrhiza								# Total daily consumption
    } else {
      daily.consumption <- 0.082179938 * RmS * storage_term_phloem + 0.096020683 * RmS * storage_term_xylem.sh + 0.821799379 * RmS * storage_term_xylem.st +
        RmN * storage_term_needles + RmR * storage_term_roots +
        (1 + common[[c("Rg.S")]]) * height.pot.growth * (0.082179938 * storage_term_phloem + 0.821799379 * storage_term_xylem.st + 0.096020683 * storage_term_xylem.sh) +
        (1 + common[[c("Rg.S")]]) * wall.pot.growth * (0.082179938 * storage_term_phloem + 0.821799379 * storage_term_xylem.st + 0.096020683 * storage_term_xylem.sh) +
        (1 + common[[c("Rg.N")]]) * storage_term_needles * needle.pot.growth +
        (1 + common[[c("Rg.R")]]) * storage_term_roots * root.pot.growth +
        (1 + common[[c("Rg.N")]]) * storage_term_needles * bud.pot.growth +
        to.mycorrhiza								# Total daily consumption
    }
    tot.growth.sofar <- cumsum(root.tot.growth + height.tot.growth + wall.tot.growth + needle.tot.growth + bud.tot.growth)
    tot.growth <- cumsum(daily.tot.growth)[365]											# Total growth so far + formation of buds
    tot.consumption <- tot.growth + tot.Rm[365] + tot.Rg[365]
    tot.resp <- tot.Rm + tot.Rg		 								# Total carbon used so far to growth and respiration


    # Measurements and output
    for (i in 1 : n.days) {
      export_daily[i + n.days.export, 1] <- NA
      export_daily[i + n.days.export, 2] <- year
      export_daily[i + n.days.export, 3] <- i
      export_daily[i + n.days.export, 4] <- bud.tot.growth[i]
      export_daily[i + n.days.export, 5] <- wall.tot.growth[i]
      export_daily[i + n.days.export, 6] <- needle.tot.growth[i]
      export_daily[i + n.days.export, 7] <- root.tot.growth[i]
      export_daily[i + n.days.export, 8] <- height.tot.growth[i]
      export_daily[i + n.days.export, 9] <- Rg.tot[i]
      export_daily[i + n.days.export, 10] <- Rm.tot[i]
      export_daily[i + n.days.export, 11] <- height.tot[i]
      export_daily[i + n.days.export, 12] <- wall.tot[i]
      export_daily[i + n.days.export, 13] <- if (sperling_model == FALSE) storage[i] else sugar.needles[i] + sugar.roots[i] + sugar.phloem[i] + sugar.xylem.st[i] + sugar.xylem.sh[i] + starch.needles[i] + starch.roots[i] + starch.phloem[i] + starch.xylem.st[i] + starch.xylem.sh[i]
      export_daily[i + n.days.export, 14] <- if (sperling_model == FALSE) sugar[i] else sugar.needles[i] + sugar.roots[i] + sugar.phloem[i] + sugar.xylem.st[i] + sugar.xylem.sh[i]
      export_daily[i + n.days.export, 15] <- if (sperling_model == FALSE) starch[i] else starch.needles[i] + starch.roots[i] + starch.phloem[i] + starch.xylem.st[i] + starch.xylem.sh[i]
      export_daily[i + n.days.export, 16] <- if (sperling_model == FALSE) storage_term[i] else storage_term_needles[i]+storage_term_phloem[i]+storage_term_roots[i] + storage_term_xylem.sh[i] + storage_term_xylem.sh[i]
      export_daily[i + n.days.export, 17] <- to.mycorrhiza[i]
      export_daily[i + n.days.export, 18] <- mycorrhiza.tot[i]
      export_daily[i + n.days.export, 19] <- P[i]
      export_daily[i + n.days.export, 20] <- if (sperling_model == FALSE) to_sugar[i] else NA
      export_daily[i + n.days.export, 21] <- if (sperling_model == FALSE) to_starch[i] else NA
      if (xylogenesis == TRUE) {
        export_daily[i + n.days.export, 22] <- daily.consumption[i]
        export_daily[i + n.days.export, 23] <- ring_width[i]
      } else {
        export_daily[i + n.days.export, 22] <- Daily.H.tot[i]
        export_daily[i + n.days.export, 23] <- Daily.N.tot[i]
      }
      export_daily[i + n.days.export, 24] <- GD.tot[i]
      if (sperling_model == T) {
        export_daily[i + n.days.export, 25] <- sugar.needles[i]
        export_daily[i + n.days.export, 26] <- sugar.phloem[i]
        export_daily[i + n.days.export, 27] <- sugar.xylem.sh[i]
        export_daily[i + n.days.export, 28] <- sugar.xylem.st[i]
        export_daily[i + n.days.export, 29] <- sugar.roots[i]
        export_daily[i + n.days.export, 30] <- starch.needles[i]
        export_daily[i + n.days.export, 31] <- starch.phloem[i]
        export_daily[i + n.days.export, 32] <- starch.xylem.sh[i]
        export_daily[i + n.days.export, 33] <- starch.xylem.st[i]
        export_daily[i + n.days.export, 34] <- starch.roots[i]
      }
      if (xylogenesis == TRUE) {
        export_daily[i + n.days.export, 25] <- n.E.tot[i]
        export_daily[i + n.days.export, 26] <- n.W.tot[i]
        export_daily[i + n.days.export, 27] <- n.M.tot[i]
      }
    }

    n.days.export <- n.days.export + n.days

    # Yearly output data
    export_yearly[n.year, 1] <- year
    export_yearly[n.year, 2] <- if (sperling_model == F) starch[n.days] else starch.needles[n.days] + starch.roots[n.days] + starch.phloem[n.days] + starch.xylem.st[n.days] + starch.xylem.sh[n.days]
    export_yearly[n.year, 3] <- if (sperling_model == F) sugar[n.days] else sugar.needles[n.days] + sugar.roots[n.days] + sugar.phloem[n.days] + sugar.xylem.st[n.days] + sugar.xylem.sh[n.days]
    export_yearly[n.year, 4] <- wall.tot[n.days]
    export_yearly[n.year, 5] <- height.tot[n.days] + B0 * CH * parameters[c("HH0"), c(site)] / 1000
    export_yearly[n.year, 6] <- needle.tot[n.days] + parameters[c("HN0"), c(site)] * repol[[c("m.N")]]
    export_yearly[n.year, 7] <- root.tot[n.days]
    export_yearly[n.year, 8] <- tot.Rm[n.days]
    export_yearly[n.year, 9] <- tot.Rg[n.days]
    export_yearly[n.year, 10] <- tot.P[n.days]
    export_yearly[n.year, 11] <- cumsum(PF)[n.days]
    export_yearly[n.year, 12] <- cum.Daily.H.tot[n.days]
    export_yearly[n.year, 13] <- cum.Daily.N.tot[n.days]
    export_yearly[n.year, 14] <- if (xylogenesis == FALSE) tot.mm[n.days] else ring_width[n.days]  # mm
    export_yearly[n.year, 15] <- needle_mass[n.year]
    export_yearly[n.year, 16] <- sum(needle_cohorts[n.year,])
    if (sperling_model == T) {
      export_yearly[n.year, 17] <- sugar.needles[n.days]
      export_yearly[n.year, 18] <- sugar.phloem[n.days]
      export_yearly[n.year, 19] <- sugar.xylem.sh[n.days]
      export_yearly[n.year, 20] <- sugar.xylem.st[n.days]
      export_yearly[n.year, 21] <- sugar.roots[n.days]
      export_yearly[n.year, 22] <- starch.needles[n.days]
      export_yearly[n.year, 23] <- starch.phloem[n.days]
      export_yearly[n.year, 24] <- starch.xylem.sh[n.days]
      export_yearly[n.year, 25] <- starch.xylem.st[n.days]
      export_yearly[n.year, 26] <- starch.roots[n.days]
    }
    if (xylogenesis == TRUE) {
      export_yearly[n.year, 17] <- ew_width[n.days]            # mm
      export_yearly[n.year, 18] <- ring_density[n.days]        # kg C m-3
      export_yearly[n.year, 19] <- ew.cells_tot[n.days]        # no
      export_yearly[n.year, 20] <- lw.cells_tot[n.days]        # no
    }

    if (sperling_model == T) {
      if (count%%2 == 0) {n.year=n.year + 1}
    } else {
      n.year=n.year + 1
    } # For actual years

    count <- count + 1 # For Sperling
  }   # loop of the years ends

export_daily[, 1] <- format(seq(as.POSIXct(as.character(paste0(years[1], "0101")), format = "%Y%m%d"), as.POSIXct(as.character(paste0(years[n.year-1], "1231")), format = "%Y%m%d"), by = "day"), "%Y-%m-%d")
out <- list(export_daily, export_yearly)
names(out) <- c("Daily", "Yearly")


return(out)

}























