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

  LH.estim = TRUE,				# LH depends on the GPP during previous july-august

  LN.estim = TRUE,				# LN depends on the GPP during previous july-august
  mN.varies = TRUE,				# needle mass (in maintenance respiration) is 2/3 of the total during period 1.10 - 31.5.

  LD.estim = TRUE,				# LD depends on the GPP during March-August
  sD.estim.T.count = FALSE,			# sD depends on the number of days when g in growing window - analogue to needles

  trees_grow = FALSE,				# can be false if mature trees are modelled and not for a very long period
  growth_decreases = FALSE,			# the height and diameter growth (alfa_S and alfaD) decrease during the simulation
  needle_mass_grows = FALSE,		# Is needle mass dynamic i.e. the modelled growth is also respiring etc and following for some years? If true, note that root mass is related to needle mass

  mychorrhiza = TRUE, 			# If allocation to mychorrhiza is taken into account
  root_as_Ding = TRUE,

  sperling_model = FALSE,       # Dynamic sugar model using Sperling's enzyme dynamics TODO: add this as an option
  xylogenesis = FALSE,    # TODO: add the new part from the Lettosuo version

  s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
  s.H0 = 1					# and for shoot grwoth
  ) {

  #####
  ## Input tests!
  #####
  # weather
  if (ncol(weather) != 7) {stop("Incomplete weather data - not enough variables")}
  if (sum(names(weather) == c("date", "T", "P", "TSA", "TSB", "MB", "Rain")) != 7) {stop("Incomplete weather data - incorrect variables, or named incorrectly")}

  # Check that the sites are within the sites allowed
  if ((site %in% c("Hyde", "Lettosuo")) == F) {stop("Unknown site: Please pick between Hyde and Lettosuo")}

  if (sperling_model == TRUE) {if (mychorrhiza == T) {
    mychorrhiza = FALSE
    warning("Mycorrhiza has been changed to mycorrhiza = false as mycorrhiza is included explicitly in the Sperling submodel")}
  }

  # TODO: add tests for the parameter inputs!

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

  export_yearly <- data.frame(matrix(ncol=15, nrow=length(years)))
  export_daily <- data.frame(matrix(ncol=22, nrow=length(years)*365))
  # TODO: decide the most sensible outputs
  # TODO: maybe make outputs depend on the original settings
  names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                            "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts")
  names(export_daily) <- c("year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                            "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage", "to.mycorrhiza", "mycorrhiza.tot",
                            "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot")
  n.days.export <- 0
  n.year <- 1
  n.days.sugar.export <- 0 # TODO: do I need this value?

  #####
  ## Non yearly dependent coefficients
  #####
  # TODO: should these go in the loop anyway for readability or are there other things that can join this?

  # height growth coefficient

  height_growth_coefficient <- diameter_growth_coefficient <- NULL
  height_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("height_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))
  diameter_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("diameter_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))



  if (growth_decreases == TRUE) {
    height_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("height_growth_coefficient_max"),c(site)], ratios[c("height_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
    diameter_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("diameter_growth_coefficient_max"),c(site)], ratios[c("diameter_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
  }

  cell.volume.ew<-(parameters[c("cell.d.ew"), c(site)])^2*parameters[c("cell.l.ew"),c(site)]
  cell.volume.lw<-(parameters[c("cell.d.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)])*parameters[c("cell.l.lw"),c(site)]		# the tangential width of the cell is the same for early and late wood

  wall.volume.ew<-cell.volume.ew-(parameters[c("cell.d.ew"), c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])^2*(parameters[c("cell.l.ew"),c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])	# m3
  wall.volume.lw<-cell.volume.lw-(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.l.lw"),c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])	# m3

  cell.density.ew<-parameters[c("cell.wall.density.ew"),c(site)]*wall.volume.ew/cell.volume.ew		# to calculate the wood density
  cell.density.lw<-parameters[c("cell.wall.density.lw"),c(site)]*wall.volume.lw/cell.volume.lw

  daily.rate.ew<-wall.volume.ew/parameters[c("tau.We"),c(site)]		#m3 cell wall day-1
  daily.rate.lw<-wall.volume.lw/parameters[c("tau.Wl"),c(site)]		#m3 cell wall day-1

  Carbon.daily.rate.ew<-daily.rate.ew*parameters[c("cell.wall.density.ew"),c(site)]	# kg C day-1	carbon to early wood wall formation
  Carbon.daily.rate.lw<-daily.rate.lw*parameters[c("cell.wall.density.lw"),c(site)]	# kg C day-1	carbon to late wood wall formation


  if (site == "Hyde") {
    stem.no <- cbind(1997 : 2020, rep(1010, length.out = length(1997 : 2020)))	# Photosynthesis calculated with SPP-model for tree class 15-20 cm. Then compared with eddy GPP, determined with which stem no. the portion of eddy GPP is same as SPP estimate.
  }

  # Carbon to height growth
  CH<-parameters[c("density_tree"),c(site)]*parameters[c("carbon_share"),c(site)]
  M.suc<-12*common[[c("M.C")]]+22*common[c("M.H")]+11*common[[c("M.O")]]


  #####
  ## Year loop
  #####

  LAI <- needle_mass <- NULL

  for (year in years) {

    n.days <- if (year %in% leap_years) 366 else 365

    if (n.year == 1) {
      rep = repola(parameters[c("D0"), c(site)], parameters[c("h0"), c(site)], n.year, needle_mas = NULL, ste = site, params = parameters, reps = repo)
      needle_mass[1] <- rep[[c("needle_mass")]]
    } else {
      rep = repola(parameters[c("D0"), c(site)], parameters[c("h0"), c(site)], n.year, needle_mas = needle_mass, ste = site, params = parameters, reps = repo)
      if (needle_mass_grows==FALSE) {
        needle_mass[n.year]=needle_mass[n.year-1] 		# constant needle mass, needlemass determined in sitespecific parameters (e.g. parameters_hyde.r)
      } else {
        needle_mass <- rep[[c("needle_mass")]] # calculated using biomass equations, reached height and diameter
      }
    }

    needle_cohorts <- matrix(ncol = parameters[c("n_age"), c(site)], nrow = length(years))				# Needle age classes (assumed three classes as in central Finland)

    # Needle age classes: oldest fall off, the others age with one year
    # The youngest needle class is the (needle_length[year-]/average_length) * 1/3*(estimated needle mass by the biomass equations)
    if (n.year > 1) {
      needle_cohorts[n.year, 1] <- cum.Daily.N.tot[365] / parameters[c("n_lenght"), c(site)] * needle_mass[n.year] / 3
      needle_cohorts[n.year, 2] <- needle_cohorts[n.year-1, 1]
      needle_cohorts[n.year, 3] <- needle_cohorts[n.year-1, 2]
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

    # Initalising the basic values for these variables
    # TODO: are they in the right place?
    B0<-pi/4*parameters[c("D0"), c(site)]^2		# basal area in the beginning
    h00<-parameters[c("h0"), c(site)]
    D00 <- parameters[c("D0"), c(site)]


    ### Photosynthesis
    # Photosynthesis of one tree per day (g C / m2 --> kg C / tree)
    P <- tot.P <- NULL

    # TODO: where does stem.no come from?
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

    # TODO: height_growth_coefficient replace this
    LH <- parameters[c("LH0"), c(site)] * height_growth_coefficient[which(height_growth_coefficient[,1] == year), 2]
    if (LH.estim == TRUE) LH <- LH * GPP_previous_sum[which(GPP_previous_sum[,1] == year), 2] / mean(GPP_previous_sum[,2])


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

    # Daily potential length growth (mm/day)  and length (mm)  of needles
    GN <- g * fN * LN
    HN <- parameters[c("HN0"), c(site)] + cumsum(GN)

    # Carbon to needle growth per day (if there's no carbon limitation) (kg C / day)
    needle.pot.growth <- rep[[c("m.N")]] * GN * (max(HH) / parameters[c("h_increment"), c(site)])		# Simulated shoot length/average shoot length

    ## Root growth
    sR <- fR <- GR <- root.pot.growth <- NULL

    if (site == "Hyde") {
      LR <- LR0 <- 2 * rep[[c("m.R.tot")]] / parameters[c("sRc"), c(site)] # root growth rate
    } else {
      LR <- LR0 / parameters[c("root.lifetime"),c(site)]
    }

    gR <- (Tsb > common[[c("TR0")]]) * (1 / (1 + exp(-common[[c("a")]] * ((Tsb - common[[c("TR0")]]) - common[[c("b")]])))) * (1 - 1 / exp(M.soil * 10)) # Temp and M driving the phase of the annual cycle of root growth

    sR <- parameters[c("sR0"), c(site)] + cumsum(gR)										# The phase of the annual cycle of root growth

    fR <- (sR < parameters[c("sRc"), c(site)]) * (sR > 0) * (sin(2 * pi / parameters[c("sRc"), c(site)] * (sR - parameters[c("sRc"), c(site)] / 4)) + 1) / 2	# A function driven by the phase of the annual cycle (annual pattern of growth) [0,1]

    if(root_as_Ding == TRUE) {
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
    }

    # Carbon to root growth per day (if there's no carbon limitation) (kg C / day)
    GR <- fR * LR * gR

    root.pot.growth <- GR								# Carbon to root growth per day (if there's no carbon limitation) (kg C / day)
    root.pot.growth[is.na(root.pot.growth)] <- 0

    ## Diameter growth
    g.sD.T <- g.sD.GPP <- NULL
    sDA <- sD <- fD <- GD <- tot.cells.pot <- NULL

    for (yy in 1:n.days) if(yy<s.D0) g.sD.T[yy] <- 0 else (g.sD.T[yy] <- g[yy])

    sD <- parameters[c("sD0.Trad"), c(site)] + cumsum(g.sD.T)						# The phase of the annual cycle of diameter growth, sD:

    if (sD.estim.T.count == TRUE){
      sDA <- sD0.T.count + cumsum(g.sD.T)
      sD <- cumsum(sDA > 0)
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

    GD <- g.sD.T * fD * LD						# Daily potential number of new cells per day (in one radial cell row), used to be called division
    GD[is.na(GD)] <- 0
    tot.cells.pot <- cumsum(GD)

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

    # n.E.pot<-round(n.E.pot, digits=0)
    # n.W.pot<-round(n.W.pot, digits=0)
    # n.M.pot<-round(cumsum(GD) - n.W.pot - n.E.pot, digits=0)

    # carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
    CE.ew <- CE.lw <- NULL
    for (i in 1 : n.days){
      CE.ew[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.ew"), c(site)] / 2)^2 * (parameters[c("cell.l.ew"), c(site)]) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
      CE.lw[i] <- common[[c("osmotic.sugar.conc")]] * pi * (parameters[c("cell.d.lw"), c(site)] / 2)^2 * (parameters[c("cell.l.lw"), c(site)]) * M.suc / (common[[c("gas.const")]] * (Temp[i] + common[[c("abs_zero")]]) * 1000) * 12 * common[[c("M.C")]] / M.suc / tau.E[i]
    }
    CE.ew <- unlist(CE.ew)
    CE.lw <- unlist(CE.lw)
    # The number of forming cell rows in the tree
    # TODO: in the orginal hyde version cell.l was one value, in the new parameters I have seperated it into lw and ew so the value of this formula is still the same, but it should be
    # combined with the Lettosuo growth method in a more sensible way
    n.rows <- ratios[1,c(site)] * parameters[c("h0"), c(site)] / parameters[c("cell.l.ew"),c(site)] * pi * parameters[c("D0"), c(site)] / parameters[c("cell.d.ew"), c(site)]

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
    CW <- rep(1.8 * 10^-11, length.out = n.days)

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

    ## New bud growth

    sB <- fB <- GB <- CB <- bud.pot.growth <- bud.pot <- NULL

    # Factor driving the phase of the annual cycle of buds (no. of days after needle growth onset).
    day.no <- 1 : n.days
    sB <- cumsum(day.no>parameters[c("sB0"), c(site)])								# The phase of the annual cycle of bud growth

    fB <- (sB < parameters[c("sBc"), c(site)]) * (sin(2 * pi / parameters[c("sBc"), c(site)] * (sB - parameters[c("sBc"), c(site)] / 4)) + 1) / 2		# A function driven by the phase of the annual cycle (annual pattern of growth) [0,1]

    bud.pot.growth <- g * fB * parameters[c("LB"), c(site)]
    bud.pot <- cumsum(bud.pot.growth)

    ############### Maintenance respiration (without carbon limitation) ##########################
    RmN.a <- RmS.a <- RmR.a <- Rm.a <- RmN <- RmS <- RmR <- NULL

    m.S.tot <- ratios[1,c(site)] * (B0) * (parameters[c("h0"), c(site)]) * parameters[c("density_tree"),c(site)] * parameters[c("carbon_share"),c(site)]		# woody carbon mass
    Ra.share <- -0.0007 * Tsa^2 + 0.0424 * Tsa + 0.0273					# share of autotrophic soil respiration
    Ra.share[Ra.share<0] = 0

    RmS.a <- parameters[c("Rm.S"),c(site)] * (exp(log(parameters[c("Q10.S"),c(site)]) / 10 * (Temp)) - 0.7) * m.S.tot * ratios[c("Q10.S"),c(site)]			# Maintenance respiration of wood
    RmR.a <- parameters[c("Rm.R"),c(site)] * (exp(log(parameters[c("Q10.R"),c(site)]) / 10 * (Tsa)) - exp(-log(parameters[c("Q10.R"),c(site)]) / 2)) * rep[[c("m.R.tot")]] * Ra.share	# Maintenance respiration of roots
    RmN.a <- parameters[c("Rm.N"),c(site)] * (exp(log(parameters[c("Q10.N"),c(site)]) / 10 * (Temp)) - 0.7) * rep[[c("m.N.tot")]]						# Maintenance respiration of needles

    if (mN.varies == TRUE){
      m.N.tot2 = NULL
      for(yy in 1 : n.days) {
        m.N.tot2[yy] = if(yy > 150 & yy < 285) rep[[c("m.N.tot")]] else rep[[c("m.N.tot")]] * 2/3
      }
      RmN.a <- parameters[c("Rm.N"),c(site)] * (exp(log(parameters[c("Q10.N"),c(site)]) / 10 * (Temp)) - 0.7) * m.N.tot2
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

    if (sperling_model == TRUE) {
      W.crit.needles <- sperling[c("starch.needles0"), c(site)] + sperling[c("sugar.needles0"), c(site)]
      W.crit.phloem <- sperling[c("starch.phloem0"), c(site)] + sperling[c("sugar.phloem0"), c(site)]
      W.crit.roots <- sperling[c("starch.roots0"), c(site)] + sperling[c("sugar.roots0"), c(site)]
      W.crit.xylem.sh <- sperling[c("starch.xylem.sh0"), c(site)] + sperling[c("sugar.xylem.sh0"), c(site)]
      W.crit.xylem.st <- sperling[c("starch.xylem.st0"), c(site)] + sperling[c("sugar.xylem.st0"), c(site)]

      a.k.needles <- 1/(1-1/exp(sperling[c("alfa.needles"),c("site")]*(sperling[c("sugar.needles0"), c(site)]+sperling[c("starch.needles0"), c(site)])))
      a.k.phloem <- 1/(1-1/exp(sperling[c("alfa.phloem"),c("site")]*(sperling[c("sugar.phloem0"), c(site)]+sperling[c("starch.phloem0"), c(site)])))
      a.k.roots <- 1/(1-1/exp(sperling[c("alfa.roots"),c("site")]*(sperling[c("sugar.roots0"), c(site)]+sperling[c("starch.roots0"), c(site)])))
      a.k.xylem.sh <- 1/(1-1/exp(sperling[c("alfa.xylem.sh"),c("site")]*(sperling[c("sugar.xylem.sh0"), c(site)]+sperling[c("starch.xylem.sh0"), c(site)])))
      a.k.xylem.st <- 1/(1-1/exp(sperling[c("alfa.xylem.st"),c("site")]*(sperling[c("sugar.xylem.st0"), c(site)]+sperling[c("starch.xylem.st0"), c(site)])))
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
      print("Add submodel")
    }

    ########### Total growth and carbon consumption

    #  Occurred growth kg C day-1 (potential growth * storage effect)
    for(i in 1 : n.days) {
      root.tot.growth[i] <- storage_term[i] * root.pot.growth[i]
      height.tot.growth[i] <- storage_term[i] * height.pot.growth[i]
      needle.tot.growth[i] <- storage_term[i] * needle.pot.growth[i]
      wall.tot.growth[i] <- storage_term[i] * wall.pot.growth[i]
      bud.tot.growth[i] <- storage_term[i] * bud.pot.growth[i]
      GD.tot[i] <- storage_term[i] * GD[i]
      Rm.tot[i] <- storage_term_Rm[i] * Rm.a[i]
      RmR.tot[i] <- storage_term_Rm[i] * RmR[i]
    }

    # Total occurred growth so far kg C

    root.tot <- cumsum(root.tot.growth)
    height.tot <- cumsum(height.tot.growth)
    needle.tot <- cumsum(needle.tot.growth)
    wall.tot <- cumsum(wall.tot.growth)
    bud.tot <- cumsum(bud.tot.growth)

    # Dimensional growth
    Daily.H.tot <- height.tot.growth / (B0 * CH / 1000 * ratios[1,c(site)])		# mm day-1
    cum.Daily.H.tot <- parameters[c("HH0"), c(site)] + cumsum(Daily.H.tot)					# mm
    Daily.N.tot <- needle.tot.growth / (rep[[c("m.N")]] * max(cum.Daily.H.tot) / parameters[c("h_increment"), c(site)])		# mm day-1
    cum.Daily.N.tot <- parameters[c("HN0"), c(site)] + cumsum(Daily.N.tot)							# mm

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

    # Carbon to mychorrhiza
    to.mycorrhiza = (Tsb > 10) * (storage > optimal.level.myco) * P * parameters[c("growth.myco"),c(site)]
    mycorrhiza.tot = cumsum(to.mycorrhiza)

    # Carbon used for respiration
    Rg.tot <- tot.Rm <- tot.Rg <- NULL
    Rg.tot <- common[[c("Rg.N")]] * needle.tot.growth + common[[c("Rg.R")]] * root.tot.growth + common[[c("Rg.S")]] * height.tot.growth + common[[c("Rg.S")]] * wall.tot.growth
    Rg.root <- common[[c("Rg.R")]] * root.tot.growth
    tot.Rm <- cumsum(Rm.tot)
    tot.Rg <- cumsum(Rg.tot)

    parameters[c("starch0"), c(site)] <- starch[n.days]		# initial value of starch for the following year
    parameters[c("sugar0"), c(site)] <- sugar[n.days]		# the initial value for labile sugars for the following year (~optimal.level)

    # New tree dimensions if trees grow
    if (trees_grow == TRUE) {
      parameters[c("D0"), c(site)] <- parameters[c("D0"), c(site)] + tot.mm[n.days] * 2 / 1000			# The diameter in the begining of next year
      B0 <- pi / 4 * parameters[c("D0"), c(site)]^2							# Basal area
      parameters[c("h0"), c(site)] <- parameters[c("h0"), c(site)] + (cum.Daily.H.tot[n.days] / 1000)	# The top shoot probably grows faster than the average of other shoots
    }

    #### Total daily growth, total growth so far, total carbon consumption and total respiration (kg C day-1 or kg C)
    daily.tot.growth <- daily.tot <- daily.consumption <- tot.consumption <- NULL
    daily.tot.growth <- height.tot.growth + wall.tot.growth + needle.tot.growth + root.tot.growth + bud.tot.growth  	# Total daily growth
    daily.consumption <- storage_term_Rm * Rm.a +
      (1 + common[[c("Rg.S")]]) * storage_term * height.pot.growth +
      (1 + common[[c("Rg.S")]]) * storage_term * wall.pot.growth +
      (1 + common[[c("Rg.N")]]) * storage_term * needle.pot.growth +
      (1 + common[[c("Rg.R")]]) * storage_term * root.pot.growth +
      (1 + common[[c("Rg.N")]]) * storage_term * bud.pot.growth +
      to.mycorrhiza								# Total daily consumption
    tot.growth.sofar <- cumsum(root.tot.growth + height.tot.growth + wall.tot.growth + needle.tot.growth + bud.tot.growth)
    tot.growth <- cumsum(daily.tot.growth)[365]											# Total growth so far + formation of buds
    tot.consumption <- tot.growth + tot.Rm[365] + tot.Rg[365]
    tot.resp <- tot.Rm + tot.Rg		 								# Total carbon used so far to growth and respiration

    # Measurements and output
    for (i in 1 : 365) {
      export_daily[i + n.days.export, 1] <- year
      export_daily[i + n.days.export, 2] <- i
      export_daily[i + n.days.export, 3] <- bud.tot.growth[i]
      export_daily[i + n.days.export, 4] <- wall.tot.growth[i]
      export_daily[i + n.days.export, 5] <- needle.tot.growth[i]
      export_daily[i + n.days.export, 6] <- root.tot.growth[i]
      export_daily[i + n.days.export, 7] <- height.tot.growth[i]
      export_daily[i + n.days.export, 8] <- Rg.tot[i]
      export_daily[i + n.days.export, 9] <- Rm.tot[i]
      export_daily[i + n.days.export, 10] <- height.tot[i]
      export_daily[i + n.days.export, 11] <- wall.tot[i]
      export_daily[i + n.days.export, 12] <- storage[i]
      export_daily[i + n.days.export, 13] <- sugar[i]
      export_daily[i + n.days.export, 14] <- starch[i]
      # TODO: what is this output?
      if (sperling_model == TRUE) export_daily[i + n.days.export, 15] <- storage_term_needles[i]+storage_term_phloem[i]+storage_term_roots[i]
      if (sperling_model == FALSE) export_daily[i + n.days.export, 15] <- storage_term[i]
      export_daily[i + n.days.export, 16] <- to.mycorrhiza[i]
      export_daily[i + n.days.export, 17] <- mycorrhiza.tot[i]
      export_daily[i + n.days.export, 18] <- P[i]
      export_daily[i + n.days.export, 19] <- to_sugar[i]
      export_daily[i + n.days.export, 20] <- to_starch[i]
      export_daily[i + n.days.export, 21] <- Daily.H.tot[i]
      export_daily[i + n.days.export, 22] <- Daily.N.tot[i]
      export_daily[i + n.days.export, 23] <- GD.tot[i]
    }
    n.days.export <- n.days.export + 365


    # Yearly output data
    export_yearly[n.year, 1] <- year
    export_yearly[n.year, 2] <- starch[365]
    export_yearly[n.year, 3] <- sugar[365]
    export_yearly[n.year, 4] <- wall.tot[365]
    export_yearly[n.year, 5] <- height.tot[365] + B0 * CH * parameters[c("HH0"), c(site)] / 1000
    export_yearly[n.year, 6] <- needle.tot[365] + parameters[c("HN0"), c(site)] * rep[[c("m.N")]]
    export_yearly[n.year, 7] <- root.tot[365]
    export_yearly[n.year, 8] <- tot.Rm[365]
    export_yearly[n.year, 9] <- tot.Rg[365]
    export_yearly[n.year, 10] <- tot.P[365]
    export_yearly[n.year, 11] <- cumsum(PF)[365]
    export_yearly[n.year, 12] <- cum.Daily.H.tot[365]
    export_yearly[n.year, 13] <- cum.Daily.N.tot[365]
    export_yearly[n.year, 14] <- tot.mm[365]
    export_yearly[n.year, 15] <- needle_mass[n.year]
    export_yearly[n.year, 16] <- sum(needle_cohorts[n.year,])

    n.year=n.year + 1

  }   # loop of the years ends

  out <- list(export_daily, export_yearly)
  names(out) <- c("Daily", "Yearly")

return(out)

}
