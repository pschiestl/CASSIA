model_xylogenesis = function(g.sD.T, LD, fD, storage_reduction) {
  # TODO: make this function operational with all inputs

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
      ew_lw[i] = ifelse(sD[division_day] < (sDc^2 / Uggla), 1, 2)
      tau.E <- ifelse(sD[division_day] < (sDc^2 / Uggla), round(tau.Ee, 0), round(tau.El, 0))
      tau.W <- ifelse(sD[division_day] < (sDc^2 / Uggla), round(tau.We, 0), round(tau.Wl, 0))
      cell.wall.density[i] <- ifelse(sD[division_day] < (sDc^2 / Uggla), cell.wall.density.ew, cell.wall.density.lw)
      cell.volume.growth.per.day = ifelse(sD[division_day] < (sDc^2 / Uggla), cell.volume.growth.per.day.ew, cell.volume.growth.per.day.lw)
      cell.wall.volume.growth.per.day = ifelse(sD[division_day] < (sDc^2 / Uggla), cell.wall.volume.growth.per.day.ew, cell.wall.volume.growth.per.day.lw)
      for (d in 1:tau.E) {    # time of enlargement
        cell_volume[(division_day + d), i] = cell_volume[(division_day + d - 1), i] + cell.volume.growth.per.day * soil_moisture_effect[division_day + d] * storage_reduction[division_day + d]
        cell_wall_volume[(division_day + d), i] = 0
        carbon_to_enlargement[(division_day + d), i] = cell.volume.growth.per.day * osmotic.sugar.conc * M.suc / (gas.const * (Temp[division_day + d] + abs_zero) * 1000) * 12 * M.C / M.suc  * storage_reduction[division_day + d]
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
    tau.E <- (sD < sDc^2 / Uggla) * (tau.Ee) + (sD >= sDc^2 / Uggla) * (tau.El)       # Duration of cell enlargement of earlywood cells and late wood cells
    tau.W <- (sD < sDc^2 / Uggla) * (tau.We) + (sD >= sDc^2 / Uggla) * (tau.Wl)       # Duration of cell wall formation of earlywood cells and late wood cells

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
      CE.ew[i] <- osmotic.sugar.conc * pi * (cell.d.ew / 2)^2 * (cell.l) * M.suc / (gas.const * (Temp[i] + abs_zero) * 1000) * 12 * M.C / M.suc / tau.E[i]
      CE.lw[i] <- osmotic.sugar.conc * pi * (cell.d.lw / 2)^2 * (cell.l) * M.suc / (gas.const * (Temp[i] + abs_zero) * 1000) * 12 * M.C / M.suc / tau.E[i]
    }

    carbon.enlargement <- carbon.wall <- en.growth <- en.release <- wall.growth <- d.growth <- NULL
    # The carbon needed to enlarge the cells  (kg C /day) (in one radial cell row)
    carbon.enlargement <- (sD < sDc^2 / Uggla) * CE.ew * n.E + (sD >= sDc^2 / Uggla) * CE.lw * n.E
    # Carbon to enlargement per day (kg C /day) (in the whole tree)
    en.growth <- n.rows * carbon.enlargement

    for(i in 1:round(tau.Ee)) en.release[i] <- 0	# the carbon used in enlargement is released after some days.
    for(i in (round(tau.Ee) + 1) : n.days) en.release[i] <- en.growth[i-round(tau.Ee)]

    # Carbon to wall formation
    # Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!
    CW <- (sD < sDc^2 / Uggla) * Carbon.daily.rate.ew + (sD >= sDc^2 / Uggla) * Carbon.daily.rate.lw
    CW <- rep(2.9 * 10^-11, length.out = n.days)

    # The use of carbon to wall growth kg C per day
    wall.growth <- n.rows * CW * n.W

    en.growth[is.na(en.growth)] <- 0
    en.release[is.na(en.release)] <- 0
    wall.growth[is.na(wall.growth)] <- 0

    # Total use of carbon to wall growth so far kg C
    wall.tot <- cumsum(wall.growth)

    ew_cells <- (sD < sDc^2 / Uggla) * tot.cells
    ew_cells[sD > sDc^2 / Uggla] = max(ew_cells)
    lw_cells <- tot.cells - ew_cells

    # Calculation of the daily annual ring width, starting from 16.12.2013.
    # Note! This does not include the size of enlarging cells, only wall forming and mature
    diameter_cells = n.M + n.W
    diameter_ew_cells = ifelse(diameter_cells <= ew_cells, diameter_cells, max(ew_cells))
    diameter_lw_cells = diameter_cells - diameter_ew_cells
    ew_width <- diameter_ew_cells * cell.d.ew * 1000 #+ ifelse(sD < sDc^2 / Uggla, n.E * cell.d.ew / 2)
    lw_width <- diameter_lw_cells * cell.d.lw * 1000 #+ ifelse(sD >= sDc^2 / Uggla, n.E * cell.d.lw / 2)
    ring_width = ew_width + lw_width
    cell_wall_thickness = c(rep(wall.thickness.ew, round(ew_cells[365],0)), rep(wall.thickness.lw, round(lw_cells[365], 0)))
    cell_d_final = c(rep(cell.d.ew, round(ew_cells[365],0)), rep(cell.d.lw, round(lw_cells[365], 0)))
    cell_density = c(rep(cell.density.ew, round(ew_cells[365],0)), rep(cell.density.lw, round(lw_cells[365], 0)))
    ring_density = cumsum(CW * n.W)[365] / (ew_cells[365] * cell.d.ew^2 * cell.l.ew + lw_cells[365] * cell.d.lw^2 * cell.l.lw)

    ew_cells = max(ew_cells)
    lw_cells = max (lw_cells)
  }

  list_xylogenesis = list("GD" = GD, "tot.cells" = tot.cells, "n.E" = n.E, "n.W" = n.W, "n.M" = n.M, "ew_cells" = ew_cells, "lw_cells" = lw_cells,
                          "ring_width" = ring_width, "ew_width" = ew_width, "lw_width" = lw_width, "en.growth" = en.growth,
                          "en.release" = en.release, "wall.growth" = wall.growth, "wall.tot" = wall.tot, "cell_d_final" = cell_d_final,
                          "cell_wall_thickness" = cell_wall_thickness, "cell_density" = cell_density, "ring_density" = ring_density)


}


