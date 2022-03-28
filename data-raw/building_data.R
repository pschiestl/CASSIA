data.direct <- "C:/Users/Käyttäjä/OneDrive - University of Helsinki/Research Assistant/CASSIA/Package/CASSIA/data"

######
## ratios_p
######

ratios_p <- data.frame(matrix(nrow = 9, ncol = 2))
names(ratios_p) <- c("Hyde", "Lettosuo")
# Lettosuo: The multiplier between a cylinder (with diam=D0, height=h0) and the total biomass of stem, coarse roots and branches
ratios_p[1,] <- c(0.6,0.6) # form_factor
# Lettosuo: Helmisaari et al. 2006 Tree physiology -> 2.0 for VT, 3.8 for MT, 5.7 for OMT. Very nice curve for needles / fine roots vs. fine root N %
# leads now to ~ 100 gC m-2 (roots < 2 mm) but result lower than in Leppälammi-Kujansuu et al. (2013, Plant Soil) where they found ca 225 gC m-2 (roots < 2 mm) in control and 300-350 in fertilized
# ca 225 gC m-2 (roots < 2 mm) in control and 300-350 in fertilized (Leppälammi-Kujansuu et al. 2013, Plant Soil)
ratios_p[2,] <- c(NA,1/2.9) # needle_fineroot_ratio
# depends on tree size, specicies and site
ratios_p[3,] <- c(0.8, 0.8) # sapwood.share
ratios_p[4,] <- c(4.3, 1) # height_growth_coefficient repeated value
ratios_p[5,] <- c(1.6, 1) # diameter_growth_coefficient repeated value
ratios_p[6,] <- c(5.5, 1.28) # height_growth_coefficient repeated value max if gorwth decreases
ratios_p[7,] <- c(3.8, 0.88) # height_growth_coefficient repeated value min if gorwth decreases
ratios_p[8,] <- c(1.9, 1.19) # diameter_growth_coefficient repeated value max if gorwth decreases
ratios_p[9,] <- c(1.5, 0.94) # diameter_growth_coefficient repeated value min if gorwth decreases
row.names(ratios_p) <- c("form_factor", "needle_fineroot_ratio", "sapwood.share",
                       "height_growth_coefficient", "diameter_growth_coefficient",
                       "height_growth_coefficient_max", "height_growth_coefficient_min",
                       "diameter_growth_coefficient_max", "diameter_growth_coefficient_min")
save(ratios_p, file = paste0(data.direct, "/ratios_p.RData"))

######
## GPP approximations - same for both sites, but GPP ref is the vairable name for Lettosuo
######

# The GPP sum in previous July-August (1996-2019 as modelling years are 1997-2020) to be used in LS-estimation. (g C / m2)(sum of the two months) NOTE! 1996 is missing, thus the value is app. the mean of all time
GPP_previous_sum <- cbind(1997:2020, c(460, 497.1, 431.3, 410.3, 451.6, 473.7, 453.4, 437.3, 448.0, 436.3, 397.0, 461.4, 420.4, 500.7, 476.2, 518.1, 529.0, 483.1, 481.3,  494.5,  463.6,  484.4, 447.9, 476.6))
# Daily GPP in 2009 (g C / m2) (1.1. - 31.12.), used in LD-estimation (new cell). E.g. daily average GPP could be used as well.
GPP_ref <- c(0.07,0.06,0.00,0.00,0.00,0.00,0.00,0.00,0.05,0.12,0.04,0.00,0.00,0.12,0.07,0.00,0.00,0.00,0.12,0.00,0.05,0.08,0.00,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.06,0.07,0.08,0.12,0.00,0.00,0.02,0.00,0.01,0.00,0.00,0.00,0.10,0.00,0.07,0.01,0.07,0.06,0.20,0.12,0.00,0.00,0.06,0.01,0.00,0.08,0.16,0.14,0.05,0.31,0.16,0.06,0.13,0.19,0.08,0.00,0.02,0.17,0.21,0.36,0.41,0.41,0.24,0.20,0.11,0.45,0.27,0.23,0.15,0.08,0.27,0.11,0.40,0.15,0.28,0.15,0.44,0.21,0.65,0.72,1.19,0.47,1.06,1.62,1.49,1.01,1.42,1.75,1.16,1.39,1.49,1.05,1.14,1.12,1.04,1.17,0.84,1.06,1.95,1.92,1.99,2.84,3.60,3.68,4.06,4.62,4.03,3.87,4.92,3.07,4.45,3.46,5.08,4.48,3.33,5.20,4.21,5.32,4.67,6.14,5.25,6.19,5.73,5.25,4.21,5.26,5.50,7.52,6.48,7.31,7.24,6.67,6.08,8.29,7.06,6.74,6.59,8.83,9.85,6.74,1.89,5.62,6.57,8.35,8.29,8.88,8.53,7.65,8.65,6.79,8.35,4.17,4.97,8.53,9.52,5.66,4.29,9.83,8.57,8.28,9.10,7.66,7.95,8.27,8.19,10.59,10.13,10.29,10.21,9.59,9.61,10.04,8.60,8.77,7.17,4.58,8.78,9.39,7.77,6.32,9.12,9.44,8.34,9.63,9.19,11.36,4.01,10.24,8.99,9.64,9.19,7.23,9.35,10.18,9.72,9.40,9.42,8.00,7.97,8.40,7.90,9.93,9.70,9.77,9.96,9.00,9.50,9.31,8.27,3.40,7.94,8.37,8.37,3.95,6.32,7.10,7.90,7.93,7.72,8.08,6.24,6.27,6.47,3.99,5.82,5.09,4.33,6.66,5.44,2.38,5.53,6.69,4.46,5.98,3.28,5.33,2.86,5.33,5.98,5.28,5.17,3.88,4.95,5.46,3.84,4.77,4.01,3.83,4.84,3.47,2.38,3.53,4.04,3.13,3.51,2.91,3.17,2.80,3.35,3.77,2.41,2.61,1.92,1.59,2.77,1.10,2.60,2.14,2.15,1.34,1.36,1.90,1.24,1.20,0.90,0.68,0.88,0.83,1.23,1.57,0.84,0.82,0.77,0.62,0.32,0.31,0.61,1.13,0.12,0.64,0.17,0.80,0.61,0.47,0.54,0.06,0.20,0.08,0.00,0.00,0.00,0.02,0.00,0.14,0.15,0.18,0.26,0.18,0.16,0.10,0.04,0.40,0.15,0.15,0.18,0.13,0.20,0.09,0.17,0.16,0.03,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.04,0.03,0.04,0.04,0.00,0.00,0.06,0.00,0.09,0.04,0.00,0.04)

save(GPP_previous_sum, file = paste0(data.direct, "/GPP_previous_sum.RData"))
save(GPP_ref, file = paste0(data.direct, "/GPP_ref.RData"))

######
## parameters_p
######
parameters_p <- data.frame(matrix(nrow = 30, ncol = 2))
names(parameters_p) <- c("Hyde", "Lettosuo")
#### Respiration
# Needles Q10
parameters_p[1,] <- c(1.898,1.898)  # Q10.N
# Needles R0
parameters_p[2,] <- c(0.00267, 0.0020) # Rm.N
# Wood (stem + brances + transport roots) Q10
parameters_p[3,] <- c(1.74788, 1.74788) # Q10.S
# Wood (stem + brances + transport roots) R0
parameters_p[4,] <- c(5.5576 * 10^-5, 5.5576 * 10^-5) # Rm.S
# Fine roots Q10
parameters_p[5,] <- c(2.5575, 2.9662) # Q10.R
# Fine roots R0
parameters_p[6,] <- c(0.00958, 0.0059) # Rm.R

#### Gorwth
# Lettosuo: (as for pine, needs to be changed!)
# root growth beginning
parameters_p[7,] <- c(-2.5, -2.5) # sR0
# root growth cessation
parameters_p[8,] <- c(30, 30) # sRc

#### Mycorrhiza
parameters_p[9,] <- c(0.1, 0.1) # growth.myco
# years
# Lettosuo: years, control site, 0.92 warmed, fertilized, irrigated - depends on nutrient availability, soil temperature etc. (Leppälammi-Kujansuu et al. 2014, Plant and soil) -affects probably also root mass!
parameters_p[10,] <- c(1.7, 2.01) # root.lifetime

#### Shoots
# Lettosuo: parameters_p from 2015 and 2016, temperature sum calculated from 1.1.
# Hyde: parameters_p from 2008, temperature sum calculated from 1.1.
parameters_p[11,] <- c(10, 10) # HH0
parameters_p[12,] <- c(-1.359200388, -4.12008) # sH0
parameters_p[13,] <- parameters_p[14,] <- c(8.226401284, NA) # LH, LH0
parameters_p[15,] <- c(14.59636279, 12.8328) # sHc

#### Needles
# Lettosuo: parameters_p from 2015 and 2016
# Hyde: parameters_p from 2008
parameters_p[16,] <- c(-8.37584, -3.56589) # sN0
# TODO: LN is not a constant in the code
parameters_p[17,] <- parameters_p[18,] <- c(1.849493, NA) # LN, LN0 n_lenght/(0.5*sNc)
parameters_p[19,] <- c(5.263883, 7.60671) # sNc
parameters_p[20,] <- c(1, 1) # HN0

#### Diameter
# Lettosuo: parameters_p from 2007, 2008 and 2009, L not correct!
parameters_p[21,] <- c(-3.724083738, -3.5) # sD0.Trad
parmaeters_p[22,] <- parameters_p[23,] <- c(1.293443902, 1.9) # LD, LD0
parameters_p[24,] <- c(5.077004992, 5.2) # sDc
parameters_p[25,] <- c(NA, 8.8) # sDc.T.count

# Lettosuo: Duration parameters_p, determined with year 2009 data C:/LocalData/schiestl/CASSIA/CASSIA_spruce/shoot_measurements/kaikki.xlsx
# duration of early wood cell enlargement,
parameters_p[26,] <- c(10.68685877, 5.5) # tau.Ee
# duration of late wood cell enlargement
parameters_p[27,] <- c(8.789131263, 4.8) # tau.El
# duration of early wood cell wall formation
parameters_p[28,] <- c(25.29448857, 17.8) # tau.We
# duration of late wood cell wall formation
parameters_p[29,] <- c(35.12148687, 19.2) # tau.Wl

# GPP effect on the daily value of parameter LD
parameters_p[30,] <- c(5, 5) # tau.GPP
# Division to early wood and late wood
# Lettosuo: Probably varies with site fertility (Lundgren et al. 2004 Cell Wall thickness and tangential and radial... Silva Fennica). Check late wood shares also from Jyske et al. 2008 Silva Fennica: Wood density within Norway spruce stems, latewood share should be 12-22 % (high in outer rings)
parameters_p[31,] <- c(1.95, 1.8) # Uggla

#### Buds
# Bud growth begins at 20.6.
parameters_p[32,] <- c(171, 171) # sB0 TODO: make this work with the sugar addition
parameters_p[33,] <- c(85, 85) # sBc
parameters_p[34,] <- c(0.005, 0.005) # LB


#### Xylogenesis

#### Cell diameters
# early wood cell diameter Havimo, M. et al. 2008, Silva Fennica (32.4 in Saren et al. 2001: Structural variation of tracheid... Journal of structural biology)
parameters_p[35,]<-c(35.7*10^-6, 32.1*10^-6) # cell.d.ew
# late wood cell diameter Havimo, M. et al. 2008, Silva Fennica
parameters_p[36,]<-c(24.2*10^-6, 27.5*10^-6) # cell.d.lw
# cell length of early wood (m) (M?kinen et al. 2008)   # Average of measured, varies between rings - may vary between provenances. (3.1 in Saren et al. 2001: Structural variation of tracheid... Journal of structural biology)
parameters_p[37,]<- c(2.59*10^-3, 2.89*10^-3) # cell.l.ew
# cell length of late wood (m) (M?kinen et al. 2008)   # Average of measured, varies between rings - may vary between provenances.
parameters_p[38,]<- c(2.73*10^-3, 2.97*10^-3) # cell.l.lw

# kg C m-3 , 1500/1114kg m-3, Relationships between the intra-ring wood density assessed by X-ray densitometry and optical anatomical measurements in conifers. Consequences for the cell wall apparent density determination
parameters_p[39,]<- c(114/2, 1140/2) # cell.wall.density.ew
# kg C m-3 , 1500/1114kg m-3, Relationships between the intra-ring wood density assessed by X-ray densitometry and optical anatomical measurements in conifers. Consequences for the cell wall apparent density determination
parameters_p[40,]<-c(114/2, 1360/2) 	# cell.wall.density.lw

# TODO: cell.wall.desity is only defined for all in pine

# Havimo, M. et al. 2008, Silva Fennica (3.1 in Saren et al. 2001: Structural variation of tracheid... Journal of structural biology)
parameters_p[41,]<-c(2.61*10^-6,3.1*10^-6) # wall.thickness.ew
# Havimo, M. et al. 2008, Silva Fennica
parameters_p[42,]<-c(5.23*10^-6,3.88*10^-6) # wall.thickness.lw

# Comes from cell volume / tauE
# TODO: this is part of the equations for pine, should check this for spruce as well
parameters_p[43,] <- c(NA, 5.49*10^-13) # cell.volume.growth.per.day.ew
# Comes from cell volume / tauE
parameters_p[44,] <- c(NA,4.62*10^-13)  # cell.volume.growth.per.day.lw

parameters_p[45,] <- c(400, 400)	# kg m-3, Scots pine, density_tree
parameters_p[46,] <- c(0.5, 0.5)	# kg kg-1, carbon_share

# Both: Diameter and height in of storage tree in 2015
parameters_p[47,] <- c(0.175, 0.175) # D0
parameters_p[48,] <- c(17.9, 17.9) #h0
# Hyde: needles stay alive for three years
# Lettosuo: Helmisaari et al. 2006 Tree physiology -> 2.0 for VT, 3.8 for MT, 5.7 for OMT. Very nice curve for needles / fine roots vs. fine root N %
# leads now to ~ 100 gC m-2 (roots < 2 mm) but result lower than in Leppälammi-Kujansuu et al. (2013, Plant Soil) where they found ca 225 gC m-2 (roots < 2 mm) in control and 300-350 in fertilized
# ca 225 gC m-2 (roots < 2 mm) in control and 300-350 in fertilized (Leppälammi-Kujansuu et al. 2013, Plant Soil)
parameters_p[49,] <- c(3, 5) # n_age
# average needle length between years, mm
parameters_p[50,] <- c(34.241, 13) # n_lenght
# average length growth (with trend removed)
parameters_p[51,] <- c(309.0938, 309.0938)# h_increment
# m2 kg-1 ICOS 2018, half of total needle area / needle mass
parameters_p[52,] <- c(13, 5.5) # SLA

row.names(parameters_p) <- c("Q10.N", "Rm.NR", "Q10.S", "Rm.S", "Q10.R", "Rm.R", "sR0", "sRc", "growth.myco", "root.lifetime",
"HH0", "sH0", "LH", "LH0", "sHc", "sN0", "LN", "LN0", "sNc", "HN0", "sD0.Trad", "LD", "LD0", "sDc", "sDc.T.count",
"tau.Ee", "tau.El", "tau.We", "tau.Wl", "tau.GPP", "Uggla", "sB0", "sBc", "LB",
"cell.d.ew", "cell.d.lw", "cell.l.ew", "cell.l.lw", "cell.wall.density.ew", "cell.wall.density.lw",
"wall.thickness.ew", "wall.thickness.lw", "cell.volume.growth.per.day.ew", "cell.volumn.growth.per.day.lw",
"density_tree", "carbon_share", "D0", "h0", "n_age", "n_lenght", "h_increment", "SLA")
save(parameters_p, file = paste0(data.direct, "/parameters_p.RData"))

######
## Splerling parameters_p
######
sperling_p <- data.frame(matrix(ncol = 2, nrow = 30))
names(sperling_p) <- c("Hyde", "Lettosuo")
sperling_p[1,] <- c(0.3246781, NA) # starch0
sperling_p[2,] <- c(0.4184208, NA) # sugar0
sperling_p[3,] <- c(0.03, NA) # starch.needles0
sperling_p[4,] <- c(0.037, NA) # starch.phloem0
sperling_p[5,] <- c(0.034, NA) # starch.xylem.sh0
sperling_p[6,] <- c(0.166, NA) # starch.xylem.st0
sperling_p[7,] <- c(0.057, NA) # starch.roots0
sperling_p[8,] <- c(0.087, NA) # sugar.needles0
sperling_p[9,] <- c(0.027, NA) # sugar.phloem0
sperling_p[10,] <- c(0.014, NA) # sugar.roots0
sperling_p[11,] <- c(0.0249, NA) # sugar.xylem.sh0
sperling_p[12,] <- c(0.021, NA) # sugar.xylem.st0
sperling_p[13,] <- sperling_p[14,] <- sperling_p[15,] <- sperling_p[16,] <- sperling_p[17,] <- sperling_p[29,] <- c(0.0, 0.0) # Wala
sperling_p[18,] <- c(0.4211, 0.4211) # carbon.sugar
sperling_p[19,] <- c(0.4444, 0.4444) # carbon.starch
sperling_p[20,] <- sperling_p[21,] <- sperling_p[22,] <- sperling_p[23,] <- sperling_p[24,] <- sperling_p[30,] <- c(3, NA) # alfa
sperling_p[25,] <- c(2, NA) # tau.s
sperling_p[26,] <- c(2, NA) # tau.t
sperling_p[27,] <- c(0.3246781, NA) # starch00
sperling_p[28,] <- c(0.4184208, NA) # sugar00
sperling_p[31,] <- c(3, NA)   # Q10s As taken from a crop study for the tree study this model is based on I am assuming that these numbers can be generalised here as well
sperling_p[32,] <- c(1.8, NA)  # Q10d As taken from a crop study for the tree study this model is based on I am assuming that these numbers can be generalised here as well
sperling_p[33,] <- c(0.23, NA) # SCb kg C, just above lowest value in starch and sugar 15_16 -- JS notes
sperling_p[34,] <- c(0.41, NA) # sugar.level, 23/9/2015, Sugar, starch and sugar 15_16 -- JS notes kg C, "assumed to equal SC at senescence"
sperling_p[35,] <- c(0.017, NA) # Ad0.needles
sperling_p[36,] <- c(0.008, NA) # Ad0.phloem
sperling_p[37,] <- c(2e-04, NA) # Ad0.roots
sperling_p[38,] <- c(2e-04, NA) # Ad0.xylem.sh
sperling_p[39,] <- c(0.047, NA) # Ad0.xylem.st
sperling_p[40,] <- c(0.197, NA) # lamda.needles
sperling_p[41,] <- c(0.05301, NA) # lamda.phloem
sperling_p[42,] <- c(0.211, NA) # lamda.roots
sperling_p[43,] <- c(0.00401, NA) # lamda.xylem.sh
sperling_p[44,] <- c(0.00401, NA) # lamda.xylem.st
sperling_p[45,] <- c(0.729, NA) # delta.needles
sperling_p[46,] <- c(0.832, NA) # delta.phloem
sperling_p[47,] <- c(0.853, NA) # delta.roots
sperling_p[48,] <- c(0.762, NA) # delta.xylem.sh
sperling_p[49,] <- c(0.294, NA) # delta.xylem.st
sperling_p[50,] <- c(0.3, NA) # k_np
sperling_p[51,] <- c(0.072, NA) # k_pr
sperling_p[52,] <- c(0.188, NA) # k_pxsh
sperling_p[53,] <- c(0.17, NA) # k_pxst
sperling_p[54,] <- c(0.025, NA) # myco.thresh

row.names(sperling_p) <- c("starch0", "sugar0", "starch.needles0", "starch.phloem0", "starch.xylem.sh0", "starch.xylem.st0", "starch.roots0",
                     "sugar.needles0", "sugar.phloem0", "sugar.roots0", "sugar.xylem.sh0", "sugar.xylem.st0",
                     "Wala.needles", "Wala.phloem", "Wala.xylem.sh", "Wala.xylem.st", "Wala.roots",
                     "carbon.sugar", "carbon.starch",
                     "alfa.needles", "alfa.phloem", "alfa.xylem.sh", "alfa.xylem.st", "alfa.roots",
                     "tau.s", "tau.t", "starch00", "sugar00", "Wala", "alfa", "Q10s", "Q10d", "SCb", "sugar.level",
                     "Ad0.needles", "Ad0.phloem", "Ad0.roots", "Ad0.xylem.sh", "Ad0.xylem.st",
                     "lamda.needles", "lamda.phloem", "lamda.roots", "lamda.xylem.sh", "lamda.xylem.st",
                     "delta.needles", "delta.phloem", "delta.roots", "delta.xylem.sh", "delta.xylem.st",
                     "k_np", "k_pr", "k_pxsh", "k_pxst", "myco.thresh")

save(sperling_p, file = paste0(data.direct, "/sperling_p.RData"))

#######
## common_p parameters_p
#######

common_p <- data.frame(matrix(ncol = 16, nrow = 1))
# parameter of function g
common_p[1,1] <- 0.185 # a
# parameter of function g
common_p[1,2] <- 18.4	# b
# root temperature factor
common_p[1,3] <- 0 # TR0
# abs_zero
common_p[1,4] <- 273.15

# TODO: can't find these in the original script
# Soil water potential and conductance
# Duursma et al 2008 Table 2: A-hirzon (Note. wrong units), Right ones: H?ltt? et al. 2009 Table 2
common_p[1,5]<-4.14	# b.s
# H?ltt? et al. 2009 Table 2
common_p[1,6]<-0.62	#theetta.FC
# MPa H?ltt? et al. 2009 Table 2, note psi.s should be psi.e
common_p[1,7]<--680/10^6 # phi.e
# mol m-1 s-1 MPa-1, H?ltt? et al 2009, Table 2
common_p[1,8]<-24.5	# K.sat
#index (m root m-2), H?ltt? et al 2009
common_p[1,9]<-5300	# R.length
#kg mol-1, H?ltt? et al 2009
common_p[1,10]<-0.018	# M.H2O
#m, H?ltt? et al 2009
common_p[1,11]<-4.25*10^-3	# r.cyl
#m, H?ltt? et al 2009
common_p[1,12]<-3.0*10^-3	# r.root

#  to prevent from dividing with zero
common_p[1,13] <- 0.00000000000001 # ypsilon

# Growth respiration, share of growth
# Needles
common_p[1,14]<-0.35 # Rg.N
# Wood
common_p[1,15]<-0.3	# Rg.S
# Fine roots
common_p[1,16]<-0.35 # Rg.R

common_p[1,17]<-8.314 # gas.const
common_p[1,18]<-12.01	# M.C g / mol
common_p[1,19]<-1.008 # M.H
common_p[1,20]<-16 # M.O
common_p[1,21]<-2000000 # osmotic.sugar.conc Osmotic sugar concentration 2 MPa H?ltt? et al. 2000

names(common_p) <- c("a", "b", "TR0", "abs_zero", "b.s", "theetta.FC", "phi.e", "K.sat", "R.length",
                     "M.H20", "r.cyl", "r.root", "ypsilon", "Rg.N", "Rg.S", "Rg.R", "gas.const", "M.C", "M.H",
                     "M.O", "osmotic.sugar.conc")
save(common_p, file = paste0(data.direct, "/common_p.RData"))

######
## Repola
######

repo_p <- data.frame(matrix(ncol = 5, nrow = 1))
repo_p[1,1]<--6.303 # b0.repo
repo_p[1,2]<-14.472 # b1.repo
repo_p[1,3]<--3.976 # b2.repo
repo_p[1,4]<-0.109 # uk.repo
repo_p[1,5]<-0.118 # eki.repo
names(repo_p) <- c("b0.repo", "b1.repo", "b2.repo", "uk.repo", "eki.repo")
save(repo_p, file = paste0(data.direct, "/repo_p.RData"))

######
## Leap year
######
leap_years <- c(1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2004, 2008, 2012, 2016)
save(leap_years, file = paste0(data.direct, "/leap_years.RData"))

######
## Weather data for Hyde data between 2010 and 2019
######

inputdirectory <- "C:/Users/Käyttäjä/OneDrive - University of Helsinki/Research Assistant/CASSIA/Input/"
files <- list.files(path = inputdirectory, pattern = "Hyde_enfact")
files <- paste0(inputdirectory, files)
Hyde_weather.list <- lapply(files, read.table, header=TRUE, sep="\t", col.names=c("T","P","TSA","TSB", "MB", "Rain"))
Hyde_weather <- do.call("rbind", Hyde_weather.list)
Hyde_weather$date <- seq(as.POSIXct(as.character("2010-01-01")), as.POSIXct(as.character("2019-12-31")), by = "day")
Hyde_weather$date <- substring(Hyde_weather$date, 1, 10)
Hyde_weather <- Hyde_weather[,c("date","T","P","TSA","TSB", "MB", "Rain")]
save(Hyde_weather, file = paste0(data.direct, "/Hyde_weather.RData"))
