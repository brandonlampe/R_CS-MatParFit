##################################
# FILE SET UP TO CHECK ZETA CALCULATION
# PARAMETERS ARE FOR COMPARISON OF CALLAHAN RESULTS IN MATHCAD FILE
##################################

# ==== DIFFERENTIAL EQUATION ====
ZETA.01 <- function(Time, State, Parm){

#   # ---- EFFECTIVE STRESS PARAMETERS FOR CRUSHED SALT FORMULATION OF MODEL ----
#   KAP0 <- 10.119
#   KAP1 <- 1.005
#   DDT  <- 0.896
#   NK   <- 1.331
#   KAP2 <- 1
#   ETA0  <- 0.102854
#   ETA1  <- 3.9387
#   ETA2   <- 1
#   NF 		<- 3.5122

  # ---- EFFECTIVE STRESS PARAMETERS FOR CRUSHED SALT FORMULATION OF MODEL ----
#   # ---- VALUES FOR MATHCAD CHECK
  KAP0 <- 10.119
  KAP1 <- 1.005
  NK   <- 1.331
  #--------------
  DDT  <- 0.8854
  KAP2 <- 1
  ETA0  <- 0.015
  ETA1  <- 1
  ETA2   <- 0.7
  NF   	<- 8.40075

#   # ---- Munson-Dawson Creep Parameters (17) ---- FOR CLEAN SALT
#   A1   	<- 8.386e22
#   A2 		<- 9.672e12
#   Q1R 	<- 12581
#   Q2R 	<- 5033
#   N1 		<- 5.5
#   N2 		<- 5.0
#   B1 		<- 6.0856e6
#   B2 		<- 3.034e-2
#   Q 		<- 5335
#   S0 		<- 20.57
#   M 		<- 3
#   K0 		<- 6.275e5
#   C 		<- 9.198e-3
#   ALPHA <- -17.37
#   BETA 	<- -7.738
#   DELTA <- 0.58
#   MU 		<- 12400

  # ---- Munson-Dawson Creep Parameters (17) ---- FOR ARGILLACEOUS SALT
  A1    <- 1.407e23
  A2 		<- 1.314e13
  Q1R 	<- 12581
  Q2R 	<- 5033
  N1 		<- 5.5
  N2 		<- 5.0
  B1 		<- 8.988e6
  B2 		<- 4.289e-2
  Q 		<- 5335
  S0 		<- 20.57
  M 		<- 3
  K0 		<- 2.47e6
  C 		<- 9.198e-3
  ALPHA <- -14.96
  BETA 	<- -7.738
  DELTA <- 0.58
  MU 		<- 12400

  #---- DATA INTERPRETED FROM TEST DATA
#   TEMP <- temp.interp(Time)
#   AS   <- as.interp(Time)
#   LS   <- ls.interp(Time)
#   D    <- d.interp(Time)

  #---- constant data values for comparison to Callahan's analysis
  TEMP <- 300
  AS   <- (-6/300.8)*Time
  LS   <- (-2/300.8)*Time
  D    <- 0.9

  # ---- calculate variables ----
# browser()
  MS  <- (2.0 * LS + AS) / 3   # MEAN STRESS
  DS 	<- LS - AS				       # STRESS DIFFERENCE
  DEN <- D                     # CURRENT FRACTIONAL DENSITY

  Z3  <- State[1] # internal variable "zeta" for the transient function (FU)

  VAR <- ifelse(DEN <= DDT, DDT, DEN)

  # ---- Equivalent Stress ----
  OMEGAA   <- ((1 - DEN) * NF / (1 - (1 - DEN)^(1/NF))^NF)^(2/(NF + 1))
  ETA		   <- ETA0 * OMEGAA^ETA1
  TERMA	   <- ((2 - DEN)/DEN)^((2 * NF)/(NF + 1))

  # ---- Eqn. 2-3 (SAND97-2601) ----
  # Equivalent stress measure for Disl. Creep
  SEQF	<- sqrt(ETA * MS^2 + ETA2 * TERMA * DS^2)

  # ==== START: equivalent inelastic strain rate form for dislocation creep ====
  # ---- Steady State Strain Rate Calc ----
  ES1 <- A1 * (SEQF / MU)^N1 * exp(-Q1R/TEMP)	# Dislocation climb - Eqn. 2-30
  ES2 <- A2 * (SEQF / MU)^N2 * exp(-Q2R/TEMP)	# Undefined Mechanism - Eqn. 2-31

  # Slip - Eqn. 2-32 (SAND98-2601)
#   browser()
  ARG <- Q * ((SEQF - S0) / MU)
  ES3 <- (B1 * exp(-Q1R / TEMP) + B2 * exp(-Q2R / TEMP)) * sinh(ARG) * Heaviside(SEQF - S0)

  ESS = ES1 + ES2 + ES3 # Steady-state strain rate, Eqn. 2-29 (SAND97-2601)

  # ---- EVALUATE TRANSIENT FUNCTION, 3 branches: work hardening, equilibrium, recovery
#   browser()
  EFT  <- K0 * exp(C * TEMP) * (SEQF / MU) ^ M  # Transient Strain Limit, Eqn. 2-28
  BIGD <- ALPHA + BETA * log10(SEQF / MU)       # Work-Hardening parameter, Eqn 2-28

  FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                      exp(-DELTA * (1 - Z3 / EFT) ^ 2)))

  MD <- FU * ESS  # equivalent inelastic strain rate form for dislocation creep, Eqn 2-23

  DZ3 <- (FU - 1) * ESS  # derivative of internal variable "ZETA"
  DZ <- list(c(DZ3), MD, FU, ESS, ES1, ES2, ES3,
             DZ3, EFT, SEQF, BIGD, AS, LS, ETA, OMEGAA)

  return(DZ)
}

# =================================================================================
PAR.TEST <- DATA.INP[which(DATA.INP$ITEST == "SC1B"),] # SUBSET OF DATA FOR ANALYSIS
# PAR.TEST <- PAR.TEST[1:28,]

# ---- linear interpolation functions to be called in "ZETA.01" ----
temp.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TEMP)
as.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$AS)
ls.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$LS)
d.interp    <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$D)

# ---- intial values for state variables ----
Z3	<- 1e-12 # internal variable "ZETA" for the transient function (FU)
# integral of Eqn 2-27, (initial values)

IC <- (c(Z3 = Z3)) # array of initial values

# TIME <- PAR.TEST$TIME
TIME <- c(0, 50, 52.5, 65, 91.61, 126.101, 174.21, 240.386, 300) #for comparison

# ---- function for Predicting the Creep Strain(E) Rates ----
 ZETA.ODE <- ode(func =ZETA.01, y = IC, times = TIME, verbose = TRUE,
                 maxsteps = 1000000, rtol = 1e-14, atol = 1e-14, hmax = 1 )

# ==== ANALYSIS OF RESULTS=============================
ODE.DT <- data.table(ZETA.ODE)
setnames(ODE.DT,c("TIME","Z3",
                  "MD", "FU", "ESS",
                  "ES1", "ES2", "ES3", "DZ3", "EFT",
                  "SEQF", "BIGD", "AS", "LS", "ETA", "OMEGAA"))

ZETA.EFT <- ODE.DT$Z3 - ODE.DT$EFT
ODE.DT2 <- ODE.DT[1:28,]

# ========================================================================
# # ---- compare to callahan results ---

ES1.CALLAHAN <- as.numeric(c(0, 0, 0, 0, 0, 1.298e-15, 1.05e-14, 8.451e-14, 3.553e-13))
ES2.CALLAHAN <- as.numeric(c(0, 3.766e-15, 5.206e-15, 1.783e-14, 1.371e-13, 9.212e-13,
                 6.346e-12, 4.351e-11, 1.638e-10))
ES3.CALLAHAN <- as.numeric(c(0, 0, 0, 0, 0, 0, 0, 0, 0))
ESS.CALLAHAN <- as.numeric(ES1.CALLAHAN + ES2.CALLAHAN + ES3.CALLAHAN)
Z3.CALLAHAN <- as.numeric(c(1e-12, 1.836e-7, 2.111e-7, 3.886e-7, 1.037e-6, 2.585e-6,
                6.504e-6, 1.628e-5, 3.057e-5))
LS.CALLAHAN <- as.numeric(-c(0, 0.323, 0.34, 0.423, 0.601, 0.831, 1.151, 1.593, 1.99))
AS.CALLAHAN <- as.numeric(-c(0, 0.99, 1.04, 1.29, 1.822, 2.512, 3.474, 4.798, 5.99))
SEQF.CALLAHAN <- as.numeric(c(0, 0.719, 0.754, 0.932, 1.309, 1.798, 2.48, 3.418, 4.264))

RATZ3 <- as.numeric(ODE.DT$Z3 / Z3.CALLAHAN)

CALLAHAN <- data.table(as.numeric(ODE.DT$TIME),
                      as.numeric(Z3.CALLAHAN), as.numeric(LS.CALLAHAN),
                      as.numeric(AS.CALLAHAN), as.numeric(SEQF.CALLAHAN),
                      as.numeric(RATZ3))

BCL <- data.table(as.numeric(ODE.DT$TIME), as.numeric(ODE.DT$Z3),
                  as.numeric(ODE.DT$LS), as.numeric(ODE.DT$AS),
                  as.numeric(ODE.DT$SEQF), as.numeric(RATZ3))

NAME <- rbind(cbind(rep("CALLAHAN", 9)),cbind(rep("BCL", 9)))

MCD.DT <- data.table(NAME, rbind(CALLAHAN, BCL))

setnames(MCD.DT, c("NAME", "TIME","ZETA", "LS",
                   "AS", "SEQF", "RATZ3"))

P.Z <- ggplot(MCD.DT, aes(x = TIME, y = ZETA, color = NAME))
P.Z <- P.Z + geom_line() + geom_point()
P.Z

P.SEQF <- ggplot(MCD.DT, aes(x = TIME, y = SEQF, color = NAME))
P.SEQF <- P.SEQF + geom_line() + geom_point()
P.SEQF


# # ---- EXPORT AND SAVE FILES AS PDF
PATH = "/Users/Lampe/GrantNo456417/MatParameterFitting/R_CS-MatParFit"
FILE.NAME = "/Zeta_MCD_Files_V02.pdf"
pdf(file = paste(PATH, FILE.NAME, sep = ""), onefile = TRUE)
P.Z
P.SEQF
dev.off()
# #=======================================================================
# # ----  PLOT DATA ---------------------------
# P.Z <- ggplot(ODE.DT, aes(x = TIME))
# P.Z <- P.Z + geom_point(aes(y = Z3, color = "ZETA"))
# P.Z <- P.Z + geom_point(aes(x = TIME, y = EFT, color = "TRANS STRAIN LIMIT (EFT)"))
# P.Z
#
# P.Z2 <- ggplot(ODE.DT2, aes(x = TIME))
# P.Z2 <- P.Z2 + geom_point(aes(y = Z3, color = "ZETA"))
# P.Z2 <- P.Z2 + geom_point(aes(y = EFT, color = "TRANS STRAIN LIMIT (EFT)"))
# P.Z2
#
# P.ZE <- ggplot(ODE.DT, aes(x = TIME))
# P.ZE <- P.ZE + geom_point(aes(y = ZETA.EFT, color = "ZETA MINUS EFT"))
# P.ZE
#
# P.ZETA <- ggplot(ODE.DT, aes(x = TIME))
# P.ZETA <- P.ZETA + geom_point(aes(y = Z3, color = "ZETA"))
# P.ZETA
#
# P.ZETA2 <- ggplot(ODE.DT2, aes(x = TIME))
# P.ZETA2<- P.ZETA2 + geom_point(aes(y = Z3, color = "ZETA"))
# P.ZETA2
#
# P.SEQF <- ggplot(ODE.DT2, aes(x = TIME))
# P.SEQF <- P.SEQF + geom_point(aes(y = SEQF, color = "Equivalent Stress"))
# P.SEQF
#
# P.FU <- ggplot(ODE.DT2, aes(x = TIME))
# P.FU <- P.FU + geom_point(aes(y = FU, color = "Transient Function"))
# P.FU
#
# P.ESS <- ggplot(ODE.DT2, aes(x = TIME))
# P.ESS <- P.ESS + geom_point(aes(y = ESS, color = "Steady-State Strain Rate"))
# P.ESS
#
# P.S <- ggplot(ODE.DT2, aes(x = TIME, STRESS))
# P.S <- P.S + geom_point(aes(y = AS, color = "axial stress"))
# P.S <- P.S + geom_point(aes(y = LS, color = "lateral stress"))
# P.S
#
# # # ---- EXPORT AND SAVE FILES AS PDF
# PATH = "/Users/Lampe/GrantNo456417/MatParameterFitting/R_CS-MatParFit"
# FILE.NAME = "/ZetaAnalysis.pdf"
# pdf(file = paste(PATH, FILE.NAME, sep = ""), onefile = TRUE)
#   P.Z2
#   P.ZETA2
#   P.SEQF
#   P.FU
#   P.ESS
#   P.S
# dev.off()
#
