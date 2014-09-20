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
  #   NF   	<- 3.5122

  # ---- EFFECTIVE STRESS PARAMETERS FOR CRUSHED SALT FORMULATION OF MODEL ----
  #   # ---- VALUES FOR MATHCAD CHECK
  KAP0 <- 10.119
  KAP1 <- 1.005
  NK   <- 1.331
  #--------------
  DDT  <- 0.8854
  KAP2 <- 1
  ETA0  <- 0.15
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
  AS   <- (-6/300.9)*Time
  LS   <- (-2/300.9)*Time
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
  EFT  <- K0 * exp(C * TEMP) * (SEQF / MU) ^ M  # Transient Strain Limit, Eqn. 2-28
  BIGD <- ALPHA + BETA * log10(SEQF / MU)       # Work-Hardening parameter, Eqn 2-28

  FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                    exp(-DELTA * (1 - Z3 / EFT) ^ 2)))

  MD <- FU * ESS  # equivalent inelastic strain rate form for dislocation creep, Eqn 2-23

  DZ3 <- (FU - 1) * ESS  # derivative of internal variable "ZETA"
  DZ <- list(c(DZ3), MD, FU, ESS, ES1, ES2, ES3,
             DZ3, EFT, SEQF, BIGD, AS, LS)

  return(DZ)
}

ZETA.EFT <- ODE.DT$Z3 - ODE.DT$EFT
ODE.DT2 <- ODE.DT[1:28,]

# ========================================================================
# ---- compare to callahan results ---
Z3.CALLAHAN <-c(1e-12, 1.836e-7, 2.111e-7, 3.886e-7, 1.037e-6, 2.585e-6,
                6.504e-6, 1.628e-5, 3.057e-5)

NAME.MCD <- rbind(cbind(rep("Callahan", 9)),cbind(rep("BCL", 9)))
TIME.MCD <- cbind(rep(ODE.DT$TIME, 2))
ZETA.MCD <- rbind(cbind(Z3.CALLAHAN), cbind(ODE.DT$Z3))
RATZ.MCD <- rbind(cbind(rep(ODE.DT$Z3 /Z3.CALLAHAN,2)))
MCD.DT <- data.table(NAME.MCD, TIME.MCD, ZETA.MCD, RATZ.MCD)

setnames(MCD.DT, c("NAME", "TIME", "ZETA", "BCLtoCAL"))

P.MCD <- ggplot(MCD.DT, aes(x = TIME, y = ZETA, color = NAME))
P.MCD <- P.MCD + geom_point() + geom_line()
P.MCD

P.RATZ <- ggplot(MCD.DT, aes(x = TIME, y = BCLtoCAL))
P.RATZ <- P.RATZ + geom_point() + geom_line()
P.RATZ

# # ---- EXPORT AND SAVE FILES AS PDF
PATH = "/Users/Lampe/GrantNo456417/MatParameterFitting/R_CS-MatParFit"
FILE.NAME = "/Zeta_MCD_Files.pdf"
pdf(file = paste(PATH, FILE.NAME, sep = ""), onefile = TRUE)
P.MCD
P.RATZ
dev.off()