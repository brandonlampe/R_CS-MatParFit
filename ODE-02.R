#############################
# solves differential equation for creep consolidation with presecribed stress
# runs an array of tests at a time
#############################

source("DoAll.R")

# ==== DIFFERENTIAL EQUATION ====
STRAINS.02 <- function(Time, State, Parm){
  #   with(as.list(c(State, Parm)),{
  # function for calculating axial and lateral strain rates
  # Input must be in vector or matrix form, no data frames
  # Eqns. referenced from:  SAND97-2601
  # CPar: EAT0, ETA1, ETA2, NF, AA1, PP, NSP, R1, R3, R4, QSR
  # FPar: KAP0, KAP1, KAP2, NK, DDT
  # TestData:


  # ============= parameters hard coded into function directly ========
  KAP0 <- 10.119
  KAP1 <- 1.005
  DDT  <- 0.896
  NK   <- 1.331
  KAP2 <- 1

  ETA0  <- 0.102854        # -
  ETA1  <- 3.9387          # -
  ETA2   <- 1               # constant -
  NF 		<- 3.5122          # -
  AA1 	<- 0.3147          # -
  PP 		<- 1.6332          # -
  NSP 	<- 0.557621        # -
  R1 		<- 0.0194          # [K/(MPa-sec)]
  R3 		<- 15.1281         # -
  R4 		<- 0.1677765       # -
  QSR 	<- 1077.46         # [K]

  # ---- Munson-Dawson Creep Parameters (17) ----
  A1 		<- 8.386e22
  A2 		<- 9.672e12
  Q1R 	<- 12581
  Q2R 	<- 5033
  N1 		<- 5.5
  N2 		<- 5.0
  B1 		<- 6.0856e6
  B2 		<- 3.034e-2
  Q 		<- 5335
  S0 		<- 20.57
  M 		<- 3
  K0 		<- 6.275e5
  C 		<- 9.198e-3
  ALPHA <- -17.37
  BETA 	<- -7.738
  DELTA <- 0.58
  MU 		<- 12400

  # ---- fitting assumptions ----
  RHOIS <- 2160.0  # ASSUMED IN SITU SALT DENSITY

  # ---- interpolated input variables ----
  TIME <- time.interp(Time)
  TEMP <- temp.interp(Time)
  AS   <- as.interp(Time)
  LS   <- ls.interp(Time)

  # ---- calculate variables ----
  MS 	<- (2.0 * LS + AS) / 3 	# MEAN STRESS
  DS 	<- LS - AS				      # STRESS DIFFERENCE
  D0 	<- 1382.4 / RHOIS			  # EMPLACED FRACTIONAL DENSITY (0.64 FRAC DENSITY)
  DI 	<- RHOI / RHOIS			    # FRACTIONAL DENSITY at the start of creep

  # ==============================================================
  Z1	<- State[1]  # Predicted axial strain (initial values)
  Z2	<- State[2] # Predicted lateral strain (initial values)
  Z3	<- State[3] # internal variable "zeta" for the transient function (FU)

  RHOI <- Parm[1]  # DENSITY AT THE START OF CREEP
  DD 	 <- Parm[2]  # AVERAGE GRAIN SIZE [MM]
  W    <- Parm[3]

  # integral of Eqn 2-27, (initial values)

  # ==== define the differential equation ====
  VOL  	<- Z1 + 2*Z2			    # VOLUMETRIC STRAIN
  VOLT	<- VOL + log(D0/DI)	  # VOLUMETRIC STRAIN + INITIAL TRUE STRAIN ESTIMATE
  DEN   <- DI/exp(VOL)        # CURRENT FRACTIONAL DENSITY

  ifelse(DEN >= 1,{
    MD <- 0  # if fractional density is 1, disclocation creep = 0
    SP <- 0},# if fractional density is 1, pressure solutioning = 0
{
  VAR <- ifelse(DEN <= DDT, DDT, DEN) # DEFINE DENSITY floor ISH
  # ==== DEBUG ====

  # ---- Equivalent Stress ----
  OMEGAA 	<- ((1 - DEN) * NF / (1 - (1 - DEN)^(1/NF))^NF)^(2/(NF + 1))
  OMEGAK 	<- ((1 - VAR) * NK / (1 - (1 - VAR)^(1/NK))^NK)^(2/(NK + 1))
  ETA		<- ETA0 * OMEGAA^ETA1
  KAP		<- KAP0 * OMEGAK^KAP1
  TERMA	<- ((2 - DEN)/DEN)^((2 * NF)/(NF + 1))
  TERMK	<- ((2 - DEN)/DEN)^((2 * NK)/(NK + 1))

  # ---- Eqn. 2-3 (SAND97-2601) ----
  # Equivalent stress measure for Disl. Creep and Press Sol'ing
  SEQF	<- sqrt(ETA * MS^2 + ETA2 * TERMA * DS^2)
  # Equivalent stress measure for Flow Potential
  SEQ		<- sqrt(KAP * MS^2 + KAP2 * TERMK * DS^2)

  # ---- Eqn. 2-17 (SAND97-2601) ----
  ALPHA2	<- KAP * MS / 3
  BETA2	<- KAP2 * TERMK * DS

  # ---- Eqn. 2-20, WithOUT dislocation creep and pressure solutioning ----
  F2A <- 	(ALPHA2 - BETA2)/SEQ        # fit to axial strains
  F2L <-	(ALPHA2 + 0.5 * BETA2)/SEQ  # fit to lateral strains
  F2V <-  3 * ALPHA2 / SEQ            # fit to volumetric strains

  # ==== START: equivalent inelastic strain rate form for dislocation creep ====
  # ---- Steady State Strain Rate Calc ----
  ES1 <- A1 * (SEQF / MU)^N1 * exp(-Q1R/TEMP)	# Dislocation climb - Eqn. 2-30
  ES2 <- A2 * (SEQF / MU)^N2 * exp(-Q2R/TEMP)	# Undefined Mechanism - Eqn. 2-31

  # Slip - Eqn. 2-32 (SAND98-2601)
  H   <- SEQF - S0 # HEAVISIDE FUNCTION
  ARG <- Q * (SEQF - S0) / MU
  ES3 <- ifelse(H >= 0, 0.5 * (B1 * exp(-Q1R / TEMP) + B2 * exp(-Q2R / TEMP)) *
                  (exp(ARG) - exp(-ARG)),0)

  ESS = ES1 + ES2 + ES3 # Steady-state strain rate, Eqn. 2-29 (SAND97-2601)

  # ---- EVALUATE TRANSIENT FUNCTION, 3 branches: work hardening, equilibrium, recovery
  EFT  <- K0 * exp(C * TEMP) * (SEQF / MU) ^ M  # Transient Strain Limit, Eqn. 2-28
  BIGD <- ALPHA + BETA * log10(SEQF / MU)       # Work-Hardening parameter, Eqn 2-28

  FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                    exp(-DELTA * (1 - Z3 / EFT) ^ 2)))

  MD <- FU * ESS  # equivalent inelastic strain rate form for dislocation creep, Eqn 2-23

  # ==== START: Equivalent Inelastic Strain Rate Form for Pressure Solutioning ====
  # ---- Calculate initial volumetric strain - Based on spherical packing ----
  CR <- abs(exp(VOLT) - 1) # USES THE DEFINITION OF ENGINEERING STRAIN

  # ---- Determine functional form - either large or small strains, Eqn 2-34 ----
  GAMMA <- ifelse(CR <= 0.15, 1,
                  (abs((D0 - exp(VOLT)) / ((1 - D0) * exp(VOLT)))) ^ NSP)
  # Small Strains (Vol Strain > - 15%)
  # Large Strains (Vol Strain < - 15%)

  # ---- component of eqn 2-35 ---
  X3 <- exp((R3 - 1) * VOLT) / (abs(1 - exp(VOLT))) ^ R4

  # ---- determine value of moisture function (w) ----
  M2 <- ifelse (W == 0, 0, W ^ AA1)

  # ---- Equivalent Inelastic Strain Rate Form for Pressure Solutioning, Eqn 2-35
  G2 <- 1 / DD ^ PP # calculate grain size function
  T2 <- exp(-QSR / TEMP) / TEMP
  SP <- R1 * M2 * G2 * T2 * X3 * GAMMA * SEQF })

DZ1 <- (MD + SP) * F2A # derivative: axial strain rate
DZ2 <- (MD + SP) * F2L # derivative: lateral strain rate
DZ3 <- (FU - 1) * ESS  # derivative of internal variable "ZETA"

DZ <- list(c(DZ1, DZ2, DZ3), MD, FU, ESS, ES1, ES2,
           ES3, SP, DZ1, DZ2, DZ3, F2A, F2L, EFT, SEQ, SEQF, BIGD,
           AS, LS, DEN, GAMMA)
return(DZ)
}

# ==== RESETS VALUES OF Z1, Z2, Z3 TO ZERO AT THE START OF EACH TEST ====
# ---- SETS STRAINS AND RATES = 0 AT TIME == 0 ----
reset <- function(Time, State){
  with (as.list(State),{

    # ---- if time == 0, derivative = 0 ----
    if(Time == 0){
      Z1  <- 1e-12 # initial axial strain (initial values)
      Z2	<- 1e-12 # initial lateral strain (initial values)
      Z3	<- 1e-12 # initial internal variable "zeta"
    }
    return(c(Z1, Z2, Z3))
  })
}

# ==== sepcify input parameters for ODE ====

PAR.TEST <- DATA.INP
# PAR.TEST <- PAR.TEST[1:100,]
# debug.out <- paste(CurrentDirectory,"debug_SC1B.csv",sep = "/")
# write.table(ODE.DT,file = debug.out, sep = ",")

# ---- linear interpolation functions to be called in "strain_Rates.01" ----
time.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TIME)
temp.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TEMP)
as.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$AS)
ls.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$LS)

RHOI   <- as.numeric(PAR.TEST$RHOI[1])	# DENSITY AT THE START OF CREEP
DD 		<- as.numeric(PAR.TEST$DD[1])	# AVERAGE GRAIN SIZE [MM]
W 		<- as.numeric(PAR.TEST$W[1])	# WATER CONENT BY PERCENT WEIGHT

PARM <- c(RHOI = RHOI, DD = DD, W = W) # CONSTANT TEST SPECIFIC PARAMETERS

# ---- intial values for state variables ----
Z1	<- 1e-12 # Predicted axial strain (initial values)
Z2	<- 1e-12 # Predicted lateral strain (initial values)
Z3	<- 1e-12 # internal variable "xi" for the transient function (FU)
# integral of Eqn 2-27, (initial values)

IC <- (c(Z1 = Z1, Z2 = Z2, Z3 = Z3)) # array of initial values

TIME <- PAR.TEST$TIME

# ---- function for Predicting the Creep Strain(E) Rates ----
P.CER <- dede(func = STRAINS.02, parms = PARM, y = IC,
              times = TIME, verbose = TRUE, events = reset)

# ==== DEBUG ====
ODE.DT <- data.table(P.CER)
setnames(ODE.DT,c("TIME", "Z1", "Z2", "Z3", "MD", "FU", "ESS",
                  "ES1", "ES2", "ES3", "SP", "DZ1",
                  "DZ2", "DZ3", "F2A", "F2L", "EFT", "SEQ",
                  "SEQF", "BIGD", "AS", "LS","DEN","GAMMA"))

# ==== GATHER MEASURED VS PREDICTED DATA AND PLOT ====
# ---- MEASURED DATA ----
M.DF <- data.frame(cbind(PAR.TEST$TIME ,PAR.TEST$AS ,PAR.TEST$LS
                         ,PAR.TEST$ELC ,PAR.TEST$EAC ,PAR.TEST$D))
# ---- PREDICTED DATA ----
P.DF <- data.frame(cbind(ODE.DT$TIME, ODE.DT$AS , ODE.DT$LS, ODE.DT$Z2
                         , ODE.DT$Z1, ODE.DT$DEN))
MVSP.DF <- rbind(M.DF, P.DF)
MVSP.DF$ITEST <- rep(TestName,nrow(MVSP.DF))
MVSP.DF$TYPE <- c(rep("Meas",nrow(M.DF)), rep("Pred",nrow(P.DF)))
colnames(MVSP.DF) <- c("TIME", "AS", "LS", "ELC", "EAC", "D", "ITEST", "TYPE")

EAC.MVSP <- ggplot(MVSP.DF, aes(x = TIME, color = TYPE))
EAC.MVSP <- EAC.MVSP + geom_point(aes(y = EAC, shape = TYPE))
EAC.MVSP <- EAC.MVSP + ylab("AXIAL CREEP STRAIN")
EAC.MVSP

ELC.MVSP <- ggplot(MVSP.DF, aes(x = TIME, color = TYPE))
ELC.MVSP <- ELC.MVSP + geom_point(aes(y = ELC, shape = TYPE))
ELC.MVSP <- ELC.MVSP + ylab("LATERAL CREEP STRAIN")
ELC.MVSP

DEN.MVSP <- ggplot(MVSP.DF, aes(x=TIME, color = TYPE))
DEN.MVSP <- DEN.MVSP + geom_point(aes(y = D))
DEN.MVSP <- DEN.MVSP + ylab("FRACTIONAL DENSITY")
DEN.MVSP

PATH = "/Users/Lampe/GrantNo456417/MatParameterFitting/R_CS-MatParFit"
FILE.NAME = "/CS_LAMPE-MEAS.pdf"
pdf(file = paste(PATH, FILE.NAME, sep = ""), onefile = TRUE)

EAC.MVSP
ELC.MVSP
DEN.MVSP

dev.off()
