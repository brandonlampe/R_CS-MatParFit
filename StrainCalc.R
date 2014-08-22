# SCRIPT USED TO CALCULATE STRAINS & TEST FUNCTIONS PRIOR TO IMPLEMENTATION
# MATERIAL PARAMETER FITS WILL BE DETERMINED USING "ShearFit.R" & "CreepFit.R"
#
# # ============= START ================
# KAP0.INP <- 2.09011272735447
# KAP1.INP <- 1.32951908075078
# DDT.INP  <- 0.900012234998625
# NK.INP   <- 2.75382483799085
# KAP2.INP <- 1

KAP0.INP <- 10.119
KAP1.INP <- 1.005
DDT.INP  <- 0.896
NK.INP   <- 1.331
KAP2.INP <- 1

FLOW.INP<- cbind(KAP0.INP, KAP1.INP, KAP2.INP,
                 DDT.INP, NK.INP) # values for fiting parameters

# ==== INPUT INITIAL VALUES FOR CREEP PARAMETERS ====
# ---- Creep parameters, results from Callahan fits - shear tests only ----
ETA0  <- 0.1029          # -
ETA1  <- 3.9387          # -
ETA2 	<- 1               # constant -
NF 		<- 3.5122          # -
AA1 	<- 0.3147          # -
PP 		<- 1.6332          # -
NSP 	<- 0.5576          # -
R1 		<- 1.041 * 10 ^ -6 # [K/(MPa-sec)]
R3 		<- 15.1281         # -
R4 		<- 0.1678          # -
QSR 	<- 1077.46         # [K]

# # ---- Creep parameters, Lampe fit ----
# ETA0  <- 0.00432927511331067  # -
# ETA1  <- 7.82850052666962     # -
# ETA2  <- 1                    # constant -
# NF 		<- 3.27379002003667     # -
# AA1 	<- 47.8455448083257     # -
# PP 		<- 27.8645951602493     # -
# NSP 	<- -22.9778664287705    # -
# R1 		<- 1.99701832841102e-05 # [K/(MPa-sec)]
# R3 		<- -42.2194105338298    # -
# R4 		<- 14.2521644139889     # -
# QSR 	<- -4337.39844125619    # [K]

CREEP.INP <- cbind(ETA0, ETA1, ETA2, NF, AA1, PP,
                   NSP, R1, R3, R4, QSR) # INITIAL CREEP PARAMETER VALUES

CPar <- CREEP.INP
FPar <- FLOW.INP

# ---- use subset of full data set for debugging ----
browser()
TestData <- DATA.INP[which(DATA.INP$ITEST == "SC1B"),] # SUBSET OF DATA FOR ANALYSIS

# ---- Flow Potential Parameters (5) *KAP2 HELD CONST. ----
KAP0	<- as.numeric(FPar[1])
KAP1 	<- as.numeric(FPar[2])
KAP2 	<- as.numeric(FPar[3])	# Constant = 1
DDT 	<- as.numeric(FPar[4])
NK 		<- as.numeric(FPar[5])

# ---- Creep Consolidation Parameters (11) *ETA2 HELD CONST
ETA0 	<- as.numeric(CPar[1])
ETA1 	<- as.numeric(CPar[2])
ETA2 	<- as.numeric(CPar[3])	# Constant = 1
NF 		<- as.numeric(CPar[4])
AA1 	<- as.numeric(CPar[5])
PP 		<- as.numeric(CPar[6])
NSP 	<- as.numeric(CPar[7])
R1 		<- as.numeric(CPar[8])
R3 		<- as.numeric(CPar[9])
R4 		<- as.numeric(CPar[10])
QSR 	<- as.numeric(CPar[11])

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
RHOIS <- 2160.0	# ASSUMED IN SITU SALT DENSITY
NTIME <- 10^6		# NORMALIZING TIME
DSP 	<- 0.64		# FRACTIONAL DENSITY OF RANDOM DENSE SPHERICAL PARTICLES

# ---- Values input into function (18)----
#ICASE	<- as.numeric(TestData[,1])	  # TEST TYPE (1:Hyd Cons, 2:Shear Cons, 3:compaction)
ITEST <- as.character(TestData[,2])	# TEST ID
TIME 	<- as.numeric(TestData[,3]) 	# TIME [SEC]
DT 		<- as.numeric(TestData[,4])	  # DELTA TIME [SEC]
#TF 		<- as.numeric(TestData[,5])	  # TOTAL TEST TIME [SEC]
TEMP 	<- as.numeric(TestData[,6])	  # TEMP [K]
AS 		<- as.numeric(TestData[,7])	  # AXIAL STRESS [MPA]
LS 		<- as.numeric(TestData[,8])	  # LATERAL STRESS [MPA]
#EVT 	<- as.numeric(TestData[,9])	  # TOTAL TRUE VOLUMETRIC STRAIN
EVC 	<- as.numeric(TestData[,10])	# CREEP TRUE VOLUMETRIC STRAIN
#EAT 	<- as.numeric(TestData[,11])	# TOTAL TRUE AXIAL STRAIN
EAC 	<- as.numeric(TestData[,12])	# CREEP TRUE AXIAL STRAIN
RHO 	<- as.numeric(TestData[,13])	# CURRENT DENSITY [KG/M3]
D 		<- as.numeric(TestData[,14])	# FRACTIONAL DENSITY
RHO0 	<- as.numeric(TestData[,15])	# DENSITY AT THE START OF CONSOLIDATION (<RHOI)
RHOI 	<- as.numeric(TestData[,16])	# DENSITY AT THE START OF CREEP
DD 		<- as.numeric(TestData[,17])	# AVERAGE GRAIN SIZE [MM]
W 		<- as.numeric(TestData[,18])	# WATER CONENT BY PERCENT WEIGHT

# ---- calculate variables ----
MS 	<- (2.0 * LS + AS) / 3 	# MEAN STRESS
DS 	<- LS - AS				      # STRESS DIFFERENCE
ELC	<- (EVC - EAC) / 2		  # CREEP TRUE LATERAL STRAIN
D0 	<- 1382.4 / RHOIS			  # EMPLACED FRACTIONAL DENSITY (0.64 FRAC DENSITY)
DI 	<- RHOI / RHOIS			    # INITIAL FRACTIONAL DENSITY

WT1 <- DT / NTIME	  # WEIGHTING FUNCTION FOR CREEP CONSOLIDATION PARAMETERS
WT 	<- 1					  # WEIGHTING FUNCTION FOR FLOW PARAMETERS
#DC 	<- DD					  # SET GRAIN SIZE FOR DCCS TESTS

Z1	<- EAC  # Predicted axial strain (initial values)
Z2	<- ELC  # Predicted lateral strain (initial values)
Z3	<- rep(0, size(ELC)[2]) # internal variable "xi" for the transient function (FU)
                            # integral of Eqn 2-27, (initial values)

# ==== define the differential equation ====
# ---- only calculate strain rates at TIME > 0 ----
# browser()
ERATE.OUT <- data.frame(ifelse(cbind(TIME > 0, TIME > 0, TIME > 0),
{
  VOL  	<- Z1 + 2*Z2			    # VOLUMETRIC STRAIN
  VOLT	<- VOL + log(DSP/DI)	# USED FOR INITIAL ESTIMATE OF VOLUMETRIC STRAIN
  #DEN		<- DI/exp(VOL)			# CURRENT FRACTIONAL DENSITY
  DEN   <- D                  # CURRENT FRACTIONAL DENSITY

  ifelse(D >= 1,
{
  MD <- 0  # if fractional density is 1, disclocation creep = 0
  SP <- 0},# if fractional density is 1, pressure solutioning = 0

{
  VAR <- ifelse(DEN <= DDT, DDT, DEN) # DEFINE DENSITY CEILING ISH

  # ---- Equivalent Stress ----
  OMEGAA 	<- ((1 - DEN) * NF / (1 - (1 - DEN)^(1/NF))^NF)^(2/(NF + 1))
  OMEGAK 	<- ((1 - VAR) * NK / (1 - (1 - VAR)^(1/NK))^NK)^(2/(NK + 1))
  ETA		<- ETA0 * OMEGAA^ETA1
  KAP		<- KAP0 * OMEGAK^KAP1
  TERMA	<- ((2 - DEN)/DEN)^((2 * NF)/(NF + 1))
  TERMK	<- ((2 - DEN)/DEN)^((2 * NK)/(NK + 1))

  # ---- Eqn. 2-3 (SAND97-2601) ----
  SEQF	<- sqrt(ETA * MS^2 + ETA2 * TERMA * DS^2)	# Equivalent stress measure for Disc. Creep and Press Sol'ing
  SEQ		<- sqrt(KAP * MS^2 + KAP2 * TERMK * DS^2)	# Equivalent stress measure for Flow Potential

  # ---- Eqn. 2-17 (SAND97-2601) ----
  ALPHA2	<- KAP * MS / 3
  BETA2	<- KAP2 * TERMK * DS

  # ---- Eqn. 2-20 divided by equivalent stress (for later calculation) ----
  F2A <- 	(ALPHA2 - BETA2)/SEQ
  F2L <-	(ALPHA2 + 0.5 * BETA2)/SEQ

  # ==== START: equivalent inelastic strain rate form for dislocation creep ====

  # ---- Steady State Strain Rate Calc ----
  ES1 <- A1 * (SEQF / MU)^N1 * exp(-Q1R/TEMP)	# Dislocation climb - Eqn. 2-30
  ES2 <- A2 * (SEQF / MU)^N2 * exp(-Q2R/TEMP)	# Undefined Mechanism - Eqn. 2-31

  # Slip - Eqn. 2-32 (SAND98-2601)
  H   <- SEQF - S0                              # HEAVISIDE FUNCTION
  ARG <- Q * (SEQF - S0) / MU
  ES3 <- ifelse(H > 0, 0.5 * (B1 * exp(-Q1R / TEMP) +
                                (B2 * exp(-Q2R / TEMP)) *
                                (exp(ARG) - exp(-ARG))),0)

  ESS = ES1 + ES2 + ES3 # Steady-state strain rate, Eqn. 2-29 (SAND97-2601)

  # ---- EVALUATE TRANSIENT FUNCTION, 3 branches: work hardening, equilibrium, recovery
  EFT  <- K0 * exp(C * TEMP) * (SEQF / MU) ^ M  # Transient Strain Limit, Eqn. 2-28
  BIGD <- ALPHA + BETA * log10(SEQF / MU)       # Work-Hardening parameter, Eqn 2-28
  # ===================================
  # ----variables initialized outside of loop for debugging ----
  Z3  <- rep(0, size(ELC)[2])
  DT.I <- DT
  POW <- rep(0, size(ELC)[2])
  ERROR.OUT <- rep(1, size(ELC)[2])
  ERROR.INN <- rep(1, size(ELC)[2])
  # FU.I calculated with vector of zeros (Z3)
  FU.I <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                     exp(-DELTA * (1 - Z3 / EFT) ^ 2)))
  DZ3.I <- (FU.I - 1) * ESS # derivative calculated

  # derivative integrated to yield orignal values (first iteration)
  Z3.I <- ifelse(TIME > 0, 0.5 * DT * (shift(DZ3.I, 1, "right") + DZ3.I), 0)

  Z3.CHECK <- Z3.I # for debugging

  IT.OUT <- 0
  IT.INN <- 0
  IMAX <- 100
  TOL <- 1e-9
  IN.TIME <- proc.time()
  ERROR.OUT <- rep(1, size(ELC)[2])

 while((rev(ERROR.OUT) > TOL) && (IT.OUT < IMAX)){
    IT.OUT <- IT.OUT + 1
    Z3 <- Z3.I
    FU.I <- ifelse(Z3 == EFT, {1},{ifelse(Z3 < EFT,{
                          exp(BIGD * (1 - Z3 / EFT) ^ 2)},{
                          exp(-DELTA * (1 - Z3 / EFT) ^ 2)})})
    DZ3.I <- (FU.I - 1) * ESS
    Z3.I <- ifelse(TIME > 0, 0.5 * DT * (shift(DZ3.I, 1, "right") + DZ3.I), 0)
    ERROR.OUT <- abs(Z3.I - Z3)}
# ---- SPLIT TIME ----
#     if(ERROR.OUT > TOL){
# #       while((rev(ERROR.INN) < TOL) && (IT.INN < IMAX)){
#         POW <- POW + 1
#         IT.INN <- IT.INN + 1
#         DT.H <- DT.I/(2^POW)
#         Z3 <- Z3.I
#         FU.I <- ifelse(Z3 == EFT, {1},{ifelse(Z3 < EFT,{
#           exp(BIGD * (1 - Z3 / EFT) ^ 2)},{
#             exp(-DELTA * (1 - Z3 / EFT) ^ 2)})})
#         Z3.I <- ifelse(TIME > 0, 0.5 * DT.H * (shift(DZ3.I, 1, "right") + DZ3.I), 0)
#         ERROR.INN <- abs(Z3.I - Z3)
#         print(c("INNER WITH", IT.INN, ERROR.INN))}
# }

# ---- END OUTER WITH
    print(c("OUTER WITH", IT.OUT, ERROR.OUT))}
#===============================================
  OUT.TIME <- proc.time()
  Z3.I - Z3

abs(sum(Z3.I - Z3))

  FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                    exp(-DELTA * (1 - Z3 / EFT) ^ 2)))
  DZ3 <- (FU - 1) * ESS
  # ==================================
    MD <- FU * ESS  # equivalent inelastic strain rate form for dislocation creep, Eqn 2-23

  # ==== START: Equivalent Inelastic Strain Rate Form for Pressure Solutioning ====

  # ---- Calculate initial volumetric strain - Based on spherical packing ----
  CR <- abs(exp(VOLT) - 1)

  # ---- Determine functional form - either large or small strains, Eqn 2-34 ----
  GAMMA <- ifelse(CR <= 0.15, 1, abs((D0 - exp(VOLT)) / ((1 - D0) * exp(VOLT))) ^ NSP)
  # Small Strains (Vol Strain > - 15%)
  # Large Strains (Vol Strain < - 15%)

  # ---- component of eqn 2-35 ---
  X3 <- exp((R3 - 1) * VOLT) / (abs(1 - exp(VOLT))) ^ R4

  # ---- determine value of moisture function (w) ----
  M2 <- ifelse (W == 0, 0, W ^ AA1)     # moisture content  = 0
  # moisture content > 0

  G2 <- 1 / DD ^ PP # calculate grain size function
  T2 <- exp(-QSR / TEMP) / TEMP

  # ---- Equivalent Inelastic Strain Rate Form for Pressure Solutioning, Eqn 2-35
  SP <- R1 * M2 * G2 * T2 * X3 * GAMMA * SEQF})  # end check for D < 1

DZ1 <- (MD + SP) * F2A # Predicted axial strain rate / derivative of strain
DZ2 <- (MD + SP) * F2L # Predicted lateral strain rate / derivative of strain
DZ3 <- (FU - 1) * ESS  # Predicted Steady-State Creep Rate

c(DZ1, DZ2, DZ3)},{c(0,0,0)}))

colnames(ERATE.OUT) <- c("FEAR", "FELR", "FEVR")         # column names

DATA.FIT <- cbind(TestData, ERATE.OUT, WT1)       # merge data

DT.DATA.FIT <- data.table(DATA.FIT)       # change format to data table
setkey(DT.DATA.FIT, ITEST)                # create a value for subsets

# ---- data table ALL DATA (INCLUDING integrated strain rates) ----
DT.FE <- DT.DATA.FIT[, c("IFEAR", "IFELR", "IFEVR"):=list(
  as.vector(cumtrapz(TIME, FEAR)),
  as.vector(cumtrapz(TIME, FELR)),
  as.vector(cumtrapz(TIME, FEVR))), by = ITEST]

# ---- check lambda calculation ----
DT.RR <- data.table(ifelse(cbind(DT.FE$TIME > 0, DT.FE$TIME > 0),{
  RR1 <- DT.FE$EAC / DT.FE$IFEAR # ratio of measured to predicted axial strains
  RR2 <- DT.FE$ELC / DT.FE$IFELR # ratio of measured to predicted lateral strains
  c(RR1, RR2)},{
  c(0, 0)}))

lambda <- 1 - (((1 - DT.RR$RR1) ^ 2 + (1 - DT.RR$RR2) ^ 2 ) * WT1) ^ (1/2)

# ---- plot fit comparison (axial strain rate)----
library(ggplot2)
ggSUB.EAR <- ggplot(data = DT.FE, aes(x=TIME, y=EAR))
ggSUB.EAR <- ggSUB.EAR + geom_line()
ggSUB.EAR <- ggSUB.EAR + geom_point(aes(y=FEAR))
# ggSUB.EAR <- ggSUB.EAR + facet_wrap(~ITEST, ncol=3, scales = "free")
ggSUB.EAR <- ggSUB.EAR + xlim(0,6e6) + ylim(-7.5e-6,0)
ggSUB.EAR <- ggSUB.EAR + ylab("Axial Strain Rate: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggSUB.EAR

# ---- plot fit comparison (axial strain )----
ggSUB.EA <- ggplot(data = DT.FE, aes(x=TIME, y=EAC))
ggSUB.EA <- ggSUB.EA + geom_line()
ggSUB.EA <- ggSUB.EA + geom_point(aes(y=IFEAR))
# ggSUB.EA <- ggSUB.EA + facet_wrap(~ITEST, ncol=3, scales = "free")
ggSUB.EA <- ggSUB.EA + xlim(0,6e6) + ylim(-0.25,0)
ggSUB.EA <- ggSUB.EA + ylab("Axial Strain: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggSUB.EA

