#==== FUNCTION FOR EVALUATING SHEAR PARAMETERS ONLY ====
shear_FIT <- function(FPar,TestData) {
  # this functions calculates the later/axial strain rate
  # from Callahan 1999, Equation 4-3 on page 31
  # VALID FOR AXIAL COMPRESSION
  # FPar: KAP0, KAP1, DDT, NK  
  # TestData: D, AS,LS  MS, DS
  
  #browser()
  
	# deconstruct matrices
	KAP0   <- as.numeric(FPar[1]) #Parameters$a #[1]
	KAP1   <- as.numeric(FPar[2]) #Parameters$b #[2]
	DDT    <- as.numeric(FPar[3]) #Parameters$c #[3]
	NK     <- as.numeric(FPar[4]) #Parameters$d #[4]

	AS     <- as.numeric(TestData[,1]) # Measured axial stress [MPa]
	LS     <- as.numeric(TestData[,2]) # Measured lateral stress [MPa]
	D      <- as.numeric(TestData[,3]) # Fractional density

	MS     <- (2*LS + AS)/3  # mean stress (MPa)
	DS     <- LS - AS 			 # stress difference (MPa)
  KAP2   <- 1	             # CONSTANT

	VAR    <- ifelse(D <= DDT,DDT,D)
	NUM    <- (1 - VAR) * NK
	DEN    <- (1 - (1 - VAR) ^ (1 / NK)) ^ NK
	OMEGAK <- (NUM / DEN) ^ (2 / (NK + 1))
	KAP    <- KAP0 * OMEGAK ^ KAP1
	TERMK  <- ((2 - D) / D) ^ (2 * NK / (NK + 1))

	ALPHA2 <- KAP * MS / 3
	BETA2  <- KAP2 * TERMK * DS

	F2A    <- ALPHA2 - BETA2
	F2L    <- ALPHA2 + 0.5 * BETA2
	FRAT   <- F2L/F2A # calculated ratio of lateral-to-axial strain rates
	return(FRAT)}

# ==== Calculate residuals in the shear parameter fit ====
shear_RESID <- function (p,Observed,xx) {
  # Response function that defines residuals FOR SHEAR PARAMETERS ONLY
  # called by nls.lm
	# p  = input values for parameters:  KAP0, KAP1, DDT, NK
  # xx = values from test data: D, MS, DS
  # Observed = measured RAT (ELR/EAR)
  
  #browser()
  Residual <- Observed - shear_FIT(p,xx)}

# ======================================================
strain_Rates <- function(CPar, FPar, TestData){
  # function for calculating axial and lateral strain rates
  # Input must be in vector or matrix form, no data frames
  # Eqns. referenced from:  SAND97-2601
  # CPar: EAT0, ETA1, ETA2, NF, AA1, PP, NSP, R1, R3, R4, QSR
  # FPar: KAP0, KAP1, KAP2, NK, DDT
  # TestData: 
  
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
	#ITEST <- as.character(TestData[,2])	# TEST ID
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
	D0 	<- 1382.4 / RHOIS			  # EMPLACED FRACTIONAL DENSITY ( NOT SURE WHERE 1382.4 CAME FROM?)
	DI 	<- RHOI / RHOIS			    # INITIAL FRACTIONAL DENSITY

	WT1 <- DT / NTIME	  # WEIGHTING FUNCTION FOR CREEP CONSOLIDATION PARAMETERS
	WT 	<- 1					  # WEIGHTING FUNCTION FOR FLOW PARAMETERS
	#DC 	<- DD					  # SET GRAIN SIZE FOR DCCS TESTS

	Z1	<- EAC  # Predicted axial strain (initial values)
	Z2	<- ELC  # Predicted lateral strain (initial values)
	Z3	<- 0

	# ==== define the differential equation ====
  # ---- only calculate strain rates at TIME > 0 ----
 # browser()
 DZ <- ifelse(cbind(TIME > 0, TIME > 0, TIME > 0),
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
  FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                    exp(-DELTA * (1 - Z3 / EFT) ^ 2)))
  
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

cbind(DZ1, DZ2, DZ3)},{cbind(0,0,0)})
return(DZ)}
# ======================================================
integrate.trap <- function(RATE,TIME){
  # performs integral with respect to time on Strain Rate
  
  TIME <- rep(seq(0, 5, 1),3)
  
}

	
