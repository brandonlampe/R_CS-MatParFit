source("DoAll.R")

# ---- compare test results for SC1B, LAMPE VS BEAN
DATA.BEAN.RAW <- read.table("/Users/Lampe/GrantNo456417/MatParameterFitting/fromJimBean/1elem_SC1B_out.txt",
                        skip = 1)
CNAME <- read.table("/Users/Lampe/GrantNo456417/MatParameterFitting/fromJimBean/1elem_SC1B_out.txt",
                    nrows = 1)
DATA.BEAN.RAW2 <- read.table("/Users/Lampe/GrantNo456417/MatParameterFitting/fromJimBean/1elem_SC1B.xmgr.out",
                             skip = 2)

colnames(DATA.BEAN.RAW) <- c("TIME", "Z3","EFT", "SEQ", "SEQF", "OMEGAK", "FU",
                         "ESS", "AS", "LS")
colnames(DATA.BEAN.RAW2) <- c("TIME", "Z2","Z1", "OMEGAK", "OMEGAA","DEN",
                             "MD", "SP")

DATA.BEAN.RAW2$TIME <- NULL  #REPEATED DATA
DATA.BEAN.RAW2$OMEGAK <- NULL#REPEATED DATA

# remove rows that extend beyond test data and combine files
DATA.BEAN <- data.frame(cbind(DATA.BEAN.RAW,DATA.BEAN.RAW2))
DATA.BEAN[,"NAME"] <- "BEAN"
DATA.BEAN[,"ITEST"] <- "SC1B"
#write.table(DATA.BEAN,file = "dataFromJimBean_2014-09-19.txt", sep = ",")

###########################################################################
# SOLVE FOR THE STRAINS USING CALLAHAN'S MODEL
# TEST SC1B
###########################################################################

# ==== "events" function, for specifying step functions ====
# ---- SETS STRAINS AND RATES = 0 AT TIME == 0 ----
eventfun <- function(Time, State, Parm){
  with (as.list(State),{

    # ---- if time == 0, derivative = 0 ----
    DZ <- ifelse(c(Time == 0, Time == 0, Time == 0), {c(0, 0, 0)}, {DZ})

    return(DZ)
  })
}

# ==== DIFFERENTIAL EQUATION ====
STRAINS.02 <- function(Time, State, Parm){
  with(as.list(c(State, Parm)),{
    # function for calculating axial and lateral strain rates
    # Input must be in vector or matrix form, no data frames
    # Eqns. referenced from:  SAND97-2601
    # CPar: EAT0, ETA1, ETA2, NF, AA1, PP, NSP, R1, R3, R4, QSR
    # FPar: KAP0, KAP1, KAP2, NK, DDT
    # TestData:

    # ============= parameters hard coded into function directly ========
    # browser()
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
    R1 		<- 1.041 * 10 ^ -6 # [K/(MPa-sec)]
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
    #   ELC	<- (EVC - EAC) / 2		  # CREEP TRUE LATERAL STRAIN
    D0 	<- 1382.4 / RHOIS			  # EMPLACED FRACTIONAL DENSITY (0.64 FRAC DENSITY)
    DI 	<- RHOI / RHOIS			    # FRACTIONAL DENSITY at the start of creep

    # ==============================================================
    # integral of Eqn 2-27, (initial values)

    # ==== define the differential equation ====
    # browser()
    VOL  	<- Z1 + 2*Z2			    # TRUE VOLUMETRIC STRAIN
    VOLT	<- VOL + log(D0/DI)	# VOLUMETRIC STRAIN + INITIAL TRUE STRAIN ESTIMATE
    DEN		<- DI/exp(VOL)			# CURRENT FRACTIONAL DENSITY

    VAR <- ifelse(DEN <= DDT, DDT, DEN) # DEFINE DENSITY floor ISH
    # ==== DEBUG ====
    DEBUG.VAR <- ifelse(DEN <= DDT, 1, -1)

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
    # ==== DEBUG ====
    DEBUG.ES3 <- ifelse(H >=0, 1, -1)

    ESS = ES1 + ES2 + ES3 # Steady-state strain rate, Eqn. 2-29 (SAND97-2601)

    # ---- EVALUATE TRANSIENT FUNCTION, 3 branches: work hardening, equilibrium, recovery
    EFT  <- K0 * exp(C * TEMP) * (SEQF / MU) ^ M  # Transient Strain Limit, Eqn. 2-28
    BIGD <- ALPHA + BETA * log10(SEQF / MU)       # Work-Hardening parameter, Eqn 2-28

    FU <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                      exp(-DELTA * (1 - Z3 / EFT) ^ 2)))

    # ==== DEBUG ====
    DEBUG.FU <- ifelse(Z3 == EFT, 0, ifelse(Z3 < EFT, -1, 1))

    # equivalent inelastic strain rate form for dislocation creep, Eqn 2-23
    MD <- FU * ESS

    # ==== START: Equivalent Inelastic Strain Rate Form for Pressure Solutioning ====
    # ---- Calculate ENGINEERING volumetric strain ----
    CR <- abs(exp(VOLT) - 1) # USES THE DEFINITION OF ENGINEERING STRAIN

    # ---- Determine functional form - either large or small strains, Eqn 2-34 ----
    GAMMA <- ifelse(CR <= 0.15, 1, (abs((D0 - exp(VOLT)) /
                                          ((1 - D0) * exp(VOLT)))) ^ NSP)
    # Small Strains (Vol Strain > - 15%)
    # Large Strains (Vol Strain < - 15%)
    # ==== DEBUG ====
    DEBUG.GAMMA <- ifelse(CR <= 0.15, 1,-1)

    # ---- component of eqn 2-35 ---
    X3 <- exp((R3 - 1) * VOLT) / (abs(1 - exp(VOLT))) ^ R4

    # ---- determine value of moisture function (w) ----
    M2 <- ifelse (W == 0, 0, W ^ AA1)

    # ---- Equivalent Inelastic Strain Rate Form for Pressure Solutioning, Eqn 2-35
    G2 <- 1 / DD ^ PP # calculate grain size function
    T2 <- exp(-QSR / TEMP) / TEMP

    R1 <- 0.0194
    SP <- R1 * M2 * G2 * T2 * X3 * GAMMA * SEQF #})

    DZ1 <- (MD + SP) * F2A # derivative: axial strain rate
    DZ2 <- (MD + SP) * F2L # derivative: lateral strain rate
    DZ3 <- (FU - 1) * ESS  # derivative of internal variable "zeta"
    #     browser()
    DZ <- list(c(DZ1, DZ2, DZ3), MD, FU, ESS, ES1, ES2, ES3, SP, DZ1, DZ2,
               DZ3, F2A, F2L, EFT, SEQ, SEQF, BIGD, DEBUG.GAMMA,
               DEBUG.ES3, AS, LS, OMEGAA, OMEGAK,VAR, DEN, GAMMA)

    return(DZ)
  })
}

TestName <- "SC1B"
PAR.TEST <- DATA.INP[which(DATA.INP$ITEST == TestName),] # SUBSET OF DATA FOR ANALYSIS
# debug.out <- paste(CurrentDirectory,"debug_SC1B.csv",sep = "/")
# write.table(ODE.DT,file = debug.out, sep = ",")

# ---- linear interpolation functions to be called in "strain_Rates.01" ----
time.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TIME)
temp.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TEMP)
# as.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$AS)
# ls.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$LS)
as.interp   <- approxfun(x = DATA.BEAN$TIME, y = DATA.BEAN$AS / 10^6)
ls.interp   <- approxfun(x = DATA.BEAN$TIME, y = DATA.BEAN$LS / 10^6)

RHOI   <- as.numeric(PAR.TEST$RHOI[1])  # DENSITY AT THE START OF CREEP
DD 		<- as.numeric(PAR.TEST$DD[1])	# AVERAGE GRAIN SIZE [MM]
W 		<- as.numeric(PAR.TEST$W[1])	# WATER CONENT BY PERCENT WEIGHT

PARM <- c(RHOI = RHOI, DD = DD, W = W) # CONSTANT TEST SPECIFIC PARAMETERS

# ---- intial values for state variables ----
Z1	<- 0 # Predicted axial strain (initial values)
Z2	<- 0 # Predicted lateral strain (initial values)
Z3	<- 0 # internal variable "xi" for the transient function (FU)
# integral of Eqn 2-27, (initial values)

IC <- (c(Z1 = Z1, Z2 = Z2, Z3 = Z3)) # array of initial values

TIME <- DATA.BEAN$TIME

# ---- function for Predicting the Creep Strain(E) Rates ----
P.CER <- ode(func = STRAINS.02, parms = PARM, y = IC,
             times = TIME, verbose = TRUE)#, hmax = 500, maxsteps=10000 )

# ==== DEBUG ====
ODE.DT <- data.table(P.CER)
setnames(ODE.DT,c("TIME", "Z1", "Z2", "Z3",
                  "MD", "FU", "ESS",
                  "ES1", "ES2", "ES3", "SP", "DZ1",
                  "DZ2", "DZ3", "F2A", "F2L", "EFT", "SEQ",
                  "SEQF", "BIGD", "DEBUG.GAMMA",
                  "DEBUG.ES3", "AS", "LS", "OMEGAA", "OMEGAK", "VAR", "DEN","GAMMA"))

DATA.LAMPE <- data.frame(ODE.DT$TIME,ODE.DT$Z3,ODE.DT$EFT, ODE.DT$SEQ * 10^6,
                         ODE.DT$SEQF * 10^6,ODE.DT$OMEGAK, ODE.DT$FU,
                         ODE.DT$ESS, ODE.DT$AS * 10^6, ODE.DT$LS * 10^6,
                         ODE.DT$Z2,ODE.DT$Z1, ODE.DT$OMEGAA, ODE.DT$DEN,
                         ODE.DT$MD, ODE.DT$SP)

colnames(DATA.LAMPE) <- c("TIME","Z3","EFT", "SEQ", "SEQF", "OMEGAK", "FU",
                         "ESS", "AS", "LS", "Z2", "Z1", "OMEGAA", "DEN",
                         "MD", "SP")
DATA.LAMPE[,"NAME"] <- "LAMPE"
DATA.LAMPE[,"ITEST"] <- TestName

DATA.ALL <- rbind(DATA.BEAN, DATA.LAMPE)

# # add ggplot here
P.DEN <- ggplot(DATA.ALL, aes(x = TIME, color = NAME))
P.DEN <- P.DEN + geom_line(aes(y = DEN))
P.DEN <- P.DEN + ylab("FRACTIONAL DENSITY")
P.DEN

P.ZETA <- ggplot(DATA.ALL, aes(x = TIME, y = Z3, color = NAME))
P.ZETA <- P.ZETA + geom_line()
P.ZETA <- P.ZETA + ylab("INTERNAL VARIABLE (ZETA)")
P.ZETA

P.Z1 <- ggplot(DATA.ALL, aes(x = TIME, y = Z1, color = NAME))
P.Z1 <- P.Z1 + geom_line()
P.Z1 <- P.Z1 + ylab("AXIAL STRAIN")
P.Z1

P.Z2 <- ggplot(DATA.ALL, aes(x = TIME, y = Z2, color = NAME))
P.Z2 <- P.Z2 + geom_line()
P.Z2 <- P.Z2 + ylab("LATERAL STRAIN")
P.Z2

P.AS <- ggplot(DATA.ALL, aes(x = TIME, y = AS, color = NAME))
P.AS <- P.AS + geom_line()
P.AS <- P.AS + ylab("AXIAL STRESS [PA]")
P.AS

P.LS <- ggplot(DATA.ALL, aes(x = TIME, y = LS, color = NAME))
P.LS <- P.LS + geom_line()
P.LS <- P.LS + ylab("LATERAL STRESS [PA]")
P.LS

P.SEQ <- ggplot(DATA.ALL, aes(x = TIME, y = SEQ, color = NAME))
P.SEQ <- P.SEQ + geom_line()
P.SEQ <- P.SEQ + ylab("EQUIVALENT STRESS FOR FLOW")
P.SEQ

P.SEQF <- ggplot(DATA.ALL, aes(x = TIME, y = SEQF, color = NAME))
P.SEQF <- P.SEQF + geom_line()
P.SEQF <- P.SEQF + ylab("EQUIVALENT STRESS FOR CREEP")
P.SEQF

P.OMEGAA <- ggplot(DATA.ALL, aes(x = TIME, y = OMEGAA, color = NAME))
P.OMEGAA <- P.OMEGAA + geom_line()
P.OMEGAA <- P.OMEGAA + ylab("FRACTION DENSITY FOR EQUIVALENT STRESS (CREEP)")
P.OMEGAA

P.OMEGAK <- ggplot(DATA.ALL, aes(x = TIME, y = OMEGAK, color = NAME))
P.OMEGAK <- P.OMEGAK + geom_line()
P.OMEGAK <- P.OMEGAK + ylab("FRACTION DENSITY FOR EQUIVALENT STRESS (FLOW)")
P.OMEGAK

P.FU <- ggplot(DATA.ALL, aes(x = TIME, y = FU, color = NAME))
P.FU <- P.FU + geom_line()
P.FU <- P.FU + ylab("TRANSIENT FUNCTION (F)")
P.FU

P.EFT <- ggplot(DATA.ALL, aes(x = TIME, y = EFT, color = NAME))
P.EFT <- P.EFT + geom_line()
P.EFT <- P.EFT + ylab("TRANSIENT STRAIN LIMIT")
P.EFT

P.ESS <- ggplot(DATA.ALL, aes(x = TIME, y = ESS, color = NAME))
P.ESS <- P.ESS + geom_line()
P.ESS <- P.ESS + ylab("STEADY STATE CREEP STRAIN RATE")
P.ESS

P.MD <- ggplot(DATA.ALL, aes(x = TIME, y = MD, color = NAME))
P.MD <- P.MD + geom_line()
P.MD <- P.MD + ylab("DISLOCATION CREEP STRAIN RATE")
P.MD

P.SP <- ggplot(DATA.ALL, aes(x = TIME, y = SP, color = NAME))
P.SP <- P.SP + geom_line()
P.SP <- P.SP + ylab("PRESSURE SOLUTIONING CREEP STRAIN RATE")
P.SP


# # # ---- EXPORT AND SAVE FILES AS PDF
PATH = "/Users/Lampe/GrantNo456417/MatParameterFitting/R_CS-MatParFit"
FILE.NAME = "/CS_LAMPE-BEAN-03.pdf"
pdf(file = paste(PATH, FILE.NAME, sep = ""), onefile = TRUE)

P.Z1
P.Z2
P.MD
P.SP
P.ESS
P.EFT
P.ZETA
P.FU
P.DEN
P.SEQ
P.SEQF
P.OMEGAA
P.OMEGAK
P.AS
P.LS

dev.off()
