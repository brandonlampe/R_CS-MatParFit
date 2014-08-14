# CALCULATE AXIAL AND LATERAL STRAIN RATES STRAIN RATES

source("Func.R")

# ==== INPUT FLOW PARAMETERS - THESE WILL REMAIN CONSTAT ====
# ---- FLOW PARAMETERS (from Callahan) ----
KAP0.INP <- 10.119
KAP1.INP <- 1.005
DDT.INP  <- 0.896
NK.INP   <- 1.331
KAP2.INP <- 1

# ---- FLOW PARAMETERS (Lampe) ----
# KAP0.INP <- 2.09011272735447
# KAP1.INP <- 1.32951908075078
# DDT.INP  <- 0.900012234998625
# NK.INP   <- 2.75382483799085
# KAP2.INP <- 1                  # CONSTANT

FLOW.INP<- cbind(KAP0.INP, KAP1.INP, KAP2.INP,
                 DDT.INP, NK.INP) # values for fiting parameters

# ==== INPUT INITIAL VALUES FOR CREEP PARAMETERS ====
# ---- Creep parameters, results from Callahan fits - shear tests only ---- 
ETA0  <- 0.102854        # -
ETA1 	<- 3.9387          # -
ETA2 	<- 1               # constant - 
NF 		<- 3.5122          # - 
AA1 	<- 0.3147          # - 
PP 		<- 1.6332          # - 
NSP 	<- 0.557621        # - 
R1 		<- 1.041 * 10 ^ -6 # [K/(MPa-sec)]
R3 		<- 15.1281         # -
R4 		<- 0.1677765       # -
QSR 	<- 1077.46         # [K]

CREEP.INP <- cbind(ETA0, ETA1, ETA2, NF, AA1, PP,
                   NSP, R1, R3, R4, QSR) # INITIAL CREEP PARAMETER VALUES

# ---- calculate strain rates ----
ERATE.OUT <- strain_Rates(CREEP.INP, FLOW.INP, DATA.INP) # STRAIN RATES
colnames(ERATE.OUT) <- c("FEAR", "FELR", "FEVR")         # column names

# ---- calculate strains (integral of strain rates) ----
DATA.FIT <- cbind(DATA.INP, ERATE.OUT)

DT.DATA.FIT <- data.table(DATA.FIT)       # change format to data table
setkey(DT.DATA.FIT, ITEST)                # create a value for subsets 

# ---- data table containing calculated strains (integrated strain rates) ----
DT.FE <- DT.DATA.FIT[, c("IFEAR", "IFELR", "IFEVR"):=list(
  as.vector(cumtrapz(TIME, FEAR)),
  as.vector(cumtrapz(TIME, FELR)),
  as.vector(cumtrapz(TIME, FEVR))), by = ITEST]

# ---- check data ----
# library(ggplot2)
ggSUB.EAC <- ggplot(data = SUB.FEAC, aes(x=TIME, y=EAC))
ggSUB.EAC <- ggSUB.EAC + geom_line()
ggSUB.EAC <- ggSUB.EAC + geom_point(aes(y=FEAC))
# ggSUB.EAC <- ggSUB.EAC + facet_wrap(~ITEST, ncol=3, scales = "free")
 ggSUB.EAC <- ggSUB.EAC + xlim(0,6e6) + ylim(-0.25,0)
ggSUB.EAC <- ggSUB.EAC + ylab("Axial Strain: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggSUB.EAC

# ---- MEASURED & FIT AXIAL STRAINS VS TIME (FACETED) ----
ggF.EAC <- ggplot(data = DT.FE, aes(x=TIME, y=EAC, color = ITEST))
ggF.EAC <- ggF.EAC + geom_point()
ggF.EAC <- ggF.EAC + geom_line(aes(y=IFEAR)) + facet_wrap(~ITEST, ncol=3, scales = "free")
ggF.EAC <- ggF.EAC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggF.EAC <- ggF.EAC + xlim(0,6e6)
ggF.EAC <- ggF.EAC + ylab("Axial Strain: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
ggF.EAC

ggF.EAR <- ggplot(data = DATA.FIT, aes(x=TIME, y=EAR, color = ITEST))
ggF.EAR <- ggF.EAR + geom_point()
ggF.EAR <- ggF.EAR + geom_line(aes(y=FEAR)) + facet_wrap(~ITEST, ncol=3, scales = "free")
ggF.EAR <- ggF.EAR + facet_wrap(~ITEST, ncol=3, scales = "free")
ggF.EAR <- ggF.EAR + xlim(0,6e6)
ggF.EAR <- ggF.EAR + ylab("Axial Strain: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
ggF.EAR

