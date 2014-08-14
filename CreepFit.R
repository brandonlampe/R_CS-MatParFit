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
ETA0.INP  <- 0.102854        # -
ETA1.INP 	<- 3.9387          # -
ETA2.INP 	<- 1               # constant -
NF.INP 		<- 3.5122          # -
AA1.INP 	<- 0.3147          # -
PP.INP 		<- 1.6332          # -
NSP.INP 	<- 0.557621        # -
R1.INP 		<- 1.041 * 10 ^ -6 # [K/(MPa-sec)]
R3.INP 		<- 15.1281         # -
R4.INP 		<- 0.1677765       # -
QSR.INP 	<- 1077.46         # [K]

# ---- INITIAL CREEP PARAMETER VALUES ----
CREEP.INP <- cbind(ETA0.INP, ETA1.INP, ETA2.INP, NF.INP, AA1.INP, PP.INP,
                   NSP.INP, R1.INP, R3.INP, R4.INP, QSR.INP)


PAR.CREEP <- nls.lm(par = CREEP.INP, fn = creep_RESID,
                    xx = DATA.INP, shear = FLOW.INP, control = nls.lm.control(
                      nprint=1, maxiter = 1000, maxfev = 10000))

str(PAR.CREEP)

print(PAR.CREEP)
ETA0.FIT  <- PAR.CREEP[[1]][1]  # -
ETA1.FIT  <- PAR.CREEP[[1]][2]  # -
ETA2.FIT 	<- PAR.CREEP[[1]][3]  # constant -
NF.FIT 		<- PAR.CREEP[[1]][4]  # -
AA1.FIT 	<- PAR.CREEP[[1]][5]  # -
PP.FIT 		<- PAR.CREEP[[1]][6]  # -
NSP.FIT 	<- PAR.CREEP[[1]][7]  # -
R1.FIT 		<- PAR.CREEP[[1]][8]  # [K/(MPa-sec)]
R3.FIT 		<- PAR.CREEP[[1]][9]  # -
R4.FIT 		<- PAR.CREEP[[1]][10] # -
QSR.FIT 	<- PAR.CREEP[[1]][11] # [K]

#---- put calculated parameters into data frame
FIT_OUT.CREEP <- data.frame(matrix(data = NA, nrow = 2, ncol = 11))
FIT_OUT.CREEP[1,] <- data.frame(ETA0.FIT, ETA1.FIT, ETA2.FIT, NF.FIT, AA1.FIT, PP.FIT,
                            NSP.FIT, R1.FIT, R3.FIT, R4.FIT, QSR.FIT)

colnames(FIT_OUT.CREEP) <- c("ETA0", "ETA1", "ETA2", "NF", "AA1",
                                 "PP","NSP", "R1", "R3", "R4", "QSR")
rownames(FIT_OUT.CREEP) <- c("Final", "Intial")

FIT_OUT.CREEP[2,] <- CREEP.INP
print(FIT_OUT.CREEP)

# ---- lost SSE into matrix ----
SSE_OUT.CREEP <- PAR.CREEP[[8]]

#==== capture output in text file ====
FIT_FILE.CREEPpar <- paste(CurrentDirectory,"ParameterFits_CREEP.OUT",sep = "/")
FIT_FILE.CREEPsse <- paste(CurrentDirectory,"SSE_CREEP.OUT",sep = "/")

write.table(FIT_OUT.CREEP,file = FIT_FILE.CREEPpar, sep = ",",
            col.names = colnames(FIT_OUT.CREEP))
write.table(SSE_OUT.CREEP,file = FIT_FILE.CREEPsse, sep = ",",
            col.names = "Sum of Squared Error")

# # ---- check data ----
# # library(ggplot2)
# ggSUB.EAC <- ggplot(data = SUB.FEAC, aes(x=TIME, y=EAC))
# ggSUB.EAC <- ggSUB.EAC + geom_line()
# ggSUB.EAC <- ggSUB.EAC + geom_point(aes(y=FEAC))
# # ggSUB.EAC <- ggSUB.EAC + facet_wrap(~ITEST, ncol=3, scales = "free")
#  ggSUB.EAC <- ggSUB.EAC + xlim(0,6e6) + ylim(-0.25,0)
# ggSUB.EAC <- ggSUB.EAC + ylab("Axial Strain: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
# ggSUB.EAC
#
# # ---- MEASURED & FIT AXIAL STRAINS VS TIME (FACETED) ----
# ggF.EAC <- ggplot(data = DT.FE, aes(x=TIME, y=EAC, color = ITEST))
# ggF.EAC <- ggF.EAC + geom_point()
# ggF.EAC <- ggF.EAC + geom_line(aes(y=IFEAR)) + facet_wrap(~ITEST, ncol=3, scales = "free")
# ggF.EAC <- ggF.EAC + facet_wrap(~ITEST, ncol=3, scales = "free")
# ggF.EAC <- ggF.EAC + xlim(0,6e6)
# ggF.EAC <- ggF.EAC + ylab("Axial Strain: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
# ggF.EAC
#
# ggF.EAR <- ggplot(data = DATA.FIT, aes(x=TIME, y=EAR, color = ITEST))
# ggF.EAR <- ggF.EAR + geom_point()
# ggF.EAR <- ggF.EAR + geom_line(aes(y=FEAR)) + facet_wrap(~ITEST, ncol=3, scales = "free")
# ggF.EAR <- ggF.EAR + facet_wrap(~ITEST, ncol=3, scales = "free")
# ggF.EAR <- ggF.EAR + xlim(0,6e6)
# ggF.EAR <- ggF.EAR + ylab("Axial Strain: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
# ggF.EAR

