source("CS_DE.R")
source("CS_DE-01.R")

PAR.TEST <- DATA.INP[which(DATA.INP$ITEST == "SC1B"),] # SUBSET OF DATA FOR ANALYSIS
# debug.out <- paste(CurrentDirectory,"debug_SC1B.csv",sep = "/")
# write.table(ODE.DT,file = debug.out, sep = ",")

# ---- linear interpolation functions to be called in "strain_Rates.01" ----
time.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TIME)
temp.interp <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$TEMP)
as.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$AS)
ls.interp   <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$LS)
d.interp    <- approxfun(x = PAR.TEST$TIME, y = PAR.TEST$D)

RHOI   <- as.numeric(PAR.TEST$RHOI[1])	# DENSITY AT THE START OF CREEP
DD 		<- as.numeric(PAR.TEST$DD[1])	# AVERAGE GRAIN SIZE [MM]
W 		<- as.numeric(PAR.TEST$W[1])	# WATER CONENT BY PERCENT WEIGHT

PARM <- c(RHOI = RHOI, DD = DD, W = W) # CONSTANT TEST SPECIFIC PARAMETERS

# ---- intial values for state variables ----
Z1	<- 0 # Predicted axial strain (initial values)
Z2	<- 0 # Predicted lateral strain (initial values)
Z3	<- 0 # internal variable "xi" for the transient function (FU)
                           # integral of Eqn 2-27, (initial values)

IC <- (c(Z1 = Z1, Z2 = Z2, Z3 = Z3)) # array of initial values

TIME <- PAR.TEST$TIME

# ---- function for Predicting the Creep Strain(E) Rates ----
P.CER <- dede(func = strain_Rates.01, parms = PARM, y = IC,
             times = TIME, hmax = 1, verbose = TRUE)#, events = list(func = eventfun, time = 0),
             #method = "impAdams")

# P.CER <- dede(func = strain_Rates.02, parms = PARM, y = IC,
#               times = TIME, verbose = TRUE, hmax = 20)

 # ==== DEBUG ====
ODE.DT <- data.table(P.CER)
setnames(ODE.DT,c("TIME", "Z1", "Z2", "Z3",
                  "MD", "FU", "ESS",
                  "ES1", "ES2", "ES3", "SP", "DZ1",
                  "DZ2", "DZ3", "F2A", "F2L", "EFT", "SEQ",
                  "SEQF", "BIGD", "DEBUG.FU", "DEBUG.GAMMA",
                  "DEBUG.ES3", "AS", "LS"))

# add ggplot here
# P.EC <- ggplot(ODE.DT, aes(x = TIME, y = Z1))
# P.EC <- P.EC + geom_point(color = "green")
# P.EC <- P.EC + geom_point(aes(x = TIME, y = Z2, color = "red"))
# P.EC
#
# P.ELC <- ggplot(ODE.DT, aes(x = TIME, y = Z2))
# P.ELC <- P.ELC + geom_point()
# P.ELC
#
P.Z <- ggplot(ODE.DT, aes(x = TIME, y = Z3))
P.Z <- P.Z + geom_point(color = "red")
P.Z
