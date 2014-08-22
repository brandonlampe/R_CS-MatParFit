source("CS_DE.R")

PAR.TEST <- DATA.INP[which(DATA.INP$ITEST == "SC1B"),] # SUBSET OF DATA FOR ANALYSIS

# ---- linear interpolation functions to be called in ODE ----
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
P.CER <- ode(func = strain_Rates.01, parms = PARM, y = IC,
             times = TIME, method = "impAdams")

