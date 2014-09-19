source("CS_DE.R")

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
             times = TIME)#, events = list(func = eventfun, time = 0),
             #method = "impAdams")

# ==== DEBUG ====
ODE.DT <- data.table(P.CER)
setnames(ODE.DT,c("TIME", "Z1", "Z2", "Z3",
                  "MD", "FU", "FU.TEST", "ESS",
                  "ES1", "ES2", "ES3", "SP", "DZ1",
                  "DZ2", "DZ3", "F2A", "F2L", "EFT", "SEQ",
                  "SEQF", "BIGD", "DEBUG.FU", "DEBUG.GAMMA",
                  "DEBUG.ES3", "AS", "LS"))

plot(P.CER[1:100,"time"], P.CER[1:100,"Z1"], ylim = c(-.15,0), col = 27,
     xlab = "Time", ylab = "Axial Strain") # green
par(new = T)
plot(PAR.TEST$TIME[1:100], PAR.TEST$EAC[1:100], ylim= c(-.15,0),xlab = "Time",
     ylab = "Axial Strain") # MEASURED AXIAL STRAIN
par(new = F)

# # ==== DIAGNOSTICS BELOW ====
# # diagnostics(P.CER)
# plot(PAR.TEST$TIME, PAR.TEST$EAR) # MEASURED AXIAL STRAIN RATE
#
# plot(P.CER[,"time"], P.CER[,"Z1"])  # axial strain
# plot(P.CER[,"time"], P.CER[,"Z2"])  # lateral strain
# plot(P.CER[,"time"], P.CER[,"Z3"])  # xi
#
# plot(P.CER[1:25,"time"], P.CER[1:25,"Z1"])
# plot(P.CER[1:25,"time"], P.CER[1:25,"Z2"])
# # plot(P.CER[1:25,"time"], P.CER[1:25,"Z3"])
#
# COMP <- cbind(TIME = PAR.TEST$TIME, EAC = PAR.TEST$EAC, PEAC = P.CER[,"Z1"],
#               REAC = (PAR.TEST$EAC / P.CER[,"Z1"]))
#
# plot(COMP[,"TIME"], COMP[,"REAC"])
#
# x <- cbind(a = 1:3, pi = pi) # simple matrix with dimnames
# attributes(x)
#
# ## strip an object's attributes:
# attributes(x) <- NULL
# x # now just a vector of length 6
#
# mostattributes(x) <- list(mycomment = "really special", dim = 3:2,
#                           dimnames = list(LETTERS[1:3], letters[1:5]), names = paste(1:6))
