# Fit Parameters shear consolidation parameters

source("Func.R")
#====================================
AS <- as.numeric(DATA.NZ$AS)          # axial stress [MPa]
LS <- as.numeric(DATA.NZ$LS)          # lateral stress [MPa]
D  <- as.numeric(DATA.NZ$D)           # fraction density

InVar    <- cbind(AS,LS,D)            # varibales input into nls.LM for fitting
RAT.MEAS <- ifelse(DATA.NZ$TIME > 0,
                   DATA.NZ$ELR /
                   DATA.NZ$EAR, 0) # measured lateral-to-axial strain rates

# ---- starting parameter values (from Callahan) ----
KAP0.INP <- 10.119
KAP1.INP <- 1.005
DDT.INP  <- 0.896
NK.INP   <- 1.331

# ---- approximations from initial run (Lampe) ----
# KAP0.INP <- 2.09011272735447
# KAP1.INP <- 1.32951908075078
# DDT.INP  <- 0.900012234998625
# NK.INP   <- 2.75382483799085

ParaStart <- cbind(KAP0.INP, KAP1.INP,
                   DDT.INP, NK.INP) # starting values for fiting parameters

PAR.SHEAR <- nls.lm(par = ParaStart, fn = shear_RESID,
                    Observed = RAT.MEAS, xx = InVar,
                    control = nls.lm.control(
                      nprint=1, maxiter = 1000, maxfev = 10000))

#---- summary of fitting results ----
str(PAR.SHEAR)

#---- Calculated parameter values
KAP0.FIT <- PAR.SHEAR[[1]][1]
KAP1.FIT <- PAR.SHEAR[[1]][2]
DDT.FIT  <- PAR.SHEAR[[1]][3]
NK.FIT   <- PAR.SHEAR[[1]][4]

#---- put calculated parameters into data frame
FIT_OUT.SC <- data.frame(KAP0.FIT, KAP1.FIT, DDT.FIT, NK.FIT)
names(FIT_OUT.SC)[1] <- paste("KAP0")
names(FIT_OUT.SC)[2] <- paste("KAP1")
names(FIT_OUT.SC)[3] <- paste("DDT")
names(FIT_OUT.SC)[4] <- paste("NK")

FIT_OUT.SC["Initial",] <- ParaStart
row.names(FIT_OUT.SC)[1] <- "Final"
print(FIT_OUT.SC)

# ---- SSE into matrix ----
SSE_OUT.SC <- PAR.SHEAR[[8]]

#==== capture output in text file ====
FIT_FILE.SCpar <- paste(CurrentDirectory,"ParameterFits_SC.OUT",sep = "/")
FIT_FILE.SCsse <- paste(CurrentDirectory,"SSE_SC.OUT",sep = "/")

write.table(FIT_OUT.SC,file = FIT_FILE.SCpar, sep = ",",
            col.names = colnames(FIT_OUT.SC))
write.table(SSE_OUT.SC,file = FIT_FILE.SCsse, sep = ",",
            col.names = "Sum of Squared Error")

#==== output summary ====
#	A list with components:
#
#	par	 The best set of parameters found.
#
#	hessian	 A symmetric matrix giving an estimate of the Hessian at the solution found.
#				fvec The result of the last fn evaluation; that is, the residuals.
#
#	info	 info is an integer code indicating the reason for termination.
#			0 Improper input parameters.
#			1 Both actual and predicted relative reductions in the sum of squares are at
#				most ftol.
#			2 Relative error between two consecutive iterates is at most ptol.
#			3 Conditions for info = 1 and info = 2 both hold.
#			4 The cosine of the angle between fvec and any column of the Jacobian is at
#				most gtol in absolute value.
#			5 Number of calls to fn has reached maxfev.
#			6 ftol is too small. No further reduction in the sum of squares is possible.
#			7 ptol is too small. No further improvement in the approximate solution par is
# 				possible.
#			8 gtol is too small. fvec is orthogonal to the columns of the Jacobian to
#				machine precision.
#			9 The number of iterations has reached maxiter. message character string
#				indicating reason for termination.
#
#	diag	 The result list of diag. See Details.
#
#	niter	 The number of iterations completed before termination.
#
#	rsstrace The residual sum of squares at each iteration. Can be used to check the
#			 progress each iteration.
#
#	deviance The sum of the squared residual vector.