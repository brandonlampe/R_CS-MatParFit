#---- define parameters----
# this routine labels the parameters loaded via "Load.R" in the data frame "DATA.ALL"
# blank rows or other misc issues with data format are dealt with in this routine

#---- define columns data frame contained measured data ----
ICASE.MEAS	<- as.numeric(DATA.MEAS$V1)	   # TEST TYPE
ITEST.MEAS <- as.character(DATA.MEAS$V2)	 # TEST ID
TIME.MEAS 	<- as.numeric(DATA.MEAS$V3)    # TIME [SEC]
DT.MEAS 		<- as.numeric(DATA.MEAS$V4)	   # DELTA TIME [SEC]
TF.MEAS 		<- as.numeric(DATA.MEAS$V5)	   # TOTAL TEST TIME [SEC]
TEMP.MEAS 	<- as.numeric(DATA.MEAS$V6)	   # TEMP [K]
AS.MEAS 		<- as.numeric(DATA.MEAS$V7)	   # AXIAL STRESS [MPA]
LS.MEAS 		<- as.numeric(DATA.MEAS$V8)	   # LATERAL STRESS [MPA]
EVT.MEAS 	<- as.numeric(DATA.MEAS$V9)	     # TOTAL TRUE VOLUMETRIC STRAIN
EVC.MEAS 	<- as.numeric(DATA.MEAS$V10)	   # CREEP TRUE VOLUMETRIC STRAIN
EAT.MEAS 	<- as.numeric(DATA.MEAS$V11)	   # TOTAL TRUE AXIAL STRAIN
EAC.MEAS 	<- as.numeric(DATA.MEAS$V12)	   # CREEP TRUE AXIAL STRAIN
RHO.MEAS 	<- as.numeric(DATA.MEAS$V13)	   # CURRENT DENSITY [KG/M3]
D.MEAS 		<- as.numeric(DATA.MEAS$V14)	   # FRACTIONAL DENSITY
RHO0.MEAS 	<- as.numeric(DATA.MEAS$V15)   # DENSITY AT THE START OF CONSOLIDATION (<RHOI)
RHOI.MEAS 	<- as.numeric(DATA.MEAS$V16)   # DENSITY AT THE START OF CREEP
DD.MEAS 		<- as.numeric(DATA.MEAS$V17)   # AVERAGE GRAIN SIZE [MM]
W.MEAS 		<- as.numeric(DATA.MEAS$V18)	   # WATER CONENT BY PERCENT WEIGHT
EVR.MEAS   <- as.numeric(DATA.MEAS$V19)    # true volumetric strain rate (/sec)
EAR.MEAS   <- as.numeric(DATA.MEAS$V20)	   # true axial strain rate (/sec)
ELR.MEAS   <- as.numeric(DATA.MEAS$V21)	   # true lateral strain rate (/sec)

# ----  # lateral-to-axial strain rates ----
RAT.MEAS   <- ifelse(DATA.MEAS$V3 > 0, {
  as.numeric(DATA.MEAS$V21) / as.numeric(DATA.MEAS$V20)},{0})

DATA.MEAS$V22 <- RAT.MEAS                 # assign RAT to column in data frame
ELC.MEAS   <- 0.5*(EVC.MEAS - EAC.MEAS)    # true lateral creep strain
DATA.INP <- cbind(DATA.MEAS, ELC.MEAS)    # create data frame of input values

#---- define columns in 'DATA.INP' -> NON ZERO TIME VALUES ----
ICASE.INP <- as.numeric(DATA.INP$ICASE)	  # TEST TYPE
ITEST.INP <- as.character(DATA.INP$V2)	# TEST ID
TIME.INP  <- as.numeric(DATA.INP$V3)    # TIME [SEC]
DT.INP 	 <- as.numeric(DATA.INP$V4)	  # DELTA TIME [SEC]
TF.INP 	 <- as.numeric(DATA.INP$V5)	  # TOTAL TEST TIME [SEC]
TEMP.INP  <- as.numeric(DATA.INP$V6)	  # TEMP [K]
AS.INP 	 <- as.numeric(DATA.INP$V7)	  # AXIAL STRESS [MPA]
LS.INP 	 <- as.numeric(DATA.INP$V8)	  # LATERAL STRESS [MPA]
EVT.INP 	 <- as.numeric(DATA.INP$V9)	  # TOTAL TRUE VOLUMETRIC STRAIN
EVC.INP 	 <- as.numeric(DATA.INP$V10)	  # CREEP TRUE VOLUMETRIC STRAIN
EAT.INP 	 <- as.numeric(DATA.INP$V11)	  # TOTAL TRUE AXIAL STRAIN
EAC.INP 	 <- as.numeric(DATA.INP$V12)	  # CREEP TRUE AXIAL STRAIN
RHO.INP 	 <- as.numeric(DATA.INP$V13)	  # CURRENT DENSITY [KG/M3]
D.INP 		 <- as.numeric(DATA.INP$V14)	  # FRACTIONAL DENSITY
RHO0.INP  <- as.numeric(DATA.INP$V15)	  # DENSITY AT THE START OF CONSOLIDATION (<RHOI)
RHOI.INP  <- as.numeric(DATA.INP$V16)	  # DENSITY AT THE START OF CREEP
DD.INP 	 <- as.numeric(DATA.INP$V17)	  # AVERAGE GRAIN SIZE [MM]
W.INP 		 <- as.numeric(DATA.INP$V18)	  # WATER CONENT BY PERCENT WEIGHT
EVR.INP   <- as.numeric(DATA.INP$V19)   # true volumetric strain rate (/sec)
EAR.INP   <- as.numeric(DATA.INP$V20)	  # true axial strain rate (/sec)
ELR.INP   <- as.numeric(DATA.INP$V21)	  # true lateral strain rate (/sec)
RAT.INP   <- as.numeric(DATA.INP$V22)   # measured lateral-to-axial strain rate ratio
ELC.INP   <- as.numeric(DATA.INP$V23)    # true lateral creep strain

colnames(DATA.INP) <- c("ICASE", "ITEST", "TIME", "DT", "TF", "TEMP", "AS",
                       "LS", "EVT", "EVC", "EAT", "EAC", "RHO", "D", "RHO0",
                       "RHOI", "DD", "W", "EVR", "EAR", "ELR", "RAT", "ELC")

# ---- CHECK IMPORTED DATA TO BE SURE IT'S ALL THERE! ----
TIME.CHECK <- sum(TIME.INP == 0)
print(c("CHECK TEST IMPORT #:", TIME.CHECK))

RAT.CHECK <- sum(is.na(ELR.INP / EAR.INP))
RAT.LOC   <- which(is.na(ELR.INP / EAR.INP))
#print(c("RAT CHECK:  CNT:", RAT.CHECK," LOCATION: ",RAT.LOC, "TEST ID", DATA.INP$ITEST[RAT.LOC]))

EAR.CHECK <- sum(EAR.INP == 0)
EAR.LOC   <- which(EAR.INP == 0)
#print(c("EAR CHECK:  CNT:", EAR.CHECK," LOCATION: ",EAR.LOC, "TEST ID", DATA.INP$ITEST[EAR.LOC]))

# ---- SUMMARY OF DATA ----
DATA.SUMMARY <- ddply(DATA.INP, c("ITEST"), function(df) c(
  mean(df$ICASE), max(df$TIME / 60 / 60 / 24), min(df$AS), max(df$AS), mean(df$AS),
  min(df$LS), max(df$LS), mean(df$LS),mean(df$LS) - mean(df$AS), min(df$EVC),
  max(df$EVC), min(df$EAC), max(df$EAC), mean(df$TEMP), min(df$D), max(df$D), 
  mean(df$W), min(df$EVR), max(df$EVR), min(df$EAR), max(df$EAR))) 

colnames(DATA.SUMMARY) <- c("ITEST","ICASE","TIME [DAY]", "S.A.MIN", "S.A.MAX",
                            "S.A.AVG", "S.L.MIN", "S.L.MAX", "S.L.AVG","S.D.AVG",
                            "E.V.MIN", "E.V.MAX", "E.A.MIN", "E.A.MAX", "T.AVG",
                            "FD.MIN","FD.MAX", "W.AVG", "ER.V.MIN", "ER.V.MAX", 
                            "ER.A.MIN","ER.A.MAX")

print(DATA.SUMMARY)
