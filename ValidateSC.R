# Routine to calculate strains based on input stresses

FEAR.NZ <- ERATE.OUT[,1] # fit of axial strain rates
FELR.NZ <- ERATE.OUT[,2] # fit of lateral strain rates
FRAT.NZ <- FELR.NZ / FEAR.NZ 

CNAMES <- colnames(DATA.NZ)#, FEAR.NZ, FELR.NZ, FRAT.NZ)
PLOT.DATA <- DATA.NZ#
colnames(PLOT.DATA) <- CNAMES

PLOT.DATA <- cbind(PLOT.DATA, FEAR.NZ, FELR.NZ, FRAT.NZ) # MERGE TEST AND FIT DATA


# ---- plot lateral-to-axial strain rate ratios ----
ggVal.SC <- ggplot(data = PLOT.DATA, aes(x=TIME.NZ, y=RAT.NZ, color = ITEST.NZ))
ggVal.SC <- ggVal.SC + geom_point()
ggVal.SC <- ggVal.SC + geom_line(aes(y=FRAT.NZ))+ facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggVal.SC <- ggVal.SC + facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggVal.SC <- ggVal.SC + xlim(0,6e6)
ggVal.SC <- ggVal.SC + ylab("Lateral-To-Axial Strain Rates: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggVal.SC

# ---- plot axial strain rate ratios ----
ggEAR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME.NZ, y=EAR.NZ, color = ITEST.NZ))
ggEAR.SC <- ggEAR.SC + geom_point()
ggEAR.SC <- ggEAR.SC + geom_line(aes(y=FEAR.NZ))+ facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggEAR.SC <- ggEAR.SC + facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggEAR.SC <- ggEAR.SC + xlim(0,6e6)
ggEAR.SC <- ggEAR.SC + ylab("Axial Strain Rates: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggEAR.SC

# ---- plot lateral strain rate ratios ----
ggELR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME.NZ, y=ELR.NZ, color = ITEST.NZ))
ggELR.SC <- ggELR.SC + geom_point()
ggELR.SC <- ggELR.SC + geom_line(aes(y=FELR.NZ))+ facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggELR.SC <- ggELR.SC + facet_wrap(~ITEST.NZ, ncol=3, scales = "free")
ggELR.SC <- ggELR.SC + xlim(0,6e6)
ggELR.SC <- ggELR.SC + ylab("Lateral Strain Rates: Calculated (dot) Vs. Measured (line)") + xlab("Time [sec]")
ggELR.SC