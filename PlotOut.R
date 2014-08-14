# Routine to calculate strains based on input stresses

FEAR <- ERATE.OUT[,1] # fit of axial strain rates
FELR <- ERATE.OUT[,2] # fit of lateral strain rates
FRAT <- FELR / FEAR 

PLOT.DATA <- cbind(DATA.NZ, ELC, FEAR, FELR, FRAT)

# ---- plot lateral-to-axial strain rate ratios ----
ggVal.SC <- ggplot(data = PLOT.DATA, aes(x=TIME, y=RAT, color = ITEST))
ggVal.SC <- ggVal.SC + geom_point()
ggVal.SC <- ggVal.SC + geom_line(aes(y=FRAT))+ facet_wrap(~ITEST, ncol=3, scales = "free")
ggVal.SC <- ggVal.SC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggVal.SC <- ggVal.SC + xlim(0,6e6)
ggVal.SC <- ggVal.SC + ylab("Lateral-To-Axial Strain Rates: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
ggVal.SC

#=======================================================
# ---- plot axial strain rates: compare measured to fit ----
ggEAR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME, y=EAR, color = ITEST))
ggEAR.SC <- ggEAR.SC + geom_point()
ggEAR.SC <- ggEAR.SC + geom_line(aes(y=FEAR))+ facet_wrap(~ITEST, ncol=3, scales = "free")
ggEAR.SC <- ggEAR.SC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggEAR.SC <- ggEAR.SC + xlim(0,6e6)
ggEAR.SC <- ggEAR.SC + ylab("Axial Strain Rates: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
ggEAR.SC

# ---- measured axial strain only -----
ggEAR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME, y=EAR, color = ITEST))
ggEAR.SC <- ggEAR.SC + geom_point()
ggEAR.SC <- ggEAR.SC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggEAR.SC <- ggEAR.SC + xlim(0,6e6)# + ylim(-1.5e-5,5.0e-6)
ggEAR.SC <- ggEAR.SC + ylab("Axial Strain Rates [/sec]: Measured (dot)") + xlab("Time [sec]")
ggEAR.SC

# ---- plot axial strain rates ----
ggFEAR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME, y=FEAR, color = ITEST))
ggFEAR.SC <- ggFEAR.SC + geom_line()
ggFEAR.SC <- ggFEAR.SC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggFEAR.SC <- ggFEAR.SC + xlim(0,6e6)# + ylim(-1.5e-5,5.0e-6)
ggFEAR.SC <- ggFEAR.SC + ylab("Axial Strain Rates [/sec]: Calculated (line)") + xlab("Time [sec]")
ggFEAR.SC

#===================================================
# ---- plot lateral strain rates ----
ggELR.SC <- ggplot(data = PLOT.DATA, aes(x=TIME, y=ELR, color = ITEST))
ggELR.SC <- ggELR.SC + geom_point()
ggELR.SC <- ggELR.SC + geom_line(aes(y=FELR))+ facet_wrap(~ITEST, ncol=3, scales = "free")
ggELR.SC <- ggELR.SC + facet_wrap(~ITEST, ncol=3, scales = "free")
ggELR.SC <- ggELR.SC + xlim(0,6e6)
ggELR.SC <- ggELR.SC + ylab("Lateral Strain Rates: Calculated (line) Vs. Measured (dot)") + xlab("Time [sec]")
ggELR.SC