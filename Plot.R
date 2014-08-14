#==== plotting data ====

#==== set working directory and load libraries ====
#store current directory
Initial.dir <- getwd()

#change to new directory
CurrentDirectory <- "~/R/WorkingDirectory/MatFit_v1-0/Do.R"
setwd(CurrentDirectory)

# load necessary libraries
library("ggplot2")
library("scales")

#---- load data set ----
source("Load.R")


#---- label/define parameters for plotting ----
source("DefPara.R")
#===== plot data ===============================

#---- load data to be plotted into a single data frame ----
PlotData.df <- data.frame(ID,T,LS,AS,DD,MS,DS,RHO,EV,EA,EL,EAR,ELR,RAT)

qplot(T,EA, data = PlotData.df, colour = ID, geom = c("point", "line"))
qplot(T,EA, data = PlotData.df, colour = ID, geom = c("point", "line"),facets = ID~.)

plot(T,EA)

gg <- ggplot(data = PlotData.df, aes(x=T, y=EA, color = ID)) + geom_point() + geom_line(aes(y=EL))

gg1 <- ggplot(data = PlotData.df, aes(x=T, y=EA, color = ID)) + geom_point()
gg1 <- gg1 + geom_line(aes(y=EL))+ facet_wrap(~ID, ncol=3, scales = "free")
gg1 <- gg1 + xlim(0,6e6)
gg1 <- gg1 + ylab("Axial and Lateral Strains") + xlab("Time [sec]")
gg1

#ggplot(PlotData.df, aes(x=T, y=EA, colour = ID)) + geom_point() + geom_line(aes(y=EL)) + facet_wrap(~ID, ncol=3, scales = "free")


p1 <- ggplot(Shear.df, aes(x=RHO, y=RAT, colour = ID)) + geom_point() + geom_line(aes(y=NewFit)) +xlim(1500,2200) + facet_wrap(~ID, ncol=3, scales = "free")

#


