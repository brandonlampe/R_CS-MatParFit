
#Calls all the routines necessary for fitting shear parameters

#==== set working directory and load libraries ====
#store current directory
Initial.dir <- getwd()

#change to new directory
CurrentDirectory <- "~/R/WorkingDirectory/MatFit_v1-2"
setwd(CurrentDirectory)

# load necessary libraries
library("minpack.lm")
library("ggplot2")
library("plyr")
library("data.table")
library("pracma")

#--------------------------

source("Load.R")
source("DefPara.R")
#source("Plot_Input.R")
#source("ShearFit.R")
#source("ValidateSC.R")
#source("Plot_Output.R")

# clear workspace
#rm(list = ls())