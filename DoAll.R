#Calls all the routines necessary for fitting shear parameters

# ---- clear memory ----
rm(list = ls())
#==== set working directory and load libraries ====
# #store current directory
# Initial.dir <- getwd()
#
# #change to new directory
CurrentDirectory <- getwd()
# CurrentDirectory <- "~/R/WorkingDirectory/MatFit_v1-2"
# setwd(CurrentDirectory)

# load libraries
library("minpack.lm")
library("ggplot2")
library("plyr")
library("data.table")
library("pracma")
library("binhf")
library("foreach")
library("iterators")
library("deSolve")
library("FME")
library(grid)
library(gridExtra)
library(fBasics)

#--------------------------

source("Load.R")
source("DefPara.R")
#source("Plot_Input.R")
#source("ShearFit.R")
#source("ValidateSC.R")
#source("Plot_Output.R")

