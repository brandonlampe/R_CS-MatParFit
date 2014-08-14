#==== routine for loading consolidation data (of '.csv' format) into R

#---- identify file path where documents are located ---
f <- list.files(path =
  "/Users/Lampe/GrantNo456417/MatParameterFitting/CompiledTestData/C-S_Database/CSV_M_SC",
  all.files = FALSE, full.names = TRUE)

#---- summary of the contents paths loaded into 'f' ---
print(str(f))

#---- create list containing data frames, with the following dim: [[16]][100,22]
M <- length(f) # number of data sets loaded
DATA.LIST <- list()
DATA.LIST <- lapply(f,function(i){
  read.csv(i,header = FALSE, sep = ",", dec = ".",stringsAsFactors = FALSE)})

#---- bind all data frames into one data frame that contains all measured data ----
DATA.MEAS <- data.frame(DATA.LIST[[1]])
for (i in 2:M){
	DATA.MEAS = rbind(DATA.MEAS,DATA.LIST[[i]])
	}

# print(summary(Data.all))
