# Run for OSMAC Sampling----
library(compiler)
library(here)
library(LaplacesDemon)
library(matrixStats)

enableJIT(1)

# Load the Scaled Data----
load(here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Init.RData"))
#load(here("Init.RData"))

# Load the OSMAC Sample----
load(here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Run_OSMAC.RData"))
#load(here("Run_OSMAC.RData"))

# Generate for Random sample of 1000 different times ---
## Pi=MMSE and Pi=MVC Proportional probabilities ---
Final_Parameter<-Run_OSMAC(Replicates = Replicates, r1 = r0, r2 = Choice,
                           Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = nrow(Original_Data),
                           alpha = rep(1/length(All_Models),length(All_Models)),
                           All_Combinations = All_Models,
                           All_Covariates = All_Models[[length(All_Models)]])

