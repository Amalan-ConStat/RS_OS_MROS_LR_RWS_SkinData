# Run for OSMAC Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the Scaled Data----
#load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Init.RData"))
load(here("Init.RData"))

# Load the OSMAC Sample----
#load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Run_OSMAC.RData"))
load(here("Run_OSMAC.RData"))

# All Models ----
for (Iterator in 1:length(All_Models)) 
{
  print(All_Models[[Iterator]])
  # Generate for Random sample of 1000 different times ---
  ## Pi=MMSE and Pi=MVC Proportional probabilities ---
  Final_Parameter<-Run_OSMAC(Replicates=Replicates,r1=r0,r2=Choice,
                             Y=as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                             X=as.matrix(Original_Data[,colnames(Original_Data) %in% All_Models[[Iterator]]]),
                             Real_Data=as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y",All_Models[[Iterator]])]),
                             N=nrow(Original_Data),Model="Real_Model",
                             Iterator=Iterator)
}

