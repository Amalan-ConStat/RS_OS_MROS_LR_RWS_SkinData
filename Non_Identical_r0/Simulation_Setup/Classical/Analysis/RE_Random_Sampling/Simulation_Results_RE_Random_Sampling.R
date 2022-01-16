# Run for Rare Event Random Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the Scaled Data----
load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","RE_Random_Sampling","Init.RData"))

# Load the Rare Event Random Sample----
load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","RE_Random_Sampling","RE_RandomSample.RData"))

# All Models ----
for (i in 1:length(All_Models)) 
{
  print(All_Models[[i]])
  # Generate for Rare Evnet Random sample of 1000 different times
  Final_Parameter<-Run_RE_RandomSample(Replicates = Replicates, 
                                       FullData = Original_Data[,colnames(Original_Data) %in% c("Y",All_Models[[i]])],
                                       N=nrow(Original_Data),
                                       Subsample_Size=Subsample_Size,Choices=Choice)
  
  Est_Param_RE_RandomSample<-Final_Parameter$EstPar
  Utility_RE_RandomSample<-Final_Parameter$Utility
  Bias_RE_RandomSample<-Final_Parameter$Bias
  SelectedData_RE_RandomSample<-Final_Parameter$Subsampled_Data
  
  # Save the Results ---
  save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,SelectedData_RE_RandomSample,
       file = here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","RE_Random_Sampling",
                   "Results",paste0("RE_Random_Sample_output_",i,".RData")))
  
  save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,SelectedData_RE_RandomSample,
       file = here("Non_Identical_r0","Outputs","Classical","RE_Random_Sampling",
                   paste0("RE_Random_Sample_output_",i,".RData")))
}

