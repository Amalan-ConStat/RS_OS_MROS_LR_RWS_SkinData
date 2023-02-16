# Run for Random Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Scaled Data----
load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","Random_Sampling","Init.RData"))

# Load the Random Sample----
load(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","Random_Sampling","Run_RandomSample.RData"))

# All Models ----
for (i in 1:length(All_Models)) 
{
  print(All_Models[[i]])
  # Generate for Random sample of 1000 different times
  Final_Parameter<-Run_RandomSample(Replicates = Replicates, 
                                    FullData = Original_Data[,colnames(Original_Data) %in% c("Y",All_Models[[i]])],
                                    N=nrow(Original_Data),
                                    Subsample_Size=Subsample_Size,Choices=Choice)
  
  Est_Param_RandomSample<-Final_Parameter$EstPar
  Utility_RandomSample<-Final_Parameter$Utility
  SelectedData_RandomSample<-Final_Parameter$Subsampled_Data
  Bias_RandomSample<-Final_Parameter$Bias
  
  # Save the Results ----
  save(Est_Param_RandomSample,Utility_RandomSample,SelectedData_RandomSample,Bias_RandomSample,
       file = here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","Random_Sampling",
                   "Results",paste0("Random_Sample_output_",i,".RData")))
  
  save(Est_Param_RandomSample,Utility_RandomSample,SelectedData_RandomSample,Bias_RandomSample,
       file = here("Non_Identical_r0","Outputs","Classical","Random_Sampling",
                   paste0("Random_Sample_output_",i,".RData")))
}

