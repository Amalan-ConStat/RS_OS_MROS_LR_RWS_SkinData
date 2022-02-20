library(here)
library(compiler)
library(LaplacesDemon)

# Load the OSMAC algorithm
source(here("Non_Identical_r0","Simulation_Setup","Classical","R_Scripts","OSMAC_Algorithm.R"))

AlgTwoStp<-cmpfun(AlgTwoStp)
getMLE<-cmpfun(getMLE)
Cordeiro<-cmpfun(Cordeiro)

# Using OSMAC to Sub-sample from Big Data ----
Run_OSMAC<-function(Replicates,r1,r2,Y,X,Real_Data,N,Model,Iterator)
{
  # From r0 and r conduct OSMAC subsampling for Logistic regression-----
  Parameter_mMSE<-list()
  Parameter_mVc<-list()

  Bias_mMSE<-list()
  Bias_mVc<-list()

  Utility_mMSE<-list()
  Utility_mVc<-list()

  #Sample_mMSE<-list()
  #Sample_mVc<-list()
  
  i=1
  while(i <= Replicates)
  {
    tryCatch({
      Results<-AlgTwoStp(r1,r2,Y,X,Real_Data,n=N,Model)
      
      Parameter_mMSE[[i]]<-Results$opt$Est_Theta_mMSE
      Parameter_mVc[[i]]<-Results$opt$Est_Theta_mVc
      
      Bias_mMSE[[i]]<-Results$opt$Bias_mMSE
      Bias_mVc[[i]]<-Results$opt$Bias_mVc
      
      Utility_mMSE[[i]]<-Results$opt$Utility_mMSE
      Utility_mVc[[i]]<-Results$opt$Utility_mVc
      
      #Sample_mMSE[[i]]<-cbind(i,Results$opt$Sample_mMSE)
      #Sample_mVc[[i]]<-cbind(i,Results$opt$Sample_mVc)
      
      print(i)
      i<-i+1
    },error=function(e){ i=i})
    
  }

  Final_param_mMSE<-do.call(rbind,Parameter_mMSE)
  Final_param_mVc<-do.call(rbind,Parameter_mVc)

  Final_bias_mMSE<-do.call(rbind,Bias_mMSE)
  Final_bias_mVc<-do.call(rbind,Bias_mVc)

  Final_utility_mMSE<-do.call(rbind,Utility_mMSE)
  Final_utility_mVc<-do.call(rbind,Utility_mVc)
  
  #Final_Sample_mMSE<-do.call(rbind,Sample_mMSE)
  #Final_Sample_mVc<-do.call(rbind,Sample_mVc)
    
  Results_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_param_mMSE),
                      "mVc_Output"=cbind.data.frame("mVc",Final_param_mVc))
  colnames(Results_OSMAC$mMSE_Output)<-colnames(Results_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                                                             paste0("Theta",0:(ncol(Real_Data[,-1])-1)))

  Bias_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_bias_mMSE),
                   "mVc_Output"=cbind.data.frame("mVc",Final_bias_mVc))
  colnames(Bias_OSMAC$mMSE_Output)<-colnames(Bias_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                                                       paste0("Theta",0:(ncol(Real_Data[,-1])-1)))

  Utility_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_utility_mMSE),
                      "mVc_Output"=cbind.data.frame("mVc",Final_utility_mVc))
  colnames(Utility_OSMAC$mMSE_Output)<-colnames(Utility_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                                                             "A_optimality","D_optimality")
  
  # Sample_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_Sample_mMSE),
  #                    "mVc_Output"=cbind.data.frame("mVc",Final_Sample_mVc))
  # colnames(Sample_OSMAC$mMSE_Output)<-colnames(Sample_OSMAC$mVc_Output)<-c("Type","Simulation",
  #                                                                          "Subsample_Size","Y","SP",
  #                                                                          paste0("X",0:(ncol(Real_Data[,-1])-1)))
  
  save(list= c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"#,"Sample_OSMAC"
               ),
       file = here(#"Non_Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC",
                   "Results",#Model,
                   paste0("OSMAC_output_",Iterator,".RData")))

  # save(list = c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"#,"Sample_OSMAC"
  #               ),
  #      file = here("Non_Identical_r0","Outputs","Classical","OSMAC",#Model,
  #                  paste0("OSMAC_output_",Iterator,".RData")))
}

Run_OSMAC<-cmpfun(Run_OSMAC)

#Save the OSMAC Sample function----
save(list =c("Run_OSMAC","AlgTwoStp","getMLE","Cordeiro"),
     file=here("Non_Identical_r0","Simulation_Setup","Classical","Analysis",
               "OSMAC","Run_OSMAC.RData"))

rm(list = ls())

# Run the OSMAC sampling method ----
source(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC",
            "Simulation_Results_OSMAC_Sampling.R"))
