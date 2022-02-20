library(here)
library(compiler)
library(LaplacesDemon)

# Load the OSMAC algorithm
source(here("Identical_r0","Simulation_Setup","Classical","R_Scripts","OSMAC_Algorithm.R"))

AlgTwoStp<-cmpfun(AlgTwoStp)
getMLE<-cmpfun(getMLE)
Cordeiro<-cmpfun(Cordeiro)

# Using OSMAC to Sub-sample from Big Data ----
Run_OSMAC<-function(Replicates,r1,r2,Y,X,N,alpha,All_Combinations,All_Covariates)
{
  # From r0 and r conduct OSMAC subsampling for Logistic regression-----
  Parameter_Single_mMSE<-Parameter_Single_mVc<-replicate(length(All_Combinations),list())
  Parameter_ModelFree_mMSE<-Parameter_ModelFree_mVc<-replicate(length(All_Combinations),list())
  
  Bias_Single_mMSE<-Bias_Single_mVc<-replicate(length(All_Combinations),list())
  Bias_ModelFree_mMSE<-Bias_ModelFree_mVc<-replicate(length(All_Combinations),list())
  
  Utility_Single_mMSE<-Utility_Single_mVc<-replicate(length(All_Combinations),list())
  Utility_ModelFree_mMSE<-Utility_ModelFree_mVc<-replicate(length(All_Combinations),list())
  
  # Sample_Single_mMSE<-Sample_Single_mVc<-replicate(length(All_Combinations),list())
  # Sample_ModelFree_mMSE<-Sample_ModelFree_mVc<-replicate(length(All_Combinations),list())
  
  # Full_SP_Single<-replicate(length(All_Combinations),list())
  # Full_SP_ModelFree<-replicate(length(All_Combinations),list())
  
  i=1
  while(i <= Replicates)
  {
    tryCatch({
      Results<-AlgTwoStp(r1,r2,Y,X,n=N,alpha,All_Combinations,All_Covariates)
      
      for (j in 1:length(All_Combinations)) 
      {
        Parameter_Single_mMSE[[j]][[i]]<-Results$opt$opt_Single[[j]]$Est_Theta_mMSE
        Parameter_Single_mVc[[j]][[i]]<-Results$opt$opt_Single[[j]]$Est_Theta_mVc
        Parameter_ModelFree_mMSE[[j]][[i]]<-Results$opt$opt_MF[[j]]$Est_Theta_mMSE
        Parameter_ModelFree_mVc[[j]][[i]]<-Results$opt$opt_MF[[j]]$Est_Theta_mVc
        
        Bias_Single_mMSE[[j]][[i]]<-Results$opt$opt_Single[[j]]$Bias_mMSE
        Bias_Single_mVc[[j]][[i]]<-Results$opt$opt_Single[[j]]$Bias_mVc
        Bias_ModelFree_mMSE[[j]][[i]]<-Results$opt$opt_MF[[j]]$Bias_mMSE
        Bias_ModelFree_mVc[[j]][[i]]<-Results$opt$opt_MF[[j]]$Bias_mVc
        
        Utility_Single_mMSE[[j]][[i]]<-Results$opt$opt_Single[[j]]$Utility_mMSE
        Utility_Single_mVc[[j]][[i]]<-Results$opt$opt_Single[[j]]$Utility_mVc
        Utility_ModelFree_mMSE[[j]][[i]]<-Results$opt$opt_MF[[j]]$Utility_mMSE
        Utility_ModelFree_mVc[[j]][[i]]<-Results$opt$opt_MF[[j]]$Utility_mVc
      }
      print(i)
      i<-i+1
    },error=function(e){ i=i; print("Failed")})
  }

  for (j in 1:length(All_Combinations)) 
  {
    # Single Model
    Final_param_mMSE<-do.call(rbind,Parameter_Single_mMSE[[j]])
    Final_param_mVc<-do.call(rbind,Parameter_Single_mVc[[j]])
    
    Final_bias_mMSE<-do.call(rbind,Bias_Single_mMSE[[j]])
    Final_bias_mVc<-do.call(rbind,Bias_Single_mVc[[j]])
    
    Final_utility_mMSE<-do.call(rbind,Utility_Single_mMSE[[j]])
    Final_utility_mVc<-do.call(rbind,Utility_Single_mVc[[j]])
    
    Results_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_param_mMSE),
                        "mVc_Output"=cbind.data.frame("mVc",Final_param_mVc))
    colnames(Results_OSMAC$mMSE_Output)<-
      colnames(Results_OSMAC$mVc_Output)<-c("Type","Subsample_Size",paste0("Theta",0:(length(All_Combinations[[j]])-1) ))
    
    Bias_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_bias_mMSE),
                     "mVc_Output"=cbind.data.frame("mVc",Final_bias_mVc))
    colnames(Bias_OSMAC$mMSE_Output)<-
      colnames(Bias_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                         paste0("Theta",0:(length(All_Combinations[[j]])-1)))
    
    Utility_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_utility_mMSE),
                        "mVc_Output"=cbind.data.frame("mVc",Final_utility_mVc))
    colnames(Utility_OSMAC$mMSE_Output)<-colnames(Utility_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                                                               "A_optimality","D_optimality")
    
    save(list= c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"),
         file = here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC",
                     "Results","Single_Model",paste0("OSMAC_output_",j,".RData")))
    
    save(list = c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"),
         file = here("Identical_r0","Outputs","Classical","OSMAC","Single_Model",
                    paste0("OSMAC_output_",j,".RData")))
    
    # Model Free
    Final_param_mMSE<-do.call(rbind,Parameter_ModelFree_mMSE[[j]])
    Final_param_mVc<-do.call(rbind,Parameter_ModelFree_mVc[[j]])
    
    Final_bias_mMSE<-do.call(rbind,Bias_ModelFree_mMSE[[j]])
    Final_bias_mVc<-do.call(rbind,Bias_ModelFree_mVc[[j]])
    
    Final_utility_mMSE<-do.call(rbind,Utility_ModelFree_mMSE[[j]])
    Final_utility_mVc<-do.call(rbind,Utility_ModelFree_mVc[[j]])
    
    Results_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_param_mMSE),
                        "mVc_Output"=cbind.data.frame("mVc",Final_param_mVc))
    colnames(Results_OSMAC$mMSE_Output)<-
      colnames(Results_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                            paste0("Theta",0:(length(All_Combinations[[j]])-1)))
    
    Bias_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_bias_mMSE),
                     "mVc_Output"=cbind.data.frame("mVc",Final_bias_mVc))
    colnames(Bias_OSMAC$mMSE_Output)<-
      colnames(Bias_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                         paste0("Theta",0:(length(All_Combinations[[j]])-1)))
    
    Utility_OSMAC<-list("mMSE_Output"=cbind.data.frame("mMSE",Final_utility_mMSE),
                        "mVc_Output"=cbind.data.frame("mVc",Final_utility_mVc))
    colnames(Utility_OSMAC$mMSE_Output)<-colnames(Utility_OSMAC$mVc_Output)<-c("Type","Subsample_Size",
                                                                               "A_optimality","D_optimality")
    
    save(list= c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"),
         file = here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC",
                     "Results","ModelFree",paste0("OSMAC_output_",j,".RData")))
    
    save(list = c("Results_OSMAC","Bias_OSMAC","Utility_OSMAC"),
         file = here("Identical_r0","Outputs","Classical","OSMAC","ModelFree",
                    paste0("OSMAC_output_",j,".RData")))
  }
  
}

Run_OSMAC<-cmpfun(Run_OSMAC)

#Save the OSMAC Sample function----
save(list =c("Run_OSMAC","AlgTwoStp","getMLE","Cordeiro"),
     file=here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Run_OSMAC.RData"))

rm(list = ls())

# Run the OSMAC sampling method ----
source(here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC",
            "Simulation_Results_OSMAC_Sampling.R"))
