library(here)
library(compiler)
library(LaplacesDemon)

# Cordeiro Bias Estimation ----
Cordeiro<-function(XData,With_bias)
{
  p <- as.vector(invlogit(XData%*%as.vector(With_bias)))
  W <- diag(p*(1-p))
  inverse_term <- solve(t(XData)%*%W%*%XData)

  Term1 <- inverse_term%*%t(XData)%*%W
  Term2 <- diag(diag(XData%*%(inverse_term)%*%t(XData))) %*% as.matrix(p-0.5)

  bias <- as.vector(Term1%*%Term2)
  return(bias)
}

# Using Random Sampling Sub-sample from Big Data ----
Run_RandomSample <- function(Replicates,FullData,Subsample_Size,N,Choices)
{
  #RandomSample<-list()
  #Sample_RS<-list()
  
  All_Parameter<-list()
  Est_Parameter <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))

  All_Bias<-list()
  Bias <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))

  All_optimality<-list()
  optimality <- matrix(NA,nrow = length(Choices),ncol = 3)

  for(i in 1:Replicates)
  {
    Sampled<-sample(1:N,size = Subsample_Size)
    Temp_Data<-as.data.frame(FullData[Sampled,])
    #RandomSample[[i]]<-Temp_Data

    for (j in 1:length(Choices))
    {
      Temp_Temp_Data<-Temp_Data[1:Choices[j],]
      Results<-glm(Y~.-1,data=Temp_Temp_Data,family = binomial)
      Est_Parameter[j,]<-c(Choices[j],Results$coefficients)
      optimality[j,]<-c(Choices[j], psych::tr(vcov(Results)), det(solve(summary(Results)$cov.unscaled)) )
      Bias[j,]<-c(Choices[j],Cordeiro(XData = as.matrix(Temp_Temp_Data[,-1]),
                                      With_bias = as.vector(Results$coefficients)))
    }

    All_Parameter[[i]]<-Est_Parameter
    All_optimality[[i]]<-optimality
    All_Bias[[i]]<-Bias
    #Sample_RS[[i]]<-cbind(i,RandomSample[[i]])
    print(i)
  }

  Final_Parameter<-do.call(rbind,All_Parameter)
  Final_optimality<-do.call(rbind,All_optimality)
  Final_Bias<-do.call(rbind,All_Bias)
  #Final_Sample_RS<-do.call(rbind,Sample_RS)
  
  return(list("EstPar"=Final_Parameter,
              "Utility"=Final_optimality,
              "Bias"=Final_Bias#,
              #"Subsampled_Data"=Final_Sample_RS
              )
         )
}

Run_RandomSample<-cmpfun(Run_RandomSample)

#Save the Random Sample function for Scaled Data ----
save(Run_RandomSample,Cordeiro,
     file=here("Non_Identical_r0","Simulation_Setup","Classical","Analysis",
               "Random_Sampling","Run_RandomSample.RData"))

rm(list = ls())

# Run the Random sampling method ----
source(here("Non_Identical_r0","Simulation_Setup","Classical","Analysis","Random_Sampling",
            "Simulation_Results_Random_Sampling.R"))
