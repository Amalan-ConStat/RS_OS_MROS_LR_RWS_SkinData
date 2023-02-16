library(here)

# Random Sampling Data ----
Base_Path<-"Random_Sampling"

## No Correlated Covariate Data -----
load(here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))
save(list = c("Original_Data","Replicates","Subsample_Size","All_Models","Choice"),
     file = here("Non_Identical_r0","Simulation_Setup","Analysis",
                 Base_Path,"Init.RData"))        

# Rare Event Random Sampling Data ----
Base_Path<-"RE_Random_Sampling"

## No Correlated Covariate Data -----
load(here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))
save(list = c("Original_Data","Replicates","Subsample_Size","All_Models","Choice"),
     file = here("Non_Identical_r0","Simulation_Setup","Analysis",
                 Base_Path,"Init.RData"))        

# OSMAC Data ----
Base_Path<-"OSMAC"

## No Correlated Covariate Data -----
load(here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))
save(list = c("Original_Data","Replicates","Subsample_Size","r0","All_Models","Choice"),
     file = here("Non_Identical_r0","Simulation_Setup","Analysis",
                 Base_Path,"Init.RData"))        

remove(Base_Path)

# OSMAC Model Free Data ----
Base_Path<-"OSMAC_Model_Free"

## No Correlated Covariate Data -----
load(here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))
save(list = c("Original_Data","Replicates","Subsample_Size","r0","All_Models","Choice"),
     file = here("Non_Identical_r0","Simulation_Setup","Analysis",
                 Base_Path,"Init.RData"))

remove(Base_Path)

rm(list = ls())
