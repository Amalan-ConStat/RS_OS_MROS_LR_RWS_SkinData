library(here)

# OSMAC Data ----
Base_Path<-"OSMAC"

## No Correlated Covariate Data -----
load(here("Identical_r0","Generate_Big_Data","Scaled.RData"))
save(list = c("Original_Data","Replicates","Subsample_Size","r0","All_Models","Choice"),
     file = here("Identical_r0","Simulation_Setup","Classical","Analysis","OSMAC","Init.RData"))        

remove(Base_Path)

rm(list = ls())
