library(here)
library(ezknitr)
library(rmarkdown)

# Random Sampling ----
ezknit(file=here("Non_Identical_r0","Rmarkdown","Classical","Random_Sampling.Rmd"),
       out_dir=here("Non_Identical_r0","htmloutputs","Classical","Random_Sampling"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = FALSE)
open_output_dir()

# Rare Event Random Sampling ----
ezknit(file=here("Non_Identical_r0","Rmarkdown","Classical","RE_Random_Sampling.Rmd"),
       out_dir=here("Non_Identical_r0","htmloutputs","Classical","RE_Random_Sampling"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = FALSE)
open_output_dir()

# OSMAC Method ----
ezknit(file=here("Non_Identical_r0","Rmarkdown","Classical","OSMAC_Method.Rmd"),
       out_dir=here("Non_Identical_r0","htmloutputs","Classical","OSMAC"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = FALSE)
open_output_dir()

# OSMAC Model Free Method ----
ezknit(file=here("Non_Identical_r0","Rmarkdown","Classical","OSMAC_Model_Free_Method.Rmd"),
       out_dir=here("Non_Identical_r0","htmloutputs","Classical","OSMAC_Model_Free"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = FALSE)
open_output_dir()
