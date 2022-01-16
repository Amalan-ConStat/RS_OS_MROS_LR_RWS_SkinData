library(here)
library(ezknitr)
library(rmarkdown)

# OSMAC Method ----
ezknit(file=here("Identical_r0","Rmarkdown","Classical","OSMAC_Method.Rmd"),
       out_dir=here("Identical_r0","htmloutputs","Classical","OSMAC"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = FALSE)
open_output_dir()
