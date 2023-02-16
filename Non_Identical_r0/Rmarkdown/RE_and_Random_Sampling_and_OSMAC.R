library(here)
library(rmarkdown)

# Random Sampling ----
render(input=here("Non_Identical_r0","Rmarkdown","Random_Sampling.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Non_Identical_r0","htmloutputs","Random_Sampling"))

# Rare Event Random Sampling ----
render(input=here("Non_Identical_r0","Rmarkdown","RE_Random_Sampling.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Non_Identical_r0","htmloutputs","RE_Random_Sampling"))

# OSMAC Method ----
render(input=here("Non_Identical_r0","Rmarkdown","OSMAC_Method.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Non_Identical_r0","htmloutputs","OSMAC"))

# OSMAC Model Free Method ----
render(input=here("Non_Identical_r0","Rmarkdown","OSMAC_Model_Free_Method.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Non_Identical_r0","htmloutputs","OSMAC_Model_Free"))
