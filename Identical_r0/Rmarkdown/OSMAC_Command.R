library(here)
library(rmarkdown)

# OSMAC Method ----
render(input =here("Identical_r0","Rmarkdown","OSMAC_Method.Rmd"),
       output_file = "OSMAC_Method",
       output_format = "html_document",
       output_dir=here("Identical_r0","htmloutputs","OSMAC"))
