library(here)
library(rmarkdown)

# Best_Subsampling_Method ----
render(input=here("Identical_r0","Rmarkdown","Summary","Best_Subsampling_Method.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Identical_r0","Summary","Best_Subsampling"))
