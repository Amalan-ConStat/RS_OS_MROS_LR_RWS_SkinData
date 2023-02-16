library(here)
library(rmarkdown)

# Best_Subsampling_Method ----
render(input=here("Non_Identical_r0","Rmarkdown","Classical","Summary","Best_Subsampling_Method.Rmd"),
       output_file = "Best_Subsampling_Method",
       output_format = "html_document",
       output_dir=here("Non_Identical_r0","Summary","Classical","Best_Subsampling"))
