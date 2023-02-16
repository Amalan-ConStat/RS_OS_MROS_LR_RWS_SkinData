library(here)
library(ezknitr)
library(rmarkdown)

# Best_Subsampling_Method ----
ezknit(file=here("Non_Identical_r0","Rmarkdown","Classical","Summary","Best_Subsampling_Method.Rmd"),
        out_dir=here("Non_Identical_r0","Summary","Classical","Best_Subsampling"),
        fig_dir = c("Plots"),
        verbose = TRUE,keep_md = FALSE)
open_output_dir()
