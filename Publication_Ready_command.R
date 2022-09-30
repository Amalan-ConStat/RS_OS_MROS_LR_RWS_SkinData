library(here)
library(ezknitr)
library(rmarkdown)

# Publication Ready ----
ezknit(file=here("Publication_Ready.Rmd"),
       out_dir=here("Publication_Ready"),
       fig_dir = c("Plots"),
       verbose = TRUE,keep_md = TRUE,keep_html = FALSE)
open_output_dir()
