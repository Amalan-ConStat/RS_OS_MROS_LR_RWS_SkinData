library(here)
library(rmarkdown)

# Publication Ready ----
render(input=here("Publication_Ready.Rmd"),
       output_file = "Publication_Ready",
       output_format = "github_document",
       output_dir=here("Publication_Ready"))
