setwd("/projet/gabi/work/ahulot/GREMLINGLASSO/2020_04_02/")
library(knitr)
library(rmarkdown)
library(parallel)

LibPath = "/projet/gabi/save/ahulot/R_lib/"

knit("2020_04_02.Rmd", tangle = TRUE)
render("2020_04_02.Rmd", output_format = "html_document", params = list(LibPath = "/projet/gabi/save/ahulot/R_lib/"))

quit("no")
