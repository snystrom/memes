# set the cran mirror
options("repos"="http://cran.rstudio.com") 

packages <- c("tidyverse",
              "devtools",
              "roxygen2",
              "testthat",
              "knitr", 
							"rmarkdown",
              "BiocManager")

packages <- setdiff(packages, installed.packages()[, "Package"])

if (length(packages) != 0){
  (install.packages(packages))
}

# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

update.packages(ask=FALSE)
