#!/usr/bin/Rscript

source("https://bioconductor.org/biocLite.R") ;
biocLite("ballgown", ask=FALSE, dependencies = TRUE) ;
install.packages("sf", repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
install.packages(c("DataCombine", "glue", "sunburstR", "docopt"), repos = "https://cloud.r-project.org/", dependencies = TRUE)
