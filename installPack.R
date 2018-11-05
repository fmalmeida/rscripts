#!/usr/bin/Rscript

source("https://bioconductor.org/biocLite.R") ;
biocLite("BiocUpgrade") ; 
biocLite("Rsamtools", ask=FALSE, dependencies = TRUE);
install.packages("sf", repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
install.packages(c("DataCombine", "glue", "sunburstR", "docopt", "dplyr"), repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
biocLite("ballgown", ask=FALSE, dependencies = TRUE) ;
