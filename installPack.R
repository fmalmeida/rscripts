#!/usr/bin/Rscript

source("https://bioconductor.org/biocLite.R") ;
#biocLite("Gviz", ask=FALSE, dependencies = TRUE) ;
biocLite("ballgown", ask=FALSE, dependencies = TRUE) ;
biocLite("ggbio", ask=FALSE, dependencies = TRUE) ;
install.packages("sf", repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
#install.packages(c("DataCombine", "optparse", "plotly", "plyr", "sunburstR"), repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages(c("DataCombine", "glue", "sunburstR"), repos = "https://cloud.r-project.org/", dependencies = TRUE)
