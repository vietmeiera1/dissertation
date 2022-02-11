
#check version of R. 
R.version

#install new version of R here:
# https://cran.revolutionanalytics.com/ 

#Install package
install.packages("installr")

#load package
library(installr)

#package help
? `installr-package`

#updateR function
updateR()


#Packages to Install after update

install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dada2")
#not available for R version 4, 3, 

install.packages("installr")
install.packages("phyloseq")

source("https://bioconductor.org/biocLite.R")
biocLite("dada2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.4")




