# BioNetStat

##BioNetStat

Package to perform network analysis. BioNetStat is able to compare two or more correlation networks.

## Installation
### Dependencies
R (>= 3.0.0), shiny (>= 0.8.0), WGCNA, igraph, shinyBS, RColorBrewer, Hmisc, pathway, psych, RJSONIO, whisker, yaml, pheatmap, preprocessCore, GO.db, AnnotationDbi, impute, and ggplot2

### Installation steps

1. If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
2. Open the R command line interface, and install all BNS dependencies (if they have not been installed yet):
```
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
biocLite("GO.db")
biocLite("AnnotationDbi")
biocLite("pathview")
install.packages(c("WGCNA", "igraph", "RColorBrewer", "Hmisc", "psych", "RJSONIO", "whisker", "yaml", "pheatmap", "ggplot2","devtools","shiny")) 
```

2. a. If pathview don't install try to ilstall this libraries in linux shell terminal.
```
$ sudo apt-get install libxml2-dev
$ sudo apt-get install libcurl4-openssl-dev
$ sudo apt-get install libssl-dev
```
2. b. And, try to install 'pathview' again
```
source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
```
3. Please, install version 0.20 for shinyBS. We are working to make the BioNetStat package compatible with the new versions of the packages as soon as possible. To install the recommended versions for shinyBS, just type the following commands on the R command-line:
```
devtools::install_version("shinyBS", "0.20")
```
4. Install the BioNetStat (BNS) package
```
library(devtools)
install_github("jardimViniciusC/BioNetStat")
```
5. To run BNS, just type the following code:
```
library(BioNetStat)
library(shiny)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
Wait for the browser page to open, and enjoy BoiNetStat!

## How to use BioNetStat
The tutorials with examples data sets are in 'Tutorial para uso do BioNetStat em linhas de comando.md' or in 'Tutorial BioNetStat para a interface.md'.
