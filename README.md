# BioNetStat

## Package to perform network analysis. 
### BioNetStat is able to compare two or more correlation networks.

The diversity of interactions that occurs in biological systems from the cells and organelles to the whole biosphere, can be assessed with tools of the networks theory. The dynamic in structure and in the interactions among systems elements is an inherent trait of those systems. Several tools have been proposed to compare networks, representing the many states assumed by a system. However, until the present, none of them is able to compare structural characteristics among more than two networks simultaneously. Due to the many states that can be assumed by a given biological system, we developed a statistical tool to compare two or more networks and point key variables in a system. BioNetStat is able to compare correlation networks using traits that are based on graph spectra (the group of eigenvalues in an adjacency matrix), such as the spectral distribution. This measure is associated with several structural characteristics of networks such as the number of walks, diameter, and clicks. In addition to the spectral distribution, BioNetStat can also compare networks by using spectral entropy, degree distribution, and nodes centralities. Until now the BioNetStat theoretical base is available only in the [master dissertation](http://www.teses.usp.br/autor.php?autor=37A71EBFAC13) in Portuguese. The paper is under production and will be posted here as soon as possible.

Authors: Jardim, V., Santos, S., Fujita, A., Buckeridge, S. (2018). BioNetStat: A tool for biological network analysis.

Maintainer: Laboratório de fisiologia ecologica de plantas [LAFIECO](http://www.lafieco.com.br/)

## Installation
### Dependencies
R (>= 3.4), shiny (>= 0.8.0), igraph, shinyBS, pathview

### Installation steps

1. If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
2. Open the R command line interface, and install all BNS dependencies (if they have not been installed yet):
```Rscript
install.packages(c("shiny","WGCNA", "igraph", "RColorBrewer", "Hmisc", "psych", "RJSONIO", "whisker", "yaml", "pheatmap", "ggplot2","devtools","matrix","mgcv"))
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
biocLite("GO.db")
biocLite("AnnotationDbi")
biocLite("pathview")
biocLite ("Biostrings")
```

2. a. If 'pathview' package don't install, try to install this libraries in linux shell terminal.
```Rscript
$ sudo apt-get install libxml2-dev
$ sudo apt-get install libcurl4-openssl-dev
$ sudo apt-get install libssl-dev
```
2. b. And, try to install 'pathview' again
```Rscript
source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
```
3. Please, install version 0.20 for shinyBS. We are working to make the BioNetStat package compatible with the new versions of the packages as soon as possible. To install the recommended versions for shinyBS, just type the following commands on the R command-line:
```Rscript
devtools::install_version("shinyBS", "0.20")
```
4. Install the BioNetStat (BNS) package
```Rscript
library(devtools)
install_github("jardimViniciusC/BioNetStat")
```
## Running BioNetStat

After installed and the next time you want to use BioNetStat, to run, just type the following code:
```Rscript
library(BioNetStat)
library(shiny)
```
If you want to use a Graphical interface type
```Rscript
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
Wait for the browser page to open, and enjoy BioNetStat!

## How to use BioNetStat?
The tutorials with examples data sets are in [Tutorial para uso do BioNetStat em linhas de comando](tutorials/tutorial_BNS_linha_de_comando.md) or in [Tutorial BioNetStat para a interface](tutorials/tutorial_BNS_interface.md ).

