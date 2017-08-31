# BioNetStat

##BioNetStat

Package to perform network analysis. BioNetStat is able to compare two or more correlation networks.

## Installation
### Dependencies
R (>= 3.0.0), shiny (>= 0.8.0), WGCNA, igraph, shinyBS, RColorBrewer, Hmisc, pathway, psych, RJSONIO, whisker, yaml, pheatmap, preprocessCore, GO.db, AnnotationDbi, impute, and ggplot2

### Installation steps

1. If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
2. Open the R command line interface, and install all BNS dependencies (if they have not been installed yet):
````
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
biocLite("GO.db")
biocLite("AnnotationDbi")
biocLite("pathview")
install.packages(c("WGCNA", "igraph", "RColorBrewer", "Hmisc", "psych", "RJSONIO", "whisker", "yaml", "pheatmap", "ggplot2","devtools")) 

````
2a. If pathview don't install try to ilstall this libraries in linux shell terminal.
````
```
3. Install the BioNetStat (BNS) package
````
library(devtools)
install_github("jardimViniciusC/BioNetStat")
````
