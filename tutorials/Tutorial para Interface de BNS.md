# Tutorial para inteface do BioNet Stat

Após a instalação do BioNetStat explicada em [README](/README.md), é necessário apenas abrir o R e rodar os seguintes comandos.
```Rscript
library(devtools)
library(shiny)
library(BioNetStat)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
Ao iniciar o BioNetStat você verá essa figura.
![Image of BioNetStat](https://octodex.github.com/images/yaktocat.png)
