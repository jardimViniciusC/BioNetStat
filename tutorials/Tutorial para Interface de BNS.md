# Tutorial para inteface do BioNet Stat

Após a instalação do BioNetStat explicada em [README](/README.md), é necessário apenas abrir o R e rodar os seguintes comandos.
```Rscript
library(devtools)
library(shiny)
library(BioNetStat)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
Ao iniciar o BioNetStat você verá essa figura.
![Image of BioNetStat](/inst/shiny/www/images/bionetstat_open_image.png)

Ao clicar no icone 'Browse..' da seção 1(Load data) uma janela se abrirá para que a tabela de valores de variáveis seja selecionada. Neste tutorial a o arquivo a ser selecionado é o 'bnsDataTest.csv'.
O arquivo será carregado e o programa irá identificar as colunas que classificadas como 'numeric' pelo R. Uma visão prévia da tabela inserida irá aparecer na sua tela. Além disso, uma seção 'Factors' reconhecerá quais as colunas são classificadas como 'factor' pelo R. Nessa seção é possível selecionar os fatores que serão usados para selecionar os estados comparados. A seleção dos estados pode ser feita em 'Choose the conditions to be compared:'

Após a seleção da tabela de valores de variáveis, é possivel carregar o arquivo 'grupo de variáveis' (opcional) em 'Variable set database'. Em nosso tutorial o arquivo selecionado pode ser 'c2.cp.v5.2.symbols.gmt' que define os grupos de variáveis de acordo com as vias gênicas as quais elas estão associadas. Caso o usuário não carregue nenhum arquivo, o programa irá comparar as redes com todas as variáveis carregadas no arquivo 'Variables values data'

Após a seleção da tabela de variáveis, a escolha dos estados comparados e do grupo de variáveis o programa ficará como na imagem a seguir:
![Image of BioNetStatFiles](/inst/shiny/www/images/bionetstat_selectedData_image.png)
