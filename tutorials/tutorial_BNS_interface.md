---
output:
  pdf_document: default
  html_document: default
---
# BioNetStat user interface tutorial

Being R and BioNetStat installed, as showed in [README](/README.md), you have just open R and run the following commands.
```Rscript
library(devtools)
library(shiny)
library(BioNetStat)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
### First screen
When you start BioNetStat you will see this image.

<img src="https://github.com/jardimViniciusC/BioNetStat/blob/master/inst/shiny/www/images/bionetstat_open_image.png?raw=true" alt="Image of BioNetStat" height="200"/>

If you are using Rstudio an alternative window will open, therefore BNS will performs well when is used in browser. So, hit "open in browser" in the upper left of the screen and the BNS will open in browser as the first picture.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_open_image_Rstudio.png" alt="Image of BioNetStat_rstudio" height="150"/>

### Loading the files

By clicking on the 'Browse ..' icon in section 1 (Load data) a window will open for the table of variables values to be selected.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_expr_selection.png" alt="data_selection" height="150"/>

In this tutorial the file to be selected is ![bnsDataTest_log2.csv](/data/bnsDataTest_log2.csv).
The file will load and the software will identify the colunms classified as "numeric" by R. A previous view of table is showed on the screen. In addition, a 'Factors' section will recognize which columns are rated as 'factor' by R. In this section you can select the factors that will be used to select the states (treatments, conditions) compared, each state will be a network to be compared. The selection of states can be done in 'Choose the conditions to be compared:'

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_factor_selection.png" alt="Image of BioNetStatFactors" height="150"/>

The tutorial to build your table of variables values is ![here](https://github.com/jardimViniciusC/BioNetStat/blob/master/inst/shiny/help/helpData.html).

After selecting the variable values table, it is possible to load the file 'variable group' (optional) into 'Variable set database'. In our tutorial the selected file is ![c2.cp.v5.2.symbols.gmt](/data/c2.cp.v5.2.symbols.gmt) that sets the groups of variables according to the metabolic pathways which they are associated. If the user does not load any files, the program will compare the networks with all variables loaded in the file 'Variables values data'

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_set_selection.png" alt="Image of BioNetStatFiles" height="150"/><img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_set_selection2.png" alt="Image of BioNetStatFiles" height="150"/>

The tutorial to build your table of variables groups is ![here](https://github.com/jardimViniciusC/BioNetStat/blob/master/inst/shiny/help/helpData.html).

### Parameters setting

After load input files, you migth set the parameters to network construction:
1. Lower and higher numbers of variables (nodes) that will be in the networks compared in 'Variable sets size range'. The program returns the number of sets of variables that are between these values (figure above).
2. Dependency measure used to infer the correlation network.
3. Association strength that will be used as the threshold for the formation of an edge between two nodes and which threshold value is used.
4. Type of network to be built, with or without weight on the edges. In the case of nets with weight, which measure of association will be used as edge weight.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net_image.png" alt="Image of BioNetStatNetParameters" height="150"/><img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net_image2.png" alt="Image of BioNetStatNetParameters" height="150"/>


5. Networks comparison method. If the spectral or degree distribution were used to compare the networks it is necessary to set the method to set the bandwith ("Silveram" or "Sturges" are aviable).
6. Number of permutations and if a seed will be used on the resampling method also could be used.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net2_image.png" alt="Image of BioNetStatNet2Parameters" height="150"/><img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net2_image2.png" alt="Image of BioNetStatNet2Parameters" height="150"/>


### Runnig the differential network analysis
After loading the files and selecting the analysis parameters, click on 'Start analysis' to compare the networks.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_runningMethod_image.png" alt="Image of BioNetStatrunning" height="150"/>

### Results
1. While the analysis is performed, this warning is shown:

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_execution_image.png" alt="Image of BioNetStatwaiting" height="150"/>

2. When the analysis is finished the program will show the parameters used to perform the analysis.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_results1_image.png" alt="Image of BioNetStatres1" height="150"/>

3. The final table is shown as in the following figure:

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_results2_image.png" alt="Image of BioNetStatres2" height="150"/>

#### Node differential analysis
4. It is also possible to compare the importance of nodes. In the 'Further analysis' tab you migth select the groups of variables according to a threshold criteria. The group of variables is chosen and the program will generate the nodes comparison by degree centrality.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert1_image.png" alt="Image of BioNetStatresVert1" height="150"/>

5. When the table appear, it is possible to select other structural properties to compare the nodes importance. The diferential node analysis table is shown as in the following figure:

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert2_image.png" alt="Image of BioNetStatresVert2" height="150"/>

6. Another available functionality is to observe the variables in the metabolic pathways of the database [KEGG] (http://www.kegg.jp/), which allow you to know more in depth the metabolic pathways that the studied variables are realted. The "KEGG pathway visualization" tab, under " ". 

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert3_image.png" alt="Image of BioNetStatresVert3" height="150"/>

you need to insert a table that contains two columns, where the first must be the names of the variables you are studying and the second the respective KEGG code of these variables. These codes can be found on the website itself ![KEGG](http://www.kegg.jp/). In addition to the table with the codes of the variables should be informed the color used in the construction of the graph, if the variables studied are genes / proteins or metabolites. The variables that will appear in the figure are the same ones analyzed in the differential analysis table of the vertex and they can be filtered or according to the value of the test or the pvalor or the value associated with it. It is necessary to choose which ![metabolic pathway](http://www.kegg.jp/kegg/pathway.html) will be used to produce the visualization and in which ![Species](http://www.kegg.jp/kegg/catalog/org_list.html).

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert4_image.png" alt="Image of BioNetStatresVert4" height="150"/>

No mapa KEGG, você observará um que os retangulos ou os círculos estão divididos em colunas, representando os tratamentos que você escolheu comparar. A intesidade das cores é relativa aos valores de centralidades escolhidos conforme a legenda. a seguir um exemplo de um mapa construido com genes da via "pathway in cancer" (05200) do kegg, para humanos ("hsa"). 

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert5_image.png" alt="Image of BioNetStatresVert5" height="150"/>

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert6_image.png" alt="Image of BioNetStatresVert6" height="150"/>

Ao clicar no botão para salver é necessário esperar para que o programa faça o donwload do mapa, salve-o no diretório "Downloads" em um arquivo compactado. A visualização se encontra dentro do arquivo compactado.

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert7_image.png" alt="Image of BioNetStatresVert7" height="150"/>

#### Visualização das redes

A visualização das redes nos ajuda a compreender as mudanças dos padrões de maneira mais facil. Essa visualização pode ser feita de diversas maneiras e algumas delas estão disponíveis no BioNetStat na aba "Network visualization plots". Como diversos estados podem ser comparados pelo BNS nesta seção é necessário escolher apenas dois deles para serem visualizados e comparados por vez.

![Image of BioNetStatclassSelection](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_class_selection.png)

A primeira forma de visualizar as redes é em forma de heatmaps das matrizes de adjacência (matriz que representa a rede). inicialmente o usuário pode escolher os parametros (cor,formato e tamanho) para a construção da figura.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_plot_settings.png)

7. E então pode visualizar as redes dos estados selecionados.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_network_vis.png)

8. É possível o usuário comparar a associação entre pares de variáveis escolhidas nos tratamentos analisados:
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_assoc.png)

9. As redes produzidas podem ser comparadas através de propriedades estruturais globais como entropia espectral da rede ou centralidade de grau média:
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_top_properties.png)

10. A matriz de diferenças permite analizar para todos os pares de variáveis como sua relação muda entre os estados.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_comparison_net.png)

11. A lista de associação entre as variáveis no permite ter uma ideia mais precisa das diferenças entre as forças de associação nos dois estados.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_association_list.png)

12. Por fim, a rede pode ser exportada para o programa de integração de níveis biológicos S.I.T. (System Integration Tool) onde a rede pode ser visualizadas e manipulada.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/blob/master/inst/shiny/www/images/bionetstat_sit_list.png)

#### Estudo o comportamento das variáveis

O BioNetStat disponibiliza também um estudo do comportamento das variáveis para que o usuário possa ter uma visão completa do seu objeto de estudo.
1. Inicialmente as variáveis podem ser exploradas por meio da construção de um heatmap onde os parametros da figura são selecionados.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_parameters_heatmap.png)

e o heatmap é produzido. Nas colunas estão as amostras com as cores indicando à qual estado cada amostra pertence. É possível agrupar tanto as colunas (amostras) quanto as linhas (variáveis) do heatmap na seção "heatmap cluster options" na figura acima.
![Image of BioNetStatplotSettings](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_heatmap.png)

2. Nos mapas do banco de dados KEGG, podemos observar também as expressões do genes, concentrações da proteínas e dos metabólitos. Em "KEGG pathway visualization" o usuário deve inserir o código das variáveis, uma tabela que contem duas colunas, onde na primeira deve estar os nomes das variáveis que você está estudando e na segunda o respectivo código KEGG dessas variáveis. Esses códigos podem ser encontrados no próprio site do [KEGG](http://www.kegg.jp/). As variáveis que irão aparecer na figura são as mesmas analizadas na tabela de análise diferencial do vértice e elas podem ser filtradas ou de acordo com o valor do teste ou o pvalor ou o qvalor associado a ele. Na figura será mostrado um valor que representará todas as amostras de cada estado, que pode ser a média, a mediana, o mínimo ou o máximo. É necessário que se escolha qual [via metabólica](http://www.kegg.jp/kegg/pathway.html) será usada para produzir a visualização e em qual [especie](http://www.kegg.jp/kegg/catalog/org_list.html).

![Image of BioNetStatresVert5](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_keeg_parameters_exp.png).
![Image of BioNetStatresVert5](https://github.com/jardimViniciusC/BioNetStat/blob/master/inst/shiny/www/images/bionetstat_keeg_parameters_exp_codes.png).

No mapa KEGG, você observará um que os retangulos ou os círculos estão divididos em colunas, representando os tratamentos que você escolheu comparar. A intesidade das cores é relativa aos valores de centralidades escolhidos conforme a legenda. a seguir um exemplo de um mapa construido com genes da via "pathway in cancer" (05200) do kegg, para humanos ("hsa"). 
![Image of BioNetStatresVert5](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_keeg_parameters_exp2.png)
![Image of BioNetStatresVert6](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert6_image.png)
Ao clicar no botão para salver é necessário esperar para que o programa faça o donwload do mapa, salve-o no diretório "Downloads" em um arquivo compactado. A visualização se encontra dentro do arquivo compactado.
![Image of BioNetStatresVert7](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert7_image.png)

3. O BioNetStat disponibiliza a realização de um teste estatístico para verificar a mudança no comportamento das variáveis entre os estados (teste t, para dois estados, e ANOVA, para mais de dois, com os respectivos testes não paramétricos). Os resultados so apresentados na tabela como a seguir:

![Image of BioNetStatresVert7](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_anova.png)

4. Por fim, a última análise disponível no BioNetStat é a visualização dos valores de uma determinada variável em um grafico boxplot. Onde inicialmente o usuário deve escolher os parametros da figura formada

![Image of BioNetStatresVert7](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_boxplot_parameters.png)

e posteriormente a variável que será visualizada.

![Image of BioNetStatresVert7](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_boxplot.png)
