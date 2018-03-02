---
output:
  html_document: default
  pdf_document: default
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

In this tutorial the file to be selected is [bnsDataTest_log2.csv](/data/bnsDataTest_log2.csv).
The file will load and the software will identify the colunms classified as "numeric" by R. A previous view of table is showed on the screen. In addition, a 'Factors' section will recognize which columns are rated as 'factor' by R. In this section you can select the factors that will be used to select the states (treatments, conditions) compared, each state will be a network to be compared. The selection of states can be done in 'Choose the conditions to be compared:'

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_factor_selection.png" alt="Image of BioNetStatFactors" height="150"/>

The tutorial to build your table of variables values is here(INSERT A LINK).

After selecting the variable values table, it is possible to load the file 'variable group' (optional) into 'Variable set database'. In our tutorial the selected file is [c2.cp.v5.2.symbols.gmt](/data/c2.cp.v5.2.symbols.gmt) that sets the groups of variables according to the metabolic pathways which they are associated. If the user does not load any files, the program will compare the networks with all variables loaded in the file 'Variables values data'

<img src="https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bns_set_selection.png" alt="Image of BioNetStatFiles" height="150"/>

The tutorial to build your table of variables groups is here(INSERT A LINK).

### Selecionando os parâmetros

Com todos os arquivos de entrada já carregados é necessário escolher os parâmetros utilizados para comparar as redes:
1. Número mínimo e máximo de variáveis (nós) que estarão nas redes comparadas em 'Variable sets size range'. O programa retorna a quantidade de conjuntos de variáveis que estão entre esses valores (figura acima).
2. Medida de dependência usada para inferir a rede de correlação.
3. Força de associação que será usada como limiar para a formação de uma aresta entre dois nós e qual o valor de limiar usado.
4. Tipo de rede que será construída, com ou sem peso nas arestas. No caso de redes com peso, qual a medida de associação será usada como peso das arestas.

![Image of BioNetStatNetParameters](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net_image.png)

5. Método de comparação das redes. Se o usuário escolher comparar as redes pelas distribuições do espectro ou de grau ele deve selecionar qual medidad de largura de banda será usada em 'Bandwidth', podendo ser 'Silverman' ou 'Sturges'.
6. Numero de permutações e se será usada uma semente para os testes de permutação aleatória.

![Image of BioNetStatNet2Parameters](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_setingParameters_net2_image.png)

### Rodando a análise diferencial de múltiplas redes
Após carregar os arquivos e selecionar os parâmetros de análise, clique em 'Start analysis' para realizar a comparação das redes.

![Image of BioNetStatrunning](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_runningMethod_image.png)

### Resultados
1. Enquanto a análise é realizada, o programa mostrará essa imagem:
![Image of BioNetStatwaiting](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_execution_image.png)

2. Quando a análise tiver terminado o programa mostrará os parâmetros usados para realizar a análise.
![Image of BioNetStatres1](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_results1_image.png)

3. A tabela final é apresentada como na figura a seguir:
![Image of BioNetStatres2](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_results2_image.png)

#### Análise diferencial de vértices
4. É possível também realizar a comparação da importância do vértices. Na aba 'Further analysis' é possível selecionar os grupos de variáveis de acordo com seu critério de corte. O grupo de variáveis é escolhido e o programa irá gerar a comparação de vértices pela centralidade de grau.

![Image of BioNetStatresVert1](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert1_image.png)

5. Quando a tabela aparece, é possível selecionar outras propriedades estruturais para comparar os vértices do grafos. A tabela resultante da análise diferencial de vértices é mostrada como na figura a seguir.

![Image of BioNetStatresVert2](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert2_image.png)

6. Outra funcionalidade disponível é observar as variáveis em vias metabólicas do banco de dados [KEGG](http://www.kegg.jp/), o que nos permitire conhecer mais a fundo as vias metabólicas que as variáveis estudadas estão inseridas. A aba "KEGG pathway visualization", em "" 
![Image of BioNetStatresVert2](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert3_image.png)

é preciso inserir uma tabela que contem duas colunas, onde na primeira deve estar os nomes das variáveis que você está estudando e na segunda o respectivo código KEGG dessas variáveis. Esses códigos podem ser encontrados no próprio site do [KEGG](http://www.kegg.jp/). Além da tabela com os códigos das variáveis deverá ser informado a cor usada na construção do gráfico, se as variáveis estudadas são genes/proteínas ou metabólitos. As variáveis que irão aparecer na figura são as mesmas analizadas na tabela de análise diferencial do vértice e elas podem ser filtradas ou de acordo com o valor do teste ou o pvalor ou o qvalor associado a ele. É necessário que se escolha qual [via metabólica](http://www.kegg.jp/kegg/pathway.html) será usada para produzir a visualização e em qual [especie](http://www.kegg.jp/kegg/catalog/org_list.html). 

![Image of BioNetStatresVert5](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert5_image.png)

No mapa KEGG, você observará um que os retangulos ou os círculos estão divididos em colunas, representando os tratamentos que você escolheu comparar. A intesidade das cores é relativa aos valores de centralidades escolhidos conforme a legenda. a seguir um exemplo de um mapa construido com genes da via "pathway in cancer" (05200) do kegg, para humanos ("hsa"). 

![Image of BioNetStatresVert5](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert5_image.png).

![Image of BioNetStatresVert6](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert6_image.png)

Ao clicar no botão para salver é necessário esperar para que o programa faça o donwload do mapa, salve-o no diretório "Downloads" em um arquivo compactado. A visualização se encontra dentro do arquivo compactado.
![Image of BioNetStatresVert7](https://github.com/jardimViniciusC/BioNetStat/raw/master/inst/shiny/www/images/bionetstat_resultsVert7_image.png)

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
