# Tutorial para inteface do BioNetStat

Após a instalação do BioNetStat explicada em [README](/README.md), é necessário apenas abrir o R e rodar os seguintes comandos.
```Rscript
library(devtools)
library(shiny)
library(BioNetStat)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
### Tela inicial
Ao iniciar o BioNetStat você verá essa imagem.
![Image of BioNetStat](/inst/shiny/www/images/bionetstat_open_image.png)

Caso você esteja usando o Rstudio uma janela alternativa será aberta, o BNS funcionará de maneira mais eficiente se for usado no browser do seu navegador, por isso clique em "open in browser" na esquerda superior da tela e o BNS será inicia no navegador como na primeira figura.
![Image of BioNetStat_rstudio](/inst/shiny/www/images/bionetstat_open_image_Rstudio.png)

### Carregando os arquivos

Ao clicar no icone 'Browse..' da seção 1(Load data) uma janela se abrirá para que a tabela de valores de variáveis seja selecionada. 
![data_selection](/inst/shiny/www/images/bns_expr_selection.png)
Neste tutorial o arquivo a ser selecionado é o [bnsDataTest_log2.csv](/data/bnsDataTest_log2.csv).
O arquivo será carregado e o programa irá identificar as colunas que são classificadas como 'numeric' pelo R. Uma visão prévia da tabela inserida irá aparecer na sua tela. Além disso, uma seção 'Factors' reconhecerá quais as colunas são classificadas como 'factor' pelo R. Nessa seção é possível selecionar os fatores que serão usados para selecionar os estados (tratamentos, condições) comparados. A seleção dos estados pode ser feita em 'Choose the conditions to be compared:'
![Image of BioNetStatFactors](/inst/shiny/www/images/bns_factor_selection.png)
Após a seleção da tabela de valores de variáveis, é possivel carregar o arquivo 'grupo de variáveis' (opcional) em 'Variable set database'. Em nosso tutorial o arquivo selecionado pode ser [c2.cp.v5.2.symbols.gmt](/data/c2.cp.v5.2.symbols.gmt) que define os grupos de variáveis de acordo com as vias gênicas as quais elas estão associadas. Caso o usuário não carregue nenhum arquivo, o programa irá comparar as redes com todas as variáveis carregadas no arquivo 'Variables values data'
![Image of BioNetStatFiles](/inst/shiny/www/images/bns_set_selection.png)

### Selecionando os parâmetros

Com todos os arquivos de entrada já carregados é necessário escolher os parâmetros utilizados para comparar as redes:
1. Número mínimo e máximo de variáveis (nós) que estarão nas redes comparadas em 'Variable sets size range'. O programa retorna a quantidade de conjuntos de variáveis que estão entre esses valores (figura acima).
2. Medida de dependência usada para inferir a rede de correlação.
3. Força de associação que será usada como limiar para a formação de uma aresta entre dois nós e qual o valor de limiar usado.
4. Tipo de rede que será construída, com ou sem peso nas arestas. No caso de redes com peso, qual a medida de associação será usada como peso das arestas.
![Image of BioNetStatNetParameters](/inst/shiny/www/images/bionetstat_setingParameters_net_image.png)
5. Método de comparação das redes. Se o usuário escolher comparar as redes pelas distribuições do espectro ou de grau ele deve selecionar qual medidad de largura de banda será usada em 'Bandwidth', podendo ser 'Silverman' ou 'Sturges'.
6. Numero de permutações e se será usada uma semente para os testes de permutação aleatória.

![Image of BioNetStatNet2Parameters](/inst/shiny/www/images/bionetstat_setingParameters_net2_image.png)

### Rodando a análise diferencial de múltiplas redes
Após carregar os arquivos e selecionar os parâmetros de análise, clique em 'Start analysis' para realizar a comparação das redes.

![Image of BioNetStatrunning](/inst/shiny/www/images/bionetstat_runningMethod_image.png)

### Resultados
1. Enquanto a análise é realizada, o programa mostrará essa imagem:
![Image of BioNetStatwaiting](/inst/shiny/www/images/bionetstat_execution_image.png)

2. Quando a análise tiver terminado o programa mostrará os parâmetros usados para realizar a análise.
![Image of BioNetStatres1](/inst/shiny/www/images/bionetstat_results1_image.png)

3. A tabela final é apresentada como na figura a seguir:
![Image of BioNetStatres2](/inst/shiny/www/images/bionetstat_results2_image.png)

#### Análise diferencial de vértices
4. É possível também realizar a comparação da importância do vértices. Na aba 'Further analysis' é possível selecionar os grupos de variáveis de acordo com seu critério de corte. O grupo de variáveis é escolhido e a o programa irá gerar a comparação de vértices pela centralidade de grau.

![Image of BioNetStatresVert1](/inst/shiny/www/images/bionetstat_resultsVert1_image.png)

5. Quando a tabela aparece é possível selecionar outros propriedades estruturais para comparar os vértices do grafos. A tabela resultante da análise diferencial de vértices é mostrada como na figura a seguir.

![Image of BioNetStatresVert2](/inst/shiny/www/images/bionetstat_resultsVert2_image.png)
