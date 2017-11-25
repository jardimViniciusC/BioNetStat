## Instalando o pacote "devtools" para baixar o pacote via GitHub
```
install.packages("devtools")
library(devtools)
```
## Instalando o BioNetStat
```
install_github("jardimViniciusC/BioNetStat")
library(BioNetStat)
```
Para usar a Interface grafica do BioNetStat basta carregar a biblioteca 'shiny' (caso ela não esteja instalada, leia o arquivo README.md) e rodar a função que carrega a interface.
```
library(shiny)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")
```
Um guia de como usar a interface gráfica está nas aba 'Help' da própria interface ou em 'Tutorial BioNetStat para a interface.md'.

## Lendo a tabela de valores das variáveis
A função 'readVarFile', lê apenas os valores numéricos da tabela que serão usados para a construção do grafo. Neste exemplo a tabela 'cancer_data.csv' é formada por valores de expressão de 134 genes.

A função 'doLabels' lê, na mesma tabela, uma coluna que definirá a qual estado cada amostra pertence, definindo quais redes serão comparadas. Neste exemplo as amostras 'cancer_data.csv' são divididas em 4 tecido cancerígenos. Quando não é especificada qual coluna definirá a classificação das amostras a função ira usar a primeira coluna da classe 'factor' (como mostrado no objeto 'labmat'). Quando se deseja especificar qual coluna deverá ser usada e quais estados (tratamentos) serão comparados deve-se usar os argumentos 'factorName' e 'classes', respectivamente (como no objeto labmat2).

Os dados usados neste exemplo 
```
matriz<-readVarFile(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv",sep = ";",dec=".")

labmat<-doLabels(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv")

labmat2<-doLabels(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv",factorName = "histologic_diagnosis",classes = c("Oligodendroglioma","Astrocytoma"))
```
## Lendo o arquivo que indica os grupos de variáveis
A função readGmtFile lê uma tabela que indica os na primeira coluna os nomes dos grupos e nas colunas seguintes as variáveis que pertencem a cada grupo. Neste arquivo, as colunas devem ser delimitadas por tabulação.
```
varSets<-readGmtFile(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/c2.cp.v5.2.symbols.gmt")
```
## Escolhendo os parâmetros para construção da matriz de adjacência

Nessa função o usuário escolhe quais os parâmetros para a construção das redes como o método estátístico usado (method), onde é possível escolher entre as correlação de Pearson, Spearman ou Kendall ou então inserir uma função que retorne uma matriz de adjacência. Se o usuário escolher uma das três correlações ele deve selecionar qual valor será usado como força de associação ("corr", "pvalue", "fdr") em 'association', a força que será usada como threshold na construção das redes ("corr", "pvalue", "fdr") e o valor numérico (entre 0 e 1) usado como threshold, em 'thr.value'. Além disso é possivel escolher se a rede terá peso ou não nas arestas em 'weighted'.
```
funAdjMat<-adjacencyMatrix(method = "pearson",association = "corr",threshold = "corr",thr.value = 0.5,weighted = T)
```

## Comparando as redes

A função 'diffNetAnalysis' realiza o teste de comparação de múltiplas redes. Para comparar as redes é necessário escolher um dos métodos de comparação em 'method'. Nos argumentos 'varFile', 'labels' e 'varSets', o usruário insere os objetos da matriz de valores, a classificação das amostras e os grupos de variáveis, respectivamente. A função que constrói a matriz de correlação é inserida em 'adjacencyMatrix'. O usuário deve definir o numero de permutações para o cálculo do p-valor ('numPermutations'), o numero mínimo de vértices usados para construir as redes ('min.vert'). Se o usuário escolher comparar as redes pelas distribuições do espectro ou de grau ele deve selecionar qual medidad de largura de banda será usada em 'options', podendo ser 'bandwidth'='Silverman' ou 'bandwidth'='Sturges'.
```
# Choose one structural property
metodos<-list(spectralDistributionTest, spectralEntropyTest, degreeDistributionTest,degreeCentralityTest,
           betweennessCentralityTest, closenessCentralityTest, eigenvectorCentralityTest,
           clusteringCoefficientTest)
           
res<-diffNetAnalysis(method = metodos[1],varFile = matriz, labels = labmat, varSets = NULL,adjacencyMatrix = funAdjMat,
                numPermutations = 1000, min.vert = 10,options = list("bandwidth"="Silverman"))
res
```
## Comparando a importância dos vértices

A função 'diffNetAnalysis' também realiza o teste de comparação de vertices em múltiplas redes. Da mesma forma, é necessário escolher um dos métodos de comparação (metodos) em 'method'. Nos argumentos 'varFile', 'labels' e 'varSets', o usruário insere os objetos da matriz de valores, a classificação das amostras e os grupos de variáveis, respectivamente. A função que constrói a matriz de correlação é inserida em 'adjacencyMatrix'. O usuário deve definir o numero de permutações para o cálculo do p-valor ('numPermutations'), o numero mínimo de vértices usados para construir as redes ('min.vert').
```
metodos<-list(degreeCentralityVertexTest,betweennessCentralityVertexTest, closenessCentralityVertexTest, eigenvectorCentralityVertexTest,
              clusteringCoefficientVertexTest)

  res<-diffNetAnalysis(method = metodos[1],options = list("bandwidth"="Silverman"),varFile = matriz,
                       labels = labmat, varSets = NULL,adjacencyMatrix = funAdjMat, numPermutations = 1000, print = T,
                       seed = F,min.vert = 10, resultsFile ="resultados.RData" )
  res$all
```

## Teste das funções que constróem mapas metabólicos

```
nomes.matriz1<-c(do.call(cbind,lapply(strsplit(x=names(matriz),split = "_"), function(x) x[1])))
nomes.matriz2<-c(do.call(cbind,lapply(strsplit(x=names(matriz),split = "_"), function(x) x[2])))
colnames(matriz)<-nomes.matriz1
colnames(matriz)<-nomes.matriz2
res

centralityPathPlot(gene.data = res$all, threshold = NULL, thr.value = 0.5 ,species = "hsa",pathway.id = "05200",file.name = "teste")
```
