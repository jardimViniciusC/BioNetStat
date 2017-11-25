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
library(shiny)
runGitHub("jardimViniciusC/BioNetStat",subdir = "inst/shiny")

## Lendo a tabela de valores das variáveis
A função 'readVarFile', lê apenas os valores numéricos da tabela que serão usados para a construção do grafo. Neste exemplo a tabela 'cancer_data.csv' é formada por valores de expressão de 134 genes.

A função 'doLabels' lê, na mesma tabela, uma coluna que definirá a qual estado cada amostra pertence, definindo quais redes serão comparadas. Neste exemplo as amostras 'cancer_data.csv' são divididas em 4 tecido cancerígenos. Quando não é especificada qual coluna definirá a classificação das amostras a função ira usar a primeira coluna da classe 'factor' (como mostrado no objeto 'labmat'). Quando se deseja especificar qual coluna deverá ser usada e quais estados (tratamentos) serão comparados deve-se usar os argumentos 'factorName' e 'classes', respectivamente (como no objeto labmat2).

Os dados usados neste exemplo 
```
matriz<-readVarFile(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv",sep = ";",dec=".")

nomes.matriz1<-c(do.call(cbind,lapply(strsplit(x=names(matriz),split = "_"), function(x) x[1])))
nomes.matriz2<-c(do.call(cbind,lapply(strsplit(x=names(matriz),split = "_"), function(x) x[2])))
colnames(matriz)<-nomes.matriz1

labmat<-doLabels(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv")

labmat2<-doLabels(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/cancer_data.csv",factorName = "histologic_diagnosis",classes = c("Oligodendroglioma","Astrocytoma"))
```
## Lendo o arquivo que indica os grupos de variáveis
A função readGmtFile lê uma tabela que indica os na primeira coluna os nomes dos grupos e nas colunas seguintes as variáveis que pertencem a cada grupo. Neste arquivo, as colunas devem ser delimitadas por tabulação.
```
varSets<-readGmtFile(fileName = "~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/dados_de_teste/c2.cp.v5.2.symbols.gmt")
```
## Escolhendo os parâmetros para construção da matriz de adjacência
```
funAdjMat<-adjacencyMatrix(method = "pearson",association = "corr",threshold = "corr",thr.value = 0.5,weighted = T)
```
# Teste das funções que comparam redes
```
metodos<-list(spectralDistributionTest, spectralEntropyTest, degreeDistributionTest,degreeCentralityTest,
           betweennessCentralityTest, closenessCentralityTest, eigenvectorCentralityTest,
           clusteringCoefficientTest)
res<-list()
for(i in 1:length(metodos)){
  res<-diffNetAnalysis(method = metodos[i],options = list("bandwidth"="Silverman"),varFile = matriz,
                labels = labmat, varSets = NULL,adjacencyMatrix = funAdjMat, numPermutations = 10, print = T,
                seed = F,min.vert = 10, resultsFile ="resultados.RData" )
  print(res)

}
```
############## Teste das funções que comparam vertices #######################

metodos<-list(degreeCentralityVertexTest,betweennessCentralityVertexTest, closenessCentralityVertexTest, eigenvectorCentralityVertexTest,
              clusteringCoefficientVertexTest)

for(i in metodos){
  res<-diffNetAnalysis(method = i,options = list("bandwidth"="Silverman"),varFile = matriz,
                       labels = labmat, varSets = NULL,adjacencyMatrix = funAdjMat, numPermutations = 10, print = T,
                       seed = F,min.vert = 10, resultsFile ="resultados.RData" )
  print(head(res$all))

}

############## Teste das funções que constróem mapas metabólicos #######################
colnames(matriz)<-nomes.matriz2
res

centralityPathPlot(gene.data = res$all, threshold = NULL, thr.value = 0.5 ,species = "hsa",pathway.id = "05200",file.name = "teste")
