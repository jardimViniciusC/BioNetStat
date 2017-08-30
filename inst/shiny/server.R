library(shiny);library(shinyBS);library(BioNetStat)
source("~/Dropbox/mestrado/pacote_r/codigo_coga_vinicius/codigo_para_implmentar_anova/creating_BioNetStat_com_git/shiny_BioNetStat/shiny/global.R")
# ------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------

# Vector that maps the names of methods for collapsing rows to functions
collapsingMethodsMatrix <- as.matrix(read.table("collapsingMethods.txt",  header=T))
collapsingMethods <- c()
for (i in 1:nrow(collapsingMethodsMatrix)) {
  if (!is.na(collapsingMethodsMatrix[i,"Source"]))
    source(collapsingMethodsMatrix[i,"Source"])
  collapsingMethods[i] <- collapsingMethodsMatrix[i,"Name"]
}
rownames(collapsingMethodsMatrix) <- collapsingMethods

# Vector that maps the names of methods for network comparison to functions
networkTestsMatrix <- as.matrix(read.table("networkTests.txt",  header=T))
networkTests <- c()
for (i in 1:nrow(networkTestsMatrix)) {
  if (!is.na(networkTestsMatrix[i,"Source"]))
    source(networkTestsMatrix[i,"Source"])
  networkTests[i] <- networkTestsMatrix[i,"Name"]
}
rownames(networkTestsMatrix) <- networkTests

# Matrix that maps dependence measure names to functions 
correlationMeasuresMatrix <- as.matrix(read.table("dependenceMeasures.txt",  header=T))
correlationMeasures <- matrix(NA, nrow(correlationMeasuresMatrix), 2)
rownames(correlationMeasures) <- correlationMeasuresMatrix[,1]
colnames(correlationMeasures) <- colnames(correlationMeasuresMatrix)[2:3]
for (i in 1:nrow(correlationMeasuresMatrix)) {
  if (!is.na(correlationMeasuresMatrix[i,4]))
    source(correlationMeasuresMatrix[i,4])
  correlationMeasures[i, ] <- correlationMeasuresMatrix[i,2:3]
}

# Vector that maps gene score methods to functions
geneScoresMatrix <- as.matrix(read.table("geneScores.txt",  header=T))
geneScores <- c()
for (i in 1:nrow(geneScoresMatrix)) {
  if (!is.na(geneScoresMatrix[i,3]))
    source(geneScoresMatrix[i,3])
  geneScores[i] <- geneScoresMatrix[i,1]
}
rownames(geneScoresMatrix) <- geneScores

# Vector that maps gene set score methods network comparison to functions
networkScoresMatrix <- as.matrix(read.table("networkScores.txt",  header=T))
networkScores <- c()
for (i in 1:nrow(networkScoresMatrix)) {
  if (!is.na(networkScoresMatrix[i,"Source"]))
    source(networkScoresMatrix[i,"Source"])
  networkScores[i] <- networkScoresMatrix[i,"Name"]
}
rownames(networkScoresMatrix) <- networkScores


# ------------------------------------------------------------------------------
# Defining the server logic
# ------------------------------------------------------------------------------

# Define server logic 
shinyServer(function(input, output, session) {
  
  values <- reactiveValues(canExecute=F, completed=F, expr=NULL, labels=NULL, 
                           filteredGeneSets=NULL, classes=NULL, 
                           associationMeasure=NULL, correlationMeasure=NULL, 
                           networkType=NULL, threshold=NULL, networkTest=NULL,
                           seed=NULL, options=NULL, exprInputFileName=NULL,
                           labelsInputFileName=NULL, geneSetsFileName=NULL)
  
  # Loading data -------------------------------------------------------------
  
  source("dataLoading.R", local=T)
  
  # Setting execution parameters ---------------------------------------------
  
  source("executionParameters.R", local=T)
  
  # Showing analysis results -------------------------------------------------
  
  source("results.R", local=T)
  
  # Further analysis ---------------------------------------------------------
  
  source("furtherAnalysis.R", local=T)
  
})