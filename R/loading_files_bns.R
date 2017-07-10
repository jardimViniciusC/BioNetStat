
readExprTxtFile <- function(fileName,path=NULL,dec=".",sep.csv=";",sep.txt="\t"){#readSampleTable
  if(is.null(path)) path<-fileName
  if(lapply(strsplit(as.character(path), '[.]'),rev)[[1]][1]=="txt"){ # Se o arquivo for TXT
    expr <- read.table(fileName, header=TRUE,dec=dec,sep=sep.txt)
    colnames(expr) <- toupper(colnames(expr))
    cols <- colnames(expr)
    if (cols[1] != "NAME")
      stop(paste("Wrong gene expression data format. The fist column",
                 "of the file should be \"Name\"."))
    
    names <- expr[, "NAME"]
    n <- ncol(expr)
    
    if (cols[2] == "DESCRIPTION")
      expr <- expr[, 3:n]
    else
      expr <- expr[, 2:n]
    
    expr <- as.matrix(expr)
    rownames(expr) <- names
    return(expr) ## Lembrar que a matrix agora está com os genes na primeira coluna
  }
  if(lapply(strsplit(as.character(path), '[.]'),rev)[[1]][1]=="csv"){ # Se o arquivo for CSV
    table <- read.table(fileName,header=T,dec=dec,sep=sep.csv) # atentar para o decimal como virgula
    expr <- table[,sapply(table,is.numeric)]
    # names <- expr[1,]
    n <- nrow(expr)
    # expr <- as.matrix(expr)
    # colnames(expr) <- names
    return(expr)
    }
}
# # chooseClass(fileName = "~/Dropbox/mestrado/Link para teste de rede 03_02/tabela_rede experimento 1.csv")
# # fileName = "~/Dropbox/mestrado/Link para Carmem/carmem.csv"
# chooseClass <- function(fileName) {
#   table <- read.csv(fileName,header=T,dec=",",sep=";") # atentar para o decimal como virgula e separador ponto e virgula
#   class <- names(which(!sapply(table,is.numeric)))
#   return(class)
# }
# # factorName="hor"
# readClass <- function(fileName, factorName) {
#   table <- read.csv(fileName,header=T,dec=",",sep=";") # atentar para o decimal como virgula
#   class <- as.factor(table[,which(names(table)==factorName)])
#   # labels <- combn(levels(class),2)
#   labels <- levels(class)
#   return(labels)
# }
# 
# labelsAux <- function(labels, class1, class2) {
#   l <- array(NA, length(labels))
#   names <- levels(labels)
#   i <- which(names == class1)
#   j <- which(names == class2)
#   symbols <- unique(labels)
#   i1 <- which(labels == symbols[i])
#   i2 <- which(labels == symbols[j])
#   l[i1] <- 0
#   l[i2] <- 1 
#   l[-(union(i1, i2))] <- -1
#   return(l)
# }
# 
# labelsAux <- function(labels, classes) {
#   l <- array(NA, length(labels))
#   names <- levels(labels)
#   i<-vector(len=length(classes))
#   for(p in 1:length(i)) i[p] <- which(names == classes[p])
#   
#   symbols <- unique(labels)
#   j<-list()
#   v<-0
#   for(p in 1:length(i)){
#     j[[p]] <- which(labels == symbols[i[p]])
#     l[j[[p]]] <- v
#     v <- v+1
#   }
#   l[is.na(l)] <- -1
#   return(l)
# }

doLabels <- function(fileName, factorName=NULL, classes=NULL,dec=",",sep=";") {
  options(stringsAsFactors = T)
  table <- read.csv(fileName,header=T,dec=dec,sep=sep) # atentar para o decimal como virgula e separador ponto e virgula
  if(is.null(factorName)) factor <- names(which(!sapply(table,is.numeric)))[1]
  else if(!any(factorName==names(which(!sapply(table,is.numeric))))) stop(paste("The factorName",factorName," doesn't exists in Data frame"))
  else factor <- factorName
  labels<-table[,factor]
  if(is.null(classes)) classes<-levels(labels)
  
  l <- array(NA, length(labels))
  names <- unique(labels)
  i<-vector(len=length(classes))
  for(p in 1:length(i)) i[p] <- which(names == classes[p])
  
  symbols <- unique(labels)
  j<-list()
  v<-0
  for(p in 1:length(i)){
    j[[p]] <- which(labels == symbols[i[p]])
    l[j[[p]]] <- v
    v <- v+1
  }
  l[is.na(l)] <- -1
  return(l)
}

#' Read a collection of gene sets (*.gmt)
#'
#' 'readGmtFile' reads a tab-delimited text file containing a collection of 
#' gene sets.
#' @param fileName a string containing the file name
#' @return a list of gene sets. Each element of the list is a character vector
#' v, where v[1] contains the gene set name, v[2] descriptions about the set,
#' v[3..length(v)] the genes that belong to the set.
#' @export
readGmtFile <- function(fileName) {
  lines <- readLines(fileName)
  geneSets <- strsplit(lines, "\t")
  n <- length(geneSets)
  for (i in 1:n) {
    if (length(geneSets[[i]]) < 3)
      stop(paste("Wrong gene sets file format. All gene sets should", 
                 "contain at least one gene."))
  }
  return(geneSets)
}
# readGmtFile("~/Dropbox/mestrado/Link para Carmem/cargmt.gmt")

 diffNetAnalysis <- function(method, options, expr, labels, geneSets=NULL, 
                            adjacencyMatrix, numPermutations=1000, print=TRUE, 
                            resultsFile=NULL, seed=NULL, min.vert=5) {
   if(is.null(geneSets)) geneSets <- list(c("all",colnames(expr)))
   # else geneSets[[length(geneSets)+1]] <- c("all",colnames(expr))
   geneSets<-lapply(geneSets, function(x){ 
     if(sum(x %in% colnames(expr))>=min.vert) x
     else NA})# 
   geneSets<-geneSets[!is.na(geneSets)]
  results <- data.frame(matrix(NA, nrow=length(geneSets), ncol=5+length(unique(labels[labels!=-1]))))
  output<-list()
  names <- array(NA, length(geneSets))
  for (i in 1:length(geneSets)) {
    names[i] <- geneSets[[i]][1]
  }
  rownames(results) <- names
  if(any(labels=="-1"))colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste("Factor",unique(labels[-which(labels=="-1")])))
  else colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste("Factor",unique(labels)))
  #temp <- tempfile("CoGA_results", fileext=".txt")
  for (i in 1:length(geneSets)) {
    setName <- geneSets[[i]][1]
    if (print)
      cat(paste("Testing ", setName, " (", i, " of ", length(geneSets), 
                ")", "\n", sep=""), append=T)
    genes <- geneSets[[i]][geneSets[[i]] %in% colnames(expr)]
    if (!is.null(seed))
      set.seed(seed)
    result <- method(expr[,genes], labels, adjacencyMatrix=
                       adjacencyMatrix1,  numPermutations=
                       numPermutations, options=options)
    if(!is.list(result)){
      output[[setName]]<-result
      results<-output
    }
    if(is.list(result)){
    results[setName, "N of Networks"] <- sum(unique(labels)!=-1)
    results[setName, "Test statistic"] <- result[[1]]
    results[setName, "Nominal p-value"] <- result$p.value
    results[setName, "Set size"] <- length(genes)
    results[setName, 6:ncol(results)] <- result$Partial*100/sum(result$Partial)
    if (!is.null(resultsFile))
      save(results, file=resultsFile)
    }
  }
  if(is.list(result)) results[, "Q-value"] <- p.adjust(results[, "Nominal p-value"], method="fdr")
  return(results)
}

