################################################################################################################
#' Read variable values matrix
#' @param fileName the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param path the path to the directory that contains the file
#' @param dec the character used in the file for decimal points.
#' @param sep the field separator character. Values on each line of the file are separated by this character. If sep = "" the separator is ‘white space’, that is one or more spaces, tabs, newlines or carriage returns, if sep=NULL (default), the function uses tabulation for .txt files or ";" for .csv files.
#' @export
#'
readVarFile <- function(fileName,path=NULL,dec=".",sep=NULL){#readSampleTable
  if(is.null(path)) path<-fileName
  if(lapply(strsplit(as.character(path), '[.]'),rev)[[1]][1]=="txt"){ # Se o arquivo for TXT
    if(is.null(sep)) sep="\t"
    expr <- read.table(fileName, header=TRUE,dec=dec,sep=sep)
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
    if(is.null(sep)) sep=";"
    table <- read.table(fileName,header=T,dec=dec,sep=sep) # atentar para o decimal como virgula
    expr <- table[,sapply(table,is.numeric)]
    # names <- expr[1,]
    n <- nrow(expr)
    # expr <- as.matrix(expr)
    # colnames(expr) <- names
    return(expr)
    }
}
################################################################################################################
# chooseClass <- function(fileName) {
#   table <- read.csv(fileName,header=T,dec=",",sep=";") # atentar para o decimal como virgula e separador ponto e virgula
#   class <- names(which(!sapply(table,is.numeric)))
#   return(class)
# }
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
################################################################################################################

#' Class vector of data table
#'
#' @param fileName the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param factorName string indicating the column name used to determine the labels of each row of matrix data. The NULL (default) indicates that the first column will be used.
#' @param classes a vector of strings indicating which labels of choosed column will be compared, the minimum are two labels. The NULL (default) indicates that all classes will be compared.
#' @param dec the character used in the file for decimal points.
#' @param sep the field separator character. Values on each line of the file are separated by this character. If sep = "" the separator is ‘white space’, that is one or more spaces, tabs, newlines or carriage returns, if sep=NULL (default), the function uses tabulation for .txt files or ";" for .csv files.
#' @export
doLabels <- function(fileName, factorName=NULL, classes=NULL,dec=".",sep=";") {
  options(stringsAsFactors = T)
  table <- read.csv(fileName,header=T,dec=dec,sep=sep)
  if(is.null(factorName)) factor <- names(which(!sapply(table,is.numeric)))[1]
  else if(!any(factorName==names(which(!sapply(table,is.numeric))))) stop(paste("The factorName",factorName," doesn't exists in Data frame"))
  else factor <- factorName
  labels<-table[,factor]
  if(is.null(classes)) classes<-levels(labels)
  else if(length(classes)<2) stop(paste("Just one class was indicated. More than one class of",factorName,"has to be indicated"))
  if(any(!c(classes %in% levels(labels)))) if(sum(!c(classes %in% levels(labels)))==1) stop(paste("Class",classes[!c(classes %in% levels(labels))],"isn't contained in",factorName))
  else stop(paste("Classes",paste(classes[!c(classes %in% levels(labels))],collapse = ", "),"aren't contained in",factorName))
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

################################################################################################################
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
      stop(paste("Wrong variable sets file format. All sets should",
                 "contain at least one variable"))
  }
  return(geneSets)
}

################################################################################################################
#' Class vector of data table
#'
#' @param method a function that receives two adjacency matrices and returns a list containing a statistic theta that measures the difference between them, and a p-value for the test H0: theta = 0 against H1: theta > 0.
#' @param options a list contaning paremeters used by 'method'.
#' @param varFile a numeric matrix contaning variables values data.
#' @param labels a vector of -1s, 0s, and 1s associating each sample with a phenotype. The value 0 corresponds to the first phenotype class of interest, 1 to the second phenotype class of interest, and -1 to the other classes, if there are more than two classes in the gene expression data.
#' @param varSets a list of gene sets. Each element of the list is a character vector v, where v[1] contains the gene set name, v[2] descriptions about the set, v[3..length(v)] the genes that belong to the set.
#' @param adjacencyMatrix a function that receives a numeric matrix containing gene expression data and returns the adjacency matrix of the inferred co-expression graph.
#' @param numPermutations the number of permutations for the permutation test.
#' @param print a logical. If true, it prints execution messages on the screen. resultsFile: path to a file where the partial results of the analysis will be saved. If NULL, then no partial results are saved.
#' @param resultsFile a ".RData" file name to be saved in tha work directory.
#' @param seed the seed for the random number generators. If it is not null then the sample permutations are the same for all the gene sets.
#' @param min.vert
#' @return a data frame containing the name, size, test statistic, nominal p-value and adjusted p-value (q-value) associated with each gene set.
#' @export
#'
 diffNetAnalysis <- function(method, options, varFile, labels, varSets=NULL,
                            adjacencyMatrix, numPermutations=1000, print=TRUE,
                            resultsFile=NULL, seed=NULL, min.vert=5) {
   if(is.null(varSets)) varSets <- list(c("all",colnames(varFile)))
   # else varSets[[length(varSets)+1]] <- c("all",colnames(varFile))
   varSets<-lapply(varSets, function(x){
     if(sum(x %in% colnames(varFile))>=min.vert) x
     else NA})#
   varSets<-varSets[!is.na(varSets)]
  results <- data.frame(matrix(NA, nrow=length(varSets), ncol=5+length(unique(labels[labels!=-1]))))
  output<-list()
  names <- array(NA, length(varSets))
  for (i in 1:length(varSets)) {
    names[i] <- varSets[[i]][1]
  }
  rownames(results) <- names
  if(any(labels=="-1"))colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste("Factor",unique(labels[-which(labels=="-1")])))
  else colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste("Factor",unique(labels)))
  #temp <- tempfile("CoGA_results", fileext=".txt")
  for (i in 1:length(varSets)) {
    setName <- varSets[[i]][1]
    if (print)
      cat(paste("Testing ", setName, " (", i, " of ", length(varSets),
                ")", "\n", sep=""), append=T)
    genes <- varSets[[i]][varSets[[i]] %in% colnames(varFile)]
    if (!is.null(seed))
      set.seed(seed)
    result <- method(varFile[,genes], labels, adjacencyMatrix=
                       adjacencyMatrix,  numPermutations=
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
    }
    # if (!is.null(resultsFile)) save(results, file=resultsFile)
  }
  if(is.list(result)) results[, "Q-value"] <- p.adjust(results[, "Nominal p-value"], method="fdr")
  return(results)
}

 #' Run BNS
 #'
 #' Run BNS on the browser user interface.
 #' @export
 runBNS <- function() {
   shiny::runApp(system.file('shiny', package='BioNetStat'))
 }
 