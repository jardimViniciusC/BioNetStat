################################################################################################################
#' Run BNS
#'
#' Run BNS on the browser user interface.
#' @return open BioNetStat user interface
#' @examples  # run  runBioNetStat() # to open user interface of BioNetStat
#' @export
runBioNetStat <- function(){
  shiny::runApp(system.file('shiny', package='BioNetStat'))
}

################################################################################################################
#' Read variable values matrix
#' @param fileName the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param path the path to the directory that contains the file. Used only by Graphical Interface.
#' @param dec the character used in the file for decimal points.
#' @param sep the field separator character. Values on each line of the file are separated by this character. If sep = "" the separator is white space, that is one or more spaces, tabs, newlines or carriage returns, if sep=NULL (default), the function uses tabulation for .txt files or ";" for .csv files.
#' @param check.names a logical value. If TRUE, the names of the data table kept as they are. Otherwise, the blank space, "-","/" and ",", are replaced by dots.
#' @return a dataframe containing only the numeric columns of selected file. Each column is considered as a variable and each row as a sample.
#' @examples 
#' # Glioma file
#' gliomaData <- system.file("extdata", "bnsData.csv", package = "BioNetStat")
#' varFile<-readVarFile(gliomaData)
#' 
#' # Random file
#' test1 <- as.data.frame(cbind(rep(LETTERS[1:4],each=10),matrix(rnorm(120),40,30)))
#' tf<-tempfile(fileext = ".csv")
#' write.table(test1, tf,sep=";",row.names=FALSE)
#' a<-readVarFile(fileName=tf)
#' @export
#'
readVarFile <- function(fileName,path=NULL,dec=".",sep=NULL,check.names=TRUE){#readSampleTable
  if(is.null(path)) path<-fileName
  if(lapply(strsplit(as.character(path), '[.]'),rev)[[1]][1]=="txt"){ # Se o arquivo for TXT
    if(is.null(sep)) sep="\t"
    expr <- read.table(fileName, header=TRUE,dec=dec,sep=sep,check.names = check.names)
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
    table <- read.table(fileName,header=TRUE,dec=dec,sep=sep,check.names = check.names) # atentar para o decimal como virgula
    expr <- table[,vapply(table,is.numeric,FUN.VALUE = vector(length = 1))]
    n <- nrow(expr)
    return(expr)
  }
}
################################################################################################################

#' Class vector of data table
#'
#' @param fileName the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param factorName string indicating the column name used to determine the labels of each row of matrix data. The NULL (default) indicates that the first column will be used.
#' @param classes a vector of strings indicating which labels of choosed column will be compared, the minimum are two labels. The NULL (default) indicates that all classes will be compared.
#' @param dec the character used in the file for decimal points.
#' @param sep the field separator character. Values on each line of the file are separated by this character. If sep = "" the separator is white space, that is one or more spaces, tabs, newlines or carriage returns, if sep=NULL (default), the function uses tabulation for .txt files or ";" for .csv files.
#' @return a vector that identify each row of the readVarFile object as a sample belonging to a state (network). 
#' @examples 
#' # Glioma file
#' gliomaData <- system.file("extdata", "bnsData.csv", package = "BioNetStat")
#' labels<-doLabels(gliomaData)
#' 
#' # Random file
#' test1 <- as.data.frame(cbind(rep(LETTERS[1:4],each=10),matrix(rnorm(120),40,30)))
#' tfl<-tempfile(fileext = ".csv")
#' write.table(test1, tfl,sep=";",row.names=FALSE)
#' labels<-doLabels(tfl)
#' @export
doLabels <- function(fileName, factorName=NULL, classes=NULL,dec=".",sep=";") {
  options(stringsAsFactors = TRUE)
  table <- read.csv(fileName,header=TRUE,dec=dec,sep=sep)
  if(is.null(factorName)) factor <- names(which(!vapply(table,is.numeric,FUN.VALUE = vector(length = 1))))[1]
  else if(!any(factorName==names(which(!vapply(table,is.numeric,FUN.VALUE = vector(length = 1)))))) stop(paste("The factorName",factorName," doesn't exists in dataframe"))
  else factor <- factorName
  labels<-table[,factor]
  if(is.null(classes)) classes<-levels(labels)
  else if(length(classes)<2) stop(paste("Just one class was indicated. More than one class of",factorName,"has to be indicated"))
  if(any(!c(classes %in% levels(labels)))) if(sum(!c(classes %in% levels(labels)))==1) stop(paste("Class",classes[!c(classes %in% levels(labels))],"isn't contained in",factorName))
  else stop(paste("Classes",paste(classes[!c(classes %in% levels(labels))],collapse = ", "),"aren't contained in",factorName))
  l <- array(NA, length(labels))
  l2 <- array(NA, length(labels))
  symbols <- unique(labels)
  i<-vector(length=length(classes))
  for(p in seq_len(length(i))) i[p] <- which(symbols == classes[p])
  
  symbols <- unique(labels)
  j<-list()
  v<-0
  for(p in seq_len(length(i))){
    j[[p]] <- which(labels == symbols[i[p]])
    l[j[[p]]] <- v
    l2[j[[p]]] <- as.character(symbols[v+1])
    v <- v+1
  }
  l[is.na(l)] <- -1
  l2[is.na(l2)] <- "Not selected"
  return(data.frame(code=l,names=l2))
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
#' @examples
#' # Read example gmt file
#' gmt_fname <- system.file("extdata", "c2.cp.v5.2.symbols.gmt", package = "BioNetStat")
#' deneSets <- readGmtFile(gmt_fname)
#' @export
readGmtFile <- function(fileName) {
  lines <- readLines(fileName)
  geneSets <- strsplit(lines, "\t")
  n <- length(geneSets)
  for (i in seq_len(n)) {
    if (length(geneSets[[i]]) < 3)
      stop(paste("Wrong variable sets file format. All sets should",
                 "contain at least one variable"))
  }
  return(geneSets)
}

################################################################################################################
#' Differential network analysis method
#'
#' @param method a function that receives two adjacency matrices and returns a list containing a statistic theta that measures the difference between them, and a p-value for the test H0: theta = 0 against H1: theta > 0.
#' @param options a list contaning parameters used by 'method'. Used only in degreeDistributionTest, spectralEntropyTest and spectralDistributionTest functions. It can be set to either \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @param varFile a numeric matrix contaning variables values data.
#' @param labels a vector of -1s, 0s, and 1s associating each sample with a phenotype. The value 0 corresponds to the first phenotype class of interest, 1 to the second phenotype class of interest, and -1 to the other classes, if there are more than two classes in the gene expression data.
#' @param varSets a list of gene sets. Each element of the list is a character vector v, where v[1] contains the gene set name, v[2] descriptions about the set, v[3..length(v)] the genes that belong to the set.
#' @param adjacencyMatrix a function that receives a numeric matrix containing gene expression data and returns the adjacency matrix of the inferred co-expression graph.
#' @param numPermutations the number of permutations for the permutation test.
#' @param print a logical. If true, it prints execution messages on the screen. resultsFile: path to a file where the partial results of the analysis will be saved. If NULL, then no partial results are saved.
#' @param resultsFile a ".RData" file name to be saved in tha work directory.
#' @param seed the seed for the random number generators. If it is not null then the sample permutations are the same for all the gene sets.
#' @param min.vert lower number of nodes (variables) that has to be to compare the networks.
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation, or a list of BiocParallelParam instances, to be applied in sequence for nested calls to BiocParallel functions. #MulticoreParam()
#' @param na.rm remove the NA values by excluding the rows ("row") or the columns ("col") that contaings it. If NULL (default) the NA values are not removed.
#' @return a data frame containing the name, size, test statistic, nominal p-value and adjusted p-value (q-value) associated with each gene set.
#' @examples 
#' # Glioma data
#' data("varFile")
#' gliomaData <- system.file("extdata", "bnsData.csv", package = "BioNetStat")
#' labels<-doLabels(gliomaData)
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue", threshold="fdr",
#'  thr.value=0.05, weighted=FALSE)
#' diffNetAnalysis(method=degreeCentralityTest, varFile=varFile, labels=labels, varSets=NULL,
#'  adjacencyMatrix=adjacencyMatrix1, numPermutations=1, print=TRUE, resultsFile=NULL,
#'   seed=NULL, min.vert=5, option=NULL)
#' # The numPermutations number is 1 to do a faster example, but we advise to use unless 1000 permutations in real analysis
#'   
#' # Random data
#' set.seed(1)
#' varFile <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-data.frame(code=rep(0:3,10),names=rep(c("A","B","C","D"),10))
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue", threshold="fdr",
#'  thr.value=0.05, weighted=FALSE)
#' diffNetAnalysis(method=degreeCentralityTest, varFile=varFile, labels=labels, varSets=NULL,
#'  adjacencyMatrix=adjacencyMatrix1, numPermutations=10, print=TRUE, resultsFile=NULL,
#'   seed=NULL, min.vert=5, option=NULL)
#' @export
#'
diffNetAnalysis <- function(method, options=list(bandwidth="Sturges"), varFile, labels, varSets=NULL,
                            adjacencyMatrix, numPermutations=1000, print=TRUE, resultsFile=NULL, 
                            seed=NULL, min.vert=5,BPPARAM=NULL, na.rm=NULL) {
  if(!is.data.frame(labels)) stop(("The object 'labels' is not a data frame. It has to be a dataframe with two columns (code and names)."))
  if(ncol(labels)!=2) stop(("Number of columns of 'labels' object is different from two. It has to be a two columns object (code and names)."))
  if(!is.null(na.rm)){ # Removing NA's
    if(na.rm=="col") {
      if(sum(apply(expr,2,function(x)any(is.na(x))))==1) print(paste(sum(apply(expr,2,function(x)any(is.na(x)))),"column was removed"))
      if(sum(apply(expr,2,function(x)any(is.na(x))))>1) print(paste(sum(apply(expr,2,function(x)any(is.na(x)))),"columns were removed"))
      varFile<-varFile[,!apply(varFile,2,function(x)any(is.na(x)))]
      }
    if(na.rm=="row") {
      if(sum(apply(expr,1,function(x)any(is.na(x))))==1) print(paste(sum(apply(expr,1,function(x)any(is.na(x)))),"row was removed"))
      if(sum(apply(expr,1,function(x)any(is.na(x))))>1) print(paste(sum(apply(expr,1,function(x)any(is.na(x)))),"rows were removed"))
      labels<-labels[!apply(varFile,1,function(x)any(is.na(x))),]
      varFile<-varFile[!apply(varFile,1,function(x)any(is.na(x))),]
    }
  }
  if(is.null(varSets)) varSets <- list(c("all",colnames(varFile)))
  if(!is.list(varSets)) stop("The varSet object is not a list. Please, insert the sets of variables as a list object")
  # else varSets[[length(varSets)+1]] <- c("all",colnames(varFile))
  varSets<-lapply(varSets, function(x){
    if(sum(x %in% colnames(varFile))>=min.vert) x
    else NA})#
  varSets<-varSets[!is.na(varSets)]
  results <- data.frame(matrix(NA, nrow=length(varSets), ncol=5+length(unique(labels$code[labels$code!=-1]))))
  output<-list()
  names <- array(NA, length(varSets))
  for (i in seq_len(length(varSets))) {
    names[i] <- varSets[[i]][1]
  }
  rownames(results) <- names
  if(any(labels$code=="-1"))colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste(unique(labels$names[-which(labels$code=="-1")])))
  else colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "Q-value",paste(unique(labels$names)))
  #temp <- tempfile("CoGA_results", fileext=".txt")
  for (i in seq_len(length(varSets))) {
    setName <- varSets[[i]][1]
    if (print)
      cat(paste("Testing ", setName, " (", i, " of ", length(varSets),
                ")", "\n", sep=""), append=TRUE)
    genes <- varSets[[i]][varSets[[i]] %in% colnames(varFile)]
    if (!is.null(seed))
      set.seed(seed)
    result <- method(varFile[,genes], labels$code, adjacencyMatrix=
                       adjacencyMatrix,  numPermutations=
                       numPermutations, options=options,BPPARAM=BPPARAM)
    if(!is.list(result)){
      for(c in colnames(result[,-c(1:3)])) colnames(result)[colnames(result)==c]<-labels$names[labels$code==c][1]
      output[[setName]]<-result
      results<-output
    }
    if(is.list(result)){
      results[setName, "N of Networks"] <- sum(unique(labels$code)!=-1)
      results[setName, "Test statistic"] <- result[[1]]
      results[setName, "Nominal p-value"] <- result$p.value
      results[setName, "Set size"] <- length(genes)
      results[setName, 6:ncol(results)] <- result$Partial
    }
    # if (!is.null(resultsFile)) save(results, file=resultsFile)
  }
  if(is.list(result)) results[, "Q-value"] <- p.adjust(results[, "Nominal p-value"], method="fdr")
  return(results)
}
