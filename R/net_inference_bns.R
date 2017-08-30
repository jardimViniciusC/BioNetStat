#' Run BNS on the browser user interface.
#' @export
runBioNetStat <- function(name='BioNetStat'){
  library(shiny)
  runApp(system.file('shiny', package='BioNetStat'))
}

#' Adjacency matrix
#' @description creates a function that infers a graph from variables values matrix
#' @param method: a function that measures the association between the variables values.
#' @param association: a charactere string indicating wich value will be used as association value. The options are "corr" for the correlation value, "pvalue" for nominal pvalue associated to correlation or "fdr" for corrected pvalue for mutiple tests.
#' @param threshold: a charactere string indicating wich value will be used as threshold value. The options are "corr" for the correlation value, "pvalue" for nominal pvalue associated to correlation or "fdr" for corrected pvalue for mutiple tests. If NULL, no edge is removed.
#' @param thr.value: a numeric value. The function removes all edges weighted by a value less
#' than or equal to 'thr.value'.
#' @param weighted:a logical value. If TRUE, then the edges of the graph are weighted by the
#' association degrees between the variables. Otherwise, the edges are are weighted by one.
#' @return a function that creates an adjacency matrix from variable values data.
#' @export


adjacencyMatrix <- function(method, association=c("corr","pvalue","fdr"), threshold=c("fdr", "pvalue","corr", NULL),
                            thr.value=0.05, weighted=T,abs.values=T) {# Lembrar que o threshold pode ser para a correlacao ou para p-valor
  return(
    function(expr) {
      if(any(method==c("pearson","spearman","kendall"))){
        if(nrow(expr)==3 | nrow(expr)==4) A <- list(cor(as.matrix(expr), method=method))
          else{
            if(any(method==c("pearson","spearman"))){
              A <- Hmisc::rcorr(as.matrix(expr), type=method) #tem que estar como matriz, não le data frames.
            }
            if(method=="kendall"){
              A <- psych::corr.test(as.matrix(expr), method="kendall", adjust="none")
              A$P <- A$p
            }
          }
        if (threshold=="corr"){
          A[[1]][A[[1]]>-(thr.value) & A[[1]]<(thr.value) | A[[1]]=="NaN" | is.na(A[[1]])]<-0
          diag(A[[1]])<-0
          A$P[A[[1]]>-(thr.value) & A[[1]]<(thr.value) | A[[1]]=="NaN"]<-1
          A$P <- 1 - A$P
        }
        if (threshold=="pvalue" | threshold=="fdr"){
          if(nrow(expr)==3 | nrow(expr)==4) stop("The method do not calculate p-value nor fdr for less than 5 samples, please choose threshold=corr ")
          if(threshold=="fdr")  A$P<-matrix(p.adjust(A$P,method = "fdr"),ncol(A[[1]]))
          A[[1]][A$P>=thr.value | is.na(A$P)]<-0
          A$P[A$P>=thr.value | is.na(A$P)]<-1
          diag(A$P) <- 1
          A$P <- 1 - A$P
        }
        if(abs.values==T) A[[1]] <- abs(A[[1]])

        if(association == "corr") A <- A[[1]]
        else
          if(nrow(expr)==3 | nrow(expr)==4) stop("The method do not calculate p-value nor fdr for less than 5 samples, please choose association=corr ")
            else A <- A$P

        if (!weighted) {
          A[which(A > thr.value)] <- 1
          A[which(A <= thr.value)] <- 0
        }
        A[is.na(A)]<-0
      }
        else{
          A<-method(expr)
          if(!is.matrix(A) | ncol(A)!=nrow(A) | isSymmetric(A)) stop("The method inserted have to output a quadratic and simetrical matrix ")
          if(!is.null(thr.value)){
            if(!weighted) A[which(A > thr.value)] <- 1
            A[which(A <= thr.value)] <- 0
          }
          A[is.na(A)]<-0
        }
      return(A)
    }
  )
}
