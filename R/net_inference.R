#' Adjacency matrix
#' @description creates a function that infers a graph from variables values matrix
#' @param method a function that measures the association between the variables values.
#' @param association a charactere string indicating wich value will be used as association value. The options are "corr" for the correlation value, "pvalue" for nominal pvalue associated to correlation or "fdr" for corrected pvalue for mutiple tests.
#' @param threshold a charactere string indicating wich value will be used as threshold value. The options are "corr" for the correlation value, "pvalue" for nominal pvalue associated to correlation or "fdr" for corrected pvalue for mutiple tests. If NULL, no edge is removed.
#' @param thr.value a numeric value. The function removes all edges weighted by a value less
#' than or equal to 'thr.value'.
#' @param weighted a logical value. If TRUE, then the edges of the graph are weighted by the
#' association degrees between the variables. Otherwise, the edges are are weighted by one.
#' @param abs.values a logical value. If TRUE, then the negatives edges of the graph are changed by its absolutes values. Otherwise, the negative edges are kept with negative weights.
#' @return a function that creates an adjacency matrix from variable values data.
#' @examples set.seed(3)
#' expr <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-rep(0:3,10)
#' functionAdjacencyMatrix <- adjacencyMatrix(method="spearman", association="pvalue",
#'  threshold="fdr", thr.value=0.05, weighted=FALSE)
#' @importFrom Hmisc rcorr
#' @importFrom psych corr.test
#' @export

adjacencyMatrix <- function(method, association="none", threshold="none",
                            thr.value=0.05, weighted=TRUE,abs.values=TRUE) {# Lembrar que o threshold pode ser para a correlacao ou para p-valor
  return(
    function(expr) {
      if(all(association!=c("fdr", "pvalue","corr","none"))) stop("Choose association argument as one of these options: corr, pvalue, fdr or none")
      if(!weighted) association<-"none"
      if(weighted & association=="none") stop("As you choose compare weigthed networks, you have to set association argument as one of these options: corr, pvalue, fdr")
      if(all(threshold!=c("fdr", "pvalue","corr","none"))) stop("Choose threshold argument as one of these options: corr, pvalue, fdr or none")
      if(association=="fdr" & (threshold=="pvalue" | threshold=="corr")) stop("In weighted networks, as you choose fdr as association value, you can not choose correlation or pvalue as threshold")
      if(any(method==c("pearson","spearman","kendall"))){
        if(nrow(expr)==3 | nrow(expr)==4) A <- list(cor(as.matrix(expr), method=method))
          else{
            if(any(method==c("pearson","spearman"))){
              A <- rcorr(as.matrix(expr), type=method)
            }
            if(method=="kendall"){
              A <- corr.test(as.matrix(expr), method="kendall", adjust="none")
              A$P <- A$p
            }
          }
        if(association=="fdr" | threshold=="fdr") A$FDR<-matrix(p.adjust(A$P,method = "fdr"),ncol(A[[1]]))
        if(association == "corr" | association=="none") AM <- A[[1]]
          else if(nrow(expr)==3 | nrow(expr)==4) stop("The method do not calculate p-value nor fdr for less than 5 samples, please choose association=corr ")
              else if(association=="pvalue") AM <- 1-A$P
                      else AM<- 1-A$FDR
        if (threshold=="corr") AM[A[[1]]>-(thr.value) & A[[1]]<(thr.value) | A[[1]]=="NaN" | is.na(A[[1]])]<-0
        if (threshold=="pvalue" | threshold=="fdr"){
          if(nrow(expr)==3 | nrow(expr)==4) stop("The method do not calculate p-value nor fdr for less than 5 samples, please choose threshold=corr ")
          if(threshold=="fdr") AM[A$FDR>=thr.value | is.na(A$FDR)]<-0
            else AM[A$P>=thr.value | is.na(A$P)]<-0
        }
        if(abs.values) AM <- abs(AM)
        diag(AM)<-0
        if (!weighted) AM[which(AM !=0)] <- 1
        AM[is.na(AM)]<-0
      }
        else{
          AM<-method(expr)
          if(!is.matrix(AM) | ncol(AM)!=nrow(AM) | isSymmetric(AM)) stop("The method inserted have to output a quadratic and simetrical matrix ")
          if(!is.null(thr.value)){
            if(!weighted) AM[which(AM > thr.value)] <- 1
            AM[which(AM <= thr.value)] <- 0
          }
          AM[is.na(AM)]<-0
        }
      return(AM)
    }
  )
}
