# ------------------------------------------------------------------------------
# Graph features
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

clusterCoef <- function(A) {
  s <- sum(A)
  results <- c()
  n <- nrow(A)
  d <- (n-1)*(n-2)
  results<-apply(A,1,function(x) (s - 2*sum(x))/d)
  return(results)
}

# Inverting Matrix Weights
invWeigthts<-function(A){
  A[which(A > 0)] <- 1/A[which(A > 0)]
  return(A)
}

# Density functions of the degrees of two graphs
#
# 'degreeDensities' estimates the density functions of the degrees for two
# graphs at the same coordinates
# @param G1 an igraph graph object
# @param G2 an igraph graph object
# @param npoints number of points used in density function estimation
# @param options a list containing parameters. It can be set to either
# \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
# @return a list containing the components, f1 (density estimate of the
# graph G1), and f2 (density estimate of the graph G2). Each component is
# a list, where the first element is the vector 'x' of 'npoints' coordinates
# of the points where the density function is estimated, and the second is
# a vector 'y' of the estimated density values.
# @examples G1<-erdos.renyi.game(30,0.6)
# G2<-barabasi.game(30,power = 1)
# d<-degreeDensities(G1, G2, npoints=1024, options=list(bandwidth="Sturges"))
# par(mfrow=c(1,2))
# plot(d$f1$x,d$f1$y,main="Erdos-Renyi\n Degree distribution",
# xlab="Degree",ylab="Frequency")
# plot(d$f2$x,d$f2$y,main="Barabasi\n Degree distribution",
# xlab="Degree",ylab="Frequency")
degreeDensities <- function(G1, G2, npoints=1024, options=list(bandwidth="Sturges")) {
  n1 <- vcount(G1)
  n2 <- vcount(G2)
  e1 <- graph.strength(G1)
  e2 <- graph.strength(G2)
  from <- min(e1, e2)
  to <- max(e1, e2)
  f1 <- gaussianDensity(e1, from=from, to=to, bandwidth=options$bandwidth, npoints=npoints)
  f2 <- gaussianDensity(e2, from=from, to=to, bandwidth=options$bandwidth, npoints=npoints)
  if (sum(is.na(f1)) > 0 || sum(is.na(f2)) > 0)
    return(NA)
  return(list("f1"=f1, "f2"=f2))
}

#' Density functions of the degrees of n graphs
#'
#' 'nDegreeDensities' estimates the density functions of the degrees for n
#' graphs at the same coordinates
#' @param Gs a list of n igraph graphs objects
#' @param npoints number of points used in density function estimation
#' @param bandwidth a parameters. It can be set to either "Sturges" or "Silverman".
#' @param from the lower value used to build the distribution
#' @param to the higher value used to build the distribution
#' @return a list containing the components 'x' and 'densities'. 
#' The first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function i estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @examples G<-list()                                   
#' G[[1]]<-erdos.renyi.game(30,0.6)
#' G[[2]]<-barabasi.game(30,power = 1)
#' G[[3]]<-watts.strogatz.game(2,30,2,0.3)
#' d<-nDegreeDensities(G, npoints=1024, bandwidth="Sturges")
#' par(mfrow=c(1,3))
#' plot(d$x,d$densities[,1],main="Erdos-Renyi\n Degree distribution",
#' xlab="Degree",ylab="Frequency")
#' plot(d$x,d$densities[,2],main="Barabasi\n Degree distribution",
#' xlab="Degree",ylab="Frequency")
#' plot(d$x,d$densities[,3],main="Watts-Strogatz\n Degree distribution",
#' xlab="Degree",ylab="Frequency")
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
nDegreeDensities <- function(Gs, npoints=1024, bandwidth="Sturges",from=NULL,to=NULL) {
  e<-lapply(Gs,graph.strength)
  densities <- matrix(NA, npoints, length(Gs))
  if(is.null(from) || is.null(to)){
    from <- min(unlist(e))
    to <- max(unlist(e))
  }
  for(i in seq_len(length(Gs))){
    f <- gaussianDensity(e[[i]], from=from, to=to, bandwidth=bandwidth, npoints=npoints)
    if (any(is.na(f))) return(NA)
    densities[,i]<-f$y
    x<-f$x
  }
  if (sum(is.na(x)) > 0 || sum(is.na(densities)) > 0)
    return(NA)
  return(list("x"=x, "densities"=densities))
}

#' Kullback-Liebler divergence among the density functions of the degrees of
#' two or more graphs
#'
#' 'KLdegree' computes the Kullback-Liebler divergence among the density
#' functions of the degrees of two or more graphs
#'
#' @param f a list containing the components 'x' and 'densities'. 
#' The first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function i estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @return returns a list containing the components 'theta' and 'partial'. 
#' 'theta' is a value representaing the Kullback-Liebler divergence among the corresponding distributions.
#' 'partial' is a vector of KL divergences between each network distribuiton and the average degree distribution.
#' @examples G<-list()                                   
#' G[[1]]<-erdos.renyi.game(30,0.6)
#' G[[2]]<-barabasi.game(30,power = 1)
#' G[[3]]<-watts.strogatz.game(2,30,2,0.3)
#' f<-nDegreeDensities(G, npoints=1024, bandwidth="Sturges")
#' KLdegree(f)
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
KLdegree<-function(f){
  if(any(is.na(f))){
    cat('Empty graph')
    return(list("measure"=NA, "p.value"=NA,"Partial"=NA))
  }
  meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
  partial <- vector(length=ncol(f$densities))
  for (j in seq_len(ncol(f$densities))) {
    f1 <- list("x"=f$x, "y"=f$densities[,j])
    partial[j] <- KL(f1, meanDensity)/ncol(f$densities)
  }
  return(list(theta=sum(partial),Partial=partial))
}

# Jensen-Shannon divergence between the density functions of the degrees of
# two graphs
#
# 'JSdegree' computes the Jensen-Shannon divergence between the density
# functions of the degrees of two graphs
#
# @param G1 an igraph graph object
# @param G2 an igraph graph object
# @param options a list containing parameters. It can be set to either
# \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
# @return a list containing the components, f1 (density estimate of the
# graph G1), and f2 (density estimate of the graph G2). Each component is
# a list, where the first element is the vector 'x' of 'npoints' coordinates
# of the points where the density function i estimated, and the second is
# a vector 'y' of the estimated density values.
JSdegree <- function(G1, G2, options=list(bandwidth="Sturges")) {
  f <- degreeDensities(G1, G2, options=options)
  if (sum(is.na(f)) > 0)
    return(NA)
  f1 <- f$f1
  f2 <- f$f2
  fm <- f1
  fm$y <- (f1$y + f2$y)/2
  return((KL(f1, fm) + KL(f2, fm))/2)
}

# Returns the spectral density for a given adjacency matrix A
spectralDensity <- function(A, bandwidth="Sturges", npoints=1024) {
  eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/sqrt(nrow(A)))
  return(gaussianDensity(eigenvalues, bandwidth=bandwidth, npoints=npoints))
}

# Returns the spectral densities for given adjacency matrices A1 and A2 at the
# same points
spectralDensities <- function(A1, A2, bandwidth="Sturges",
                              npoints=1024) {
  n1 <- nrow(A1)
  n2 <- nrow(A2)
  e1 <- (as.numeric(eigen(A1, only.values = TRUE)$values)/sqrt(n1))
  e2 <- (as.numeric(eigen(A2, only.values = TRUE)$values)/sqrt(n2))
  from <- min(e1, e2)
  to <- max(e1, e2)
  f1 <- gaussianDensity(e1, from=from, to=to, bandwidth=bandwidth, npoints=npoints)
  f2 <- gaussianDensity(e2, from=from, to=to, bandwidth=bandwidth, npoints=npoints)
  if (sum(is.na(f1)) > 0 || sum(is.na(f2)) > 0)
    return(NA)
  return(list("f1"=f1, "f2"=f2))
}

#' Spectral Density functions of n graphs
#'
#' Returns the spectral densities for a list of adjacency matrices at the
#' same points
#' @param A a list of adjacency matrices
#' @param from the lower value used to build the distribution
#' @param to the higher value used to build the distribution
#' @param bandwidth a parameters. It can be set to either "Sturges" or "Silverman".
#' @return a list containing the components 'x' and 'densities'. 
#' The first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function i estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @examples A<-list()
#' A[[1]]<-as.matrix(as_adj(erdos.renyi.game(30,0.6,directed = FALSE)))
#' A[[2]]<-as.matrix(as_adj(barabasi.game(30,power = 1,directed = FALSE)))
#' A[[3]]<-as.matrix(as_adj(watts.strogatz.game(1,30,2,0.3)))
#' d<-nSpectralDensities(A, bandwidth="Sturges")
#' par(mfrow=c(1,3))
#' plot(d$x,d$densities[,1],main="Erdos-Renyi\n Spectral distribution",
#' xlab="Eigenvalue",ylab="Frequency")
#' plot(d$x,d$densities[,2],main="Barabasi\n Spectral distribution",
#' xlab="Eigenvalue",ylab="Frequency")
#' plot(d$x,d$densities[,3],main="Watts-Strogatz\n Spectral distribution",
#' xlab="Eigenvalue",ylab="Frequency")
#' @seealso \code{KLdegree}
#' @seealso \code{density}
#' @import igraph
#' @export
nSpectralDensities <- function (A, from=NULL, to=NULL, bandwidth="Silverman") {
  npoints <- 1024
  ngraphs <- length(A)
  n <- ncol(A[[1]])
  spectra <- matrix(NA, n, ngraphs)
  for (i in seq_len(ngraphs)) {
    Adj <- A[[i]]
    eigenvalues <- (as.numeric(eigen(Adj, only.values = TRUE)$values)/
                      sqrt(nrow(Adj)))
    spectra[,i] <- eigenvalues
  }
  densities <- matrix(NA, npoints, ngraphs)
  minimum <- min(spectra)
  maximum <- max(spectra)
  if (!is.null(from) && !is.null(to)) {
    minimum <- from
    maximum <- to
  }
  for (i in seq_len(ngraphs)) {
    f <- gaussianDensity(spectra[,i], bandwidth=bandwidth,
                         from=minimum, to=maximum,
                         npoints=npoints)

    densities[,i] <- f$y
    x <- f$x
  }
  return(list("x"=x, "densities"=densities))
}

#' Kullback-Liebler divergence among the spectral density functions of 
#' two or more graphs
#'
#' 'KLspectrum' computes the Kullback-Liebler divergence among the spectral density
#' functions of two or more graphs
#'
#' @param f a list containing the components 'x' and 'densities'. 
#' The first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function i estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @return returns a list containing the components 'theta' and 'partial'. 
#' 'theta' is a value representaing the Kullback-Liebler divergence among the corresponding distributions.
#' 'partial' is a vector of KL divergences between each network distribuiton and the average spectral distribution.
#' @examples A<-list()
#' A[[1]]<-as.matrix(as_adj(erdos.renyi.game(30,0.6,directed = FALSE)))
#' A[[2]]<-as.matrix(as_adj(barabasi.game(30,power = 1,directed = FALSE)))
#' A[[3]]<-as.matrix(as_adj(watts.strogatz.game(1,30,2,0.3)))
#' f<-nSpectralDensities(A, bandwidth="Sturges")
#' KLspectrum(f)
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
KLspectrum<-function(f){
  meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
  partial <- vector(length = ncol(f$densities))
  for (j in seq_len(ncol(f$densities))) {
    f1 <- list("x"=f$x, "y"=f$densities[,j])
    partial[j] <- KL(f1, meanDensity)/ncol(f$densities)
  }
  return(list(theta=sum(partial),Partial=partial))
}

# Given two adjacency matrices, returns the Jensen-Shannon divergence between
# the corresponding graphs
#
# 'JSspectrum' computes the Jensen-Shannon divergence between the spectral density
# functions of two graphs
#
# @param G1 an adjacency matrix
# @param G2 an adjacency matrix
# @param bandwidth a parameters. It can be set to either "Sturges" or "Silverman".
# @return returns a value representaing the Jensen-Shannon divergence between the corresponding graphs
JSspectrum <- function(A1, A2, bandwidth="Sturges") {
  f <- spectralDensities(A1, A2, bandwidth=bandwidth)
  if (sum(is.na(f)) > 0)
    return(NA)
  f1 <- f$f1
  f2 <- f$f2
  fm <- f1
  fm$y <- (f1$y + f2$y)/2
  return((KL(f1, fm) + KL(f2, fm))/2)
}

# Given two spectral densities, returns the Jensen-Shannon divergence between
# the corresponding graphs
#
# 'JS' computes the Jensen-Shannon divergence between the spectral density
# functions of two graphs
# 
JS <- function(f1, f2) {
  fm <- f1
  fm$y <- (f1$y + f2$y)/2
  return((KL(f1, fm) + KL(f2, fm))/2)
}


# Returns the absolute difference of the adjacency matrix A1 and A2 spectral
# entropies
absDiffSpectralEntropy <- function(A1, A2, bandwidth="Sturges") {
  fs <- spectralDensities(A1, A2, bandwidth=bandwidth)
  if (sum(is.na(fs)) > 0)
    return(NA)
  H1 <- entropy(fs$f1)
  H2 <- entropy(fs$f2)
  return(abs(H1 - H2))
}

# ------------------------------------------------------------------------------
# Node scores
# ------------------------------------------------------------------------------

#' Node score
#' @description Node score (degree, betweenness, closenness, eigenvector centralities or clustering coefficient) for each network analysed.
#' @param expr Matrix of variables (columns) vs samples (rows).
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state.
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix.
#' @return a list of vector containing the node scores (degree, betweenness, closenness, eigenvector centralities or clustering coefficient) of each network.
#' @examples 
#' set.seed(1)
#' expr <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-rep(0:3,10)
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue",
#'  threshold="fdr", thr.value=0.05, weighted=FALSE)
#' @name nodeScore

#' @rdname nodeScore
#' @examples
#' 
#' # Degree centrality
#' degreeCentrality(expr, labels, adjacencyMatrix1)
#' @export
degreeCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, graph.strength)
  return(result)
}

#' @rdname nodeScore
#' @examples 
#' 
#' # Betweenness Centrality
#' betweennessCentrality(expr, labels, adjacencyMatrix1)
#' @export
betweennessCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, betweenness)
  return(result)
}

#' @rdname nodeScore
#' @examples 
#' 
#' # Closenness Caentrality
#' closenessCentrality(expr, labels, adjacencyMatrix1)
#' @export
closenessCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, closeness)
  return(result)
}

#' @rdname nodeScore
#' @examples 
#' 
#' # Eigenvector centrality
#' eigenvectorCentrality(expr, labels, adjacencyMatrix1)
#' @export
eigenvectorCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, function(x) evcent(x)$vector)
  return(result)
}

#' @rdname nodeScore
#' @examples
#' 
#' # Clustering coefficient
#' clusteringCoefficient(expr, labels, adjacencyMatrix1)
#' @export
clusteringCoefficient <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) result <- lapply(A, clusterCoef)
  else {
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    result <- lapply(G, transitivity, type="local", isolates="zero")
  }
  return(result)
}

# ------------------------------------------------------------------------------
# Network features
# ------------------------------------------------------------------------------

#' Network features
#' @name networkFeature
#' @description Network feature average nodes scores (degree, betweenness, closenness, eigenvector centralities or clustering coefficient) or spectral entropies for each network analysed.
#' @param expr Matrix of variables (columns) vs samples (rows)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the spectral entropie or average node score of each network.
#' @examples set.seed(1)
#' expr <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-rep(0:3,10)
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue",
#'  threshold="fdr", thr.value=0.05, weighted=FALSE)


#' @rdname networkFeature
#' @examples 
#' 
#' # Average degree centrality
#' averageDegreeCentrality(expr, labels, adjacencyMatrix1)
#' @export
averageDegreeCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- degreeCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' @rdname networkFeature
#' @examples 
#' 
#' # Average betweenness centrality
#' averageBetweennessCentrality(expr, labels, adjacencyMatrix1)
#' @export
averageBetweennessCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- betweennessCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' @rdname networkFeature
#' @examples 
#' 
#' # Average closeness centrality
#' averageClosenessCentrality(expr, labels, adjacencyMatrix1)
#' @export
averageClosenessCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- closenessCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' @rdname networkFeature
#' @examples 
#' 
#' # Average eigenvector centrality
#' averageEigenvectorCentrality(expr, labels, adjacencyMatrix1)
#' @export
averageEigenvectorCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- eigenvectorCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' @rdname networkFeature
#' @examples 
#' 
#' # Average clustering coefficient
#' averageClusteringCoefficient(expr, labels, adjacencyMatrix1)
#' @export
averageClusteringCoefficient <- function(expr, labels, adjacencyMatrix) {
  result <- clusteringCoefficient(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' @rdname networkFeature
#' @examples 
#' 
#' # Average shortest path
#' averageShortestPath(expr, labels, adjacencyMatrix1)
#' @export
averageShortestPath <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, average.path.length)
  return(result)
}

# Returns the spectral entropy for a given adjacency matrix A
spectralEntropy <- function(A, bandwidth="Sturges") {
  f <- spectralDensity(A, bandwidth=bandwidth)
  if (sum(is.na(f)) > 0)
    return(NA)
  y <- f$y
  n <- length(y)
  i <- which(y != 0)
  y[i] <- y[i]*log(y[i])
  return(-trapezoidSum(f$x, y))
}

#' @rdname networkFeature
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return spectralEntropies. A list of values containing the spectral entropy of each network.
#' @examples 
#' 
#' # Spectral entropies
#' spectralEntropies(expr, labels, adjacencyMatrix1, options=list(bandwidth="Sturges"))
#' @export
spectralEntropies <- function(expr, labels, adjacencyMatrix, options=list(bandwidth="Sturges")) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in seq_len(length(unique(labels)))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
  }
  fs<- nSpectralDensities(A,bandwidth=options$bandwidth)
  if (sum(is.na(fs)) > 0)
    return(NA)
  entropies<-list()
  for(j in seq_len(length(A))){ # uma entropia para cada grafo
    entropies[[j]]<-entropy(list("x"=fs$x, "y"=fs$densities[,j]))
  }
  return(entropies)
}

# ------------------------------------------------------------------------------
# Test of equality between network properties
# ------------------------------------------------------------------------------

resInt <- function(A,expr,weighted,fun){
  n <- ncol(expr)# numero de genes
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
  s<-lapply(G, fun)# guarda a centralidade
  s<-do.call(rbind,s) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
  partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))# calcula a distancia euclidiana entre cada vetor de betweenness e o vetor medio
  return(c(sum(partial),partial)) # A estatítica é a soma das distancias
}

#' Network equality test
#' @name networkTest
#' @description Test of equality between network properties
#' @param expr Matrix of variables (columns) vs samples (rows)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @param options a list containing parameters. Used only in degreeDistributionTest, spectralEntropyTest and spectralDistributionTest functions. It can be set to either \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation, or a list of BiocParallelParam instances, to be applied in sequence for nested calls to BiocParallel functions.
#' @return A list containing:
#' "measure" - difference among two or more networks associated with each phenotype. To compare networks by centralities and clustering coefficient, 
#' one uses euclidian distance. In spectral entropy comparison, one uses the absolute difference. In distributions (spectral and degree) comparison, 
#' one uses Kulback-Liebler divergence.
#' "p.value" - the Nominal p-value of the test.
#' "Partial" - a vector with the weigths of each network in a measure value.
#' @examples 
#' set.seed(1)
#' expr <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-rep(0:3,10)
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue",
#'  threshold="fdr", thr.value=0.05, weighted=FALSE)
#' # The numPermutations number is 1 to do a faster example, but we advise to use unless 1000 permutations in real analysis

#' @rdname networkTest
#' @examples 
#' 
#' # Degree centrality test
#' degreeCentralityTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
degreeCentralityTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  output<-resInt(A,expr,weighted,graph.strength)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE)
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    return(resInt(A,expr,weighted,graph.strength)[1])
  }, BPPARAM=BPPARAM)
  results<-do.call(c,results)
  pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
  return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Betweenness centrality test
#' betweennessCentralityTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
betweennessCentralityTest <- function(expr, labels, adjacencyMatrix,numPermutations=1000, options=NULL,BPPARAM=MulticoreParam()) {
  # Betweenness centrality test for many graphs

    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
    if(any(!v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    output<-resInt(A,expr,weighted,betweenness)
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE)
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      if (!is.null(weighted)) A<-lapply(A,invWeigthts)
      return(resInt(A,expr,weighted,betweenness)[1])
    }, BPPARAM=BPPARAM)
    results<-do.call(c,results)
    pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
    return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Closeness centrality test
#' closenessCentralityTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
closenessCentralityTest <- function(expr, labels, adjacencyMatrix,numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  # Closeness centrality test for many graphs
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
    if(any(!v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    output<-resInt(A,expr,weighted,closeness)
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE)
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      if (!is.null(weighted)) A<-lapply(A,invWeigthts)
      return(resInt(A,expr,weighted,closeness)[1])
    }, BPPARAM=BPPARAM)
    results<-do.call(c,results)
    pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
    return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Eigenvector centrality test
#' eigenvectorCentralityTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
eigenvectorCentralityTest <- function(expr, labels, adjacencyMatrix,numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  # Eigenvector centrality test for many graphs
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  output<-resInt(A,expr,weighted,function(x) evcent(x)$vector)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE)
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    return(resInt(A,expr,weighted,function(x) evcent(x)$vector)[1])
  }, BPPARAM=BPPARAM)
  results<-do.call(c,results)
  pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
  return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Clustering coefficient test
#' clusteringCoefficientTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
clusteringCoefficientTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
    if(any(!v)) weighted <- TRUE
    if (!is.null(weighted)) { 
      n <- ncol(expr)
      s<-lapply(A, clusterCoef)
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
      partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
      output<-c(sum(partial),partial)
    }
    else{
      output<-resInt(A,expr,weighted,function(x){transitivity(x,type="local", isolates="zero")})
    }
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE)
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      if (!is.null(weighted)){ 
        s<-lapply(A, clusterCoef)
        s<-do.call(rbind,s)
        s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
        res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
        return(sum(res))
      }
      else return(resInt(A,expr,weighted,function(x){transitivity(x,type="local", isolates="zero")}))
    }, BPPARAM=BPPARAM)
    results<-do.call(c,results)
    pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
    return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
}

# #' @rdname networkTest
# #' @examples
# #'
# #' # Shortest path test
# #' shortestPathTest(expr, labels, adjacencyMatrix1,numPermutations=1)
# #' @export
# shortestPathTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
#   # Shortest path test for many graphs
#   lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
#   if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
#   A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
#   weighted <- NULL # Define o weighted como NULL, assim como na função original
#   v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
#   if(any(!v)) weighted <- TRUE
#   if (!is.null(weighted)) A<-lapply(A,invWeigthts)
#   # if(is.null(weighted)) output<-resInt(A,expr,weighted,function(x){average.path.length(x,directed=FALSE)})
#   # else output<-resInt(A,expr,weighted,function(y){apply(distances(y), 1, function(x){ min(x[x!=0])})})
#   output<-resInt(A,expr,weighted,function(y){apply(distances(y), 1, function(x){ min(x[x!=0])})})
#   results<-bplapply(seq_len(numPermutations),function(i){
#     l <- sample(labels, replace = FALSE)
#     A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
#     if(is.null(weighted)){
#       A<-lapply(A,invWeigthts)
#       return(resInt(A,expr,weighted,function(x){average.path.length(x,directed=FALSE)})[1])
#     }
#     else return(resInt(A,expr,weighted,function(y){apply(distances(y), 1, function(x){ min(x[x!=0])})})[1])
#   }, BPPARAM=BPPARAM)
#   results<-do.call(c,results)
#   pvalue <- (1 + sum(results >= output[1]))/(numPermutations + 1)
#   return(list("measure"=output[1], "p.value"=pvalue,"Partial"=output[-1]))
# }

#' @rdname networkTest
#' @examples
#' 
#' # Degree distribution test
#' degreeDistributionTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
degreeDistributionTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=list(bandwidth="Silverman"), BPPARAM=MulticoreParam()) {

    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
    if(any(!v)) weighted <- TRUE
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
    f<-nDegreeDensities(Gs=G, bandwidth=options$bandwidth)
    result<-KLdegree(f)
    if(length(result$Partial)==1 & all(is.na(result$Partial))) result$Partial<-rep(NA,length(G))
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE)
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
      f<-nDegreeDensities(Gs=G, bandwidth=options$bandwidth)
      return(KLdegree(f)$theta)
    }, BPPARAM=BPPARAM)
    results<-do.call(c,results)
  pvalue <- (1 + sum(results >= result$theta))/(numPermutations + 1)
  return(list("measure"=result$theta, "p.value"=pvalue,"Partial"=result$Partial))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Spectral entropy test
#' spectralEntropyTest(expr, labels, adjacencyMatrix1,numPermutations=1,
#'  options=list(bandwidth="Sturges"))
#' @export
spectralEntropyTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=list(bandwidth="Sturges"), BPPARAM=MulticoreParam()) {
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
    entropies<-vector(length=length(A)) # vetor para guardar entropias
    for(j in seq_len(length(A))){ # uma entropia para cada grafo
      entropies[j]<-entropy(list("x"=f$x, "y"=f$densities[,j]))
    }
    meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities)) # Calcula a entropia média a partir de uma distribuicao media
    result<-sqrt((sum((entropies-entropy(meanDensity))^2))/length(entropies)) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
    partial<-sqrt(((entropies-entropy(meanDensity))^2)/length(entropies)) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE) # Reamostra os labels sem reposicao
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
      entropies<-vector(length=length(A)) # vetor para guardar entropias
      for(j in seq_len(length(A))){ # uma entropia para cada grafo
        entropies[j]<-entropy(list("x"=f$x, "y"=f$densities[,j]))
      }
      meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities)) # Calcula a entropia média a partir de uma distribuicao media
      return(sqrt((sum((entropies-entropy(meanDensity))^2))/length(entropies))) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
    }, BPPARAM=BPPARAM)
    results<-do.call(c,results)
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1) # calculo do pvalor
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' @rdname networkTest
#' @examples 
#' 
#' # Spectral distribution test
#' spectralDistributionTest(expr, labels, adjacencyMatrix1,numPermutations=1,
#'  options=list(bandwidth="Sturges"))
#' @export
spectralDistributionTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=list(bandwidth="Sturges"),BPPARAM=MulticoreParam()) {
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
    f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
    result<-KLspectrum(f)
    results<-bplapply(seq_len(numPermutations),function(i){
      l <- sample(labels, replace = FALSE)
      A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
      f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
      return(KLspectrum(f)$theta)
    }, BPPARAM = BPPARAM)
    results<-do.call(c,results)
  pvalue <- (1 + sum(results >= result$theta))/(numPermutations + 1)
  return(list("measure"=result$theta, "p.value"=pvalue,"Partial"=result$Partial))
}

# ------------------------------------------------------------------------------
# Test of equality of vertices between network properties
# ------------------------------------------------------------------------------

resVertexInt <- function(A,expr,weighted,fun){
  n <- ncol(expr)# numero de genes
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
  s<-lapply(G, fun)# guarda o betweenness
  s<-do.call(rbind,s) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
  sp<-s[-(nrow(s)),]
  result<-apply(s[-(nrow(s)),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,sum)
  return(cbind(result,t(sp))) # A estatítica é a soma das distancias
}

retTable <- function(results,output,expr,numPermutations,lab){
  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= output[,1]),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(output[,1],3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(output[,-1],3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}

#' Node score equality test
#' @name nodeTest
#' @description Nodes scores equality test between network
#' @param expr Matrix of variables (columns) vs samples (rows)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @param options argument non used in this function
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation, or a list of BiocParallelParam instances, to be applied in sequence for nested calls to BiocParallel functions.
#' @return A table, containing on the columns, the following informations for each variable (rows):
#' "Test Statistic" - difference among the degree centrality of a node in two or more networks associated with each phenotype
#' "Nominal p-value" - the Nominal p-value of the test
#' "Q-value" - the q-value of the test, correction of p-value by FDR to many tests
#' "Factor n" - the node degree centrality in each network compared
#' @examples 
#' set.seed(1)
#' expr <- as.data.frame(matrix(rnorm(120),40,30))
#' labels<-rep(0:3,10)
#' adjacencyMatrix1 <- adjacencyMatrix(method="spearman", association="pvalue",
#'  threshold="fdr", thr.value=0.05, weighted=FALSE)
#' # The numPermutations number is 1 to do a faster example, but we advise to use unless 1000 permutations in real analysis

#' @rdname nodeTest
#' @examples 
#' 
#' # Degree centrality test
#' degreeCentralityVertexTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
degreeCentralityVertexTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL,BPPARAM=MulticoreParam()) {
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  output<-resVertexInt(A,expr,weighted,graph.strength)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE)
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    return(resVertexInt(A,expr,weighted,graph.strength)[,1])
  }, BPPARAM=BPPARAM)
  return(retTable(results,output,expr,numPermutations,lab))
}

#' @rdname nodeTest
#' @examples 
#' 
#' # Betweenness centrality test
#' betweennessCentralityVertexTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
betweennessCentralityVertexTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  # Betweenness centrality test for many graphs
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  output<-resVertexInt(A,expr,weighted,betweenness)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    return(resVertexInt(A,expr,weighted,betweenness)[,1])
  }, BPPARAM=BPPARAM)
  return(retTable(results,output,expr,numPermutations,lab))
}

#' @rdname nodeTest
#' @examples 
#' 
#' # Closeness centrality test
#' closenessCentralityVertexTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
closenessCentralityVertexTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  output<-resVertexInt(A,expr,weighted,closeness)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    return(resVertexInt(A,expr,weighted,closeness)[,1])
  }, BPPARAM=BPPARAM)
  return(retTable(results,output,expr,numPermutations,lab))
}

#' @rdname nodeTest
#' @examples 
#' 
#' # Eigenvector centrality test
#' eigenvectorCentralityVertexTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
eigenvectorCentralityVertexTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  output<-resVertexInt(A,expr,weighted,function(x) evcent(x)$vector)
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE)
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    return(resVertexInt(A,expr,weighted,function(x) evcent(x)$vector)[,1])
  }, BPPARAM=BPPARAM)
  return(retTable(results,output,expr,numPermutations,lab))
}

#' @rdname nodeTest
#' @examples 
#' 
#' # Clustering coefficient test
#' clusteringCoefficientVertexTest(expr, labels, adjacencyMatrix1,numPermutations=1)
#' @export
clusteringCoefficientVertexTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000, options=NULL, BPPARAM=MulticoreParam()) {
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  A<-lapply(lab, function(x) adjacencyMatrix(expr[labels==x,]))
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  v<-vapply(A,FUN = function(x) sum(x==0) + sum(x==1) == length(x),FUN.VALUE = vector(length = 1))
  if(any(!v)) weighted <- TRUE
  if (!is.null(weighted)) { 
    n <- ncol(expr)
    s<-lapply(A, clusterCoef)
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    sp<-s[-(nrow(s)),]
    result<-apply(s[-(nrow(s)),],1,function(x) abs(x-s[dim(s)[1],]))
    result<-apply(result,1,sum)
    output<-cbind(result,t(sp))
  }
  else output<-resVertexInt(A,expr,weighted,function(x){transitivity(x,type="local", isolates="zero")})
  results<-bplapply(seq_len(numPermutations),function(i){
    l <- sample(labels, replace = FALSE)
    A<-lapply(lab, function(x) adjacencyMatrix(expr[l==x,]))
    if (!is.null(weighted)) { 
      n <- ncol(expr)
      s<-lapply(A, clusterCoef)
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
      res<-apply(s[seq_len(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
      return(apply(res,1,sum))
    }
    else return(resVertexInt(A,expr,weighted,function(x){transitivity(x,type="local", isolates="zero")})[,1])
  }, BPPARAM=BPPARAM)
  return(retTable(results,output,expr,numPermutations,lab))
}
