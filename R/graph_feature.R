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

#' Density functions of the degrees of two graphs
#'
#' 'degreeDensities' estimates the density functions of the degrees for two
#' graphs at the same coordinates
#' @param G1 an igraph graph object
#' @param G2 an igraph graph object
#' @param npoints number of points used in density function estimation
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return a list containing the components, f1 (density estimate of the
#' graph G1), and f2 (density estimate of the graph G2). Each component is
#' a list, where the first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function is estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
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
#' 'degreeDensities' estimates the density functions of the degrees for n
#' graphs at the same coordinates
#' @param Gs a list of n igraph graphs objects
#' @param npoints number of points used in density function estimation
#' @param bandwidth a parameters. It can be set to either "Sturges" or "Silverman".
#' @param from the lower value used to build the distribution
#' @param to the higher value used to build the distribution
#' @return a list containing the components, f1 (density estimate of the
#' graph G1), and f2 (density estimate of the graph G2). Each component is
#' a list, where the first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function is estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
nDegreeDensities <- function(Gs, npoints=1024, bandwidth="Sturges",from=NULL,to=NULL) {
  e<-list()
  for(i in 1:length(Gs)){
    e[[i]] <- graph.strength(Gs[[i]])
  }
  densities <- matrix(NA, npoints, length(Gs))
  if(is.null(from) || is.null(to)){
    from <- min(unlist(e))
    to <- max(unlist(e))
  }
  for(i in 1:length(Gs)){
    f <- gaussianDensity(e[[i]], from=from, to=to, bandwidth=bandwidth, npoints=npoints)
    if (any(is.na(f))) return(NA)
    densities[,i]<-f$y
    x<-f$x
  }
  if (sum(is.na(x)) > 0 || sum(is.na(densities)) > 0)
    return(NA)
  return(list("x"=x, "densities"=densities))
}

#' Jensen-Shannon divergence between the density functions of the degrees of
#' two graphs
#'
#' 'JSdegree' computes the Jensen-Shannon divergence between the density
#' functions of the degrees of two graphs
#'
#' @param G1 an igraph graph object
#' @param G2 an igraph graph object
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return a list containing the components, f1 (density estimate of the
#' graph G1), and f2 (density estimate of the graph G2). Each component is
#' a list, where the first element is the vector 'x' of 'npoints' coordinates
#' of the points where the density function i estimated, and the second is
#' a vector 'y' of the estimated density values.
#' @seealso \code{graph.strength}
#' @seealso \code{density}
#' @import igraph
#' @export
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

# Returns the spectral densities for a list of adjacency matrices at the
# same points
nSpectralDensities <- function (graphs, from=NULL, to=NULL, bandwidth="Silverman") {
  npoints <- 1024
  ngraphs <- length(graphs)
  n <- ncol(graphs[[1]])
  spectra <- matrix(NA, n, ngraphs)
  for (i in 1:ngraphs) {
    A <- graphs[[i]]
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                      sqrt(nrow(A)))
    spectra[,i] <- eigenvalues
  }
  densities <- matrix(NA, npoints, ngraphs)
  minimum <- min(spectra)
  maximum <- max(spectra)
  if (!is.null(from) && !is.null(to)) {
    minimum <- from
    maximum <- to
  }
  for (i in 1:ngraphs) {
    f <- gaussianDensity(spectra[,i], bandwidth=bandwidth,
                         from=minimum, to=maximum,
                         npoints=npoints)

    densities[,i] <- f$y
    x <- f$x
  }
  return(list("x"=x, "densities"=densities))
}


# Given two adjacency matrices, returns the Jensen-Shannon divergence between
# the corresponding graphs
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

#' Given two spectral densities, returns the Jensen-Shannon divergence between
#' the corresponding graphs
#'
#' 'JS' computes the Jensen-Shannon divergence between the spectral density
#' functions of two graphs
#'
#' @param f1 an igraph graph object
#' @param f2 an igraph graph object
#' @return returns the Jensen-Shannon divergence between the corresponding graphs
#' @seealso \code{density}
#' @import igraph
#' @export
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
# Gene scores
# ------------------------------------------------------------------------------

#' Degree centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of vector containing the degree centrality of node of each network
#' @export
degreeCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, graph.strength)
  return(result)
}

#' Betweenness centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of vector containing the betweenness centrality of node of each network
#' @export
betweennessCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
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

#' Closeness centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of vector containing the closeness centrality of node of each network
#' @export
closenessCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
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

#' Eigenvector centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of vector containing the eigenvector centrality of node of each network
#' @export
eigenvectorCentrality <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  result <- lapply(G, function(x) evcent(x)$vector)
  return(result)
}

#' Clustering coefficient
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of vector containing the clustering coefficient of node of each network
#' @export
clusteringCoefficient <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
    v[a]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
  }
  weighted <- NULL
  if(any(v)) weighted <- TRUE
  if (weighted)
    result <- lapply(A, clusterCoef)
  else {
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    result <- lapply(G, transitivity, type="local", isolates="zero")
  }
  return(result)
}

# ------------------------------------------------------------------------------
# Network features
# ------------------------------------------------------------------------------

#' Average degree centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average degree centrality of each network.
#' @export
averageDegreeCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- degreeCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' Average betweenness centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average betweenness centrality of each network.
#' @export
averageBetweennessCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- betweennessCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' Average closeness centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average closeness centrality of each network.
#' @export
averageClosenessCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- closenessCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' Average eigenvector centrality
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average eigenvector centrality of each network.
#' @export
averageEigenvectorCentrality <- function(expr, labels, adjacencyMatrix) {
  result <- eigenvectorCentrality(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' Average clustering coefficient
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average clustering coefficient of each network.
#' @export
averageClusteringCoefficient <- function(expr, labels, adjacencyMatrix) {
  result <- clusteringCoefficient(expr, labels, adjacencyMatrix)
  return(lapply(result,mean))
}

#' Average shortest path
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @return a list of values containing the average shortest path of each network.
#' @export
averageShortestPath <- function(expr, labels, adjacencyMatrix) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
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

#' Spectral entropies
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return a list of values containing the spectral entropy of each network.
#' @export
spectralEntropies <- function(expr, labels, adjacencyMatrix, options=list(bandwidth="Sturges")) {
  A<-list()
  v<-vector(length=length(unique(labels)))
  for (a in 1:length(unique(labels))){
    A[[a]]<-adjacencyMatrix(expr[labels==unique(labels)[a],])
  }
  fs<- nSpectralDensities(A,bandwidth=options$bandwidth)
  if (sum(is.na(fs)) > 0)
    return(NA)
  entropies<-list()
  for(j in 1:length(A)){ # uma entropia para cada grafo
    entropies[[j]]<-entropy(list("x"=fs$x, "y"=fs$densities[,j]))
  }
  return(entropies)
}

# ------------------------------------------------------------------------------
# Test of equality between network properties
# ------------------------------------------------------------------------------

#' Degree centrality test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
degreeCentralityTest <- function(expr, labels, adjacencyMatrix,
                                 numPermutations=1000) {
    A<-list()
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    v<-vector(length=length(lab))
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    n <- ncol(expr)
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    s<-lapply(G, graph.strength)
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
    result<-mean(partial)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
        d=d+1
      }
      weighted <- NULL # Define o weighted como NULL, assim como na função original
      if(any(v)) weighted <- TRUE
      n <- ncol(expr)
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
      s<-lapply(G, graph.strength)
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
      res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
      results[i]<-mean(res)
    }
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Betweenness centrality test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
betweennessCentralityTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {
  # Betweenness centrality test for many graphs

    A<-list() # Lista onde serão colocadas as matrizes de adjacência
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1 # inicia o contador do indice da matriz dentro de "A"
    v<-vector(length=length(lab))
    for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
      A[[d]] <- adjacencyMatrix(expr[labels==k,]) # formação da matriz de adjacência
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    n <- ncol(expr) # numero de genes

    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
    s<-lapply(G, betweenness)# guarda o betweenness
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
    partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n)) # calcula a distancia euclidiana entre cada vetor de betweenness e o vetor medio
    result<-mean(partial) # A estatítica é a média das distancias
    results <- vector(length=numPermutations) # Vetor onde serão guardados todos os resultados das N permutações

    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
      A<-list() # DAQUI EM DIANTE OS PASSOS SÃO REPETIDOS PARA AS MATRZES PERMUTADAS
      d=1 # inicia o contador do indice da matriz dentro de "A"
      v<-vector(length=length(lab))
      for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
        A[[d]] <- adjacencyMatrix(expr[l==k,]) # formação da matriz de adjacência
        v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
        d=d+1
      }
      weighted <- NULL # Define o weighted como NULL, assim como na função original
      if(any(v)) weighted <- TRUE
      if (!is.null(weighted)) A<-lapply(A,invWeigthts)
      n <- ncol(expr) # numero de genes
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
      s<-lapply(G, betweenness)# guarda o betweenness
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
      res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n)) # calcula a distancia euclidiana entre cada vetor de betweenness e o vetor medio
      results[i]<-mean(res) # A estatítica é a média das distancias
    }

  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Closeness centrality test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
closenessCentralityTest <- function(expr, labels, adjacencyMatrix,
                                    numPermutations=1000) {
  # Closeness centrality test for many graphs
    A<-list() # Lista onde serão colocadas as matrizes de adjacência
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1 # inicia o contador do indice da matriz dentro de "A"
    v<-vector(length=length(lab))
    for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
      A[[d]] <- adjacencyMatrix(expr[labels==k,]) # formação da matriz de adjacência
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    n <- ncol(expr) # numero de genes
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
    s<-lapply(G, closeness)# guarda o closeness
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
    partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n)) # calcula a distancia euclidiana entre cada vetor de betweenness e o vetor medio
    result<-mean(partial) # A estatítica é a média das distancias
    results <- vector(length=numPermutations) # Vetor onde serão guardados todos os resultados das N permutações

    for (i in 1:numPermutations) {

      l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
      A<-list() # DAQUI EM DIANTE OS PASSOS SÃO REPETIDOS PARA AS MATRZES PERMUTADAS
      d=1 # inicia o contador do indice da matriz dentro de "A"
      v<-vector(length=length(lab))
      for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
        A[[d]] <- adjacencyMatrix(expr[l==k,]) # formação da matriz de adjacência
        v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
        d=d+1
      }
      weighted <- NULL # Define o weighted como NULL, assim como na função original
      if(any(v)) weighted <- TRUE
      if (!is.null(weighted)) A<-lapply(A,invWeigthts)
      n <- ncol(expr) # numero de genes
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
      s<-lapply(G, closeness)# guarda o closeness
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
      res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n)) # calcula a distancia euclidiana entre cada vetor de betweenness e o vetor medio
      results[i]<-mean(res) # A estatítica é a média das distancias
    }
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Eigenvector centrality test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
eigenvectorCentralityTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {
  # Eigenvector centrality test for many graphs
    A<-list()
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    v<-vector(length=length(lab))
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    # if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    n <- ncol(expr)
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    s<-lapply(G, function(x) evcent(x)$vector)
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
    result<-mean(partial)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        d=d+1
      }
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
      s<-lapply(G, function(x) evcent(x)$vector)
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
      res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
      results[i]<-mean(res) # A estatítica é a média das distancias
    }
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Clustering coefficient test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
clusteringCoefficientTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {
    A<-list()
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    v<-vector(length=length(lab))
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    n <- ncol(expr)
    if (!is.null(weighted)) s<-lapply(A, clusterCoef)
    else{
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
      s<-lapply(G, transitivity, type="local", isolates="zero")
    }
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    partial<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
    result<-mean(partial)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        d=d+1
      }
      G<-list()
      s<-matrix(NA,length(A)+1,n)
      if (!is.null(weighted)) s<-lapply(A, clusterCoef)
      else{
        G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
        s<-lapply(G, transitivity, type="local", isolates="zero")
      }
      s<-do.call(rbind,s)
      s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
      res<-apply(s[-dim(s)[1],],1, function(x) dist(rbind(x,s[dim(s)[1],]))/sqrt(n))
      results[i]<-mean(res)
    }
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Shortest path test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
shortestPathTest <- function(expr, labels, adjacencyMatrix,
                             numPermutations=1000) {
  # Shortest path test for many graphs
    A<-list()
    weighted <- NULL
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      if (sum(!(A[[d]] %in% c(1,0))) != 0) {
        weighted <- TRUE
        i1 <- which(A[[d]] > 0)
        A[[d]][i1] <- 1/A[[d]][i1]
      }
      d=d+1
    }
    G<-list()
    s<-vector(length=length(A))
    for(g in 1:length(A)){
      G[[g]] <- graph.adjacency(A[[g]], mode="undirected", weighted=weighted)
      if(is.null(weighted)) s[g]<-average.path.length(G[[g]], directed=FALSE)
      else {
        dist<-distances(G[[g]])
        s[g]<-mean(apply(dist, 1, function(x){ min(x[x!=0])}))
      }
    }
    s<-abs(s-mean(s))
    partial <- s
    result<-sum(s)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        if (sum(!(A[[d]] %in% c(1,0))) != 0) {
          weighted <- TRUE
          i1 <- which(A[[d]] > 0)
          A[[d]][i1] <- 1/A[[d]][i1]
        }
        d=d+1
      }
      G<-list()
      s<-vector(length=length(A))
      for(g in 1:length(A)){
        G[[g]] <- graph.adjacency(A[[g]], mode="undirected", weighted=weighted)
        if(is.null(weighted)) s[g]<-average.path.length(G[[g]], directed=FALSE)
        else {
          dist<-distances(G[[g]])
          s[g]<-mean(apply(dist, 1, function(x){ min(x[x!=0])}))
        }
      }
      s<-abs(s-mean(s))
      results[i]<-sum(s)
    }

  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

#' Degree distribution test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @export
degreeDistributionTest <- function(expr, labels, adjacencyMatrix,
                                   numPermutations=1000, options=list(bandwidth="Sturges")) {

    A<-list()
    weighted <- NULL
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      if (sum(!(A[[d]] %in% c(1,0))) != 0) weighted <- TRUE
      d=d+1
    }
    G<-list()
    for(g in 1:length(A)) G[[g]]<-graph.adjacency(A[[g]], mode="undirected", weighted=weighted)
    f<-nDegreeDensities(Gs=G, bandwidth=options$bandwidth)
    if(any(is.na(f))){
      cat('Empty graph')
      return(list("measure"=NA, "p.value"=NA))
    }
    meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
    partial <- vector(length=length(G))
    for (j in 1:length(G)) {
      f1 <- list("x"=f$x, "y"=f$densities[,j])
      partial[j] <- KL(f1, meanDensity)/length(G)
    }
    result <- sum(partial)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        d=d+1
      }
      G<-list()
      for(g in 1:length(A)) G[[g]]<-graph.adjacency(A[[g]], mode="undirected", weighted=weighted)
      f<-nDegreeDensities(Gs=G, bandwidth=options$bandwidth)
      if(any(is.na(f))){
        res<-NA
      }
      else{
        meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
        res<-0
        for (j in 1:length(G)) {
          f1 <- list("x"=f$x, "y"=f$densities[,j])
          res <- res + KL(f1, meanDensity)/length(G)
      }
      }
      results[i]<-res
    }

  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("measure"=result, "p.value"=pvalue,"Partial"=partial))
}

# ------------------------------------------------------------------------------

# Test for the absolute difference between the entropies of two graphs (genes
# networks), described by:
# H0: |H(A1) - H(A2)| = 0
# H1: |H(A1) - H(A2)| > 0
# where A1 and A2 are the adjacency matrices of the graphs
#
# -- Input: --
# "expr" - Matrix of genes (rows) vs microarrays (columns)
# "labels" - binary vector in which a position indicates the phenotype (0 or 1)
# of the corresponding microarray
# "adjacencyMatrix" - a function that returns the adjacency matrix for a given
# variables values matrix
# "numPermutations" - number of permutations that will be carried out in the
# permutation test
#
# -- Output: --
# A list containing:
# "abs.diff" - the absolute difference between the entropies of the two gene
# networks associated with each phenotype (0 or 1)
# "p.value" - the Nominal p-value of the test for the Jensen-Shannon divergence

#' Spectral entropy test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return A list containing:
# "abs.diff" - the absolute difference between the entropies of the two gene
# networks associated with each phenotype (0 or 1)
# "p.value" - the Nominal p-value of the test for the Jensen-Shannon divergence
#' @export
spectralEntropyTest <- function(expr, labels, adjacencyMatrix, numPermutations=1000,
                                options=list(bandwidth="Sturges")) {
    A<-list()
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1 # inicia contador da lista de matrizes "A"
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,]) #constrói a matrix de adjacencia selecionando determinado fator de labels
      d=d+1
    }
    f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
    entropies<-vector(length=length(A)) # vetor para guardar entropias
    for(j in 1:length(A)){ # uma entropia para cada grafo
      entropies[j]<-entropy(list("x"=f$x, "y"=f$densities[,j]))
    }
    meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities)) # Calcula a entropia média a partir de uma distribuicao media
    result<-sqrt((sum((entropies-entropy(meanDensity))^2))/length(entropies)) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
    partial<-sqrt(((entropies-entropy(meanDensity))^2)/length(entropies)) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
    results <- vector(length=numPermutations) # vetor para os resultados das permutacoes, para calculo do pvalor.
    for (i in 1:numPermutations) {
      d=1 # inicia contador da lista de matrizes "A"
      l <- sample(labels, replace = FALSE) # Reamostra os labels sem reposicao
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,]) # Monta as novas matrizes
        d=d+1
      }
      f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
      entropies<-vector(length=length(A)) # vetor para guardar entropias
      for(j in 1:length(A)){ # uma entropia para cada grafo
        entropies[j]<-entropy(list("x"=f$x, "y"=f$densities[,j]))
      }
      meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities)) # Calcula a entropia média a partir de uma distribuicao media
      results[i]<-sqrt((sum((entropies-entropy(meanDensity))^2))/length(entropies)) # idem e Calcula a raiz da soma dos quadrados das diferenças entre as entropias e a média
      }

  pvalue <- (1 + sum(results >= result))/(numPermutations + 1) # calculo do pvalor
  return(list("abs.diff"=result, "p.value"=pvalue,"Partial"=partial))
}

# Test for the Jensen-Shannon divergence between two graphs (genes networks),
# described by:
# H0: JS(A1, A2) = 0
# H1: JS(A1, A2) > 0
# where A1 and A2 are the adjacency matrices of the graphs
#
# -- Input: --
# "expr" - Matrix of genes (rows) vs microarrays (columns)
# "labels" - binary vector in which a position indicates the phenotype (0 or 1)
# of the corresponding microarray
# "adjacencyMatrix" - a function that returns the adjacency matrix for a given
# variables values matrix (the "expr" matrix)
# "numPermutations" - number of permutations that will be carried out in the
# permutation test
#
# -- Output: --
# A list containing:
# "JS" - The Jensen-Shannon divergence between the gene networks associated with
# each phenotype (0 or 1)
# "p.value" - the Nominal p-value of the test for the Jensen-Shannon divergence

#' Spectral distribution test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @param options a list containing parameters. It can be set to either
#' \code{list(bandwidth="Sturges")} or \code{list(bandwidth="Silverman")}.
#' @return A list containing:
#' "JS" - The Jensen-Shannon divergence between the gene networks associated with
#' each phenotype (0 or 1)
#' "p.value" - the Nominal p-value of the test for the Jensen-Shannon divergence
#' @export
spectralDistributionTest <- function(expr, labels, adjacencyMatrix,
                                     numPermutations=1000,
                                     options=list(bandwidth="Sturges")) {
    A<-list()
    lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
    if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[labels==k,])
      if (sum(is.na(A[[d]]))>0 | sum(is.infinite(as.matrix(A[[d]])))>0) stop(paste(sum(is.na(A[[d]]))+sum(is.infinite(as.matrix(A[[d]]))),"encontrados na matrix A[[",d ,"]] na medida da estatistica "))
      d=d+1
    }
    f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
    meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
    partial <- vector(length = length(A))
    for (j in 1:length(A)) {
      f1 <- list("x"=f$x, "y"=f$densities[,j])
      partial[j] <- KL(f1, meanDensity)/length(A)
    }
    result<-sum(partial)
    results <- vector(length=numPermutations)
    for (i in 1:numPermutations) {
      l <- sample(labels, replace = FALSE)
      A<-list()
      d=1
      for(k in lab) {
        A[[d]] <- adjacencyMatrix(expr[l==k,])
        if (sum(is.infinite(as.matrix(A[[d]])))>0) stop(print(expr[l==k,]))
        d=d+1
      }

      f<-nSpectralDensities(A, bandwidth=options$bandwidth) # Lista com as coordenadas (x,y) da dist. espectral dos grafos de "A"
      meanDensity <- list("x"=f$x, "y"=rowMeans(f$densities))
      res <- 0
      for (j in 1:length(A)) {
        f1 <- list("x"=f$x, "y"=f$densities[,j])
        res <- res + KL(f1, meanDensity)/length(A)
      }
      results[i]<-res
    }
  pvalue <- (1 + sum(results >= result))/(numPermutations + 1)
  return(list("KL"=result, "p.value"=pvalue,"Partial"=partial))
}

# ------------------------------------------------------------------------------
# Test of equality of vertices between network properties
# ------------------------------------------------------------------------------

#' Degree centrality test for Vector
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
degreeCentralityVertexTest <- function(expr, labels, adjacencyMatrix,
                                       numPermutations=1000) {
  A<-list()
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  d=1
  v<-vector(length=length(lab))
  for(k in lab) {
    A[[d]] <- adjacencyMatrix(expr[labels==k,])
    v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
    d=d+1
  }
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  if(any(v)) weighted <- TRUE
  n <- ncol(expr)
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  s<-lapply(G, graph.strength)
  s<-do.call(rbind,s)
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
  sp<-s
  result<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,mean)
  results <- list()
  for (i in 1:numPermutations) {
    l <- sample(labels, replace = FALSE)
    A<-list()
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[l==k,])
      d=d+1
    }
    n <- ncol(expr)
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    s<-do.call(rbind,lapply(G, graph.strength))
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    res<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
    results[[i]]<-apply(res,1,mean)
  }
  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= result),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(result,3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(t(sp[1:length(G),]),3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}

#' Betweenness centrality test
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
betweennessCentralityVertexTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {
  # Betweenness centrality test for many graphs

  A<-list() # Lista onde serão colocadas as matrizes de adjacência
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  d=1 # inicia o contador do indice da matriz dentro de "A"
  v<-vector(length=length(lab))
  for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
    A[[d]] <- adjacencyMatrix(expr[labels==k,]) # formação da matriz de adjacência
    v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
    d=d+1
  }
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  n <- ncol(expr) # numero de genes

  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
  s<-lapply(G, betweenness)# guarda o betweenness
  s<-do.call(rbind,s)
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
  sp<-s
  result<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,mean)
  results <- list() # Lista onde serão guardados todos os resultados das N permutações
  for (i in 1:numPermutations) {
    l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
    A<-list() # DAQUI EM DIANTE OS PASSOS SÃO REPETIDOS PARA AS MATRZES PERMUTADAS
    d=1 # inicia o contador do indice da matriz dentro de "A"
    v<-vector(length=length(lab))
    for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
      A[[d]] <- adjacencyMatrix(expr[l==k,]) # formação da matriz de adjacência
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    n <- ncol(expr) # numero de genes
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
    s<-lapply(G, betweenness)# guarda o betweenness
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
    res<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
    results[[i]]<-apply(res,1,mean)
  }

  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= result),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(result,3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(t(sp[1:length(G),]),3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}

#' Closeness centrality test for vectors
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
closenessCentralityVertexTest <- function(expr, labels, adjacencyMatrix,
                                    numPermutations=1000) {
  A<-list() # Lista onde serão colocadas as matrizes de adjacência
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  d=1 # inicia o contador do indice da matriz dentro de "A"
  v<-vector(length=length(lab))
  for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
    A[[d]] <- adjacencyMatrix(expr[labels==k,]) # formação da matriz de adjacência
    v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
    d=d+1
  }
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  if(any(v)) weighted <- TRUE
  if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  n <- ncol(expr) # numero de genes
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
  s<-lapply(G, closeness)# guarda o closeness
  s<-do.call(rbind,s)
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
  sp<-s
  result<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,mean)
  results <- list()
  for (i in 1:numPermutations) {

    l <- sample(labels, replace = FALSE) # Faz a permutação dos labels
    A<-list() # DAQUI EM DIANTE OS PASSOS SÃO REPETIDOS PARA AS MATRZES PERMUTADAS
    d=1 # inicia o contador do indice da matriz dentro de "A"
    v<-vector(length=length(lab))
    for(k in lab) { # Looping para criar uma matrix de adjacência para cada rede
      A[[d]] <- adjacencyMatrix(expr[l==k,]) # formação da matriz de adjacência
      v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
      d=d+1
    }
    weighted <- NULL # Define o weighted como NULL, assim como na função original
    if(any(v)) weighted <- TRUE
    if (!is.null(weighted)) A<-lapply(A,invWeigthts)
    n <- ncol(expr) # numero de genes
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)# guarda a rede
    s<-lapply(G, closeness)# guarda o closeness
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean)) # "s" é uma matrix onde serão guardados os vetores dos betweenness e um vetor de média deles
    res<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
    results[[i]]<-apply(res,1,mean)
  }
  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= result),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(result,3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(t(sp[1:length(G),]),3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}

#' Eigenvector centrality test for vectors
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
eigenvectorCentralityVertexTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {
  A<-list()
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  d=1
  v<-vector(length=length(lab))
  for(k in lab) {
    A[[d]] <- adjacencyMatrix(expr[labels==k,])
    v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
    d=d+1
  }
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  if(any(v)) weighted <- TRUE
  # if (!is.null(weighted)) A<-lapply(A,invWeigthts)
  n <- ncol(expr)
  G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
  s<-lapply(G, function(x) evcent(x)$vector)
  s<-do.call(rbind,s)
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
  sp<-s
  result<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,mean)
  results <- list()
  for (i in 1:numPermutations) {
    l <- sample(labels, replace = FALSE)
    A<-list()
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[l==k,])
      d=d+1
    }
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    s<-lapply(G, function(x) evcent(x)$vector)
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    res<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
    results[[i]]<-apply(res,1,mean)
  }
  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= result),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(result,3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(t(sp[1:length(G),]),3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}

#' Clustering coefficient test for vectors
#' @param expr Matrix of genes (rows) vs microarrays (columns)
#' @param labels a vector in which a position indicates the phenotype of the corresponding sample or state
#' @param adjacencyMatrix a function that returns the adjacency matrix for a given variables values matrix
#' @param numPermutations number of permutations that will be carried out in the permutation test
#' @export
clusteringCoefficientVertexTest <- function(expr, labels, adjacencyMatrix,
                                      numPermutations=1000) {

  A<-list()
  lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
  if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
  d=1
  v<-vector(length=length(lab))
  for(k in lab) {
    A[[d]] <- adjacencyMatrix(expr[labels==k,])
    v[d]<-(sum(!(A[[1]] %in% c(1,0))) != 0)
    d=d+1
  }
  weighted <- NULL # Define o weighted como NULL, assim como na função original
  if(any(v)) weighted <- TRUE
  n <- ncol(expr)
  if (!is.null(weighted)) s<-lapply(A, clusterCoef)
  else{
    G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
    s<-lapply(G, transitivity, type="local", isolates="zero")
  }
  s<-do.call(rbind,s)
  s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
  sp<-s
  result<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
  result<-apply(result,1,mean)
  results <- list()
  for (i in 1:numPermutations) {
    l <- sample(labels, replace = FALSE)
    A<-list()
    d=1
    for(k in lab) {
      A[[d]] <- adjacencyMatrix(expr[l==k,])
      d=d+1
    }
    if (!is.null(weighted)) s<-lapply(A, clusterCoef)
    else{
      G<-lapply(A,graph.adjacency, mode="undirected", weighted=weighted)
      s<-lapply(G, transitivity, type="local", isolates="zero")
    }
    s<-do.call(rbind,s)
    s<-rbind(s,apply(s,MARGIN=2,FUN=mean))
    res<-apply(s[1:(dim(s)[1]-1),],1,function(x) abs(x-s[dim(s)[1],]))
    results[[i]]<-apply(res,1,mean)
  }
  results<-do.call(rbind,results)
  pvalue <- (1 + apply(t(t(results) >= result),2,sum))/(numPermutations + 1)
  saida<-cbind(stat=round(result,3),pvalue=round(pvalue,4),qvalue=round(p.adjust(pvalue,method="fdr"),4),round(t(sp[1:length(G),]),3))
  rownames(saida)<-colnames(expr)
  saida<-saida[order(saida[,"pvalue"]),]
  colnames(saida)<-c("Test Statistic","Nominal p-value","Q-value",paste("Factor",0:max(as.numeric(lab))))
  return(saida)
}
