% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_feature.R
\name{nDegreeDensities}
\alias{nDegreeDensities}
\title{Density functions of the degrees of n graphs}
\usage{
nDegreeDensities(
  Gs,
  npoints = 1024,
  bandwidth = "Sturges",
  from = NULL,
  to = NULL
)
}
\arguments{
\item{Gs}{a list of n igraph graphs objects}

\item{npoints}{number of points used in density function estimation}

\item{bandwidth}{a parameters. It can be set to either "Sturges" or "Silverman".}

\item{from}{the lower value used to build the distribution}

\item{to}{the higher value used to build the distribution}
}
\value{
a list containing the components 'x' and 'densities'. 
The first element is the vector 'x' of 'npoints' coordinates
of the points where the density function i estimated, and the second is
a vector 'y' of the estimated density values.
}
\description{
'nDegreeDensities' estimates the density functions of the degrees for n
graphs at the same coordinates
}
\examples{
G<-list()                                   
G[[1]]<-erdos.renyi.game(30,0.6)
G[[2]]<-barabasi.game(30,power = 1)
G[[3]]<-watts.strogatz.game(2,30,2,0.3)
d<-nDegreeDensities(G, npoints=1024, bandwidth="Sturges")
par(mfrow=c(1,3))
plot(d$x,d$densities[,1],main="Erdos-Renyi\n Degree distribution",
xlab="Degree",ylab="Frequency")
plot(d$x,d$densities[,2],main="Barabasi\n Degree distribution",
xlab="Degree",ylab="Frequency")
plot(d$x,d$densities[,3],main="Watts-Strogatz\n Degree distribution",
xlab="Degree",ylab="Frequency")
}
\seealso{
\code{graph.strength}

\code{density}
}
