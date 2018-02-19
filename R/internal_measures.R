# ------------------------------------------------------------------------------
# Shannon entropy and Kullback Leibler divergence
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Given a partition x[1]...x[n] and y[i] = f(x[i]), returns the trapezoid sum
# approximation for int_{x[1]}^{x[n]}{f(x)dx}
trapezoidSum <- function (x, y) {
    n <- length(x)
    delta <- (x[2] - x[1])
    area <- sum(y[2:(n-1)])
    area <- (area + (y[1] + y[n])/2)*delta
    return(area)
}

# Returns the kernel bandwidth for a given sample x
kernelBandwidth <- function(x) {
    n <- length(x)
    # Sturges' criterion
    nbins <- ceiling(log2(n) + 1)
    return(abs(max(x) - min(x))/nbins)
}

# Returns the density function for a given sample x at n points in the interval
# [form, to]
gaussianDensity <- function(x, from=NULL, to=NULL, bandwidth="Sturges", npoints=1024) {
    if (bandwidth == "Sturges")
        bw<- kernelBandwidth(x)
    else if (bandwidth == "Silverman")
        bw <- bw.nrd0(x)
    if (bw == 0)
       return(NA)
    if (is.null(from) || is.null(to))
        f <- density(x, bw=bw, n=npoints)

    else
        f <- density(x, bw=bw, from=from, to=to, n=npoints)

    area <- trapezoidSum(f$x, f$y)
    return(list("x"=f$x, "y"=f$y/area))
}



# Returns the spectral entropy for a given density function f
entropy <- function(f) {
    y <- f$y
    n <- length(y)
    i <- which(y != 0)
    y[i] <- y[i]*log(y[i])
    return(-trapezoidSum(f$x, y))
}

# ------------------------------------------------------------------------------
# Kullback-Leibler divergence between graphs
# ------------------------------------------------------------------------------

# Given two spectral densities, returns the Kullback-Leibler divergence between
# the corresponding graphs
KL <- function(f1, f2) {
    y <- f1$y
    n <- length(y)
    for (i in 1:n) {
        if (y[i] != 0 && f2$y[i] == 0)
            return(Inf)
        if (y[i] != 0)
            y[i] <- y[i]*log(y[i]/f2$y[i])
    }
    return(trapezoidSum(f1$x, y))
}

# Given two adjacency matrices, returns the Kullback-Leibler divergence between
# the corresponding graphs
kullbackLeibler <- function(A1, A2, bandwidth="Sturges") {
    f <- spectralDensities(A1, A2, bandwidth=bandwidth)
    if (sum(is.na(f)) > 0)
        return(NA)
    return(KL(f$f1, f$f2))
}
