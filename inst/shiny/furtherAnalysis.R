# Ruturns the selected gene set
searchGeneSet <- reactive({
    geneSets <- geneSetsInput()
    if (is.null(input$filterGeneSets))
        return(NULL)
    if (input$filterGeneSets %in% c("tested", "pvalueThreshold", "qvalueThreshold")) {
        geneSets <- values$filteredGeneSets
    }
    if (is.null(geneSets))
        return(NULL)
    geneSet <- input$selectGeneSet
    if (is.null(geneSet))
        return(NULL)
    for (i in 1:length(geneSets)) {
        if (geneSet == geneSets[[i]][1]) {
            return(geneSets[[i]][2:length(geneSets[[i]])])
        }
    }
    return(NULL)
})

# Network properties -----------------------------------------------------------

plotSelectedData <- reactive({
    expr <- exprInput()
    labels <- labels()
    geneSets <- geneSetsInput()
    classes <- input$classes
    geneSet <- input$selectGeneSet
    filterGeneSets <- input$filterGeneSets
    if (!is.null(classes))
        classes <- strsplit(classes, " ")
    if (is.null(filterGeneSets))
      return(NULL)
    if (filterGeneSets %in% c("tested", "pvalueThreshold", "qvalueThreshold")) {
        expr <- values$expr
        labels <- values$labels
        geneSets <- values$filteredGeneSets
        classes <- values$classes
    }
    if (is.null(expr) || is.null(labels) || is.null(classes) || is.null(geneSet))
        return(NULL)
    genes <- searchGeneSet()
    if (is.null(genes))
        return(NULL)
    i <- which(genes %in% colnames(expr))
    if (length(i) == 0)
        genes <- NULL
    else
        genes <- genes[i]
    if (is.null(genes))
        return(NULL)
    return(list("expr"=expr[,genes], "labels"=labels, "classes"=classes))
})

plotAdjacencyMatrix <- reactive({
    correlationMeasure <- input$correlationMeasure
    associationMeasure <- input$associationMeasure
    networkType <- input$networkType
    threshold <- input$correlationValue
    edgeWeight <- input$edgeWeight
    signedCorrelation <- input$signedCorrelation
    if(is.null(correlationMeasure) ||
       is.null(associationMeasure) || is.null(networkType))
      return(NULL)
    if (is.null(signedCorrelation))
      signedCorrelation <- F
    # if (associationMeasure == "correlation")
    #   col <- 1
    # else
    #   col <- 2
    associationEdge<-ifelse(associationMeasure=="correlation",
                            "corr", ifelse(associationMeasure=="qvalue",
                                           "fdr", "pvalue"))
    adjacencyMatrix <- adjacencyMatrix(correlationMeasures[correlationMeasure, 1],
                                       associationEdge,
                                       ifelse(networkType=="weighted", ifelse(edgeWeight=="correlation",
                                                                              "corr", ifelse(edgeWeight=="qvalue",
                                                                                             "fdr", "pvalue")), associationEdge),
                                       threshold,
                                       ifelse(networkType=="weighted", T, F),abs.values =!signedCorrelation )
})

# Returns a adjacency matrix whose first half columns belongs to the class1 gene
# network and second hallf columns belongs to the class2
# gene network
adjacencyMatrices <- reactive({
    data <- plotSelectedData()
    classes<-c(input$selectClassNetwork1,input$selectClassNetwork2)
    class <- input$factorsinput
    cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
    if (is.null(data))
        return(NULL)
    adjMatrix <- plotAdjacencyMatrix()
    if (is.null(adjMatrix))
        return(NULL)
    r1 <- adjMatrix(data$expr[data$labels==cla[cla[,1]==classes[1],2],])
    r2 <- adjMatrix(data$expr[data$labels==cla[cla[,1]==classes[2],2],])
    #diag(r1) <- diag(r2) <- 1
    genes <- colnames(data$expr)
    colnames(r1) <- colnames(r2) <- rownames(r1) <- rownames(r2) <- genes
    r <- cbind(r1, r2)
    colnames(r) <- c(genes, genes)
    return(r)
})

# Rendering --------------------------------------------------------------------

# _____Further analysis tab

# Render radio buttons to select a collection of gene sets
output$filterGeneSets <- renderUI({
    if (values$completed) {
        return(
            radioButtons("filterGeneSets", paste("Choose a collection of",
                   "gene sets:"),
                    c("All gene sets with p-values less than the threshold"=
                      "pvalueThreshold",
                      "All gene sets with q-values less than the threshold"=
                      "qvalueThreshold",
                      "All tested gene sets"="tested",
                      "All loaded gene sets"="all")
            )
        )
    }
    else if (!is.null(filteredGeneSets()) && !is.null(exprInput()) &&
             !is.null(labelsInput())) {
        return(
            radioButtons("filterGeneSets", paste("Choose a collection of",
                         "gene sets:"), c("All filtered gene sets"=
                                          "filtered",
                                          "All loaded gene sets"="all")
            )
        )
    }
    else return(NULL)
})

# Render a numeric input for the p-value threshold to filter the gene sets
# that will be displayed
output$geneSetThreshold <- renderUI({
    if (is.null(input$filterGeneSets))
        return(NULL)
    if (values$completed) {
        if (input$filterGeneSets %in% c("pvalueThreshold", "qvalueThreshold"))
            return(numericInput("geneSetThreshold", paste("Enter a",
                          "threshold:"), 0.05, min=0, max=1, step=0.05))
    }
    return(NULL)
})

# Render a select input to choose a genes set
output$selectGeneSet <- renderUI({
    if (is.null(filteredGeneSets()) || is.null(exprInput()) ||
         is.null(labelsInput()))
        return(NULL)

    if (is.null(input$filterGeneSets))
        return(NULL)

    n <- 0
    results <- results()

    if (input$filterGeneSets %in% c("pvalueThreshold", "qvalueThreshold")) {
        if (is.null(results))
            return(NULL)
        if (is.null(input$geneSetThreshold))
            return(NULL)
        if (input$filterGeneSets == "pvalueThreshold")
            i <- which(results[, "Nominal p-value"] <= input$geneSetThreshold)
        else
            i <- which(results[, "q-value"] <= input$geneSetThreshold)
        n <- length(i)
        if (n != 0)
            geneSets <- results[i, "Set name"]
    }

    else if (input$filterGeneSets == "tested") {
        geneSets <- values$filteredGeneSets
        n <- length(geneSets)
        names <- vector()
        if (n != 0)
            for (i in 1:n)
                names[i] <- geneSets[[i]][1]
        geneSets <- names
    }

    else if (input$filterGeneSets == "filtered") {
        i <- filteredGeneSets()
        geneSets <- geneSetsInput()
        n <- length(i)
        if (n != 0) {
            names <- vector()
            for (j in 1:n)
                names[j] <- geneSets[[i[j]]][1]
        }
        geneSets <- names
    }

    else if (input$filterGeneSets == "all") {
        geneSets <- geneSetsInput()
        n <- length(geneSets)
        names <- vector()
        if (n != 0)
            for (i in 1:n)
                names[i] <- geneSets[[i]][1]
        geneSets <- names
    }

    else
        return(NULL)

    if (n == 0)
        return("No filtered genes set.")

    geneSets <- sort(geneSets)

    selectInput("selectGeneSet", "Select a genes set:", choices=geneSets)
})

# Render selected genes set information
output$geneSetInfo <- renderUI ({
    if (!canPlotHeatmaps())
        return(NULL)
    if(is.null(input$filterGeneSets))
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    expr <- data$expr
    geneSet <- input$selectGeneSet
    if (is.null(geneSet))
        return(NULL)
    genes <- searchGeneSet()
    if (is.null(genes))
        return(NULL)
    n1 <- length(genes)
    i <- which(genes %in% colnames(expr))
    n2 <- length(i)
    msg <- paste("You have selected the ", geneSet, ". ", n2, " of the ",
                 n1,
                 " genes in this set were found in your expression data.")
})

output$selectedGeneSet <- renderUI({
    geneSet <- input$selectGeneSet
    if (is.null(geneSet))
        return(NULL)
    return(h4(paste(geneSet, "set analyses")))
})

output$networkScore <- renderUI({
    if (is.null(input$networkType))
        return(NULL)
    i <- which(networkScoresMatrix[,"Type"] == "both")
    if (input$networkType == "weighted")
        i <- union(i, which(networkScoresMatrix[,"Type"] == "weighted"))
    else if (input$networkType == "unweighted")
        i <- union(i, which(networkScoresMatrix[,"Type"] == "unweighted"))
    names <- rownames(networkScoresMatrix)[i]
    selectInput("networkScore", "Select a method to measure a network feature:",
                        names)
})

output$networkScoreOptions <- renderUI({
    if (is.null(input$networkTest))
        return(NULL)
    if (input$networkType == "")
        return(NULL)
    if(is.null(input$networkScore))
        return(NULL)
    options <- networkScoresMatrix[input$networkScore, "Options"]
    if (options == "")
        return(NULL)
    options <- strsplit(options, "=")
    name <- options[[1]][1]
    options <- strsplit(options[[1]][2], ",")
    options <- options[[1]]
    selectInput("networkScoreOptions", paste(name, ":", sep=""), options)
 })

output$networkScoresComparison  <- renderUI({
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)

    adjMatrix <- plotAdjacencyMatrix()
    networkScore <- input$networkScore
    signedCorrelation <- input$signedCorrelation
    if (is.null(data) || is.null(adjMatrix) || is.null(networkScore))
        return(NULL)

    options <- input$networkScoreOptions
    name <- networkScoresMatrix[networkScore, "Options"]
    if (is.null(options) && name != "")
        return(NULL)
    if (!is.null(options)) {
        name <- strsplit(name, "=")[[1]][1]
        ops <- list()
        ops[tolower(name)] <- options
        options <- ops
    }

    adjacencyMatrix <- function(expr) {
        M <- adjMatrix(expr)
        return(abs(M))
    }

    r <- match.fun(networkScoresMatrix[networkScore, 2])(data$expr, data$labels,
                                                adjacencyMatrix, options=options)
    r1 <- round(r[[1]], 6)
    r2 <- round(r[[2]], 6)
    diff <- round(r[[1]]-r[[2]], 6)
    classes <- list(c(input$selectClassNetwork1,input$selectClassNetwork2))
    result <- div(class="row-fluid",
                  div(class="span4",
                      p(h5(paste(classes[[1]][1], " score:", sep="")), r1)),
                  div(class="span4",
                      p(h5(paste(classes[[1]][2], " score:", sep="")), r2)),
                  div(class="span4",
                      p(h5(paste("Difference between ", classes[[1]][1], " and ",
                          classes[[1]][2], " scores:", sep="")), diff)))
    return(result)
})

# Gene scores ------------------------------------------------------------------

source("differentialVertexAnalysis.R", local=T)
source("keggPathway.R", local=T)

# Network visualization plots --------------------------------------------------

source("networkVisualization.R", local=T)

# Gene scores ------------------------------------------------------------------

source("geneScores.R", local=T)

# Gene expression analysis -----------------------------------------------------

source("geneExpression.R", local=T)
