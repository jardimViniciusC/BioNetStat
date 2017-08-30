geneScoresComparison <- reactive({
    data <- plotSelectedData()
    if (is.null(data)) {
        c1 <- "Class 1"
        c2 <- "Class 2"
    }
    else {
        c1 <- data$classes[[1]][1]
        c2 <- data$classes[[1]][2]
    }
    result <- data.frame(matrix(NA, nrow=1, ncol=4))
    colnames(result) <- c("Gene symbol", paste(c1, "score"), paste(c2, "score"),
                           paste("Difference between", c1, "and", c2, "scores"))
    if (is.null(data))
        return(result)
    adjMatrix <- plotAdjacencyMatrix()
    geneScore <- input$geneScore
    signedCorrelation <- input$signedCorrelation
    if (is.null(data) || is.null(adjMatrix) || is.null(geneScore))
        return(result)
    
    adjacencyMatrix <- function(expr) {
        M <- adjMatrix(expr)
        return(abs(M))
    }

    r <- match.fun(geneScoresMatrix[geneScore, 2])(data$expr, data$labels, 
                                          adjacencyMatrix)

    genes <- rownames(data$expr)

    result <- data.frame(matrix(NA, nrow=length(genes), ncol=4))
    colnames(result) <- c("Gene symbol", paste(c1, "score"), paste(c2, "score"),
                           paste("Difference between", c1, "and", c2, "scores"))

    result[, 1] <- genes
    result[, 2] <- round(r[[1]], 6)
    result[, 3] <- round(r[[2]], 6)
    result[, 4] <- round(r[[1]] - r[[2]], 6)

    return(result)
})

# Rendering --------------------------------------------------------------------

#_____Gene scores tab

# Render radio buttons to choose gene scores comparison file format. 
output$geneScoresType <- renderUI({
    if (is.null(geneScoresComparison())) {
        return(NULL)
    }
    radioButtons("geneScoresType", paste("Select a file format to save the",
                 "gene scores table:"),
                  c("CSV", "R data"))
})

# Render button to download the gene scores comparison table
output$downloadGeneScoresButton <- renderUI({
    if (is.null(input$geneScoresType))
        return(NULL)
    downloadButton("downloadGeneScores", "Save gene scores")
})

# Prepare file with the statistics of the absolute differences between 
# correlations for download 
output$downloadGeneScores <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        classes <- data$classes
        c1 <- classes[[1]][1]
        c2 <- classes[[1]][2]
        geneScore <- input$geneScore
        geneScore <-  gsub(" ", "_", tolower(geneScore))
        name <- paste(input$selectGeneSet, "_gene_scores_", c1, "_vs_", c2, "_",
                      input$networkType , "_", input$plotCorrelationMeasure, 
                      "_", input$plotAssociationMeasure,
                      ifelse(input$networkType == "unweighted", 
                      paste("_threshold=", input$plotEdgeThreshold, sep=""),
                      ""),  sep="")
        if (input$geneScoresType == "R data")
            name <- paste(name, ".RData", sep="")
        else 
            name <- paste(name, ".csv", sep="")
        return(name)
    },
    content = function(filename) {
        geneScores <- geneScoresComparison()
        if (input$geneScoresType == "R data")
            save(geneScores, file=filename)
        else
            write.csv(geneScores, filename, row.names=F)
    }
)

output$geneScore <- renderUI({
    if (is.null(input$networkType))
        return(NULL)
    i <- which(geneScoresMatrix[,"Type"] == "both")
    if (input$networkType == "weighted")
        i <- union(i, which(geneScoresMatrix[,"Type"] == "weighted"))
    else if (input$networkType == "unweighted")
        i <- union(i, which(geneScoresMatrix[,"Type"] == "unweighted")) 
    names <- rownames(geneScoresMatrix)[i]
    selectInput("geneScore", "Select a method to compute the gene scores:",
                        names)
})

#output$geneScoresComparison  <- renderDataTable({
#    return(geneScoresComparison())
#})

output$geneScoresComparison  <- renderChart2({
    table <- geneScoresComparison()
    return(dTable(table))
})