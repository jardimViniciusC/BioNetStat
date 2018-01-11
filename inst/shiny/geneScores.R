geneScoresComparison <- reactive({
    data <- plotSelectedData()
    if (is.null(data)) {
      data<-list(classes=c("Class 1","Class 2"))
    }
    result <- data.frame(matrix(NA, nrow=1, ncol=1+length(data$classes)))
    colnames(result) <- c("Nodes", paste(data$classes,"score",sep = " "))
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
    result <- do.call(cbind,r)
    result<-cbind(rownames(result),result)
    colnames(result) <- c("Nodes", paste(data$classes,"score",sep = " "))
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
        classes <- list(c(input$selectClassNetwork1,input$selectClassNetwork2))
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

output$geneScoresComparison  <- renderDataTable({
   return(geneScoresComparison())
})

# output$geneScoresComparison  <- renderChart2({
#     table <- geneScoresComparison()
#     return(dTable(table))
# })
