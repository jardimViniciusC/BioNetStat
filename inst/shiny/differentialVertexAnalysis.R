# Returns the analysis results matrix
vertexAnalysisTable <- reactive ({
  data<-plotSelectedData()
  classes <- data$classes
  results <- data.frame(matrix(NA, 1, ncol=4+length(classes)))#
  colnames(results) <- ifelse(is.null(classes),c("Node","Test statistic", "Nominal p-value", "q-value"),c("Node","Test statistic", "Nominal p-value", "q-value",paste(data$classes,"score",sep = " ")))#

  if (!values$canExecute || input$start==0)
    return(results)
  # createAlert(session, inputId = "resultsWarning",
  # message = paste("The analysis is running..."),
  # type = "info",
  # dismiss = TRUE,
  # block = FALSE,
  # append = FALSE)
  isolate({
    expr <- data$expr
    labels <- data$labels
    classes <- input$factorsinput
    geneSets <- NULL
    numPermutations <- values$numPermutations
    correlationMeasure <- values$correlationMeasure
    associationMeasure <- values$associationMeasure
    edgeWeight <- input$edgeWeight
    networkType <- values$networkType
    threshold <- input$correlationValue

    vertexFunc <- input$vertexFunction
    options <- NULL
    seed <- values$seed
    printParameters <- function(){print(values$parameters)}
    if (is.null(expr) || is.null(labels) ||# is.null(geneSets) ||
        is.null(numPermutations) || is.null(correlationMeasure))
      return(NULL)
    # if (associationMeasure=="correlation")
    #     col <- 1
    # else
    #     col <- 2
    # if (!is.null(options)) {
    #   name <- networkTestsMatrix[Function, "Options"]
    #   name <- strsplit(name, "=")[[1]][1]
    #   ops <- list()
    #   ops[tolower(name)] <- options
    #   options <- ops
    # }
    # networkInference <-
    #              match.fun(correlationMeasures[correlationMeasure, col])
    associationEdge<-ifelse(associationMeasure=="correlation",
                            "corr", ifelse(associationMeasure=="qvalue",
                                           "fdr", "pvalue"))
    print <- F
    adjacencyMatrix <- adjacencyMatrix(correlationMeasures[correlationMeasure, 1],
                                       associationEdge,
                                       ifelse(networkType=="weighted", ifelse(edgeWeight=="correlation",
                                                                              "corr", ifelse(edgeWeight=="qvalue",
                                                                                             "fdr", "pvalue")), associationEdge),
                                       threshold,
                                       ifelse(networkType=="weighted", T, F))
    logFile=stdout()
    saida<-list()
    method <- match.fun(differentialVertexAnalysis[vertexFunc, 2])
    if(is.null(geneSets)) geneSets <- list(c("all",colnames(expr)))
    results <- data.frame(matrix(NA, nrow=length(geneSets), ncol=5+length(classes)))#
    names <- array(NA, length(geneSets))
    for (i in 1:length(geneSets)) {
      names[i] <- geneSets[[i]][1]
    }
    rownames(results) <- names
    colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "q-value",classes)#
    temp <- tempfile(paste(paste(classes,collapse = ", "), "are being compared", sep=""),
                     fileext=".txt")
    withProgress(session, min=1, max=length(geneSets), {
      setProgress(message = 'Analysis in progress',
                  detail = 'This may take a while...')
      print(3)
      for (i in 1:length(geneSets)) {
        setName <- geneSets[[i]][1]
        setProgress(value = i)
        msg <- paste("Testing ", setName, " (", i, " of ", length(geneSets),
                     " sets)", sep="")
        setProgress(detail = msg)
        if (print)
          cat(msg, file=logFile, append=T)
        genes <- geneSets[[i]][geneSets[[i]] %in% colnames(expr)]
        if (!is.null(seed))
          set.seed(seed)
        result <- method(expr[,genes], labels, adjacencyMatrix=
                           adjacencyMatrix,  numPermutations=
                           numPermutations, options=options)
        if(!is.list(result)){
          saida[[setName]]<-result
          results<<-saida
        }
        if(is.list(result)){
          results[setName, "N of Networks"] <<- sum(unique(labels)!=-1)
          results[setName, "Test statistic"] <<- round(result[[1]],4)
          results[setName, "Nominal p-value"] <<- round(result$p.value,4)
          results[setName, "Set size"] <<- length(genes)
          results[setName, 6:ncol(results)] <<- round(result$Partial*100/sum(result$Partial),1)
        }
      }
      if(is.list(result)) results[, "q-value"] <<- round(p.adjust(results[, "Nominal p-value"], method="fdr"),4)#
    })
    print(4)
    results <- results[[1]]
    results <- cbind(rownames(results), results)
    colnames(results) <- c("Node","Test statistic", "Nominal p-value", "q-value",paste(data$classes,"score",sep = " "))#
    print(5)
    return(results)
  })
})

# Rendering --------------------------------------------------------------------

#_____Gene scores tab

# Render radio buttons to choose gene scores comparison file format.
output$vertexScoresType <- renderUI({
    if (is.null(geneScoresComparison())) {
        return(NULL)
    }
    radioButtons("vertexScoresType", paste("Select a file format to save the",
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
output$downloadVertexAnalysisTable <- downloadHandler(
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

output$vertexFunction <- renderUI({
  if (is.null(input$networkType))
    return(NULL)
  i <- which(differentialVertexAnalysis[,"Type"] == "both")
  if (input$networkType == "weighted")
    i <- union(i, which(differentialVertexAnalysis[,"Type"] == "weighted"))
  else if (input$networkType == "unweighted")
    i <- union(i, which(differentialVertexAnalysis[,"Type"] == "unweighted"))
  names <- rownames(differentialVertexAnalysis)[i]
    selectInput("vertexFunction", "Select a method to perform the nodes differential analysis:",
                        names)
})

output$vertexAnalysisTable  <- renderDataTable({
   return(vertexAnalysisTable())
})

# output$geneScoresComparison  <- renderChart2({
#     table <- geneScoresComparison()
#     return(dTable(table))
# })

# Select the application messages tab when the user clicks on the Start GNEA
# button
# observe({
#   input$startVertex
#   if (input$startVertex != 0)
#     updateTabsetPanel(sessionVertex, "FurtherTabSelected", selected =
#                         "Analysis results")
# })
