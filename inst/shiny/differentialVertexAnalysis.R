# Returns the analysis results matrix
vertexAnalysisTable <- reactive ({
  data<-plotSelectedData()
  vertexFunc <- input$vertexFunction
  classes <- data$classes
  results <- data.frame(matrix(NA, 1, ncol=4+length(classes)))#
  colnames(results) <- ifelse(is.null(classes),c("Node","Test statistic", "Nominal p-value", "q-value"),c("Node","Test statistic", "Nominal p-value", "q-value",paste(data$classes,"score",sep = " ")))#
  if (!values$canExecute || input$start==0)
    return(NULL)
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
    numPermutations <- values$numPermutations
    correlationMeasure <- input$correlationMeasure
    thrMeasure <- values$thrMeasure
    edgeWeight <- input$edgeWeight
    networkType <- values$networkType
    threshold <- input$thrValue
    options <- NULL
    seed <- values$seed
    printParameters <- function(){print(values$parameters)}
    if (is.null(expr) | is.null(labels) | is.null(numPermutations) | is.null(correlationMeasure))# is.null(geneSets) ||
      return(NULL)
    # if (thrMeasure=="correlation")
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
    thrEdge<-ifelse(thrMeasure=="none","none",
                    ifelse(thrMeasure=="correlation", "corr", 
                           ifelse(thrMeasure=="qvalue", "fdr", "pvalue")))
    print <- F
    funAdjMat <- adjacencyMatrix(method = correlationMeasures[correlationMeasure, 1],
                                       association = ifelse(edgeWeight=="correlation","corr", ifelse(edgeWeight=="qvalue","fdr", "pvalue")),
                                       threshold = thrEdge,
                                       thr.value = ifelse(thrEdge=="corr",threshold,1-threshold),
                                       weighted = ifelse(networkType=="weighted", T, F))
    # logFile=stdout()
    # saida<-list()
    method <- match.fun(differentialVertexAnalysis[vertexFunc, 2])

    results <- data.frame(matrix(NA, nrow=ncol(expr), ncol=5+length(classes)))#
    # names <- array(NA, length(geneSets))
    # for (i in 1:length(geneSets)) {
    #   names[i] <- geneSets[[i]][1]
    # }

    rownames(results) <- names(expr)
    colnames(results) <- c("N of Networks","Set size", "Test statistic", "Nominal p-value", "q-value",classes)#
    # temp <- tempfile(paste(paste(classes,collapse = ", "), "are being compared", sep=""),
                     # fileext=".txt")
    # withProgress(session, min=1, max=length(geneSets), {
      # setProgress(message = 'Analysis in progress',
                  # detail = 'This may take a while...')
      # for (i in 1:length(geneSets)) {
        # setName <- geneSets[[i]][1]
        # setProgress(value = i)
        # msg <- paste("Testing ", setName, " (", i, " of ", length(geneSets),
        #              " sets)", sep="")
        # setProgress(detail = msg)
        # if (print)
          # cat(msg, file=logFile, append=T)
        # genes <- geneSets[[i]][geneSets[[i]] %in% colnames(expr)]
        if (!is.null(seed))
          set.seed(seed)
        # result <- method(expr[,genes], labels, adjacencyMatrix=
        #                    adjacencyMatrix,  numPermutations=
        #                    numPermutations, options=options)
        results <- method(expr, labels, adjacencyMatrix=
                            funAdjMat,  numPermutations=
                           numPermutations, options=options)
        
        # if(!is.list(result)){
        #   saida[[setName]]<-result
        #   results<<-saida
        # }
        # if(is.list(result)){
        #   results[setName, "N of Networks"] <<- sum(unique(labels)!=-1)
        #   results[setName, "Test statistic"] <<- round(result[[1]],4)
        #   results[setName, "Nominal p-value"] <<- round(result$p.value,4)
        #   results[setName, "Set size"] <<- length(genes)
        #   results[setName, 6:ncol(results)] <<- round(result$Partial*100/sum(result$Partial),1)
        # }
      # }
      # if(is.list(result)) results[, "q-value"] <<- round(p.adjust(results[, "Nominal p-value"], method="fdr"),4)#
    # })
    results <- as.data.frame(results)
    results <- as.data.frame(cbind(rownames(results), results))
    colnames(results) <- c("Node","Test statistic", "Nominal p-value", "q-value",paste(data$classes,"score",sep = " "))#
    return(results)
  })
})

# Rendering --------------------------------------------------------------------

#_____Gene scores tab

# Render radio buttons to choose gene scores comparison file format.
output$vertexScoresType <- renderUI({
    if (is.null(vertexAnalysisTable())) {
        return(NULL)
    }
    radioButtons("vertexScoresType", paste("Select a file format to save the",
                 "gene scores table:"),
                  c("CSV", "R data"))
})

# Render button to download the gene scores comparison table
output$downloadVertexAnalysisButton <- renderUI({
    if (is.null(vertexAnalysisTable()))
        return(NULL)
    downloadButton("downloadVertexAnalysisTable", "Save gene scores")
})

# Prepare file with the statistics of the absolute differences between
# correlations for download
output$downloadVertexAnalysisTable <- downloadHandler(
    filename = function() {
      classes <- values$classes
      if (input$resultsType == "R data")
        name <- paste("Vertex_analysis", paste(classes,collapse = "_"),".RData", sep="")
      else
        name <- paste("Vertex_analysis", paste(classes,collapse = "_"),".csv" , sep="")
      return(name)
    },
    content = function(filename) {
        results <- vertexAnalysisTable()

        if (input$resultsType == "R data") {
          save(results, file=filename)
        }
        else {
          write.table(results, filename, append=T, row.names=F, col.names=T,
                      sep=";", dec=".",quote = F)
        }
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

output$vertexAnalysisTable  <- DT::renderDataTable({
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
