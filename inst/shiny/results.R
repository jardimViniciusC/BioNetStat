# Results computation ------------------------------------------------------

# Returns the analysis results matrix
results <- reactive ({
    classes <- input$factorsinput
    results <- data.frame(matrix(NA, 1, ncol=6+length(classes)))#
    colnames(results) <- c("Variables set name","N of Networks","Set size", "Test statistic", "Nominal p-value", "q-value",classes)#
    if (!values$canExecute || input$start==0)
        return(results)
    createAlert(session, inputId = "resultsWarning",
                message = paste("The analysis is running..."),
                type = "info",
                dismiss = TRUE,
                block = FALSE,
                append = FALSE)
    isolate({
        expr <- values$expr
        labels <- values$labels
        classes <- input$factorsinput
        geneSets <- values$filteredGeneSets
        numPermutations <- values$numPermutations
        correlationMeasure <- values$correlationMeasure
        associationMeasure <- values$associationMeasure
        networkType <- values$networkType
        threshold <- input$correlationValue
        # threshold <- ifelse(associationMeasure=="correlation", 
        #                     input$correlationValue, ifelse(associationMeasure=="qvalue", 
        #                                                    input$qvalueThreshold, input$pvalueThreshold))
        networkTest <- values$networkTest
        options <- values$networkTestOptions
        seed <- values$seed
        printParameters <- function(){print(values$parameters)}
        if (is.null(expr) || is.null(labels) ||# is.null(geneSets) ||
            is.null(numPermutations) || is.null(correlationMeasure))
            return(NULL)
        # if (associationMeasure=="correlation")
        #     col <- 1
        # else 
        #     col <- 2
        if (!is.null(options)) {
            name <- networkTestsMatrix[networkTest, "Options"]
            name <- strsplit(name, "=")[[1]][1]
            ops <- list()
            ops[tolower(name)] <- options
            options <- ops
        }
        # networkInference <- 
        #              match.fun(correlationMeasures[correlationMeasure, col])
        print <- F
        adjacencyMatrix <- adjacencyMatrix(correlationMeasures[correlationMeasure, 1],
                                           ifelse(associationMeasure=="correlation", 
                                                  "corr", ifelse(associationMeasure=="qvalue", 
                                                                                "fdr", "pvalue")),
                                           ifelse(associationMeasure=="correlation", 
                                                  "corr", ifelse(associationMeasure=="qvalue", 
                                                                 "fdr", "pvalue")),
                                           threshold,
                                           ifelse(networkType=="weighted", T, F))
        logFile=stdout()
        method <- match.fun(networkTestsMatrix[networkTest, 2])
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
                  output[[setName]]<<-result
                  results<<-output
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
        results <- cbind(rownames(results), results)
        colnames(results)[1] <- "Set name"
        return(results) 
    })
})

# Rendering --------------------------------------------------------------------

# _____Results tab

# Render application messages 
output$appMessages <- renderUI({
    if (input$start == 0) {
        return("")
    }
    isolate({
        expr <- exprInput()
        labels <- labels()
        geneSets <-geneSetsInput()
        min <- input$minSize
        max <- input$maxSize
        numPermutations <- input$numPermutations
        correlationMeasure <- input$correlationMeasure
        associationMeasure <- input$associationMeasure
        networkTest <- input$networkTest
        networkType <- input$networkType
        threshold <- input$edgeThreshold
        values$canExecute <- F
        seed <- input$seed
        options <- input$networkTestOptions
        exprInputFileName <- input$expr$name
        geneSetsInputFileName <- input$geneSets$name
        labelsInputFileName <- input$labels$name
        if (is.null(expr)) {
            msg <- paste("To run the analysis, you must load the variables",
                         "data. Please, load the data on the side panel and",
                         "try again.")
            stop(msg)
        }
        if (is.null(labels)) {
            msg <- paste("To run the analysis, you must load the variables data with conditions information.",
                         "Please, load the variables data containing a column with conditions labels on the side panel and",
                         "try again.")
            stop(msg)
        }
        # if (is.null(geneSets)) {
        #     msg <- paste("To run the analysis, you must load the gene set data.",
        #                  "Please, load the data on the side panel and",
        #                  "try again.")
        #     stop(msg)
        # }
        if (is.na(min)) {
            msg <- "Please, enter a valid minimum set size."
            stop(msg)
        }
        if (is.na(max)) {
            msg <- "Please, enter a valid maximum set size."
            stop(msg)
        }
        if (min > max) {
            msg <- paste("Minimum set size should be less than or",
                         "equal to the maximum set size.")
            stop(msg)
        }
        if (min < minSize()) {
            min <- minSize()
            updateNumericInput(session, "minSize", min)
        }
        if (max > maxSize()) {
            max <- maxSize()
            updateNumericInput(session, "maxSize", max) 
        }
        if (length(filteredGeneSets) == 0) {
            msg <- paste("There are no gene sets with sizes between", 
                          min, "and", max, ". Please, choose a new size",
                          "interval on the side panel and try to run",
                          "again.")
            stop(msg)
        }
        if (minSize() < 2 && min < 2) {
            msg <- paste("All gene sets should contain at least 2 variables.", 
                         "Please, choose a new minimum size on the side", 
                         "panel and try again.")
            stop(msg)
        }
        if (numPermutations < 1) {
            msg <- paste("The number of permutations should be positive.")
            stop(msg)
        }
        geneSetDatabase <- paste("\"", geneSetsInputFileName, "\"", sep="")
        classes <- input$factorsinput
        
        # classes <- strsplit(classes, " ")
        l <- labels()
        names <- input$factorsinput
        
        phen <- paste(paste(names, " (", table(l)[-1],
                      " samples)", 
                      sep=""),collapse = ", ")
        associationMsg <- c("correlation"="absolute correlation",
                            "pvalue"="1 - pvalue", "qvalue"="1 - qvalue")
        if (is.null(seed))
          seedMsg <- "created from the current time and the process ID"
        else 
          seedMsg <- seed

        if (networkTestsMatrix[networkTest, "Options"] == "")
          options <- NULL
        if (is.null(options))
          optionsMsg <- ""
        else { 
          name <- networkTestsMatrix[networkTest, "Options"]
          # name <- strsplit(name, "=")[[1]][1]
          optionsMsg <- paste("*", name, "-", options, "\n")
        }
        thresholdMsg <- ""
        if (networkType == "unweighted")
          thresholdMsg <- paste("* Association degree threshold for edge", 
                                "selection -", threshold, "\n")

        parameters <- paste("* Variables values and conditions data from \"", exprInputFileName, "\"", 
                  " - ", ncol(expr)," variables and ", nrow(expr), " samples",
                  "\n", "* Conditions data - ", phen, "\n", 
                  "* Variable sets from ", geneSetDatabase, " - ", 
                  length(filteredGeneSets()), " filtered gene sets.", 
                  " The set sizes vary between ", min, " and ", 
                  max, "\n",
                  "* Network type - ", networkType, "\n",
                  "* Association measure - ", correlationMeasure, " (", 
                  associationMsg[associationMeasure], ")", "\n", 
                  thresholdMsg,
                  "* Method for networks comparison - ", networkTest, "\n",
                  optionsMsg, 
                  "* Number of label permutations - ", input$numPermutations,
                  "\n",  
                  "* Seed used to generate random permutations - ", 
                  seedMsg, "\n\n",  sep="")
        msg <- paste("The analysis is running with the following parameters:", 
                  "\n\n", parameters, sep="")
        cat(msg)
        if (optionsMsg != "")
          optionsMsg <- div(optionsMsg, br())
        if (thresholdMsg != "")
          thresholdMsg <- div(thresholdMsg, br())
        msg <- p(h5("Execution parameters"),
                  "* Variables values data from \"", exprInputFileName, "\"", 
                  " - ", ncol(expr)," variables and ", nrow(expr), " samples",
                  br(), "* Conditions data from \"", exprInputFileName, 
                  "\" - ", phen, br(), 
                  "* Variable sets from ", geneSetDatabase, " - ", 
                  length(filteredGeneSets()), " filtered variable sets.", 
                  " The set sizes vary between ", min, " and ", 
                  max, br(),
                  "* Network type - ", networkType, br(),
                  "* Association measure - ", correlationMeasure, " (", 
                  associationMsg[associationMeasure], ")", br(),
                  thresholdMsg,
                  "* Method for networks comparison - ", networkTest, br(),
                  optionsMsg,
                  "* Number of label permutations - ", input$numPermutations,
                  br(),
                  "* Seed used to generate random permutations - ", 
                  seedMsg)
        values$expr <- expr
        values$filteredGeneSets <- geneSets[filteredGeneSets()]
        values$labels <- l
        values$numPermutations <- numPermutations
        values$completed <- F
        values$canExecute <- T
        values$classes <- classes
        values$correlationMeasure <- correlationMeasure
        values$associationMeasure <- associationMeasure
        values$networkTest <- networkTest
        values$networkType <- networkType
        values$threshold <- threshold
        values$seed <- seed
        values$networkTestOptions <- options
        values$parameters <- parameters
        values$exprInputFileName <- exprInputFileName
        values$labelsInputFileName <- labelsInputFileName
        values$geneSetsInputFileName <- geneSetsInputFileName
        return(msg)
    })
})

# # Render table results
# output$results <- renderDataTable({
#    if (!values$completed)
#      return(NULL)
#    isolate({
#      return(results())
#    })}, options=list(aoColumns = list(list(bSearchable = FALSE),
#                                       list(bSearchable = FALSE),
#                                       list(bSearchable = FALSE),
#                                       list(bSearchable = FALSE),
#                                       list(bSearchable = FALSE))))



output$results <- renderDataTable({
      table <- results()
      #
      # colnames(table)[5] <- paste("Q-value <img src=\"images/info.png\" ",
      #                             "title=\"Adjusted p-value for ",
      #                             "multiple comparisons (Benjamin and ",
      #                             "Hochberg FDR method)\" />", sep="")
      #
      return(table)
})

# output$results <- renderChart2({
#   table <- results()
#   colnames(table)[5] <- paste("Q-value <img src=\"images/info.png\" ", 
#                               "title=\"Adjusted p-value for ",
#                               "multiple comparisons (Benjamin and ",
#                               "Hochberg FDR method)\" />", sep="")
#   return(dTable(table))
# })


# Render radio buttons that show the result file format options. 
output$resultsType <- renderUI({
    if (!values$completed) {
        return(NULL)
    }
    radioButtons("resultsType", "Select a file format to save the results:",
                  c("CSV", "R data"))

})

# Render download results button
output$downloadResultsButton <- renderUI({
    if (!values$completed)
        return(NULL)
    downloadButton("downloadResults", "Creat table of results")
})

# Prepare results file for download 
output$downloadResults <- downloadHandler(
    filename = function() {
        classes <- values$classes
        if (input$resultsType == "R data")
            name <- paste("BioNetStat_res_", paste(classes,collapse = "_"),".RData", sep="")
        else
            name <- paste("BioNetStat_res_", paste(classes,collapse = "_"),".csv" , sep="")
        return(name)
    },
    content = function(filename) {
        results <- results()
        rownames(results) <- NULL
        parameters <- matrix(NA, 18, ncol(results))
        expr <- values$expr
        classes <- input$factorsinput
        c1 <- classes[1]
        c2 <- classes[2]
        labels <- values$labels
        associationMsg <- c("correlation"="absolute correlation",
                            "pvalue"="1 - pvalue", "qvalue"="1 - qvalue")
        options <- values$networkTestOptions
        if (is.null(values$seed))
          seedMsg <- "created from the current time and the process ID"
        else 
          seedMsg <- values$seed

        if (networkTestsMatrix[values$networkTest, "Options"] == "")
          options <- NULL
        if (is.null(options))
          optionsMsg <- ""
        else { 
          name <- networkTestsMatrix[values$networkTest, "Options"]
          # name <- strsplit(name, "=")[[1]][1]
          optionsMsg <- paste(tolower(name), "-", options)
        }
        thresholdMsg <- "none (full graph)"
        if (values$networkType == "unweighted")
          thresholdMsg <- values$threshold
        parameters[1, 1] <- "BioNetStat differential network analysis"
        parameters[2, 1] <- "Date:"
        parameters[2, 2] <- date()
        parameters[3, 1] <- "Variables values input file:"
        parameters[3, 2] <- values$exprInputFileName
        parameters[4, 1] <- "Number of variables:" 
        parameters[4, 2] <- ncol(expr)
        parameters[5, 1] <- "Total number of samples:" 
        parameters[5, 2] <- nrow(expr)
        parameters[6, 1] <- "Sample labels input file:"
        # parameters[6, 2] <- values$labelsInputFileName
        parameters[7, 1] <- paste("Classes", 
                                  sep="")
        parameters[7, 2:(length(classes)+1)] <- classes
        parameters[8, 1] <- paste("Number of samples from each class:", 
                                  sep="")
        parameters[8, 2:(length(classes)+1)] <- table(labels[which(labels!=-1)])
        parameters[9, 1] <- "Variable sets collection:"
        parameters[9, 2] <- ifelse(is.null(values$geneSetsInputFileName),"All viriables",values$geneSetsInputFileName)
        parameters[10, 1] <- "Number of tested variable sets:"
        parameters[10, 2] <- ifelse(is.null(values$filteredGeneSets),1,length(values$filteredGeneSets))
        if(!is.null(values$filteredGeneSets)){
          parameters[11, 1] <- "Minimum variable set size:"
          parameters[11, 2] <-  min(as.numeric(lapply(values$filteredGeneSets, 
                                                      length)))
          parameters[12, 1] <- "Maximum variable set size:"
          parameters[12, 2] <-  max(as.numeric(lapply(values$filteredGeneSets, 
                                                      length)))
        }
        parameters[13, 1] <- "Network type:"
        parameters[13, 2] <- values$networkType
        parameters[14, 1] <- "Association measure:"
        parameters[14, 2] <- paste(values$correlationMeasure, " (", 
                                   associationMsg[values$associationMeasure], 
                                   ")", sep="")
        parameters[15, 1] <- "Association degree threshold for edge selection:"
        parameters[15, 2] <- thresholdMsg
        parameters[16, 1] <- "Method for networks comparison:"
        parameters[16, 2] <- values$networkTest
        if (!is.null(options))
          parameters[16, 2] <- paste(parameters[16, 2], " (", optionsMsg, ")",
                                     sep="")
        parameters[17, 1] <- "Number of label permutations:"
        parameters[17, 2] <- values$numPermutations
        parameters[18, 1] <- "Seed used to generate random permutations:"
        parameters[18, 2] <- seedMsg

        if (input$resultsType == "R data") {
          parameters <- parameters[,1:2]
          parameters <- as.table(parameters)
          rownames(parameters) <- NULL
          colnames(parameters) <- NULL
          save(results, parameters, file=filename)
        }
        else {
          parameters <- rbind(parameters, rep(NA, ncol(results)))
          parameters <- rbind(parameters, rep(NA, ncol(results)))
          parameters <- rbind(parameters, colnames(results))
          write.table(parameters, filename, na="", col.names=F, row.names=F, 
                      sep=",",quote = F)
          write.table(results, filename, append=T, row.names=F, col.names=F, 
                      sep=",", dec=".",quote = F)
        }
    }
)

# Render download results button
output$downloadNetworksButton <- renderUI({
  if (!values$completed)
    return(NULL)
  downloadButton("downloadNetworks", "Creat files to vizualize networks in S.I.T. software")
})

# Prepare results file for download 
output$downloadNetworks <- downloadHandler(
  filename = function() {
    classes <- values$classes
    name <- paste("SIT_networks_", paste(classes,collapse = "_"),".tar", sep="")
    return(name)
  },
  content = function(filename) {
    results <- results()
    rownames(results) <- NULL
    expr <- values$expr
    labels <- values$labels
    classes <- input$factorsinput
    geneSets <- values$filteredGeneSets
    numPermutations <- values$numPermutations
    correlationMeasure <- values$correlationMeasure
    associationMeasure <- values$associationMeasure
    networkType <- values$networkType
    threshold <- input$correlationValue
    networkTest <- values$networkTest
    options <- values$networkTestOptions
    seed <- values$seed
    adjacencyMatrix <- adjacencyMatrix(correlationMeasures[correlationMeasure, 1],
                                       ifelse(associationMeasure=="correlation", 
                                              "corr", ifelse(associationMeasure=="qvalue", 
                                                             "fdr", "pvalue")),
                                       ifelse(associationMeasure=="correlation", 
                                              "corr", ifelse(associationMeasure=="qvalue", 
                                                             "fdr", "pvalue")),
                                       threshold,
                                       ifelse(networkType=="weighted", T, F),abs.values = F)
    files<-vector(length = length(geneSets))
    if(!dir.exists(file.path("/tmp/", "BNS_to_SIT"))) dir.create("/tmp/BNS_to_SIT")
    setwd("/tmp/BNS_to_SIT")
    if(is.null(geneSets)) geneSets <- list(c("all",colnames(expr)))
    for (i in 1:length(geneSets)) {
      setName <- geneSets[[i]][1]
      genes <- geneSets[[i]][geneSets[[i]] %in% colnames(expr)]
      if (!is.null(seed))
        set.seed(seed)
      g<-make_empty_graph(n=0,directed = F)
      lab<-levels(as.factor(labels)) # salva os fatores de labels em lab.
      if(any(lab=="-1")) lab<-lab[-which(lab=="-1")] # se houver o fator "-1" ele é retirado dos fatores.
      for(l in lab){
        ma <- adjacencyMatrix(expr[labels==l,genes])
        g1<-graph.adjacency(ma,weighted=TRUE,mode="undirected")
        V(g1)$name<-paste(V(g1)$name,classes[as.numeric(l)+1],sep = "_")
        g<-g+g1
      }
      df <- get.data.frame(g)
      weight<-vector(length = dim(df)[1])
      for(w in 3:dim(df)[2]){
        weight[which(!is.na(df[,w]))]<-df[which(!is.na(df[,w])),w]
      }
      df<-cbind(df[1:2],weight)
      write.table(df,file = paste("SIT_", setName,".txt",sep=""), append=T, 
                  row.names=F, col.names=F, sep=",", dec=".")
    }
    tar(tarfile = filename, files ="/tmp/BNS_to_SIT")
    file.remove(list.files())
  }
)


# Render message when GNEA analysis completes
output$isCompleted <- renderUI({
    if (!values$completed || input$start==0)
      return(NULL)
    isolate({
        msg <- "The analysis completed successfully."
        if (!is.null(results)) {
          cat(paste("\n", msg, "\n", sep=""))
          updateTabsetPanel(session, "tabSelected", 
                            selected="Analysis results")
          return(h5("Results"))
        }
        else
          return(NULL)
      })
})

# Observing user inputs --------------------------------------------------------

# Verifies when it must run GNEA
observe({
    values$canExecute
    if (values$canExecute) {
        results()
        if (!is.null(results)) {
            values$completed <- T
        }
    }
})


# Select the application messages tab when the user clicks on the Start GNEA
# button
observe({
    input$start
    if (input$start != 0)
    updateTabsetPanel(session, "tabSelected", selected = 
                      "Analysis results")
})

# Alerts -----------------------------------------------------------------------

createAlert(session, inputId = "resultsWarning",
            message = paste("The analysis is not running. To start, load the",
                            "data, set the execution parameters, and then",
                            "click on the \"Start analysis\" button on the",
                            "sidebar."),
            type = "info",
            dismiss = TRUE,
            block = FALSE,
            append = FALSE)

observe({
    if (values$completed) {
        r <- results()
        if (sum(is.na(r)) > 0) {
            createAlert(session, inputId = "resultsWarning",
                message = paste("There are missing results in your analysis.",
                                "Please, check out the \"Help\" section on",
                                "\"Interpreting results\" to know",
                                "more about missing p-values."),
                type = "warning",
                dismiss = TRUE,
                block = FALSE,
                append = FALSE
            )
        }
        else {
            createAlert(session, inputId = "resultsWarning",
                  message = paste("The analysis completed successfully."),
                  type = "success",
                  dismiss = TRUE,
                  block = FALSE,
                  append = FALSE
              )
        }
    }
    #else {
    #    closeAlert(session, "resultsWarning")
    #}
})

# createAlert(session, anchorId = "resultsWarning",
#             message = paste("The analysis is not running. To start, load the",
#                             "data, set the execution parameters, and then",
#                             "click on the \"Start analysis\" button on the",
#                             "sidebar."),
#             type = "info",
#             dismiss = TRUE,
#             block = FALSE,
#             append = FALSE)
# 
# observe({
#   if (values$completed) {
#     r <- results()
#     if (sum(is.na(r)) > 0) {
#       createAlert(session, anchorId = "resultsWarning",
#                   message = paste("There are missing results in your analysis.",
#                                   "Please, check out the \"Help\" section on", 
#                                   "\"Interpreting results\" to know", 
#                                   "more about missing p-values."),
#                   type = "warning",
#                   dismiss = TRUE,
#                   block = FALSE,
#                   append = FALSE
#       )
#     }
#     else {
#       createAlert(session, anchorId = "resultsWarning",
#                   message = paste("The analysis completed successfully."),
#                   type = "success",
#                   dismiss = TRUE,
#                   block = FALSE,
#                   append = FALSE
#       )
#     }
#   }
#   #else {
#   #    closeAlert(session, "resultsWarning")
#   #}
# })