# Gene set properties ----------------------------------------------------------

# Vector of the gene set ids (the position of the set in the list of gene
# sets) that were filtered.
filteredGeneSets <- reactive({
    if (is.null(exprInput()))
        geneSetSizes <- geneSetSizes()
    else
        geneSetSizes <- geneSetSizesInChip()
    if (is.null(geneSetSizes))
        return(1)
    a <- as.numeric(input$minSize)
    b <- as.numeric(input$maxSize)
    ge <- which(geneSetSizes >= a)
    le <- which(geneSetSizes <= b)
    return(intersect(ge, le))
})

# Number of gene sets that were filtered.
geneSetsCount <- reactive({
  filteredGenes<-filteredGeneSets()
  if(is.null(filteredGenes)) return(1)
  else return(length(filteredGenes))
})

# Rendering --------------------------------------------------------------------

# _____Slide bar

# Render filtered gene sets information:
# - Number of filtered gene sets
output$geneSetsCount <- renderUI({
    n <- geneSetsCount()
    if (n < 2)
        paste(n, "variable set was filtered.")
    else
        paste(n, "variable sets were filtered.")
})
# - Minimum gene set size
output$minSize <- renderUI({
    min <- as.numeric(ifelse(is.null(minSize()), 16, minSize()))
    max <- as.numeric(ifelse(is.null(maxSize()), 1024, maxSize()))
    numericInput("minSize", "Minimum variable set size", ifelse(min<10,min,10), min=min, max=max)
})
# - Maximum gene set size
output$maxSize <- renderUI({
    min <- as.numeric(ifelse(is.null(minSize()), 16, minSize()))
    max <- as.numeric(ifelse(is.null(maxSize()), 1024, maxSize()))
    numericInput("maxSize", "Maximum variable set size", max, min=min, max=max)
})


# Render a select input of the classes that will be tested
# output$classes <- renderUI({
#     if (is.null(labelsInput()))
#         return(NULL)
#     labels <- labelsInput()
#     names <- labels[[2]][2:length(labels[[2]])]
#     classes <- combn(names, 2)
#     options <- vector()
#     for (i in 1:ncol(classes)) {
#         options[i] <- paste(classes[1, i], classes[2, i])
#     }
#     selectInput("classes", "Choose the phenotype classes you want to test:",
#                 options)
# })

# Select predefined gene sets
output$predefGeneSets <- renderUI({
    selectInput("predefGeneSets", "",
                        names(predefGeneSets))
})

output$correlationMeasure <- renderUI({
    selectInput("correlationMeasure", h5("Select an association measure:"),
                rownames(correlationMeasures))
})

output$linkFormation <- renderUI({
  if(input$thrMeasure=="none")
    return()
  switch(input$thrMeasure, "correlation"= sliderInput("thrValue",h5("Minimum value (threshold) for link formation"),min=0, max=1, value = 0.7),
         "pvalue" = sliderInput("thrValue", h5("Minimum value (threshold) for link formation"),min = 0,max = 1,value = 0.95),
         "qvalue" = sliderInput("thrValue", h5("Minimum value (threshold) for link formation"),min = 0,max = 1,value = 0.95))
})

# Select predefined gene sets
output$networkTest <- renderUI({
    if (is.null(input$networkType))
        return(NULL)
    i <- which(networkTestsMatrix[,"Type"] == "both")
    if (input$networkType == "weighted")
        i <- union(i, which(networkTestsMatrix[,"Type"] == "weighted"))
    else if (input$networkType == "unweighted")
        i <- union(i, which(networkTestsMatrix[,"Type"] == "unweighted"))
    names <- rownames(networkTestsMatrix)[i]
    selectInput("networkTest", "Select a method to compare networks:",
                 names)
})

output$networkTestDescription <- renderUI({
    if (is.null(input$networkTest))
        return(NULL)
    description <- networkTestsMatrix[input$networkTest,"Description"]
    if (description == "")
        return(NULL)
    return(helpText(description))
})

output$networkTestOptions <- renderUI({
    if (is.null(input$networkTest))
        return(NULL)
    if (input$networkType == "")
        return(NULL)
    options <- networkTestsMatrix[input$networkTest, "Options"]
    if (options == "")
        return(NULL)
    options <- strsplit(options, "=")
    name <- options[[1]][1]
    options <- strsplit(options[[1]][2], "#")
    description <- options[[1]][2]
    options <- strsplit(options[[1]][1], ",")
    options <- options[[1]]
    name <- paste(name, ":", sep="")
    if (!is.na(description))
        name <- p(name, img(src="images/info.png", title=description))
    selectInput("networkTestOptions", name, options)
 })
