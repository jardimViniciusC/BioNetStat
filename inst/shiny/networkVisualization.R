# Network visualization properties ---------------------------------------------

genesOrder <- function(adjacencyMatrix) {
    ord <- order.dendrogram(as.dendrogram(hclust(as.dist(1 - adjacencyMatrix))))
}

# Returns correlation matrix colors
heatmapColors <- reactive({
    col <- input$heatmapColors
    if (col == "Green-Black-Red")
        return(colorRampPalette(c("green","black", "red"),space="rgb")(41))
    name <- switch(col,
                "Blue-White-Red"="RdBu",
                "Green-Yellow-Red"="RdYlGn",
                "Blue-Yellow-Red"="RdYlBu"
            )
    return(colorRampPalette(rev(brewer.pal(n=7, name=name)))(41))
})

# Verifies if it can plot the correlation heatmaps
canPlotHeatmaps <- reactive({
    if (is.null(filteredGeneSets()) ||
         is.null(labelsInput()))
        return(F)

    if (is.null(input$filterGeneSets))
        return(F)

    if (is.null(input$selectGeneSet))
        return(F)

    if (is.null(input$heatmapColors))
        return(F)

    expr <- exprInput()
    labels <- labelsInput()
    geneSets <- geneSetsInput()
    classes <- input$classes

    if (input$filterGeneSets %in% c("tested", "pvalueThreshold", "qvalueThreshold")) {
        expr <- values$expr
        labels <- values$labels
        geneSets <- values$filteredGeneSets
        classes <- values$classes
    }

    if (is.null(expr) || is.null(labels) || is.null(geneSets) ||
        is.null(classes))
        return(F)

    return(T)
})

# Returns a matrix of the absolute differences between the gene correlations
corAbsDiff <- reactive({
    data <- plotSelectedData()

    c1 <- "Class 1"
    c2 <- "Class 2"
    legend <- "Difference between associations"
    result <- data.frame(matrix(NA, nrow=1, ncol=5))

    colnames(result) <- c("Gene 1", "Gene 2", paste(c1, "association"),
                          paste(c2, "association"), legend)

    if (is.null(data))
        return(result)
    classes <- data$classes
    expr <- data$expr
    if (!canPlotHeatmaps())
        return(result)
    option <- input$heatmapDiffOptions
    if (is.null(option))
        return(result)
    c1 <- input$selectClassNetwork1
    c2 <- input$selectClassNetwork2
    genes <- colnames(expr)
    names <- combn(genes, 2)
    r <- adjacencyMatrices()
    n <- length(genes)
    r1 <- r[, 1:n]
    r2 <- r[, (n+1):(2*n)]
    diff <- r1 - r2
    legend <- paste("Difference between gene associations (", option, ")",
                    sep="")
    if (option ==  paste(c2, "-", c1))
        diff <- -diff
    else if (option == "abs") {
        diff <- abs(diff)
        legend <- "Absolute difference between gene associations"
    }

    result <- data.frame(matrix(NA, nrow=ncol(names), ncol=5))

    colnames(result) <- c("Gene 1", "Gene 2", paste(c1, "association"),
                          paste(c2, "association"), legend)
    for (i in 1:ncol(names)) {
        g1 <- names[1, i]
        g2 <- names[2, i]
        result[i, "Gene 1"] <- g1
        result[i, "Gene 2"] <- g2
        result[i, paste(c1, "association")] <- round(r1[g1, g2], 5)
        result[i, paste(c2, "association")] <- round(r2[g1, g2], 5)
        result[i, legend] <- round(diff[g1, g2], 5)
    }

    return(result)
})

# Returns a matrix of the absolute differences between the gene correlations
sitTable <- reactive({
  data <- plotSelectedData()

  c1 <- "Class 1"
  c2 <- "Class 2"
  legend <- "Difference between associations"
  result <- data.frame(matrix(NA, nrow=1, ncol=5))

  colnames(result) <- c("Gene 1", "Gene 2", paste(c1, "association"),
                        paste(c2, "association"), legend)

  if (is.null(data))
    return(result)
  classes <- data$classes
  expr <- data$expr
  if (!canPlotHeatmaps())
    return(result)
  option <- input$heatmapDiffOptions
  if (is.null(option))
    return(result)
  c1 <- input$selectClassNetwork1
  c2 <- input$selectClassNetwork2
  genes <- colnames(expr)
  names <- combn(genes, 2)
  r <- adjacencyMatrices()
  n <- length(genes)
  r1 <- r[, 1:n]
  r2 <- r[, (n+1):(2*n)]
  diff <- r1 - r2
  legend <- paste("Difference between gene associations (", option, ")",
                  sep="")
  if (option ==  paste(c2, "-", c1))
    diff <- -diff
  else if (option == "abs") {
    diff <- abs(diff)
    legend <- "Absolute difference between gene associations"
  }

  result <- data.frame(matrix(NA, nrow=ncol(names), ncol=3))
  result2 <- data.frame(matrix(NA, nrow=ncol(names), ncol=3))

  colnames(result) <- c("Gene 1", "Gene 2","association")
  colnames(result2) <- c("Gene 1", "Gene 2","association")
  for (i in 1:ncol(names)) {
    g1 <- names[1, i]
    g2 <- names[2, i]
    result[i, "Gene 1"] <- paste(g1,c1,sep = "_")
    result[i, "Gene 2"] <- paste(g2,c1,sep = "_")
    result[i, "association"] <- round(r1[g1, g2], 5)
  }
  for (i in 1:ncol(names)) {
    g1 <- names[1, i]
    g2 <- names[2, i]
    result2[i, "Gene 1"] <- paste(g1,c2,sep = "_")
    result2[i, "Gene 2"] <- paste(g2,c2,sep = "_")
    result2[i, "association"] <- round(r2[g1, g2], 5)
  }
  res<-rbind(result,result2)
  return(res)
})

# Rendering --------------------------------------------------------------------

# _____Network visualization plots tab

# Render select input for heatmap colors

# Render a select input of the classes that will be tested
output$factorsToNetViz1 <- renderUI({
  if (is.null(input$factorsinput))
    return(NULL)
  classes <- input$factorsinput
  # options <- vector()
  # for (i in 1:ncol(classes)) {
  #   options[i] <- paste(classes[1, i], classes[2, i])
  # }
  selectizeInput("selectClassNetwork1", p(h4(strong("Class 1\n")),h5("Choose the first condition to be visualized:")),
                           choices = classes, multiple=F)#c(options)
})
output$factorsToNetViz2 <- renderUI({
  if (is.null(input$factorsinput))
    return(NULL)
  classes <- input$factorsinput
  claN<-input$selectClassNetwork1
  # options <- vector()
  # for (i in 1:ncol(classes)) {
  #   options[i] <- paste(classes[1, i], classes[2, i])
  # }
  selectizeInput("selectClassNetwork2", p(h4(strong("Class 2\n")),h5("Choose the second condition to be visualized:")),
                 choices = classes[!(classes==claN)], multiple=F)#c(options)
})

output$heatmapColors  <- renderUI({
    selectInput("heatmapColors", "Select a color scheme:",
                 c( "Blue-White-Red","Green-Black-Red",
                    "Green-Yellow-Red", "Blue-Yellow-Red"))
})

output$networkPlotDimensions <- renderUI({
    if (is.null(input$networkPlotFormat))
        return(NULL)
    unit <- "pixels"
    default <- 480
    min <- 10
    max <- 10000
    if (input$networkPlotFormat == "PDF") {
        unit <- "inches"
        default <- 7
        min <- 1
        max <- 100
    }
    div(
        p(paste("Enter the plot dimensions (in ", unit, "):", sep="")),
        numericInput("networkPlotWidth", "Width:", default, min=min,
            max=max),
        numericInput("networkPlotHeight", "Height:", default, min=min,
            max=max)
    )
})

# Render button to download the class 1 network plot
output$downloadNetworkPlot1Button <- renderUI({
    if (!canPlotHeatmaps())
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    c1 <- input$selectClassNetwork1
    w <- input$networkPlotWidth
    h <- input$networkPlotHeight
    if (is.na(w) || is.na(h) || is.null(w) || is.null(h) || w < 1 || h < 1
        || !is.numeric(w) || !is.numeric(h))
        return(NULL)
    if (is.null(input$networkPlotFormat))
        return(NULL)
    downloadButton("downloadNetworkPlot1",
                   paste("Save", c1, "network plot"))
})

# Render button to download the class 2 network plot
output$downloadNetworkPlot2Button <- renderUI({
    if (!canPlotHeatmaps())
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    c2 <-input$selectClassNetwork2
    w <- input$networkPlotWidth
    h <- input$networkPlotHeight
    if (is.na(w) || is.na(h) || is.null(w) || is.null(h) || w < 1 || h < 1
        || !is.numeric(w) || !is.numeric(h))
        return(NULL)
    if (is.null(input$networkPlotFormat))
        return(NULL)
    downloadButton("downloadNetworkPlot2",
                   paste("Save", c2, "network plot"))
})

# Render button to download the plot of the differences between the networks
output$downloadNetworkDiffPlotButton <- renderUI({
    if (!canPlotHeatmaps())
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    c2 <- data$classes[[1]][2]
    w <- input$networkPlotWidth
    h <- input$networkPlotHeight
    if (is.na(w) || is.na(h) || is.null(w) || is.null(h) || w < 1 || h < 1
        || !is.numeric(w) || !is.numeric(h))
        return(NULL)
    if (is.null(input$networkPlotFormat))
        return(NULL)
    downloadButton("downloadNetworkDiffPlot", "Save plot")
})

# Prepare class 1 network plot for download
output$downloadNetworkPlot1 <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        classes <- data$classes
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        format <- input$networkPlotFormat
        if (format == "PNG")
            ext <- ".png"
        else if (format == "JPG")
            ext <- ".jpg"
        else
            ext <- ".pdf"
        name <- paste(input$selectGeneSet, "_", c1, "_network_",
                      input$networkType , "_", input$correlationMeasure,
                      "_", input$associationMeasure,
                      ifelse(input$networkType == "unweighted",
                             paste("_threshold=", input$plotEdgeThreshold, "_",
                                   sep=""), ""), ext, sep="")
    },
    content = function(filename) {
        data <- plotSelectedData()
        r <- adjacencyMatrices()
        breaks <- seq(-1, 1, 0.05)
        n <- nrow(r)
        ord <- genesOrder(r[, 1:n])
        col <- heatmapColors()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        format <- input$networkPlotFormat
        if (format == "PNG")
           saveFunction <- png
        else if (format == "JPG")
            saveFunction <- jpeg
        else
            saveFunction <- pdf
        saveFunction(filename, width=input$networkPlotWidth,
                     height=input$networkPlotHeight)
        pheatmap::pheatmap(r[ord, ord], col=col, cluster_rows=F, cluster_cols=F,
                 border_color=F, scale="none", breaks=breaks,
                 main=paste(c1, "network"))
        dev.off()
    }
)

# Prepare class 1 network plot for download
output$downloadNetworkPlot2 <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        format <- input$networkPlotFormat
        if (format == "PNG")
            ext <- ".png"
        else if (format == "JPG")
            ext <- ".jpg"
        else
            ext <- ".pdf"
        name <- paste(input$selectGeneSet, "_", c2, "_network_",
                      input$networkType , "_", input$correlationMeasure,
                      "_", input$associationMeasure,
                      ifelse(input$networkType == "unweighted",
                             paste("_threshold=", input$plotEdgeThreshold, "_",
                                   sep=""), ""), ext, sep="")
    },
    content = function(filename) {
        data <- plotSelectedData()
        r <- adjacencyMatrices()
        breaks <- seq(-1, 1, 0.05)
        n <- nrow(r)
        ord <- genesOrder(r[, 1:n])
        col <- heatmapColors()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        format <- input$networkPlotFormat
        if (format == "PNG")
           saveFunction <- png
        else if (format == "JPG")
            saveFunction <- jpeg
        else
            saveFunction <- pdf
        saveFunction(filename, width=input$networkPlotWidth,
                     height=input$networkPlotHeight)
        pheatmap::pheatmap(r[ord, n + ord], col=col, cluster_rows=F, cluster_cols=F,
                 border_color=F, scale="none", breaks=breaks,
                 main=paste(c2, "network"))
        dev.off()
    }
)


# Prepare the plot of the differences between the networks for download
output$downloadNetworkDiffPlot <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        if (input$heatmapDiffOptions == "abs")
            main <- paste("_abs_diff_between_", c1, "_and_", "_networks_",
                          sep="")
        else if (input$heatmapDiffOptions == paste(c2, "-", c1))
            main <- paste("_diff_between_", c2, "_and_", c1, "_networks_",
                         sep="")
        else
            main <- paste("_diff_between_", c1, "_and_", c2, "_networks_",
                          sep="")
        format <- input$networkPlotFormat
        if (format == "PNG")
            ext <- ".png"
        else if (format == "JPG")
            ext <- ".jpg"
        else
            ext <- ".pdf"
        name <- paste(input$selectGeneSet, main,
                      input$networkType , "_", input$correlationMeasure,
                      "_", input$associationMeasure,
                      ifelse(input$networkType == "unweighted",
                             paste("_threshold=", input$plotEdgeThreshold, "_",
                                   sep=""), ""), ext, sep="")
    },
    content = function(filename) {
        data <- plotSelectedData()
        r <- adjacencyMatrices()
        breaks <- seq(-1, 1, 0.05)
        n <- nrow(r)
        ord <- genesOrder(r[, 1:n])
        col <- heatmapColors()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2
        diff <- r[, 1:n] - r[, (n+1):(2*n)]
        format <- input$networkPlotFormat
        if (input$heatmapDiffOptions == "abs") {
            diff <- abs(diff)
            main <- "Absolute differences between association degrees"
        }
        else if (input$heatmapDiffOptions == paste(c2, "-", c1)) {
            diff <- -diff
            main <- paste("Differences between", c2, "and", c1,
                          "association degrees")
        }
        else {
            main <- paste("Differences between", c1, "and", c2,
                          "association degrees")
        }

        if (format == "PNG")
           saveFunction <- png
        else if (format == "JPG")
            saveFunction <- jpeg
        else
            saveFunction <- pdf
        saveFunction(filename, width=input$networkPlotWidth,
                     height=input$networkPlotHeight)
        #pheatmap(r[ord, ord], col=col, cluster_rows=F, cluster_cols=F,
        #         border_color=F, scale="none", breaks=breaks,
        #         main=paste(c1, "network"))
        #pheatmap(r[ord, n + ord], col=col, cluster_rows=F, cluster_cols=F,
        #         border_color=F, scale="none", breaks=breaks,
        #         main=paste(c2, "network"))
        pheatmap::pheatmap(diff[ord, ord], col=col, cluster_rows=F,
                 cluster_cols=F, border_color=F, scale="none",
                 main=main)
        dev.off()
    }
)

# Render gene correlation matrix for class 1
output$heatmapClass1 <- renderPlot({
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    classes <- input$selectClassNetwork1
    if (!canPlotHeatmaps())
        return(NULL)
    r <- adjacencyMatrices()
    if (is.null(r))
        return(NULL)
    breaks <- seq(-1, 1, 0.05)
    n <- nrow(r)
    col <- heatmapColors()
    ord <- genesOrder(r[, 1:n])
    pheatmap::pheatmap(r[ord, ord], col=col, cluster_rows=F, cluster_cols=F,
             border_color=F, scale="none", breaks=breaks,
             main=paste(classes, "network"))
})

# Render gene correlation matrix for class 2
output$heatmapClass2 <- renderPlot({
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    classes <- input$selectClassNetwork2
    if (!canPlotHeatmaps())
        return(NULL)
    r <- adjacencyMatrices()
    if (is.null(r))
        return(NULL)
    breaks <- seq(-1, 1, 0.05)
    n <- nrow(r)
    col <- heatmapColors()
    ord <- genesOrder(r[, 1:n])
    pheatmap::pheatmap(r[ord, n + ord], col=col, cluster_rows=F,
             cluster_cols=F, border_color=F, scale="none", breaks=breaks,
             main=paste(classes, "network"))
})

# Render matrix of differences options
output$heatmapDiffOptions <- renderUI({
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    classes <- list(c(input$selectClassNetwork1,input$selectClassNetwork2))
    if (!canPlotHeatmaps())
        return(NULL)
    c1 <- classes[[1]][1]
    c2 <- classes[[1]][2]
    options <- c(paste(c1, "-", c2),
                 paste(c2, "-", c1),
                 "Absolute differences between the association degrees"="abs")
    radioButtons("heatmapDiffOptions", h4(strong("Choose a matrix of differences:")),
                 options)
})

# Render difference correlation matrix
output$heatmapDiff <- renderPlot({
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    classes <- list(c(input$selectClassNetwork1,input$selectClassNetwork2))
    if (is.null(input$heatmapDiffOptions))
        return(NULL)
    r <- adjacencyMatrices()
    if (is.null(r))
        return(NULL)
    n <- nrow(r)
    col <- heatmapColors()

    c1 <- classes[[1]][1]
    c2 <- classes[[1]][2]

    diff <- r[, 1:n] - r[, (n+1):(2*n)]
    if (input$heatmapDiffOptions == "abs") {
        diff <- abs(diff)
        main <- "Absolute differences between association degrees"
    }
    else if (input$heatmapDiffOptions == paste(c2, "-", c1)) {
        diff <- -diff
        main <- paste("Differences between", c2, "and", c1,
                      "association degrees")
    }
    else {
        main <- paste("Differences between", c1, "and", c2,
                      "association degrees")
    }
    ord <- genesOrder(r[, 1:n])
    pheatmap::pheatmap(diff[ord, ord], col=col, cluster_rows=F,
             cluster_cols=F, border_color=F, scale="none",
             main=main)
})

# Render select inputs of two genes
output$selectGenes <- renderUI({
    if (!canPlotHeatmaps())
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    expr <- data$expr
    genes <- colnames(expr)
    if (is.null(genes))
        return(NULL)
    genes <- sort(genes)
    div(class="row-fluid",
        div(class="span4",
            selectInput("gene1", "Select a variable:", genes)),
        div(class="span4",
            selectInput("gene2", "Select another variable:", genes)))
})


# Render the correlation between two genes
output$corr <- renderDataTable({
    if (!canPlotHeatmaps())
        return(NULL)
    if (is.null(input$gene1) || is.null(input$gene2))
        return(NULL)
    data <- plotSelectedData()
    if (is.null(data))
        return(NULL)
    expr <- data$expr
    labels <- data$expr
    classes <- c(input$selectClassNetwork1,input$selectClassNetwork2)
    genes <- colnames(expr)
    i1 <- which(genes == input$gene1)
    i2 <- which(genes == input$gene2)
    if(length(i1) == 0 || length(i2) == 0)
        return(NULL)
    n <- length(genes)
    r <- adjacencyMatrices()
    if(is.null(r))
        return(NULL)
    #p1 <- cor.test(expr[input$gene1, labels==0], expr[input$gene2,
    #               labels==0], method=correlationMeasure)$p.value
    #p2 <- cor.test(expr[input$gene1, labels==1], expr[input$gene2,
    #               labels==1], method=correlationMeasure)$p.value
    #m <- matrix(NA, 2, 2)
    m <- matrix(NA, 1, 2)
    rownames(m) <- c("Association degree")
    colnames(m) <- c(classes[1], classes[2])
    m[1,1] <- r[i1, i2]
    m[1,2] <- r[i1, i2+n]
    #m[2,1] <- p1
    #m[2,2] <- p2
    return(m)
    #p(h5(classes[[1]][1], " correlation: "), r[i1, i2], br(),
     # h5(classes[[1]][2], " correlation: "), r[i1, i2+n])
})

# Render radio buttons that show the file with the statistics of the
# absolute differences between the gene correlations format options.
output$absDiffType <- renderUI({
    if (is.null(corAbsDiff())) {
        return(NULL)
    }
    radioButtons("absDiffType", h4(strong(paste("Select a file format to save the",
                 "list of gene association degrees:"))),
                  c("CSV", "R data"))
})

# Render button to download the statistics of the absolute differences
# between correlations
output$downloadAbsDiffButton <- renderUI({
    if (is.null(input$absDiffType))
        return(NULL)
    downloadButton("downloadAbsDiff", "Save list of association degrees")
})

# Prepare file with the statistics of the absolute differences between
# correlations for download
output$downloadAbsDiff <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        c1 <- input$selectClassNetwork1
        c2 <- input$selectClassNetwork2

        name <- paste(input$selectGeneSet, "_gene_association_degrees_", c1,
                      "_vs_", c2, "_",
                      input$networkType , "_", input$correlationMeasure,
                      "_", input$associationMeasure,
                      ifelse(input$networkType == "unweighted",
                             paste("_threshold=", input$plotEdgeThreshold, "_",
                             sep=""), ""), sep="")
        if (input$absDiffType == "R data")
            name <- paste(name, ".RData", sep="")
        else
            name <- paste(name, ".csv", sep="")
        return(name)
    },
    content = function(filename) {
        associationDegrees <- corAbsDiff()
        if (input$absDiffType == "R data")
            save(associationDegrees, file=filename)
        else
            write.csv(associationDegrees, filename, row.names=F)
    }
)

# Render table containing the average absolute difference of the gene
#correlations
output$corAbsDiff <- renderDataTable({
   corAbsDiff()
})

# output$corAbsDiff <- renderChart2({
#     table <- corAbsDiff()
#     return(dTable(table))
# })


# Render radio buttons that show the file with the statistics of the
# absolute differences between the gene correlations format options.
output$sitTableType <- renderUI({
  if (is.null(sitTable())) {
    return(NULL)
  }
  radioButtons("sitTableType", h4(strong(paste("Select a file format to save the",
                                               "list of gene association degrees:"))),
               c("CSV", "R data","TXT"))
})

# Render button to download the statistics of the absolute differences
# between correlations
output$downloadsitTableButton <- renderUI({
  if (is.null(input$sitTableType))
    return(NULL)
  downloadButton("downloadsitTable", "Save file to be opened in S.I.T. and to visualize and manipulate the network")
})

# Prepare file with the statistics of the absolute differences between
# correlations for download
output$downloadsitTable <- downloadHandler(
  filename = function() {
    data <- plotSelectedData()
    c1 <- input$selectClassNetwork1
    c2 <- input$selectClassNetwork2

    name <- paste(input$selectGeneSet, "_gene_association_degrees_", c1,
                  "_vs_", c2, "_",
                  input$networkType , "_", input$correlationMeasure,
                  "_", input$associationMeasure,
                  ifelse(input$networkType == "unweighted",
                         paste("_threshold=", input$plotEdgeThreshold, "_",
                               sep=""), ""), sep="")
    if (input$sitTableType == "R data") name <- paste(name, ".RData", sep="")
    if (input$sitTableType == "CSV") name <- paste(name, ".csv", sep="")
    if(input$sitTableType == "TXT") name <- paste(name, ".txt", sep="")
    return(name)
  },
  content = function(filename) {
    associationDegrees <- sitTable()

    if (input$sitTableType == "R data") save(associationDegrees, file=filename)
    if (input$sitTableType == "CSV")  write.csv(associationDegrees, filename, row.names=F)
    if (input$sitTableType == "TXT")  write.table(associationDegrees, filename, row.names=F,sep="\t",dec=",")
  }
)

output$sitTable <- renderDataTable({
  sitTable()
})

