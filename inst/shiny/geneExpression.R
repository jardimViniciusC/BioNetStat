exprBoxplot <- function(expr, labels, title="",
                        ylab="Expression level") {
    data <- data.frame(labels=labels, expr=expr)
    ggplot2::qplot(labels, expr, data=data, geom=c("boxplot", "jitter"),
       color=labels, main=title, xlab="", ylab=ylab)
}

exprHeatmap <- function(expr, labels, col, clusterCols, clusterRows) {
    expr<-t(expr)
    annotation <- data.frame(Class=labels)
    rownames(annotation) <- colnames(expr)
    # c <- c("#FF0000FF", "#33FF00FF")
    c <- 1:length(unique(labels))
    names(c) <- unique(labels)
    colors <- list(Group=unique(labels))
    pheatmap::pheatmap(expr, annotation = annotation, col=col, annotation_cols=colors,
        cluster_cols=clusterCols, cluster_rows=clusterRows, scale="none", border_color=F)
}

# Returns correlation matrix colors
exprHeatmapColors <- reactive({
    col <- input$exprHeatmapColors
    if (col == "Green-Black-Red")
        return(colorRampPalette(c("green","black", "red"),space="rgb")(41))
    name <- switch(col,
                "Blue-White-Red"= "RdBu",
                "Green-Yellow-Red"="RdYlGn",
                "Blue-Yellow-Red"="RdYlBu"
            )
    return(colorRampPalette(rev(brewer.pal(n=7, name=name)))(41))
})

# Verifies if it can plot the expression heatmap
canPlotExprHeatmap <- reactive({
    if (is.null(input$exprHeatmapColors))
        return(F)
    data <- plotSelectedData()
    if (is.null(data))
        return(F)
    return(T)
})

diffExprTests <- reactive({
    data <- plotSelectedData()
    result <- data.frame(matrix(NA, nrow=1, ncol=5))
    if(length(input$factorsinput)==2){
      colnames(result) <- c("Variable symbol", "Difference between means",
                            "Difference between means test p-value",
                            "Difference in location",
                            "Wilcoxon-Mann-Whitney test p-value")
      if (is.null(data))
          return(result)
      expr <- data$expr
      labels <- data$labels
      classes <- input$factorsinput
      c1 <- classes[1]
      c2 <- classes[2]
      class <- input$factorsinput
      cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
      genes <- names(expr)
      result <- data.frame(matrix(NA, nrow=length(genes), ncol=5))
      colnames(result) <- c("Variable symbol", "Difference between means",
                            "Difference between means test p-value", "Difference in location",
                            "Wilcoxon-Mann-Whitney test p-value")
      for (i in 1:length(genes)) {
          result[i,"Variable symbol"] <- genes[i]
          expr1 <- expr[labels==cla[cla[,1]==c1,2], i]
          expr2 <- expr[labels==cla[cla[,1]==c2,2], i]
          t <-  t.test(expr1, expr2)
          result[i,"Difference between means"] <- round(t$estimate[1] -
                                                        t$estimate[2], 4)
          result[i,"Difference between means test p-value"] <- t$p.value
          w <- wilcox.test(expr1, expr2, conf.int=T)
          result[i,"Difference in location"] <- round(w$estimate, 4)
          result[i, "Wilcoxon-Mann-Whitney test p-value"] <- w$p.value
      }
    }
    else{
      colnames(result) <- c("Variable symbol", "F-value",
                            "Analisys of veriance test p-value",
                            "F-value",
                            "Kruskal-Wallis test p-value")
      if (is.null(data))
        return(result)
      expr <- data$expr
      labels <- data$labels
      class <- input$factorsinput
      cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
      genes <- names(expr)
      result <- data.frame(matrix(NA, nrow=length(genes), ncol=5))
      colnames(result) <- c("Variable symbol", "F-value",
                            "Analisys of veriance test p-value",
                            "F-value KS",
                            "Kruskal-Wallis test p-value")
      for (i in 1:length(genes)) {
        result[i,"Variable symbol"] <- genes[i]
        expr<-list()
        for(j in 1:length(class)){
          expr[[j]]<-data$expr[labels==cla[cla[,1]==class[j],2],]
        }
        expr <- do.call(rbind,expr)
        l<-list()
        for(j in 1:length(class)){
          l[[j]]<-rep(class[j], sum(labels==cla[cla[,1]==class[j],2]))
        }
        l<-do.call(c,l)
        t <-anova(lm(expr[,i]~as.factor(l)))
        result[i,"F-value"] <- round(t$`F value`[1], 4)
        result[i,"Analisys of veriance test p-value"] <- t$`Pr(>F)`[1]
        k<-kruskal.test(expr[,i],as.factor(l))
        result[i,"F-value KS"] <- round(k$statistic, 4)
        result[i, "Kruskal-Wallis test p-value"] <- k$p.value
      }
    }
    return(result)
})


# Rendering --------------------------------------------------------------------

#_____Gene expression tab

# Render select input for heatmap colors
output$exprHeatmapColors  <- renderUI({
    data <- plotSelectedData()
    if (is.null(data)) {
        return(NULL)
    }
    selectInput("exprHeatmapColors", "Select a color scheme:",
                 c("Green-Black-Red", "Blue-White-Red",
                    "Green-Yellow-Red", "Blue-Yellow-Red"))
})

output$exprHeatmapClustering  <- renderUI({
    data <- plotSelectedData()
    if (is.null(data)) {
        return(NULL)
    }
    checkboxGroupInput("exprHeatmapClustering", "",
                        c("Rows (genes)" = "row",
                          "Columns (samples)" = "col"))
})


output$exprHeatmap <- renderPlot({
    if (!canPlotExprHeatmap())
        return(NULL)
    col <- exprHeatmapColors()
    data <- plotSelectedData()
    clustering <- input$exprHeatmapClustering
    class <- input$factorsinput
    cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
    labels <- data$labels
    expr<-list()
    for(i in 1:length(class)){
      expr[[i]]<-data$expr[labels==cla[cla[,1]==class[i],2],]
    }
    expr <- do.call(rbind,expr)
    l<-list()
    for(i in 1:length(class)){
      l[[i]]<-rep(class[i], sum(labels==cla[cla[,1]==class[i],2]))
    }
    l<-do.call(c,l)
    clusterCols <- F
    clusterRows <- F
    if (!is.null(clustering)) {
        if ("col" %in% clustering)
            clusterCols <- T
        if ("row" %in% clustering)
            clusterRows <- T
    }
    exprHeatmap(expr, l, col, clusterCols, clusterRows)
})

# Render numeric inputs for setting the expression heatmap dimensions
output$exprHeatmapDimensions <- renderUI({
    format <- input$exprHeatmapFormat
    min <- 1
    max <- 100
    defaultWidth <- 755
    defaultHeight <- 480
    unit <- "pixels"
    if (is.null(input$exprHeatmapFormat))
        return(NULL)
    if (input$exprHeatmapFormat == "PDF") {
        min <- 10
        max <- 10000
        defaultWidth <- 11
        defaultHeight <- 7
        unit <- "inches"
    }
    div(p(paste("Enter the plot dimensions (in ", unit, "):", sep=""),
        numericInput("exprWidth", "Width:",
                     defaultWidth, min=min, max=max),
        numericInput("exprHeight",
                     "Height:", defaultHeight, min=min,
                     max=max)))
})

# Render button to download plots
output$downloadExprHeatmapButton <- renderUI({
    if (!canPlotExprHeatmap())
        return(NULL)
    w <- input$exprWidth
    h <- input$exprHeight
    if (is.na(w) || is.na(h) || is.null(w) || is.null(h) || w < 1 || h < 1
        || !is.numeric(w) || !is.numeric(h))
        return(NULL)
    if (is.null(input$exprHeatmapFormat))
        return(NULL)
    downloadButton("downloadExprHeatmap", "Save gene expression heatmap")
})

# Prepare expression heatmap for download
output$downloadExprHeatmap <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        classes <- input$factorsinput
        format <- input$exprHeatmapFormat
        if (format == "PNG")
            ext <- ".png"
        else if (format == "JPG")
            ext <- ".jpg"
        else
            ext <- ".pdf"
        name <- paste(input$selectGeneSet, "_expression_heatmap_", paste(classes,collapse = "_"), ext, sep="")
    },
    content = function(filename) {
        data <- plotSelectedData()
        class <- input$factorsinput
        cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
        col <- exprHeatmapColors()
        labels <- data$labels
        expr<-list()
        for(i in 1:length(class)){
          expr[[i]]<-data$expr[labels==cla[cla[,1]==class[i],2],]
        }
        expr <- do.call(rbind,expr)
        l<-list()
        for(i in 1:length(class)){
          l[[i]]<-rep(class[i], sum(labels==cla[cla[,1]==class[i],2]))
        }
        l<-do.call(c,l)
        format <- input$exprHeatmapFormat
        if (format == "PNG")
            saveFunction <- png
        else if (format == "JPG")
            saveFunction <- jpeg
        else
            saveFunction <- pdf
        saveFunction(filename, width=input$exprWidth, height=input$exprHeight)
        clusterCols <- F
        clusterRows <- F
        clustering <- input$exprHeatmapClustering
        if (!is.null(clustering)) {
            if ("col" %in% clustering)
                clusterCols <- T
            if ("row" %in% clustering)
                clusterRows <- T
        }
        exprHeatmap(expr, l, col, clusterCols, clusterRows)
        dev.off()
    }
)

# Render radio buttons to choose gene scores comparison file format.
output$diffExprTestsType <- renderUI({
    if (is.null(diffExprTests())) {
        return(NULL)
    }
    radioButtons("diffExprTestsType", paste("Select a file format to save the",
                 "results:"),
                  c("CSV", "R data"))
})

# Render button to download plots
output$downloadDiffExprTestsButton <- renderUI({
    if (is.null(diffExprTests()))
        return(NULL)
    downloadButton("downloadDiffExprTests", "Save results")
})

# Prepare plots for download
output$downloadDiffExprTests <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        classes <- input$factorsinput
        name <- paste(input$selectGeneSet, "_diff_expr_tests_", paste(classes,collapse = "_"),
                      sep="")
        if (input$diffExprTestsType == "CSV")
            name <- paste(name, ".csv", sep="")
        else if (input$diffExprTestsType == "R Data")
            name <- paste(name, ".RData", sep="")
    },
    content = function(filename) {
        diffExpressionAnalysis <- diffExprTests()
        if (input$diffExprTestsType == "R data")
            save(diffExpressionAnalysis, file=filename)
        else
            write.csv(diffExpressionAnalysis, filename, row.names=F)
    }
)

output$diffExprTests <- DT::renderDataTable({
   table <- diffExprTests()
if (!is.null(table))
# colnames(table)[1] <- "Gene symbol <img src=\"images/info.png\" title=\"teste\" />"
   return(table)
}#, options=list(aoColumns = list(list(bSearchable = FALSE),
                                # list(bSearchable = FALSE),
                                # list(bSearchable = FALSE),
                                # list(bSearchable = FALSE),
                                # list(bSearchable = FALSE)))
)

# output$diffExprTests <- renderChart2({
#     table <- diffExprTests()
#     data <- plotSelectedData()
#     name2 <- "Difference between means"
#     name4 <- "Difference in location"
#     c1 <- "population 1"
#     c2 <- "population 2"
#     if (!is.null(data)) {
#         classes <- list(c(input$selectClassNetwork1,input$selectClassNetwork2))
#         c1 <- classes[[1]][1]
#         c2 <- classes[[1]][2]
#         name2 <- paste(name2, " (", c1, " - ", c2, ")", sep="")
#         name4 <- paste(name4, " (", c1, " - ", c2, ")", sep="")
#     }
#
#     colnames(table)[2] <- name2
#     colnames(table)[3] <- paste("Difference between means test p-value <img ",
#                                 "src=\"images/info.png\" title=\"P-value of the ",
#                                 "t-test for the null hypothesis that ",
#                                 "the means of the two populations are equal.\"",
#                                 "/>", sep="")
#     colnames(table)[4] <- paste(name4, " <img src=\"images/info.png\" title=\"Median ",
#                                 "of the difference between a sample from ", c1,
#                                 " and a sample from ", c2, ".\"/>", sep="")
#     colnames(table)[5] <- paste("Difference between meadians test p-value <img",
#                                 " src=\"images/info.png\" title=",
#                                 "\"Wilcoxon-Mann-Whitney test p-value. The ",
#                                 "null hypothesis is that the medians are ",
#                                 "equal.\"/>", sep="")
#
#     return(dTable(table))
# })

# Render select inputs of a gene
output$selectGene <- renderUI({
    if (!canPlotExprHeatmap())
        return(NULL)
    expr <- exprInput()
    if (input$filterGeneSets %in% c("tested", "pvalueThreshold", "qvalueThreshold")) {
        expr <- values$expr
    }
    genes <- searchGeneSet()
    if (is.null(genes))
        return(NULL)
    i <- which(genes %in% names(expr))
    if (length(i) == 0)
        return(NULL)
    genes <- genes[i]
    genes <- sort(genes)
    selectInput("selectGene", "Select a gene:", genes)
})

output$exprBoxplot <- renderPlot({
    if (!canPlotExprHeatmap())
        return(NULL)
    gene <- input$selectGene
    if (is.null(gene))
        return(NULL)
    data <- plotSelectedData()
    if (!(gene %in% names(data$expr)))
        return(NULL)
    class <- input$factorsinput
    cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
    labels <- data$labels
    expr<-list()
    for(j in 1:length(class)){
      expr[[j]]<-data$expr[labels==cla[cla[,1]==class[j],2],]
    }
    expr <- do.call(rbind,expr)
    l<-list()
    for(j in 1:length(class)){
      l[[j]]<-rep(class[j], sum(labels==cla[cla[,1]==class[j],2]))
    }
    l<-do.call(c,l)
    plot(exprBoxplot(expr[,gene], l, gene))
})

output$exprBoxplotDimensions <- renderUI({
    if (is.null(input$exprBoxplotFormat))
        return(NULL)
    unit <- "pixels"
    default <- 480
    min <- 10
    max <- 10000
    if (input$exprBoxplotFormat == "PDF") {
        unit <- "inches"
        default <- 7
        min <- 1
        max <- 100
    }
    div(
        p(paste("Enter the plot dimensions (in ", unit, "):", sep="")),
        numericInput("boxplotWidth", "Width:", default, min=min,
            max=max),
        numericInput("boxplotHeight", "Height:", default, min=min,
            max=max)
    )
})

# Render button to download plots
output$downloadExprBoxplotButton <- renderUI({
    if (!canPlotExprHeatmap())
        return(NULL)
    if (is.null(input$selectGene))
        return(NULL)
    if (is.null(input$exprBoxplotFormat))
        return(NULL)
    w <- input$boxplotWidth
    h <- input$boxplotHeight
    if (is.na(w) || is.na(h) || is.null(w) || is.null(h) || w < 1 || h < 1
        || !is.numeric(w) || !is.numeric(h))
        return(NULL)
    downloadButton("downloadExprBoxplot", "Save boxplot")
})

# Prepare plots for download
output$downloadExprBoxplot <- downloadHandler(
    filename = function() {
        data <- plotSelectedData()
        classes <- input$factorsinput
        format <- input$exprBoxplotFormat
        if (format == "PNG")
            ext <- ".png"
        else if (format == "JPG")
            ext <- ".jpg"
        else
            ext <- ".pdf"
        name <- paste(input$selectGene, "_boxplot_", paste(classes,collapse = "_"), ext,
                      sep="")
    },
    content = function(filename) {
        format <- input$exprBoxplotFormat
        if (format == "PNG")
            saveFunction <- png
        else if (format == "JPG")
            saveFunction <- jpeg
        else
            saveFunction <- pdf
        data <- plotSelectedData()
        gene <- input$selectGene
        class <- input$factorsinput
        cla<-cbind(levels(as.factor(class)),c(0:(length(class)-1)))
        labels <- data$labels
        expr<-list()
        for(j in 1:length(class)){
          expr[[j]]<-data$expr[labels==cla[cla[,1]==class[j],2],]
        }
        expr <- do.call(rbind,expr)
        l<-list()
        for(j in 1:length(class)){
          l[[j]]<-rep(class[j], sum(labels==cla[cla[,1]==class[j],2]))
        }
        l<-do.call(c,l)
        saveFunction(filename, width=input$boxplotWidth,
                     height=input$boxplotHeight)
        plot(exprBoxplot(expr[,gene], l, gene))
        dev.off()
    }
)
