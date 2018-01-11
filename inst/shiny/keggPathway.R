# Returns correlation matrix colors
exprKeegMapColors <- reactive({
  col <- input$exprKeegMapColors
  if (col == "Green-Black-Red")
    return(colorRampPalette(c("green","black", "red"),space="rgb")(41))
  name <- switch(col,
                 "Blue-White-Red"= "RdBu",
                 "Green-Yellow-Red"="RdYlGn",
                 "Blue-Yellow-Red"="RdYlBu"
  )
  return(colorRampPalette(rev(brewer.pal(n=7, name=name)))(41))
})

# Rendering --------------------------------------------------------------------
output$exprKeegMapColors  <- renderUI({
  data <- plotSelectedData()
  if (is.null(data)) {
    return(NULL)
  }
  selectInput("exprKeegMapColors", "Select a color scheme:",
              c("Green-Black-Red", "Blue-White-Red",
                "Green-Yellow-Red", "Blue-Yellow-Red"))
})

output$thrValue <- renderUI({
  if(is.null(input$thresholdSelected))
    return(NULL)
  table<-vertexAnalysisTable()
  if(is.null(table))
    return(NULL)
  switch(input$thresholdSelected, "statistic"= sliderInput("thrValue",h5("Minimum value (threshold) for link construction"),min=min(as.numeric(table[,2])), max=max(as.numeric(table[,2])), value = min(as.numeric(table[,2]))),
         "pvalue" = sliderInput("thrValue", h5("Minimum value (threshold) for link construction"),min = 0,max = 1,value = 0.05),
         "qvalue" = sliderInput("thrValue", h5("Minimum value (threshold) for link construction"),min = 0,max = 1,value = 0.05))
})

# Prepare file with the statistics of the absolute differences between
# correlations for download
output$downloadKeggMap <- downloadHandler(
  filename = paste("data-", Sys.Date(), ".csv", sep="")
    ,
  content = function(file) {
    results <- vertexAnalysisTable()
    fileCodes<-input$keggCodes
    codes<-read.csv(fileCodes$datapath)
    if(input$selectingDataType=="gene")
      if(!dir.exists(file.path("/tmp/", "pathMap"))) dir.create("/tmp/pathMap")
      setwd("/tmp/pathMap")
      l<-which(results[,1] %in% codes[,1])
      print(results[,1])
      print(codes[,1])
      print(which(codes[,1] %in% results[,1]))
      print(codes[l,2])
      results[,1]<-codes[l,2]
      centralityPathPlot(gene.data=results, cpd.data=NULL, threshold=input$thresholdSelected, thr.value=input$thrValue, species=input$speciesID , pathway.id=input$pathID, kegg.native=input$keggNative, file.name="file",
                       limit = NULL, bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                       mid =list(gene = "white", cpd = "white"),high = list(gene = "red",cpd = "red"))
      tar(tarfile = file, files ="/tmp/pathMap")
      file.remove(list.files())
  }
)

output$exprKeegMapDimensions <- renderUI({
  format <- input$exprKeegMapFormat
  min <- 1
  max <- 100
  defaultWidth <- 755
  defaultHeight <- 480
  unit <- "pixels"
  if (is.null(input$exprKeegMapFormat))
    return(NULL)
  if (input$exprKeegMapFormat == "PDF") {
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

output$exprKeegMapClustering  <- renderUI({
  data <- plotSelectedData()
  if (is.null(data)) {
    return(NULL)
  }
  checkboxGroupInput("exprKeegMapClustering", "",
                     c("Rows (genes)" = "row",
                       "Columns (samples)" = "col"))
})
