
# Rendering --------------------------------------------------------------------
output$exprKeegMapColors  <- renderUI({
  data <- plotSelectedData()
  if (is.null(data)) {
    return(NULL)
  }
  selectInput("exprKeegMapColors", "Select a color scheme:",
              c("Green"="green3", "Blue"="blue","Yellow"="yellow", "Red"="red", "Gray"="gray"))
})
output$exprKeegMapColors2  <- renderUI({
  data <- plotSelectedData()
  if (is.null(data)) {
    return(NULL)
  }
  selectInput("exprKeegMapColors2", "Select a color scheme:",
              c("Green"="green3", "Blue"="blue","Yellow"="yellow", "Red"="red", "Gray"="gray"))
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
      results[,1]<-codes[l,2]
      centralityPathPlot(gene.data=results, cpd.data=NULL, threshold=input$thresholdSelected, thr.value=input$thrValue, species=input$speciesID , pathway.id=input$pathID, kegg.native=input$keggNative, file.name="file",
                       limit = NULL, bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                       mid =list(gene = "white", cpd = "white"),high = list(gene = input$exprKeegMapColors,cpd = input$exprKeegMapColors))
      tar(tarfile = file, files ="/tmp/pathMap")
      file.remove(list.files())
  }
)

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
    results[,1]<-codes[l,2]
    centralityPathPlot(gene.data=results, cpd.data=NULL, threshold=input$thresholdSelected, thr.value=input$thrValue, species=input$speciesID , pathway.id=input$pathID, kegg.native=input$keggNative, file.name="file",
                       limit = NULL, bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                       mid =list(gene = "white", cpd = "white"),high = list(gene = input$exprKeegMapColors,cpd = input$exprKeegMapColors))
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
