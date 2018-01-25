
# Rendering --------------------------------------------------------------------
output$exprKeegMapColors  <- renderUI({
  data <- plotSelectedData()
  if (is.null(data)) {
    return(NULL)
  }
  selectInput("exprKeegMapColors", "Select a color scheme:",
              c("Green"="green3", "Blue"="blue","Yellow"="yellow", "Red"="red", "Black"="black"))
})
# output$exprKeegMapColors2  <- renderUI({
#   data <- plotSelectedData()
#   if (is.null(data)) {
#     return(NULL)
#   }
#   selectInput("exprKeegMapColors2", "Select a color scheme:",
#               c("Green"="green3", "Blue"="blue","Yellow"="yellow", "Red"="red", "Gray"="gray"))
# })

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

output$thrValue2 <- renderUI({
  if(is.null(input$thresholdSelected2))
    return(NULL)
  table<-vertexAnalysisTable()
  if(is.null(table))
    return(NULL)
  switch(input$thresholdSelected2, "statistic"= sliderInput("thrValue2",h5("Minimum value (threshold) for link construction"),min=min(as.numeric(table[,2])), max=max(as.numeric(table[,2])), value = min(as.numeric(table[,2]))),
         "pvalue" = sliderInput("thrValue2", h5("Minimum value (threshold) for link construction"),min = 0,max = 1,value = 0.05),
         "qvalue" = sliderInput("thrValue2", h5("Minimum value (threshold) for link construction"),min = 0,max = 1,value = 0.05))
})

# Prepare file with the statistics of the absolute differences between
# correlations for download
output$downloadKeggMap <- downloadHandler(
  filename = paste("data", Sys.Date(),input$speciesID, input$pathID, "centrality_bns.tar", sep="_")
    ,
  content = function(file) {
    results <- vertexAnalysisTable()
    fileCodes<-input$keggCodes
    codes<-read.csv(fileCodes$datapath)
    if(!dir.exists(file.path("/tmp/", "pathMap"))) dir.create("/tmp/pathMap")
    setwd("/tmp/pathMap")
    if(input$selectingDataType=="gene"){
        l<-which(results[,1] %in% codes[,1])
        results[,1]<-codes[l,2]
        centralityPathPlot(gene.data=results, cpd.data=NULL, threshold=input$thresholdSelected, thr.value=input$thrValue, species=input$speciesID , pathway.id=input$pathID, kegg.native=input$keggNative, file.name="file",
                         limit = NULL, bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                         mid =list(gene = "white", cpd = "white"),high = list(gene = input$exprKeegMapColors,cpd = input$exprKeegMapColors))
    }
    if(input$selectingDataType=="cpd"){
      l<-which(results[,1] %in% codes[,1])
      results[,1]<-codes[l,2]
      centralityPathPlot(gene.data=NULL, cpd.data=results, threshold=input$thresholdSelected, thr.value=input$thrValue, species=input$speciesID , pathway.id=input$pathID, kegg.native=input$keggNative, file.name="file",
                         limit = NULL, bins = list(gene = 15,cpd = 15), both.dirs= list(gene = F,cpd = F),
                         mid =list(gene = "white", cpd = "white"),high = list(gene = input$exprKeegMapColors,cpd = input$exprKeegMapColors))
    }
      tar(tarfile = file, files ="/tmp/pathMap")
      file.remove(list.files())
  }
)

output$downloadKeggExprMap <- downloadHandler(
  filename = paste("data", Sys.Date(),input$speciesID2, input$pathID2, "expression_bns.tar", sep="_")
  ,
  content = function(file) {
    data<-plotSelectedData()
    labels <- data$labels
    expr<-data$expr
    results <- vertexAnalysisTable()
    fileCodes<-input$keggCodes2
    codes<-read.csv(fileCodes$datapath)
    fun<-input$funSelected
    if(!dir.exists(file.path("/tmp/", "pathMap"))) dir.create("/tmp/pathMap")
      setwd("/tmp/pathMap")
    if(input$selectingDataType=="gene"){
      l<-match(names(expr),codes[,1])
      s<-match(results[,1],codes[,1])
      names(expr)<-codes[l,2]
      rownames(results)<-results[,1]<-codes[s,2]
      pathPlot(gene.data=t(expr), cpd.data=NULL, labels, varr.diff.list=results, threshold=input$thresholdSelected2, thr.value=input$thrValue2, FUN=fun,species=input$speciesID2 , pathway.id=input$pathID2, kegg.native=input$keggNative2, file.name="pathExpr")
    }
      if(input$selectingDataType=="cpd"){
        l<-which(results[,1] %in% codes[,1])
        results[,1]<-codes[l,2]
        pathPlot(gene.data=NULL, cpd.data=expr, labels, varr.diff.list=results, threshold=input$thresholdSelected2, thr.value=input$thrValue2, FUN=fun,species=input$speciesID2 , pathway.id=input$pathID2, kegg.native=input$keggNative2, file.name="pathExpr")
      }
    tar(tarfile = file, files ="/tmp/pathMap")
    file.remove(list.files())
  }
)

output$codeFile <- renderText({
  fileCodes<-input$keggCodes
  if(is.null(fileCodes)) return (NULL)
  codes<-read.csv(fileCodes$datapath)
  if(ncol(codes)==2) return(NULL)
  if(ncol(codes)!=2) return("The kegg code file has to has to have just two columns. The fisrt one with the variables names and the second one with the respective kegg codes")
})
output$codeFile2 <- renderText({
  fileCodes<-input$keggCodes2
  if(is.null(fileCodes)) return (NULL)
  codes<-read.csv(fileCodes$datapath)
  if(ncol(codes)==2) return(NULL)
  if(ncol(codes)!=2) return("The kegg code file has to has to have just two columns. The fisrt one with the variables names and the second one with the respective kegg codes")
})
