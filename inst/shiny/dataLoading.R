# Reading file inputs ----------------------------------------------------------

# Read gene expression data from file
loadExprFile <- reactive({
  inFile <- input$expr
  if (is.null(inFile))
    return(NULL)
  #n <- nchar(inFile$name)
  #ext <- substr(inFile$name, n-3, n)D
    expr <- readVarFile(inFile$datapath,inFile$name,check.names=F)
  return(expr)
})

# Read gene expression data from file
exprInput <- reactive({
  expr <- loadExprFile()
  if (is.null(input$annotationOptions))
    return(expr)
  if (input$annotationOptions == "nocollapse")
    return(expr)
  if (is.null(input$collapsingMethod))
    return(expr)
  name <- collapsingMethodsMatrix[input$collapsingMethod, "Options"]
  annotation <- annotationInput()
  if (is.null(annotation))
    return(expr)
  if (name == "")
    options=NULL
  else {
    if (is.null(input$collapsingMethodOptions))
      return(expr)
    name <- strsplit(name, "=")
    name <- name[[1]][1]
    name <- strsplit(name, " ")
    name <- name[[1]]
    newName <- tolower(name[1])
    n <- length(name)
    if (n > 1) {
      for (i in 2:n) {
        str <- name[i]
        first <- toupper(substr(str, 1, 1))
        tail <- ""
        if (nchar(str) > 1)
          tail <- tolower(substr(str, 2, nchar(str)))
        newName <- paste(newName, first, tail, sep="")
      }
    }
    options <- list()
    options[newName] <- input$collapsingMethodOptions
  }
  method <- match.fun(collapsingMethodsMatrix[input$collapsingMethod,
                                              "Function"])
  return(collapseExprData(expr, annotation, method, options))
})

# Read file and return the names of categorical columns
chooseClass <- function(fileName) {
  table <- read.csv(fileName,header=T,dec=".",sep=";") # atentar para o decimal como virgula e separador ponto e virgula
  class <- names(which(!sapply(table,is.numeric)))
  return(class)
}

# Read categorical class file
labelsInput <- reactive({
  inFile <- input$expr
  if (is.null(inFile))
    return(NULL)
  labels <- chooseClass(inFile$datapath)
  return(labels)
})

# Read the data and returns the classes of choosed column
readClass <- function(fileName, factorName) {
  table <- read.csv(fileName,header=T,dec=".",sep=";") # atentar para o decimal como virgula
  class <- as.factor(table[,which(names(table)==factorName)])
  # labels <- combn(levels(class),2)
  labels <- levels(class)
  return(labels)
}
################################################
# Read categorical class file
classesInput <- reactive({
  inFile <- input$expr
  if (is.null(inFile) | is.null(input$classes))
    return(NULL)
  labels <- readClass(inFile$datapath, input$classes)
  return(labels)
})
############################################################
# Read gene set data file
geneSetsInput <- reactive ({
  inFile <- input$geneSets
  expr <- exprInput()
  if (is.null(inFile)) geneSets<-list(c("All",colnames(expr)))
  else geneSets <- readSetFile(inFile$datapath)
  return(geneSets)
})

# Phenotype properties ---------------------------------------------------------

# Vector of gene set sizes
geneSetSizes <- reactive({
  geneSets <- geneSetsInput()
  if (is.null(geneSets))
    return(NULL)
  sizes <- array(NA, length(geneSets))
  for (i in 1:length(geneSets)) {
    sizes[i] <- length(geneSets[[i]]) - 1
  }
  return(sizes)
})

# Vector of gene set sizes, tanking into account only the genes that are
# present in the gene expression data
geneSetSizesInChip <-  reactive({
  geneSets <- geneSetsInput()
  expr <- exprInput()
  if (is.null(geneSets) || is.null(expr))
    return(NULL)
  sizes <- array(NA, length(geneSets))
  for (i in 1:length(geneSets)) {
    genes <- which(geneSets[[i]] %in% colnames(expr))
    sizes[i] <- length(genes)
  }
  return(sizes)
})

# Min gene set size. If gene expression data is loaded, it takes into
# account only the genes that are present in the gene expression data.
# Otherwise, it considers all genes in the gene set file.
minSize <- reactive({
  if (is.null(exprInput())) {
    if (is.null(geneSetSizes()))
      return(dim(exprInput())[2])
    return(min(geneSetSizes()))
  }
  if (is.null(geneSetSizes()))
    return(dim(exprInput())[2])
  if(is.null(geneSetSizesInChip()))
    return(NULL)
  return(min(geneSetSizesInChip()))
})

# Max gene set size. If gene expression data is loaded, it takes into
# account only the genes that are present in the gene expression data.
# Otherwise, it considers all genes in the gene set file.
maxSize <- reactive({
  if (is.null(exprInput())) {
    if (is.null(geneSetSizes()))
      return(dim(exprInput())[2])
    return(max(geneSetSizes()))
  }
  if (is.null(geneSetSizes()))
    return(dim(exprInput())[2])
  if(is.null(geneSetSizesInChip()))
    return(NULL)
  return(max(geneSetSizesInChip()))
})

# Returns a string containg information about the phenotype classes
# classes <- reactive({
#   labels <- labelsInput()
#   symbols <- unique(labels[[3]])
#   classes <- paste(labels[[2]][2], " (", length(which(
#     labels[[3]] == symbols[1])), " samples)", sep="")
#   if (as.numeric(labels[[1]][2]) > 2)
#     for (i in 3:as.numeric(labels[[1]][2])) {
#       classes <- paste(classes, ", ", labels[[2]][i], " (",
#                        length(which(labels[[3]] == symbols[i-1])),
#                        " samples)", sep="")
#     }
#   i <- as.numeric(labels[[1]][2]) + 1
#   classes <- paste(classes, " and ", labels[[2]][i], " (",
#                    length(which(labels[[3]] == symbols[i-1])),
#                    " samples).", sep="")
#   return(classes)
# })

# Returns a vector containing 0 for each class1 sample, 1 for each class2
# sample and -1 for the remaining samples.
labels <- reactive({
  inFile<-input$expr
  classes <- input$classes
  if (is.null(classes))
    return(NULL)
  # classes <- strsplit(classes, " ")
  # labels <- labelsInput()
  l <- doLabels(fileName = inFile$datapath,factorName = classes,classes = input$factorsinput)
  return(l[,1])
})

# _____Data tab

# Render gene expression data information
output$expr <- renderUI({
  expr <- exprInput()
  if (is.null(input$expr) || is.null(expr))#
    return("No variable values data file was loaded.")
  else {
    dim <- dim(expr)
    filename <- paste("\"", input$expr$name, "\"", sep="")
    msg <- paste(filename, "was loaded successfully. The",
                 "matrix has", dim[2], "columns (variables) and", dim[1],
                 "rows (samples).")
    return(msg)
  }
})

# Render phenotype data information
output$labels <- renderUI({
  labels <- labelsInput()
  if (is.null(labels))
    return("No categorical class file was loaded.")
  else {
    filename <- paste("\"", input$classes, "\"", sep="")
    # factors <- paste("\"", input$factors, "\"", sep="")
    # classes <- classes()
    if (is.null(input$factorsinput)) msg <- paste("Factor",filename, "was choosed.", "Please, choose the classes to be compared")
      else{ if(length(input$factorsinput)==1) msg <- paste("Factor ",filename, " was choosed. ",
                                                        "The first class is ", input$factorsinput,
                                                        ". Please, choose unless one more class to be compared.",sep = "")
            else{ if(length(input$factorsinput)==2) msg <- paste("Factor",filename, "was choosed.",
                 "The factors are: ", paste(input$factorsinput[-c(length(input$factorsinput))],collapse = ", ")
                 , "and", input$factorsinput[c(length(input$factorsinput))])
              else msg <- paste("Factor ",filename, " was choosed. ",
                            "The factors are: ", paste(input$factorsinput[-c(length(input$factorsinput))],collapse = ", ")
                            , ", and ", input$factorsinput[c(length(input$factorsinput))], sep="")
            }
      }
    return(msg)
  }
})

# Render phenotype data information
# output$labels2 <- renderUI({
#   labels <- labels()
#   if (is.null(labels))
#     return("No categorical class file was loaded.")
#   else {
#     msg2 <- paste(levels(as.factor(labels)),collapse = ";")
#     print(levels(as.factor(labels)))
#     return(msg2)
#   }
# })

# Render a select input of the classes that will be tested
output$classes <- renderUI({
  if (is.null(labelsInput()))
    return(NULL)
  labels <- labelsInput()
  wellPanel(selectInput("classes", p(h4(strong("Column classes name\n")),h5("Choose the column of classes you want to test:",img(src="images/info.png", title="Select the column of your data table that you want to use to set the samples (rows) states."))),
              choices = c(labels)))
})

# Render a select input of the classes that will be tested
output$factors <- renderUI({
  if (is.null(classesInput()))
    return(NULL)
  classes <- classesInput()
  # options <- vector()
  # for (i in 1:ncol(classes)) {
  #   options[i] <- paste(classes[1, i], classes[2, i])
  # }
  wellPanel(selectizeInput("factorsinput", p(h4(strong("Classes\n")),h5("Choose the classes to be compared:",img(src="images/info.png", title="Select which states (samples) you want to bild a network to be compared. "))),
              choices = classes, multiple=T))#c(options)
})

# output$factor2 <- renderUI({
#   if (is.null(classesInput()))
#     return(NULL)
#   classes <- classesInput()
#   classes <- classes[-which(classes==input$factors1)]
#   # options <- vector()
#   # for (i in 1:ncol(classes)) {
#   #   options[i] <- paste(classes[1, i], classes[2, i])
#   # }
#   selectInput("factors2", "Choose the second factor:",
#               choices = classes)#c(options)
# })

# Render gene sets information
output$geneSets <- renderUI({
  geneSets <- geneSetsInput()
  expr <- exprInput()
  if (length(geneSets)==1 & geneSets[[1]][1]=="All")
    return(paste("No variables sets collection file was loaded. \n All", dim(expr)[2], "variables will be compared. \n (Depending on the number of variables, it may take a while) "))
  if (minSize()==0 & maxSize()==0)
    return(paste("Something went wrong. There is no variable name in Variable set database corresponding to collumn name of Variables values file"))
  else {
    filename <- paste("\"", input$geneSets$name, "\"", sep="")
    msg <- paste(filename, "was loaded successfully.",
                 "The collection contains", length(geneSets),
                 "variables sets ranging from")
    if (is.null(expr)) {
      msg <- paste(msg, min(geneSetSizes()), "to", max(geneSetSizes()),
                   "variables.")
    }
    else {
      msg <- paste(msg, " ", min(geneSetSizes()), " (", minSize(),
                   " in the variables values data) to ", max(geneSetSizes()),
                   " (", maxSize(), " in the variables values data) variables.",
                   sep="")
    }
    return(msg)
  }
})
output$table<-renderDataTable({
  if (is.null(exprInput())){
    expr <- exprInput()
  }
  else{
    expr <- round(exprInput(),3) 
  }
  expr
})
# Observing user inputs --------------------------------------------------------

# Select the data tab when a file is loaded
observe({
  input$expr
  input$geneSets
  input$labels
  input$classes
  input$factor
  updateTabsetPanel(session, "tabSelected", selected = "Loaded data")
  # updateTabsetPanel(sessionVertex, "tabSelected", selected = "Loaded data")
})
