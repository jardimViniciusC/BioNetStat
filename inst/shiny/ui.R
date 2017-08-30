# ------------------------------------------------------------------------------
# CoGA User Interface (UI) definition VINICIUS
# ------------------------------------------------------------------------------
library(shiny)
# library(shinyBs)
# library(BioNetStat)
# Define UI for application 
shinyUI(fluidPage(
  
#   # Application title --------------------------------------------------------
  headerPanel(
    p(img(src="images/figuraBNS.png"), align="center"),
    windowTitle="CoGA"
  ),

  # Sidebar to load data and set application parameters ----------------------
  sidebarPanel(
    # Input files
    h3("1. Load data"),
    wellPanel(
      h5("Variables values data"),
      fileInput("expr", 
                paste("Please, select the text file (*.csv) containing",
                      "the variables values data"),
                accept=c("Comma-Seperated Values", "text/csv", ".csv"))
    ),
    h3("Factors"),
    
    uiOutput("classes"),
    
    uiOutput("factors"),
    # uiOutput("factor2"),
    
    wellPanel(
      h5("Variable set database"),
      fileInput("geneSets", 
                paste("Load the Set of variables file (*.gmt)", 
                      "describing the variables sets"),
                accept=c('.gmt'))
    ),
    # Parameters
    # h3("2. Set parameters"),
    # wellPanel(
    #   h5("Classes (conditions) being compared"),
    #   uiOutput("classes")
    # ),
    wellPanel(
      h5("Variable sets size range"),
      p("Enter the minimum and maximum variable set sizes."),
      uiOutput("minSize"),
      uiOutput("maxSize"),
      uiOutput("geneSetsCount")
    ),
    wellPanel(
      h5("Method for network inference"),
      uiOutput("correlationMeasure"),
      radioButtons(
        "associationMeasure", 
        "",
        c("Absolute correlation"="correlation", "1 - p-value"="pvalue",
          "1 - q-value"="qvalue")
      )
    ),
    wellPanel(
      h5("Minimum value threshold for link construction"),
      uiOutput("linkFormation")
    ),
    wellPanel(
      h5("Network type"),
      radioButtons(
        "networkType",
        "",
        c("Weighted"="weighted", "Unweighted"="unweighted")
      )#,
      # conditionalPanel(
      #   "input.networkType=='unweighted'",
      #   numericInput("edgeThreshold", 
      #                paste("Enter a association degree threshold to", 
      #                      "define the network edges:"), 
      #                0.95, min=0, max=1, step=0.01)
      # )
    ),
    wellPanel(
      h5("Method for networks comparison"),
      uiOutput("networkTest"),
      uiOutput("networkTestDescription"),
      uiOutput("networkTestOptions")
    ), 
    wellPanel(
      h5("Permutation test settings"),
      numericInput("numPermutations", 
                   "Enter the number of label permutations:", 
                   1000, min=0),
      numericInput("seed", 
                   "Enter a seed to generate random permutations:", 0)
    ),
    h3("3. Run differential network analysis"),
    wellPanel(
      actionButton("start", "Start analysis")
    )
    ),
  # Main panel ---------------------------------------------------------------
  mainPanel(
    h3("Reports"),
      
    # Tabset
    tabsetPanel(
      id="tabSelected",
      # Data tab
      tabPanel(
        h5("Loaded data"),
        h5("Variable values data"),
        uiOutput("expr"), 
        h5("Sample labels"),
        uiOutput("labels"),
        h5("Variable sets collection"), 
        uiOutput("geneSets"),
        value = "Loaded data",
        dataTableOutput("table")
      ),
      # Results tab
      tabPanel(
        h5("Analysis results"),
        bsAlert(inputId = "resultsWarning"),
        # bsAlert(anchorId =  "resultsWarning"),
        # progressInit(),
        uiOutput("appMessages"),
        uiOutput("isCompleted"),
        wellPanel(
          uiOutput("resultsType"),
          uiOutput("downloadResultsButton")
        ),
        dataTableOutput("results"),# 'datatables'
        # chartOutput('results', 'datatables'),
        wellPanel(
          uiOutput("downloadNetworksButton")
        ),
        value = "Analysis results"
      ),
      # Help tab
      tabPanel(
        h5("Help"), 
        bsCollapse(
          multiple=TRUE, 
          id="help",
          bsCollapsePanel(
            "1. Loading the dataset", 
            includeHTML("help/helpData.html"), 
            id="helpData"
          ),
          bsCollapsePanel(
            "2. Setting parameters and running the analysis", 
            includeHTML("help/helpParameters.html"), 
            id="helpParameters"
          ),
          bsCollapsePanel(
            "3. Interpreting the results", 
            includeHTML("help/helpResults.html"), 
            id="helpResults"
          ),
          bsCollapsePanel(
            "4. Further Analyses", 
            includeHTML("help/helpFurtherAnalyses.html"), 
            id="helpFurtherAnalyses"
          ),
          bsCollapsePanel(
            "5. Running BNS from the command-line interface", 
            includeHTML("help/helpCommandLine.html"), 
            id="helpCommandLine"
          ),
          bsCollapsePanel(
            "6. References", 
            includeHTML("help/helpReferences.html"), 
            id="helpReferences"
          )
        ),
        value = "Help"            
      )
      ))
))
  