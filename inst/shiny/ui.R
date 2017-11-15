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
      fileInput("expr", p(paste("Please, select the text file (*.csv) containing the variables values data"),
                          img(src="images/info.png", title="Choose a file containing the table with the variables in the colums and the samples in the rows. The conditions of measure rows have to be indicated by a factor column")),
                accept=c("Comma-Seperated Values", "text/csv", ".csv"))
    ),
    h3("Factors"),

    uiOutput("classes"),

    uiOutput("factors"),
    # uiOutput("factor2"),

    wellPanel(
      h5("Variable set database"),
      fileInput("geneSets",
                p(paste("Load the Set of variables file (*.gmt)",
                      "describing the variables sets"),img(src="images/info.png", title="Choose a file containing the groups of variables. The format of the file that has to be inserted is explained in help section.")),
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
      ),
      conditionalPanel(
        "input.networkType=='weighted'",
        radioButtons(
          "edgeWeight",
          "Enter a association value to define the network edges weights:",
          c("Absolute correlation"="correlation", "1 - p-value"="pvalue",
            "1 - q-value"="qvalue")
        )
      )
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
      # Further analysis tab
      tabPanel(
        h5("Further analyses"),
        wellPanel(
          h5("Gene set selection"),
          div(
            class="row-fluid",
            div(class="span6", uiOutput("filterGeneSets"),
                uiOutput("geneSetThreshold")),
            div(class="span6", uiOutput("selectGeneSet"))
          ),
          uiOutput("geneSetInfo")
        ),

        br(),
        uiOutput("selectedGeneSet"),
        br(),

        tabsetPanel(
          id="FurtherTabSelected",
          tabPanel(
            "Nodes differential analysis",
            wellPanel(
              h5("Nodes Test"),
              uiOutput("vertexFunction"),
              h5("Run differential node analysis"),
              # actionButton("startVertex", "Start analysis"),
              # bsAlert(inputId = "vertexResultsWarning"),
              uiOutput("vertexScoresType"),
              uiOutput("downloadVertexAnalysisButton")
            ),
            br(),
            dataTableOutput("vertexAnalysisTable")
          # ),
          # tabPanel(
          #   "Network visualization plots"#,
            # bsCollapsePanel(
            #   "Plot settings",
            #   div(
            #     class="row-fluid",
            #     div(
            #       class="span4",
            #       h5("Colors selection"),
            #       uiOutput("heatmapColors")
            #     ),
            #     div(
            #       class="span4",
            #       h5("Plot format"),
            #       radioButtons(
            #         "networkPlotFormat",
            #         "Select a format to save the plots:",
            #         c("PDF","PNG", "JPG"))
            #     ),
            #     div(
            #       class="span4",
            #       h5("Plot dimensions"),
            #       uiOutput("networkPlotDimensions")
            #     )
            #   )
            # ),
            # bsCollapsePanel(
            #   "Network visualization",
            #   conditionalPanel(
            #     "input.associationMeasure=='correlation'",
            #     checkboxInput("signedCorrelation",
            #                   "Show negative correlations")
            #   ),
            #   br(),
            #   div(class="row-fluid",
            #       div(class="span6", plotOutput("heatmapClass1"),
            #           uiOutput("downloadNetworkPlot1Button")),
            #       div(class="span6", plotOutput("heatmapClass2"),
            #           uiOutput("downloadNetworkPlot2Button"))),
            #   br(),
            #   wellPanel(
            #     h5("Association between two gene products"),
            #     uiOutput("selectGenes"),
            #     tableOutput("corr")
            #   )
            # ),
            # bsCollapsePanel(
            #   "Differences between the gene networks",
            #   uiOutput("heatmapDiffOptions"),
            #   br(),
            #   bsCollapsePanel(
            #     "Matrix of differences",
            #     plotOutput("heatmapDiff", width="50%"),
            #     uiOutput("downloadNetworkDiffPlotButton")),
            #   bsCollapsePanel(
            #     "List of gene association degrees",
            #     wellPanel(
            #       uiOutput("absDiffType"),
            #       uiOutput("downloadAbsDiffButton")
            #     ),
            #     br(),
            #     chartOutput("corAbsDiff", "datatables")
            #   )
            # )
          # ),
          # tabPanel(
          #   "Gene set properties"#,
          #   # wellPanel(
          #   #   h5("Gene set network topological properties"),
          #   #   uiOutput("networkScore"),
          #   #   uiOutput("networkScoreOptions"),
          #   #   uiOutput("networkScoresComparison")
          #   # )
          # ),
          # tabPanel(
          #   "Gene scores",
          #   wellPanel(
          #     h5("Gene scores"),
          #     uiOutput("geneScore"),
          #     uiOutput("geneScoresType"),
          #     uiOutput("downloadGeneScoresButton")
          #   ),
          #   br(),
          #   dataTableOutput("geneScoresComparison")
          #   # chartOutput("geneScoresComparison", "datatables")
          # ),
          # tabPanel(
          #   "Gene expression analysis"#,
            # bsCollapsePanel(
            #   "Gene expression heatmap",
            #   wellPanel(
            #     div(
            #       class="row-fluid",
            #       div(class="span4",
            #           h5("Colors selection"),
            #           uiOutput("exprHeatmapColors"),
            #           h5("Heatmap clustering options"),
            #           uiOutput("exprHeatmapClustering")),
            #       div(
            #         class="span4",
            #         h5("Plot format"),
            #         radioButtons(
            #           "exprHeatmapFormat",
            #           "Select a format to save the plots:",
            #           c("PDF","PNG", "JPG"))
            #       ),
            #       div(class="span4",
            #           h5("Plot dimensions"),
            #           uiOutput("exprHeatmapDimensions"),
            #           uiOutput(paste("downloadExpr",
            #                          "HeatmapButton",
            #                          sep="")))
            #     )
            #   ),
            #   plotOutput("exprHeatmap", width="100%")
            # ),
            # bsCollapsePanel(
            #   "Tests for differential expression",
            #   wellPanel(
            #     uiOutput("diffExprTestsType"),
            #     uiOutput("downloadDiffExprTestsButton")
            #   ),
            #   br(),
            #   chartOutput("diffExprTests", "datatables")
            # ),
            # bsCollapsePanel(
            #   "Gene expression boxplot",
            #   wellPanel(
            #     h5("Gene selection"),
            #     uiOutput("selectGene")
            #   ),
            #   wellPanel(
            #     div(class="row-fluid",
            #         div(class="span6",
            #             h5("Plot format"),
            #             radioButtons(
            #               "exprBoxplotFormat",
            #               "Select a format to save the plots:",
            #               c("PDF","PNG", "JPG"))),
            #         div(class="span6",
            #             h5("Plot dimensions"),
            #             uiOutput("exprBoxplotDimensions"),
            #             uiOutput("downloadExprBoxplotButton")))
            #   ),
            #   br(),
            #   plotOutput("exprBoxplot")
            # )
          )
        ),
        value = "Further analyses"


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
