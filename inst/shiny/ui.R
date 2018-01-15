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
      h4(strong("Variables values data")),
      fileInput("expr", h5(paste("Please, select the text file (*.csv) containing the variables values data"),
                          img(src="images/info.png", title="Choose a file containing the table with the variables in the colums and the samples in the rows. The conditions of measure rows have to be indicated by a factor column")),
                accept=c("Comma-Seperated Values", "text/csv", ".csv"),placeholder = "Select a variable values file")
    ),
    uiOutput("classes"),
    uiOutput("factors"),
    wellPanel(
      h4(strong("Variable set database")),
      fileInput("geneSets",
                h5(paste("Load the Set of variables file (*.gmt)",
                      "describing the variables sets"),img(src="images/info.png", title="Choose a file containing the groups of variables. The format of the file that has to be inserted is explained in help section.")),
                accept=c('.gmt'),placeholder = "Select a variables groups file")
    ),
    # Parameters
    h3("2. Set parameters"),
    wellPanel(
      h4(strong("Variable sets size range")),
      h5("Enter the minimum and maximum variable set sizes."),
      uiOutput("minSize"),
      uiOutput("maxSize"),
      uiOutput("geneSetsCount")
    ),
    wellPanel(
      h3("Network construction"),
      h4(strong("a. Network type")),
      radioButtons(
        "networkType",
        h5("Select what network type will be compared"),
        c("Weighted"="weighted", "Unweighted"="unweighted")
      ),
      h4(strong("b. Method for network inference")),
      uiOutput("correlationMeasure"),
      h4(strong("c. Method for selecting the network links")),
      selectInput(
        "associationMeasure",
        h5("Association measure"),
        c("Absolute correlation"="correlation", "1 - p-value"="pvalue",
          "1 - q-value"="qvalue")
      ),
      uiOutput("linkFormation"),
      conditionalPanel(
        "input.networkType=='weighted'",
        radioButtons(
          "edgeWeight",
          h4(strong("d. Enter a association value to define the network edges weights:")),
          c("Absolute correlation"="correlation", "1 - p-value"="pvalue",
            "1 - q-value"="qvalue")
        )
      )
    ),
    wellPanel(
      h4(strong("Method for networks comparison")),
      uiOutput("networkTest"),
      uiOutput("networkTestDescription"),
      uiOutput("networkTestOptions")
    ),
    wellPanel(
      h4(strong("Permutation test settings")),
      numericInput("numPermutations",
                   "Enter the number of label permutations:",
                   1000, min=0),
      numericInput("seed",
                   "Enter a seed to generate random permutations:", 0)
    ),
    h3("3. Run differential network analysis"),
    wellPanel(
      actionButton("start", "Start analysis", class = "btn-primary")
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
        div(
          class="span4",
          h5("Classes selection"),
          uiOutput("factorsToNetViz1"),
          uiOutput("factorsToNetViz2")
        ),
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
            dataTableOutput("vertexAnalysisTable"),
            bsCollapsePanel(
              "KEGG pathway visualization",
                fileInput("keggCodes", h5(paste("Please, select the codes table (*.csv) containing the variables names and coresponding Kegg code")),
                          accept=c("Comma-Seperated Values", "text/csv", ".csv"),placeholder = "Select a file"),
              div(
                class="span4",
                h5("Colors selection"),
                uiOutput("exprKeegMapColors"),
                radioButtons(
                  "selectingDataType",
                  "Select a what kind of data you are analysing:",
                  c("Genes or/and Proteins"="gene", "Metabolite"="cpd"))
              ),
              selectInput(
                "thresholdSelected",
                h5("Association measure"),
                c("p-value"="pvalue",
                  "q-value"="qvalue","Test statistics"="statistic")
              ),
              uiOutput("thrValue"),
              textInput("speciesID", label = h5("Write the code of the species that you want to analyze"), value = "hsa"),
              textInput("pathID", label = h5("Write the code of the pathway that you want to analyze"), value = "05200"),
              checkboxInput("keggNative", label = "Kegg Native plot", value = TRUE),
              downloadButton("downloadKeggMap","Save the Kegg Map")
            )
          ),
          ## NETWORK VISUALIZATION
          tabPanel(
            "Network visualization plots",
          bsCollapsePanel(
            "Plot settings",
            div(
              class="row-fluid",
              div(
                class="span4",
                h5("Colors selection"),
                uiOutput("heatmapColors")
              ),
              div(
                class="span4",
                h5("Plot format"),
                radioButtons(
                  "networkPlotFormat",
                  "Select a format to save the plots:",
                  c("PDF","PNG", "JPG"))
              ),
              div(
                class="span4",
                h5("Plot dimensions"),
                uiOutput("networkPlotDimensions")
              )
            )
          ),
          bsCollapsePanel(
            "Network visualization",
            conditionalPanel(
              "input.associationMeasure=='correlation'",
              checkboxInput("signedCorrelation",
                            "Show negative correlations")
            ),
            br(),
            div(class="row-fluid",
                div(class="span6", plotOutput("heatmapClass1"),
                    uiOutput("downloadNetworkPlot1Button")),
                div(class="span6", plotOutput("heatmapClass2"),
                    uiOutput("downloadNetworkPlot2Button"))),
            br(),
            wellPanel(
              h5("Association between two gene products"),
              uiOutput("selectGenes"),

              dataTableOutput("corr")
            )
          ),
          bsCollapsePanel(
            "Differences between the gene networks",
            uiOutput("heatmapDiffOptions"),
            br(),
            bsCollapsePanel(
              "Matrix of differences",
              plotOutput("heatmapDiff", width="50%"),
              uiOutput("downloadNetworkDiffPlotButton")),
            bsCollapsePanel(
              "List of gene association degrees",
              wellPanel(
                uiOutput("absDiffType"),
                uiOutput("downloadAbsDiffButton")
              ),
              br(),
              dataTableOutput("corAbsDiff")#, "datatables")
              # chartOutput("corAbsDiff")#, "datatables")
            )
          )
          ),
          tabPanel(
            "Gene set properties",
            wellPanel(
              h5("Gene set network topological properties"),
              uiOutput("networkScore"),
              uiOutput("networkScoreOptions"),
              uiOutput("networkScoresComparison")
            )
          ),
          tabPanel(
            "Gene scores",
            wellPanel(
              h5("Gene scores"),
              uiOutput("geneScore"),
              uiOutput("geneScoresType"),
              uiOutput("downloadGeneScoresButton")
            ),
            br(),
            dataTableOutput("geneScoresComparison")
            # chartOutput("geneScoresComparison", "datatables")
          ),
          tabPanel(
            "Gene expression analysis",
            bsCollapsePanel(
              "Gene expression heatmap",
              wellPanel(
                div(
                  class="row-fluid",
                  div(class="span4",
                      h5("Colors selection"),
                      uiOutput("exprHeatmapColors"),
                      h5("Heatmap clustering options"),
                      uiOutput("exprHeatmapClustering")),
                  div(
                    class="span4",
                    h5("Plot format"),
                    radioButtons(
                      "exprHeatmapFormat",
                      "Select a format to save the plots:",
                      c("PDF","PNG", "JPG"))
                  ),
                  div(class="span4",
                      h5("Plot dimensions"),
                      uiOutput("exprHeatmapDimensions"),
                      uiOutput(paste("downloadExpr",
                                     "HeatmapButton",
                                     sep="")))
                )
              ),
              plotOutput("exprHeatmap", width="100%")
            ),
            bsCollapsePanel(
              "KEGG pathway visualization",
              fileInput("keggCodes2", h5(paste("Please, select the codes table (*.csv) containing the variables names and coresponding Kegg code")),
                        accept=c("Comma-Seperated Values", "text/csv", ".csv"),placeholder = "Select a file"),
              div(
                class="span4",
                h5("Colors selection"),
                uiOutput("exprKeegMapColors2"),
                radioButtons(
                  "selectingDataType2",
                  "Select a what kind of data you are analysing:",
                  c("Genes or/and Proteins"="gene", "Metabolite"="cpd"))
              ),
              selectInput(
                "thresholdSelected2",
                h5("Association measure"),
                c("p-value"="pvalue",
                  "q-value"="qvalue","Test statistics"="statistic")
              ),
              uiOutput("thrValue2"),
              textInput("speciesID2", label = h5("Write the code of the species that you want to analyze"), value = "hsa"),
              textInput("pathID2", label = h5("Write the code of the pathway that you want to analyze"), value = "05200"),
              checkboxInput("keggNative2", label = "Kegg Native plot", value = TRUE),
              downloadButton("downloadKeggMap2","Save the Kegg Map")
            ),
            bsCollapsePanel(
              "Tests for differential expression",
              wellPanel(
                uiOutput("diffExprTestsType"),
                uiOutput("downloadDiffExprTestsButton")
              ),
              br(),
              dataTableOutput("diffExprTests")
              # chartOutput("diffExprTests", "datatables")
            ),
            bsCollapsePanel(
              "Gene expression boxplot",
              wellPanel(
                h5("Gene selection"),
                uiOutput("selectGene")
              ),
              wellPanel(
                div(class="row-fluid",
                    div(class="span6",
                        h5("Plot format"),
                        radioButtons(
                          "exprBoxplotFormat",
                          "Select a format to save the plots:",
                          c("PDF","PNG", "JPG"))),
                    div(class="span6",
                        h5("Plot dimensions"),
                        uiOutput("exprBoxplotDimensions"),
                        uiOutput("downloadExprBoxplotButton"))),
                plotOutput("exprBoxplot", width="100%")
              )
            )
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
