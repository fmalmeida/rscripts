#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)

# Shiny options
options(
    shiny.maxRequestSize = 50*1024^2 # Input size to 50Mb
)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # App theme
    theme = shinytheme("cosmo"),

    # Application title
    titlePanel("A shiny app to automatize karyoploteR"),

    sidebarLayout(
        sidebarPanel(
            
            tabsetPanel(
                tabPanel("Input files",
                         
                         br(), # newline
                         
                         ####################
                         ### Loading data ###
                         ####################
                         
                         # Genome in BED format
                         fileInput(
                             inputId = "genomeBED",
                             label   = "BED file of input genome (Chr, Start, End)."
                         ),
                         
                         # Genome window
                         uiOutput("genomeWindow"),
                         
                         # Complete GFF
                         fileInput(
                             inputId = "mainGFF",
                             label   = "Genome GFF that will be used to calculate feature density"
                         ),
                         
                         # Select feature
                         textInput(
                             inputId = "feature", 
                             label   = "Feature type to plot density in the plot",
                             value   = "gene"
                         ),
                         
                         # Regions GFF
                         fileInput(
                             inputId = "regionGFFs",
                             label   = "GFFs (subsets) to be drawn as regions (Multiple inputs allowed)", multiple = TRUE
                         )), 
                
                tabPanel("Plot parameters", 
                         
                         br(), # newline
                         
                         # plot title
                         textInput(
                             inputId = "plot-title",
                             label = "Please enter the desired plot title",
                             value = "karyoploteR template"
                         ))
            )
            
        ),
        mainPanel(
            uiOutput("Overview"),
            br(),
            verbatimTextOutput("features"),
            uiOutput("delimiter"),
            plotOutput("karyotype")
        )
    )
))
