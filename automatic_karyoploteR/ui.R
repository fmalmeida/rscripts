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
library(shinycssloaders)

# Shiny options
options(
    repos = BiocManager::repositories(),
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
            
            width = "2",
            
            tabsetPanel(
                tabPanel("Input files",
                         
                         br(), # newline
                         
                         ####################
                         ### Loading data ###
                         ####################
                         
                         # Genome in BED format
                         fileInput(
                             inputId = "genomeBED",
                             label   = "BED file of input genome (Chr, Start, End).",
                             accept  = ".bed"
                         ),
                         
                         # Genome window
                         uiOutput("genomeWindow"),
                         
                         # Complete GFF
                         fileInput(
                             inputId = "mainGFF",
                             label   = "Genome GFF that will be used to calculate feature density",
                             accept  = c(".gff", ".gff3")
                         ),
                         
                         # Select feature
                         uiOutput("feature_select"),
                         
                         # Regions GFF
                         fileInput(
                             inputId = "regionGFFs",
                             label   = "GFFs (subsets) to be drawn as regions (Multiple inputs allowed)",
                             multiple = TRUE
                         )),
                
                tabPanel("Plot parameters",
                         
                         br(), # newline
                         
                         tabsetPanel(
                          
                          # Params for title adjustment
                          tabPanel("Plot title",
                                    
                                    br(),
                                    
                                    # plot title
                                    textInput(
                                        inputId = "plot-title",
                                        label = "Please enter the desired plot title",
                                        value = "karyoploteR template"
                                    ),
                                    sliderInput(
                                        inputId = "titlecex",
                                        label   = "Title font cex (size)",
                                        min     = 0.1,
                                        max     = 2,
                                        value   = 0.9,
                                    )),
                           
                           # Params for chr ideogram adjustment
                           tabPanel("Chr ideogram",
                                    
                                    br(),
                                    
                                    sliderInput(
                                        inputId = "chrcex",
                                        label   = "Chr label font cex (size)",
                                        min     = 0.1,
                                        max     = 1,
                                        value   = 0.5,
                                    ),
                                    sliderInput(
                                        inputId = "chroffset",
                                        label   = "Chr label positioning (x axis)",
                                        min     = -10,
                                        max     = 10,
                                        value   = -0.1,
                                        step    = 0.1
                                    ), 
                                    sliderInput(
                                        inputId = "ideoheight",
                                        label   = "Ideogram height",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 300
                                    )),
                           
                          # Params for tick controlling
                          tabPanel("Handling ticks",
                                    
                                    br(),
                                    
                                    sliderInput(
                                        inputId = "tickdistance",
                                        label   = "Distance between ticks",
                                        min     = 0,
                                        max     = 1*10^10,
                                        value   = 2.5*10^6,
                                        step    = 500
                                    )),
                          
                          # Params for title adjustment
                          tabPanel("Plot colors",
                                   
                                   br(),
                                   
                                   # plot title
                                   textInput(
                                       inputId = "data_panel_colors",
                                       label = "Which colors (in #HEX or not) to use in data panels?"
                                   )),
                           
                           # KaryoploteR data.panels
                           tabPanel("karyoploteR panels",
                                    
                                    br(),
                                    
                                    sliderInput(
                                        inputId = "plot-type",
                                        label   = "Plot type",
                                        min     = 1,
                                        max     = 2,
                                        value   = 1,
                                        step = 1
                                    ),
                                    
                                    sliderInput(
                                        inputId = "bottommargin",
                                        label   = "Bottom margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 500
                                    ),
                                    sliderInput(
                                        inputId = "topmargin",
                                        label   = "Top margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 500
                                    ),
                                    sliderInput(
                                        inputId = "data1height",
                                        label   = "Data panel 1 height",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 300
                                    ),
                                    sliderInput(
                                        inputId = "data1outmargin",
                                        label   = "Data panel 1 out margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 500
                                    ),
                                    sliderInput(
                                        inputId = "data1inmargin",
                                        label   = "Data panel 1 in margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 200
                                    ),
                                    sliderInput(
                                        inputId = "data2height",
                                        label   = "Data panel 2 height",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 300
                                    ),
                                    sliderInput(
                                        inputId = "data2outmargin",
                                        label   = "Data panel 2 out margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 500
                                    ),
                                    sliderInput(
                                        inputId = "data2inmargin",
                                        label   = "Data panel 2 in margin",
                                        min     = 1,
                                        max     = 2000,
                                        value   = 200
                                    ))
                           
                         ))
            )
            
        ),
        mainPanel(
            width = 9,
            verbatimTextOutput("test"),
            uiOutput("delimiter"),
            withSpinner(plotOutput("karyotype", width = "100%", height = 1000), color="#0dc5c1")
        )
    )
))
