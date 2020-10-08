# A shiny app to execute karyoploteR packages with user's gff/bed inputs
#
# Author: Felipe Marques de Almeida

#################
### Tutorials ###
#################

# https://bernatgel.github.io/karyoploter_tutorial/#Tutorial
# https://www.bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
# https://www.rapidtables.com/web/color/RGB_Color.html

######################
### Load libraries ###
######################
suppressMessages(library(shiny))
suppressMessages(library(karyoploteR))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(plyranges))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(randomcoloR))

########################
### Useful functions ###
########################

# Load genome chr from bed (3 cols only)
# It is possible to filter chrs based on min and max sizes
# Default: All chrs
grFromBed <- function(bed, minSize=1, maxSize=NULL) {
    
    if (!is.null(maxSize)) {
        # Open bed
        genome.bed <- read.csv(sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE,
                               file = bed,
                               col.names = c('Chr', 'Start', 'End')) %>%
            dplyr::filter((End - Start) >= minSize & (End - Start) <= maxSize)
        
        # Set BED start to one
        genome.bed$Start <- 1
        
        # Return GRanges
        return(
            with(genome.bed, GRanges(Chr, IRanges(Start, End)))
        )
        
    } else {
        
        # Open bed
        genome.bed <- read.csv(sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE,
                               file = bed,
                               col.names = c('Chr', 'Start', 'End'))
        
        # Set BED start to one
        genome.bed$Start <- 1
        
        # Return GRanges
        return(
            with(genome.bed, GRanges(Chr, IRanges(Start, End)))
        )
        
    }
}
readTmpGff <- function(input) {
    f <- file.path(tempdir(), "tmp.gff")
    
    write.table(
        read.csv(sep = "\t", file = input,
                 quote = "", header = F, comment.char = "#"),
        file = f, quote = F, sep = "\t", row.names = F, col.names = F
    )
    
    return(f)
}

# Define server logic
shinyServer(function(input, output) {
    
    ##################################
    ### List GFF possible features ###
    ##################################
    output$feature_select <- renderUI({
        
        # Wait for GFF
        req(input$mainGFF)
        
        # Generate
        gff.gr  <- import.gff(input$mainGFF$datapath)
        ft_list <- unique(gff.gr$type)
        selectInput("features", label = h3("Select feature for density plot"), 
                    choices = ft_list,
                    selected = 1)
    })
    
    #############################################
    ### Create window slider for density plot ###
    #############################################
    output$density_window <- renderUI({
        
        # Wait for GFF
        req(input$mainGFF)
        
        # Window
        textInput(
            inputId = "density_window",
            label   = "Set the window size for feature density plot",
            value   = 500
        )
    })
    
    ###########################################
    ### Create color input for density plot ###
    ###########################################
    output$density_color <- renderUI({
        
        # Wait for GFF
        req(input$mainGFF)
        
        # Color
        textInput(
            inputId = "density_color",
            label   = "Set desired color for density plot",
            value   = "#313D7C"
        )
        
    })
    
    ######################
    ### Results header ###
    ######################
    output$delimiter <- renderUI({
        
        HTML(paste0(
            "<b><h3><span style='font-weight: bold;'>Resulting karyotypes will appear here</span></h3></b>",
            sep = ""))
    })
    
    ################################
    ### Checking genome BED file ###
    ################################
    output$genomeWindow <- renderUI({
        req(input$genomeBED) # This makes the function wait for the input
        genome_bed <- read.csv(input$genomeBED$datapath, sep = "\t", header = FALSE)
        max_val <- max(genome_bed$V3)
        
        # UI
        # Genome window
        sliderInput(
            inputId = "genomeWindow",
            label = "Select chr size threasholds for plot (min and max sizes)",
            min = 0,
            max = max_val,
            value = c(0, max_val),
            step = 1000
        )
    })
    
    ################################
    ### Render ideo color option ###
    ################################
    output$datapanel_color <- renderUI({
        req(input$genomeBED) # This makes the function wait for the input
        
        # UI
        textInput(
            inputId = "datapanel_color",
            label   = "Set color for data panels background",
            value   = "#F7F7F7"
        )
    })
    
    #####################################
    ### Generate the body of the plot ###
    ### already from the bed file     ###
    #####################################
    custom_plot <- reactive({
        # Wait for regions to plot
        req(input$genomeBED)
        
        # Plot parameters
        pp <- getDefaultPlotParams(plot.type = input$`plot-type`)
        pp$bottommargin   <- input$bottommargin
        pp$topmargin      <- input$topmargin
        pp$ideogramheight <- input$ideoheight
        pp$data1height    <- input$data1height
        pp$data1outmargin <- input$data1outmargin
        pp$data1inmargin  <- input$data1inmargin
        pp$data2height    <- input$data2height
        pp$data2outmargin <- input$data2outmargin
        pp$data2inmargin  <- input$data2inmargin
        pp$data1max = 1
        
        # Load genome bed as GRanges
        genome.gr <- grFromBed(input$genomeBED$datapath, 
                               minSize = input$genomeWindow[1], 
                               maxSize = input$genomeWindow[2])
        
        # Create karyotype body
        kp <- plotKaryotype(plot.type = input$`plot-type`, main = input$`plot-title`, cex = input$titlecex, 
                            genome = genome.gr, plot.params = pp, labels.plotter = NULL)
        kpAddChromosomeNames(kp, cex = input$chrcex, xoffset = input$chroffset, yoffset = -0.1) # The options yoffset and xoffset control the position where the names will appear
        
        if (is.null(input$datapanel_color)) {
            color <- "#F7F7F7"
        } else {
            color <- as.character(input$datapanel_color)
        }
        kpDataBackground(kp, color = color)
        kpAddBaseNumbers(kp, tick.dist = as.integer(input$tickdistance), tick.len = input$ticklen, cex = input$tickcex)
        
        # Plot features density
        if (!is.null(input$mainGFF) & !is.null(input$features)) {
            
            # Load main gff as GRanges
            gff.gr <- import.gff(input$mainGFF$datapath)
            
            # Filter GFF features for density plot
            filtered.gff.gr <- gff.gr[gff.gr$type == as.character(input$features)]
            
            # Gene density
            if (is.null(input$density_color)) {
                color <- randomColor()
            } else {
                color <- "#313D7C"
            }
            
            if (is.null(input$density_window)) {
                d_window <- 500
            } else {
                d_window <- as.integer(input$density_window)
            }
            
            kpPlotDensity(kp, filtered.gff.gr, window.size = d_window, 
                          data.panel = "ideogram", col=color, border=color,
                          r0 = 0.1, r1 = 0.85)
        }
        
        # Plot regions data
        if (!is.null(input$regionGFFs)) {
            
            # Calculate step for plotting data in data.panel
            step     = 1 / length(input$regionGFFs[,1])
            r0_val   = 0
            r1_val   = step - 0.05
            
            # Plot
            for (i in 1:length(input$regionGFFs[,1])) {
                
                if (input$data_panel_colors == "" || length(as.list(strsplit(input$data_panel_colors, ",")[[1]])) != length(input$regionGFFs[,1])) {
                    color = randomColor()
                } else {
                    color = as.list(strsplit(input$data_panel_colors, ",")[[1]])[[i]]
                }
                
                region <- import.gff(input$regionGFFs[[i, 'datapath']])
                label  <- input$regionGFFs[[i, 'name']]
                kpPlotRegions(kp, data=region, col=color, border=color, r0=r0_val, r1=r1_val)
                kpAddLabels(kp, labels = label, r0=r0_val, r1=r1_val, cex = .5)
                
                # Sum
                r0_val = r0_val + step
                r1_val = r1_val + step
            }
            
        }
        
        print(kp)
    })
    
    output$karyotype <- renderPlot({
        custom_plot()
    })
    
    # output$downloadPlot <- downloadHandler(
    #     filename = function() { paste("teste", '.png', sep='') },
    #     content = function(file) {
    #         png(file)
    #         custom_plot()
    #         dev.off()
    #     })
})
