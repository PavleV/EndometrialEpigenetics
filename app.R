#
# Shiny web application for visualising scATAC-seq data
#
#

library(shiny)
source("functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Endometrial Epigenetics"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(

            # Choose how to deal with gene names and genomic locations
            radioButtons("link", label = h4("Link genes and genomic coordinates"),
                         choices = list("By gene name" = "bygene", "By genomic location" = "bylocation", "Separate genes and genomic coordinates" = "separate"),
                         selected = "bygene"),

            # select Biopsies
            checkboxGroupInput(inputId="which", label=h4("Select samples for RNA plot:"),
                               choices = c("Sample 1"="S169","Sample 2"="S170","Sample 3"="S506", "Sample 4"="S508"),
                               selected = c("S169","S170","S506","S508")),

            # Choose gene name
            conditionalPanel(
                condition = "input.link != 'bylocation'",
                textInput("genename", label = h4("Gene Name(s)"), value = "SCARA5", placeholder="Enter one or more gene names..."),
            ),

            conditionalPanel(
                condition = "input.link == 'bygene'",
                # select gene regions
                radioButtons("gene_region", label = h4("Gene regions"),
                             choices = list("TSS" = "TSS", "Whole gene" = "Whole", "TTS" = "TTS"),
                             selected = "Whole"),
                #numeric input of up and downstream
                numericInput("add_up", label = h4("Add upstream bases"), value = 0),
                numericInput("add_down", label = h4("Add downstream bases"), value = 0)
            ),

            # choose a genomic region

            conditionalPanel(
                condition = "input.link != 'bygene'",
                textInput("genome_coord", label = h4("Genomic region (hg19)"), value = "chr8:27696656-27881112", placeholder="Enter chromosomal location...")

            )


        ),

        # Show the generated plots
        mainPanel(
            plotOutput("rnaPlot"),
            plotOutput("atacPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$rnaPlot <- renderPlot({

        if(!is.null(input$genename) & input$link != 'bylocation')
        {
            # split multiple entries into character vector
            geneNames.to.plot <- str_split_1(input$genename, pattern=regex("[:,_-[:space:]]"))
        }
        if(!is.null(input$genome_coord) & input$link == 'bylocation')
        {
            genomic_coordinates <- extractCoord(input$genome_coord)
            geneNames.to.plot <- geneKey.ranges$mcols.GeneName[subjectHits(findOverlaps(genomic_coordinates,geneKey.ranges))]
        }

        # Check if geneNames.to.plot is NULL or empty
        #if (is.null(geneNames.to.plot) || length(geneNames.to.plot) == 0) {
        # Show an error message
        #    return(plot_error_message("Please enter a valid gene name."))
        #}

        req(geneNames.to.plot)
        # draw the line plot of RNA expression based on genename/ensemblID and selected biopsies
        plotRNA.FUN(mydata = AH_EL_RNA_ALLREPS, geneName=geneNames.to.plot, ensemblID=NULL, biopsies = input$which)
    })


    output$atacPlot <- renderPlot({

        if(!is.null(input$genename) & input$link == 'bygene')
        {
            gene_names <- str_split_1(input$genename, pattern=regex("[:,_-[:space:]]"))
            genomic_coordinates <- extractCoordfromGene(genes=gene_names, geneCoordKey=geneKey.ranges, add.upstream = input$add_up, add.downstream = input$add_down, gene.region = input$gene_region)

        }
        if(!is.null(input$genome_coord) & input$link != 'bygene')
        {
            genomic_coordinates <- extractCoord(input$genome_coord)

        }


        # Check if genomic_coordinates is NULL or empty
        if (is.null(genomic_coordinates) || length(genomic_coordinates) == 0) {
            # Show an error message
            return("Invalid genomic coordinates. Please enter valid coordinates.")
        }

        req(genomic_coordinates)

        # draw the line plot based on peaks within the specified genomic coordinates

        plotATAC.FUN(mydata=ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates = genomic_coordinates)
    })

}

# Run the application
shinyApp(ui = ui, server = server)
