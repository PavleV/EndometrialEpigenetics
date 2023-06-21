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
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application
shinyApp(ui = ui, server = server)
