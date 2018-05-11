library(shiny)
#library(shinyWidgets)
library(plotly)
library(reshape2)

ui <- fluidPage(
    
    #sidebarLayout(
        #sidebarPanel(
    
    fileInput(inputId='data', 'choose a data file', multiple=FALSE, accept=NULL),
    div(style="display:inline-block", selectInput(inputId="what", label="Choose what to plot",
                        c("Countings" = "counts", "Frequencies"="freqs"),
                        width="200px")),
    div(style="display:inline-block", numericInput(inputId="species", label="Species number ...", value=1, min=1, max=NA, step=1, width="200px")),
    div(style="display:inline-block", actionButton("start", "Start")),
    
        #),
        #mainPanel(
            plotlyOutput('barplot')
        #)
    #)
)

server <- function (input, output) {
    
    # Data
    data <- reactive({
        if (is.null(input$data)) {return(NULL)}
        dataframe <-read.table(input$data$datapath, header=TRUE, dec=".", sep=",", row.names=1)
        
        # Exclude transitions X -> X (plot only X -> Y)
        skip <- c()
        for (i in seq(1, ncol(dataframe),by=(ncol(dataframe)/sqrt(ncol(dataframe))+1))) {
            skip <- c(skip, i)
        }
        
        # values of skip to iterate ('skip'+1 without it's last value)
        starts <- skip[1:(length(skip)-1)]+1
        
        # Columns reordering, allows to have X -> Y and Y -> X bars at the same positions on the graphs
        # Patterns for iteration had not been easy to find
        reordered <- c()
        l <- length(starts) - 1
        
        for (i in 1:length(starts)) {
            step <- length(starts)
            substep <- 1
            for (j in seq(starts[i],(starts[i]+l))) {
                couple <- c(j, j+step)
                reordered <- c(reordered, couple)
                substep <- substep +1
                step <- length(starts) * substep
            }
            l <- l-1
        }
        
        final_data <- dataframe[,reordered]
        rm(list=c('skip','starts','i','j','l','couple','step','substep','reordered'))
        return(final_data)
    })
    
    observe({
        
        if(input$start != 0) {
 
            if (input$what == "counts") {
                data_plots <- data()[seq(1, nrow(data()), 3),]
            } else if (input$what == "freqs") {
                data_plots <- data()[seq(2, nrow(data()), 3),]
            }
            
            pvalues <- data()[seq(3, nrow(data()), 3),]
            
            row.names(data_plots) <- substrLeft(row.names(data_plots),5)
            row.names(pvalues) <- substrLeft(row.names(pvalues),5)
            
            # Select the couple of species
            spec1 <- data_plots[input$species,]
            spec1p <- pvalues[input$species,]
            
            # Melting dataframe
            m <- melt(spec1, variable="Type")
            mp <- melt(spec1p, variable="Type")
            
            # Split odd and even lines, one being the bars on top and the others the bars below.
            mAbove <- m[seq(1,nrow(m), by=2),]
            mAbove <- droplevels(mAbove)
            mBelow <- m[seq(2,nrow(m), by=2),]
            mBelow <- droplevels(mBelow)
            
            # Bind the splitted liens in the case we want to filter easily 0 values
            bind <- cbind(mAbove,mBelow)
            # exclude lines where both transi = 0
            bind.filt <- bind[rowSums(bind[,c(2,4)])!=0,]
            bind.filt <- droplevels(bind.filt)
            
            # Find max y axis value
            max_y <- max(mBelow[,2], mAbove[,2])
            
            output$barplot <- renderPlotly({
                pa4 <- plot_ly() %>%
                    add_bars(x = bind.filt[,1], 
                             y = bind.filt[,2], 
                             name='From X to Y',
                             showlegend = TRUE,
                             marker=list(color = 'rgb(102, 178, 1255)')) %>%
                    layout(yaxis=list(range=c(0,max_y)))
                
                pb4 <- plot_ly() %>%
                    add_bars(x = bind.filt[,3], 
                             y = bind.filt[,4], 
                             name='From Y to X',
                             showlegend = TRUE,
                             marker=list(color = 'rgb(102, 255, 102)')) %>%
                    layout(title=paste("Transitions between", strsplit(rownames(data_plots)[input$species],">")[[1]][1], "to", strsplit(rownames(data_plots)[input$species],">")[[1]][2], sep=" "),
                           yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
                           xaxis = list(side='top'))
                
                subplot(pa4, pb4, nrows = 2, margin=0.1)
            })
        }
        
    })
    
}

# Simplify row names
substrLeft <- function(x, n){
    sapply(x, function(xx)
        substr(xx, 0, n)
    )
}

shinyApp(ui=ui, server=server)