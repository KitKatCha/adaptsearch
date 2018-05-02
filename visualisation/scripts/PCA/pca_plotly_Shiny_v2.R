library(shiny)
library(shinythemes)
library(plotly)
library(FactoMineR)
library(DT)

ui <- navbarPage(theme=shinytheme('sandstone'), 'PCA analysis',
                 
    tabPanel('Main',
             sidebarLayout(
                 sidebarPanel(
                     fileInput("input_file", "choose a data file", multiple=FALSE, accept=NULL, width="300px"),
                     actionButton("start", "Start"),
                     div(style="display:inline-block; width:150px", selectInput(inputId="what", label="Choose what to plot", c("Countings" = "counts", "Frequencies"="freqs"))),
                     numericInput(inputId="axisX", label="Axis X", value=1, min=1, max=NA, step=1, width="75px"),
                     div(style="display:inline-block; width:75px", numericInput(inputId="axisY", label="Axis Y", value=2, min=1, max=NA, step=1)),
                     selectInput(inputId="what_inds", label="Individuals color range", c("cos2" = "cos2", "Contribution" = "contrib"), width="150px"),
                     selectInput(inputId="what_vars", label="Variables color range", c("Correlation"="cor", "cos2" = "cos2", "Contribution" = "contrib"), width="150px"),
                     numericInput(inputId="ind_colr", label="Axis for individuals color range", value=1, min=1, max=NA, step=1, width="210"),
                     numericInput(inputId="var_colr", label="Axis for variables color range", value=1, min=1, max=NA, step=1, width="210")
                ),
                mainPanel(
                    div(style="display:inline-block; width:150px", checkboxInput("show_labels", label="Show labels", value=TRUE)),
                    div(style="display:inline-block; width:200px", checkboxInput("hide_data", label="Hide unrepresentative data", value=FALSE)),
                    conditionalPanel(condition ="input.hide_data",
                                     selectInput(inputId="hide_what", label="Reference", c("cos2" = "cos2", "Contribution" = "contrib")),
                                     sliderInput("hide_range", "Hide values under (%of the max value)...", step=5, min=0, max=100, value=0)
                                     ),
                    h4('Individuals'),
                    plotlyOutput('PCA_ind'),
                    h4('Variables'),
                    plotlyOutput('PCA_var'),
                    h4('Biplot (on scaled coordinates)'),
                    plotlyOutput('PCA_biplot'),
                    h4('Explained variance per dimension'),
                    plotlyOutput('PCA_eigen')
                )
             )
    ),
    
    tabPanel('Individuals data',
             mainPanel(
                 h4('Cos2'),
                 dataTableOutput('indcos2'),
                 h4('Contribution'),
                 dataTableOutput('indcontrib')
             )
    ),
    
    tabPanel('Variables data',
             mainPanel(
                 h4('Cos2'),
                 dataTableOutput('varcos2'),
                 h4('Contribution'),
                 dataTableOutput('varcontrib'),
                 h4('Correlation'),
                 dataTableOutput('varcor')
             )
    )
)

server <- function (input, output) {
    
    # get data and compute PCA
    data <- reactive({
        if (is.null(input$input_file)) {return(NULL)}
        data <-read.table(input$input_file$datapath, header=TRUE, dec=".", sep=",", row.names=1)
        
        if (input$input_file$name != "gc_and_others.csv") {
            if (input$what == "counts") {selected_data <- data[seq(1, nrow(data), 2),]} 
            else if (input$what == "freqs") {selected_data <- data[seq(2, nrow(data), 2),]}
            
            if (grepl('transitions', input$input_file$name)) {
                row.names(selected_data) <- substrLeft(row.names(selected_data),3) # faire plutôt un strsplit ?
            } else {
                row.names(selected_data) <- substrLeft(row.names(selected_data),1)
            }
            
        } else {
            selected_data <- data
        }
        
        res.pca <- PCA(selected_data, scale.unit=TRUE, graph=F, axes=c(1,2))
        # indx <- (apply(res.pca$cos2, 2, function(x) any(is.na(x))))
        return(res.pca)
    })
    
    observe({
        if(input$start != 0) {
            
            # Plotting individuals
            output$PCA_ind <- renderPlotly({
                ind_coord <- as.data.frame(data()$ind$coord)
                if (input$what_inds == "cos2") {colors <- as.data.frame(data()$ind$cos2)}
                else if (input$what_inds == "contrib") {colors <- as.data.frame(data()$ind$contrib)}
                if (input$show_labels) {mod="markers+text"}
                else {mod="markers"}
                
                # Essayer avec which pour correspondance entre tables
                # if (input$hide_data) {
                #     max <- max(colors)
                #     ind_coord <- ind_coord[colors[,input$ind_colr]>(input$range*max/100)]
                # }
                
                p <- plot_ly(ind_coord,
                             x=ind_coord[,input$axisX],
                             y=ind_coord[,input$axisY],
                             type = 'scatter',
                             text=rownames(ind_coord),
                             textposition='top',
                             mode=mod,
                             color=colors[,input$ind_colr],
                             colors="OrRd",
                             marker=list(symbol=27, size=11))
                
                p <- layout(p, title = "PCA on individuals",
                            xaxis = list(title = input$axisX),
                            yaxis = list(title = input$axisY))
                
                p
            })
            
            # Plotting variables
            output$PCA_var <- renderPlotly({
                var_coord <- as.data.frame(data()$var$coord)
                if (input$what_vars == "cor") {colors <- as.data.frame(data()$var$cor)}
                else if (input$what_vars == "cos2") {colors <- as.data.frame(data()$var$cos2)}
                else if (input$what_vars == "contrib") {colors <- as.data.frame(data()$var$contrib)}
                if (input$show_labels) {mod="markers+text"}
                else {mod="markers"}
                
                p2 <- plot_ly(var_coord,
                              x=var_coord[,input$axisX],
                              y=var_coord[,input$axisY],
                              type = 'scatter',
                              text=rownames(var_coord),
                              textposition='top',
                              mode=mod,
                              color=colors[,input$var_colr],
                              colors="Blues",
                              marker=list(symbol=4, size=11))
                
                p2 <- layout(p2, title = "PCA on variables",
                             xaxis = list(title = input$axisX, range=c(-1.1,1.1)),
                             yaxis = list(title = input$axisY, range=c(-1.1,1.1)))
                
                p2
            })
            
            # Plotting both individuals and variables
            output$PCA_biplot <- renderPlotly({
                
                ind_coord <- scale(as.data.frame(data()$ind$coord), center=TRUE, scale=TRUE)
                if (input$what_inds == "cos2") {colors <- as.data.frame(data()$ind$cos2)}
                else if (input$what_inds == "contrib") {colors <- as.data.frame(data()$ind$contrib)}
                
                var_coord <- scale(as.data.frame(data()$var$coord), center=TRUE, scale=TRUE)
                if (input$what_vars == "cor") {colors2 <- as.data.frame(data()$var$cor)}
                else if (input$what_vars == "cos2") {colors2 <- as.data.frame(data()$var$cos2)}
                else if (input$what_vars == "contrib") {colors2 <- as.data.frame(data()$var$contrib)}
                
                if (input$show_labels) {mod="markers+text"}
                else {mod="markers"}
                
                biplot <- plot_ly(
                    x=ind_coord[,input$axisX],
                    y=ind_coord[,input$axisY],
                    type = 'scatter',
                    text=rownames(ind_coord),
                    textposition='top',
                    mode=mod,
                    #color=colors[,input$ind_colr],
                    #colosr="OrRd
                    marker=list(color="rgb(225, 50, 50)", symbol=27, size=11)) %>%
                    
                    add_trace(x=var_coord[,input$axisX],
                              y=var_coord[,input$axisY],
                              type = 'scatter',
                              text=rownames(var_coord),
                              textposition='top',
                              mode=mod,
                              #color=colors2[,input$var_colr],
                              #colors="BuGn",
                              marker=list(color="rgb(86,135,212)", symbol=4, size=11))
                
                biplot <- layout(biplot, title="PCA - bi-plot",
                                 showlegend=FALSE,
                                 xaxis = list(title=input$axisX),
                                 yaxis = list(title=input$axisY))
                
                # p6 <- plot_ly(showscale = TRUE) %>%
                #     add_surface(z = ~z1, cmin = min(z1), cmax = max(z2), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
                #     add_surface(z = ~z2, cmin = min(z1), cmax = max(z2), colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
                #     layout(title="SURFACE 1 + SURFACE 2\n(Distinct colorscales as defined)", scene = list(zaxis = list(range = c(min(z1),max(z2)))))
                
                biplot
            })
            
            output$PCA_eigen <- renderPlotly({
                eig <- as.data.frame(data()$eig)
                eig_names <- factor(rownames(eig), levels = rownames(eig))
                eigen <- plot_ly(eig) %>%
                    add_trace(x=eig_names,
                              y=eig[,2],
                              type='bar',
                              name="Explained variation") %>%
                    add_trace(x=eig_names,
                              y=eig[,3],
                              type='scatter',
                              mode='lines',
                              line=list(color="rgb(225, 50, 50)"),
                              name="Cumulative explained variation") %>%
                    layout(legend = list(orientation = 'h'))
                
                eigen
            })
            
            output$indcos2 <- renderDataTable(data()$ind$cos2)
            output$indcontrib <- renderDataTable(data()$ind$contrib)
            output$varcos2 <- renderDataTable(data()$var$cos2)
            output$varcor <- renderDataTable(data()$var$cor)
            output$varcontrib <- renderDataTable(data()$var$contrib)
        }
    })
}

substrLeft <- function(x, n){
    sapply(x, function(xx)
        substr(xx, 0, n)
    )
}

shinyApp(ui=ui, server=server)