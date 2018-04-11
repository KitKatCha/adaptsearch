#load packages
library(shiny)
library(shinyWidgets)
library(plotly)
library(FactoMineR)

ui <- bootstrapPage(


    # Parameters
    fileInput("input_file", "choose a data file", multiple=FALSE, accept=NULL),
    actionButton("start", "Start"),
    div(style="display:inline-block", selectInput(inputId="what", label="Choose what to plot",
        c("Countings" = "counts", "Frequencies"="freqs"),
        width="150px")),
    div(style="display:inline-block", numericInput(inputId="axisX", label="Axis X", value=1, min=1, max=NA, step=1, width="50px")),
    div(style="display:inline-block", numericInput(inputId="axisY", label="Axis Y", value=2, min=1, max=NA, step=1, width="50px")),
    div(style="display:inline-block", selectInput(inputId="what_inds", label="Individuals color range", 
        c("cos2" = "cos2", "Contribution" = "contrib"),
        width="160px")),
    div(style="display:inline-block", selectInput(inputId="what_vars", label="Variables color range", 
        c("Correlation"="cor", "cos2" = "cos2", "Contribution" = "contrib"),
        width="150px")),
    div(style="display:inline-block", numericInput(inputId="ind_colr", label="Axis for individuals color range", value=1, min=1, max=NA, step=1, width="230px")),
    div(style="display:inline-block", numericInput(inputId="var_colr", label="Axis for variables color range", value=1, min=1, max=NA, step=1, width="200px")),
    # switchInput(
    #     inputId = "plot_details",
    #     label = strong("Show details"),
    #     value = FALSE,
    #     labelWidth = "120px"
    # ),

    # UI
    fluidRow(
        column(
            6,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('PCA_ind')
            )
        ),
        column(
            6,
            div(style = "border-style: solid; border-width: thin;border-color: #000000",
                plotlyOutput('PCA_var')
            )
        )
    ),
    fluidRow(
        column(
            6,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('PCA_biplot')
            )
        ),
        column(
            6,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('PCA_eigen')
            )
        )
    ),
    fluidRow(
        column(
            6,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('indcos2')
            )
        ),
        column(
            6,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('indcontrib')
            )
        )
    ),
    fluidRow(
        column(
            4,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('varcos2')
            )
        ),
        column(
            4,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('varcontrib')
            )
        ),
        column(
            4,
            div(style = "border-style: solid; border-width: thin; border-color: #000000",
                plotlyOutput('varcor')
            )
        )
    )
)

server <- function (input, output) {
    
    # get data and compute PCA
    data <- reactive({
        if (is.null(input$input_file)) {return(NULL)}
        data <-read.table(input$input_file$datapath, header=TRUE, dec=".", sep=",", row.names=1)
        # if (input$data$name != "gc_and_others.csv") {
        #     toDelete <- seq(1, nrow(dat), 3)
        #     dat_good <- dat[toDelete,]
        # }
        if (input$what == "counts") {selected_data <- data[seq(1, nrow(data), 3),]} 
        else if (input$what == "freqs") {selected_data <- data[seq(2, nrow(data), 3),]}
        row.names(selected_data) <- substrLeft(row.names(selected_data),2)
        res.pca <- PCA(selected_data, scale.unit=TRUE, graph=F, axes=c(1,2))
        return(res.pca)
    })
    
    observe({
        if(input$start != 0) {

            # Plotting individuals
            output$PCA_ind <- renderPlotly({
                ind_coord <- as.data.frame(data()$ind$coord)
                if (input$what_inds == "cos2") {colors <- as.data.frame(data()$ind$cos2)}
                else if (input$what_inds == "contrib") {colors <- as.data.frame(data()$ind$contrib)}
                
                p <- plot_ly(ind_coord,
                             x=ind_coord[,input$axisX],
                             y=ind_coord[,input$axisY],
                             type = 'scatter',
                             text=rownames(ind_coord),
                             textposition='top',
                             mode="markers+text",
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
                else if (input$what_inds == "cos2") {colors <- as.data.frame(data()$var$cos2)}
                else if (input$what_vars == "contrib") {colors <- as.data.frame(data()$var$contrib)}

                p2 <- plot_ly(var_coord,
                              x=var_coord[,input$axisX],
                              y=var_coord[,input$axisY],
                              type = 'scatter',
                              text=rownames(var_coord),
                              textposition='top',
                              mode="markers+text",
                              color=colors[,input$var_colr],
                              colors="BuGn",
                              marker=list(symbol=4, size=11))

                p2 <- layout(p2, title = "PCA on variables",
                             xaxis = list(title = input$axisX),
                             yaxis = list(title = input$axisY))

                p2
            })

            # Plotting both individuals and variables
            output$PCA_biplot <- renderPlotly({
                
                ind_coord <- as.data.frame(data()$ind$coord)
                if (input$what_inds == "cos2") {colors <- as.data.frame(data()$ind$cos2)}
                else if (input$what_inds == "contrib") {colors <- as.data.frame(data()$ind$contrib)}
                
                var_coord <- as.data.frame(data()$var$coord)
                if (input$what_vars == "cor") {colors2 <- as.data.frame(data()$var$cor)}
                else if (input$what_inds == "cos2") {colors2 <- as.data.frame(data()$var$cos2)}
                else if (input$what_vars == "contrib") {colors2 <- as.data.frame(data()$var$contrib)}
                
                biplot <- plot_ly(
                    x=ind_coord[,input$axisX],
                    y=ind_coord[,input$axisY],
                    type = 'scatter',
                    text=rownames(ind_coord),
                    textposition='top',
                    mode="markers+text",
                    color=colors[,input$ind_colr],
                    colors="OrRd",
                    marker=list(symbol=27, size=11)) %>%
                    add_trace(x=var_coord[,input$axisX],
                              y=var_coord[,input$axisY],
                              type = 'scatter',
                              text=rownames(var_coord),
                              textposition='top',
                              mode="markers+text",
                              color=colors2[,input$var_colr],
                              colors="BuGn",
                              marker=list(symbol=4, size=11))

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
                eig <- data()$eig
                eigen <- plot_ly(eig) %>%
                    add_trace(x=row.names(eig),
                              y=eig[,2],
                              type='bar',
                              name="Explained variation") %>%
                    add_trace(eig,
                              x=row.names(eig),
                              y=eig[,3],
                              type='scatter',
                              mode='lines',
                              name="Cumulative explained variation") %>%
                    layout(legend = list(orientation = 'h'))

                eigen
            })
            
            output$indcos2 <- renderPlotly({
                i2 <- plot_ly(x=colnames(data()$ind$cos2),
                              y=row.names(data()$ind$cos2),
                              z=data()$ind$cos2,
                              colors="OrRd",
                              type="heatmap")
                
                i2 <- layout(i2, title = "Individuals - cos2")
                
                i2
            })
            
            output$indcontrib <- renderPlotly({
                i1 <- plot_ly(x=colnames(data()$ind$contrib),
                              y=row.names(data()$ind$contrib),
                              z=data()$ind$contrib,
                              colors="OrRd",
                              type="heatmap")
                
                i1 <- layout(i1, title = "Individuals - contribution")
                
                i1
            })
            
            output$varcos2 <- renderPlotly({
                v3 <- plot_ly(x=colnames(data()$var$cos2),
                              y=row.names(data()$var$cos2),
                              z=data()$var$cos2,
                              colors="BuGn",
                              type="heatmap")
                
                v3 <- layout(v3, title = "Variables - cos2")
                
                v3
            })
            
            output$varcontrib <- renderPlotly({
                v2 <- plot_ly(x=colnames(data()$var$contrib),
                              y=row.names(data()$var$contrib),
                              z=data()$var$contrib,
                              colors="BuGn",
                              type="heatmap")
                
                v2 <- layout(v2, title = "Variables - contribution")
                
                v2
            })
            
            output$varcor <- renderPlotly({
                v1 <- plot_ly(x=colnames(data()$var$cor),
                              y=row.names(data()$var$cor),
                              z=data()$var$cor,
                              colors="BuGn",
                              type="heatmap")
                
                v1 <- layout(v1, title = "Variables - correlation")
                
                v1
            })
        }
    })
}

substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}

shinyApp(ui=ui, server=server)