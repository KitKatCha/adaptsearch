#load packages
library(shiny)
library(shinyWidgets)
library(plotly)

ui <- bootsrapPage(


    # Parameters
    selectInput(inputId="what", label="Choose what to plot",
        c("Countings" = "counts", "Frequencies"="freqs")),
    numericInput(inputId="axisX", label="Axis X", value=1, min=1, max=NA, step=1),
    numericInput(inputId="axisY", label="Axis Y", value=2, min=1, max=NA, step=1),
    selectInput(inputId="what_inds", label="Choose parameter for individuals color range", 
        c("cos2" = "cos2", "Contribution" = "contrib")),
    selectInput(inputId="what_vars", label="Choose parameter for variables color range", 
        c("Correaltion"="cor", "cos2" = "cos2", "Contribution" = "contrib")),
    numericInput(inputId="ind_colr", label="Choose axis for individuals color range", value=1, min=1, max=NA, step=1),
    numericInput(inputId="var_colr", label="Choose axis for variables color range", value=2, min=1, max=NA, step=1),
    switchInput(
        inputId = "plot_details",
        label = strong("Show indivuals and variables details"),
        value = FALSE
    ),

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
    )
    conditionalPanel(
        condition = "plot_details",
        fluidRow(
            colum(
                6,
                div(style = "border-style: solid; border-width: thin; border-color: #000000",
                    plotlyOutput('indcos2')
                )
            ),
            colum(
                6,
                div(style = "border-style: solid; border-width: thin; border-color: #000000",
                    plotlyOutput('indcontrib')
                )
            )            
        )
        fluidRow(
            colum(
                4,
                div(style = "border-style: solid; border-width: thin; border-color: #000000",
                    plotlyOutput('varcos2')
                )
            ),
            colum(
                4,
                div(style = "border-style: solid; border-width: thin; border-color: #000000",
                    plotlyOutput('varcontrib')
                )
            ),
            colum(
                4,
                div(style = "border-style: solid; border-width: thin; border-color: #000000",
                    plotlyOutput('varcor')
                )
            )
        )
    )
)

server <- function (input, output) {

    data <- read.table("aa_freqs.csv", header=TRUE, dec=".", sep=",", row.names=1)
    counts <- data[seq(1, nrow(data), 3),]
    freqs <- data[seq(2, nrow(data), 3),]
    row.names(freqs) <- substrLeft(row.names(freqs),2)
    row.names(counts) <- substrLeft(row.names(counts),2)
    res.pca = PCA(freqs, scale.unit=TRUE, graph=F, axes=c(1,2))
    ind <- as.data.frame(res.pca$ind$coord)
    cos2 <- as.data.frame(res.pca$ind$cos2)
    var <- as.data.frame(res.pca$var$coord)

    output$PCA_ind <- renderPlotly({
        p <- plot_ly(ind, 
             x=ind[,input$axisX], 
             y=ind[,input$axisY],
             type = 'scatter',
             text=rownames(input$what),
             textposition='top',
             mode="markers+text", 
             color=input$what_inds[,input$ind_colr],
             colors="OrRd",
             marker=list(size=11))

        p <- layout(p, title = "PCA on individuals", 
                    xaxis = list(title = input$axisX),
                    yaxis = list(title = input$axisY))
        p
    })

    output$PCA_var <- renderPlotly({
        p2 <- plot_ly(var, 
             x=var[,input$axisX], 
             y=var[,input$axisY],
             type = 'scatter',
             text=colnames(input$what),
             textposition='top',
             mode="markers+text", 
             color=input$what_vars[,input$var_colr],
             colors="OrRd",
             marker=list(size=11))

        p2 <- layout(p2, title = "PCA on variables", 
                    xaxis = list(title = input$axisX),
                    yaxis = list(title = input$axisY))
        p2
    })

    Output$PCA_biplot <- renderPlotly({

    })

    Output$PCA_eigen <- renderPlotly({
        eigen <- plot_ly(res.pca$eig) %>%
          add_trace(x=row.names(res.pca$eig),
                    y=res.pca$eig[,2],
                    type='bar') %>%
          add_trace(res.pca$eig,
                    x=row.names(res.pca$eig),
                    y=res.pca$eig[,3],
                    type='scatter',
                    mode='lines') %>%
          layout(yaxis=list(type='log'))

eigen
    })

}

substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}