library(shiny)
library(shinythemes)
library(ggplot2)
library(FactoMineR)
library(factoextra)

ui <- navbarPage(theme=shinytheme("sandstone"), "Pretty PCA analysis",
    tabPanel("Data for PCA",
        sidebarLayout(
            sidebarPanel(
                fileInput("data", "choose a data file", multiple=FALSE, accept=NULL),
                selectInput(inputId="axisX", label="Axis X", value=1, min=1, max=10),
                selectInput(inputId="axisY", label="Axis Y", value=2, min=1, max=10),
                radioButton(inputId="want_biplot", label="Bi-plot ?", choices=c("No","Yes"))
                actionButton("startpcs", "Start")
            ),
            mainPanel(
                dataTableOutput("initial_data")
            )
        )
    ),

    tabPanel("PCA plot",
        mainPanel(
            plotOutput("plot_inds"),
            plotOutput("plot_vars"),
            conditionalPanel(
                condition="input.want_biplot == Yes",
                plotOutput("biplot")
            )           
        )
    ),

    tabPanel("Eigen",
        mainpanel(
            plotOutput("axes"),

        )
    )

    tabPanel("Individuals data",
        mainpanel(            
            dataTableOutput("ind_coord"),
            dataTableOutput("ind_cor"),     
            dataTableOutput("ind_cos2"),
            dataTableOutput("ind_contrib")
        )
    )

    tabPanel("Variables data",
        mainpanel(
            dataTableOutput("var_coord"),
            dataTableOutput("var_cor"),     
            dataTableOutput("var_cos2"),
            dataTableOutput("var_contrib")
        )
    )
)

server <- function(input, output, session) {
    # Data
    data <- reactive({
        if (is.null(input$data)) {return(NULL)}
        read.table(input$data$datapath, header=TRUE, dec=".", sep=" ", row.names=1)
        res.pca = PCA(data, scale.unit=TRUE, graph=F, axes=c(input$axisX,input$axisY))
    })

    # plot_ind
    output$plot_inds <- renderPlot({
        fviz_pca_ind(res.pca, col.ind="cos2", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE)
    })

    # plot var
    output$plot_vars <- renderPlot({
        fviz_pca_var(res.pca, col.var="contrib", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE)
    })

    # bi-plot
    output$biplot <- renderPlot({
        fviz_pca_biplot(res.pca, repel=TRUE)
    })

    # Axes
    output$axes <- renderPlot({
        barplot(res.pca$eig[,2], xlab="Axes",ylab="Percentage of variance")
    })

    # Individus
    output$ind_coord <- renderDataTable({
        res.pca$ind$coord
    })    
    output$ind_cor <- renderDataTable({
        res.pca$ind$cor
    })    
    output$ind_cos2 <- renderDataTable({
        res.pca$ind$cos2
    })    
    output$ind_contrib <- renderDataTable({
        res.pca$ind$contrib
    })

    # Variables
    output$var_coord <- renderDataTable({
        res.pca$var$coord
    })    
    output$ind_cor <- renderDataTable({
        res.pca$var$cor
    })    
    output$ind_cos2 <- renderDataTable({
        res.pca$var$cos2
    })    
    output$ind_contrib <- renderDataTable({
        res.pca$var$contrib
    })

}

shinyApp(ui=ui, server=server)