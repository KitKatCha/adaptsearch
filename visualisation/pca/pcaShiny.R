library(shiny)
library(shinythemes)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(ggcorrplot)

ui <- navbarPage(theme=shinytheme("sandstone"), "PCA analysis",
                 tabPanel("Data for PCA",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("data", "choose a data file", multiple=FALSE, accept=NULL)
                            ),
                            mainPanel(
                              dataTableOutput("initial_data")
                            )
                          )
                 ),
                 
                 tabPanel("PCA plot",
                          sidebarLayout(
                            sidebarPanel(
                              numericInput(inputId="axisX", label="Axis X", value=1, min=1, max=NA, step=1),
                              numericInput(inputId="axisY", label="Axis Y", value=2, min=1, max=NA, step=1),
                              actionButton("start", "Start")
                            ),
                            mainPanel(
                              downloadButton(outputId="savedPCA", label="Download the plots"),
                              
                              plotOutput("plot_inds", height = 500, height=500,
                                         dblclick = "plot1_dblclick",
                                         brush = brushOpts(
                                           id = "plot1_brush",
                                           resetOnNew = TRUE
                                         )
                              ),
                              
                              plotOutput("plot_vars", height = 500, height=500,
                                         dblclick = "plot2_dblclick",
                                         brush = brushOpts(
                                           id = "plot2_brush",
                                           resetOnNew = TRUE
                                         )
                              ),
                              plotOutput("biplot", height = 500, height=500,
                                         dblclick = "plot3_dblclick",
                                         brush = brushOpts(
                                           id = "plot3_brush",
                                           resetOnNew = TRUE
                                         )
                              ),
                              plotOutput("axes")
                            )
                          )
                 ),    
                 
                 tabPanel("Individuals data",
                          mainPanel(
                            plotOutput("indcos2", width=1000, height=250,
                                       dblclick = "plot4_dblclick",
                                       brush = brushOpts(
                                         id = "plot4_brush",
                                         resetOnNew = TRUE
                                       )
                            ),
                            plotOutput("indcontrib", width=1000, height=250,
                                       dblclick = "plot5_dblclick",
                                       brush = brushOpts(
                                         id = "plot5_brush",
                                         resetOnNew = TRUE
                                       )
                            ),
                            dataTableOutput("ind_coord"),
                            dataTableOutput("ind_cos2"),
                            dataTableOutput("ind_contrib")
                          )
                 ),
                 
                 tabPanel("Variables data",
                          mainPanel(
                            downloadButton(outputId="savedCor", label="Download the plots"),
                            plotOutput("varcor", width=1000, height=250,
                                       dblclick = "plot6_dblclick",
                                       brush = brushOpts(
                                         id = "plot6_brush",
                                         resetOnNew = TRUE
                                       )
                            ),
                            plotOutput("varcos2", width=1000, height=250,
                                       dblclick = "plot7_dblclick",
                                       brush = brushOpts(
                                         id = "plot7_brush",
                                         resetOnNew = TRUE
                                       )
                            ),
                            plotOutput("varcontrib", width=1000, height=250,
                                       dblclick = "plot8_dblclick",
                                       brush = brushOpts(
                                         id = "plot8_brush",
                                         resetOnNew = TRUE
                                       )
                            ),
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
    dat <-read.table(input$data$datapath, header=TRUE, dec=".", sep=",", row.names=1)
    if (input$data$name != "gc_counts.csv") {
      toDelete <- seq(1, nrow(dat), 2)
      dat_good <- dat[toDelete,]
    }
  })
  
  # data_table
  output$initial_data <- renderDataTable(
    data(), options=list(pageLength=10)
  )
  
  observe({
    if(input$start != 0) {
      # pca object
      pca <- reactive({
        if (is.null(data())) {return(NULL)}
        res.pca = PCA(data(), scale.unit=TRUE, graph=F)
      })
      
      output$savedPCA <- downloadHandler(
        filename = function() {
          paste("PCA_Plots", "pdf", sep=".")
        },
        # do not save zoomed plot
        content = function(file) {
          pdf(file) # Create a new pdf device
          print(fviz_pca_ind(pca(), col.ind="cos2", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE))
          print(fviz_pca_var(pca(), col.var="contrib", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel=TRUE))
          print(fviz_pca_biplot(pca(), repel=TRUE))
          print(fviz_screeplot(pca(), addlabels = TRUE))
          dev.off() # Close the pdf device
        }
      )
      
      #https://shiny.rstudio.com/articles/selecting-rows-of-data.html
      isolate({
        ranges1 <- reactiveValues(x1 = NULL, y1 = NULL)
        # plot_ind
        output$plot_inds <- renderPlot({
          if (is.null(data())) {return(NULL)}
            fviz_pca_ind(pca(), axes=c(input$axisX, input$axisY), 
                         col.ind="cos2", 
                         gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), 
                         repel=TRUE) +
              coord_cartesian(xlim = ranges1$x1, ylim = ranges1$y1, expand = FALSE)
        })
        
        # When a double-click happens, check if there's a brush on the plot.
        # If so, zoom to the brush bounds; if not, reset the zoom.
        observeEvent(input$plot1_dblclick, {
          brush <- input$plot1_brush
          if (!is.null(brush)) {
            ranges1$x1 <- c(brush$xmin, brush$xmax)
            ranges1$y1 <- c(brush$ymin, brush$ymax)
            
          } else {
            ranges1$x1 <- NULL
            ranges1$y1 <- NULL
          }
        })
        
        # plot var
        ranges2 <- reactiveValues(x2 = NULL, y2 = NULL)
        output$plot_vars <- renderPlot({
          if (is.null(data())) {return(NULL)}
          fviz_pca_var(pca(), axes=c(input$axisX, input$axisY), 
                       col.var="contrib", 
                       gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), 
                       repel=TRUE) +
            coord_cartesian(xlim = ranges2$x2, ylim = ranges2$y2, expand = FALSE)
        })
        
        # When a double-click happens, check if there's a brush on the plot.
        # If so, zoom to the brush bounds; if not, reset the zoom.
        observeEvent(input$plot2_dblclick, {
          brush <- input$plot2_brush
          if (!is.null(brush)) {
            ranges2$x2 <- c(brush$xmin, brush$xmax)
            ranges2$y2 <- c(brush$ymin, brush$ymax)
            
          } else {
            ranges2$x2 <- NULL
            ranges2$y2 <- NULL
          }
        })
        
        # bi-plot
        ranges3 <- reactiveValues(x3 = NULL, y3 = NULL)
        output$biplot <- renderPlot({
          if (is.null(data())) {return(NULL)}
          fviz_pca_biplot(pca(), axes=c(input$axisX, input$axisY), repel=TRUE) +
            coord_cartesian(xlim = ranges3$x3, ylim = ranges3$y3, expand = FALSE)
        })
        
        # When a double-click happens, check if there's a brush on the plot.
        # If so, zoom to the brush bounds; if not, reset the zoom.
        observeEvent(input$plot3_dblclick, {
          brush <- input$plot3_brush
          if (!is.null(brush)) {
            ranges3$x3 <- c(brush$xmin, brush$xmax)
            ranges3$y3 <- c(brush$ymin, brush$ymax)
            
          } else {
            ranges3$x3 <- NULL
            ranges3$y3 <- NULL
          }
        })
      })
      
      # Axes
      output$axes <- renderPlot({
        if (is.null(data())) {return(NULL)}
        fviz_screeplot(pca(), addlabels = TRUE)
      })
      
      # Individus
      ind <- reactive({
        get_pca_ind(pca())
      })
      
      # Individuals : plots cos2 and contrib
      ranges4 <- reactiveValues(x4 = NULL, y4 = NULL)
      output$indcos2 <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggcorrplot(ind()$cos2, method="square", lab=TRUE, title="cos2 for the individuals") +
          coord_cartesian(xlim = ranges4$x4, ylim = ranges4$y4, expand = FALSE)
      })
      
      # When a double-click happens, check if there's a brush on the plot.
      # If so, zoom to the brush bounds; if not, reset the zoom.
      observeEvent(input$plot4_dblclick, {
        brush <- input$plot4_brush
        if (!is.null(brush)) {
          ranges4$x4 <- c(brush$xmin, brush$xmax)
          ranges4$y4 <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges4$x4 <- NULL
          ranges4$y4 <- NULL
        }
      })
      
      ranges5 <- reactiveValues(x5 = NULL, y5 = NULL)
      output$indcontrib <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggcorrplot(ind()$contrib, method="square", lab=TRUE, title="contributions of the individuals") +
          coord_cartesian(xlim = ranges5$x5, ylim = ranges5$y5, expand = FALSE)
      })
      
      observeEvent(input$plot5_dblclick, {
        brush <- input$plot5_brush
        if (!is.null(brush)) {
          ranges5$x5 <- c(brush$xmin, brush$xmax)
          ranges5$y5 <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges5$x5 <- NULL
          ranges5$y5 <- NULL
        }
      })
      
      # Individuals : raw tables cos2 and contrib
      output$ind_coord <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        ind()$coord
      }, options=list(pageLength=10))    
      output$ind_cos2 <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        ind()$cos2
      })    
      output$ind_contrib <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        ind()$contrib
      })
      
      # Variables
      var <- reactive({
        get_pca_var(pca())
      })
      
      output$savedCor <- downloadHandler(
        filename = function() {
          paste("variables_Data", "pdf", sep=".")
        },
        content = function(file) {
          pdf(file) # Create a new pdf device
          plot1 <- ggcorrplot(var()$cor, method="square", lab=TRUE, lab_size=2, title="correlations variables - dimensions")
          plot2 <- ggcorrplot(var()$cos2, method="square", lab=TRUE, lab_size=2, title="cos2 for the variables")
          plot3 <- ggcorrplot(var()$contrib, method="square", lab=TRUE, lab_size=2, title="contributions of the variables")
          print(plot1)
          print(plot2)
          print(plot3)
          dev.off() # Close the pdf device
        }
      )
      
      # plots variables cor, contrib, cos2
      ranges6 <- reactiveValues(x = NULL, y = NULL)
      output$varcor <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggcorrplot(var()$cor, method="square", lab=TRUE, title="correlations variables - dimensions") +
          coord_cartesian(xlim = ranges6$x, ylim = ranges6$y, expand = FALSE)
      })
      
      # When a double-click happens, check if there's a brush on the plot.
      # If so, zoom to the brush bounds; if not, reset the zoom.
      observeEvent(input$plot6_dblclick, {
        brush <- input$plot6_brush
        if (!is.null(brush)) {
          ranges6$x <- c(brush$xmin, brush$xmax)
          ranges6$y <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges6$x <- NULL
          ranges6$y <- NULL
        }
      })
      
      ranges7 <- reactiveValues(x = NULL, y = NULL)
      output$varcos2 <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggcorrplot(var()$cos2, method="square", lab=TRUE, title="cos2 for the variables") +
          coord_cartesian(xlim = ranges7$x, ylim = ranges7$y, expand = FALSE)
      })
      
      # When a double-click happens, check if there's a brush on the plot.
      # If so, zoom to the brush bounds; if not, reset the zoom.
      observeEvent(input$plot7_dblclick, {
        brush <- input$plot7_brush
        if (!is.null(brush)) {
          ranges7$x <- c(brush$xmin, brush$xmax)
          ranges7$y <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges7$x <- NULL
          ranges7$y <- NULL
        }
      })
      
      ranges8 <- reactiveValues(x = NULL, y = NULL)
      output$varcontrib <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggcorrplot(var()$contrib, method="square", lab=TRUE, title="contributions of the variables") +
          coord_cartesian(xlim = ranges8$x, ylim = ranges5$y, expand = FALSE)
      })
      
      # When a double-click happens, check if there's a brush on the plot.
      # If so, zoom to the brush bounds; if not, reset the zoom.
      observeEvent(input$plot8_dblclick, {
        brush <- input$plot8_brush
        if (!is.null(brush)) {
          ranges8$x <- c(brush$xmin, brush$xmax)
          ranges8$y <- c(brush$ymin, brush$ymax)
          
        } else {
          ranges8$x <- NULL
          ranges8$y <- NULL
        }
      })
      
      # raw tables : variables cor, contrib, cos2
      output$var_coord <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        var()$coord
      }, options=list(pageLength=10))

      output$var_cor <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        var()$cor
      }, options=list(pageLength=10))

      output$var_cos2 <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        var()$cos2
      }, options=list(pageLength=10))

      output$var_contrib <- renderDataTable({
        if (is.null(data())) {return(NULL)}
        var()$contrib
      }, options=list(pageLength=10))
    }
  })
}

shinyApp(ui=ui, server=server)
