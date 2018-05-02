library(shiny)
library(shinyWidgets)
library(plotly)
library(ggplot2)
library(ggdendro)

ui <- fluidPage(

    mainPanel(
        plotlyOutput('heatmapDendro')
    )
)

server <- function (input, output) {

    # Code mainly taken from https://plot.ly/ggplot2/ggdendro-dendrograms/

    counts <- data[seq(1, nrow(data), 3),]
    freqs <- data[seq(2, nrow(data), 3),]

    dd.col <- as.dendrogram(hclust(dist(freqs)))
    dd.row <- as.dendrogram(hclust(dist(t(freqs))))
    dx <- dendro_data(dd.row)
    dy <- dendro_data(dd.col)

    # x/y dendograms
    px1 <- ggplot(segment(dx)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend))      
    py2 <- ggplot(segment(dy)) + 
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
      coord_flip() + 
      scale_y_reverse(expand=c(0.2, 0))

    # heatmap
    col.ord <- order.dendrogram(dd.col)
    row.ord <- order.dendrogram(dd.row)
    xx <- scale(freqs)[col.ord, row.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$spec <- xx_names[[1]]
    df$spec <- with(df, factor(spec, levels=spec, ordered=TRUE))
    mdf <- reshape2::melt(df, id.vars="spec")
    p <- ggplot(mdf, aes(x = variable, y = spec)) + geom_tile(aes(fill = value))

    output$heatmapDendro <- renderPlotly({
        subplot(py2, p, margin=0.03, widths=c(0.25,0.75))
    })
}

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}