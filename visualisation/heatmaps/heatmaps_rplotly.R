library(plotly)
library(ggplot2)
library(ggdendro)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
data <- read.table("aa_freqs.csv", header=TRUE, dec=".", sep=",", row.names=1)

counts <- data[seq(1, nrow(data), 3),]
freqs <- data[seq(2, nrow(data), 3),]

substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}

row.names(freqs) <- substrLeft(row.names(freqs),2)
row.names(counts) <- substrLeft(row.names(counts),2)

dd.col <- as.dendrogram(hclust(dist(freqs)))
dd.row <- as.dendrogram(hclust(dist(t(freqs))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments) 
py <- ggdend(dy$segments) + coord_flip()

px1 <- ggplot(segment(dx)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) 
  #+ geom_text(data=label(dx), aes(label=label, x=x, y=0))
py2 <- ggplot(segment(dy)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + coord_flip()
  #+ geom_text(data=label(dy), aes(label=label, x=x, y=0)) 

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

# hide axis ticks and grid lines
eaxis <- list(
  showticklabels = FALSE,
  showgrid = FALSE,
  zeroline = FALSE
)

p_empty <- plot_ly(filename="r-docs/dendrogram") %>%
  # note that margin applies to entire plot, so we can
  # add it here to make tick labels more readable
  layout(margin = list(l = 200),
         xaxis = eaxis,
         yaxis = eaxis)

subplot(px1, p_empty, p, py2, nrows = 2, margin = 0.02, heights = c(0.25,0.75), widths=c(0.75,0.25))
