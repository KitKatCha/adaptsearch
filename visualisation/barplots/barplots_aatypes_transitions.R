library(plotly)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
dataframe <- read.table("aatypes_transitions_freqs.csv", header=TRUE, dec=".", sep=",", row.names=1)

data <- dataframe[,c(2,5,3,9,4,13,7,10,12,15,8,14)] # 1,6,11,16 : X -> X
d2 <- dataframe[,c(1,6,11,16)]

# Split countings, frequencies and pvalues
counts <- data[seq(1, nrow(data), 3),]
freqs <- data[seq(2, nrow(data), 3),]
pvalues <- data[seq(3, nrow(data), 3),]

# Simplify row names
substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}

row.names(counts) <- substrLeft(row.names(counts),5)
row.names(freqs) <- substrLeft(row.names(freqs),5)
row.names(pvalues) <- substrLeft(row.names(pvalues),5)

# Split lines according to their plotting orientation
ab <- c(2,3,4,7,8,12)
bl <- c(5,9,13,10,14,15)
rm <- c(1,6,11,16)

spec1 <- freqs[1,]
spec1p <- pvalues[1,]
m <- melt(spec1, variable="Type")
mp <- melt(spec1p, variable="Type")

# Split lines for freqs
mAbove <- m[seq(1,12, by=2),]
mAbove <- droplevels(mAbove)
mBelow <- m[seq(2,12, by=2),]
mBelow <- droplevels(mBelow)

# Split lines for pvalues
mpAbove <- mp[seq(1,12, by=2),]
mpAbove <- droplevels(mpAbove)
mpBelow <- mp[seq(2,12, by=2),]
mpBelow <- droplevels(mpBelow)

# Find max y axis value
max_y <- max(mBelow[,2], mAbove[,2])

#####################################################################

# Try with gplot and ggplotly

g1 <- ggplot(mAbove, aes(x=Type, y=value)) + 
  theme_hc()+ scale_colour_hc() +
  geom_bar(stat="identity", position="identity", color="tan1", fill="tan1") +
  theme(axis.text.x=element_text(angle=20, vjust=0.5), axis.title.x=element_blank()) +
  ylim(0, max_y)

g2 <- ggplot(mBelow, aes(x=Type, y=value)) +
  theme_hc()+ scale_colour_hc() +
  geom_bar(stat="identity", position="identity", color="steelblue1", fill="steelblue1") +
  theme(axis.text.x=element_text(angle=20, vjust=-2), axis.title.x=element_blank()) +
  scale_x_discrete(position = "top") +
  scale_y_reverse(limits=c(max_y, 0))

pg1 <- ggplotly(g1)
pg2 <- ggplotly(g2)
subplot(pg1,pg2, nrows=2, margin=0.075)

grid.arrange(g1, g2, nrow=2)

###############################################

# Only Plotly, without pvalues

pa <- plot_ly() %>%
  add_bars(x = mAbove[,1], 
           y = mAbove[,2], 
           name='Transitions X to Y', 
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  layout(yaxis=list(range=c(0,max_y)))

pb <- plot_ly() %>%
  add_bars(x = mBelow[,1], 
           y = mBelow[,2], 
           name='Transitions Y to X', 
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
         xaxis = list(side='top'))

subplot(pa, pb, nrows=2, margin=0.1)

###############################################

# Try with plotly

ay <- list(
  tickfont = list(color = "red"),
  overlaying = "y",
  side = "right",
  title = "pvalues",
  range=c(0, 1)
)

p1 <- plot_ly() %>%
  add_bars(x = mAbove[,1], 
           y = mAbove[,2], 
           name='Transitions X to Y',
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  add_markers(x=mpAbove[,1], 
              y=mpAbove[,2], 
              name='pvalues', 
              yaxis='y2', 
              marker=list(color = 'red')) %>%
  layout(yaxis=list(range=c(0,max_y)), 
         yaxis2=ay,
         title=row.names(freqs[1,]))

ay2 <- list(
  tickfont = list(color = "red"),
  overlaying = "y3",
  side = "right",
  title = "pvalues",
  range=c(1, 0),
  autorange=F,
  autorange = "reversed"
)

p2 <- plot_ly() %>%
  add_bars(x = mBelow[,1], 
           y = mBelow[,2], 
           name='Transitions Y to X', 
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  add_markers(x=mpBelow[,1], 
              y=mpBelow[,2], 
              name='pvalues',
              showlegend = FALSE,
              yaxis='y2', 
              marker=list(color = 'red')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"), 
         yaxis2 = ay2,
         xaxis = list(side='top'))

subplot(p1, p2, nrows=2, margin=0.1, titleX=T)
