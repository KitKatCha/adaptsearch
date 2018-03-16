library(plotly)
library(reshape2)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
dataframe <- read.table("aa_transitions_freqs.csv", header=TRUE, dec=".", sep=",", row.names=1)

skip <- c()
for (i in seq(1,400,by=21)) {
  skip <- c(skip, i)
}

starts <- skip[1:(length(skip)-1)]+1

reordered <- c()
l <- 18

for (i in 1:length(starts)) {
  step <- 19
  substep <- 1
  for (j in seq(starts[i],(starts[i]+l))) {
    couple <- c(j, j+step)
    print(couple)
    reordered <- c(reordered, couple)
    substep <- substep +1
    step <- 19 * substep
  }
  l <- l-1
}
reordered

data <- dataframe[,reordered]

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

spec1 <- freqs[1,]
spec1p <- pvalues[1,]
m <- melt(spec1, variable="Type")
mp <- melt(spec1p, variable="Type")

mAbove <- m[seq(1,380, by=2),]
mAbove <- droplevels(mAbove)
mBelow <- m[seq(2,380, by=2),]
mBelow <- droplevels(mBelow)

bind <- cbind(mAbove,mBelow)
# exclude lines where both transi = 0
bind.filt <- bind[rowSums(bind[,c(2,4)])!=0,]
bind.filt <- droplevels(bind.filt)

# Find max y axis value
max_y <- max(mBelow[,2], mAbove[,2])

# Split in subTables, otherwise graphes are too larges
mAbove.1 <- mAbove[1:64,]
mAbove.2 <- mAbove[65:128,]
mAbove.3 <- mAbove[129:190,]
mBelow.1 <- mBelow[1:64,]
mBelow.2 <- mBelow[64:128,]
mBelow.3 <- mBelow[129:190,]


# without pvalues
pa1 <- plot_ly() %>%
  add_bars(x = mAbove.1[,1], 
           y = mAbove.1[,2], 
           name='Transitions X to Y', 
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  layout(yaxis=list(range=c(0,max_y)))

pb1 <- plot_ly() %>%
  add_bars(x = mBelow.1[,1], 
           y = mBelow.1[,2], 
           name='Transitions Y to X', 
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
         xaxis = list(side='top'))

pa2 <- plot_ly() %>%
  add_bars(x = mAbove.2[,1], 
           y = mAbove.2[,2], 
           name='Transitions X to Y',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  layout(yaxis=list(range=c(0,max_y)))

pb2 <- plot_ly() %>%
  add_bars(x = mBelow.2[,1], 
           y = mBelow.2[,2], 
           name='Transitions Y to X',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
         xaxis = list(side='top'))

pa3 <- plot_ly() %>%
  add_bars(x = mAbove.3[,1], 
           y = mAbove.3[,2], 
           name='Transitions X to Y',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  layout(yaxis=list(range=c(0,max_y)))

pb3 <- plot_ly() %>%
  add_bars(x = mBelow.3[,1], 
           y = mBelow.3[,2], 
           name='Transitions Y to X',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
         xaxis = list(side='top'))

subplot(pa1, pb1, nrows = 2, margin=0.1)
subplot(pa2, pb2, nrows = 2, margin=0.1)
subplot(pa3, pb3, nrows = 2, margin=0.1)

####################"

# Without null transis

pa4 <- plot_ly() %>%
  add_bars(x = bind.filt[,1], 
           y = bind.filt[,2], 
           name='Transitions X to Y',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 178, 1255)')) %>%
  layout(yaxis=list(range=c(0,max_y)))

pb4 <- plot_ly() %>%
  add_bars(x = bind.filt[,3], 
           y = bind.filt[,4], 
           name='Transitions Y to X',
           showlegend = FALSE,
           marker=list(color = 'rgb(102, 255, 102)')) %>%
  layout(yaxis = list(range=c(max_y,0), autorange=F, autorange = "reversed"),
         xaxis = list(side='top'))

subplot(pa4, pb4, nrows = 2, margin=0.1)

