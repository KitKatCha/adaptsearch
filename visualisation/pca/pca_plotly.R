library(plotly)
library(FactoMineR)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
data <- read.table("aatypes_freqs.csv", header=TRUE, dec=".", sep=",", row.names=1)

counts <- data[seq(1, nrow(data), 3),]
freqs <- data[seq(2, nrow(data), 3),]

substrLeft <- function(x, n){
  sapply(x, function(xx)
    substr(xx, 0, n)
  )
}

row.names(freqs) <- substrLeft(row.names(freqs),2)
row.names(counts) <- substrLeft(row.names(counts),2)

res.pca = PCA(freqs, scale.unit=TRUE, graph=F, axes=c(1,2))
ind <- as.data.frame(res.pca$ind$coord)
cos2 <- as.data.frame(res.pca$ind$cos2)
var <- as.data.frame(res.pca$var$coord)
cos2_v <- as.data.frame(res.pca$var$cos2)

biplot <- plot_ly(ind) %>%
  add_trace(x=ind[,1],
            y=ind[,2],
            type='scatter',
            text=rownames(freqs),
            textposition='top',
            mode="markers+text", 
            color=cos2[,1],
            colors="OrRd",
            marker=list(symbol=27, size=11)) %>%
  add_trace(var, 
          x=var[,1], 
          y=var[,2],
          type = 'scatter',
          text=colnames(freqs),
          textposition='top',
          mode="markers+text", 
          color=cos2_v[,1],
          colors="BuGn",
          marker=list(symbol=4, size=11))

biplot <- layout(biplot, title="PCA - bi-plot")

p <- plot_ly(ind, 
             x=ind[,1], 
             y=ind[,2],
             type = 'scatter',
             text=rownames(freqs),
             textposition='top',
             mode="markers+text", 
             color=cos2[,1],
             colors="OrRd",
             marker=list(size=11))

p <- layout(p, title = "PCA on individuals", 
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"))
p

i1 <- plot_ly(x=colnames(res.pca$ind$contrib),
             y=row.names(res.pca$ind$contrib),
             z=res.pca$ind$contrib,
             type="heatmap")

i1 <- layout(i1, title = "Individuals - contribution") 

i2 <- plot_ly(x=colnames(res.pca$ind$cos2),
              y=row.names(res.pca$ind$cos2),
              z=res.pca$ind$cos2,
              type="heatmap")

i2 <- layout(i2, title = "Individuals - cos2")

v1 <- plot_ly(x=colnames(res.pca$var$cor),
              y=row.names(res.pca$var$cor),
              z=res.pca$var$cor,
              colors="RdYlBu",
              type="heatmap")

v1 <- layout(v1, title = "Variables - correlation")

v2 <- plot_ly(x=colnames(res.pca$var$contrib),
              y=row.names(res.pca$var$contrib),
              z=res.pca$var$contrib,
              type="heatmap")

v2 <- layout(v2, title = "Variables - contribution")

v3 <- plot_ly(x=colnames(res.pca$var$cos2),
              y=row.names(res.pca$var$cos2),
              z=res.pca$var$cos2,
              type="heatmap")

v3 <- layout(v3, title = "Variables - cos2")

# Ajouter un mode coloration selon aa types, ou bien un mode ave différents marqueur

p2 <- plot_ly(var, 
             x=var[,1], 
             y=var[,2],
             type = 'scatter',
             text=colnames(freqs),
             textposition='top',
             mode="markers+text", 
             color=cos2_v[,1],
             colors="BuGn",
             marker=list(size=11))

p2 <- layout(p2, title = "PCA on variables", 
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"))
p2

plot_ly(var, x=var[,1], y=var[,2], type='scatter', mode='lines')

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