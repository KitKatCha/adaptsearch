library(plotly)

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

res.pca = PCA(freqs, scale.unit=FALSE, graph=F, axes=c(1,2))
ind <- as.data.frame(res.pca$ind$coord)
cos2 <- as.data.frame(res.pca$ind$cos2)
p <- plot_ly(ind, 
             x=ind[,2], 
             y=ind[,3],
             type = 'scatter',
             text=rownames(freqs),
             mode="markers", 
             color=cos2[,3], 
             marker=list(size=11))

p <- layout(p, title = "PCA on individuals", 
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"))
p

m <- plot_ly(x=colnames(res.pca$ind$contrib),
             y=row.names(res.pca$ind$contrib),
             z=res.pca$ind$contrib,
             type="heatmap")

####################

# cor = TRUE indicates that PCA is performed on 
# standardized data (mean = 0, variance = 1)
pcaCars <- princomp(mtcars, cor = TRUE)
names(pcaCars) # view objects stored in pcaCars
summary(pcaCars) # proportion of variance explained
plot(pcaCars, type = "l") # scree plot
carsHC <- hclust(dist(pcaCars$scores), method = "ward.D2") # cluster cars
plot(carsHC) # dendrogram
carsClusters <- cutree(carsHC, k = 3) # cut the dendrogram into 3 clusters
# add cluster to data frame of scores
carsDf <- data.frame(pcaCars$scores, "cluster" = factor(carsClusters))
carsDf <- transform(carsDf, cluster_name = paste("Cluster",carsClusters))

p <- plot_ly(carsDf, x=Comp.1, y=Comp.2, text=rownames(carsDf),
             mode="markers", color=cluster_name, marker=list(size=11)) 

p <- layout(p, title = "PCA Clusters from Hierarchical Clustering of Cars Data", 
            xaxis = list(title = "PC 1"),
            yaxis = list(title = "PC 2"))

p
