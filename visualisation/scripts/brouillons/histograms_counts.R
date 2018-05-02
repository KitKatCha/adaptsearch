library(ggplot2)
library(reshape2)
library(gridExtra)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
data <- read.table("aatypes_freqs.csv", header=TRUE, dec=".", sep=",")
data_1 <- data[seq(2, nrow(data), 3),]
data_2 <- data[seq(3, nrow(data), 3),]

ggplot(melt(data_1, id="Species"), aes(x=Species, y=value, color=variable)) +
                   geom_bar(stat="identity", fill="white")

plot1 <- ggplot() +
  geom_bar(aes(variable, value, color=variable), melt(data_1[1,], id="Species"), stat="identity", fill="white") +
  geom_point(aes(variable, value), melt(data_2[1,], id="Species"))

plot2 <- ggplot() +
  geom_bar(aes(variable, value, color=variable), melt(data_1[2,], id="Species"), stat="identity", fill="white") +
  geom_point(aes(variable, value), melt(data_2[2,], id="Species"))

plot3 <- ggplot() +
  geom_bar(aes(variable, value, color=variable), melt(data_1[3,], id="Species"), stat="identity", fill="white") +
  geom_point(aes(variable, value), melt(data_2[3,], id="Species"))

plot4 <- ggplot() +
  geom_bar(aes(variable, value, color=variable), melt(data_1[4,], id="Species"), stat="identity", fill="white") +
  geom_point(aes(variable, value), melt(data_2[4,], id="Species"))

# Faire un grid extra avec une espèce par ligne
grid.arrange(plot1, plot2, plot3, plot4, nrow=4)
