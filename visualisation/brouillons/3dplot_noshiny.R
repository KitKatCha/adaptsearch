library(scatterplot3d)
library(reshape2)
setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_nuc")
nuc <- read.table("nuc_compositions.csv", header=TRUE, dec=".", sep=",", row.names=1)
#nuc <- nuc[,seq(1, ncol(nuc)-1)]

groups <- c()
for (i in 1:nrow(nuc)) {
  groups <- c(groups, rep(i, ncol(nuc)))
}

species <- c()
for (i in 1:(ncol(nuc)/4)) {
  species <- c(species, rep(i, (ncol(nuc)/10)))
}
species <- rep(species, nrow(nuc))

values=c()
for (i in (1:nrow(nuc))) {
  for (j in (1:ncol(nuc))) {
    values <- c(values, nuc[i,j])
  }
}

nucleotide <- c("A","T","C","G")
nucletodide <- rep(nucleotide, (ncol(nuc)*nrow(nuc))/4)


test <- data.frame(groups, species, values, nucleotide)

# scatterplot3d(test, type = "h", pch=16, main = "3D plot",
#               angle=55, lab=c(3,10), lab.z=4,
#               xlab = "Orthogroups", ylab="Species", zlab="codons frequencies",
#               x.ticklabs = c("groupe1","groupe2","groupe3"), 
#               y.ticklabs = c("Ac", "Pu", "Am", "Ap", "Pf", "Pg", "Ph", "Ps", "Pp", "Pa"),
#               grid = TRUE)

colors <- c("#999999", "#E69F00", "#56B4E9", "#000000")
colors <- colors[as.numeric(test$nucleotide)]
s3d<-scatterplot3d(test[,1:3], color=colors, 
                                pch = 19, type="h",
                                lab=c(3,10), lab.z=4,
                                xlab = "Orthogroups", ylab="Species", zlab="codons frequencies",
                                x.ticklabs = c("groupe1","groupe2","groupe3"), 
                                y.ticklabs = c("Ac", "Pu", "Am", "Ap", "Pf", "Pg", "Ph", "Ps", "Pp", "Pa"),
                                grid = TRUE)
legend(s3d$xyz.convert(7.5, 3, 4.5), legend = levels(test$nucleotide),
       col =  c("#999999", "#E69F00", "#56B4E9", "#000000"), pch = 16)

s3d <- with(test, scatterplot3d(groups, species, values, color=colors, 
                                pch = 19, type="h",
                                lab=c(3,10), lab.z=4,
                                xlab = "Orthogroups", ylab="Species", zlab="codons frequencies",
                                x.ticklabs = c("groupe1","groupe2","groupe3"), 
                                y.ticklabs = c("Ac", "Pu", "Am", "Ap", "Pf", "Pg", "Ph", "Ps", "Pp", "Pa"),
                                grid = TRUE))

legend(s3d$xyz.convert(0, 0, 0), pch = 19, yjust=0,
       legend=levels(test$nucleotide), col=seq_along(levels(test$nucleotide)) )
