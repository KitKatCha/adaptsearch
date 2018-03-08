library(ggplot2)
library(easyGgplot2)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(ggcorrplot)
library(gridExtra)
dat <- read.table("aatypes_transitions.csv", header=TRUE, dec=".", sep=",", row.names=1)
toDelete <- seq(1, nrow(dat), 2)
dat_good <- dat[toDelete,]
res.pca = PCA(dat_good, scale.unit=TRUE, graph=F, axes=c(1,2))
par(mfrow=c(1,3))
corrplot(get_pca_var(res.pca)$cor, method="number", is.corr=FALSE, title="correlations variables - dimensions")
corrplot(get_pca_var(res.pca)$cos2, method="number", is.corr=FALSE, title="cos2 for the variables")
corrplot(get_pca_var(res.pca)$contrib, method="number", is.corr=FALSE, title="contributions of the variables")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)

plot1 <- ggcorrplot(get_pca_var(res.pca)$cor, method="square", lab = TRUE, lab_size=3, 
                    title="correlations variables - dimensions")
plot2 <- ggcorrplot(get_pca_var(res.pca)$cos2, method="square", lab = TRUE, lab_size=3, 
                    title="cos2 for the variables")
plot3 <- ggcorrplot(get_pca_var(res.pca)$contrib, method="square", lab = TRUE, lab_size=3, 
                    title="contributions of the variables")

grid.arrange(plot1, plot2, plot3, nrow=3)

ggplot2.multiplot(plot1,plot2,plot3, cols=2)

pdf("toto.pdf")
print(plot1)
print(plot2)
print(plot3)
dev.off()
