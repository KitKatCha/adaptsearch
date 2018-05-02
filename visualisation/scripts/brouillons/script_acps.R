# Embryon de script
# Les fichiers de données ont été crées manuellement, il faut automatiser

library(FactoMineR)
library(factoextra)
# Fréquences codons
codons <- read.table("pca_codons", dec=".", header=TRUE, sep=" ", row.names=1)
res.pca_codons <- PCA(codons, scale.unit=TRUE, axes=c(1,2), graph=F)
#plot(res.pca_codons, choix="ind")
#plot(res.pca_codons, choix="var", main="Codons frequences - Variables")
fviz_pca_var(res.pca_codons, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
fviz_pca_ind(res.pca_codons, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca_codons, repel = TRUE)

# Fréquences AA
aa <- read.table("pca_aa", dec=".", header=TRUE, sep=" ", row.names=1)
res.pca_aa = PCA(aa, scale.unit=TRUE, graph=T)
fviz_pca_var(res.pca_aa, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
fviz_pca_ind(res.pca_aa, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca_aa, repel = TRUE)

# Fréquences types AA
aatypes <- read.table("pca_aatypes", dec=".", header=TRUE, sep=" ", row.names=1)
res.pca_aatypes = PCA(aatypes, scale.unit=TRUE, graph=T)
fviz_pca_var(res.pca_aatypes, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
fviz_pca_ind(res.pca_aatypes, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca_aatypes, repel = TRUE)

# GC and various stuff
various <- read.table("pca_various", dec=".", header=TRUE, sep=" ", row.names=1)
res.pca_various = PCA(various, scale.unit=TRUE, graph=T)
fviz_pca_var(res.pca_various, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
fviz_pca_ind(res.pca_various, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca_various, repel = TRUE)

###########

