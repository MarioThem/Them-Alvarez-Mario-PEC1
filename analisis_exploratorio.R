library(ggplot2)
library(stats)
library(SummarizedExperiment)
library(factoextra)
load("ST003674.rda")
st003674_assay<-assay(st003674)
# Se creaN nuevos SE con los valores de cada metabolito para cada una de las condiciones experimentales
hombres_INT<-st003674[,st003674$ST003674.Pneumonia=="INT" & 
                        st003674$ST003674.Sex=="M "]
hombres_LOB<-st003674[,st003674$ST003674.Pneumonia=="LOB" & 
                        st003674$ST003674.Sex=="M "]
mujeres_INT<-st003674[,st003674$ST003674.Pneumonia=="INT" & 
                        st003674$ST003674.Sex=="F "]
mujeres_LOB<-st003674[,st003674$ST003674.Pneumonia=="LOB" & 
                        st003674$ST003674.Sex=="F "]
# Se calcula la media para cada metabolito en cada grupo experimental
medias_1 <- rowMeans(assay(hombres_INT), na.rm = TRUE)
medias_2 <- rowMeans(assay(hombres_LOB), na.rm = TRUE)
medias_3 <- rowMeans(assay(mujeres_INT), na.rm = TRUE)
medias_4 <- rowMeans(assay(mujeres_LOB), na.rm = TRUE)

# Se crea un data frame con las medias por metabolito y grupo experimental
medias <- data.frame(medias_1,medias_2,medias_3,medias_4)
colnames(medias)<- c("hombres-INT", "hombres-LOB", "mujeres_INT", "mujeres-LOB")
# Se estandariza dividiendo por la media de cada grupo experimental
medias$`hombres-INT`<-medias$`hombres-INT`/mean(medias$`hombres-INT`)
medias$`hombres-LOB`<-medias$`hombres-LOB`/mean(medias$`hombres-LOB`)
medias$mujeres_INT<-medias$mujeres_INT/mean(medias$mujeres_INT)
medias$`mujeres-LOB`<-medias$`mujeres-LOB`/mean(medias$`mujeres-LOB`)
# Boxplot para ver los valores generales de los metabolitos
boxplot(medias)
# PCA para determinar posibles grupos y outliers
pc <- prcomp(t(na.omit(st003674_assay)), scale. = TRUE)
# Contribuciones a la variabilidad de las componentes del PCA(screeplot)
fviz_eig(pc, addlabels = TRUE)
# Visualización de las 50 variables que más contribuyen en el círculo de correlación de las 2 primeras PC
fviz_pca_var(pc, col.var = "contrib",
             repel = TRUE, select.var = list(contrib=50))
# Visualización de la contribución de los individuos a los 2 primeros PC
grupos <- factor(paste(st003674$ST003674.Sex,st003674$ST003674.Pneumonia))
grupos_neumonia <- st003674$ST003674.Pneumonia
fviz_pca_ind(pc, repel = TRUE, col.ind = grupos,
             addEllipses = TRUE, ellipse.type = "confidence")
fviz_pca_ind(pc, repel = TRUE, col.ind = grupos_neumonia,
             addEllipses = TRUE, ellipse.type = "confidence")
# Lo mismo pero en bar plot
fviz_contrib(pc, choice = "var", top = 10)
fviz_contrib(pc, choice = "var", top = 10, axes = 2)
fviz_contrib(pc, choice = "var", top = 10, axes = 3)
fviz_contrib(pc, choice = "var", top = 10, axes = 4)
fviz_contrib(pc, choice = "ind", top = 10, axes = 1)
# Cluster analysis
dist <-dist(medias, method = "euclidean")
hc <- hclust(dist, method = "ward.D2")
fviz_dend(hc)
res.coph <- cophenetic(hc)
cor(dist, res.coph)
grp <- cutree(hc, k=4)
table(grp)
rownames(st003674_assay)[grp==3|grp==4|grp==1]
fviz_dend(hc,k=4, repel = TRUE)
