# Se importan las librerias correspondientes
library(ggplot2)
library(stats)
library(SummarizedExperiment)
library(factoextra)
library(SummarizedExperiment)
library(metabolomicsWorkbenchR)
library(car)
library(knitr)
load("ST003674.rda")
# Se inspeccionan las dimensiones del slot assay
dim(assay(st003674))
# Se muestran los nombres de las columnas y filas del assay
head(rownames(assay(st003674)))
head(colnames(assay(st003674)))
# Se muestra un resumen numérico de los primeros 5 registros
summary(st003674_assay[1:5])
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
medias$`hombres-INT`<-(medias$`hombres-INT`-mean(medias$`hombres-INT`))/sd(medias$`hombres-INT`)
medias$`hombres-LOB`<-(medias$`hombres-LOB`-mean(medias$`hombres-LOB`))/sd(medias$`hombres-LOB`)
medias$mujeres_INT<-(medias$mujeres_INT-mean(medias$mujeres_INT))/sd(medias$mujeres_INT)
medias$`mujeres-LOB`<-(medias$`mujeres-LOB`-mean(medias$`mujeres-LOB`))/sd(medias$`mujeres-LOB`)
# Boxplot para ver los valores generales de los metabolitos
Boxplot(medias, 
        main="Boxplot de las medias de los metabolitos estandarizadas por grupo", 
        id = list(cex = 0.4, n = 4, location = "l"))
# Estandarización de los datos 
assay<-na.omit(assay(st003674))
estandarizado <- sapply(assay, function(assay) (assay-mean(assay)/sd(assay)))
rownames(estandarizado)<-rownames(assay)
# PCA para determinar posibles grupos y outliers
pc <- prcomp(t(estandarizado), scale. = FALSE)
# Contribuciones a la variabilidad de las componentes del PCA(screeplot)
fviz_eig(pc, addlabels = TRUE, main = "Scree plot de los componentes del PCA")
# Visualización de la contribución de los individuos a los 2 primeros PC
grupos_sexo <- st003674$ST003674.Sex
grupos_neumonia <- st003674$ST003674.Pneumonia
fviz_pca_ind(pc, repel = TRUE, col.ind = grupos, label = "none",
             addEllipses = TRUE, ellipse.type = "confidence",
             title="Contribución de los individuos según el sexo")
fviz_pca_ind(pc, repel = TRUE, col.ind = grupos_neumonia, label="none",
             addEllipses = TRUE, ellipse.type = "confidence",
             title="Contribución de los individuos según el tipo de neumonia")
# Grafico de contribución de las variables a las dos primeras PCs
fviz_pca_var(pc)

