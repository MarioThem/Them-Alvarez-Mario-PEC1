library(ggplot2)
library(stats)
library(SummarizedExperiment)
load("ST003674.rda")
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

# Boxplot para ver los valores generales de los metabolitos
boxplot(na.omit(st003570_normalizado), 
        main="Valores de los metabolitos normalizados",
        col=c(rep("red",3), rep("blue",3)),
        ylab="Peak values",
        las=2
        )
legend(x="topleft", 
       legend = c("CaSki CM", "Hepatocyte CM"), 
       fill = c("red", "blue"))
# PCA para determinar posibles grupos 
pc <- prcomp(t(na.omit(st003570_normalizado)), scale=FALSE)
loads <- round(pc$sdev^2/sum(pc$sdev)*100,1)
xlab<-paste("PC1", loads[1], "%")
ylab<-paste("PC2", loads[2], "%")
plot(pc)
plot(pc$x, 
     xlab=xlab, 
     ylab=ylab,
     col=c(rep("red",3), rep("blue",3))
     )
text(pc$x, row.names(pc$x), cex = 0.6, pos = 3)
print(pc)
