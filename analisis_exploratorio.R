library(MSPrep)
library(ggplot2)
library(stats)
load("ST003570.rda")
datos <- assay(st003570)
# Normalizamos en logaritmo de base 10 los valores de los compuestos                      
st003570_normalizado <- log10(datos) 
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
plot(pc$x, 
     xlab=xlab, 
     ylab=ylab,
     col=c(rep("red",3), rep("blue",3))
     )
text(pc$x, row.names(pc$x), cex = 0.6, pos = 3)

