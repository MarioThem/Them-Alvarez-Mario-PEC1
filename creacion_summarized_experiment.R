# Se importan las librerias correspondientes
library(SummarizedExperiment)
library(metabolomicsWorkbenchR)
# Se descarga el mwtab del estudio seleccionado
mwtab<-do_query(context = "study", input_item = "study_id", input_value = "ST003674", 
                output_item = "mwtab")
# Se descargan los datos experimentales del estudio
datos<-as.data.frame(do_query(context = "study", input_item = "study_id", input_value = "ST003674", 
                output_item = "data"))
# Se descargan los datos de las muestras
datos_muestras<-as.data.frame(do_query(context = "study", input_item = "study_id", input_value = "ST003674", 
                              output_item = "factors"))
# Se descargan los datos de los metabolitos
datos_metabolitos<-as.data.frame(do_query(context = "study", input_item = "study_id", input_value = "ST003674", 
                                          output_item = "metabolites"))
# Se construye el slot assay
assay<-datos[ ,8:52]
for(i in seq(1,45)){
  colnames<- c(colnames, mwtab$SUBJECT_SAMPLE_FACTORS[[i]][[2]])
}
colnames(assay)<-colnames[2:46]
for(i in seq(1,139)){
  rownames<-c(rownames,
              mwtab[["NMR_METABOLITE_DATA"]][["Metabolites"]][[i]][["Metabolite"]])
}
rownames(assay)<-rownames[2:140]
assay <- assay[order(rownames(assay)),]
# Se construye el slot colData
colData <- datos_muestras[ ,c(2:4,8:9)]
rownames(colData)<-colData$ST003674.local_sample_id
colData <- colData[order(colData$ST003674.local_sample_id),]
# Se construye el slot rowData
rowData<-data.frame(datos_metabolitos[ ,1:5])
rownames(rowData)<-datos_metabolitos$metabolite_name
rowData <- rowData[order(rownames(rowData)),]
# Se construye el objeto de tipo SummarizedExperiment
st003674<-SummarizedExperiment(
  assay = assay, 
  rowData = rowData,
  colData = colData,
  metadata = mwtab$STUDY)
# Se muestran los distintos slots del objeto
head(assay(st003674))
head(rowData(st003674))
head(colData(st003674))
head(metadata(st003674))
save(st003674, file = "ST003674.rda")
