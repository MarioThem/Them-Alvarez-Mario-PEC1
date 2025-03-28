# Se importan las librerias correspondientes
library(SummarizedExperiment)
library(metabolomicsWorkbenchR)
# Se descarga el mwtab del estudio seleccionado
mwtab<-do_query(context = "study", input_item = "study_id", input_value = "ST003570", 
                output_item = "mwtab")
# Se descargan los datos del estudio
datos<-as.data.frame(do_query(context = "study", input_item = "study_id", input_value = "ST003570", 
                output_item = "data"))
# Se descargan los datos de las muestras
datos_muestras<-as.data.frame(do_query(context = "study", input_item = "study_id", input_value = "ST003570", 
                              output_item = "factors"))
# Se construye el slot assay
assay<-datos[ ,8:13]
for(i in seq(1,6)){
  colnames<- c(colnames, mwtab$SUBJECT_SAMPLE_FACTORS[[i]][[2]])
}
colnames(assay)<-colnames[2:7]
for(i in seq(1,18)){
  rownames<-c(rownames,
              mwtab[["MS_METABOLITE_DATA"]][["Metabolites"]][[i]][["KEGG identifiers"]])
}
rownames(assay)<-rownames[2:19]
# Se construye el slot colData
colData <- datos_muestras[ ,2:8]
# Se construye el slot rowData
rowData<-
  
# Se construye el objeto de tipo SummarizedExperiment
st003570<-SummarizedExperiment(
  assay = assay, 
  rowData = rowData,
  colData = colData, 
  metadata = mwtab$STUDY)
assay(st003570)
rowData(st003570)
colData(st003570)
