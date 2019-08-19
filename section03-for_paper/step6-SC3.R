####### SC3 Heatmaps #########
library(scater)
library(SC3)

anno<-data.frame(Plate = unlist(lapply(strsplit(colnames(RPKM),"_"),function(x) x[3])), Population = CAFgroups_full$cluster)
rownames(anno)<-colnames(RPKM)
pd <- new("AnnotatedDataFrame", data = anno)
sceset <-  SingleCellExperiment(assays = list(normcounts = as.matrix(RPKM)) ,colData =anno)

counts(sceset)<-normcounts(sceset)
logcounts(sceset)<- log2(normcounts(sceset)+1)
rowData(sceset)$feature_symbol <- rownames(sceset)
sceset<-sceset[!duplicated(rownames(sceset)),]
sceset <- sc3(sceset, ks = 3:5, biology = TRUE)

pdf("sc3_heatmap.pdf", useDingbats = FALSE, width=6, height = 4.5)

sc3_plot_consensus(
  sceset, k=4,
  show_pdata = c(
    "Plate",
    "Population"
  )
)
dev.off()



matrisome<-as.character(read.table("GitHub/170210_naba_matrisome.txt")[,1])
matrisome<-lookuptable[lookuptable[,1] %in% matrisome,4]

sce.matrisome <- SingleCellExperiment(assays = list(normcounts = as.matrix(RPKM.full[matrisome,])), colData = anno)
counts(sce.matrisome)<-normcounts(sce.matrisome)
logcounts(sce.matrisome)<- log2(normcounts(sce.matrisome)+1)
rowData(sce.matrisome)$feature_symbol <- rownames(sce.matrisome)
sceset<-sce.matrisome[!duplicated(rownames(sce.matrisome)),]


sce.matrisome <- sc3(sce.matrisome, ks = 4, biology = TRUE)

pdf("sce_matrisome_expression.pdf",useDingbats = FALSE, width=6, height = 4.5)
sc3_plot_expression(
  sce.matrisome, k=4,
  show_pdata = c(
    "Plate",
    "Population"
  )
)
dev.off()
