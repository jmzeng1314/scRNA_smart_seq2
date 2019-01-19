rm(list = ls()) 
load(file='../input.Rdata')
Sys.setenv(R_MAX_NUM_DLLS=999)
counts=a
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=df
head(meta) 

library(Seurat)
# https://satijalab.org/seurat/mca.html
sce <- CreateSeuratObject(raw.data = counts, 
                          meta.data =meta,
                          min.cells = 5, min.genes = 2000, 
                          project = "sce")
head(sce@meta.data)
dim(sce@data)

VlnPlot(object = sce, features.plot = c("nGene", "nUMI" ), nCol = 3)
VlnPlot(sce,c("Gapdh","Bmp3","Brca1","Brca2","nGene"))
GenePlot(object = sce, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = sce, gene1 = "Brca1", gene2 = "Brca2")
CellPlot(sce,sce@cell.names[3],sce@cell.names[4],do.ident = FALSE)

sce <- NormalizeData(object = sce, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
sce <- FindVariableGenes(object = sce, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length( sce@var.genes)
sce <- ScaleData(object = sce, vars.to.regress = c("nUMI" ))
sce <- RunPCA(object = sce, pc.genes = sce@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
VizPCA( sce, pcs.use = 1:2)
PCAPlot(sce, dim.1 = 1, dim.2 = 2)
sce <- ProjectPCA(sce, do.print = FALSE)
PCHeatmap(object = sce, pc.use = 1, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = sce, pc.use = 1:10, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)
sce <- FindClusters(object = sce, reduction.type = "pca", dims.use = 1:10, force.recalc = T,
                    resolution = 0.8, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(sce)
table(sce@meta.data$res.0.8)


sce <- RunTSNE(object = sce, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = sce)
## firstly find marker one by one by change ident.1 
markers_df <- FindMarkers(object = sce, ident.1 = 1, min.pct = 0.25)
print(x = head(markers_df))
markers_genes =  rownames(head(x = markers_df, n = 5))
VlnPlot(object = sce, features.plot =markers_genes, use.raw = TRUE, y.log = TRUE)
FeaturePlot(object = sce, 
            features.plot =markers_genes, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
# Then find markers for every cluster compared to all remaining cells, report # only the positive ones
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr)
sce.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of# every cell name
DoHeatmap(object = sce, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
save(sce.markers,sce,file='sce_seurat.Rdata')
load(file='sce_seurat.Rdata')
table(sce@meta.data[,5:6])
FeaturePlot(object = sce, 
            features.plot =top10$gene, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

