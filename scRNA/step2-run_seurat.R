rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
## 首先载入文章的数据
load(file='../input.Rdata')
counts=a
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=df
head(meta) 
# 下面的基因是文章作者给出的
gs=read.table('top18-genes-in-4-subgroup.txt')[,1]
gs
library(pheatmap)


fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
# 上面检测了 counts 和 meta 两个变量，后面需要使用

library(Seurat)
# https://satijalab.org/seurat/mca.html
# 构建 Seurat 需要的对象
# 其中 min.cells 和 min.genes 两个参数是经验值
sce <- CreateSeuratObject(raw.data = counts, 
                          meta.data =meta,
                          min.cells = 5, 
                          min.genes = 2000, 
                          project = "sce")
# 参考：https://github.com/satijalab/seurat/issues/668

sce
?seurat
table(apply(counts,2,function(x) sum(x>0) )>2000)
table(apply(counts,1,function(x) sum(x>0) )>4)
## 可以看到上面的过滤参数是如何发挥作用的
head(sce@meta.data) 
dim(sce@data)

## 默认使用细胞名字字符串的切割第一列来对细胞进行分组
# 所以在这里只有一个SS2分组信息, 我们就增加一个 group.by 参数
VlnPlot(object = sce, 
        features.plot = c("nGene", "nUMI" ), 
        group.by = 'plate',
        nCol = 2)
VlnPlot(object = sce, 
        features.plot = c("nGene", "nUMI" ), 
        group.by = 'g',
        nCol = 2)
### 同样的发现，普通的层次聚类得到的4组，很明显是检测到的基因数量的差异造成的。

## 可以给sce对象增加一个属性，供QC使用
ercc.genes <- grep(pattern = "^ERCC-", x = rownames(x = sce@raw.data), value = TRUE)
percent.ercc <- Matrix::colSums(sce@raw.data[ercc.genes, ]) / Matrix::colSums(sce@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
sce <- AddMetaData(object = sce, metadata = percent.ercc,
                   col.name = "percent.ercc")
VlnPlot(object = sce, 
        features.plot = c("nGene", "nUMI", "percent.ercc" ), 
        group.by = 'g',
        nCol = 3)
## 发现一个很有趣的现象，细胞能检测到的基因数量与其含有的ERCC序列反相关。


## 下面是一些可视化函数
VlnPlot(sce,group.by = 'plate',c("Gapdh","Bmp3","Brca1","Brca2","nGene"))

GenePlot(object = sce,  gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = sce,  gene1 = "Brca1", gene2 = "Brca2")
CellPlot(sce,sce@cell.names[3], 
         sce@cell.names[4],
         do.ident = FALSE)

sce <- NormalizeData(object = sce, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)
sce@data[1:4,1:4]
pheatmap(as.matrix(sce@data[gs,]))
# 这里需要理解  dispersion 值
sce <- FindVariableGenes(object = sce, 
                         mean.function = ExpMean, 
                         dispersion.function = LogVMR )
## 默认值是：x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1
# 需要认真读说明书：https://satijalab.org/seurat/pbmc3k_tutorial.html

# This function is unchanged from (Macosko et al.), 
# but new methods for variable gene expression identification are coming soon. 
# We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. 


## 根据经验阈值挑选的变化基因个数。
length( sce@var.genes)

# Scaling the data and removing unwanted sources of variation

## 去除一些技术误差，比如 nUMI或者ERCC
head(sce@meta.data) 
sce <- ScaleData(object = sce, 
                 vars.to.regress = c("nUMI",'nGene',"percent.ercc" ))
# 后面就不需要考虑ERCC序列了。
sce@scale.data[1:4,1:4]
pheatmap(as.matrix(sce@scale.data[gs,])) 

## 普通的PCA降维
sce <- RunPCA(object = sce, pc.genes = sce@var.genes, 
              do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
tmp <- sce@dr$pca@gene.loadings
head(tmp)
VizPCA( sce, pcs.use = 1:2)
PCAPlot(sce, dim.1 = 1, dim.2 = 2,group.by = 'plate')
PCAPlot(sce, dim.1 = 1, dim.2 = 2,group.by = 'g')
sce <- ProjectPCA(sce, do.print = FALSE)

PCHeatmap(object = sce, pc.use = 1, cells.use = 100, 
          do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = sce, pc.use = 1:10, cells.use = 100, 
          do.balanced = TRUE, label.columns = FALSE)

### 根据参数来调整最后的分组个数
sce1 <- FindClusters(object = sce, reduction.type = "pca", 
                    dims.use = 1:20, force.recalc = T,
                    resolution = 0.4, print.output = 0, 
                    save.SNN = TRUE)
PrintFindClustersParams(sce1)
table(sce1@meta.data$res.0.4) 
## resolution 是最关键的参数
sce=sce1

sce <- RunTSNE(object = sce, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = sce)
## firstly find marker one by one by change ident.1 
markers_df <- FindMarkers(object = sce, ident.1 = 1, 
                          min.pct = 0.25)
print(x = head(markers_df))
markers_genes =  rownames(head(x = markers_df, n = 5))
# 可视化最后找到的marker基因
VlnPlot(object = sce, features.plot =markers_genes, 
        use.raw = TRUE, y.log = TRUE)
# 首先可视化找到的marker基因
FeaturePlot(object = sce, 
            features.plot =markers_genes, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
# 然后可视化文献作者给出的基因
FeaturePlot(object = sce, 
            features.plot =gs[1:18], 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
FeaturePlot(object = sce, 
            features.plot =gs[19:36], 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
FeaturePlot(object = sce, 
            features.plot =gs[37:54], 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
FeaturePlot(object = sce, 
            features.plot =gs[55:72], 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


## 对每个类别细胞都找到自己的marker基因
# Then find markers for every cluster compared to all remaining cells, report # only the positive ones
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr)
sce.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
intersect(top10$gene,gs)
# setting slim.col.label to TRUE will print just the cluster IDS instead of# every cell name
DoHeatmap(object = sce, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
# save(sce.markers,sce,file='sce_seurat.Rdata')
# load(file='sce_seurat.Rdata') 
FeaturePlot(object = sce, 
            features.plot =top10$gene, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

top20 <- sce.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
intersect(top20$gene,gs)
DoHeatmap(object = sce, genes.use = top20$gene, 
          slim.col.label = TRUE, remove.key = TRUE)

 





