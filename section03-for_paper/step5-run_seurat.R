rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
## 首先counts数据
load(file='expr_raw_counts.Rdata')
counts=expr_raw
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))

dat=log2(edgeR::cpm(counts)+1) 
hc=hclust(dist(t(dat)))  
plot(hc,labels = FALSE)  
clus = cutree(hc, 4) 
group_list= as.factor(clus) 
table(group_list) 

#提取批次信息
colnames(dat) #取列名
library(stringr)
plate=str_split(colnames(dat),'_',simplify = T)[,3] #取列名，以'_'号分割，提取第三列。
#str_split()函数可以分割字符串
table(plate)

n_g = apply(counts,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）
# 这里我们定义， reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
# 而且大部分单细胞转录组技术很烂，通常超过75%的基因都没办法检测到。

df=data.frame(g=group_list,plate=plate,n_g=n_g) #新建数据框(细胞的属性信息)
library(stringr) 
meta=df
head(meta)  

# 上面得到的 counts 和 meta 两个变量，Seurat 需要使用

library(Seurat)
# https://satijalab.org/seurat/mca.html
# 构建 Seurat 需要的对象
# 其中 min.cells 和 min.genes 两个参数是经验值
sce <- CreateSeuratObject(raw.data = counts, 
                          meta.data =meta,
                          min.cells = 5, 
                          min.genes = 1000, 
                          project = "sce")
# 参考：https://github.com/satijalab/seurat/issues/668

sce
?seurat
table(apply(counts,2,function(x) sum(x>0) )>1000)
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


sce <- NormalizeData(object = sce, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)
sce@data[1:4,1:4] 
# 这里需要理解  dispersion 值
sce <- FindVariableGenes(object = sce, 
                         mean.function = ExpMean, 
                         dispersion.function = LogVMR )
## 根据经验阈值挑选的变化基因个数。
length( sce@var.genes)
## 去除一些技术误差，比如 nUMI或者ERCC
head(sce@meta.data) 
sce <- ScaleData(object = sce, 
                 vars.to.regress = c("nUMI",'nGene' ))
# 后面就不需要考虑ERCC序列了。
sce@scale.data[1:4,1:4] 
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

### 根据参数来调整最后的分组个数
sce1 <- FindClusters(object = sce, reduction.type = "pca", 
                    dims.use = 1:20, force.recalc = T,
                    resolution = 0.4, print.output = 0, 
                    save.SNN = TRUE)
PrintFindClustersParams(sce1)
table(sce1@meta.data$res.0.4) 
table(sce1@meta.data$res.0.4,sce1@meta.data$g) 
## resolution 是最关键的参数
sce=sce1

sce <- RunTSNE(object = sce, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = sce) 
## 对每个类别细胞都找到自己的marker基因
# Then find markers for every cluster compared to all remaining cells, report # only the positive ones
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr)
sce.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

head(meta) 
# 下面的基因是文章作者给出的
gs=read.table('top18-genes-in-4-subgroup.txt')[,1]
gs
library(pheatmap)
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

 





