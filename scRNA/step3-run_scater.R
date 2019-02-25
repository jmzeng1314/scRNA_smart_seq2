rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
## 首先载入文章的数据
load(file='../input.Rdata')
counts=a
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=df
head(meta) 

options(warn=-1) # turn off warning message globally
suppressMessages(library(scater))
## 创建 scater 要求的对象
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(counts)), 
  colData = meta
)
sce
exprs(sce) <- log2(
  calculateCPM(sce ) + 1)
## 只有运行了下面的函数后才有各式各样的过滤指标
genes=rownames(rowData(sce))
genes[grepl('^MT-',genes)]
genes[grepl('^ERCC-',genes)]
sce <- calculateQCMetrics(sce, 
                          feature_controls = list(ERCC = grep('^ERCC',genes)))
keep_feature <- rowSums(exprs(sce) > 0) > 5
table(keep_feature)
sce <- sce[keep_feature,]
tf=sce$total_features_by_counts 
boxplot(tf)
fivenum(tf)
table(tf>2000)
sce=sce[,tf > 2000 ]
sce
head(meta)
## 基因表达，理论上应该是跟384孔板 这个变量无关
plotExpression(sce, rownames(sce)[1:6],
               x = "plate", 
               exprs_values = "logcounts") 

# 展示高表达量基因, 绘图有点耗时
plotHighestExprs(sce, exprs_values = "counts")
plotExprsFreqVsMean(sce)

sce <- runPCA(sce)
plotPCA(sce)
reducedDimNames(sce)
# colnames(as.data.frame(colData(sce)))
head(colData(sce))
## PCA分布图上面添加临床信息--------------
plotReducedDim(sce, use_dimred = "PCA", 
                shape_by= "plate", 
                colour_by= "g")
## 考虑 ERCC 影响后继续PCA
sce2 <- runPCA(sce, 
               feature_set = rowData(sce)$is_feature_control)
plotPCA(sce2)
## PCA分布图上面添加临床信息--------------
plotReducedDim(sce2, use_dimred = "PCA", 
               shape_by= "plate", 
               colour_by= "g")

## 运行 tSNE 降维算法
set.seed(1000)
sce <- runTSNE(sce, perplexity=10)
plotTSNE(sce, 
         shape_by= "plate", 
         colour_by= "g")
## 对tSNE降维后结果进行不同的聚类
colData(sce)$tSNE_kmeans <- as.character(kmeans(sce@reducedDims$TSNE,
                                                centers = 4)$clust)
head(sce@reducedDims$TSNE)
hc=hclust(dist( sce@reducedDims$TSNE ))
clus = cutree(hc, 4) 
colData(sce)$tSNE_hc <-  as.character(clus)
plotTSNE(sce,  colour_by = "tSNE_kmeans")
plotTSNE(sce,  colour_by = "tSNE_hc")
table(colData(sce)$tSNE_hc , colData(sce)$tSNE_kmeans)

## 同样是一直降维方式，不同的算法
sce <- runDiffusionMap(sce)
plotDiffusionMap(sce,  
                 shape_by= "plate", 
                 colour_by= "g")


library(SC3) # BiocManager::install('SC3')
sce <- sc3_estimate_k(sce)
metadata(sce)$sc3$k_estimation
rowData(sce)$feature_symbol=rownames(rowData(sce))
# 耗费时间
kn=4
sc3_cluster="sc3_4_clusters"
# 非常耗时
sce <- sc3(sce, ks = kn, biology = TRUE)

sc3_plot_consensus(sce, k = kn, show_pdata = c("g",sc3_cluster))
sc3_plot_expression(sce, k = kn, show_pdata =  c("g",sc3_cluster))
sc3_plot_markers(sce, k = kn, show_pdata =  c("g",sc3_cluster))
plotPCA(sce, shape_by= "g" , colour_by =  sc3_cluster )
sc3_interactive(sce)





