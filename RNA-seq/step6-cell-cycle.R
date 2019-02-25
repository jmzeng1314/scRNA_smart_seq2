## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-12-29 23:24:48
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-12-29  First version
###
### ---------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input.Rdata')
a[1:4,1:4]
head(df) 

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。

group_list=df$g
plate=df$plate
table(plate)
 
a[1:4,1:4]
library(scran)
# https://mp.weixin.qq.com/s/nFSa5hXuKHrGu_othopbWQ
sce <- SingleCellExperiment(list(counts=dat)) 
#list() 创建列表

library(org.Mm.eg.db)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))

ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), 
                  keytype="SYMBOL", column="ENSEMBL")
#取探针名创建一个向量
#rownames(sce) 取行名（即实验检测到的基因）

if(F){
  assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
  save(assigned,file = 'cell_cycle_assigned.Rdata')
}
load(file = 'cell_cycle_assigned.Rdata')
head(assigned$scores)
table(assigned$phases)
draw=cbind(assigned$score,assigned$phases) #合并assigned$score列和assigned$phases列
colnames(draw)
attach(draw)
library(scatterplot3d)
scatterplot3d(G1, S, G2M, angle=20,
              color = rainbow(3)[as.numeric(as.factor(assigned$phases))],
              grid=TRUE, box=FALSE)
detach(draw)

library(pheatmap)
cg=names(tail(sort(apply(dat,1,sd)),100))
n=t(scale(t(dat[cg,])))
# pheatmap(n,show_colnames =F,show_rownames = F)
library(pheatmap)
df$cellcycle=assigned$phases
ac=df
rownames(ac)=colnames(n)
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,
         filename = 'all_cells_top_100_sd_all_infor.png')
dev.off()
head(ac)
table(ac[,c(1,5)])





