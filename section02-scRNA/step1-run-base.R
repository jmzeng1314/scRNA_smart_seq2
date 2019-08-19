rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
## 首先载入文章的数据
load(file='../input.Rdata')
load(file='../input_rpkm.Rdata')
counts=a
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
dat[1:4,1:4];dim(dat)
dat=log2(dat+1)
library(stringr) 
meta=df
head(meta) 
gs=read.table('top18-genes-in-4-subgroup.txt')[,1]
gs
library(pheatmap)
pheatmap(dat[gs,])
pheatmap(dat[gs,],cluster_rows = F)
tmp=data.frame(g=meta$g)
rownames(tmp)=colnames(dat)
pheatmap::pheatmap(dat[gs,],annotation_col = tmp)


table(apply(dat,2,function(x) sum(x>0) )>2000)
table(apply(dat,1,function(x) sum(x>0) )>4)
dat=dat[apply(dat,1,function(x) sum(x>0) )>4,
        apply(dat,2,function(x) sum(x>0) )>2000]



