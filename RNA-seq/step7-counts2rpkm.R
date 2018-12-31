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
## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞数量）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
## 下面是定义基因长度为 非冗余exon长度之和
if(F){
  exon_txdb=exons(txdb)
  genes_txdb=genes(txdb)
  
  o = findOverlaps(exon_txdb,genes_txdb)
  o
  t1=exon_txdb[queryHits(o)]
  t2=genes_txdb[subjectHits(o)]
  t1=as.data.frame(t1)
  t1$geneid=mcols(t2)[,1]
  # 如果觉得速度不够，就参考R语言实现并行计算
  # http://www.bio-info-trainee.com/956.html
  g_l = lapply(split(t1,t1$geneid),function(x){
    # x=split(t1,t1$geneid)[[1]]
    head(x)
    tmp=apply(x,1,function(y){
      y[2]:y[3]
    })
    length(unique(unlist(tmp)))
  })
  head(g_l)
  g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
  
  save(g_l,file = 'step7-g_l.Rdata')
}
load(file = 'step7-g_l.Rdata')
## 下面是定义基因长度为 最长转录本长度
if(F){
  
  t_l=transcriptLengths(txdb)
  head(t_l)
  t_l=na.omit(t_l)
  t_l=t_l[order(t_l$gene_id,t_l$tx_len,decreasing = T),]
  str(t_l)
  t_l=t_l[!duplicated(t_l$gene_id),]
  head(t_l)
  g_l=t_l[,c(3,5)]
  
}

head(g_l)
library(org.Mm.eg.db)
s2g=toTable(org.Mm.egSYMBOL)
head(s2g)
g_l=merge(g_l,s2g,by='gene_id')

# 参考counts2rpkm，定义基因长度为非冗余CDS之和
# http://www.bio-info-trainee.com/3298.html 
a[1:4,1:4]
ng=intersect(rownames(a),g_l$symbol)
# 有了counts矩阵和对应的基因长度信息，就很容易进行各种计算了：
exprSet=a[ng,]
lengths=g_l[match(ng,g_l$symbol),'length']
# http://www.biotrainee.com/thread-1791-1-1.html
exprSet[1:4,1:4]
total_count<- colSums(exprSet)
head(total_count)
head(lengths)
10^9/(1122*121297)
rpkm <- t(do.call( rbind,
                   lapply(1:length(total_count),
                          function(i){
  10^9*exprSet[,i]/lengths/total_count[i]
}) ))
rpkm[1:4,1:4]
a=read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz',
             header = T ,sep = '\t')
# 每次都要检测数据
a[1:4,1:4]
rpkm_paper=a[ng,] 
rpkm_paper[1:4,1:4]




