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
table(df$g)
## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞数量）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。
exprSet=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] 
exprSet=exprSet[!grepl('ERCC',rownames(exprSet)),]
group_list=ifelse(df$g==1,'me','other')
table(group_list)
library(edgeR)
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('me-other'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='me-other', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
}
boxplot(dat['Shank1',]~df$g)


group_list=ifelse(df$g==1,'me','other')
table(group_list)
library(edgeR)
do_limma_RNAseq <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('me-other'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='me-other', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
  return(DEG_limma_voom)
}
deg1=do_limma_RNAseq(exprSet,group_list)
boxplot(dat['Shank1',]~df$g)
boxplot(dat['Prelid1',]~df$g)
boxplot(dat[rownames(deg1)[1],]~df$g)

group_list=ifelse(df$g==2,'me','other');table(group_list)
deg2=do_limma_RNAseq(exprSet,group_list)
boxplot(dat[rownames(deg2)[1],]~df$g)
boxplot(dat[rownames(deg2)[2],]~df$g)

group_list=ifelse(df$g==3,'me','other');table(group_list)
deg3=do_limma_RNAseq(exprSet,group_list)

group_list=ifelse(df$g==4,'me','other');table(group_list)
deg4=do_limma_RNAseq(exprSet,group_list)

deg1=deg1[order(deg1$logFC,decreasing = T),]
deg2=deg2[order(deg2$logFC,decreasing = T),]
deg3=deg3[order(deg3$logFC,decreasing = T),]
deg4=deg4[order(deg4$logFC,decreasing = T),]

cg=c(head(rownames(deg1),18),
     head(rownames(deg2),18),
     head(rownames(deg3),18),
     head(rownames(deg4),18)
     )
library(pheatmap)
g=df$g
mat=dat[cg,]
mat=mat[,order(g)]
ac=data.frame(group=g)
rownames(ac)=colnames(dat)
mat=mat[head(rownames(deg1),18),]

n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]

pheatmap(n,show_rownames = T,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         annotation_col = ac)



plot(deg1$logFC,-log10(deg1$adj.P.Val))
with(deg1,plot( logFC,-log10( adj.P.Val)))
with(deg2,plot( logFC,-log10( adj.P.Val)))
with(deg3,plot( logFC,-log10( adj.P.Val)))
with(deg4,plot( logFC,-log10( adj.P.Val)))

diff1=rownames(deg1[abs(deg1$logFC)>3,])
diff2=rownames(deg2[abs(deg2$logFC)>3,])
diff3=rownames(deg3[abs(deg3$logFC)>3,])
diff4=rownames(deg4[abs(deg4$logFC)>3,])

library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
df <- bitr(diff1, fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
kk  <- enrichKEGG(gene         = df$ENTREZID,
                    organism     = 'mmu', 
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk)[,1:6]
dotplot(kk)























