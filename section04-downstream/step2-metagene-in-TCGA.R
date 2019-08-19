## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-08-16 22:43:41
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-08-16  First version
###
### ---------------


library(stringr)
vCAF='Esam, Gng11, Higd1b, Cox4i2, Cygb, Gja4, Eng'
vCAF=unlist(str_split(vCAF,', ')[[1]])
mCAF='Dcn, Col12a1, Mmp2, Lum, Mrc2, Bicc1, Lrrc15, Mfap5, Col3A1, Mmp14, Spon1, Pdgfrl, Serpinf1, Lrp1, Gfpt2, Ctsk, Cdh11, Itgbl1, Col6a2, Postn, Ccdc80, Lox, Vcan, Col1a1, Fbn1, Col1a2, Pdpn, Col6a1, Fstl1, Col5a2, Aebp1'
mCAF=unlist(str_split(mCAF,', ')[[1]])

ECM=c('COL1A1', 'COL1A2','COL3A1')
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=DIPK2B
# CXorf36  --> DIPK2B
endothelial=c('CDH5', 'DIPK2B','TIE1')

if(F){
  library(data.table)
  # 文件BRCA.htseq_counts.tsv.gz从UCSC的XENA数据库下载，大于100M所以不提供在这里。
  # a=fread('../../../TCGA_BRCA/UCSC_xena/TCGA-BRCA.htseq_counts.tsv.gz',data.table=F)
  a[1:4,1:4]
  # Ensembl_ID TCGA-3C-AAAU-01A TCGA-3C-AALI-01A TCGA-3C-AALJ-01A
  # 1 ENSG00000000003.13         9.348728         8.714246        10.356452
  # 2  ENSG00000000005.5         1.584963         1.584963         5.727920
  # 3 ENSG00000000419.11        10.874981        10.834471        10.329796
  # 4 ENSG00000000457.12        10.121534        11.512247         8.867279
  library(org.Hs.eg.db)
  library(stringr)
  esid=str_split(a$Ensembl_ID,
                 '[.]',simplify = T)[,1]
  columns(org.Hs.eg.db)
  rownames(a)=esid
  a=a[,-1]
  e2s=select(org.Hs.eg.db,keys = esid,columns = c( "ENSEMBL" ,  "SYMBOL" ),keytype = 'ENSEMBL')
  vCAF=toupper(vCAF);vCAF=vCAF[vCAF %in% e2s$SYMBOL,]
  mCAF=toupper(mCAF);mCAF=mCAF[mCAF %in% e2s$SYMBOL,]
  ng=e2s[match(vCAF,e2s$SYMBOL),1]
  vCAF_value=colMeans(a[ng,])
  ng=e2s[match(mCAF,e2s$SYMBOL),1]
  mCAF_value=colMeans(a[ng,])
  
  ng=e2s[match(ECM,e2s$SYMBOL),1]
  ECM_value=colMeans(a[ng,])
  ng=e2s[match(endothelial,e2s$SYMBOL),1]
  endothelial_value=colMeans(a[ng,])
  dat=data.frame(vCAF_value=vCAF_value,
                 mCAF_value=mCAF_value,
                 ECM_value=ECM_value,
                 endothelial_value=endothelial_value )
  save(dat,file = 'step2.Rdata')
}

load(file = 'step2.Rdata')
## 两种方法批量绘制散点图（多图组合）
## https://mp.weixin.qq.com/s/SG1VvfL8lVCniGCjW2DqfA 
## 画boxplot 
library(ggpubr)
colnames(dat)
ggscatter(dat, x = "vCAF_value", y = "endothelial_value",
          color = 'black', shape = 21, size = 0.5, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,  
          cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)


