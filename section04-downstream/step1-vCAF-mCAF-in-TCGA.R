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
vCAF=unlist(str_split(vCAF,', ') )
mCAF='Dcn, Col12a1, Mmp2, Lum, Mrc2, Bicc1, Lrrc15, Mfap5, Col3A1, Mmp14, Spon1, Pdgfrl, Serpinf1, Lrp1, Gfpt2, Ctsk, Cdh11, Itgbl1, Col6a2, Postn, Ccdc80, Lox, Vcan, Col1a1, Fbn1, Col1a2, Pdpn, Col6a1, Fstl1, Col5a2, Aebp1'
mCAF=unlist(str_split(mCAF,', ') )

if(F){
  library(data.table)
  # 文件BRCA.htseq_counts.tsv.gz从UCSC的XENA数据库下载，大于100M所以不提供在这里。
  # a=fread('TCGA_BRCA/UCSC_xena/TCGA-BRCA.htseq_counts.tsv.gz',data.table=F)
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
  head(esid)
  rownames(a)=esid
  e2s=select(org.Hs.eg.db,keys = esid,columns = c( "ENSEMBL" ,  "SYMBOL" ),keytype = 'ENSEMBL')
  vCAF=toupper(vCAF);vCAF=vCAF[vCAF %in% e2s$SYMBOL]
  mCAF=toupper(mCAF);mCAF=mCAF[mCAF %in% e2s$SYMBOL]
  ng=e2s[match(c(vCAF,mCAF),e2s$SYMBOL),1]
  mat=a[ng,]
  mat=mat[,-1]
  save(mat,file = 'TCGA-BRCA-vCAF-and-mCAF-genes-expression.Rdata')
}

load(file = 'TCGA-BRCA-vCAF-and-mCAF-genes-expression.Rdata')
mat[1:4,1:4] 
M=cor(t(mat))
colnames(M)=c(vCAF,mCAF)
rownames(M)=c(vCAF,mCAF)
pheatmap::pheatmap(M)
pheatmap::pheatmap(M,cluster_rows = F,cluster_cols = F)
n=t(scale(t( M )))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap::pheatmap(n,cluster_rows = F,cluster_cols = F)







