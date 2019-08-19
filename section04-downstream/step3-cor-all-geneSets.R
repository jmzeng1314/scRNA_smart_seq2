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
# for 参考文献ref32：proliferation signature
proliferation <- c("BAG1","ESR1","FOXA1","GPR160","NAT1","NAT1","MAPT","MLPH","PGR")
# for ref31: stroma-related signatures
stroma <- c('DCN', 'CSPG2', 'CDH11', 'COL3A1', 'FAP', 'PEDF', 'FBN1', 'PDGFRL', 'CTSK', 'HTRA1', 'ASPN', 'SPARC', 'COL5A2', 'LOXL1', 'MMP2', 'SPON1', 'SFRP4', 'ITGBL1', 'CALD1', 'COPZ2', 'MFAP2', 'ANGPTL2', 'PLAU', 'COL1A2', 'LRRC17', 'C1QTNF3', 'SNAI2', 'PCOLCE', 'POSTN', 'ECM2', 'FBLN1', 'ADAM12', 'MMP11', 'AEBP1', 'PDGFRB', 'GAS1', 'COL6A3', 'RARRES2', 'COL6A1', 'TGFB3', 'NDN', 'C1R', 'LRP1', 'COL10A1', 'DPYSL3', 'OLFML2B', 'MMP14', 'DACT1', 'MGC3047', 'THBS2')
# for ref29: endothelial/microvasculature signatures
microvasculature <- c('ARAP3','ADCY4','ESAM','ERG','SLC43A3','SOX7','PTPRB','PTPRM','AFAP1L1','MMRN2','TENC1','STARD9','COL4A3','LRRK1','PALD1','NPR3','ROBO4','NOTCH4','TIE1','RASIP1','ACVRL1','RAMP2','FAM110D','EGFL7','SMAD6','FGD5','ENG','CASKIN2','ACKR2','SLC9A3R2','CALCRL','HSPA12B','EPAS1','EHD4','LATS2','ICAM1','HBEGF','PLTP','C1orf54','CTTNBP2NL','MYO1B','SLCO2A1','KIFC1','EPHB4','SOX13','DRAM1','PECAM1','ENTPD1','ICAM2','CLDN5','SDPR','CDH5','GPR116','ELTD1','KDR','HILPDA','NPNT')



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
  vCAF=toupper(vCAF);vCAF=vCAF[vCAF %in% e2s$SYMBOL]
  mCAF=toupper(mCAF);mCAF=mCAF[mCAF %in% e2s$SYMBOL]
  ng=e2s[match(vCAF,e2s$SYMBOL),1]
  vCAF_value=colMeans(a[ng,])
  ng=e2s[match(mCAF,e2s$SYMBOL),1]
  mCAF_value=colMeans(a[ng,])
  ng=e2s[match(ECM,e2s$SYMBOL),1]
  ECM_value=colMeans(a[ng,])
  ng=e2s[match(endothelial,e2s$SYMBOL),1];ng
  endothelial_value=colMeans(a[ng,])
  ng=e2s[match(proliferation,e2s$SYMBOL),1];ng=ng[!is.na(ng)];ng
  
  
  proliferation_value=colMeans(a[ng,])
  ng=e2s[match(stroma,e2s$SYMBOL),1];ng=ng[!is.na(ng)];ng
  stroma_value=colMeans(a[ng,])
  ng=e2s[match(microvasculature,e2s$SYMBOL),1];ng=ng[!is.na(ng)];ng
  microvasculature_value=colMeans(a[ng,])
  
  dat=data.frame(vCAF=vCAF_value,
                 mCAF=mCAF_value,
                 ECM=ECM_value,
                 endothelial=endothelial_value,
                 microvasculature=microvasculature_value,
                 stroma=stroma_value,
                 proliferation=proliferation_value)
  save(dat,file = 'step3.Rdata')
}


load(file = 'step3.Rdata')
dat[1:4,1:4] 
M=cor( dat) 
pheatmap::pheatmap(M)


 