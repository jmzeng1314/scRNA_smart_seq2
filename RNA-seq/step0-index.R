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
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

# http://www.bio-info-trainee.com/3727.html 周末生物信息学培训班准备工作

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
BiocManager::install(c('org.Mm.eg.db','scran'))

options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror
BiocManager::install(c( 'scran'))


if(F){
  a=read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz',
               header = T ,sep = '\t')
  # 每次都要检测数据
  a[1:4,1:4]
  ## 读取RNA-seq定量结果，表达矩阵需要进行简单的过滤
  dat=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] 
  ## 2 % 的细胞有表达量
  dat[1:4,1:4]
  dat=log2(edgeR::cpm(dat)+1)
  dat[1:4,1:4] 
  # STATQUEST 
  hc=hclust(dist(t(dat))) 
  clus = cutree(hc, 4)
  group_list= as.factor(clus) 
  table(group_list)
  
  colnames(dat)
  library(stringr)
  plate=str_split(colnames(dat),'_',simplify = T)[,3]
  table(plate)
  n_g = apply(a,2,function(x) sum(x>1)) 
  
  df=data.frame(g=group_list,plate=plate,n_g=n_g)
  df$all='all'
  
  save(a,dat,df,file = '../input.Rdata')
} 

load(file = '../input.Rdata')
a[1:4,1:4]


