###
### ---------------

rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor


options()$repos  ## 查看使用install.packages安装时的默认镜像
options()$BioC_mirror ##查看使用bioconductor的默认镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") ##指定镜像，这个是中国科技大学镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##指定install.packages安装镜像，这个是清华镜像
options()$repos 
options()$BioC_mirror

# http://www.bio-info-trainee.com/3727.html 周末生物信息学培训班准备工作
 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") ##判断是否存在BiocManager包，不存在的话安装

if(!require('Seurat')){
  BiocManager::install('Seurat',ask = F,update = F)
}
if(!require('scran')){
  BiocManager::install(c( 'scran'),ask = F,update = F)
}
if(!require('monocle')){
  BiocManager::install(c( 'monocle'),ask = F,update = F)
}


## http://cole-trapnell-lab.github.io/monocle-release/monocle3/ 
BiocManager::install('destiny')
BiocManager::install(c( 'flexclust','mcclust'),ask = F,update = F)

BiocManager::install(c( 'scRNAseq'),ask = F,update = F)
BiocManager::install(c( 'dbscan'),ask = F,update = F)
BiocManager::install(c( 'M3Drop'),ask = F,update = F)


rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(reticulate))
suppressMessages(library(devtools))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))

