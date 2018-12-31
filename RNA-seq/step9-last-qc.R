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

## 需要确定最后的QC指标


