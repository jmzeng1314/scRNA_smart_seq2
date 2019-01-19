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
# 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
group_list=df$g
plate=df$plate
table(plate)
# n_g = apply(a,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）

### 绘图必须强推 ggpubr 
## 搜索到代码，然后修改即可出美图。
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input_rpkm.Rdata')
#   n_g = apply(a,2,function(x) sum(x>0)) #统计每个样本有表达的有多少行（基因）

dat[1:4,1:4]
df=metadata
head(df) 
library(ggpubr)
ggviolin(df, x = "all", y = "n_g", fill = "all", 
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "plate", y = "n_g", fill = "plate",
         #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "g", y = "n_g", fill = "g", 
         add = "boxplot", add.params = list(fill = "white"))  + stat_compare_means()



