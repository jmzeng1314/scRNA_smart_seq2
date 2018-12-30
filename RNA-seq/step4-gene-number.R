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

group_list=df$g
plate=df$plate
table(plate)
 
library(ggpubr)
ggviolin(df, x = "all", y = "n_g", fill = "all", 
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "plate", y = "n_g", fill = "plate",
         #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white")) 
library(ggpubr)
ggviolin(df, x = "g", y = "n_g", fill = "g", 
         add = "boxplot", add.params = list(fill = "white")) 
