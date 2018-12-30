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
table(group_list)

cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F,
         filename = 'all_cells_top_1000_sd.png')

if(T){

  
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,
           filename = 'all_cells_top_1000_sd_cutree1.png')
  
}

## 针对top1000的sd的基因集的表达矩阵 进行重新 聚类并且分组。
if(T){

  
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  
  hc=hclust(dist(t(n))) 
  clus = cutree(hc, 4)
  group_list=as.factor(clus)
  table(group_list)
   
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,
           filename = 'all_cells_top_1000_sd_cutree_2.png')
  dev.off()
  
}



## 下面是画PCA的必须操作，需要看说明书。
dat_back=dat
dat=dat_back
dat[1:4,1:4]
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list )
dat[1:4,1:4]
dat[,ncol(dat)]
table(dat$group_list)
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,repel =T,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
            # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
## 事实上还是有很多基因dropout非常严重。
ggsave('all_cells_PCA.png')

