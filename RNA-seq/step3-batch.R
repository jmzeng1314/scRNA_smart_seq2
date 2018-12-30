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

# plate=group_list

## 下面是画PCA的必须操作，需要看说明书。
dat_back=dat

dat=dat_back
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,plate )
dat[1:4,1:4]
table(dat$plate)

# 'princomp' can only be used with more units than variables
# Principal component analysis is underspecified if you have fewer samples than data point. 
# pca_dat =  princomp(t(dat[,-ncol(dat)]))$scores[,1:2]
pca_dat =  prcomp(t(dat[,-ncol(dat)])) 
plot(pca_dat$rotation[,1:2], t='n')
colors = rainbow(length(unique(dat$plate)))
names(colors) = unique(dat$plate)
text(pca_dat$rotation[ , 1:2], labels=dat$plate,col=colors[dat$plate])



library("FactoMineR")
library("factoextra") 
# The variable plate (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,repel =T,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df$g, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
## 事实上还是有很多基因dropout非常严重。
ggsave('all_cells_PCA_by_plate.png')

library(Rtsne) 
dat_matrix <- as.matrix(dat[,-ncol(dat)])
dat_matrix[1:4,1:4]
# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(dat_matrix,pca=FALSE,perplexity=30,theta=0.0) # Run TSNE

# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=dat$plate, asp=1)
# https://distill.pub/2016/misread-tsne/

library(ggpubr)
# Add marginal rug
head(tsne_out$Y)
df=as.data.frame(tsne_out$Y)
colnames(df)=c("X",'Y')
df$plate=dat$plate
head(df)
df$g=group_list
ggscatter(df, x = "X", y = "Y", color = "g" 
          # palette = c("#00AFBB", "#E7B800" ) 
          )


if(F){
  library(tsne)
  ## Not run:
  colors = rainbow(length(unique(iris$Species)))
  names(colors) = unique(iris$Species)
  ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }
  tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)
  head(iris[,1:4])
  
  library(Rtsne)
  iris_unique <- unique(iris) # Remove duplicates
  iris_matrix <- as.matrix(iris_unique[,1:4])
  
  # Set a seed if you want reproducible results
  set.seed(42)
  tsne_out <- Rtsne(iris_matrix,pca=FALSE,perplexity=30,theta=0.0) # Run TSNE
  
  # Show the objects in the 2D tsne representation
  plot(tsne_out$Y,col=iris_unique$Species, asp=1)
  
  
}


