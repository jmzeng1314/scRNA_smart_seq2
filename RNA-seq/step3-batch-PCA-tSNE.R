## 
### ---------------
###
### Create: Jianming Zeng
### Date: nc18-12-29 23:24:48
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: nc18-12-29  First version
###
### ---------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input.Rdata')
a[1:4,1:4]
head(df) 
## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
# 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
group_list=df$g #所有数据的聚类分组信息
plate=df$plate #批次信息
table(plate) 
## 这个时候需要从多个维度来探索两个不同的plate的单细胞群体是否有明显的差别。
# plate=group_list

## 最流行的细胞群体是否有明显的差别，肯定是hclust分群，热图展现，PCA,tSNE 等等

## 如果想了解PCA分析原理，需要阅读：https://mp.weixin.qq.com/s/Kw05PWD2m65TZu2Blhnl4w
## 首先我们使用简单的 prcomp 函数来了解 PCA分析
if(F){
  set.seed(123456789)
  #set.seed()产生随机数
  #用于设定随机数种子，一个特定的种子可以产生一个特定的伪随机序列，这个函数的主要目的，
  #是让你的模拟能够可重复出现，因为很多时候我们需要取随机数，但这段代码再跑一次的时候，
  #结果就不一样了，如果需要重复出现同样的模拟结果的话，就可以用set.seed()。
  library(pheatmap)
  library(Rtsne)
  library(ggfortify)
  library(mvtnorm)
  
  ## 同样的正态分布随机表达矩阵，是无法区分开来。
  if(T){
    ng=500 
    nc=20
    a1=rnorm(ng*nc);dim(a1)=c(ng,nc) #创建正态分布随机矩阵500行，20列
    #dim()检索或设置对象的维度
    a2=rnorm(ng*nc);dim(a2)=c(ng,nc) #因为是随机创建，这两个矩阵不一样
    a3=cbind(a1,a2)
    colnames(a3)=c(paste0('cell_01_',1:nc),
                   paste0('cell_02_',1:nc)) #添加列名
    #paste()粘贴，
    rownames(a3)=paste('gene_',1:ng,sep = '') #添加行名
    pheatmap(a3)
    dist(a3)
    a3=t(a3);dim(a3) ## PCA分析 
    pca_dat <- prcomp(a3, scale. = TRUE) #prcomp()主成分分析
    p=autoplot(pca_dat) + theme_classic() + ggtitle('PCA plot')
    # df=cbind(as.data.frame(a3),group=c(rep('b1',20),rep('b2',20)))
    # autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
    print(p)
    # 可以看到细胞无法被区分开来。
    set.seed(42)
    tsne_out <- Rtsne(a3,pca=FALSE,perplexity=10,theta=0.0) # Run TSNE
    tsnes=tsne_out$Y
    colnames(tsnes) <- c("tSNE1", "tSNE2") #添加列名
   # group=c(rep('b1',20),rep('b2',20))
   # tsnes=as.data.frame(tsnes)
   #  tsnes$group=group
   #  ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=group))
    ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point()
  }
  
  ## 同样的正态分布随机表达矩阵，但是其中部分细胞+3，可以区分开来。
  if(T){
    ng=500
    nc=20
    a1=rnorm(ng*nc);dim(a1)=c(ng,nc)
    a2=rnorm(ng*nc)+3;dim(a2)=c(ng,nc) 
    a3=cbind(a1,a2)
    colnames(a3)=c(paste0('cell_01_',1:nc),paste0('cell_02_',1:nc))
    rownames(a3)=paste('gene_',1:ng,sep = '')
    pheatmap(a3)
    a3=t(a3);dim(a3) ## PCA分析 
    
    pca_dat <- prcomp(a3, scale. = TRUE)
    p=autoplot(pca_dat) + theme_classic() + ggtitle('PCA plot')
    print(p)
    # 这个时候细胞被区分开，而且是很明显的一个主成分。
    
    set.seed(42)
    tsne_out <- Rtsne(a3,pca=FALSE,perplexity=10,theta=0.0) # Run TSNE
    tsnes=tsne_out$Y
    colnames(tsnes) <- c("tSNE1", "tSNE2")
    ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point()
  }
  
  ## 不同的正态分布随机表达矩阵，可以区分。
  if(T){
    ng=600
    nc=200
    mu1  = rnorm(ng, mean = 1)
    #rnorm()函数产生一系列的随机数，随机数个数，均值和标准差都可以设定
    mu2  = rnorm(ng, mean = 5)
    a1=rmvnorm(nc,mu1);dim(a1) #rmvnorm随机生成多元正太分布数
    a2=rmvnorm(nc,mu2) ;dim(a2)
    
    a3=rbind(a1,a2);dim(a3) #rbind()行进行合并，就是行的叠加
    rownames(a3)=c(paste0('cell_01_',1:nc),paste0('cell_02_',1:nc))
    colnames(a3)=paste('gene_',1:ng,sep = '')
    pheatmap(a3)
    
    pca_dat <- prcomp(a3, scale. = TRUE)
    p=autoplot(pca_dat) + theme_classic() + ggtitle('PCA plot')
    print(p)
    set.seed(42)
    tsne_out <- Rtsne(a3,pca=FALSE,perplexity=30,theta=0.0) # Run TSNE
    tsnes=tsne_out$Y
    colnames(tsnes) <- c("tSNE1", "tSNE2")
    ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point()
  }
  
  
}

## 然后我们 可以使用高级R包做真实的分析。
## 下面是画PCA的必须操作，需要看不同做PCA的包的说明书。
dat_back=dat

dat=dat_back
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,plate ) #cbind根据列进行合并，即叠加所有列 #矩阵添加批次信息
dat[1:4,1:4]
table(dat$plate)
# 
# # 'princomp' can only be used with more units than variables
# # Principal component analysis is underspecified if you have fewer samples than data point. 
# # pca_dat =  princomp(t(dat[,-ncol(dat)]))$scores[,1:2]
# pca_dat =  prcomp(t(dat[,-ncol(dat)])) 
# 
# plot(pca_dat$rotation[,1:2], t='n')
# colors = rainbow(length(unique(dat$plate)))
# names(colors) = unique(dat$plate)
# text(pca_dat$rotation[ , 1:2], labels=dat$plate,col=colors[dat$plate])
# 

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input_rpkm.Rdata')
dat[1:4,1:4]
head(metadata) 
plate=metadata$plate
## 下面是画PCA的必须操作，需要看不同做PCA的包的说明书。
dat_back=dat

dat=dat_back
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,plate ) #cbind根据列进行合并，即叠加所有列 #矩阵添加批次信息
dat[1:4,1:4]
table(dat$plate)
library("FactoMineR")
library("factoextra") 
# The variable plate (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,#repel =T,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$plate, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
) 
ggsave('all_cells_PCA_by_plate.png') 
## 保存你的画布的图到本地图片文件。

## scater, monocle3, seurat 

library(Rtsne) 
load(file = '../input_rpkm.Rdata')
dat[1:4,1:4]
dat_matrix <- t(dat)

if(F){
  load(file = '../input.Rdata')
  # dat is log2(cpm+1)
  dat_matrix=t(dat)
  dat_matrix[1:4,1:4]
  
}

dat_matrix =log2(dat_matrix+0.01)
# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(dat_matrix,pca=FALSE,perplexity=30,theta=0.0) # Run TSNE
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col= plate) # asp=1
# https://distill.pub/nc16/misread-tsne/



library(ggpubr)
# Add marginal rug
head(tsne_out$Y)
df=as.data.frame(tsne_out$Y)
colnames(df)=c("X",'Y')
df$plate= plate
head(df)
df$g=metadata$g
ggscatter(df, x = "X", y = "Y", color = "g" 
          # palette = c("#00AFBB", "#E7B800" ) 
          )

## 下面是想说明一个错误的例子。
if(F){
  library(tsne)
  ## Not run:
  data(iris)
  iris
  colors = rainbow(length(unique(iris$Species)))
  names(colors) = unique(iris$Species)
  colors
  ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, 
                                           col=colors[iris$Species]) }
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


