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

if(F){
  
  load(file = '../input.Rdata')
  a[1:4,1:4]
  head(metadata) #head()函数显示操作前面的信息，默认前6行
  
  ## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
  # 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
  
  group_list=metadata$g #'$'符，取列，取metadata矩阵的g列,取出层级聚类信息
  table(group_list) ##这是全部基因集的聚类分组信息
  
  ## 要完全掌握本代码，需要理解什么是热图，什么是PCA图，后面会单独讲解。
  
  cg=names(tail(sort(apply(dat,1,sd)),100)) ##取表达量标准差最大的100行的行名
  # 这个前面演示过一次，  dat=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] #筛选表达量合格的行,列数不变
  # 大家可以对照着理解 apply，需要一定的R语言功底
  
  #对dat矩阵每行求标准差，排序，取最后的100行，并获得后100行的行名（探针名）
  #sd()求标准差，对dat矩阵每一行的counts求标准差
  #sort()函数,排序;
  #tail()函数，显示操作对象后面的信息，默认后6行，这里设定取后100行
  #names()函数，获取或设置对象的名称
  library(pheatmap)
  
  ##画热图,针对top100的sd的基因集的表达矩阵,没有聚类分组
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F,
           filename = 'all_cells_top_100_sd.png')
  
  ##针对top100的sd的基因集的表达矩阵，归一化，分组画图
  if(T){
    
    x=1:10;plot((x))
    scale(x);plot(scale(x))
    
    n=t(scale(t(dat[cg,]))) #scale()函数去中心化和标准化
    #对每个探针的表达量进行去中心化和标准化
    n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
    n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
    n[1:4,1:4]
    
    # pheatmap(n,show_colnames =F,show_rownames = F)
    
    ac=data.frame(g=group_list) #制作细胞（样本）分组矩阵
    rownames(ac)=colnames(n) ##ac的行名（样本名）等于n的列名（样本名）
    ##判断分组矩阵的行（样本数）和表达矩阵的列（样本数）是否相等
    pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac,
             filename = 'all_cells_top_100_sd_cutree1.png')
    
  }
  
  ## 针对top100的sd的基因集的表达矩阵 进行重新 聚类并且分组。
  if(T){
    
    
    n=t(scale(t(dat[cg,])))
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    
    ##这个聚类分组只是对top100的sd的基因集
    hc=hclust(dist(t(n))) 
    clus = cutree(hc, 4)
    group_list=as.factor(clus)
    table(group_list) ##这个聚类分组信息是针对top100的sd的基因集的，和针对全部基因集的分组结果不一样
    table(group_list,metadata$g) ## 其中 metadata$g 是前面步骤针对全部表达矩阵的层次聚类结果。
    
    ## 下面针对本次挑选100个基因的表达矩阵的层次聚类结果进行热图展示。
    ac=data.frame(g=group_list)
    rownames(ac)=colnames(n)
    pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac,
             filename = 'all_cells_top_100_sd_cutree_2.png')
    dev.off() ##关闭画板
    
  }
  
  ###先对整个基因集聚类分组，再对top100的sd的基因集画图
  ###和选取top100的sd的基因集，聚类分组再画图，结果聚类分组信息不同
  ###两次的聚类分组信息不同，画出的图不同
  
  ## 下面是画PCA的必须操作，需要看不同做PCA的包的说明书。
  ## 如果想了解PCA分析原理，需要阅读：https://mp.weixin.qq.com/s/Kw05PWD2m65TZu2Blhnl4w
  
  dat_back=dat ##防止下面操作把数值搞坏的一个备份
  dat=dat_back  ##表达矩阵数据
  dat[1:4,1:4]
  dat=t(dat)
  dat=as.data.frame(dat) ##转换为数据框
  dat=cbind(dat,group_list ) ##cbind()合并列（横向追加）;添加分组信息
  dat[1:4,1:4]
  ## 表达矩阵可以随心所欲的取行列，基础知识需要打牢。
  dat[1:4,12197:12199]
  dat[,ncol(dat)] #ncol()列，返回列长值
  table(dat$group_list)
  library("FactoMineR")
  library("factoextra") 
  # The variable group_list (index = ) is removed
  # before PCA analysis
  ## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #'-'表示“非”
  fviz_pca_ind(dat.pca,repel =T,
               geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
               col.ind = dat$group_list, # color by groups 颜色组
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses 集中成椭圆
               legend.title = "Groups"
  )
  ## 事实上还是有很多基因dropout非常严重。
  ggsave('all_cells_PCA.png')
  
  
  
  
}

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
if(F){
  
  load(file = '../input_rpkm.Rdata')
  # a[1:4,1:4]
  dat[1:4,1:4]
  head(metadata) #head()函数显示操作前面的信息，默认前6行
  
  ## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
  # 注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。
  
  group_list=metadata$g #'$'符，取列，取metadata矩阵的g列,取出层级聚类信息
  table(group_list) ##这是全部基因集的聚类分组信息
  
  ## 要完全掌握本代码，需要理解什么是热图，什么是PCA图，后面会单独讲解。
  
  cg=names(tail(sort(apply(dat,1,sd)),100)) ##取表达量标准差最大的100行的行名
  # 这个前面演示过一次，  dat=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] #筛选表达量合格的行,列数不变
  # 大家可以对照着理解 apply，需要一定的R语言功底
  
  #对dat矩阵每行求标准差，排序，取最后的100行，并获得后100行的行名（探针名）
  #sd()求标准差，对dat矩阵每一行的counts求标准差
  #sort()函数,排序;
  #tail()函数，显示操作对象后面的信息，默认后6行，这里设定取后100行
  #names()函数，获取或设置对象的名称
  library(pheatmap)
  
  mat=log2(dat[cg,]+0.01)
  ##画热图,针对top100的sd的基因集的表达矩阵,没有聚类分组
  pheatmap(mat,show_colnames =F,show_rownames = F,
           filename = 'all_cells_top_100_sd.png')
  
  ##针对top100的sd的基因集的表达矩阵，归一化，分组画图
  if(T){
    
    x=1:10;plot((x))
    scale(x);plot(scale(x))
    
    n=t(scale(t( mat ))) #scale()函数去中心化和标准化
    #对每个探针的表达量进行去中心化和标准化
    n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
    n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
    n[1:4,1:4]
    
    # pheatmap(n,show_colnames =F,show_rownames = F)
    
    ac=data.frame(g=group_list) #制作细胞（样本）分组矩阵
    rownames(ac)=colnames(n) ##ac的行名（样本名）等于n的列名（样本名）
    ##判断分组矩阵的行（样本数）和表达矩阵的列（样本数）是否相等
    pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac,
             filename = 'all_cells_top_100_sd_cutree1.png')
    
  }
  
  ## 针对top100的sd的基因集的表达矩阵 进行重新 聚类并且分组。
  if(T){
    
    
    n=t(scale(t( mat )))
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    
    ##这个聚类分组只是对top100的sd的基因集
    hc=hclust(dist(t(n))) 
    clus = cutree(hc, 4)
    group_list=as.factor(clus)
    table(group_list) ##这个聚类分组信息是针对top100的sd的基因集的，和针对全部基因集的分组结果不一样
    table(group_list,metadata$g) ## 其中 metadata$g 是前面步骤针对全部表达矩阵的层次聚类结果。
    
    ## 下面针对本次挑选100个基因的表达矩阵的层次聚类结果进行热图展示。
    ac=data.frame(g=group_list)
    rownames(ac)=colnames(n)
    pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac,
             filename = 'all_cells_top_100_sd_cutree_2.png')
    dev.off() ##关闭画板
    
  }
  
  ###先对整个基因集聚类分组，再对top100的sd的基因集画图
  ###和选取top100的sd的基因集，聚类分组再画图，结果聚类分组信息不同
  ###两次的聚类分组信息不同，画出的图不同
  
  ## 下面是画PCA的必须操作，需要看不同做PCA的包的说明书。
  ## 如果想了解PCA分析原理，需要阅读：https://mp.weixin.qq.com/s/Kw05PWD2m65TZu2Blhnl4w
  
  dat_back=dat ##防止下面操作把数值搞坏的一个备份
  dat=dat_back  ##表达矩阵数据
  dat[1:4,1:4]
  dat=t(dat)
  dat=as.data.frame(dat) ##转换为数据框
  dat=cbind(dat,group_list ) ##cbind()合并列（横向追加）;添加分组信息
  dat[1:4,1:4]
  ## 表达矩阵可以随心所欲的取行列，基础知识需要打牢。
  dat[1:4,12197:12199]
  dat[,ncol(dat)] #ncol()列，返回列长值
  table(dat$group_list)
  library("FactoMineR")
  library("factoextra") 
  # The variable group_list (index = ) is removed
  # before PCA analysis
  ## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #'-'表示“非”
  fviz_pca_ind(dat.pca,repel =T,
               geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
               col.ind = dat$group_list, # color by groups 颜色组
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses 集中成椭圆
               legend.title = "Groups"
  )
  ## 事实上还是有很多基因dropout非常严重。
  ggsave('all_cells_PCA.png')
  
  
  
  
}

