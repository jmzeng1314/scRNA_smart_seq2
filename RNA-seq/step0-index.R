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
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor

# http://www.bio-info-trainee.com/3727.html 周末生物信息学培训班准备工作


options()$repos  ## 查看使用install.packages安装时的默认镜像
options()$BioC_mirror ##查看使用bioconductor的默认镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") ##指定镜像，这个是中国科技大学镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##指定install.packages安装镜像，这个是清华镜像
options()$repos 
options()$BioC_mirror

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") ##判断是否存在BiocManager包，不存在的话安装

library(BiocManager)


BiocManager::install(c( 'scran'),ask = F,update = F)
 
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene",ask = F,update = F)
BiocManager::install("org.Mm.eg.db",ask = F,update = F)
BiocManager::install("genefu",ask = F,update = F)
BiocManager::install("org.Hs.eg.db",ask = F,update = F)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene",ask = F,update = F)
install.packages("ggfortify")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("ggfortify")
 
##  里面代码注释较好，请务必花 4个小时以上时间理解。
if(T){
  a=read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz',
               header = T ,sep = '\t')  ##把表达矩阵文件载入R，header=T :保留文件头部信息，seq='\t'以tap为分隔符
  # 每次都要检测数据
  a[1:6,1:4] #对于a矩阵取第1~6行，第1~4列
  ## 读取RNA-seq的 counts 定量结果，表达矩阵需要进行简单的过滤
  dat=a[apply(a,1, function(x) sum(x>1) > floor(ncol(a)/50)),] 
  #筛选表达量合格的行(基因), 列（细胞）数不变
  
  ##ncol()返回矩阵的列数值；floor()四舍五入取整数；
  ##function() 定义一个函数；sum()求和
  #上面的apply()指令代表对矩阵a进行行计算，判断每行表达量>1的样本总个数，并筛选出细胞表达量合格的基因（行）
  #第一个参数是指要参与计算的矩阵——a
  #第二个参数是指按行计算还是按列计算，1——表示按行计算，2——按列计算；
  #第三个参数是指具体的运算参数,定义一个函数x（即表达量为x）
  
  #对每行中x>1的列（即样本数）求和，即得出的是每行中表达量大于1的样本数，
  #然后再筛选出大于floor(ncol(a)/50)的行，这样的行（基因）的细胞表达量才算合格
  #因为2%的细胞有表达量，所以对于768个细胞样本，每个行(基因)在细胞中的表达至少要有15.36（约等于15）个样本表达才算合格
  ## 2 % 的细胞有表达量
  
  dat[1:4,1:4]
  sum(dat[,3])
  #  0610007P14Rik in  SS2_15_0048_A5 
  log2(18*1000000/sum(dat[,3])+1)
  ## 18 -- > 6.459884  ## SS2_15_0048_A5
  dat=log2(edgeR::cpm(dat)+1) 
  ##归一化的一种选择，这里是CPM(count-per-million，每百万碱基中每个转录本的count值)
  ###CPM只对read count相对总reads数做了数量的均一化，去除文库大小差异。
 
  dat[1:4,1:4] 
  
  ## 下面介绍一下 dist 函数
  x=1:10
  y=2*x
  z=rnorm(10)
  tmp=data.frame(x,y,z)
  dist(tmp)
  head(tmp)
  dist(t(tmp))
  cor(tmp)
  dist(t(scale(tmp)))
  # 可以看到dist函数计算样本直接距离和cor函数计算样本直接相关性，是完全不同的概念。虽然我都没有调它们两个函数的默认的参数。
  # 
  # 总结：
  # 
  # - dist函数计算行与行（样本）之间的距离
  # - cor函数计算列与列（样本）之间的相关性
  # - scale函数默认对每一列（样本）内部归一化
  # - 计算dist之前，最好是对每一个样本（列）进行scale一下
  # 
  
  
  #层次聚类，因为近 800细胞，非常耗时。
  hc=hclust(dist(t(dat))) ##样本间层次聚类
  # 原始表达矩阵转置后，细胞在行，所以计算的是细胞与细胞之间的距离。
  ## statquest 有详细讲解背后的统计学原理。
  class(hc)
  ?plot.hclust
  ## 查看说明书。
  plot(hc,labels = FALSE)
  
  #t:矩阵转置，行转列，列转行
  #分类时常常需要估算不同样本之间的相似性(Similarity Measurement)
  # 这时通常采用的方法就是计算样本间”距离”(Distance)。
  
  #dist函数是R语言计算距离的主要函数。dist函数可以计算行与行两两间的距离。
  # 所以之前的矩阵里面行是基因，转置后行是样本，因为我们要计算样本与样本之间的距离。
  # dist()函数计算变量间距离
  #hclust函数用来层次聚类
  
  clus = cutree(hc, 4) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
  group_list= as.factor(clus) ##转换为因子属性
  table(group_list) ##统计频数
  
  #提取批次信息
  colnames(dat) #取列名
  library(stringr)
  plate=str_split(colnames(dat),'_',simplify = T)[,3] #取列名，以'_'号分割，提取第三列。
  #str_split()函数可以分割字符串
  table(plate)
  
  n_g = apply(a,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）
  # 这里我们定义， reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
  # 而且大部分单细胞转录组技术很烂，通常超过75%的基因都没办法检测到。

  df=data.frame(g=group_list,plate=plate,n_g=n_g) #新建数据框(细胞的属性信息)
  
  ##(样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)
  
  df$all='all' #添加列，列名为"all"，没事意思，就是后面有需要
  metadata=df
  save(a,dat,df,file = '../input.Rdata') #保存a,dat,df这变量到上级目录的input.Rdata
  # 因为另外一个项目也需要使用这个数据集，所以保存到了上级目录。
} 

if(T){
  a=read.table('../GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz',
               header = T ,sep = '\t')  ##把表达矩阵文件载入R，header=T :保留文件头部信息，seq='\t'以tap为分隔符
  # 每次都要检测数据
  a[1:6,1:4] #对于a矩阵取第1~6行，第1~4列
  ## 读取RNA-seq的 counts 定量结果，表达矩阵需要进行简单的过滤
  dat=a[apply(a,1, function(x) sum(x>0) > floor(ncol(a)/50)),] #筛选表达量合格的行,列数不变
 
  dat[1:4,1:4]   
  #层次聚类
  hc=hclust(dist(t(log(dat+0.1)))) ##样本间层次聚类
  # 如果是基因聚类，可以选择 wgcna 等算法 
  ## statquest 
  plot(hc,labels = F)
  clus = cutree(hc, 4) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
  group_list= as.factor(clus) ##转换为因子属性
  table(group_list) ##统计频数
  
  #提取批次信息
  colnames(dat) #取列名
  library(stringr)
  plate=str_split(colnames(dat),'_',simplify = T)[,3] #取列名，以'_'号分割，提取第三列。
  #str_split()函数可以分割字符串
  table(plate)
  
  n_g = apply(a,2,function(x) sum(x>0)) #统计每个样本有表达的有多少行（基因）
  # 这里我们定义， reads数量大于1的那些基因为有表达，一般来说单细胞转录组过半数的基因是不会表达的。
  # 而且大部分单细胞转录组技术很烂，通常超过75%的基因都没办法检测到。
  
  df=data.frame(g=group_list,plate=plate,n_g=n_g) #新建数据框(细胞的属性信息)
  
  ##(样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)
  
  df$all='all' #添加列，列名为"all"，没事意思，就是后面有需要
  metadata=df
  save(dat,metadata,file = '../input_rpkm.Rdata') #保存a,dat,df这变量到上级目录的input.Rdata
  # 因为另外一个项目也需要使用这个数据集，所以保存到了上级目录。
} 


load(file = '../input.Rdata') ##从上级目录载入input.Rdata

## 每次载入以前的变量，都是可以简单检查一下。
a[1:4,1:4] 
dat[1:4,1:4] 
head(df)

load(file = '../input_rpkm.Rdata') ##从上级目录载入 input_rpkm.Rdata

## 每次载入以前的变量，都是可以简单检查一下。
#a[1:4,1:4] 
dat[1:4,1:4] 
head(metadata)



