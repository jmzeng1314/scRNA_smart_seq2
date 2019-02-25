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

## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的细胞基因）
# 注意 变量a是原始的counts矩阵，变量 dat是logCPM后的表达量矩阵。

group_list=df$g
plate=df$plate
table(plate)

rownames(dat)=toupper(rownames(dat)) ##toupper()函数，把小写字符转换成大写字符
dat[1:4,1:4]

library(genefu)

if(T){
  ddata=t(dat)
  ddata[1:4,1:4]
  s=colnames(ddata);head(s);tail(s) ##把实验检测到的基因赋值给S
  library(org.Hs.eg.db) ##人类基因信息的包
  s2g=toTable(org.Hs.egSYMBOL)
  g=s2g[match(s,s2g$symbol),1];head(g) ##取出实验检测到的基因所对应的基因名
 
  # match（x, y）返回的是vector x中每个元素在vector y中对映的位置（positions in y），
  # 如果vector x中存在不在vector y中的元素，该元素处返回的是NA
  # probe Gene.symbol Gene.ID
  
  dannot=data.frame(probe=s,
                    "Gene.Symbol" =s, 
                    "EntrezGene.ID"=g)
  #s向量是实验检测基因的基因名字，g向量是标准基因ID
  # 这里s应该和g是一一对应的，制作一个数据框
  
  ddata=ddata[,!is.na(dannot$EntrezGene.ID)] #ID转换
  dim(ddata)
  #制作行为样本，列为实验检测基因（这里的剩下的实验检测基因都有标准基因ID对应）的矩阵。
  #即剔除无基因ID对应的列
  
  # !is.na去除dannot数据框EntrezGene.ID列为NA的行（去除NA值即去除没有标准基因ID对应的实验检测基因名））
  
  dannot=dannot[!is.na(dannot$EntrezGene.ID),] #去除有NA的行，即剔除无对应的基因
  head(dannot)
  ddata[1:4,1:4]
  library(genefu)
  # c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow")
  
  s<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                         annot=dannot,do.mapping=TRUE)
  table(s$subtype)
  tmp=as.data.frame(s$subtype)
  subtypes=as.character(s$subtype)
}
head(df)
df$subtypes=subtypes
table(df[,c(1,5)])

library(genefu)
pam50genes=pam50$centroids.map[c(1,3)]
pam50genes[pam50genes$probe=='CDCA1',1]='NUF2'
pam50genes[pam50genes$probe=='KNTC2',1]='NDC80'
pam50genes[pam50genes$probe=='ORC6L',1]='ORC6'
x=dat
dim(x)
x=x[pam50genes$probe[pam50genes$probe %in% rownames(x)] ,]
table(group_list)
tmp=data.frame(group=group_list,
               subtypes=subtypes)
rownames(tmp)=colnames(x)

library(pheatmap)


pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_raw.png') 


x=t(scale(t(x)))
x[x>1.6]=1.6
x[x< -1.6]= -1.6

pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_scale.png') 








