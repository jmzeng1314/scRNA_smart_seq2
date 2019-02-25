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

# 在单细胞测序中起始RNA含量越低，技术重复样本的基因表达相关性也越低，
# 因此生成标准化的基因表达谱之后，非常重要的一步是评估技术噪音（technical variability）。
# 比较常见的一种方法是计算基因表达值的变异系数的平方（CV^2）

#标准差与平均数的比值称为变异系数，记为C.V(Coefficient of Variance)。 
# 变异系数又称“标准差率”，是衡量资料中各观测值变异程度的另一个统计量。 
#当进行两个或多个资料变异程度的比较时，如果度量单位与平均数相同，可以直接利用标准差来比较。

# 平均绝对误差（Mean Absolute Deviation），又叫平均绝对离差。
# 它是是所有单个观测值与算术平均值的偏差的绝对值的平均。



rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '../input_rpkm.Rdata')
#a[1:4,1:4]
head(df)
## 载入第0步准备好的表达矩阵，及细胞的一些属性（hclust分群，plate批次，检测到的基因数量）
#注意 变量a是原始的counts矩阵，变量 dat是log2CPM后的表达量矩阵。

dat[1:4,1:4]
exprSet=dat

mean_per_gene <- apply(exprSet, 1, mean, na.rm = TRUE) #对表达矩阵每行求均值
sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE) #对表达矩阵每行求标准差
mad_per_gene <-   apply(exprSet, 1, mad, na.rm = TRUE) #对表达矩阵每行求绝对中位差
## 同样的apply函数，多次出现，请务必学透它！！！

# 构造一个数据框来存放结果。
cv_per_gene <- data.frame(mean = mean_per_gene,
  sd = sd_per_gene,
  mad=mad_per_gene,
  cv = sd_per_gene/mean_per_gene)
rownames(cv_per_gene) <- rownames(exprSet)
head(cv_per_gene)
# pairs(cv_per_gene)
with(cv_per_gene,plot(log10(mean),log10(cv)))
with(cv_per_gene,plot(log10(mean),log10(cv^2)))
cv_per_gene$log10cv2=log10(cv_per_gene$cv^2)
cv_per_gene$log10mean=log10(cv_per_gene$mean)
library(ggpubr)
cv_per_gene=cv_per_gene[cv_per_gene$log10mean < 4, ]
cv_per_gene=cv_per_gene[cv_per_gene$log10mean > 0, ]
ggscatter(cv_per_gene, x = 'log10mean', y = 'log10cv2',
          color = "black", shape = 16, size = 1, # Points color, shape and size
          xlab = 'log10(mean)RPKM', ylab = "log10(cv^2)",
          add = "loess", #添加拟合曲线
          add.params = list(color = "red",fill = "lightgray"),
          cor.coeff.args = list(method = "spearman"), 
          label.x = 3,label.sep = "\n",
          dot.size = 2,
          ylim=c(-0.5, 3),
          xlim=c(0,4) 
)



