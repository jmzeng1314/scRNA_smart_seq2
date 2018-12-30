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
load(file = '../input.Rdata')
a[1:4,1:4]
head(df)

dat[1:4,1:4]
exprSet=dat

mean_per_gene <- apply(exprSet, 1, mean, na.rm = TRUE)
sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE)
mad_perl_gene <-   apply(exprSet, 1, mad, na.rm = TRUE)
cv_per_gene <- data.frame(mean = mean_per_gene,
  sd = sd_per_gene,
  mad=mad_perl_gene,
  cv = sd_per_gene/mean_per_gene)
rownames(cv_per_gene) <- rownames(exprSet)
head(cv_per_gene)
# pairs(cv_per_gene)
with(cv_per_gene,plot(log10(mean),log10(cv)))



