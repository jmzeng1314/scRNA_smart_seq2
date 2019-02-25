# 单细胞转录组数据处理视频课程

课程说明在：https://mp.weixin.qq.com/s/AV2uTbsvJGBRq_zv7yDmNg

### 表达矩阵获取

这里的例子是2018年12月的NC文章：[Spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell RNA sequencing](https://www.nature.com/articles/s41467-018-07582-3)  使用成熟的单细胞转录组( Smart-seq2 )手段探索了癌相关的成纤维细胞 CAFs的功能和空间异质性。

在文章搜索到作者上传的数据的GEO链接： [GSE111229](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111229) 就可以找到作者处理好的表达矩阵（counts和RPKM格式的都有）

| **Supplementary file**                                       | **Size** | **Download**                                                 | **File type/resource** |
| ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ---------------------- |
| GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz | 5.3 Mb   | [(ftp)](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111229/suppl/GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz)[(http)](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111229&format=file&file=GSE111229%5FMammary%5FTumor%5Ffibroblasts%5F768samples%5FrawCounts%2Etxt%2Egz) | TXT                    |
| GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz | 23.8 Mb  | [(ftp)](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111229/suppl/GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz)[(http)](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111229&format=file&file=GSE111229%5FMammary%5FTumor%5Ffibroblasts%5F768samples%5FrpkmNormalized%2Etxt%2Egz) | TXT                    |

如果想下载作者测序的原始测序数据：[SRP133642](https://www.ncbi.nlm.nih.gov/sra?term=SRP133642)  来走一波RNA-seq上游流程就需要：

shell脚本处理RNA-seq数据上游分析全部代码在：[code](./shell.txt)

- 代码参考：https://www.jianshu.com/p/a84cd44bac67

- 视频教程见：https://www.bilibili.com/video/av28453557

### 转录组分析回顾

主要是考虑到完全复现这篇文章数据的全部处理过程，需要掌握linux，r，转录组，考虑到不少人会基础知识有点薄弱，所以通过引入`常规转录组数据分析`的演示来提醒大家巩固基础知识！

这一单元代码都在  `RNA-seq` 文件夹，进入打开后缀是 `Rproj` 的文件就会自动调用你系统的Rstudio软件，从而定位到项目。

- step0-index.R 
  - 读取 作者的counts文件，简单过滤，并且logCPM转换
- step1-check.R
  - 检测表达矩阵里面细胞的相关性，hclust结果，热图展现，PCA图展现
- step2-cv2.R
  - 检测基因的变异系数及基因表达量的均值的相关性，探索其它统计学指标
- step3-batch-PCA-tSNE.R
  - 查看表达量是否受批次效应影响，这里简短介绍了PCA和tSNE作用
- step4-gene-number.R
  - 查看不同分组条件下检测到的基因数量的多少分布
- step5-pam50.R
  - 生物学背景知识，对乳腺癌研究不感兴趣的请不要学习
- step6-cell-cycle.R
  - 生物学背景知识，最好是自行搜索了解
- step7-counts2rpkm.R
  - 完全的R代码技巧，能力不够着请不要浪费时间学习，非常辛苦的。

- step8-DEG.R
  - 根据文章

### 单细胞转录组数据分析

#### 首先学习3个R包



#### 然后尝试把这3个R包应用到该文章的数据

