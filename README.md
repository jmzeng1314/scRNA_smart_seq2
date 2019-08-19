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

这一单元代码都在  `section01-RNA-seq` 文件夹，进入打开后缀是 `Rproj` 的文件就会自动调用你系统的Rstudio软件，从而定位到项目。

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

见 section02-scRNA 文件夹代码，因为文章其实并没有使用主流单细胞转录组R包，这里仅仅是根据scRNAseq的示例数据来讲解，进入打开后缀是 `Rproj` 的文件就会自动调用你系统的Rstudio软件，从而定位到项目。

#### 然后尝试把这3个R包应用到该文章的数据

见 section03-for_paper 文件夹代码，进入打开后缀是 `Rproj` 的文件就会自动调用你系统的Rstudio软件，从而定位到项目。可以完全复现文章图表。 

#### 最后是公共数据库挖掘

在 section04-downstream 文件夹，进入打开后缀是 `Rproj` 的文件就会自动调用你系统的Rstudio软件，从而定位到项目。

### 学习笔记在单细胞天地持续更新

- [单细胞转录组学习笔记-1](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484464&idx=1&sn=cec394967fcd20a25c258bed892f2457&chksm=ea1f02b2dd688ba4c2bb0659f7f66349c49d876ec340c84763c29c7bf080dba5cfea90d4b68e&scene=21#wechat_redirect)
- [单细胞转录组学习笔记-2](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484484&idx=1&sn=eb28c69645f8dbfda6c8a4afaf931883&chksm=ea1f02c6dd688bd0ca7ab1caa3261d9c0f59dc9c83002ca3780bf88e11a94f0893bff3849900&scene=21#wechat_redirect)
- [单细胞转录组上游分析之shell回顾](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484498&idx=1&sn=45a7aa58d726b68dac1c80d4b686f386&chksm=ea1f02d0dd688bc64fd14e5f0370ec1889215fc5ef3837e34db84f1319f517598d3e69832fc2&scene=21#wechat_redirect)
- [获取Github代码包以及准备工作](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484551&idx=1&sn=f58786c646cad53047b28cdc4f9bcd9e&chksm=ea1f0205dd688b1326aa5309e506452aacaa9edadf46e44d8cc9547f1126abc553362ed96800&scene=21#wechat_redirect)
- [常说的表达矩阵，那得到之后呢？](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484603&idx=1&sn=6c13a90be81285a71b915c77caa30407&chksm=ea1f0239dd688b2fce040852e8e799d410517e0924c45bb4cd00f749cafbe3972d8def3c2ead&scene=21#wechat_redirect)
- [由表达矩阵看内部异质性](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484644&idx=1&sn=eb86415a435da0ef81c1c0ae3c2767fc&chksm=ea1f0266dd688b7030db7bb57782490f0f23d91c79bb9b89a3e6135657544dd131a9cc26a503&scene=21#wechat_redirect)
- [重复平均表达量和变异系数相关性散点图](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484658&idx=1&sn=c0c42921ad194abb76ebdd1bc800a7f7&chksm=ea1f0270dd688b66eb730b2459eb65b885fe85a2f7478c26f2aafc3b3ae6a5054f4de75d7c2c&scene=21#wechat_redirect)
- [聚类算法之PCA与tSNE](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484689&idx=1&sn=6689a065d287b0308354b00c0a93bb13&chksm=ea1f0393dd688a85faa81c61c8c2528290a00aa8e379da1f939f5d4f5cfa2127581228961835&scene=21#wechat_redirect)
- [统计细胞检测的基因数量](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484705&idx=1&sn=69b50b56286cd09bc1e9462b3218722d&chksm=ea1f03a3dd688ab57971f11b5e6433e0951cc68fff4de007d01f2ee3237413832420b7ad771b&scene=21#wechat_redirect)
- [乳腺癌领域之PAM50分类](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484714&idx=1&sn=c5a232af5b0be717b40af56f9f10f9e2&chksm=ea1f03a8dd688abe0c61699db9d9fb1248c5e749abc9354677662a8fa216f3c78ac029f0fd90&scene=21#wechat_redirect)
- [生物学背景知识之细胞周期推断](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484719&idx=1&sn=c2e28e616261078936e7d34505b03120&chksm=ea1f03addd688abb0e67f1416aef16167584de2b2231b242e1c612e6b873b02eccf1f9575d92&scene=21#wechat_redirect)
- [RPKM概念及计算方法](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484749&idx=1&sn=f734117fc9fcb604071b9b6df1739aeb&chksm=ea1f03cfdd688ad97a021f34d0b54fba7c253b7ea7c49cd00a6dc59ead3f2359b794b29576c7&scene=21#wechat_redirect)
- [差异分析及KEGG注释简介](http://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484838&idx=1&sn=2823e03598d88fd4d219988e8bd69397&chksm=ea1f0324dd688a328ad4fe11320815121461ca6b9452cf8e98aa069280931515f8638a776694&scene=21#wechat_redirect)
- 更新中