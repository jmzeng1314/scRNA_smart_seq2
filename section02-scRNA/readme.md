# 学习使用各种单细胞R包来处理数据

根据大家对NGS数据处理上游分析的掌握或需求，自行选择是否学习linux，但是需要完全理解我们的单细胞转录组表达矩阵是如何得到的，以及其生物学意义，每个基因比对到的reads数量的counts矩阵，以及去除了每个细胞测序数据量（文库大小）差异后的 rpm 矩阵，以及去除了基因长度效应的 rpkm矩阵，以及最近比较流行的 tpm 矩阵。

表达矩阵是单细胞转录组课程的开始！

### 主要数据分析要点分类

完整表单见：https://omictools.com/single-cell-rna-seq-category  我还在生信技能树写过推文介绍如何爬去工具列表，并且制作成为思维导图

#### normalization

[Linnorm](https://doi.org/10.1093/nar/gkx828),NODES, SAMstrt, SCnorm, scran, DESeq and TMM

#### feature Selection

- Detecting highly variable genes
- correlated gene pairs
- cell cycle phase
- tissue specific gene signatures

#### Dimension Reduction

MDS,PCA,t-SNE

#### clustering 

- K-means clustering 

- Mixture models 

- Hierarchical clustering

#### DEG analysis methods

monocle,MAST,SCDE, BASiCS, NODES, SAMstrt, [Seurat](http://satijalab.org/seurat/) and DESeq2

#### [Pseudotime](https://github.com/agitter/single-cell-pseudotime)

- Monocle / Monocle 2 / Census
- Wanderlust / Cycler / Wishbone
- SCUBA
- [Slingshot](https://www.biorxiv.org/content/early/2017/04/19/128843)

由于课程时间限制，以及我们所介绍的文章的数据限制，这里只能挑选`最出名的3个R包`来介绍， 它们这些R包或多或少涵盖了上面提到的部分分析内容。

值得注意的是：这里并不是说其它R包就不重要， 其实我在单细胞天地公众号也介绍过不少实用R包，请自行搜索学习，比如 `SC3, pcaReduce,SINCERA,M3Drop` 

学习下面的R包，需要掌握一些对象

### 关于测试数据

这里我们选择的是scRNAseq R包中的数据集

这个包内置的是 Pollen et al. 2014 数据集，人类单细胞细胞，分成`4`类，分别是 pluripotent stem cells 分化而成的 neural progenitor cells (“NPC”) ，还有 “GW16” and “GW21” ，“GW21+3” 这3种孕期细胞。 

首先我写了一个探索这个数据集的教程：[study_scRNAseq.html](http://bio-info-trainee.com/tmp/scRNA/study_scRNAseq.html)

### 关于seurat

学习seurat用法，当然是以官网为主，不过看英文笔记有挑战，简略带领大家一起学习咯： https://satijalab.org/seurat/get_started.html    主要学习：https://satijalab.org/seurat/pbmc3k_tutorial.html 

我这里主要演示使用 seurat包来处理 scRNAseq 这个R包内置的是 `Pollen et al. 2014` 单细胞转录组数据集 。

教程见：[study_seurat.html](http://bio-info-trainee.com/tmp/scRNA/study_seurat.html)

- `counts`矩阵进来后被包装为对象，方便操作。

- 然后一定要经过 `NormalizeData` 和 `ScaleData` 的操作

- 函数 `FindVariableGenes` 可以挑选适合进行下游分析的基因集。

- 函数 `RunPCA` 和 `RunTSNE` 进行降维

- 函数 `FindClusters` 直接就分群了，非常方便 函数 `FindAllMarkers` 可以对分群后各个亚群找标志基因。

- 函数 `FeaturePlot` 可以展示不同基因在所有细胞的表达量 

- 函数 `VlnPlot` 可以展示不同基因在不同分群的表达量差异情况 函数 `DoHeatmap` 可以选定基因集后绘制热图

### 关于scater

学习scater用法，当然是以官网为主，不过看英文笔记有挑战，简略带领大家一起学习咯：  https://bioconductor.org/packages/release/bioc/html/scater.html 

值得提醒的是 `2017年 11 月` 这个 scater 包经过了重大变革，所以如果大家看到比较旧的教程需要注意一下，通常是无法成功的。

其GitHub的教程：http://hemberg-lab.github.io/scRNA.seq.course/

我这里主要演示使用 scater 包来处理 scRNAseq 这个R包内置的是 `Pollen et al. 2014` 单细胞转录组数据集 。

教程见：[study_scater.html](http://bio-info-trainee.com/tmp/scRNA/study_scater.html)

### 关于monocle

学习monocle用法，当然是以官网为主，不过看英文笔记有挑战，简略带领大家一起学习咯：http://cole-trapnell-lab.github.io/monocle-release/monocle3/ 

我这里主要演示使用 monocle 包来处理 scRNAseq 这个R包内置的是 `Pollen et al. 2014` 单细胞转录组数据集 。

教程见：[study_monocle.html](http://bio-info-trainee.com/tmp/scRNA/study_monocle.html)

### 最后把学习的4个R包应用到文章的数据



### 扩展阅读

首先需要了解`bioconductor`：https://bioconductor.github.io/BiocWorkshops/

然后需要了解`scRNA课程`：https://hemberg-lab.github.io/scRNA.seq.course/index.html

另外一个还在持续制作中的课程：https://osca.bioconductor.org/ 

