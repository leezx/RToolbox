# commonly used data

```
hs.TFs <- read.csv("https://github.com/leezx/RToolbox/raw/master/data/Homo_sapiens_TF.txt", sep = "\t")
hs.coTFs <- read.csv("https://github.com/leezx/RToolbox/raw/master/data/Homo_sapiens_TF_cofactors.txt", sep = "\t")
hs.TFs <- unique(c(hs.TFs$Symbol, hs.TFs$Symbol))
length(hs.TFs)
 
hs.epi.reg <- read.csv("https://github.com/leezx/RToolbox/raw/master/data/KAT6A-chromatin-regulators.csv", sep = ",")
hs.epi.reg <- unique(hs.epi.reg$Gene)
hs.epi.reg <- hs.epi.reg[!grepl("NonTargetingControlGuide",hs.epi.reg)]
length(hs.epi.reg)
 
gplots::venn(list(hs.TFs=hs.TFs, hs.epi.reg=hs.epi.reg))
```

## all mouse genes
- gene.anno.mm10.3.0.0.csv - 提取自10x的gtf文件

## all human genes
- gene.anno.GRCh38.3.0.0.csv - 提取自10x的gtf文件
- gene.anno.hg19.3.0.0.csv - 提取自10x的gtf文件

## all mouse and human TFs
- method 1: [The Human Transcription Factors](http://humantfs.ccbr.utoronto.ca/), This website contains the catalog of 1639 known and likely human TFs and their motifs (version 1.01).
- method 2: http://211.67.31.242/HumanTFDB/#!/download
- method 3: Transcription factors were identified using GO terms transcription factor activity (GO:0000989), and regulation of transcription, DNA dependent (GO:0006355).
- 文章的数据，pdf转Excel，mmc13-Transcription Factors.csv
- 问题：human和mouse严格意义mapping的TF只有1000个，各自有600个左右的unique TF，查一下为什么map不上！！！

```r
# human TF v1
hTFs1 <- read.csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv", header = T)
hTFs1 <- hTFs1[!duplicated(hTFs1$HGNC.symbol),]
dim(hTFs1)
# human TF v2
hTFs2 <- read.csv("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF", header = T, sep="\t")
hTFs2 <- hTFs2[!duplicated(hTFs2$Symbol),]
dim(hTFs2)
# human TF coFactors
hTF_coFs <- read.csv("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF_cofactors", header = T, sep="\t")
hTF_coFs <- hTF_coFs[!duplicated(hTF_coFs$Symbol),]
dim(hTF_coFs)
```

## all mouse and human signaling factors
- Signaling factors were identified from the terms growth factor (Panther MF00019) and secreted (SP_PIR_KEYWORDS). Semaphorins and Slit proteins were identified manually and added to the lists.
- 文章的数据，pdf转Excel，mmc14-Signaling ligands.csv

## all mouse and human receptors
- Receptors were identified from the annotation terms receptor (SP_PIR_KEYWORDS) and receptor activity (GO:0004872). The combined lists of genes were manually screened to remove wrongly annotated genes.
- 文章的数据，pdf转Excel，mmc15-Signaling Receptors.csv


参考：
- [hg19 | GRCh38 | mm10 | 基因组 | 功能区域 | 位置提取](https://www.cnblogs.com/leezx/p/11889925.html)
- Transcription and Signaling Regulators in Developing Neuronal Subtypes of Mouse and Human Enteric Nervous System
