# commonly used data

## all mouse genes
- gene.anno.mm10.3.0.0.csv - 提取自10x的gtf文件

## all human genes
- gene.anno.GRCh38.3.0.0.csv - 提取自10x的gtf文件
- gene.anno.hg19.3.0.0.csv - 提取自10x的gtf文件

## all mouse and human TFs
- method 1: [The Human Transcription Factors](http://humantfs.ccbr.utoronto.ca/), This website contains the catalog of 1639 known and likely human TFs and their motifs (version 1.01).
- method 2: http://211.67.31.242/HumanTFDB/#!/download
- method 3: Transcription factors were identified using GO terms transcription factor activity (GO:0000989), and regulation of transcription, DNA dependent (GO:0006355).

```r
#
hTFs <- read.csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv", header = T)
#
hTFs <- read.csv("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF", header = T, sep="\t")
hTF_coFs <- read.csv("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF_cofactors", header = T, sep="\t")
```

## all mouse and human signaling factors
- Signaling factors were identified from the terms growth factor (Panther MF00019) and secreted (SP_PIR_KEYWORDS). Semaphorins and Slit proteins were identified manually and added to the lists.

## all mouse and human receptors
- Receptors were identified from the annotation terms receptor (SP_PIR_KEYWORDS) and receptor activity (GO:0004872). The combined lists of genes were manually screened to remove wrongly annotated genes.



参考：
- [hg19 | GRCh38 | mm10 | 基因组 | 功能区域 | 位置提取](https://www.cnblogs.com/leezx/p/11889925.html)
