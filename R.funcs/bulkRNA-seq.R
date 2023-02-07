

# bulk PCA
# path: EllyLab/human/bulkRNA/bulk%20RNA-seq.ipynb
library(DESeq2)
library(ggplot2)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ group)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

options(repr.plot.width=8, repr.plot.height=5)
plotPCA(vsd, intgroup=c("group")) +
  theme_bw() +
  theme(axis.text.x  = element_text(face="plain", angle=90, size = 14, color = "black"), # , vjust=0.6
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 14))

