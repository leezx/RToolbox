
# annotate peaks to TSS and genes
annotatePeak.zx <- function(peak.df.gr=norm.count, species="hs", TSS.range=1000) {
    #
    library(ChIPseeker)
    library(clusterProfiler)
    options(warn=-1)
    #
    if(class(peak.df.gr)=="data.frame") {
        tmp.gr <- GRanges(norm.count[,1:3])
        tmp.df <- norm.count
    } else if (class(peak.df.gr)=="GRanges") {
        tmp.gr <- peak.df.gr
        tmp.df <- as.data.frame(peak.df.gr)
    } else {
        message("Input is not data.frame nor GRanges, please double check!")
    }
    if (species=="hs") {
        # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        annoDb <- "org.Hs.eg.db"
    } else if (species=="mm") {
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        annoDb <- "org.Mm.eg.db"
    }
    #
    tmp.peaks.anno <- ChIPseeker::annotatePeak(tmp.gr, tssRegion=c(-TSS.range, TSS.range),
                         TxDb=txdb, annoDb=annoDb) # org.Hs.eg.db
    #
    tmp.df$geneName <- as.data.frame(tmp.peaks.anno@anno)[rownames(tmp.df),]$SYMBOL
    tmp.df$distanceToTSS <- as.data.frame(tmp.peaks.anno@anno)[rownames(tmp.df),]$distanceToTSS
    tmp.df <- subset(tmp.df, !is.na(distanceToTSS))
    tmp.df$type <- "Enhancer"
    tmp.df[abs(tmp.df$distanceToTSS)<=TSS.range,]$type <- "Promoter"
    return(tmp.df)
}

