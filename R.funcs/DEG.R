# source("https://github.com/leezx/bt2m/raw/main/notebooks/unpackaged-code/DEG.R")

#' A general function to identify DEGs between case and control
cluster_DEG_twoGroups <- function(seuratObj, cluster, group_by, ident.ref, ident.query, pooled = T, 
                                  assay = "RNA", slot = "data") {
    options(warn=-1)
    DEGs <- list()
    seuratObj[["cluster"]] <- seuratObj[[cluster]]
    seuratObj[["group_by"]] <- seuratObj[[group_by]]
    for (i in unique(seuratObj$cluster)) {
        message(sprintf("identifying DEGs in %s between %s (%s and %s)...", i, group_by, ident.query, ident.ref))
        tmp.seuratObj <- subset(seuratObj, subset = cluster == i)
        tmp.seuratObj@active.ident <- tmp.seuratObj$group_by
        #
        # tmp.DEGs <- FindMarkers(tmp.seuratObj, ident.1 = ident.query, ident.2 = ident.ref, only.pos = F,
        #                        min.pct = 0, min.diff.pct = "-Inf", logfc.threshold = 0, assay = assay)
        #
        tmp.DEGs <- presto:::wilcoxauc.Seurat(X = tmp.seuratObj, group_by = 'group_by',
                                                 assay = slot, seurat_assay = assay)
        # tmp.DEGs <- tmp.DEGs[1:(nrow(tmp.DEGs)/2),] # only take 1 group
        tmp.DEGs <- subset(tmp.DEGs, group == ident.query)
        rownames(tmp.DEGs) <- tmp.DEGs$feature
        # return(tmp.DEGs) # testing
        if(dim(tmp.DEGs)[1]<1) {
          message(sprintf("skipping %s, no DEGs found...", i))
          next
        }
        # DEGs[[i]] <- add.missing.DEGs(tmp.DEGs, rownames(tmp.seuratObj@assays$RNA@counts))
        DEGs[[i]] <- tmp.DEGs
        # fill emplty genes
        #tmp.empty <- rowSums(DEGs[[i]])==0
        #DEGs[[i]][tmp.empty,]$pval <- 1
        #DEGs[[i]][tmp.empty,]$padj <- 1
        # add gene column
        DEGs[[i]]$gene <- rownames(DEGs[[i]])
        # add log2FC and correlation
        cells_1 <- rownames(subset(seuratObj@meta.data, cluster==i & group_by==ident.ref))
        cells_2 <- rownames(subset(seuratObj@meta.data, cluster==i & group_by==ident.query))
        log2fc <- log2(rowMeans(seuratObj@assays$RNA@data[,cells_2])+1) - log2(rowMeans(seuratObj@assays$RNA@data[,cells_1])+1)
        corM <- cor(t(as.matrix(seuratObj@assays$RNA@data[,c(cells_1,cells_2)])), c(rep(0, length.out = length(cells_1)),
                                                                       rep(1, length.out = length(cells_2))
                                                                       ))
        corM[is.na(corM)] <- 0
        DEGs[[i]]$log2FC <- log2fc[DEGs[[i]]$gene]
        DEGs[[i]]$correlation <- as.data.frame(corM)[DEGs[[i]]$gene,]
        #
        # break
    }
    if (pooled) {
        i <- "pooled"
        tmp.DEGs <- presto:::wilcoxauc.Seurat(X = seuratObj, group_by = 'group_by',
                                                 assay = slot, seurat_assay = assay)
        # tmp.DEGs <- tmp.DEGs[1:(nrow(tmp.DEGs)/2),] # only take 1 group
        tmp.DEGs <- subset(tmp.DEGs, group == ident.query)
        rownames(tmp.DEGs) <- tmp.DEGs$feature
        DEGs[[i]] <- tmp.DEGs
        DEGs[[i]]$gene <- rownames(DEGs[[i]])
        cells_1 <- rownames(subset(seuratObj@meta.data, group_by==ident.ref))
        cells_2 <- rownames(subset(seuratObj@meta.data, group_by==ident.query))
        log2fc <- log2(rowMeans(seuratObj@assays$RNA@data[,cells_2])+1) - log2(rowMeans(seuratObj@assays$RNA@data[,cells_1])+1)
        corM <- cor(t(as.matrix(seuratObj@assays$RNA@data[,c(cells_1,cells_2)])), c(rep(0, length.out = length(cells_1)),
                                                                       rep(1, length.out = length(cells_2))
                                                                       ))
        corM[is.na(corM)] <- 0
        DEGs[[i]]$log2FC <- log2fc[DEGs[[i]]$gene]
        DEGs[[i]]$correlation <- as.data.frame(corM)[DEGs[[i]]$gene,]
    }
    return(DEGs)
    # check DEGs
    # HBSO.DEG.1by1.sig <- lapply(HBSO.DEG.1by1, function(x) {
    # y <- subset(x, padj<0.05 & abs(correlation)>0.1)
    # y[order(y$padj, decreasing = F),c("auc","pval","padj","logFC","correlation")]
    # options(repr.plot.width=10, repr.plot.height=5)
    # VlnPlot(HBSO.combined, features = c("CHCHD2", "C1QL4"), group.by = "subtype", split.by = "sample", pt.size = 0, combine = F)
})
}


scde_func <- function(count_matrix, cell_group, coreN = 3, fileName="scde.o.ifm.Rdata", providemodel=F, pval=1){
  options(stringsAsFactors = FALSE)
  library(scde)
  #data(es.mef.small)
  #sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
  #names(sg) <- colnames(es.mef.small)
  #table(sg)
  cd <- clean.counts(count_matrix, min.lib.size=1000, min.reads = 5, min.detected = 1)
  sg <- cell_group
  if (providemodel){
    load(fileName)
  } else {
    o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = coreN, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
    save(o.ifm, file = fileName)
  }
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
  # estimate gene expression prior
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
  groups <- cell_group
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  coreN, verbose  =  1)
  p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
  p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
  significant.genes <- which(p.values.adj<=pval)
  ord <- order(p.values.adj[significant.genes]) # order by p-value
  de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
  colnames(de) <- c("Lower_bound","log2FC","Upper_bound","p_value")
  de$log2FC <- -(de$log2FC)
  return(de)
}

## volcano plot
draw_vovalno_plot <- function(data, log2FC_thred=0.5, markedGene, my_title="") {
  colnames(data) <- c("p_value", "log2FC")
  data$sig <- "non-sig"
  data[data$p_value<0.05 & abs(data$log2FC)>=log2FC_thred,]$sig <- "sig"
  # see /Users/surgery/Project/HOME/1-projects/0.bulk_RNA-Seq/9-BACE2/BACE2.R
  library(ggrepel) #Avoid overlapping labels
  # mark_data <- data[!data$sig=="non-sig" | rownames(data) %in% whiteList,]
  mark_data <- data[markedGene,]
  if (dim(mark_data)[1] > 30) {
    mark_data <- mark_data[order(abs(mark_data$log2FC), decreasing = T),][1:30,]
  }
  mark_data$gene <- rownames(mark_data)
  #
  data$color <- data$sig
  data[data$log2FC>0 & data$sig=="sig",]$color <- "up"
  data[data$log2FC<0 & data$sig=="sig",]$color <- "down"
  #
  ggplot(data, aes(x=log2FC, y=-log10(p_value))) +
    geom_hline(aes(yintercept=1), colour="grey50", linetype="dashed", size=0.2) +
    geom_vline(aes(xintercept=0.5), colour="red", linetype="dashed", size=0.2) +
    geom_vline(aes(xintercept=-0.5), colour="blue", linetype="dashed", size=0.2) +
    geom_point(aes(color=color), stroke = 0.3, size=1) +
    scale_color_manual(values=c("blue", "grey50",  "red")) +
    labs(x = "Log2 fold change",y = "-Log10(P-value)", title = my_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(size=0.8, colour = "black")) +
    theme(axis.title =element_text(size = 11),axis.text =element_text(size = 10, color = "black"),
          plot.title =element_text(hjust = 0.5, size = 16)) +
    #scale_x_continuous(position="bottom", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    #scale_y_continuous(position="left", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    theme(legend.title=element_blank(), legend.key.size = unit(0.8, 'lines')) +
    geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic",
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  # for SAG
  # geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.15, point.padding = 0.4, segment.color = 'black', segment.size = 0.3, force = 10, max.iter = 3e3)
}

