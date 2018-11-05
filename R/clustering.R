

#' clustering by seurat
#'
#' @param rawCountM raw expression matrix from 10x (can't apply to smart-seq data)
#' @param clusterNum cluster number
#' @param colors.use colors used for clustering (should be more than the cluster num)
#' @param genes genes wanted to used in staining plot (FeaturePlot)
#' @param nCol column set in staining plot (FeaturePlot)
#' @return a seurat object and some plots
#' @export
#' @examples
#' load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Aug_30/analysis_4/gbm.4.pos.neg.ctrl.vcl.Rdata")
#' exprM_ctrl <- exprs(gbm_ctrl)
#' mm10 <- read.table("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Reeson_mar/mm10.feature.info", row.names = 1, header = F)
#' rownames(exprM_ctrl) <- mm10[rownames(exprM_ctrl),]$V2
#' rawCountM <- exprM_ctrl; mycolors <- brewer.pal(8,"Set2")
#' ENSgene <- c("Top2a","Cdk1","Aurkb","Myc","Snai2","Sox9","Sox10","Foxd3","Phox2b","Isl1","Plp1","Tubb3","Tubb4a","Elavl4","Ret","Fabp7")
#' new_seuset <- clustering_by_seurat(rawCountM, clusterNum=5, colors.use=mycolors, genes=ENSgene)
clustering_by_seurat <- function(rawCountM, clusterNum=5, colors.use=NULL, genes=NULL, nCol = 4) {
  ## reduce dimension
  seuset <- CreateSeuratObject(raw.data = rawCountM, min.cells = 2, min.genes = 100)
  # VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 2)
  # GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")
  seuset <- FilterCells(object = seuset, subset.names = c("nUMI"))
  # FilterCells(object, subset.names, low.thresholds, high.thresholds,cells.use = NULL)
  seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)
  seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  # length(x = seuset@var.genes)
  seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI"))
  seuset <- RunPCA(object = seuset, pc.genes = seuset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  # PrintPCA(object = seuset, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  # VizPCA(object = seuset, pcs.use = 1:2)
  # PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
  # PCHeatmap(object = seuset, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  # seuset <- JackStraw(object = seuset, num.replicate = 100)
  # JackStrawPlot(object = seuset, PCs = 1:9)
  # PCElbowPlot(object = seuset)
  ###
  ## Binary search the best resolution for clustering
  minResolution <- 0.01; maxResolution <- 0.99
  seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = minResolution, print.output = 0, save.SNN = TRUE)
  # PrintFindClustersParams(object = seuset)
  min_clusterNum <- length(unique(seuset@ident))
  print(paste("minResolution: ",minResolution,"; min_clusterNum: ", min_clusterNum))
  #
  seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = maxResolution, print.output = 0, save.SNN = TRUE)
  # PrintFindClustersParams(object = seuset)
  max_clusterNum <- length(unique(seuset@ident))
  print(paste("maxResolution: ",maxResolution,"; max_clusterNum: ", max_clusterNum))
  #
  tempResolution <- (minResolution+maxResolution)/2
  seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = tempResolution, print.output = 0, save.SNN = TRUE)
  # PrintFindClustersParams(object = seuset)
  mid_clusterNum <- length(unique(seuset@ident))
  #
  while(mid_clusterNum != clusterNum) {
    if (mid_clusterNum < clusterNum) {
      minResolution <- tempResolution
      tempResolution <- (tempResolution+maxResolution)/2
      seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = tempResolution, print.output = 0, save.SNN = TRUE)
      # PrintFindClustersParams(object = seuset)
      mid_clusterNum <- length(unique(seuset@ident))
      print(paste("tempResolution: ",tempResolution,"; mid_clusterNum: ", mid_clusterNum))
      #
    } else if (mid_clusterNum > clusterNum) {
      maxResolution <- tempResolution
      tempResolution <- (tempResolution+minResolution)/2
      seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = tempResolution, print.output = 0, save.SNN = TRUE)
      # PrintFindClustersParams(object = seuset)
      mid_clusterNum <- length(unique(seuset@ident))
      print(paste("tempResolution: ",tempResolution,"; mid_clusterNum: ", mid_clusterNum))
    }
  }
  # adjustedRandIndex(colData(as.matrix(deng))[seuset@cell.names, ]$cell_type2, seuset@ident)
  seuset <- RunTSNE(object = seuset, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
  # TSNEPlot(object = seuset)
  ## Assigning cell type identity to clusters
  current.cluster.ids <- sort(unique(seuset@ident))
  new.cluster.ids <- paste("Cluster", as.integer(current.cluster.ids), sep = "")
  # new.cluster.ids <- factor(new.cluster.ids, levels = new.cluster.ids)
  seuset@ident <- plyr::mapvalues(x = seuset@ident, from = current.cluster.ids, to = new.cluster.ids)
  seuset@ident <- factor(seuset@ident, levels = new.cluster.ids)
  # additional plot
  if(is.null(colors.use)) {
    TSNEPlot(object = seuset, do.label = TRUE, pt.size = 2)
  } else {
    TSNEPlot(object = seuset, do.label = TRUE, pt.size = 2, colors.use=mycolors[1:length(current.cluster.ids)])
  }
  if(!is.null(genes)) {
    FeaturePlot(seuset, genes, cols.use = c("lightgrey", "blue"), nCol = nCol)
  }
  # return
  seuset
  # end
}

#' marker identification by seurat
#'
#' @param seuset a seurat object after clustering
#' @param topNum top X markers
#' @param plotHeatmap whether to plot heatmap or not (default yes)
#' @param group.order specify the cluster order (after pseudotime)
#' @return a marker dataframe and the heatmap
#' @export
#' @examples
#' markers <- marker_identification_by_seurat(seuset, topNum=10, plotHeatmap=T, group.order=NULL)
marker_identification_by_seurat <- function(seuset, topNum=10, plotHeatmap=T, group.order=NULL) {
  markers <- FindAllMarkers(object = seuset, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  if (!is.null(group.order)) {
    markers$cluster <- factor(markers$cluster, group.order)
    markers <- markers[order(markers$cluster),]
  }
  top10 <- markers %>% group_by(cluster) %>% top_n(topNum, avg_logFC)
  # seuset@ident <- factor(seuset@ident, levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5"))
  #if (plotHeatmap) {
  DoHeatmap(object = seuset, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, group.label.rot=T, group.order=group.order)
  #}
  # return
  markers
  # end
}
