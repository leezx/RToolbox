
#' prepare the annotation for the pheatmap
#' @param sce sce object
#' @param FeatureList FeatureList
#' @param sortBy sortBy
#' @param colors colors
#'
#' @return anno list
#' @export
#' @examples
#' group.order <- c("Control", "S-HSCR", "L-HSCR")
#' FeatureList <- list("Sample"=sample.order[1:8], "Group"=group.order)
#' anno <- pre.pheatmap.anno.general(sce_HSCR_c1, FeatureList, sortBy = "Sample")
#'
pre.pheatmap.anno.general <- function(sce=sce, FeatureList=FeatureList, sortBy="impute_cluster",
                                      colors=list("1"=brewer.pal(8,"Set2"), "2"=brewer.pal(12,"Set3"))) {
  # group.order=group.order, sample.order=sample.order, Cluster="impute_cluster", Sample="cellGroup"
  anno <- list()
  annotation_colors <- FeatureList
  annotation_col <- as.data.frame(colData(sce)[,c("cellName", names(FeatureList))])
  for (j in 1:length(names(FeatureList))) {
    # print(j)
    i <- names(FeatureList)[j]
    # print(i)
    tmpNames <- FeatureList[[i]]
    # print(tmpNames)
    annotation_col[,i] <- factor(annotation_col[,i], levels = tmpNames)
    tmpcolor <- colors[[j]][1:length(tmpNames)] # mistake [] and [[]]
    if (length(tmpNames)>8) {
      tmpcolor <- brewer.pal(12,"Set3")[1:length(tmpNames)]
    }
    names(tmpcolor) <- tmpNames
    # print(annotation_colors)
    annotation_colors[[i]]  <- tmpcolor
  }
  annotation_col <- annotation_col[order(annotation_col[,sortBy], decreasing = F),]
  annotation_col$cellName <- NULL
  #######
  anno[["col"]] <- annotation_col
  anno[["colors"]] <- annotation_colors
  return(anno)
}

#' prepare the annotation for the pheatmap
#' @param sce sce object
#' @param group.order group.order
#' @param sample.order sample.order
#' @param Cluster the feature take as Cluster
#' @param Sample the feature take as Sample
#'
#' @return anno list
#' @export
#' @examples
#' options(repr.plot.width=8, repr.plot.height=7)
#' pheatmap(scale.data[sc3_marker$name, rownames(anno$col)], cluster_rows = F, cluster_cols = F, border_color = NA,
#'          color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(100),
#'          show_colnames=F, show_rownames=F, annotation_col = anno$col, annotation_colors = anno$colors)
#'
pre.pheatmap.anno <- function(sce=sce, group.order=group.order, sample.order=sample.order,
                              Cluster="impute_cluster", Sample="cellGroup") {
  anno <- list()
  #######
  annotation_col <- as.data.frame(colData(sce)[,c(Cluster, Sample)])
  colnames(annotation_col) <- c("Cluster", "Sample")
  annotation_col$Cluster <- factor(annotation_col$Cluster, levels = group.order)
  annotation_col$Sample <- factor(annotation_col$Sample, levels = sample.order)
  annotation_col <- annotation_col[order(annotation_col$Cluster, decreasing = F),]
  #######
  annotation_colors <- list()
  annotation_colors[["Cluster"]] <- brewer.pal(8,"Set2")[1:length(levels(annotation_col$Cluster))]
  annotation_colors[["Sample"]] <- brewer.pal(12,"Set3")[1:length(levels(annotation_col$Sample))]
  names(annotation_colors$Cluster) <- group.order
  names(annotation_colors$Sample) <- sample.order
  #######
  anno[["col"]] <- annotation_col
  anno[["colors"]] <- annotation_colors
  return(anno)
}

