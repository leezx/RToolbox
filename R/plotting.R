



#' plotting smooth curve of different modules across pseudotime (in multiple samples)
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param moduleDf the marker genes (modules) dataframe
#' @param sampleCol sample column (for comparison)
#' @return a ggplot object
#' @export
#' @examples
#' load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Aug_30/analysis_4/Ctrl_YFP_ENCC.Rdata")
#' exprM <- seuset@scale.data; annoDf <- all_tsne; moduleDf <- markers
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' result <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#' annoDf <- result
#' plotting_modules_smooth_curve_across_pseudotime()
plotting_modules_smooth_curve_across_pseudotime <- function(exprM, annoDf, moduleDf, sampleCol=NULL) {
  # exprM <- exprM[,rownames(annoDf)]
  all_modules_exprM <- exprM[markers$gene,]
  module_mean_expr <- data.frame()
  for (cluster in unique(annoDf$cluster)) {
    print(cluster)
    temp_genes <- markers[markers$cluster==cluster,]$gene
    temp_modules_exprM_mean <- apply(all_modules_exprM[temp_genes,], 2, mean)
    module_mean_expr <- rbind(module_mean_expr, as.data.frame(t(temp_modules_exprM_mean)))
  }
  rownames(module_mean_expr) <- unique(annoDf$cluster)
  module_mean_expr <- as.data.frame(t(module_mean_expr))
  module_mean_expr$pseudotime <- annoDf[rownames(module_mean_expr),]$pseudotime
  #library(reshape)
  module_mean_expr_melt <- melt(module_mean_expr, id.vars=c("pseudotime"))
  # smooth line
  ggplot(data=module_mean_expr_melt, aes(x=pseudotime, y=value)) +
    # geom_point(size=0.1, alpha=0.1) +
    facet_wrap( ~ variable, ncol=1) +
    labs(x = "Pseudotime", y = "Log2(expr+1)") +
    geom_smooth(method = 'loess', se=T, alpha=0.2, size=0.5, weight=1, span = 0.1)
  # end
}
