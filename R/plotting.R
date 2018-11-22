



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
#' exprM <- seuset@scale.data; annoDf <- all_tsne
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' new_annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#' annoDf <- new_annoDf
#' moduleDf <- markers
#' plotting_modules_smooth_curve_across_pseudotime(exprM, annoDf, moduleDf)
#'
#' exprM <- exprM_all; annoDf <- new_annoDf; moduleDf <- markers; annoDf$cellGroup <- substr(rownames(annoDf), 1, 4)
#' 1. seuset_3 <- SubsetData(object = seuset, ident.use = c("Cluster1", "Cluster2", "Cluster6"))
#' exprM <- seuset_3@scale.data
#' all_tsne <- TSNEPlot(object = seuset_3)$data[[1]]
#' rownames(all_tsne) <- colnames(seuset_3@data)
#' all_tsne$cluster <- as.character(seuset_3@ident[rownames(all_tsne)])
#' all_tsne$cellGroup <- substr(rownames(all_tsne), 1, 4)
#' annoDf <- all_tsne
#' startCluster <- "Cluster1"; endCluster <- "Cluster6"
#' annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#'
#' 2. exprM <- seuset_2@scale.data
#' all_tsne <- TSNEPlot(object = seuset_2)$data[[1]]
#' rownames(all_tsne) <- colnames(seuset_2@data)
#' all_tsne$cluster <- as.character(seuset_2@ident[rownames(all_tsne)])
#' all_tsne$cellGroup <- substr(rownames(all_tsne), 1, 4)
#' annoDf <- all_tsne
#' startCluster <- "Cluster1"; endCluster <- "Cluster3"
#' annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#'
plotting_modules_smooth_curve_across_pseudotime <- function(exprM, annoDf, moduleList) {
  # detect boundary
  icpt <- list()
  for (i in 1:length(unique(annoDf$cluster))){
    if (i == length(unique(annoDf$cluster))) {break}
    temp_c <- unique(annoDf$cluster)[i]
    last_c <- unique(annoDf$cluster)[i+1]
    temp_p <- max(annoDf[annoDf$cluster==temp_c,]$pseudotime)
    last_p <- min(annoDf[annoDf$cluster==last_c,]$pseudotime)
    #print(temp_p)
    #print(last_p)
    print((temp_p+last_p)/2)
    icpt[[i]] <- (temp_p+last_p)/2
  }
  # exprM <- exprM[,rownames(annoDf)]
  all_modules_exprM <- exprM[unlist(moduleList, use.names=FALSE),]
  module_mean_expr <- data.frame()
  for (i in 1:length(moduleList)) {
    temp_genes <- moduleList[[i]]
    temp_modules_exprM_mean <- apply(all_modules_exprM[temp_genes,], 2, median) # median
    module_mean_expr <- rbind(module_mean_expr, as.data.frame(t(temp_modules_exprM_mean)))
  }
  rownames(module_mean_expr) <- names(moduleList)
  module_mean_expr <- as.data.frame(t(module_mean_expr))
  module_mean_expr$pseudotime <- annoDf[rownames(module_mean_expr),]$pseudotime
  #
  #module_mean_expr$cellGroup <- annoDf[rownames(module_mean_expr),]$cellGroup
  module_mean_expr$gender <- annoDf[rownames(module_mean_expr),]$gender
  module_mean_expr2 <- module_mean_expr
  #module_mean_expr2$cellGroup <- "mean"
  module_mean_expr2$gender <- "mean"
  #
  module_mean_expr_melt <- melt(module_mean_expr, id.vars=c("pseudotime","gender"))
  module_mean_expr_melt2 <- melt(module_mean_expr2, id.vars=c("pseudotime","gender"))
  #module_mean_expr_melt <- melt(module_mean_expr, id.vars=c("pseudotime","cellGroup"))
  #module_mean_expr_melt2 <- melt(module_mean_expr2, id.vars=c("pseudotime","cellGroup"))
  #
  module_mean_expr_melt3 <- rbind(module_mean_expr_melt, module_mean_expr_melt2)
  #subset_cells <- sample(1:dim(module_mean_expr_melt3)[1], 8000)
  # smooth line
  # subset1 <- subset(module_mean_expr_melt, cellGroup%in%c("ctrl", "post", "negt"))
  # subset2 <- subset(module_mean_expr_melt, cellGroup%in%c("ctrl", "kif7"))
  ggplot(data=module_mean_expr_melt3, aes(x=pseudotime, y=value)) +
    # geom_point(size=0.1, alpha=0.1) +
    facet_wrap( ~ variable, ncol=1) +
    labs(x = "Pseudotime", y = "Log2(expr+1)") +
    geom_smooth(method = 'loess', se=T, alpha=0.2, size=0.5, weight=1, span = 0.4, aes(color=gender, fill=gender)) +
    geom_vline(aes(xintercept=icpt[[1]]), colour="gray50", linetype="dashed") +
    geom_vline(aes(xintercept=icpt[[2]]), colour="gray50", linetype="dashed") +
    geom_vline(aes(xintercept=icpt[[3]]), colour="gray50", linetype="dashed") +
    scale_color_manual(values=c("deepskyblue","red","gray50")) +
    scale_fill_manual(values=c("deepskyblue","red","gray50"))
  # end
  # , aes(color=cellGroup, fill=cellGroup)
  # save 300X600   gender: 400X600  big:500X700
  # plot(annoDf$pseudotime, factor(annoDf$cluster)) # 500X350
}

#' plotting smooth curve of different modules across pseudotime (in multiple samples)
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param moduleDf the marker genes (modules) dataframe
#' @param sampleCol sample column (for comparison)
#' @return a ggplot object
#' @export
#' @examples
#' HH_genes <- c("Shh","Dhh","Ptch1","Smo","Cdon","Gli1","Gli2","Gli3","Hhip","Spopl","Btrc","Cul1","Cul3","Sufu","Kif7")
#' HH_genes2 <- c( "Arrb2","Smurf1","Csnk1g3","Evc2","Csnk1g1","Ccnd2","Spop","Prkaca","Prkacb","Kif2a", "Gpr161","Csnk1g2","Csnk1a1","Fbxw11","Csnk1d","Csnk1e","Gsk3b","Boc","Evc","Bcl2","Arrb1","Gas1",  "Ccnd1","Smurf2")
#'
#' genes <- HH_genes
plotting_genes_smooth_curve_across_pseudotime <- function(exprM, annoDf, genes) {
  # exprM <- exprM[,rownames(annoDf)]
  all_genes_exprM <- exprM[genes,]
  gene_expr <- as.data.frame(t(all_genes_exprM))
  gene_expr$pseudotime <- annoDf[rownames(gene_expr),]$pseudotime
  gene_expr$cellGroup <- annoDf[rownames(gene_expr),]$cellGroup
  # module_mean_expr$gender <- annoDf[rownames(module_mean_expr),]$gender
  #
  genes_expr_melt <- melt(gene_expr, id.vars=c("pseudotime", "cellGroup"))
  # smooth line
  #subset1 <- subset(genes_expr_melt, cellGroup%in%c("post", "negt"))
  #subset2 <- subset(genes_expr_melt, cellGroup%in%c("ctrl", "kif7"))
  ggplot(data=genes_expr_melt, aes(x=pseudotime, y=value)) +
    scale_x_continuous(limits=c(0,1),breaks=round(seq(0,1,length.out=5),2)) +
    # geom_point(size=0.1, alpha=0.1) +
    facet_wrap( ~ variable, ncol=5) +
    labs(x = "Pseudotime", y = "Log2(expr+1)") +
    geom_smooth(method = 'loess', se=T, alpha=0.2, size=0.5, weight=1, span = 0.7, aes(color=cellGroup, fill=cellGroup)) +
    theme(axis.text.x  = element_text(face="plain", angle=30, size = 8, color = "black"))
  # end
}

#' plotting violin plot of markers across clusters
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param moduleDf the marker genes (modules) dataframe
#' @param sampleCol sample column (for comparison)
#' @return a ggplot object
#' @export
#' @examples
#' cols.use <- mycolors
#' plotting_violin_plot_of_markers_across_clusters_by_seurat(seuset, markers, topNum=3, cols.use=brewer.pal(8,"Set2"))
plotting_violin_plot_of_markers_across_clusters_by_seurat <- function(seuset, markers, topNum=3, cols.use=NULL) {
  topMarkers <- markers %>% group_by(cluster) %>% top_n(topNum, avg_logFC) #
  VlnPlot(object = seuset, features.plot = topMarkers$gene, nCol=topNum, point.size.use=0, cols.use=cols.use, size.x.use=0, size.y.use=0, size.title.use=12, single.legend=T, x.lab.rot=T)
}

