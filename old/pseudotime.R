

#' pseudotime analysis based on clustering.
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param startCluster starting Cluster
#' @param endCluster ending Cluster
#' @return a new annoDf with pseudotime value for each cell
#' @export
#' @examples
#' load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Aug_30/analysis_4/Ctrl_YFP_ENCC.Rdata")
#' exprM <- seuset@scale.data
#' all_tsne <- TSNEPlot(object = seuset)$data[[1]]
#' rownames(all_tsne) <- colnames(seuset@data)
#' all_tsne$cluster <- as.character(seuset@ident[rownames(all_tsne)])
#' all_tsne$cellGroup <- substr(rownames(all_tsne), 1, 4)
#' annoDf <- all_tsne
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' new_annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
pseudotime_based_on_clustering <- function(exprM, annoDf, startCluster, endCluster, threads=3) {
  options(stringsAsFactors = F)
  print("1. Get the center of each cluster...")
  # 1. Get the center of each cluster
  centerM <- data.frame()
  for (cluster in unique(annoDf$cluster)) {
    cells <- rownames(annoDf[annoDf$cluster == cluster,])
    centerM <- rbind(centerM, t(data.frame(apply(exprM[,cells],1,mean))))
  }
  rownames(centerM) <- unique(annoDf$cluster)
  print("2. Sort each center by Insert sorting...")
  # 2. Sort each center by Insert sorting
  orderList <- c(startCluster, endCluster)
  objectList <- unique(annoDf$cluster)[!unique(annoDf$cluster) %in% orderList]
  library(parallelDist)
  distM <- parDist(x = as.matrix(centerM), method = "euclidean", threads=threads)
  distM <- as.matrix(distM)
  distM[is.na(distM)] <- 0
  rownames(distM) <- rownames(centerM)
  colnames(distM) <- rownames(centerM)
  for (cluster in objectList) {
    len_orderList <- length(orderList)
    temp_dis_df <- data.frame()
    for (i in 1:(len_orderList-1)) {
      temp_pairs <- orderList[i:(1+i)]
      temp_dis_df <- rbind(temp_dis_df, temp_pairs)
    }
    temp_close_cluster_list <- c()
    for (j in 1:dim(temp_dis_df)[1]) {
      temp_close_cluster <- distM[cluster,temp_dis_df[j,1]] + distM[cluster,temp_dis_df[j,2]]
      temp_close_cluster_list <- c(temp_close_cluster_list, temp_close_cluster)
    }
    temp_dis_df <- cbind(as.data.frame(temp_dis_df), temp_close_cluster_list)
    temp_dis_df <- temp_dis_df[order(temp_dis_df[,3], decreasing = F),]
    # print(temp_dis_df)
    c1 <- temp_dis_df[1,1]; c2 <- temp_dis_df[1,2]
    orderList <- c(orderList[1:which(c1 == orderList)], cluster, orderList[which(c2 == orderList):length(orderList)])
  }
  # 3. Output the order of cluster
  print(c("3. The order of the cluster is: ", paste(orderList,collapse=" -> ")))
  print("4. Assign the pseudotime for each cell...")
  # 4. Assign the pseudotime for each cell
  # exprM_merged <- cbind(exprM, t(centerM))
  library(pdist)
  distM_merged <- pdist(t(exprM), centerM)
  distM_merged <- as.matrix(distM_merged)
  #distM_merged <- parDist(x = as.matrix(t(exprM_merged)), method = "euclidean", threads=3)
  rownames(distM_merged) <- colnames(exprM)
  colnames(distM_merged) <- rownames(centerM)
  #
  ordered_cells <- c()
  for (i in 1:length(orderList)) {
    direction <- T
    tempCluster <- orderList[i]
    nearestCluster <- orderList[i+1]
    if (i == length(orderList)) {nearestCluster <- orderList[i-1]; direction <- F}
    # print(nearestCluster)
    temp_cells <- rownames(annoDf[annoDf$cluster==tempCluster,])
    temp_ordered_cells <- temp_cells[order(distM_merged[temp_cells,nearestCluster], decreasing = direction)]
    ordered_cells <- c(ordered_cells, temp_ordered_cells)
  }
  annoDf <- annoDf[ordered_cells,]
  annoDf$pseudotime <- (1:dim(annoDf)[1])/dim(annoDf)[1]
  annoDf
  # end
}


#' pseudotime analysis based on control clustering result.
#'
#' @param exprM_ctrl normalized expression matrix
#' @param exprM_all normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param startCluster starting Cluster
#' @param endCluster ending Cluster
#' @return a new annoDf with pseudotime value for each cell
#' @export
#' @examples
#' load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/result1/control_all.Rdata")
#' exprM_all <- scale_data_by_seurat(merged_five_mouse_ENCC)
#' exprM_ctrl <- seuset@scale.data; annoDf <- new_annoDf
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' new_annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#' new_annoDf <- pseudotime_based_on_clustering2(exprM_ctrl, exprM_all, annoDf)
pseudotime_based_on_clustering2 <- function(exprM_ctrl, exprM_all, annoDf, threads=3) {
  options(stringsAsFactors = F)
  exprM_ctrl <- exprM_all[,rownames(annoDf)]
  print("1. Get the center of each cluster in control...")
  # 1. Get the center of each cluster
  centerM <- data.frame()
  for (cluster in unique(annoDf$cluster)) {
    cells <- rownames(annoDf[annoDf$cluster == cluster,])
    centerM <- rbind(centerM, t(data.frame(apply(exprM_ctrl[,cells],1,mean))))
  }
  rownames(centerM) <- unique(annoDf$cluster)
  print("2. Sort each center by Insert sorting...")
  # 2. Sort each center by Insert sorting
  orderList <- unique(annoDf$cluster)
  # 3. Output the order of cluster
  print(c("3. The order of the cluster is: ", paste(orderList,collapse=" -> ")))
  print("4. Assign the pseudotime for each cell...")
  # 4. Assign the pseudotime for each cell
  # exprM_merged <- cbind(exprM, t(centerM))
  library(pdist)
  # better to do a dimention reduction
  distM_merged <- pdist(t(exprM_all), centerM)
  distM_merged <- as.matrix(distM_merged)
  #distM_merged <- parDist(x = as.matrix(t(exprM_merged)), method = "euclidean", threads=3)
  rownames(distM_merged) <- colnames(exprM_all)
  colnames(distM_merged) <- rownames(centerM)
  # assign cluster for the new cells
  new_annoDf <- as.data.frame(apply(distM_merged,1,which.min))
  new_annoDf$cluster <- colnames(distM_merged)[new_annoDf[,1]]
  #
  ordered_cells <- c()
  for (i in 1:length(orderList)) {
    direction <- T
    tempCluster <- orderList[i]
    nearestCluster <- orderList[i+1]
    if (i == length(orderList)) {nearestCluster <- orderList[i-1]; direction <- F}
    # print(nearestCluster)
    temp_cells <- rownames(new_annoDf[new_annoDf$cluster==tempCluster,])
    temp_ordered_cells <- temp_cells[order(distM_merged[temp_cells,nearestCluster], decreasing = direction)]
    ordered_cells <- c(ordered_cells, temp_ordered_cells)
  }
  new_annoDf <- new_annoDf[ordered_cells,]
  new_annoDf$pseudotime <- (1:dim(new_annoDf)[1])/dim(new_annoDf)[1]
  new_annoDf[,1] <- NULL
  new_annoDf
  # end
}

#' normalize two pseudotime analysis, make sure the "common" are continuous, must at the beginning
#'
#' @param annoDf control annoDf
#' @param annoDf2 to be aligned
#' @param common starting Cluster
#' @return a new normalized annoDf2
#' @export
#' @examples
#' annoDf <- annoDf_glia; annoDf2 <- annoDf_neuron; common <- c("Cluster1", "Cluster4")
normalize_two_pseudotime <- function(annoDf, annoDf2, common) {
  #ctrl_min <- min(annoDf[annoDf$cluster %in% common,]$pseudotime)
  ctrl_max <- max(annoDf[annoDf$cluster %in% common,]$pseudotime)
  #obj_min <- min(annoDf2[annoDf2$cluster %in% common,]$pseudotime)
  obj_max <- max(annoDf2[annoDf2$cluster %in% common,]$pseudotime)
  # in
  #bias <- obj_min - ctrl_min
  #fold <- (obj_max)/(ctrl_max)
  annoDf2[annoDf2$cluster %in% common,]$pseudotime <- (annoDf2[annoDf2$cluster %in% common,]$pseudotime)/obj_max*ctrl_max
  #
  ctrl_max_rev <- max(1-annoDf[!annoDf$cluster %in% common,]$pseudotime)
  obj_max_rev <- max(1-annoDf2[!annoDf2$cluster %in% common,]$pseudotime)
  annoDf2[!annoDf2$cluster %in% common,]$pseudotime <- 1-(1-annoDf2[!annoDf2$cluster %in% common,]$pseudotime)/obj_max_rev*ctrl_max_rev
  annoDf2
}


#' pseudotime analysis based on clustering.
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @param startCluster starting Cluster
#' @param endCluster ending Cluster
#' @return a new annoDf with pseudotime value for each cell
#' @export
#' @examples
#' load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Aug_30/analysis_4/Ctrl_YFP_ENCC.Rdata")
#' exprM <- seuset@scale.data; annoDf <- all_tsne
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' new_annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
#' annoDf <- new_annoDf
module_detection_along_pseudotime <- function(exprM, annoDf, minNum=100, maxNum=200) {
  moduleList <- list()
  # smooth the exprM to alleviate the effect of dropout
  exprM <- exprM[,rownames(annoDf)]
  exprM_smooth <- t(apply(exprM, 1, smooth))
  colnames(exprM_smooth) <- colnames(exprM)
  # Iteratively get high expression module between two clusters along pseudotime
  clusterList <- unique(annoDf$cluster)
  for (i in 1:length(clusterList)) {
    formerCluster <- clusterList[i]
    latterCluster <- clusterList[i+1]
    if (i == length(clusterList)) {break}
    sub_annoDf <- annoDf[annoDf$cluster %in% c(formerCluster, latterCluster),]
    simulate_ref <- (sub_annoDf$cluster == latterCluster) * 10
    corM <- cor(sub_annoDf$pseudotime, t(exprM_smooth[,rownames(sub_annoDf)]))
    if (i == 1) {
      corM_1 <- -corM
      temp_module <- colnames(corM_1)[corM_1>0.4]
      # perform not well
      # plot(annoDf$pseudotime, exprM_smooth["Xkr4",rownames(annoDf)])
    }
  }
  # end
}

#' pseudotime analysis based on clustering.
#'
#' @param exprM normalized expression matrix contain all cells
#' @param annoDf the annotation of the cells (columns) of the exprM (must have a column "pseudotime")
#' @param sampleNum how many cells for sampling
#' @param knn how many cells for smoothing
#' @return a dataframe to show the gene, log2(sumA/sumB)
#' @export
#' @examples
#' save(exprM, annoDf, file="/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/result1/exprM_3_anno.Rdata")
#' annoDf$cellGroup <- substr(rownames(annoDf), 1, 4)
DEG_pseudotime <- function(exprM, annoDf, sampleNum=200, knn=10, ctrl="negt", comp="post", topNum=200) {
  # end
  ctrl_annoDf <- annoDf[annoDf$cellGroup==ctrl,]
  comp_annoDf <- annoDf[annoDf$cellGroup==comp,]
  # how to sample
  samplePonint <- seq(min(annoDf$pseudotime), max(annoDf$pseudotime), by = (max(annoDf$pseudotime)-min(annoDf$pseudotime))/(sampleNum-1) )
  # get nearest knn cells
  ctrl_cells <- c()
  comp_cells <- c()
  for (i in samplePonint) {
    ctrl_annoDf$dis <- abs(ctrl_annoDf$pseudotime-i)
    temp_ctrl_cells <- rownames(ctrl_annoDf[order(ctrl_annoDf$dis),][1:knn,])
    ctrl_cells <- c(ctrl_cells, temp_ctrl_cells)
    #
    comp_annoDf$dis <- abs(comp_annoDf$pseudotime-i)
    temp_comp_cells <- rownames(comp_annoDf[order(comp_annoDf$dis),][1:knn,])
    comp_cells <- c(comp_cells, temp_comp_cells)
  }
  diff <- exprM[,comp_cells] - exprM[,ctrl_cells]
  # col_anno <-
  diff_genes <- as.data.frame(rowSums(diff)/(sampleNum*knn)) # can't use log2FC, because minus value
  #diff_genes <- as.data.frame(rowSums(diff)/(sampleNum*knn))
  colnames(diff_genes) <- "log2FC"
  diff_genes$gene <- rownames(diff_genes)
  #diff_genes$log2FC <- log2(abs(diff_genes[,1]))
  #diff_genes$log2FC <- sign(diff_genes[,1]) * diff_genes$log2FC
  diff_genes_up_list <- rownames(diff_genes[order(diff_genes$log2FC, decreasing = T),][1:topNum,])
  diff_genes_down_list <- rownames(diff_genes[order(diff_genes$log2FC, decreasing = F),][1:topNum,])
  diff_gene_list <- list()
  diff_gene_list[["up"]] <- diff_genes_up_list
  diff_gene_list[["down"]] <- diff_genes_down_list
  # source("/Users/surgery/Project/HOME/myScript/zxli_lib.R")
  # result_anno <- go_pathway_by_clusterProfiler_mouse(geneList = diff_gene_list)
  # draw_go_kegg_barplot_yshu(result_anno$go_list)
  # draw_go_kegg_barplot_yshu(result_anno$kegg_list)
  diff_gene_list
}

align_two_pseudotime_based_on_clustering <- function(annoDf1, annoDf2) {
  # end
}
