

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
  print("# 4. Assign the pseudotime for each cell...")
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
#' exprM_ctrl <- seuset@scale.data; exprM_all <- scale_data_by_seurat(merged_five_mouse_ENCC); annoDf <- new_annoDf
#' startCluster<-"Cluster1"; endCluster<-"Cluster6"
#' new_annoDf <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
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
  print("# 4. Assign the pseudotime for each cell...")
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

scale_data_by_seurat <- function(rawCountM) {
  seuset <- CreateSeuratObject(raw.data = rawCountM, min.cells = 2, min.genes = 100)
  # VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 2)
  # GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")
  seuset <- FilterCells(object = seuset, subset.names = c("nUMI"))
  # FilterCells(object, subset.names, low.thresholds, high.thresholds,cells.use = NULL)
  seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)
  seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  # length(x = seuset@var.genes)
  seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI"))
  seuset@scale.data
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


align_two_pseudotime_based_on_clustering <- function(annoDf1, annoDf2) {
  # end
}
