

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
#' result <- pseudotime_based_on_clustering(exprM, annoDf, startCluster, endCluster)
pseudotime_based_on_clustering <- function(exprM, annoDf, startCluster, endCluster, threads=3) {
  options(stringsAsFactors = F)
  print("1. Get the center of each clustering...")
  # 1. Get the center of each clustering
  centerM <- data.frame()
  for (cluster in unique(annoDf$cluster)) {
    cells <- rownames(annoDf[annoDf$cluster == cluster,])
    centerM <- rbind(centerM, t(data.frame(apply(exprM[,cells],1,mean))))
  }
  rownames(centerM) <- unique(annoDf$cluster)
  print("2. Iteratively inserting each center...")
  # 2. Iteratively inserting each center
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
    tempCluster <- orderList[i]
    nearestCluster <- orderList[i+1]
    if (i == length(orderList)) {nearestCluster <- orderList[i-1]}
    # print(nearestCluster)
    temp_cells <- rownames(annoDf[annoDf$cluster==tempCluster,])
    temp_ordered_cells <- temp_cells[order(distM_merged[temp_cells,nearestCluster], decreasing = T)]
    ordered_cells <- c(ordered_cells, temp_ordered_cells)
  }
  annoDf <- annoDf[ordered_cells,]
  annoDf$pseudotime <- (1:dim(annoDf)[1])/dim(annoDf)[1]
  annoDf
  # end
}

align_two_pseudotime_based_on_clustering <- function(annoDf1, annoDf2) {
  # end
}

module_detection_along_pseudotime <- function(annoDf1, annoDf2) {
  # end
}
