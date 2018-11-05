
prepare_raw_count_matrix_from_10x_data <- function(cellranger_pipestance_path) {
  # library(cellrangerRkit)
  gbm <- load_cellranger_matrix(cellranger_pipestance_path)
  gbm <- gbm[!duplicated(gbm@featureData$symbol),]
  gbm
  # end
}
