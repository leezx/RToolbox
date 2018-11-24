
#' prepare_raw_count_matrix_from_10x_data
#'
#' @param cellranger_pipestance_path cellranger_pipestance_path
#' @return a gbm object
#' @export
#' @examples
prepare_raw_count_matrix_from_10x_data <- function(cellranger_pipestance_path) {
  # library(cellrangerRkit) # not support by official 10x
  # pbmc33k.data <- Read10X("~/Projects/datasets/pbmc33k/filtered_gene_bc_matrices/hg19/")
  pbmc33k  <- new("seurat", raw.data = pbmc33k.data)
  gbm <- load_cellranger_matrix(cellranger_pipestance_path)
  gbm <- gbm[!duplicated(gbm@featureData$symbol),]
  gbm
  # end
}


#' add_gender_info_to_annoDf
#'
#' @param exprM normalized expression matrix
#' @param annoDf the annotation of the cells (columns) of the exprM
#' @return a ggplot object
#' @export
#' @examples
add_gender_info_to_annoDf <- function(exprM, annoDf) {
  y_chr_genes <- c("Kdm5d","Eif2s3y","Gm29650","Uty","Ddx3y","Erdr1")
  gender_expr <- as.matrix(exprM[y_chr_genes,])
  annoDf$gender <- "female"
  annoDf[colSums(gender_expr)>1,]$gender <- "male"
  #return
  annoDf
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
