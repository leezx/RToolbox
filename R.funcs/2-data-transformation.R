# source("https://github.com/leezx/bt2m/raw/main/notebooks/unpackaged-code/2-data-transformation.R")

# quickly read large txt file to data.frame, matrix, or dgCMatrix
fast.read.txt <- function(fileName, sep=",", format="dataframe") {
  tmp.raw <- data.table::fread(fileName, sep=sep)
  # get rowname
  tmp.raw <- as.data.frame(tmp.raw)
  tmp.raw <- tmp.raw[!duplicated(tmp.raw[,1]),]
  rownames(tmp.raw) <- tmp.raw[,1]
  tmp.raw[,1] <- NULL
  # output format
  if (format=="dataframe") {
    return(tmp.raw)
    } else if (format=="matrix") {
      return(as.matrix(tmp.raw))
    } else if (format=="dgCMatrix") {
      library(Matrix)
      return(as(as.matrix(tmp.raw), "dgCMatrix"))
    } else {
      stop("please input correct format!!!")
    }
}

# compress gene expression for seurat
compress.expression.seurat <- function(seuratObj, compress.group) {
  # too slow
  # compressed.exprM <- aggregate(t(as.matrix(seuratObj@assays$RNA@data)),
  #                              list(seuratObj@meta.data[,compress.group]), mean)
  #
  compressed.exprM <- data.frame()
  tmp.exprM <- as.matrix(seuratObj@assays$RNA@data)
  for (i in unique(seuratObj@meta.data[,compress.group])) {
    # print(i)
    tmp.expr <- rowMeans(tmp.exprM[,seuratObj@meta.data[,compress.group]==i])
    compressed.exprM <- rbind(compressed.exprM, tmp.expr)
  }
  compressed.exprM <- as.data.frame(t(compressed.exprM))
  colnames(compressed.exprM) <- unique(seuratObj@meta.data[,compress.group])
  rownames(compressed.exprM) <- rownames(tmp.exprM)
  return(compressed.exprM)
}


#' scale the matrix
#' @param raw.data raw.data
#' @param max.value max.value after scale
#'
#' @return scale.data
#' @export
#' @examples
#' scale.data <- pre.scale.data(logcounts(sce_HSCR))
#'
pre.scale.data <- function(raw.data=logcounts(sce_HSCR), max.value=2) {
  # max.value <- 2
  scale.data <- t(x = scale(x = t(x = as.matrix(x = raw.data)), center = T, scale = T))
  scale.data[is.na(scale.data)] <- 0
  scale.data[scale.data > max.value] <- max.value
  scale.data[scale.data < -max.value] <- -max.value
  return(scale.data)
}

## 5. float matrix to interger matrix
df_float_to_int <- function(counts){
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  return(counts)
}

## 13. read csv in transform
read.tcsv <- function(file, header=TRUE, sep=",", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}
