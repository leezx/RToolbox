# VlnPlot has no output, cannot resolve the problem. 
# So, create a simple version by GPT.
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

VlnPlot_simple <- function(
  object,
  features,
  group.by = NULL,
  assay = NULL,
  layer = NULL,
  log = FALSE,
  cols = NULL,
  pt.size = 0,
  alpha = 1
) {
  assay <- assay %||% DefaultAssay(object)
  # 取分组信息
  groups <- if (!is.null(group.by) && group.by %in% colnames(object[[]])) {
    object[[group.by]][, 1, drop = TRUE]
  } else {
    Idents(object)
  }
  groups <- droplevels(as.factor(groups))

  dat_list <- lapply(features, function(f) {
    vals <- tryCatch(
      FetchData(object, vars = f, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (is.null(vals)) stop(sprintf("Feature %s not found.", f))
    y <- vals[[1]]
    if (log) y <- log1p(y)
    data.frame(
      feature = f,
      value = y,
      group = groups
    )
  })
  dat <- bind_rows(dat_list)

  p <- ggplot(dat, aes(x = group, y = value, fill = feature)) +
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = alpha) +
    stat_summary(fun = median, geom = "point", size = 0.5, color = "black") +
    facet_wrap(~feature, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  if (pt.size > 0) {
    p <- p + geom_jitter(width = 0.15, size = pt.size, alpha = 0.5)
  }
  if (!is.null(cols)) {
    p <- p + scale_fill_manual(values = cols)
  } else {
    p <- p + scale_fill_brewer(palette = "Set2")
  }

  p
}

