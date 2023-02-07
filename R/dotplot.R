# source("https://github.com/leezx/bt2m/raw/main/notebooks/unpackaged-code/dotplot.R")

#' origninal DotPlot in Seurat cannot order features

source("https://github.com/satijalab/seurat/raw/8da35ba8dc2d60f504f7cf276efd684dc19a418c/R/utilities.R")
DotPlot_order <- function (object, assay = NULL, features, cols = c("lightgrey",
    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
    group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius",
    scale.min = NA, scale.max = NA)
{
    if(is.null(assay)) assay <- DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot,
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    #
    # print(head(data.plot))
    data.plot$features.plot <- factor(data.plot$features.plot, levels = features)
    #
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by)) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}


#' Draw dotplot for human SC3 object
#' @param sce SC3 object
#' @param genes.plot genes
#' @param use group used in SC3 object
#' @param xAngle angle of x title
#' @param cols.use colors used
#' @param col.min col.min
#' @param col.max col.max
#' @param dot.min dot.min
#' @param dot.scale whether to scale or not
#' @param group.order group.order
#' @param scale.by scale.by
#' @param scale.min scale.min
#' @param scale.max scale.max
#' @param group.by group.by
#' @param title title
#' @param plot.legend plot legend or not
#' @param do.return return object or not
#' @param x.lab.rot rotate x title or not
#'
#' @return a ggplot dotplot
#' @export
#' @examples
#' options(repr.plot.width=4, repr.plot.height=6)
#' plot.dotplot.human(sce = sce_HSCR_pure, genes.plot = uniquegenes2, use="cellGroup", group.order = group.order,
#'           scale.min=0, scale.max=100, title="", plot.legend = T, xAngle=90)
#'
plot.dotplot.SC3 <- function (sce, genes.plot, use="cellGroup2",xAngle=60, cols.use = c("lightgrey", "blue"),
                              col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.order = c(),
                              scale.by = "radius", scale.min = NA, scale.max = NA, group.by, title="",
                              plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE)
{
  #cols.use = c("lightgrey", "blue"); genes.plot = tmp$human; sce <- sce_comp; use="cellGroup2";
  #col.min = -2.5; col.max = 2.5; dot.min = 0; dot.scale = 6; group.order = c();xAngle=60;
  #scale.by = "radius"; scale.min = 0; scale.max = 100; title="";
  #plot.legend = T; do.return = FALSE; x.lab.rot = FALSE
  #
  library(dplyr)
  library(tidyr)
  PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  #if (!missing(x = group.by)) {
  #    object <- SetAllIdent(object = object, id = group.by)
  #}
  #data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot <- as.data.frame(t(logcounts(sce)[genes.plot,]))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  #data.to.plot$id <- object@ident
  data.to.plot$id <- colData(sce)[,use]
  data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression,
                                                                            threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    dplyr::mutate(avg.exp.scale = scale(x = avg.exp)) %>% dplyr::mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                                                                               max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  # group.order
  if (length(group.order)>0) {
    data.to.plot$id <- factor(data.to.plot$id, levels = group.order)
  }
  # plot
  p <- ggplot(data = data.to.plot, mapping = aes(y = genes.plot, x = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min,scale.max)) +
    theme_bw(base_line_size = 0.1, base_rect_size = 0.1) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.border = element_rect(colour = 'black'), panel.grid.minor = element_blank(),
          plot.background = element_blank(), panel.grid.major = element_blank(),
          axis.text.x  = element_text(face="plain", angle=xAngle, size = 10, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="italic", size = 10, color = "black")) +
    labs(title = title) +
    guides(colour = guide_colourbar(title = "Relative expression", title.position = "left",
                                    direction="vertical", title.theme = element_text(angle=90), title.hjust=-4),
           size = guide_legend(title = "Expressed percentage (%)", title.position = "left",
                               direction="vertical", title.theme = element_text(angle=90), title.hjust=-15))
  ####################
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  x.lab.rot <- F
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
