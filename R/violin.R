
################################################################
#' a stacked version of VlnPlot in Seurat, from public website
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
StackedVlnPlot.rowCluster <- function(obj, features, angle.gene = "30",
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
  # addition funcitons
  modify_vlnplot <- function(obj,
                           feature,
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
    p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
      xlab("") + ylab(feature) + ggtitle("") +
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            # axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = 13, angle = angle.gene, face = "italic", vjust = 0.5),
            axis.title.y = element_text(size = 12, angle = 0, face = "plain", vjust = 0.5),
            axis.text.x = element_text(size = 8, face = "plain"),
            plot.margin = plot.margin ) & coord_flip()
    return(p)
  }

  ## extract the max value of the y axis
  extract_max <- function(p){
    ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  # 
  # main code
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[1]] <- plot_list[[1]] +
    theme(axis.text.y = element_text(angle = 0, size = 15, face = "bold"), axis.ticks.y = element_line())

  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x +
                             scale_y_continuous(breaks = c(y)) +
                             expand_limits(y = y))

  p <- patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  # p <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, labels = features) # not aligned
  return(p)
}

StackedVlnPlot.rowGene <- function(obj, features,
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
  # addition funcitons
  modify_vlnplot <- function(obj,
                           feature,
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
    p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
      xlab("") + ylab(feature) + ggtitle("") +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 20, angle = 0, face = "bold.italic", vjust = 0.5),
            axis.text.y = element_text(size = 10),
            plot.margin = plot.margin )
    return(p)
  }

  ## extract the max value of the y axis
  extract_max <- function(p){
    ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  # 
  # main code
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 75, size = 15), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x +
                             scale_y_continuous(breaks = c(y)) +
                             expand_limits(y = y))

  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  # p <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, labels = features) # not aligned
  return(p)
}

# example: EllyLab/mouse/singleCell/case/Kif7_ENCC/Kif7-integration/1-integration_public_and_Kif7.ipynb
# example: EllyLab/mouse/singleCell/case/Vcl_ENCC/Vcl_ENCCs_aggregate_analysis.ipynb#
#' Draw violin plot
#' @param exprData the expression matrix to be used
#' @param cellAnno the annotation of the cells
#' @param genes genes wanted to show
#' @param use the group selected in the cellAnno
#' @param direction direction of the violin plot
#' @param showGroup whether to show the group names
#' @param groupOrder the order of the group
#'
#' @return a ggplot violin object
#' @export
#' @examples
#' options(repr.plot.width=4, repr.plot.height=9)
#' plot.violin(exprData=seuset@data, cellAnno=all_tsne_2, genes=key_genes, use="final", direction = "h")
#'
#' all_tsne_3 <- subset(all_tsne_2, final %in% c('ctrl:c4','kif7:c4','ctrl:c6','kif7:c6'))
#' options(repr.plot.width=4, repr.plot.height=4)
#' plot.violin(exprData=seuset@data, cellAnno=all_tsne_3, genes=key_genes, use="final", direction = "h",
#'             groupOrder = c('ctrl:c4','kif7:c4','ctrl:c6','kif7:c6'))
#'
plot.violin <- function(exprData, cellAnno, genes, use, direction="h", showGroup=T, groupOrder=NULL) {
  # exprData=seuset@data; cellAnno=all_tsne_2; genes=key_genes; use="final"; showGroup=T; direction="h"
  ##############################
  if (showGroup==T) {
    groupSize <- 12
  } else {groupSize <- 0}
  ##############################
  # input balanced-absolute value expression matrix
  # prepare data for violin plot
  options(stringsAsFactors = F)
  exprM <- cbind(as.matrix(t(exprData[genes, rownames(cellAnno)])), group=cellAnno[,use])
  exprM_melt <- reshape::melt(as.data.frame(exprM), id.vars = c("group"))
  exprM_melt$value <- as.double(exprM_melt$value)
  if (!is.null(groupOrder)) {exprM_melt$group <- factor(exprM_melt$group, levels=groupOrder)}
  #print(exprM_melt[1:2,])
  #print(max(exprM_melt$value))
  ###############################
  # plot violin
  if (direction=="h") {
    #options(repr.plot.width=4, repr.plot.height=9)
    plot <- ggplot(exprM_melt, aes(x=variable, y=value)) +
      facet_wrap( ~ group, ncol=1, scales = "fix",strip.position = "left") + # free_y
      labs(x = "", y = "Scaled UMIs\n") +
      geom_violin(trim = T, na.rm = T, aes(fill = group), colour = "black", scale="width", size=0.3, bw=0.1) +
      scale_y_continuous(position="right") + # , limits=c(0, 3), breaks = seq(0, 3, length.out = 4)
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(legend.title=element_blank()) +
      theme(legend.position = "none") +
      theme(strip.background = element_rect(fill = NA, color = NA)) + # strip background color
      theme(strip.placement = "outside", strip.text.x = element_text(face="plain", size = 0), #italic
            strip.text.y = element_text(face="plain", size = groupSize, angle=-90)) +
      theme(panel.spacing=unit(.3, "lines"),panel.border = element_rect(color = "black", fill = NA,
                                                                        size = 0.2,colour = "black")) + #line size
      theme(axis.text.x  = element_text(face="italic", angle=90, size = 13, color = "black", vjust=0.5),
            axis.text.y  = element_text(face="plain", size = 8, color = "black"),
            axis.title =element_text(size = 12)) +
      geom_hline(yintercept=max(exprM_melt$value)/2, linetype="dashed", color = "gray10", size=0.08)
    #scale_color_manual(values=myColors)

    plot
  } else if (direction=="v") {
    plot <- ggplot(exprM_melt, aes(x=group, y=value)) +
      facet_wrap( ~ variable, ncol=1, scales = "fix",strip.position = "left") + # free_y
      labs(x = "", y = "Scaled UMIs\n") +
      geom_violin(trim = T, na.rm = T, aes(fill = variable), colour = "black", scale="width", size=0.3, bw=0.1) +
      scale_y_continuous(position="right") + # , limits=c(0, 3), breaks = seq(0, 3, length.out = 4)
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(legend.title=element_blank()) +
      theme(legend.position = "none") +
      theme(strip.background = element_rect(fill = NA, color = NA)) + # strip background color
      theme(strip.placement = "outside", strip.text.x = element_text(face="plain", size = 0), #italic
            strip.text.y = element_text(face="plain", size = groupSize, angle=-90)) +
      theme(panel.spacing=unit(.3, "lines"),panel.border = element_rect(color = "black", fill = NA,
                                                                        size = 0.2,colour = "black")) + #line size
      theme(axis.text.x  = element_text(face="italic", angle=90, size = 13, color = "black", vjust=0.5),
            axis.text.y  = element_text(face="plain", size = 8, color = "black"),
            axis.title =element_text(size = 12)) +
      geom_hline(yintercept=max(exprM_melt$value)/2, linetype="dashed", color = "gray10", size=0.08)
    #scale_color_manual(values=myColors)

    plot
  } else (strop("direction can only be h or v!"))
}

