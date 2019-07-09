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
plot.dotplot.human <- function (sce, genes.plot, use="cellGroup2",xAngle=60, cols.use = c("lightgrey", "blue"), 
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
