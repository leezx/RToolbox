#' Draw GSEA GO KEGG barplot for a gene list
#' @param anno_list a gene list
#' @param type type of the input data
#'
#' @return a dataframe
#' @export
#' @examples
#' options(repr.plot.width=10, repr.plot.height=6)
#' go_list <- plot.gsea.GO.KEGG.barplot.batch(gsea_list$go_list, type = "GO")
#' kegg_list <- plot.gsea.GO.KEGG.barplot.batch(gsea_list$kegg_list, type = "KEGG")
#'
plot.gsea.GO.KEGG.barplot.batch <- function(anno_list=gsea_list$go_list, type, f.length=40) {
  library(ggplot2)
  library(RColorBrewer)
  tmpcolors <- brewer.pal(12,"Set3")[1:length(anno_list)]
  nameList <- names(anno_list)
  if (is.null(nameList)) {
    nameList <- 1:length(anno_list)
    names(tmpcolors) <- nameList
    }
  j <- 0
  merged_list <- list()
  for (i in nameList) {
    j <- j + 1
    barplot_df <- anno_list[[i]]@result
    if (length(barplot_df)<2 | is.null(barplot_df)) {next}
    barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F),]
    barplot_df$Description <- factor(barplot_df$Description, levels=rev(barplot_df$Description))
    maxpvalue <- max(-log10(barplot_df$pvalue))
    # filter duplication
    left_go <- c()
    barplot_df <- barplot_df[!duplicated(barplot_df$core_enrichment),]
    # remove duplicate GO or KEGG terms
    for (one in 1:length(barplot_df$ID)) {
      if (one ==1) {left_go <- c(left_go, one); next}
      leftGenes <- strsplit(paste(barplot_df[left_go,]$core_enrichment, collapse = "/"), split = "/")
      tmpGenes <- strsplit(paste(barplot_df[one,]$core_enrichment, collapse = "/"), split = "/")
      if (length(intersect(leftGenes, tmpGenes))/length(tmpGenes) > 0.7) {next}
      left_go <- c(left_go, one)
    }
    barplot_df <- barplot_df[left_go,]
    # remove terms with too long description
    barplot_df <- barplot_df[sapply(as.character(barplot_df$Description), nchar) < f.length,]
    if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
    ###
    g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
    geom_bar(stat="identity", fill = tmpcolors[j]) +
    geom_text(aes(label=round(NES, digits=2)), color=ifelse(barplot_df$NES > 0,'red','blue'),
              vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
    ylim(0, maxpvalue*1.1) +
    coord_flip() +
    labs(x = "", y = "-Log10(P-value)", title=paste(type, i, sep=":")) + 
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 16, color = "black", face = "bold"), 
          axis.text.x = element_text(size = 11, color = "black", face = "plain")) +
    # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),panel.border = element_blank(),
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
          plot.margin=unit(c(0,0,0,0), "cm"), panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black'))
    # + scale_fill_manual(values=colorRampPalette(c("#3176e0", "#8831e0", "#e031b4"))(15))
    plot(g)
    merged_list[[i]] <- barplot_df
  }
  merged_list
}

#' Draw GO KEGG barplot for a gene list
#' @param anno_list a gene list
#' @param type type of the input data
#'
#' @return a dataframe
#' @export
#' @examples
#' options(repr.plot.width=4, repr.plot.height=9)
#' go_list <- plot.GO.KEGG.barplot.batch(result$go_list, type="GO")
#'
plot.ora.GO.KEGG.barplot.batch <- function(anno_list=go_list, type, f.length=40) {
  # anno_list=kegg_list
  # colors <- c("#FB8072", "#B3DE69", "#BC80BD")
  library(ggplot2)
  library(RColorBrewer)
  tmpcolors <- brewer.pal(12,"Set3")[1:length(anno_list)]
  nameList <- names(anno_list)
  if (is.null(nameList)) {
    nameList <- 1:length(anno_list)
    names(tmpcolors) <- nameList
    }
  j <- 0
  merged_list <- list()
  for (i in nameList) {
    j <- j + 1
    barplot_df <- anno_list[[i]]@result
    if (length(barplot_df)<2 | is.null(barplot_df)) {next}
    barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F),]
    barplot_df$Description <- factor(barplot_df$Description, levels=rev(barplot_df$Description))
    maxpvalue <- max(-log10(barplot_df$pvalue))
    # filter duplication
    left_go <- c()
    barplot_df <- barplot_df[!duplicated(barplot_df$geneID),]
    for (one in 1:length(barplot_df$ID)) {
      if (one ==1) {left_go <- c(left_go, one); next}
      leftGenes <- strsplit(paste(barplot_df[left_go,]$geneID, collapse = "/"), split = "/")
      tmpGenes <- strsplit(paste(barplot_df[one,]$geneID, collapse = "/"), split = "/")
      if (length(intersect(leftGenes, tmpGenes))/length(tmpGenes) > 0.7) {next}
      left_go <- c(left_go, one)
    }
    barplot_df <- barplot_df[left_go,]
    # remove terms with too long description
    barplot_df <- barplot_df[sapply(as.character(barplot_df$Description), nchar) < f.length,]
    #
    if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
    g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
    geom_bar(stat="identity", fill = tmpcolors[j]) +
    geom_text(aes(label=Count),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
    ylim(0, maxpvalue*1.1) +
    coord_flip() +
    labs(x = "", y = "-Log10(P-value)", title=paste(type, i, sep=":")) + 
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 16, color = "black", face = "bold"), axis.text.x = element_text(size = 11, color = "black", face = "plain")) +
    # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),panel.border = element_blank(),
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin=unit(c(0,0,0,0), "cm"), panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black'))
    # + scale_fill_manual(values=colorRampPalette(c("#3176e0", "#8831e0", "#e031b4"))(15))
    plot(g)
    # save, 700, 350
    merged_list[[i]] <- barplot_df
  }
  merged_list
}

#' Draw barplot for GO dataframe
#' @param barplot_df GO dataframe from clusterProfiler
#'
#' @return a ggplot barplot object
#' @export
#' @examples
#' options(repr.plot.width=4, repr.plot.height=9)
#' plot.GO.barplot(df)
#'
plot.GO.barplot2 <- function(barplot_df, color="random") {
  library(Hmisc)
  library(stringr)
  library(RColorBrewer)
  if (color=="random") {
    color <- sample(brewer.pal(12, "Set3"), 1)
  }
  # colors <- brewer.pal(10,"Paired")
  for (i in 1:dim(barplot_df)[1]) {
    barplot_df[i,]$Description <- capitalize(as.character(barplot_df[i,]$Description))
  }
  # if (length(barplot_df)<2 | is.null(barplot_df)) {next}
  barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F),]
  barplot_df$Description <- factor(barplot_df$Description, levels=rev(barplot_df$Description))
  maxpvalue <- max(-log10(barplot_df$pvalue))
  # if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
  g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
  geom_bar(stat="identity", color=color, fill=color, alpha=0.8) +
  geom_text(aes(label=Count),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
  ylim(0, maxpvalue*1.1) +
  coord_flip() +
  labs(x = "", y = "-Log10(P-value)", title="") + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 14, color = "black", face = "plain"), 
    axis.text.x = element_text(size = 12, color = "black", face = "plain"),
    axis.title =element_text(size = 15)) +
  # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), 
  # axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),
  # panel.border = element_blank(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin=unit(c(0,0,0,0), "cm"), 
        panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  scale_x_discrete(labels=function(x) str_wrap(x, width=25))
  # scale_fill_manual(values=rev(colors))
  g
}

plot.GO.barplot <- function (barplot_df, color = "random") 
{
    #barplot_df <- tmp.GO.sub
    #color = brewer.pal(8,"Set2")[2]

    #
    library(Hmisc)
    library(stringr)
    library(RColorBrewer)
    if (color == "random") {
        color <- sample(brewer.pal(12, "Set3"), 1)
    }
    for (i in 1:dim(barplot_df)[1]) {
        barplot_df[i, ]$Description <- capitalize(as.character(barplot_df[i, 
            ]$Description))
    }

    barplot_df$Description2 <- paste(barplot_df$Description," (", barplot_df$Count, ")", sep = "")

        barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F), 
            ]
        barplot_df$Description2 <- factor(barplot_df$Description2, 
            levels = rev(barplot_df$Description2))
        maxpvalue <- max(-log10(barplot_df$pvalue))

    maxpvalue*1.1

    #options(repr.plot.width=7, repr.plot.height=4)

    #color = brewer.pal(8,"Set2")[2]

    g <- ggplot(data = barplot_df, aes(x = Description2, y = -log10(pvalue))) + 
        geom_bar(stat = "identity", color = color, fill = color, alpha = 0.8, width=.5) + 
        #geom_text(aes(label = Count), color = "black", 
        #    vjust = 0.4, hjust = -0.5, size = 3, fontface = "bold") + 
        #ylim(0, maxpvalue * 1.1) + 
        coord_flip() + labs(x = "", 
        y = "-Log10(P-value)", title = "") + theme_bw() + 
        theme(legend.position = "none") + 
        theme(axis.text.y = element_text(size = 14, color = "black", 
            face = "plain"), axis.text.x = element_text(size = 12, 
            color = "black", face = "plain"), axis.title = element_text(size = 15)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.border = element_blank()) + 
        theme(axis.line = element_line(color = "black")) + 
        scale_x_discrete(labels = function(x) str_wrap(x, width = 55), expand = c(0.07, 0)) +
        scale_y_continuous(limits = c(0, maxpvalue*1.1), breaks = seq(0, maxpvalue*1.2, by = 3), expand = c(0, 0))
    g
}
                   
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
