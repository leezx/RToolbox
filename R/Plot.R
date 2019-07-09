#' draw violin plot
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

