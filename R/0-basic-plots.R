# source("https://github.com/leezx/bt2m/raw/main/notebooks/unpackaged-code/violin.R")

# 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


# get color theme
# https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
library(RColorBrewer)
example_col <- rev(brewer.pal(10,"Set3"))
example_col <- RColorBrewer::brewer.pal(5,"Set2")

library(ggplot2)
DoHeatmap(object = pbmc_small) + scale_fill_gradientn(colors = c("blue", "white", "red"))
library(viridis)
DoHeatmap(object = pbmc_small) + scale_fill_viridis()

# geom_point with geom_label_repel
# path: EllyLab/human/bulkRNA/bulk%20RNA-seq.ipynb

# facet_wrap
# path:
# facet_grid
# path: EllyLab/mouse/singleCell/case/GBS_ENCC/Tg%20(GBS-GFP)%20and%20metabolic.ipynb

# venn plot
# path: EllyLab/human/singleCell/HSCR/5-AS_and_metabolism_CGS.ipynb
require(yyplot)
require(ggplot2)
# require(ggforce)
options(repr.plot.width=5, repr.plot.height=4)
g2 <- ggplot(d) +
  geom_circle(aes_(x0 = ~x, y0 = ~y, r = ~diameters/2, fill = ~labels, color = ~labels), size=1.5) +
  # geom_text(aes_(x = ~x, y = ~y, label = ~labels)) +
  coord_fixed() +
  theme_void() +
  scale_color_manual(values = c("darkred","red")) +
  scale_fill_manual(values = alpha(c("darkred","red"), .2)) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 15), legend.position = c(0.65, 0.35))+
  geom_text(x=0.20, y=0.5, label="843", size=6) +
  geom_text(x=0.50, y=0.5, label="548", size=8, fontface="bold") +
  geom_text(x=0.80, y=0.5, label="1080", size=6)
# scale_color_manual(values = c("#E41A1C", "#984EA3", 'grey50')) +
# scale_fill_manual(values = alpha(c("#E41A1C", "#984EA3", 'grey50'), .2))
# scale_color_manual(values = c("#E41A1C", "#984EA3", 'grey50',"blue")) +
# scale_fill_manual(values = alpha(c("#E41A1C", "#984EA3", 'grey50',"blue"), .2))
g2

# scatter with fit line
# path: project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
options(repr.plot.width=4, repr.plot.height=4)
p <- ggplot(data = all.umap.peaks, aes(x = pseudotime, y = peakNum)) +
  geom_point(aes(color=cluster), size=1) +
  # geom_smooth(method = "lm", se=F, color="black", formula = y ~ x) +
  geom_smooth(method = "loess", se=F, color="black") +
  theme_bw() +
  labs(x = "Expression of \n37 RNA splicing regulators", y = "AS frequency", title = "") +
  theme(legend.position = "none") +
  theme(axis.text  = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_color_manual(values=tmp.colors) #+
# geom_text(x = 0.2, y = 1.47, label = "P-value < 2.2E-16", size = 4) +
# annotate("text", label = "P-value < 2.2E-16", x = 0.2, y = 1.47, size = 4, colour = "black") +
# scale_x_continuous(limits = c(-0.5, 0.4)) +
# scale_y_continuous(limits = c(0.5, 1.5))
p

tmp <- lm(pseudotime ~ peakNum, all.umap.peaks)
summary(tmp)

# density ridges
# path: project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
library(ggridges)
options(repr.plot.width=4, repr.plot.height=3)
p1 <- ggplot(all.umap) +
  aes(x= pseudotime, y = sample, color = sample, fill=sample) +
  geom_density_ridges(alpha = 0.75, span=0.2) +
  # theme_ridges() +
  theme_bw() +
  #scale_color_manual(values=rev(sample.colors)) +
  #scale_fill_manual(values=rev(sample.colors)) +
  labs(x = "Pseudotime", y = "", title = "") +
  theme(legend.position = "none") +
  # empty border, ticks, text
  theme(panel.border = element_blank(), # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 0, size = 10, vjust = 5, face = "plain"),
        axis.title = element_text(size = 13, face = "plain"),
        axis.text.y = element_text(size = 11, face = "bold")
  ) #+
# guides(y.sec = guide_axis_label_trans(~paste("Ha!", .y))) #+
# scale_y_discrete(sec.axis = sec_axis(name = "Density")) # not support
# theme(axis.title=element_blank(), axis.text=element_blank(), )
p1

# boxplot
draw_boxplot <- function(data=pca_data, use=c("colour_by", "gene_num"), xlab="", ylab="", title=""){
  data$x <- data[,use[1]]
  data$y <- data[,use[2]]
  n <- length(unique(data$x))
  ggplot(data, aes(x=x, y=y, fill=x)) + geom_boxplot() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = xlab, y = ylab, title = title) +
    theme(panel.spacing=unit(.5, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x  = element_text(face="plain", angle=90, size = 12, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="plain", size = 12, color = "black"),
          axis.title =element_text(size = 15)) #+
  #scale_fill_manual(values=brewer.pal(6,"Set2"))
}

draw_violin_plot <- function(){
  library(easyGgplot2)
  plot <- ggplot(MetaData.F, aes(factor(cellGroup), prolif_rate)) + geom_violin(trim = FALSE, aes(fill = cellGroup), colour = "white",scale = "width") + scale_y_continuous(position="left", breaks = seq(0, 1, by = 0.2), limits=c(-0.1, 1.1)) + geom_jitter(height = 0, width = 0.2, size=0.3) + geom_boxplot(width=0.1)

  ggplot2.customize(plot, removePanelBorder=TRUE,removePanelGrid=TRUE,backgroundColor="white",showLegend=FALSE) + theme(panel.spacing=unit(.5, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    labs(x = "", y = "proliferation index\n", title = "")  +
    theme(axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.text.x  = element_text(face="bold", angle=30, size = 14, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="plain", size = 14, color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title =element_text(size = 20)
    ) + scale_fill_manual(values=brewer.pal(6,"Set2"))
}

draw_scatter_plot <- function(){
  # with out border
  ggplot(df_simlr_5, aes(x=x, y=y, color=new_SIMLR)) +
    geom_point(size=3.5, alpha=0.55) +
    labs(x = " ",y = " ", title = " ") +
    theme_bw() +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
    theme(axis.title = element_blank() ,axis.text = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank()) +
    theme(legend.title=element_blank()) +
    scale_color_manual(values=brewer.pal(5,"Set2"))
  # with border
  ggplot(SAG_PCA[sce$cellName[(sce$cellGroup!="SAG0_4")],], aes(x=PC1, y=PC2, color=colour_by)) +
    geom_point(size=1.5, alpha=1) +
    labs(x = "PC1",y = "PC2", title = "SAG0_10") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
    theme(legend.title=element_blank()) +
    scale_color_manual(values=brewer.pal(6,"Set2"))
}

draw_barplot <- function() {
  # position_dodge
  a <- table(sHSCR_tsne$data[sHSCR_tsne$data$shape_by=="HSCR_20c7",]$colour_by)
  b <- table(sHSCR_tsne$data[sHSCR_tsne$data$shape_by=="HSCR_5c3",]$colour_by)
  a1 <- data.frame(celltype="HSCR_20c7", cluster=names(a), count=as.vector(a))
  b1 <- data.frame(celltype="HSCR_5c3", cluster=names(b), count=as.vector(b))
  countTable <- rbind(a1, b1, c("HSCR_5c3", "Cluster5", 0))
  ggplot(data=countTable, aes(x=cluster, y=as.integer(count), fill=celltype)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_brewer(palette="Blues") +
    labs(x = "",y = "cell count", title = " ") +
    theme(axis.text.x  = element_text(face="plain", angle=30, size = 15, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="plain", size = 15, color = "black"),
          axis.title =element_text(size = 15))
  # fill/stack
  SAG4 <- table(tmp_group$my_sc3_4)
  SAG10 <- table(tmp_group$my_sc3_4)
  SAG41 <- data.frame(celltype="SAG0_4", cluster=names(SAG4), count=as.vector(SAG4))
  SAG101 <- data.frame(celltype="SAG0_10", cluster=names(SAG10), count=as.vector(SAG10))
  countTable <- rbind(SAG41, SAG101)
  ggplot(data=countTable, aes(x=celltype, y=as.integer(count), fill=cluster)) +
    geom_bar(stat="identity", position="fill") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "",y = "cell percentage", title = " ") +
    theme(axis.text.x  = element_text(face="plain", angle=30, size = 15, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="plain", size = 15, color = "black"),
          axis.title =element_text(size = 15)) +
    scale_fill_manual(values=brewer.pal(6,"Set2")[2:5])
}

