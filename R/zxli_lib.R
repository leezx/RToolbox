## purpose: wrap all the frequently-used functions or commands.
# setwd("/Users/surgery/Project/HOME/myScript")

go_pathway_by_WebGestaltR <- function(geneList=markerList) {
  library(WebGestaltR)
  go_list <- list()
  kegg_list <- list()
  for (i in names(geneList)) {
    genes <- geneList[[i]]
    projectName <- i
    enrichResult_BP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                            enrichDatabase="geneontology_Biological_Process_noRedundant", 
                            enrichDatabaseFile=NULL, 
                            enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,
                            interestGeneFile = NULL, interestGene=genes, 
                            interestGeneType="genesymbol", collapseMethod="mean", 
                            referenceGeneFile=NULL, referenceGene=genes,
                            referenceGeneType="genesymbol", referenceSet=NULL, minNum=5, 
                            maxNum=2000, fdrMethod="BH", sigMethod="top", 
                            fdrThr=0.05, topThr=20, dNum=20,
                            perNum=1000, lNum=20,is.output=TRUE,outputDirectory=getwd(),
                            projectName=paste(projectName,"BP", sep = "."), 
                            keepGSEAFolder=FALSE,
                            methodType="R",dagColor="binary",
                            hostName="http://www.webgestalt.org/")
    go_list[[i]] <- enrichResult_BP
    enrichResult_pathway <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                            enrichDatabase="pathway_KEGG", 
                            enrichDatabaseFile=NULL, 
                            enrichDatabaseType=NULL, enrichDatabaseDescriptionFile=NULL,
                            interestGeneFile = NULL, interestGene=genes, 
                            interestGeneType="genesymbol", collapseMethod="mean", 
                            referenceGeneFile=NULL, referenceGene=genes,
                            referenceGeneType="genesymbol", referenceSet=NULL, minNum=5, 
                            maxNum=2000, fdrMethod="BH", sigMethod="top", 
                            fdrThr=0.05, topThr=20, dNum=20,
                            perNum=1000, lNum=20,is.output=TRUE,outputDirectory=getwd(),
                            projectName=paste(projectName,"pathway", sep = "."), 
                            keepGSEAFolder=FALSE,
                            methodType="R",dagColor="binary",
                            hostName="http://www.webgestalt.org/")
    kegg_list[[i]] <- enrichResult_pathway
  }
  # save(go_list, kegg_list, file="web_go_kegg.Rdata")
}

draw_go_kegg_barplot_web <- function(anno_list=go_list) {
  # names(kegg_list)
  for (i in names(anno_list)) {
    barplot_df <- anno_list[[i]]
    if (length(barplot_df)<2) {next}
    barplot_df <- barplot_df[order(barplot_df$PValue, decreasing = F),]
    barplot_df$description <- factor(barplot_df$description, levels=rev(barplot_df$description))
    maxpvalue <- max(-log10(barplot_df$PValue))
    if (length(barplot_df$description)>8) {barplot_df <- barplot_df[1:8,]}
    g <- ggplot(data=barplot_df, aes(x=description, y=-log10(PValue))) +
    geom_bar(stat="identity", aes(fill = description)) +
    geom_text(aes(label=C),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
    ylim(0, maxpvalue*1.1) +
    coord_flip() +
    labs(x = "", y = "-Log10(P-value)", title=i) + 
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 13, color = "black", face = "bold"), axis.text.x = element_text(size = 11, color = "black", face = "plain")) +
    # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),panel.border = element_blank(),
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin=unit(c(0,0,0,0), "cm")) + scale_fill_manual(values=colorRampPalette(c("#3176e0", "#8831e0", "#e031b4"))(10))
    plot(g)
  }
}

go_pathway_by_clusterProfiler_human_bak <- function(geneList=markerList) {
  library(clusterProfiler)
  library(org.Hs.eg.db) # human
  # library(org.Mm.eg.db) # mouse
  go_list <- list()
  kegg_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {nameList <- 1:length(geneList)}
  for (i in nameList) {
    genes <- geneList[[i]]
    projectName <- i
    gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    # gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene      = gene.df$ENTREZID,
                #universe      = genes,
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
    go_list[[i]] <- ego@result
    kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
    # kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
    kegg_list[[i]] <- kk@result
    }
    return(list("go_list"=go_list, "kegg_list"=kegg_list))
  }

ID2gene <- function(ID="4171/4175/5422/4172") {
  ID <- unlist(lapply(strsplit(ID, "/"), as.integer))
  gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = org.Hs.eg.db)
  gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
  rownames(gene.df) <- gene.df$ENTREZID
  genes <- paste(gene.df$SYMBOL, collapse ="/")
  genes
}

go_pathway_by_clusterProfiler_human <- function(geneList=markerList) {
  library(clusterProfiler)
  library(org.Hs.eg.db) # human
  # library(org.Mm.eg.db) # mouse
  go_list <- list()
  kegg_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {nameList <- 1:length(geneList)}
  for (i in nameList) {
    genes <- geneList[[i]]
    projectName <- i
    gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    # gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene      = gene.df$ENTREZID,
                #universe      = genes, # SYMBOL
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
    # remove duplications
    # too slow
    # ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    # drop top 3 level
    # ego <- dropGO(ego, level = c(1,2))
    go_list[[i]] <- ego
    kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
    # kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
    # kk$genes <- unlist(lapply(kk$geneID, ID2gene))
    kegg_list[[i]] <- kk
    }
    return(list("go_list"=go_list, "kegg_list"=kegg_list))
  }

go_pathway_by_clusterProfiler_mouse <- function(geneList=markerList) {
  library(clusterProfiler)
  # library(org.Hs.eg.db) # human
  library(org.Mm.eg.db) # mouse
  go_list <- list()
  kegg_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {nameList <- 1:length(geneList)}
  for (i in nameList) {
    genes <- geneList[[i]]
    projectName <- i
    # gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene      = gene.df$ENTREZID,
                #universe      = genes, # SYMBOL
                keyType       = "ENTREZID",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
    # remove duplications
    # too slow
    # ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    # drop top 3 level
    # ego <- dropGO(ego, level = c(1,2))
    go_list[[i]] <- ego
    # kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
    kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
    kegg_list[[i]] <- kk
    }
    return(list("go_list"=go_list, "kegg_list"=kegg_list))
  }

gsea_by_clusterProfiler_mouse <- function(geneList=DEGs_list_full, pAdjustMethod = "BH", pvalueCutoff = 0.05) {
#   geneList=list(DEGs_list_full[[1]])
  library(clusterProfiler)
  # library(org.Hs.eg.db) # human
  library(org.Mm.eg.db) # mouse
  go_list <- list()
  kegg_list <- list()
  gsea_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {nameList <- 1:length(geneList)}
  for (i in nameList) {
    genes <- geneList[[i]]$gene
    # projectName <- i
    gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
    # prepare geneList
    geneList2 = geneList[[i]][gene.df$SYMBOL,"log2FC"]
    names(geneList2) <- gene.df$ENTREZID
    geneList2 = sort(geneList2, decreasing = TRUE)
    print(length(geneList2))
    # no result, no matter how I try
    ego <- gseGO(geneList     = geneList2,
              OrgDb        = org.Mm.eg.db,
              keyType = "ENTREZID",
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              by = "fgsea", #fgsea, DOSE
              verbose      = F)
    kk <- gseKEGG(geneList     = geneList2,
               organism     = 'mmu', #hsa
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = pvalueCutoff,
               pAdjustMethod = pAdjustMethod,
               by = "fgsea",
               verbose      = F)
    if (nrow(ego@result) > 0) { go_list[[i]] <- ego }
    if (nrow(kk@result) > 0) { kegg_list[[i]] <- kk }
    }
    gsea_list[["go_list"]] <- go_list
    gsea_list[["kegg_list"]] <- kegg_list
    gsea_list
  }  

gsea_by_clusterProfiler_human <- function(geneList=DEGs_list_full, pAdjustMethod = "none", pvalueCutoff = 0.05) {
  # geneList=list(DEG_list_full[[1]])
  library(clusterProfiler)
  library(org.Hs.eg.db) # human
  # library(org.Mm.eg.db) # mouse
  go_list <- list()
  kegg_list <- list()
  gsea_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {nameList <- 1:length(geneList)}
  for (i in nameList) {
    genes <- geneList[[i]]$gene
    if (length(intersect(c("gene", "log2FC"), colnames(geneList[[i]]))) != 2) {
      stop("\ncolnames error, must contain gene and log2FC!!!")
    }
    # projectName <- i
    gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
    # prepare geneList
    geneList2 = geneList[[i]][gene.df$SYMBOL,"log2FC"]
    names(geneList2) <- gene.df$ENTREZID
    geneList2 = sort(geneList2, decreasing = TRUE)
    print(length(geneList2))
    # no result, no matter how I try
    ego <- gseGO(geneList     = geneList2,
              OrgDb        = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              by = "fgsea", #fgsea, DOSE
              verbose      = F)
    kk <- gseKEGG(geneList     = geneList2,
               organism     = 'hsa', #hsa
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = pvalueCutoff,
               pAdjustMethod = pAdjustMethod,
               by = "fgsea",
               verbose      = F)
    if (nrow(ego@result) > 0) { go_list[[i]] <- ego }
    if (nrow(kk@result) > 0) { kegg_list[[i]] <- kk }
    }
    gsea_list[["go_list"]] <- go_list
    gsea_list[["kegg_list"]] <- kegg_list
    gsea_list
  }

draw_go_kegg_barplot_yshu <- function(anno_list=go_list) {
  # anno_list=kegg_list
  # colors <- c("#FB8072", "#B3DE69", "#BC80BD")
  library(ggplot2)
  library(RColorBrewer)
  tmpcolors <- brewer.pal(9,"Set3")[1:length(anno_list)]
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
    #
    if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
    g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
    geom_bar(stat="identity", fill = tmpcolors[j]) +
    geom_text(aes(label=Count),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
    ylim(0, maxpvalue*1.1) +
    coord_flip() +
    labs(x = "", y = "-Log10(P-value)", title=i) + 
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

draw_single_merged_barplot <- function(barplot_df) {
  library(Hmisc)
  # barplot_df <- up_neurongenesis
  # barplot_df <- intermediate
  # barplot_df <- down_stemness
  # barplot_df <- c5
  colors <- brewer.pal(10,"Paired")[5:6]
  for (i in 1:dim(barplot_df)[1]) {
    barplot_df[i,]$Description <- capitalize(as.character(barplot_df[i,]$Description))
  }
  # if (length(barplot_df)<2 | is.null(barplot_df)) {next}
  barplot_df <- barplot_df[order(barplot_df$pvalue, decreasing = F),]
  barplot_df$Description <- factor(barplot_df$Description, levels=rev(barplot_df$Description))
  maxpvalue <- max(-log10(barplot_df$pvalue))
  # if (length(barplot_df$Description)>20) {barplot_df <- barplot_df[1:20,]}
  g <- ggplot(data=barplot_df, aes(x=Description, y=-log10(pvalue))) +
  geom_bar(stat="identity", aes(fill = ptype)) +
  geom_text(aes(label=Count),color="black",vjust=0.4,hjust=-0.5,size=3,fontface="bold") +
  ylim(0, maxpvalue*1.1) +
  coord_flip() +
  labs(x = "", y = "-Log10(P-value)", title="") + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 16, color = "black", face = "plain"), axis.text.x = element_text(size = 12, color = "black", face = "plain"),axis.title =element_text(size = 15)) +
  # axis.title = element_blank() ,plot.title = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.line = element_blank(),panel.border = element_blank(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  plot.margin=unit(c(0,0,0,0), "cm"), panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  scale_fill_manual(values=rev(colors))
  plot(g)
  # save, 7*3 for 4 bars; 4.5 for 6 bars
}

my_get_marker_genes_SC3 <- function(sce=sceTPM, cluster=sceTPM$cellType3, annotation_colors=color3, minNum=100, maxNum=200) {
  markerList <- list()
  markerList2 <- c()
  for (i in names(table(cluster))) {
    simulate_marker <- (cluster == i) * 10
    corM <- cor(simulate_marker, t(logcounts(sce)))
    # corM <- abs(corM)
    corM[corM<0] <- 0
    marker1 <- colnames(corM)[corM>0.4]
    if (length(marker1) < minNum) {
      marker1 <- colnames(corM)[order(corM, decreasing = T)][1:minNum]
    } else if (length(marker1) > maxNum) {
      marker1 <- colnames(corM)[order(corM, decreasing = T)][1:maxNum]
    }
    markerList[[i]] <- marker1
    markerList2 <- c(markerList2, marker1)
  }
  # annotation_colors <-color3
  names(annotation_colors) <- levels(cluster)
  annotation_col <- data.frame(Cluster = factor(cluster), row.names = sce$cellName)
  pheatmap(logcounts(sce)[unique(markerList2),order(cluster)], cluster_cols = F, cluster_rows = T, annotation_col=annotation_col, annotation_colors = list(Cluster=annotation_colors), show_colnames = F, show_rownames=F)
  # simulate_marker <- (sce$allCluster3 %in% c("Cluster1","Cluster2","Cluster3","Cluster4")) * 10
  # markerList[["Cluster1_4"]] <- marker1
  # save(markerList, file = "markerList.Rdata")
  return(markerList)
}

my_get_marker_genes_monocle2 <- function(sce=HSMM, cluster=HSMM$Cluster, minNum=100, maxNum=200) {
  library(WGCNA)
  library(RColorBrewer)
  if (length(table(cluster)) > 8) {
    tmpcolors <- brewer.pal(12,"Set3")[1:length(table(cluster))]
  } else {
    tmpcolors <- brewer.pal(8,"Set2")[1:length(table(cluster))]
  }
  names(tmpcolors) <- levels(cluster)
  markerList <- list()
  markerList2 <- c()
  markerdf <- data.frame()
  exprM <- log2(exprs(sce)+1)
  # for 10x data
  mm10 <- read.table("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/3-10x/Reeson_mar/mm10.feature.info", row.names = 1, header = F)
  rownames(exprM) <- mm10[rownames(exprM),]$V2

  for (i in names(table(cluster))) {
    simulate_marker <- (cluster == i) * 10
    corM <- WGCNA::cor(simulate_marker, t(exprM))
    # corM <- abs(corM)
    corM[corM<0] <- 0
    marker1 <- colnames(corM)[corM>0.4]
    if (length(marker1) < minNum) {
      marker1 <- colnames(corM)[order(corM, decreasing = T)][1:minNum]
    } else if (length(marker1) > maxNum) {
      marker1 <- colnames(corM)[order(corM, decreasing = T)][1:maxNum]
    }
    markerList[[i]] <- marker1
    markerList2 <- c(markerList2, marker1)
  }
  for (i in names(markerList)) {
    for (j in markerList[[i]]) {
      markerdf <- rbind(markerdf, data.frame(marker=j, cluster=i))
    }
  }
  # rownames(markerdf) <- markerdf$marker
  # markerdf[markerdf$marker %in% c("Sox10", "Fabp7", "Plp1", "Elavl4", "Phox2b", "Tubb3"),]
  # annotation_colors <-color3
  annotation_col <- data.frame(Cluster = factor(cluster), row.names = colnames(exprs(HSMM)))
  plotM <- exprM[unique(markerList2),order(cluster)]
  # scale
  plotM2 <- t(scale(t(plotM)))
  plotM2[plotM2<0] <- 0
  plotM2[plotM2>6] <- 6
  #apply(plotM, 1, function(x) {sort(x, decreasing=T)[as.integer(dim(plotM)[2]/2)]})
  pheatmap(plotM2, cluster_cols = F, cluster_rows = F, annotation_col=annotation_col, annotation_colors = list(Cluster=tmpcolors), show_colnames = F, show_rownames=F)
  # simulate_marker <- (sce$allCluster3 %in% c("Cluster1","Cluster2","Cluster3","Cluster4")) * 10
  # markerList[["Cluster1_4"]] <- marker1
  # save(markerList, file = "markerList.Rdata")
  return(markerList)
}

my_get_marker_genes2 <- function(sce, cluster=sce$allCluster3) {
  sce <- sce[,sce$allCluster3!="Cluster9"]
  markerList <- list()
  #for (i in names(table(sce$allCluster3))) {
  simulate_marker <- as.integer(as.integer(sce$allCluster3)<0)
  simulate_marker[(sce$allCluster3 %in% c("Cluster1", "Cluster2", "Cluster3", "Cluster4"))] <- 0
  simulate_marker[(sce$allCluster3 %in% c("Cluster5","Cluster6","Cluster7"))] <- 10
  simulate_marker[(sce$allCluster3 %in% c("Cluster8"))] <- 0
  corM <- cor(simulate_marker, t(logcounts(sce)), method="spearman")
  # corM <- abs(corM)
  corM[corM<0] <- 0
  marker1 <- colnames(corM)[corM>0.4]
  if (length(marker1) < 300) {
    marker1 <- colnames(corM)[order(corM, decreasing = T)][1:300]
  }
  markerList[["up"]] <- marker1
  markerList[["down"]] <- marker1
  markerList[["inter"]] <- marker1
  # }
  # simulate_marker <- (sce$allCluster3 %in% c("Cluster1","Cluster2","Cluster3","Cluster4")) * 10
  # markerList[["Cluster1_4"]] <- marker1
  # save(markerList, file = "up_down_inter_markerList.Rdata")
  allmarkers <- c()
  for (i in 1:8) {allmarkers <- c(allmarkers, markerList[[i]][1:50])}
  annotation_col <- data.frame(Cluster = factor(sce$allCluster3), row.names = sce$cellName)
  # annotation_colors <- brewer.pal(9,"Set3")
  # annotation_colors <- brewer.pal(12,"Set3")[c(1,3:8,10:11)]
  annotation_colors <- c("#a2bf00", "#4fbf00", "#ff8300", "#b800bf","#ef2a33", "#1d9ff9",  "#747ad3", "#e537a2")
  #names(annotation_colors) <- c("Ctrl_ENCC", "SAG4_ENCC","SAG10_ENCC", "ENCC-derived_neurons")
  names(annotation_colors) <- levels(sce$allCluster2)
  pheatmap(logcounts(sce)[unique(allmarkers),order(sce$allCluster3)], cluster_cols = F, cluster_rows = F, annotation_col=annotation_col, annotation_colors = list(Cluster=annotation_colors), show_colnames = F, show_rownames=F)
  # save 600 X 600
  # smmoth curve
  HSMM <- make_CellDataSet(counts(sce)[unique(allmarkers),order(sce$allCluster3)])
  genSmoothCurvesSCE <- genSmoothCurves(HSMM[, ], cores = 1, 
        trend_formula = "~sm.ns(Pseudotime, df=3)", relative_expr = T, new_data = data.frame(Pseudotime=sce$pseudotime, row.names=colnames(sce)))
  genSmoothCurvesSCE <- log10(1+genSmoothCurvesSCE)
  genSmoothCurvesSCE = Matrix::t(scale(Matrix::t(genSmoothCurvesSCE), 
        center = TRUE))
  genSmoothCurvesSCE[genSmoothCurvesSCE>3] <- 3
  genSmoothCurvesSCE[genSmoothCurvesSCE< -3] <- -3
  genSmoothCurvesSCE[is.na(genSmoothCurvesSCE)] <- 0
  colnames(genSmoothCurvesSCE) <- colnames(sce)
  row_dist <- as.dist((1 - cor(Matrix::t(genSmoothCurvesSCE)))/2)
  row_dist[is.na(row_dist)] <- 1
  hclust_method <- "ward.D2"
  pheatmap(genSmoothCurvesSCE, cluster_cols = F, cluster_rows = T, annotation_col=annotation_col, annotation_colors = list(Cluster=annotation_colors), show_colnames = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), clustering_distance_rows = row_dist, clustering_method = hclust_method, useRaster = T, border_color = NA, silent = F)
  # save 600 X 600
}

make_CellDataSet <- function(count_matrix=counts(sce)[unique(allmarkers),order(sce$allCluster3)]) {
  library(monocle)
  pd <- colnames(count_matrix)
  fd <- rownames(count_matrix)
  pd <- data.frame(pd, Pseudotime=sce$pseudotime)
  fd <- data.frame(gene_short_name=fd)
  rownames(pd) <- pd$pd
  rownames(fd) <- fd$gene_short_name
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  # see detail in monocle help doc.
  HSMM <- newCellDataSet(as.matrix(count_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  # Next, use it to estimate RNA counts
  rpc_matrix <- relative2abs(HSMM, method = "num_genes")
  # Now, make a new CellDataSet using the RNA counts
  HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
  # feature selection
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  return(HSMM)
}

sort_cell_of_aCluster_by_distance <- function(exprMatrix=logcounts(sce[rowData(sce)$gene_type=="protein_coding",]), ref=c("Cluster1", "Cluster8"), pcNum=100, threads=3) {
  prin_comp <- prcomp(t(exprMatrix), scale. = T, center = T, rank. = 100)
  pcMatrix <- prin_comp$x
  # something like outlier
  library(parallelDist)
  # dis_matrix <- parDist(x = as.matrix(pcMatrix), method = "euclidean", threads=threads)
  dis_matrix <- parDist(x = as.matrix(pcMatrix), method = "mahalanobis", threads=threads)
  dis_matrix <- as.matrix(dis_matrix)
  dis_matrix[is.na(dis_matrix)] <- 0
  rownames(dis_matrix) <- rownames(pcMatrix)
  colnames(dis_matrix) <- rownames(dis_matrix)
  # sort
  all_order <- c()
  for (i in names(table(sce$allCluster3)[table(sce$allCluster3)!=0])) {
    refc <- ref[!ref%in%i][1]
    obj_ref_matrix <- dis_matrix[sce[,sce$allCluster3==i]$cellName, sce[,sce$allCluster3==refc]$cellName]
    obj_order <- names(sort(rowSums(obj_ref_matrix), decreasing=T))
    all_order <- c(all_order, obj_order)
  }
  return(all_order)
}

deg_analysis <- function(sce) {
  deg_list <- list()
  ctrlGroup <- c("Ctrl_ENCC")
  treatGroup <- c("SAG4_ENCC", "SAG10_ENCC")
  # clusters <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
  clusters <- c("Cluster5")
  ctrl <- apply(logcounts(sce[,sce$cellGroup%in%ctrlGroup & sce$allCluster3 %in% clusters]),1,mean)
  treat <- apply(logcounts(sce[,sce$cellGroup%in%treatGroup & sce$allCluster3 %in% clusters]),1,mean)
  logFC <- sort(treat - ctrl)
  deg_list[["ENCC_up"]] <- names(logFC[logFC>=1])
  deg_list[["ENCC_down"]] <- names(logFC[logFC < -1])
  deg_list[["neuron_up"]] <- names(logFC[logFC>=1])
  deg_list[["neuron_down"]] <- names(logFC[logFC < -1])
  # save(deg_list, file="deg_list.Rdata")
}

## 1. get the gene list of a specific KEGG pathway 
# library(clusterProfiler)
# search_kegg_organism('hsa', by='kegg_code')
# library(KEGGREST)
# query <- keggGet("ko01100")
# TogoWS REST service
get_kegg_genelist <- function(){
  library(EnrichmentBrowser)
  gs <- get.kegg.genesets("hsa")
}

## 2. get the gene list of a specific GO term


## 3. do gene list enrichment analysis


## 4. single cell DEG: SCDE
# cell_group format
# ESC_10 ESC_11
# ESC    ESC
scde_func <- function(count_matrix, cell_group, coreN = 3, fileName="scde.o.ifm.Rdata", providemodel=F, pval=1){
 options(stringsAsFactors = FALSE)
 library(scde)
 #data(es.mef.small)
 #sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
 #names(sg) <- colnames(es.mef.small)  
 #table(sg)
 cd <- clean.counts(count_matrix, min.lib.size=1000, min.reads = 5, min.detected = 1)
 sg <- cell_group
 if (providemodel){
   load(fileName)
 } else {
   o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = coreN, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
   save(o.ifm, file = fileName)
 }
 valid.cells <- o.ifm$corr.a > 0
 o.ifm <- o.ifm[valid.cells, ]
 # estimate gene expression prior
 o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
 groups <- cell_group
 ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  coreN, verbose  =  1)
 p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
 p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
 significant.genes <- which(p.values.adj<=pval)
 ord <- order(p.values.adj[significant.genes]) # order by p-value
 de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
 colnames(de) <- c("Lower_bound","log2FC","Upper_bound","p_value")
 de$log2FC <- -(de$log2FC)
 return(de)
}

## 5. float matrix to interger matrix
df_float_to_int <- function(counts){
  counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  return(counts)
}

## volcano plot
draw_vovalno_plot <- function(data, log2FC_thred=0.5, markedGene, my_title="") {
  colnames(data) <- c("p_value", "log2FC")
  data$sig <- "non-sig"
  data[data$p_value<0.05 & abs(data$log2FC)>=log2FC_thred,]$sig <- "sig"
  # see /Users/surgery/Project/HOME/1-projects/0.bulk_RNA-Seq/9-BACE2/BACE2.R
  library(ggrepel) #Avoid overlapping labels
  # mark_data <- data[!data$sig=="non-sig" | rownames(data) %in% whiteList,]
  mark_data <- data[markedGene,]
  if (dim(mark_data)[1] > 30) {
    mark_data <- mark_data[order(abs(mark_data$log2FC), decreasing = T),][1:30,]
  }
  mark_data$gene <- rownames(mark_data)
  #
  data$color <- data$sig
  data[data$log2FC>0 & data$sig=="sig",]$color <- "up"
  data[data$log2FC<0 & data$sig=="sig",]$color <- "down"
  #
  ggplot(data, aes(x=log2FC, y=-log10(p_value))) +
    geom_hline(aes(yintercept=1), colour="grey50", linetype="dashed", size=0.2) + 
    geom_vline(aes(xintercept=0.5), colour="red", linetype="dashed", size=0.2) + 
    geom_vline(aes(xintercept=-0.5), colour="blue", linetype="dashed", size=0.2) + 
    geom_point(aes(color=color), stroke = 0.3, size=1) +
    scale_color_manual(values=c("blue", "grey50",  "red")) +
    labs(x = "Log2 fold change",y = "-Log10(P-value)", title = my_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size=0.8, colour = "black")) +
    theme(axis.title =element_text(size = 11),axis.text =element_text(size = 10, color = "black"),
          plot.title =element_text(hjust = 0.5, size = 16)) +
    #scale_x_continuous(position="bottom", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    #scale_y_continuous(position="left", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    theme(legend.title=element_blank(), legend.key.size = unit(0.8, 'lines')) +
    geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2, 
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  # for SAG
  # geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.15, point.padding = 0.4, segment.color = 'black', segment.size = 0.3, force = 10, max.iter = 3e3)
}


## 6. mouse gene to human gene, or reverse
# load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/mouse_to_human_homo_gene.Rdata")

## 7. human and mouse TFs
# load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/human_TF_id_name.Rdata")

## 8. draw boxplot
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

## 9. draw violin plot
# source("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/6-redo1934/normalization.R")
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

## 10. draw scatter plot
# source("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/IMR90/IMR90.R)
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

## 11. barplot
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

## 10. get color theme
# https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
# library(RColorBrewer)
# example_col <- rev(brewer.pal(10,"Set3"))
# example_col <- brewer.pal(5,"Set2")

## 11. draw a heatmap for a SC3 object
draw_heatmap_SC3 <- function(tmp_group, cluster=tmp_group$rename_sc3_3, genelist, cluster_rows=F, annotation_colors){
  #tmp_group <- tmp_group[genelist,]
  annotation_col <- data.frame(Cluster = factor(cluster), row.names = tmp_group$cellName)
  order_cells <- colnames(tmp_group)[order(cluster)]
  genelist <- genelist[genelist%in%rownames(tmp_group)]
  # pheatmap(logcounts(deng)[sc3_marker[sc3_marker$cluster=="2cell",]$name, order(deng$cell_type1)], cluster_rows = T, show_rownames = F, show_colnames = F, cluster_cols = F)
  pheatmap(logcounts(tmp_group)[genelist,order_cells],annotation_col = annotation_col, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = cluster_rows, show_rownames = T, annotation_colors=annotation_colors)
}

## 12. get target gene of a TF
get_TF_target_gene <- function(gene){
  TRED[["HMGB2"]]
  ITFP[["HMGB2"]]
  ENCODE[["HMGB2"]]
  Neph2012[["AG10803-DS12374"]][["HMGB2"]]
  TRRUST[["HMGB2"]]
  Marbach2016[["HMGB2"]]
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

## 14. a simple function for mark
qc_mark <- function(x) {if (x) {TRUE} else {FALSE}} 

## 15. calculate correlation index
get_cor_index <- function(){
  library(WGCNA)
  y <- cor(t(x3), method="pearson", nThreads = 3)
}

## 16. clustering by SC3

## 16. clustering by SIMLR
# source("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/IMR90/IMR90.R")
cluster_by_SIMLR <- function(){
  # load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/IMR90/final_plot.Rdata")
  # saveRDS(IMR90_anno_col, file = "/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/IMR90_anno_col.SIMLR.rds")
  library(SIMLR)
  # /Users/surgery/mySIMLR/R/SIMLR_Estimate_Number_of_Clusters.R
  # required external packages for SIMLR
  library(Matrix)
  library(parallel)
  # load the SIMLR R package
  source("/Users/surgery/mySIMLR/R/SIMLR_Estimate_Number_of_Clusters.R")
  source("/Users/surgery/mySIMLR/R/compute.multiple.kernel.R")
  source("/Users/surgery/mySIMLR/R/network.diffusion.R")
  IMR90_simlr = SIMLR(X = logcounts(tmp_group)[know_markers[know_markers%in%rownames(tmp_group)],], c = 6, cores.ratio = 0)
  plot(IMR90_simlr$ydata, col = IMR90_simlr$y$cluster,xlab = "SIMLR component 1", ylab = "SIMLR component 2", pch = 20, main="SIMILR 2D visualization for BuettnerFlorian")
}

## 17. simulate_phenotype
simulate_phenotype <- function(num_cau_SNP=20, num_SNP=500, samplesize=20, h_squared=0.5){
  # generate genotype in Binomial distribution
  pj <- runif(num_SNP, 0.01, 0.5) # probability of success on each trial
  xij_star <- matrix(0, samplesize, num_SNP)
  #for every SNP, select from 0, 1, 2
  for (j in 1: num_SNP) {
    xij_star[,j] <- rbinom(samplesize, 2, pj[j]) } # 2, number of trials
  
  #position of causal SNPs
  CauSNP <- sample(1:num_SNP, num_cau_SNP, replace = F)
  Ord_CauSNP <- sort(CauSNP, decreasing = F)
  
  # generate beta, which is the best predictor
  beta <- rep(0,num_SNP)
  dim(beta) <- c(num_SNP,1)
  # non-null betas follow standard normal distribution
  beta[Ord_CauSNP] <- rnorm(num_cau_SNP,0,1)
  
  # epsilon
  var_e <- sum((xij_star %*% beta)^2) # Multiplies two matrices
  # var_e <- t(beta)%*%t(xij_star)%*%xij_star%*%beta/samplesize*(1-h_squared)/h_squared
  e <- rnorm(samplesize, 0,sqrt(var_e))
  dim(e) <- c(samplesize, 1)
  
  # generate phenotype
  pheno <- xij_star %*% beta + e
  # scale(genotype matrix) # make sure var(u*k) = 1
  return(pheno)
}

## 18. parallel_lapply
parallel_lapply <- function(){
  library(parallel)
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  ####
  pVal_matrix <- do.call(cbind, parallel::parLapply(cl, expr_log, get_pVal, arg1=expr_log))
  stopCluster(cl)
}

# cell_midbrain_markers <- c("TBX1", "ERG", "SOX18", "BCL6B", "SOX17", "TFEB", "FOXS1", "TBX2", "FOXD1", "TBX15", "TBX18", "FHL5", "SPI1", "FOXO3", "IRF8", "IKZF1", "STAT1", "REL", "SOX10", "OLIG1", "ZMAT3", "ETV4", "PTF1A", "PRDM13", "PAX3", "PAX7", "ZSCAN1", "LHX3", "MSX1", "GABPB2", "SAMD13", "ENO1", "GTF2F2", "TOX", "JADE1", "DMBX1", "ZEB2", "TCF3", "HMGB1", "HMGB3", "NR2F6", "MBD4", "GTF2H2", "NKX2-3", "HES6", "TFDP2", "NEUROG1", "DLL1", "ZBTB18", "NEUROD6", "ZNF33A", "PHTF1", "NKX6-2 ", "FOXD2", "SIM1", "ZNF521", "POU4F1", "ZBED4 ", "ZNF91", "ZNF124", "LMO1", "PROX1", "KDM5B ", "ZNF441", "TERF2 ", "NFE2L3", "ZNF57", "ZFP69B", "HMX2", "ZNF300", "ZNF25", "ZNF821", "LMO3", "DEAF1 ", "ZNFX1,ZNF813", "MEIS2", "PRDM6", "TFAP2B", "NFIL3", "ZNF292", "ZNF362", "FOXD4L5", "IKZF3", "ATRX", "CITED2", "PRDM2", "CAMTA1", "FEV", "LHX4", "ISL1", "ESRRB", "ISL2", "PHOX2A", "PHOX2B", "EPAS1", "FOXF2", "ETS1", "FOXC1", "FOXQ1", "ZIC2", "OLIG2", "ETV5 IRX4", "DBX2", "FHL1", "HEY1 HMGA2 NHLH1", "NEUROD4", "NEUROD1", "NEUROG2", "EMX2 PITX3", "EN1", "PBX1", "ZNF264,ZNF248", "MEF2C", "JUNB", "BCL3", "NR4A1", "IRF1", "CSRNP1", "GLIS3", "ID4", "KLF15", "NACC2", "HMGB2", "HMGA1 SOX4 TUB", "LMO4", "ETV1 MKX", "SALL2 OTX2 NHLH2 PAX5", "BCL11B", "NR4A3", "STAT3", "IFI35", "WWTR1 POU3F1 GATA3", "HES5", "SOX9 TEAD2 ONECUT2", "ZNF608", "FOXA2")

## 19. draw vocalco plot
draw_vovalno_plot <- function(data, log2FC_thred=0.5, markedGene, my_title="") {
  colnames(data) <- c("p_value", "log2FC")
  data$sig <- "non-sig"
  data[data$p_value<0.05 & abs(data$log2FC)>=log2FC_thred,]$sig <- "sig"
  # see /Users/surgery/Project/HOME/1-projects/0.bulk_RNA-Seq/9-BACE2/BACE2.R
  library(ggrepel) #Avoid overlapping labels
  # mark_data <- data[!data$sig=="non-sig" | rownames(data) %in% whiteList,]
  mark_data <- data[markedGene,]
  if (dim(mark_data)[1] > 30) {
    mark_data <- mark_data[order(abs(mark_data$log2FC), decreasing = T),][1:30,]
  }
  mark_data$gene <- rownames(mark_data)
  #
  data$color <- data$sig
  data[data$log2FC>0 & data$sig=="sig",]$color <- "up"
  data[data$log2FC<0 & data$sig=="sig",]$color <- "down"
  #
  ggplot(data, aes(x=log2FC, y=-log10(p_value))) +
    geom_hline(aes(yintercept=1), colour="grey50", linetype="dashed", size=0.2) + 
    geom_vline(aes(xintercept=0.5), colour="red", linetype="dashed", size=0.2) + 
    geom_vline(aes(xintercept=-0.5), colour="blue", linetype="dashed", size=0.2) + 
    geom_point(aes(color=color), stroke = 0.3, size=1) +
    scale_color_manual(values=c("blue", "grey50",  "red")) +
    labs(x = "Log2 fold change",y = "-Log10(P-value)", title = my_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(size=0.8, colour = "black")) +
    theme(axis.title =element_text(size = 11),axis.text =element_text(size = 10, color = "black"),
          plot.title =element_text(hjust = 0.5, size = 16)) +
    #scale_x_continuous(position="bottom", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    #scale_y_continuous(position="left", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    theme(legend.title=element_blank(), legend.key.size = unit(0.8, 'lines')) +
    geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2, 
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  # for SAG
  # geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.15, point.padding = 0.4, segment.color = 'black', segment.size = 0.3, force = 10, max.iter = 3e3)
}


## 20. fit maturation trajectory
maturation.trajectory <- function(cm, md, expr, pricu.f=1/3, ref_gene="NANOG") {
  cat('Fitting maturation trajectory\n')
  genes <- apply(cm[rownames(expr), ] > 0, 1, mean) >= 0.02 & apply(cm[rownames(expr), ] > 0, 1, sum) >= 3
  rd <- dim.red(expr[genes, ], max.dim=50, ev.red.th=0.04)
  # for a consisten look use ref_gene expression to orient each axis
  for (i in 1:ncol(rd)) {
    if (cor(as.numeric(expr[ref_gene, ]), rd[, i]) > 0) {
      rd[, i] <- -rd[, i]
    }
  }
  
  md <- md[, !grepl('^DMC', colnames(md))]
  md <- cbind(md, rd)
  
  pricu <- principal.curve(rd, smoother='lowess', trace=TRUE, f=pricu.f, stretch=333)
  # two DMCs
  pc.line <- as.data.frame(pricu$s[order(pricu$lambda), ])
  # lambda, for each point, its arc-length from the beginning of the curve. The curve is parametrized approximately by arc-length, and hence is unit-speed.
  md$maturation.score <- pricu$lambda/max(pricu$lambda)
  
  # orient maturation score using ref_gene expression
  if (cor(md$maturation.score, as.numeric(expr[ref_gene, ])) > 0) {
    md$maturation.score <- -(md$maturation.score - max(md$maturation.score))
  }
  
  # use 1% of neighbor cells to smooth maturation score
  md$maturation.score.smooth <- nn.smooth(md$maturation.score, rd[, 1:2], round(ncol(expr)*0.01, 0))
  
  # pick maturation score cutoff to separate mitotic from post-mitotic cells
  md$in.cc.phase <- md$cc.phase != 0
  fit <- loess(as.numeric(md$in.cc.phase) ~ md$maturation.score.smooth, span=0.5, degree=2)
  md$cc.phase.fit <- fit$fitted
  # pick MT threshold based on drop in cc.phase cells
  # ignore edges of MT because of potential outliers
  mt.th <- max(subset(md, md$cc.phase.fit > mean(md$in.cc.phase)/2 & md$maturation.score.smooth >= 0.2 & md$maturation.score.smooth <= 0.8)$maturation.score.smooth)
  
  md$postmitotic <- md$maturation.score.smooth > mt.th
  return(list(md=md, pricu=pricu, pc.line=pc.line, mt.th=mt.th))
}


# for smoothing maturation score
nn.smooth <- function(y, coords, k) {
  knn.out <- FNN::get.knn(coords, k)
  w <- 1 / (knn.out$nn.dist+.Machine$double.eps)
  w <- w / apply(w, 1, sum)
  v <- apply(knn.out$nn.index, 2, function(i) y[i])
  return(apply(v*w, 1, sum))
}

# maturation score colors
my.cols.RYG <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b",
                                  "#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"))(11)
##
scRNA.seq.course <- function() {
  setwd("/Users/surgery/Project/HOME/github/scRNA.seq.course")
  #devtools::install_github("hemberg-lab/scRNA.seq.funcs")
  library(scRNA.seq.funcs)
  library(edgeR)
  library(monocle)
  library(MAST)
  library(ROCR)
  set.seed(1)

}

#
myDotPlotM <- function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
    col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.order = c(),
    scale.by = "radius", scale.min = NA, scale.max = NA, group.by, title="",
    plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
    library(dplyr)
    library(tidyr)
    PercentAbove <- function(x, threshold){
      return(length(x = x[x > threshold]) / length(x = x))
    }
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    if (!missing(x = group.by)) {
        object <- SetAllIdent(object = object, id = group.by)
    }
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
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
          axis.text.x  = element_text(face="plain", angle=0, size = 12, color = "black", vjust=0.5),
          axis.text.y  = element_text(face="italic", size = 12, color = "black")) +
    labs(title = title) +
    guides(colour = guide_colourbar(title = "Relative expression", title.position = "left",
                              direction="vertical", title.theme = element_text(angle=90), title.hjust=-4),
           size = guide_legend(title = "Expressed percentage (%)", title.position = "left",
                              direction="vertical", title.theme = element_text(angle=90), title.hjust=-15))
    #
    if (length(x = cols.use) == 1) {
        p <- p + scale_color_distiller(palette = cols.use)
    }
    else {
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

myDotPlotH <- function (sce, genes.plot, use="cellGroup2",xAngle=60, cols.use = c("lightgrey", "blue"), 
    col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.order = c(),
    scale.by = "radius", scale.min = NA, scale.max = NA, group.by, title="",
    plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
#     cols.use = c("lightgrey", "blue"); genes.plot = tmp$human; sce <- sce_comp; use="cellGroup2";
#     col.min = -2.5; col.max = 2.5; dot.min = 0; dot.scale = 6; group.order = c();xAngle=60;
#     scale.by = "radius"; scale.min = 0; scale.max = 100; title="";
#     plot.legend = T; do.return = FALSE; x.lab.rot = FALSE
    #
    library(dplyr)
    library(tidyr)
    PercentAbove <- function(x, threshold){
      return(length(x = x[x > threshold]) / length(x = x))
    }
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
#     if (!missing(x = group.by)) {
#         object <- SetAllIdent(object = object, id = group.by)
#     }

#     data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    data.to.plot <- as.data.frame(t(logcounts(sce)[genes.plot,]))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
#     data.to.plot$id <- object@ident
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
    #
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

scale_count_seurat <- function(seuset, median=F) {
    mouse_exprM <- as.matrix(t(seuset@raw.data))
    if (median==F) {median_c <- 10000}
    else {median_c <- median(apply(seuset@raw.data,2,sum))}
    mouse_exprM_scale <- mouse_exprM/rowSums(mouse_exprM)*median_c
    mouse_exprM_scale <- log2(1+mouse_exprM_scale)
    mouse_exprM_scale <- as.data.frame(as.matrix(t(mouse_exprM_scale)))
    mouse_exprM_scale
}

get_unique_markers <- function(countM=counts(sce_comp), colanno=pca_data, use="colour_by", max_pct=0.5, min_pct=0.1) {
    #countM=seuset@raw.data; colanno=all_tsne; use="cluster"; max_pct=0.5; min_pct=0.1
    # prepare
    library(gtools)
    countM <- countM[, rownames(colanno)]
    colanno$c <- colanno[,use]
    clusters <- unique(as.character(colanno$c))
    n <- length(clusters)
    # get pct of each cluster first
    # will simplify the calculation
    obj_df <- data.frame()
    for (i in clusters) {
        obj_pct <- rowSums(countM[,colanno$c %in% i]>0)/sum(colanno$c %in% i)
        obj_df <- rbind(obj_df, obj_pct)
        }
    rownames(obj_df) <- clusters
    colnames(obj_df) <- rownames(countM)

    obj_df <- obj_df[,!(colSums(obj_df>max_pct)==length(clusters) | colSums(obj_df<min_pct)==length(clusters))]

    uniq_markers=data.frame()
    # only need to compare paired cluster
    # will do further merge
    comb_df <- permutations(n, 2)
    for (j in 1:dim(comb_df)[1]) {
        obj_c <- clusters[comb_df[j,1]]
        bcg_c <- clusters[comb_df[j,2]]
        obj_pct <- obj_df[obj_c,]
        bcg_pct <- obj_df[bcg_c,]
        tmp_uniq_markers <- names(obj_pct)[obj_pct>max_pct & bcg_pct<min_pct]
        if(length(tmp_uniq_markers)>=1) {
            uniq_markers <- rbind(uniq_markers, data.frame(marker=tmp_uniq_markers, obj=obj_c, bcg=bcg_c,
                               obj_pct=as.vector(t(obj_pct[obj_c,tmp_uniq_markers])), 
                               bcg_pct=as.vector(t(bcg_pct[bcg_c,tmp_uniq_markers]))))
        }
    }
    list(uniq_markers=uniq_markers, obj_df=obj_df)
}

staing_plot <- function(sce=sce_dev, colanno=pca_full$data, use=c("X","Y"), show_genes=ENSgene, maxV=4, minV=0) {
    colanno_expr <- cbind(colanno[,use], t(logcounts(sce)[show_genes,rownames(colanno)]))
    colanno_expr_melt <- melt(colanno_expr, id.vars = use)
    # print(colnames(colanno_expr_melt))
    colnames(colanno_expr_melt) <- c("x","y","gene","expression")
    # colanno_expr_melt[colanno_expr_melt$expression < minV,]$expression <- minV
    colanno_expr_melt[colanno_expr_melt$expression > maxV,]$expression <- maxV
    #
    ggplot(colanno_expr_melt, aes(x=x, y=y, color=expression)) + 
        geom_point(size=0.2, alpha=0.7) + 
      facet_wrap( ~ gene, ncol=4) +
      #geom_point(size=2, alpha=0.5) + 
      theme_bw() + labs(x='Dimension1', y='Dimension2') + 
      geom_density_2d(color='black', size=0.1, alpha=0.2) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
      theme(strip.background = element_blank(),
            panel.border = element_rect(colour = "black"))+
      guides(colour = guide_colourbar(title = "expression", title.position = "left",
                                  direction="vertical", title.theme = element_text(angle=90))) +
      # theme(strip.text = element_text(colour = 'black')) +
      scale_colour_gradient(low = "lightgrey", high = "blue")
}

draw_cell_scatter_plot <- function(colanno=pca_data, use=c("X","Y"), colour_by="colour_by") {
    colanno$x <- colanno[,use[1]]
    colanno$y <- colanno[,use[2]]
    colanno$color <- colanno[,colour_by]
    n <- length(unique(colanno$color))
    ggplot(colanno, aes(x=x, y=y, color=color)) + 
      #facet_wrap( ~ shape_by, ncol=3) +
      geom_point(size=1, alpha=1) + 
      labs(x = "Dimension1",y = "Dimension2", title = "") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.title=element_blank()) +
        theme(strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.2))+
        theme(#axis.ticks.x = element_blank(),
            axis.text.x  = element_text(face="plain", size = 7, color = "black"),
            axis.text.y  = element_text(face="plain", size = 7, color = "black"),
            axis.title =element_text(size = 12)) +
        theme(strip.background = element_rect(fill = "gray97", color = NA, size = 0.2)) +
        theme(strip.placement = "outside", strip.text = element_text(face="bold", size = 11))+
      scale_color_manual(values=brewer.pal(n,"Set2"))
}

scale_count_seurat <- function(seuset, median=F) {
    seuset@raw.data <- seuset@raw.data[,colnames(seuset@scale.data)]
    mouse_exprM <- as.matrix(t(seuset@raw.data))
    if (median==F) {median_c <- 10000}
    else {median_c <- median(apply(seuset@raw.data,2,sum))}
    mouse_exprM_scale <- mouse_exprM/rowSums(mouse_exprM)*median_c
    mouse_exprM_scale <- log2(1+mouse_exprM_scale)
    mouse_exprM_scale <- as.data.frame(as.matrix(t(mouse_exprM_scale)))
    mouse_exprM_scale
}

supervised_clustering <- function(exprM_to_be_c, exprM_ctrl, annoDf) {
  cells <- intersect(colnames(exprM_ctrl), rownames(annoDf))
  annoDf <- annoDf[cells,]
  options(stringsAsFactors = F)
  exprM_ctrl <- exprM_ctrl[,cells]
  #
  centerM <- data.frame()
  for (cluster in unique(annoDf$cluster)) {
    cells <- rownames(annoDf[annoDf$cluster == cluster,])
    centerM <- rbind(centerM, t(data.frame(apply(exprM_ctrl[,cells],1,mean))))
  }
  rownames(centerM) <- unique(annoDf$cluster)
  #
  library(pdist)
  # better to do a dimention reduction
  distM_merged <- pdist(t(exprM_to_be_c), centerM)
  distM_merged <- as.matrix(distM_merged)
  #distM_merged <- parDist(x = as.matrix(t(exprM_merged)), method = "euclidean", threads=3)
  rownames(distM_merged) <- colnames(exprM_to_be_c)
  colnames(distM_merged) <- rownames(centerM)
  # assign cluster for the new cells
  new_annoDf <- as.data.frame(apply(distM_merged,1,which.min))
  new_annoDf$cluster <- colnames(distM_merged)[new_annoDf[,1]]
  new_annoDf$cellGroup <- substr(rownames(new_annoDf),1,3)
  new_annoDf[,1] <- NULL
  new_annoDf
}

# input: genestr list
# output: unique gene list
str2genes <- function(genestrlist) {
    unique(unlist(lapply(genestrlist, function(x) {unlist(strsplit(x, "/"))})))
}

# transfer ENTREZID ID to gene name, for KEGG and GSEA analysis
# input: ID <- "4171/4175/5422/4172"
# output: gene <- "PFN1/ITGA9/EGFR/CFL1
ID2gene <- function(ID) {
    #ID <- ="4171/4175/5422/4172"
    ID <- unlist(lapply(strsplit(ID, "/"), as.integer))
    gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = org.Mm.eg.db)
    gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
    rownames(gene.df) <- gene.df$ENTREZID
    #genes <- paste(gene.df$SYMBOL, collapse ="/")
    #genes
    gene.df$SYMBOL
}

# transfer gene ID to postion in pseudotime, case-control study
gene2pos <- function(genestr, geneorder=rownames(raw_log2fc_df)) {
    #genestr="MCM2/MCM6/POLA1/MCM3"; geneorder=rownames(raw_log2fc_df)
    genes <- unlist(strsplit(genestr, "/"))
    geneorderNum <- which(geneorder %in% genes)
    geneorderstr <- paste(round(sort(geneorderNum)/length(geneorder), 2), collapse = "/")
    geneorderstr
}

# get GO term gene list
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#
get_GO_data <- function(OrgDb, ont, keytype) {
    GO_Env <- get_GO_Env()
    use_cached <- FALSE

    if (exists("organism", envir=GO_Env, inherits=FALSE) &&
        exists("keytype", envir=GO_Env, inherits=FALSE)) {

        org <- get("organism", envir=GO_Env)
        kt <- get("keytype", envir=GO_Env)

        if (org == DOSE:::get_organism(OrgDb) &&
            keytype == kt &&
            exists("goAnno", envir=GO_Env, inherits=FALSE)) {
            ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
            ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){

            use_cached <- TRUE
        }
    }

    if (use_cached) {
        goAnno <- get("goAnno", envir=GO_Env)
    } else {
        OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
        kt <- keytypes(OrgDb)
        if (! keytype %in% kt) {
            stop("keytype is not supported...")
        }

        kk <- keys(OrgDb, keytype=keytype)
        goAnno <- suppressMessages(
            select(OrgDb, keys=kk, keytype=keytype,
                   columns=c("GOALL", "ONTOLOGYALL")))

        goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])

        assign("goAnno", goAnno, envir=GO_Env)
        assign("keytype", keytype, envir=GO_Env)
        assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
    }

    if (ont == "ALL") {
        GO2GENE <- unique(goAnno[, c(2,1)])
    } else {
        GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
    }

    GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())

    goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
    goOnt <- goOnt.df[,2]
    names(goOnt) <- goOnt.df[,1]
    assign("GO2ONT", goOnt, envir=GO_DATA)
    return(GO_DATA)
}

get_GO_Env <- function () {
    if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos)
        assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
    }
    get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
    pos <- 1
    envir <- as.environment(pos)
    if (!exists(".GOTERM_Env", envir=envir)) {
        assign(".GOTERM_Env", new.env(), envir)
    }
    GOTERM_Env <- get(".GOTERM_Env", envir = envir)
    if (exists("GOTERM.df", envir = GOTERM_Env)) {
        GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
    } else {
        GOTERM.df <- toTable(GOTERM)
        assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
    }
    return(GOTERM.df)
}
