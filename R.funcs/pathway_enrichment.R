# source("https://github.com/leezx/bt2m/raw/main/notebooks/unpackaged-code/pathway_enrichment.R")

#' GO KEGG ORA analysis
#' @param geneList a gene list
#' @param organism organism (hs and mm)
#'
#' @return a list with GO and KEGG annotation result (clusterProfiler format)
#' @export
#' @examples
#' result <- ora.go.kegg.clusterProfiler(geneList = new_moduleList, organism="mm")
#'
ora.go.kegg.clusterProfiler <- function(geneList=markerList, organism="hs") {
  library(clusterProfiler)
  go_list <- list()
  kegg_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {
    print("no name for the gene list!!!\n")
    nameList <- 1:length(geneList)
  }
  for (i in nameList) {
    genes <- geneList[[i]]
    projectName <- i
    if (organism=="mm") {
      library(org.Mm.eg.db) # mouse
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
      # remove duplications # too slow
      # ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      # drop top 3 level
      # ego <- dropGO(ego, level = c(1,2))
      go_list[[i]] <- ego
      kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
      # kk$genes <- unlist(lapply(kk$geneID, ID2gene))
      kegg_list[[i]] <- kk
    }
    else if (organism=="hs") {
      library(org.Hs.eg.db)
      gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
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
      # kk$genes <- unlist(lapply(kk$geneID, ID2gene))
      kegg_list[[i]] <- kk
    }
    else {stop("only support hs and mm now!")}
  }
  return(list("go_list"=go_list, "kegg_list"=kegg_list))
}


#' GSEA analysis by newest fgsea
#' @param geneList a gene list
#' @param use.score the score used for the GSEA analysis
#' @param organism organism (hs and mm)
#' @param pvalueCutoff pvalueCutoff
#'
#' @return a list with GO and KEGG GSEA annotation result (clusterProfiler format)
#' @export
#' @examples
#' gsea_list <- gsea.go.kegg.clusterProfiler(geneList = pheno_DEGs, use.score = "cor", organism="mm")
#'
gsea.go.kegg.clusterProfiler <- function(geneList=DEGs_list_full, use.score="cor", organism="hs", pvalueCutoff = 1) {
  pAdjustMethod = "BH";
  library(clusterProfiler)
  library(ReactomePA)
  go_list <- list()
  kegg_list <- list()
  gsea_list <- list()
  reactome_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {
    print("No name for the gene list!!! Will use 1:n\n")
    nameList <- 1:length(geneList)
  }
  for (i in nameList) {
    # genes <- geneList[[i]]$gene
    genes <- rownames(geneList[[i]])
    # projectName <- i
    if (organism=="mm") {
      library(org.Mm.eg.db) # mouse
      gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
      gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
      # prepare geneList
      geneList2 = geneList[[i]][gene.df$SYMBOL,use.score]
      names(geneList2) <- gene.df$ENTREZID
      geneList2 = sort(geneList2, decreasing = TRUE)
      print(length(geneList2))
      # no result, no matter how I try
      ego <- gseGO(geneList     = geneList2,
                   OrgDb        = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont          = "BP",
                   # nPerm        = 1000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = pvalueCutoff,
                   # pAdjustMethod = pAdjustMethod,
                   # by = "fgsea", #fgsea, DOSE
                   verbose      = F)
      # gseGO(geneRank, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10,
      # maxGSSize = 1000, pvalueCutoff=1)
      kk <- gseKEGG(geneList     = geneList2,
                    organism     = 'mmu', #hsa
                    # nPerm        = 1000,
                    minGSSize    = 10,
                    maxGSSize = 1000,
                    pvalueCutoff = pvalueCutoff,
                    # pAdjustMethod = pAdjustMethod,
                    # by = "fgsea",
                    verbose      = F)
      # gseKEGG(geneRank, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
      Reactomep <- ReactomePA::gsePathway(geneList2,
                              organism = "mouse",
                              # nPerm = 1000,
                              minGSSize = 10,
                              maxGSSize = 1000,
                              pvalueCutoff = pvalueCutoff)
    } else if (organism=="hs") {
      library(org.Hs.eg.db) # human
      gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
      gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
      # prepare geneList
      geneList2 = geneList[[i]][gene.df$SYMBOL,use.score]
      names(geneList2) <- gene.df$ENTREZID
      geneList2 = sort(geneList2, decreasing = TRUE)
      print(length(geneList2))
      # no result, no matter how I try
      ego <- gseGO(geneList     = geneList2,
                   OrgDb        = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont          = "BP",
                   # nPerm        = 1000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = pvalueCutoff,
                   # pAdjustMethod = pAdjustMethod,
                   # by = "fgsea", #fgsea, DOSE
                   verbose      = F)
      kk <- gseKEGG(geneList     = geneList2,
                    organism     = 'hsa', #hsa
                    # nPerm        = 1000,
                    minGSSize    = 10,
                    maxGSSize = 1000,
                    pvalueCutoff = pvalueCutoff,
                    # pAdjustMethod = pAdjustMethod,
                    # by = "fgsea",
                    verbose      = F)
      Reactomep <- ReactomePA::gsePathway(geneList2,
                              organism = "human",
                              # nPerm = 1000,
                              minGSSize = 10,
                              maxGSSize = 1000,
                              pvalueCutoff = pvalueCutoff)
    } else {
      stop("only support hs and mm now!!!")
    }
    if (nrow(ego@result) > 0) { go_list[[i]] <- ego }
    if (nrow(kk@result) > 0) { kegg_list[[i]] <- kk }
    if (nrow(Reactomep@result) > 0) { reactome_list[[i]] <- Reactomep }
  }
  gsea_list[["go_list"]] <- go_list
  gsea_list[["kegg_list"]] <- kegg_list
  gsea_list[["reactome_list"]] <- reactome_list
  gsea_list
}


# get genes from specific GO terms
get_GO_data <- function(OrgDb, ont, keytype) {
  library(GO.db)
  library(DOSE)
  library(GOSemSim)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  #
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
    GOTERM.df <- toTable(GO.db::GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}

# do enrichment analysis manually
simple.enrichment.analysis <- function(k,n,M,N) {
  # see clusterProfiler
  # GeneRatio：k/n [n is the number of query genes, k is the genes in this pathways]
  # BgRatio: M/N [M is the number of this pathway, N is total background genes]
  args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                        numW=M,        ## White balls
                        numB=N-M,      ## Black balls
                        numDrawn=n)    ## balls drawn
  ##
  ## calcute pvalues based on hypergeometric model
  pvalues <- apply(args.df, 1, function(n)
    phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
  )
  pvalues
}

# 1. give a gene list
# 2. see which GO it enriched most in a dataframe from clusterProfiler
genes_enrich_GO.result <- function(querys=querys, tmp.enrich=tmp.enrich) {
  # get the genes
  tmp.enrich.genes <- unique(unlist(lapply(tmp.enrich$geneID, function(x) {
    strsplit(x, split = "/")
  })))

  total.bcg <- length(tmp.enrich.genes)
  # total.bcg
  querys <- querys[querys %in% tmp.enrich.genes]
  total.query <- length(querys)
  # total.query

  tmp.enrich$enrich.p <- 1
  tmp.enrich$k <- 0
  tmp.enrich$query <- ""

  # see the manual in cnblog
  for (i in 1:nrow(tmp.enrich)) {
    genes.m <- tmp.enrich[i,"geneID"]
    genes <- unique(unlist(lapply(genes.m, function(x) {
      strsplit(x, split = "/")
    })))
    my.genes <- genes[genes %in% querys]
    k <- length(my.genes)
    M <- total.query
    n <- length(genes)
    N <- total.bcg
    tmp.enrich[i,]$k <- k
    tmp.enrich[i,]$enrich.p <- phyper(k-1, M, N-M, n, lower.tail = F)
    tmp.enrich[i,]$query <- paste(my.genes, collapse = "/")
    # break
  }

  # sort and filter
  tmp.enrich <- tmp.enrich[order(tmp.enrich$enrich.p, decreasing = F),]
  tmp.enrich <- subset(tmp.enrich, k>0)
  return(tmp.enrich)
}

# get genes in dataframe from clusterProfiler result
geneIDtoGeneList <- function(p1=p1) {
  core_gene_set <- unique(unlist(lapply(p1$geneID, function(x) {
    strsplit(x, split = "/")[[1]]
  })))
  print(length(core_gene_set))
  return(core_gene_set)
}

# a simple function to get the genes in clusterProfiler result
# get.gene.list(neuron.pathways)
get.gene.list <- function(clusterProfiler.result) {
  unique(unlist(lapply(clusterProfiler.result$geneID, function(x) {
    strsplit(x, split = "/")
  })))
}

#' transfer the ID string to gene vector or string
#' @param ID the ID string, like 4171/4175/5422/4172
#' @param organism organism (hs and mm)
#' @param returnVector return vector or string
#'
#' @return a transfromed gene name list or string
#' @export
#' @examples
#' tmpgenes <- gsea.ID2gene("4171/4175/5422/4172", organism="hs")
#'
gsea.ID2gene <- function(ID, organism="hs", returnVector=T) {
  ID <- unlist(lapply(strsplit(ID, "/"), as.integer))
  if (organism=="hs") {
    library(org.Hs.eg.db)
    gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = org.Hs.eg.db)
  } else if (organism=="mm") {
    library(org.Mm.eg.db)
    gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = org.Mm.eg.db)
  } else {stop("only support hs and mm now!")}

  gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
  rownames(gene.df) <- gene.df$ENTREZID
  if (returnVector) {
    return(gene.df$SYMBOL)
  } else {
    genes <- paste(gene.df$SYMBOL, collapse ="/")
    genes
  }
}


#' Get the genes of interested GO terms by some key words
#' @param organism organism (hs and mm)
#' @param keyWords keyWords (for grey function, e.g. cell cycle|DNA replication)
#' @param returnAll if return all the GO terms
#'
#' @return a list with the GO terms and genes
#' @export
#' @examples
#' tmpList <- toC_get_genes_of_GO_by_keyWords(organism = "hs", keyWords = "axon", returnAll = F)
#' tmpList <- toC_get_genes_of_GO_by_keyWords(organism = "hs",
#'                                            keyWords = "cell cycle|DNA replication|cell division|segregation")
#'
gsea.GO.genes.keyWords <- function(organism="hs", keyWords="axon", returnAll=F) {
  ######################
  # functions from clusterProfiler
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
    library(dplyr)
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Term")] %>% unique
  }

  get_GOTERM <- function() {
    library(GO.db) # GO.db::GOTERM
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
  #######################
  # get_GO_data <- rvcheck::get_fun_from_pkg("clusterProfiler", "get_GO_data") ## for easy installation
  if (organism=="hs") {
    GO_DATA <- get_GO_data("org.Hs.eg.db", "BP", "SYMBOL") # org.Mm.eg.db
  } else if (organism=="mm") {
    GO_DATA <- get_GO_data("org.Mm.eg.db", "BP", "SYMBOL") # org.Mm.eg.db
  } else {
    stop("currently only support hs and mm!")
  }
  # EXTID2PATHID: gene->GO terms
  # GO2ONT: type of GO terms
  # PATHID2EXTID: GO term->genes
  # PATHID2NAME: names of GO terms
  # weird thing: PATHID2EXTID != PATHID2NAME
  #######################
  returnList <- list()
  if (returnAll==T) {
    GO_DATA
  } else {
    # keyWords <- "cell cycle|DNA replication|cell division|segregation"
    tmpGO <- GO_DATA$PATHID2NAME[grep(keyWords, GO_DATA$PATHID2NAME)]
    returnList[["PATHID2NAME"]] <- tmpGO[names(tmpGO) %in% names(GO_DATA$PATHID2EXTID)]
    # the length not matched
    returnList[["PATHID2EXTID"]] <- GO_DATA$PATHID2EXTID[names(returnList[["PATHID2NAME"]])]
    returnList
  }
}

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
  tmpcolors <- c(brewer.pal(12,"Set3"), brewer.pal(9,"Set1"))[1:length(anno_list)]
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
plot.ora.GO.KEGG.barplot.batch <- function(anno_list=go_list, type, f.length=40, rmDup=F, rmPattern="viral") {
  # anno_list=kegg_list
  # colors <- c("#FB8072", "#B3DE69", "#BC80BD")
  library(ggplot2)
  library(RColorBrewer)
  tmpcolors <- c(brewer.pal(12,"Set3"), brewer.pal(9,"Set1"))[1:length(anno_list)]
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
    #####################
    # filter duplication
    if (rmDup) {
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
    }
    #####################
    # remove terms with too long description
    barplot_df <- barplot_df[sapply(as.character(barplot_df$Description), nchar) < f.length,]
    #####################
    # remove patterns
    barplot_df <- subset(barplot_df, !grepl(rmPattern, Description))
    # plotting
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

plot.GO.barplot.pair <- function(barplot_df, colors=1:2, width=35) {
  library(Hmisc)
  library(stringr)
  library(RColorBrewer)
  #
  if (length(colors)!=2) stop("Please input 1:2 or 3:4, etc..")
  tmp.colors <- brewer.pal(12, "Paired")[colors]
  #
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
    geom_bar(stat="identity", aes(color=type, fill=type), alpha=0.8) +
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
    scale_x_discrete(labels=function(x) str_wrap(x, width=width)) +
    scale_fill_manual(values=rev(tmp.colors)) +
    scale_color_manual(values=rev(tmp.colors))
  g
}

ID2gene <- function(ID="4171/4175/5422/4172") {
  ID <- unlist(lapply(strsplit(ID, "/"), as.integer))
  gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = org.Hs.eg.db)
  gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
  rownames(gene.df) <- gene.df$ENTREZID
  genes <- paste(gene.df$SYMBOL, collapse ="/")
  genes
}

