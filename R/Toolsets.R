#' description of the example function
#' @param param1 param1
#'
#' @return the return object
#' @export
#' @examples
#' example.function()
#'
example.function <- function(param1=NULL) {
	NULL
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

#' prepare the annotation for the pheatmap
#' @param sce sce object
#' @param group.order group.order
#' @param sample.order sample.order
#' @param Cluster the feature take as Cluster
#' @param Sample the feature take as Sample
#'
#' @return anno list
#' @export
#' @examples
#' options(repr.plot.width=8, repr.plot.height=7)
#' pheatmap(scale.data[sc3_marker$name, rownames(anno$col)], cluster_rows = F, cluster_cols = F, border_color = NA,
#'          color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(100),
#'          show_colnames=F, show_rownames=F, annotation_col = anno$col, annotation_colors = anno$colors)
#'
pre.pheatmap.anno <- function(sce=sce, group.order=group.order, sample.order=sample.order, 
                              Cluster="impute_cluster", Sample="cellGroup") {
    anno <- list()
    #######
    annotation_col <- as.data.frame(colData(sce)[,c(Cluster, Sample)])
    colnames(annotation_col) <- c("Cluster", "Sample")
    annotation_col$Cluster <- factor(annotation_col$Cluster, levels = group.order)
    annotation_col$Sample <- factor(annotation_col$Sample, levels = sample.order)
    annotation_col <- annotation_col[order(annotation_col$Cluster, decreasing = F),]
    #######
    annotation_colors <- list()
    annotation_colors[["Cluster"]] <- brewer.pal(8,"Set2")[1:length(levels(annotation_col$Cluster))]
    annotation_colors[["Sample"]] <- brewer.pal(12,"Set3")[1:length(levels(annotation_col$Sample))]
    names(annotation_colors$Cluster) <- group.order
    names(annotation_colors$Sample) <- sample.order
    #######
    anno[["col"]] <- annotation_col
    anno[["colors"]] <- annotation_colors
    return(anno)
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
gsea.go.kegg.clusterProfiler <- function(geneList=DEGs_list_full, use.score="cor", organism="hs", pvalueCutoff = 0.05) {
  pAdjustMethod = "BH"; 
  library(clusterProfiler)
  go_list <- list()
  kegg_list <- list()
  gsea_list <- list()
  nameList <- names(geneList)
  if (is.null(nameList)) {
		print("no name for the gene list!!!\n") 
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
			}
		else if (organism=="hs") {
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
		}
		else {
			stop("only support hs and mm now!")
		}
    if (nrow(ego@result) > 0) { go_list[[i]] <- ego }
    if (nrow(kk@result) > 0) { kegg_list[[i]] <- kk }
    }
    gsea_list[["go_list"]] <- go_list
    gsea_list[["kegg_list"]] <- kegg_list
    gsea_list
} 

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
