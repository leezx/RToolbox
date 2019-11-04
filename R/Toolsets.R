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

# get interaction of a gene list by STRING PPI database
get_interactions_by_PPI <- function(genes, species="human") {
    all.df <- data.frame()
    #
    if (species=="human") {
        species <- "Homo%20sapiens"
    } else if (species=="mouse") {
        species <- "Mus%20musculus"
    }
    for (i in genes) {
        tmpgene <- i
        url <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=",
                              tmpgene,"&species=",species, sep="")
        # 404 error
        options(warn=-1)
        options(show.error.messages= FALSE)
        options(stringsAsFactors = F)
        webDf <- try(read.table(url, header=T))
        # head(webDf)
        if (is.na(webDf)) next
        all.df <- rbind(all.df, webDf)
    }
    all.df
}

# get the candidate regulators in previous PPI dataframe
# stringId_A	STRING identifier (protein A)
# stringId_B	STRING identifier (protein B)
# preferredName_A	common protein name (protein A)
# preferredName_B	common protein name (protein B)
# ncbiTaxonId	NCBI taxon identifier
# score	combined score
# nscore	gene neighborhood score
# fscore	gene fusion score
# pscore	phylogenetic profile score
# ascore	coexpression score
# escore	experimental score
# dscore	database score
# tscore	textmining score
get_regulators_by_PPI <- function(all.df=interaction.df, min.num=5, filter=F) {
    # filter by score
    if (filter==T) {
        all.df <- subset(all.df, escore>0 | dscore>0 | tscore>0)
    }
    # 
    sort.regulators <- sort(table(interaction.df$preferredName_B), decreasing = T)
    sort.regulators <- names(sort.regulators)[sort.regulators >= min.num]
    #
    regulators.df <- data.frame()
    for (i in sort.regulators) {
        targets <- paste(subset(all.df, preferredName_B==i)$preferredName_A, collapse = ",")
        regulators.df <- rbind(regulators.df, data.frame(regulator=i, targets=targets))
    }
    regulators.df
}

# remove useless items from GO/KEGG annotation
# new.result <- filter.ora.items(result, filter = "mitochondrial|translational|riboso|Riboso")
filter.ora.items <- function(result, filter="") {
    for (i in names(result)) {
        for (j in names(result[[i]])) {
            tmpdf <- result[[i]][[j]]@result
            result[[i]][[j]]@result <- tmpdf[!grepl(filter, tmpdf$Description, ignore.case=T),]
        }
    }
    result
}

# scater inner function, for runPCA.detail using
get_mat_for_reddim <- function(object, exprs_values="logcounts", feature_set=NULL, ntop=500, scale=FALSE) {
# Picking the 'ntop' most highly variable features or just using a pre-specified set of features.
# Also removing zero-variance columns and scaling the variance of each column.
# Finally, transposing for downstream use (cells are now rows).

    exprs_mat <- assay(object, exprs_values, withDimnames=FALSE)
    rv <- rowVars(as.matrix(DelayedArray(as.matrix(exprs_mat))))

    if (is.null(feature_set)) {
        o <- order(rv, decreasing = TRUE)
        feature_set <- head(o, ntop)
    } else if (is.character(feature_set)) {
        feature_set <- .subset2index(feature_set, object, byrow=TRUE)
    }

    exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]
    rv <- rv[feature_set]

    exprs_to_plot <- t(exprs_to_plot)
    if (scale) {
        exprs_to_plot <- scater:::.scale_columns(exprs_to_plot, rv)
        rv <- rep(1, ncol(exprs_to_plot))
    }

    exprs_to_plot
}

# return the weight of each gene on each PCs
runPCA.detail <- function (object, ncomponents = 2, method = c("prcomp", "irlba"), 
    ntop = 500, exprs_values = "logcounts", feature_set = NULL, 
    scale_features = TRUE, use_coldata = FALSE, selected_variables = NULL, 
    detect_outliers = FALSE, rand_seed = NULL, ...) 
{
    if (use_coldata) {
        if (is.null(selected_variables)) {
            selected_variables <- list()
            it <- 1L
            for (field in c("pct_counts_in_top_100_features", 
                "total_features_by_counts", "pct_counts_feature_control", 
                "total_features_by_counts_feature_control", "log10_total_counts_endogenous", 
                "log10_total_counts_feature_control")) {
                out <- .qc_hunter(object, field, mode = "column", 
                  error = FALSE)
                if (!is.null(out)) {
                  selected_variables[[it]] <- out
                  it <- it + 1L
                }
            }
        }
        exprs_to_plot <- matrix(0, ncol(object), length(selected_variables))
        for (it in seq_along(selected_variables)) {
            exprs_to_plot[, it] <- .choose_vis_values(object, 
                selected_variables[[it]], mode = "column", search = "metadata")$val
        }
        if (scale_features) {
            exprs_to_plot <- scater:::.scale_columns(exprs_to_plot)
        }
    }
    else {
        exprs_to_plot <- get_mat_for_reddim(object, exprs_values = exprs_values, 
            ntop = ntop, feature_set = feature_set, scale = scale_features)
    }
    if (detect_outliers && use_coldata) {
        outliers <- mvoutlier::pcout(exprs_to_plot, makeplot = FALSE, 
            explvar = 0.5, crit.M1 = 0.9, crit.c1 = 5, crit.M2 = 0.9, 
            crit.c2 = 0.99, cs = 0.25, outbound = 0.05)
        outlier <- !as.logical(outliers$wfinal01)
        object$outlier <- outlier
    }
    method <- match.arg(method)
    if (method == "prcomp") {
        print(paste("using ", method, "...", sep=""))
        exprs_to_plot <- as.matrix(exprs_to_plot)
        ncomponents <- min(c(ncomponents, dim(exprs_to_plot)))
        pca <- prcomp(exprs_to_plot, rank. = ncomponents)
        percentVar <- pca$sdev^2
        percentVar <- percentVar/sum(percentVar)
    }
    else if (method == "irlba") {
        print(paste("using ", method, "...", sep=""))
        if (!is.null(rand_seed)) {
            .Deprecated(msg = "'rand.seed=' is deprecated.\nUse 'set.seed' externally instead.")
            set.seed(rand_seed)
        }
        ncomponents <- min(c(ncomponents, dim(exprs_to_plot) - 
            1L))
        pca <- irlba::prcomp_irlba(exprs_to_plot, n = ncomponents, 
            ...)
        percentVar <- pca$sdev^2/sum(colVars(DelayedArray(exprs_to_plot)))
    }
    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar
    if (use_coldata) {
        reducedDim(object, "PCA_coldata") <- pcs
    }
    else {
        reducedDim(object, "PCA") <- pcs
        # object$PCA_rotation <- pca$rotation
    }
    return(list(object=object, PCA_rotation=pca$rotation))
}

#' repeat the col or row for n times
#' @param x x
#' @param n n
#'
#' @return df
#' @export
#' @examples
#' rep.df <- rep.col(all_tsne[,cluster], length(groups)) 
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

#' repeat the col or row for n times
#' @param x x
#' @param n n
#'
#' @return df
#' @export
#' @examples
#' rep.df <- rep.col(all_tsne[,cluster], length(groups)) 
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' for facet plot, add background point
#' @param all_tsne all_tsne
#' @param cluster cluster
#' @param group group
#' @param cor1 x
#' @param cor2 y
#' @param sample.order sample.order
#' @param group.order group.order
#'
#' @return a merged dataframe
#' @export
#' @examples
#' merged_df <- add.background.point(all_tsne, sample.order = sample.order, group.order = group.order)
#'
add.background.point <- function(all_tsne, cluster="cluster", group="group", cor1="X", cor2="Y", 
                                 sample.order, group.order) {
    # bind_rows(replicate(3, all_tsne, simplify = F))
    groups <- unique(as.character(all_tsne[,group]))
    rep.df <- rep.col(all_tsne[,cluster], length(groups))
    colnames(rep.df) <- groups
    all_tsne_rep <- cbind(all_tsne, rep.df)
    # print(c(cor1, cor2, groups))
    # assign to others
    for (i in groups) {
        all_tsne_rep[all_tsne_rep$group!=i, i] <- "others"
        next
    }
    all_tsne_rep <- all_tsne_rep[,c(cor1, cor2, groups)]
    merged_df <- melt(all_tsne_rep, id.vars = c(cor1, cor2))
    merged_df$variable <- factor(merged_df$variable, levels = sample.order)
    merged_df$value <- factor(merged_df$value, levels = c(group.order, "others"))
    merged_df <- merged_df[order(merged_df$value, decreasing = T),]
    merged_df
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
#' @param FeatureList FeatureList
#' @param sortBy sortBy
#' @param colors colors
#'
#' @return anno list
#' @export
#' @examples
#' group.order <- c("Control", "S-HSCR", "L-HSCR")
#' FeatureList <- list("Sample"=sample.order[1:8], "Group"=group.order)
#' anno <- pre.pheatmap.anno.general(sce_HSCR_c1, FeatureList, sortBy = "Sample")
#'
pre.pheatmap.anno.general <- function(sce=sce, FeatureList=FeatureList, sortBy="impute_cluster", 
                                      colors=list("1"=brewer.pal(8,"Set2"), "2"=brewer.pal(12,"Set3"))) {
    # group.order=group.order, sample.order=sample.order, Cluster="impute_cluster", Sample="cellGroup"
    anno <- list()
    annotation_colors <- FeatureList
    annotation_col <- as.data.frame(colData(sce)[,c("cellName", names(FeatureList))])
    for (j in 1:length(names(FeatureList))) {
        # print(j)
        i <- names(FeatureList)[j]
        # print(i)
        tmpNames <- FeatureList[[i]]
        # print(tmpNames)
        annotation_col[,i] <- factor(annotation_col[,i], levels = tmpNames)
        tmpcolor <- colors[[j]][1:length(tmpNames)] # mistake [] and [[]]
        if (length(tmpNames)>8) {
            tmpcolor <- brewer.pal(12,"Set3")[1:length(tmpNames)]
        }
        names(tmpcolor) <- tmpNames
        # print(annotation_colors)
        annotation_colors[[i]]  <- tmpcolor
    }
    annotation_col <- annotation_col[order(annotation_col[,sortBy], decreasing = F),]
    annotation_col$cellName <- NULL
    #######
    anno[["col"]] <- annotation_col
    anno[["colors"]] <- annotation_colors
    return(anno)
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
