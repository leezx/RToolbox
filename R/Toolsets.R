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

# a revised version of plot_cell_trajectory in monocle package
# add use.cells, you can remove any cells in the plot, eg.sampling cells
plot_cell_trajectory2 <- function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, use.cells,
    show_backbone = TRUE, backbone_color = "black", markers = NULL, 
    use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
    show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
    cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
    theta = 0, ...) 
{
    library(tibble)
    source("https://raw.githubusercontent.com/cole-trapnell-lab/monocle-release/d8940700e70689ebc2dc6c7c4d5929a27f2de5e2/R/plotting.R")
    requireNamespace("igraph")
    gene_short_name <- NA
    sample_name <- NA
    sample_state <- pData(cds)$State
    data_dim_1 <- NA
    data_dim_2 <- NA
    lib_info_with_pseudo <- pData(cds)
    if (is.null(cds@dim_reduce_type)) {
        stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }
    if (cds@dim_reduce_type == "ICA") {
        reduced_dim_coords <- reducedDimS(cds)
    }
    else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
        reduced_dim_coords <- reducedDimK(cds)
    }
    else {
        stop("Error: unrecognized dimensionality reduction method.")
    }
    ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
        select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>% 
        mutate(sample_name = rownames(.), sample_state = rownames(.))
    dp_mst <- minSpanningTree(cds)
    if (is.null(dp_mst)) {
        stop("You must first call orderCells() before using this function")
    }
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from", 
        target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name", 
        source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
        by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name", 
        target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), 
        by = "target")
    data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
        select_(data_dim_1 = x, data_dim_2 = y) %>% rownames_to_column("sample_name") %>% 
        mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
        rownames_to_column("sample_name"), by = "sample_name")
    return_rotation_mat <- function(theta) {
        theta <- theta/180 * pi
        matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
            nrow = 2)
    }
    rot_mat <- return_rotation_mat(theta)
    cn1 <- c("data_dim_1", "data_dim_2")
    cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
    cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
    data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
    edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
    edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% 
            markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), 
                ])))
            colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
            markers_exprs <- merge(markers_exprs, markers_fData, 
                by.x = "feature_id", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
            by.y = "cell_id")
        if (use_color_gradient) {
            if (markers_linear) {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2)) + geom_point(aes(color = value), 
                  size = I(cell_size), na.rm = TRUE) + scale_color_viridis(name = paste0("value"), 
                  ...) + facet_wrap(~feature_label)
            }
            else {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2)) + geom_point(aes(color = log10(value + 
                  0.1)), size = I(cell_size), na.rm = TRUE) + 
                  scale_color_viridis(name = paste0("log10(value + 0.1)"), 
                    ...) + facet_wrap(~feature_label)
            }
        }
        else {
            if (markers_linear) {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2, size = (value * 0.1))) + facet_wrap(~feature_label)
            }
            else {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2, size = log10(value + 0.1))) + 
                  facet_wrap(~feature_label)
            }
        }
    }
    else {
        g <- ggplot(data = subset(data_df, sample_name %in% use.cells), aes(x = data_dim_1, y = data_dim_2))
    }
    if (show_tree) {
        g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
            y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
            yend = "target_prin_graph_dim_2"), size = cell_link_size, 
            linetype = "solid", na.rm = TRUE, data = edge_df)
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        if (use_color_gradient) {
        }
        else {
            g <- g + geom_point(aes_string(color = color_by), 
                na.rm = TRUE)
        }
    }
    else {
        if (use_color_gradient) {
        }
        else {
            g <- g + geom_point(aes_string(color = color_by), 
                size = I(cell_size), na.rm = TRUE)
        }
    }
    if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
        mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
        branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes, 
            sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
        g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
            branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
            y = "prin_graph_dim_2", label = "branch_point_idx"), 
            size = 4, color = "white", na.rm = TRUE, branch_point_df)
    }
    if (show_cell_names) {
        g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
    }
    if (show_state_number) {
        g <- g + geom_text(aes(label = sample_state), size = state_number_size)
    }
    g <- g + monocle_theme_opts() + xlab(paste("Component", x)) + 
        ylab(paste("Component", y)) + theme(legend.position = "top", 
        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = "white"))
    g
    # head(data_df)
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

# read the SNP annotation file from annovar software
readInfile.annovar <- function(path, sample) {
    IMR90 <- read.table(path, header = F, sep = "\t", stringsAsFactors = F)
    IMR90$id <- paste(IMR90$V4, IMR90$V5, IMR90$V6, IMR90$V8, sep = ":")
    IMR90 <- subset(IMR90, !V2 %in% c("synonymous SNV", "unknown"))
    IMR90$sample <- sample
    IMR90$type <- as.character(IMR90$V2)
    IMR90$ref <- IMR90$V7
    IMR90$alt <- IMR90$V8
    IMR90$allele <- IMR90$V9
    IMR90$gene <- unlist(lapply(IMR90$V3, function(x) {
        split.list <- strsplit(x, split = ":")
        split.list[[1]][1]
        #length(x)
    }))
    IMR90 <- IMR90[,c("id", "sample", "type", "gene", "ref", "alt", "allele")]
    IMR90
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
