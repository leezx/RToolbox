# metacell
# remove some bugs
library(hdWGCNA)
MetacellsByGroups.zx <- function (seurat_obj, group.by = c("seurat_clusters"), ident.group = "seurat_clusters", 
    k = 25, reduction = "pca", assay = NULL, cells.use = NULL, 
    slot = "counts", mode = "average", min_cells = 100, max_shared = 15, 
    target_metacells = 1000, max_iter = 5000, verbose = FALSE, 
    wgcna_name = NULL) 
{
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    if (any(grepl("#", group.by))) {
        stop("Invalid character # found in group.by, please re-name the group.")
    }
    if (!(ident.group %in% group.by)) {
        stop("ident.group must be in group.by")
    }
    if (!(mode %in% c("sum", "average"))) {
        stop("Invalid choice for mode. Mode can be either sum or average.")
    }
    if (!(reduction %in% names(seurat_obj@reductions))) {
        stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", 
            paste(names(seurat_obj@reductions), collapse = ", ")))
    }
    if (is.null(assay)) {
        assay <- DefaultAssay(seurat_obj)
    }
    else if (!(assay %in% names(seurat_obj@assays))) {
        stop(paste0("Assay ", assay, " not found in seurat_obj. Select a valid assay: ", 
            paste0(names(seurat_obj@assays), collapse = ", ")))
    }
    if (!(slot %in% c("counts", "data", "scale.data"))) {
        stop("Invalid input for slot. Valid choices are counts, data, scale.data.")
    }
    else {
        slot_dim <- dim(GetAssayData(seurat_obj, assay = assay, 
            slot = slot))
        if (any(slot_dim) == 0) {
            stop(paste(c("Selected slot ", slot, " not found in this assay.")))
        }
    }
    if (!is.null(cells.use)) {
        seurat_full <- seurat_obj
        seurat_obj <- seurat_obj[, cells.use]
    }
    if (length(group.by) > 1) {
        seurat_meta <- seurat_obj@meta.data[, group.by]
        for (col in colnames(seurat_meta)) {
            seurat_meta[[col]] <- as.character(seurat_meta[[col]])
        }
        seurat_obj$metacell_grouping <- apply(seurat_meta, 1, 
            paste, collapse = "#")
    }
    else {
        seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
    }
    groupings <- unique(seurat_obj$metacell_grouping)
    groupings <- groupings[order(groupings)]
    group_counts <- table(seurat_obj$metacell_grouping) < min_cells
    if (any(group_counts)) {
        warning(paste0("Removing the following groups that did not meet min_cells: ", 
            paste(names(group_counts)[group_counts], collapse = ", ")))
    }
    groupings <- groupings[table(seurat_obj$metacell_grouping) >= 
        min_cells]
    if (length(groupings) == 0) {
        stop("No groups met the min_cells requirement.")
    }
    meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, 
        "#")))
    colnames(meta_df) <- group.by
    meta_list <- lapply(1:nrow(meta_df), function(i) {
        x <- list(as.character(meta_df[i, ]))[[1]]
        names(x) <- colnames(meta_df)
        x
    })
    seurat_list <- lapply(groupings, function(x) {
        seurat_obj[, seurat_obj$metacell_grouping == x]
    })
    names(seurat_list) <- groupings
    metacell_list <- mapply(ConstructMetacells, seurat_obj = seurat_list, 
        name = groupings, meta = meta_list, MoreArgs = list(k = k, 
            reduction = reduction, assay = assay, slot = slot, 
            return_metacell = TRUE, mode = mode, max_shared = max_shared, 
            max_iter = max_iter, target_metacells = target_metacells, 
            verbose = verbose, wgcna_name = wgcna_name))
    names(metacell_list) <- groupings
    remove <- which(sapply(metacell_list, is.null))
    if (length(remove) > 1) {
        metacell_list <- metacell_list[-remove]
    }
    run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, 
        function(x) {
            x@misc$run_stats
        })))
    rownames(run_stats) <- 1:nrow(run_stats)
    for (i in 1:length(group.by)) {
        run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), 
            "#"))[, i]
    }
    if (length(metacell_list) > 1) {
        metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
    }
    else {
        metacell_obj <- metacell_list[[1]]
    }
    Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]
    if (!is.null(cells.use)) {
        seurat_obj <- seurat_full
    }
    # seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, 
    #     wgcna_name)
    # seurat_obj <- SetWGCNAParams(seurat_obj, params = list(metacell_k = k, 
    #     metacell_reduction = reduction, metacell_slot = slot, 
    #     metacell_assay = assay, metacell_stats = run_stats), 
    #     wgcna_name)
    # seurat_obj
    output.list <- list()
    output.list[["seurat_obj"]] <- seurat_obj
    output.list[["metacell_obj"]] <- metacell_obj
    output.list
}


