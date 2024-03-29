
# fit a principal_curve
# path: project/Data_center/analysis/ApcKO_multiomics/ApcKO-seurat.ipynb
pricu1 <- princurve::principal_curve(as.matrix(all.umap[,c("UMAP_1","UMAP_2")]),
                                     smoother='lowess', trace=F, plot_iterations=F, stretch=100, maxit=10, thresh=-1000)
pc.line1 <- as.data.frame(pricu1$s[order(pricu1$lambda), ])
all.umap$pseudotime <- pricu1$lambda/max(pricu1$lambda)
centers <- all.umap %>% dplyr::group_by(cluster) %>% summarize(UMAP_1 = median(UMAP_1),
                                                               UMAP_2 = median(UMAP_2))
options(repr.plot.width=5, repr.plot.height=5)
p <- ggplot(all.umap, aes(x=UMAP_1, y=UMAP_2, color=pseudotime)) +
  # facet_grid(cols = vars(variable)) +
  # facet_wrap( ~ variable, ncol=2) + # error in border
  geom_point(size=0.3, alpha=0.8) +
  geom_density_2d(color='black', size=0.05, alpha=0.15) +
  geom_text(data = centers, mapping = aes(label = cluster), size = 4.5, color="black") +
  geom_line(data=pc.line1, color='red', size=0.5) +
  labs(x = "UMAP_1",y = "UMAP_2", title = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) +
  # theme(legend.position = "none") +
  theme(strip.background = element_rect(fill = "gray97", color = NA)) + # strip background color
  theme(strip.placement = "outside", strip.text.x = element_text(face="plain", size = 14), #italic
        strip.text.y = element_text(face="plain", size = 11)) +
  theme(panel.spacing=unit(.3, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.2,
                                                                    colour = "black")) #+ #line size
#scale_color_manual(values=tmp.colors)
p

# a modified version of plot_cell_trajectory function in monocle package
# add use.cells, you can remove any cells in the plot, eg.sampling cells
# error history: a safer solution might be to use dplyr::slice in case you end up needing S4Vectors
# for something else. Perhaps another package relies on S4Vectors or even S4Vectors::slice explicitly.
plot_cell_trajectory2 <- function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE,
                                   use.cells, show_backbone = TRUE, backbone_color = "black",
                                   markers = NULL, use_color_gradient = FALSE, markers_linear = FALSE,
                                   show_cell_names = FALSE, show_state_number = FALSE, cell_size = 1.5,
                                   cell_link_size = 0.75, cell_name_size = 2, state_number_size = 2.9,
                                   show_branch_points = TRUE, theta = 0, ...)
{
  library(tibble)
  library(dplyr)
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
    g <- ggplot(data = subset(data_df, sample_name %in% use.cells),
                aes(x = data_dim_1, y = data_dim_2))
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
    branch_point_df <- ica_space_df %>% dplyr::slice(match(mst_branch_nodes,
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
}

# monocle2
# add API: use.cells = NULL
# function: subset the cells in the point
# example: scPipeline/human/Human-Organoid/Comparison_public_all/imputation_PHOX2B_pos.ipynb
plot_cell_trajectory.zx <- function (cds, x = 1, y = 2, color_by = "State", use.cells = NULL,
                                     show_backbone = TRUE, backbone_color = "black", markers = NULL, show_tree = TRUE,
                                     use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE,
                                     show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75,
                                     cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE,
                                     theta = 0, ...)
{
  # new packages
  require(monocle)
  require(dplyr)
  require(tibble)
  source("https://raw.githubusercontent.com/cole-trapnell-lab/monocle-release/master/R/plotting.R")
  #
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
  # subset cells
  if (!is.null(use.cells)) {
    data_df <- subset(data_df, sample_name %in% use.cells)
  }
  # plotting
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
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
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
