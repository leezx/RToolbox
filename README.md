# singleCellTools
my personal single cell RNA-seq analysis script R package

# script
prepare_smart_seq_data

prepare_10x_data

filtering_matrix

normalizing_matrix

feature_selection

Highly_variable_genes_selection

batch_effect_correction

dimension_reduction_by_PCA

dimension_reduction_by_diffusionMap



clustering_by_SC3

clustering_by_SIMLR

clustering_by_seurat

clustering_by_monocle

marker_identification_by_SC3

pseudotime_by_monocle1

pseudotime_by_monocle2

pseudotime_by_SLICER

module_identification_by_WGCNA

module_identification_by_MSIC



plotting_boxplot

plotting_violin_plot

plotting_barplot

plotting_volcano_plot

plotting_heatmap



DEG_by_scde

DEG_by_edgeR

scoring_model_by_LR



dataset_human_TFs

dataset_mouse_TFs

dataset_human_Y_chromsome_genes

dataset_mouse_Y_chromsome_genes

dataset_human_to_mouse_genes

dataset_cell_cycle_genes

dataset_cell_cycle_phase_genes

dataset_migration_genes

dataset_ENS_markers



show_notes

tools_study_conclusion


# Learning material

[R packages](http://r-pkgs.had.co.nz/)

Coursera - [Building R Packages](https://www.coursera.org/learn/r-packages/home/welcome)

From course:

- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
- [Building R Packages Pre-Flight Check List](https://github.com/rdpeng/daprocedures/blob/master/lists/Rpackage_preflight.md)
- [devtools-cheatsheet.pdf](https://www.rstudio.com/wp-content/uploads/2015/06/devtools-cheatsheet.pdf)
- [Common roxygen2 tags](https://bookdown.org/rdpeng/RProgDA/documentation.html#common-roxygen2-tags)
- [testthat: Get Started with Testing](https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf)
- 

## 以规范的R包标准来写这个包

- 将单细胞的分析按功能分类（smart-seq and 10x，clustering，pseudotime，绘图模块）；
- 标准化输入输出；
- 写好使用文档；

