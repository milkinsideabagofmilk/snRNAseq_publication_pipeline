if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown before running render_all.R", call. = FALSE)
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else "snRNAseq_publication_pipeline"
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

# 从集中参数文件读取 run_optional_cellchat 等设置。
source(file.path(script_dir, "config", "pipeline_params.R"))

rmd_files <- c(
  "snRNAseq_publication_pipeline/00_data_loading_soupx_qc.Rmd",
  "snRNAseq_publication_pipeline/01_spleen_clustering_annotation.Rmd",
  "snRNAseq_publication_pipeline/02_bm_clustering_annotation.Rmd",
  "snRNAseq_publication_pipeline/03_manual_marker_validation.Rmd",
  "snRNAseq_publication_pipeline/04_mouse_level_cell_proportion.Rmd",
  "snRNAseq_publication_pipeline/05_pseudobulk_deg_edgeR_DESeq2.Rmd",
  "snRNAseq_publication_pipeline/06_pseudobulk_go_kegg.Rmd",
  "snRNAseq_publication_pipeline/07_publication_figures_tables.Rmd"
)

if (isTRUE(run_optional_cellchat)) {
  rmd_files <- c(rmd_files, "snRNAseq_publication_pipeline/08_optional_cellchat.Rmd")
}

for (rmd in rmd_files) {
  message("Rendering: ", rmd)
  rmarkdown::render(rmd, output_dir = "snRNAseq_publication_pipeline/outputs/reports", clean = TRUE)
}
