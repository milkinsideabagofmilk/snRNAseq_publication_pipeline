if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown before running render_all.R", call. = FALSE)
}

# render_all.R — full-pipeline driver (renders 00–07 in order).
#
# Usage:
#   Rscript snRNAseq_publication_pipeline/render_all.R           # resume (default)
#   Rscript snRNAseq_publication_pipeline/render_all.R --fresh   # archive old outputs/ then rerun everything
#
# 08 was split out (2026-08-06): CellChat has a huge, unstable peak memory
# footprint, so it no longer runs with the full pipeline. Launch it separately
# via snRNAseq_publication_pipeline/run_08_cellchat.sh (tmux + cgroup sandbox).
#
# Resume mode: for each step, checks whether the report HTML is newer than all
# of its inputs (the Rmd itself + config + the R/ function library + the map
# tables in annotation_maps); if so, the step is skipped. Once a step is
# rerun, all downstream steps are rerun in cascade (downstream steps consume
# upstream rds files). The report is written at the very last moment of
# rendering, so its existence means the step completed in full.
# gate_decisions.csv is a decision log, not an input to any step, and does not
# participate in this check.
#
# Exit-code contract: 0 = all done; 10 = a manual-annotation gate was hit (not
# an error — fill in the map and rerun this command to resume from the break
# point); 1 = a real error.

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else "snRNAseq_publication_pipeline"
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

fresh <- "--fresh" %in% commandArgs(trailingOnly = TRUE)

# With --fresh, archive the old outputs/ for version isolation (rename is
# instantaneous on the same filesystem and uses no extra disk). Archiving only
# happens on an explicit --fresh; neither resume nor rerunning a single Rmd
# triggers it. After the run, a human decides whether to keep or delete the
# old archive based on the quality of the new results (delete or move back
# manually).
archive_outputs <- function(outputs_dir, archive_root) {
  if (!dir.exists(outputs_dir)) return(invisible(NULL))
  n_files <- length(list.files(outputs_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE))
  if (n_files == 0) return(invisible(NULL))
  dir.create(archive_root, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(archive_root, format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (dir.exists(dest) || !file.rename(outputs_dir, dest)) {
    stop("Failed to archive previous outputs: ", outputs_dir, " -> ", dest, call. = FALSE)
  }
  dir.create(outputs_dir, showWarnings = FALSE)
  message("Archived previous outputs -> ", dest)
  invisible(dest)
}

if (fresh) {
  archive_outputs(
    outputs_dir = file.path(script_dir, "outputs"),
    archive_root = file.path(script_dir, "outputs_archive")
  )
}

rmd_files <- c(
  "snRNAseq_publication_pipeline/00_data_loading_soupx_qc.Rmd",
  "snRNAseq_publication_pipeline/01_spleen_clustering_annotation.Rmd",
  "snRNAseq_publication_pipeline/02_bm_clustering_annotation.Rmd",
  "snRNAseq_publication_pipeline/03_manual_marker_validation.Rmd",
  "snRNAseq_publication_pipeline/04_mouse_level_cell_proportion.Rmd",
  "snRNAseq_publication_pipeline/05_pseudobulk_deg_edgeR_DESeq2.Rmd",
  "snRNAseq_publication_pipeline/06_pseudobulk_go_kegg.Rmd",
  "snRNAseq_publication_pipeline/07_publication_figures_tables.Rmd"
  # 08_optional_cellchat.Rmd is not listed here: run it separately with
  # run_08_cellchat.sh (see the header comment above).
)

# Inputs for the resume check: shared inputs = the Rmd itself + config
# parameters + the R/ function library; annotation maps are inputs only for
# 01/02 (00 does not read maps; 03+ is covered indirectly by the 01/02 cascade
# through rds files). gate_decisions.csv is a decision log, not an input to
# any step, and does not participate in the check. Raw data
# (dnbc4tools_results) is assumed static; use --fresh when the data changes.
inputs_global <- c(
  file.path(script_dir, "config", "pipeline_params.R"),
  list.files(file.path(script_dir, "R"), pattern = "\\.R$", full.names = TRUE)
)
inputs_global <- inputs_global[file.exists(inputs_global)]
map_inputs <- list.files(file.path(script_dir, "annotation_maps"), pattern = "manual_map\\.csv$", full.names = TRUE)

step_needs_render <- function(rmd) {
  report <- file.path(script_dir, "outputs", "reports", sub("\\.Rmd$", ".html", basename(rmd)))
  if (!file.exists(report)) return(TRUE)
  inputs <- c(file.path(script_dir, basename(rmd)), inputs_global)
  if (grepl("^0[12]_", basename(rmd))) inputs <- c(inputs, map_inputs)
  any(file.mtime(inputs) > file.mtime(report))
}

rendered_upstream <- FALSE
for (rmd in rmd_files) {
  if (!fresh && !rendered_upstream && step_needs_render(rmd) == FALSE) {
    message("Skip (up to date): ", rmd)
    next
  }
  message("Rendering: ", rmd)
  err <- tryCatch(
    {
      rmarkdown::render(rmd, output_dir = "snRNAseq_publication_pipeline/outputs/reports", clean = TRUE)
      NULL
    },
    error = function(e) e
  )
  if (!is.null(err)) {
    msg <- conditionMessage(err)
    gate_hit <- regmatches(msg, regexpr("GATE\\[[A-Za-z0-9_]+\\]", msg))
    if (length(gate_hit) > 0) {
      gate_id <- sub("^GATE\\[|\\]$", "", gate_hit)
      status_file <- file.path(
        script_dir, "outputs", "gate_evidence", paste0("gate_status_", gate_id, ".csv")
      )
      message("\n==== Manual-annotation gate triggered: ", gate_id, " (not an error; waiting for an annotation decision) ====")
      message("Status file (clusters to label / cell counts / map and template paths): ", status_file)
      message("Standard procedure:")
      message("  1) Rscript snRNAseq_publication_pipeline/run_gate_evidence.R ", gate_id)
      message("  2) Fill in the map based on the evidence in outputs/gate_evidence/", gate_id, "/ (label ambiguous clusters as Unknown; reuse the existing label vocabulary)")
      message("  3) Record the decision in annotation_maps/gate_decisions.csv")
      message("  4) Rerun this command — resume skips completed steps and continues from this gate")
      if (!interactive()) quit(save = "no", status = 10) else stop(msg, call. = FALSE)
    }
    message("Render failed: ", rmd, "\n", msg)
    if (!interactive()) quit(save = "no", status = 1) else stop(msg, call. = FALSE)
  }
  # Cascade: this step was rerun, so downstream steps must not be skipped as
  # "report still up to date".
  rendered_upstream <- TRUE
}
message("render_all finished: all steps rendered.")
