# run_gate_evidence.R — Evidence generation for the manual-annotation gates (persistent script, under version control).
#
# Usage (run from /home/abagofmilk/scDATA):
#   Rscript snRNAseq_publication_pipeline/run_gate_evidence.R                # process all pending gates (those with a gate_status_*.csv)
#   Rscript snRNAseq_publication_pipeline/run_gate_evidence.R 01_bcell       # process only the specified gate(s) (multiple allowed)
#   Rscript snRNAseq_publication_pipeline/run_gate_evidence.R --save-rds     # additionally save compartment/main RDS for follow-up evidence
#
# Shares the same set of functions from R/gate_cluster_helpers.R, the same
# set.seed(seed) anchor, and the same parallel settings as the 01/02 Rmds,
# ensuring cluster numbering in the evidence matches the official render
# (replaces the handwritten mirror scripts in agent_scratch from before
# 2026-08-03, eliminating a structural source of numbering drift).
# To keep the RNG sequence identical: within a tissue, compartments that
# precede the requested gate in registry order are also recomputed in order
# (without exporting evidence).
# Outputs: outputs/gate_evidence/<gate_id>/ (07 only aggregates outputs/tables
# and outputs/figures; this directory is not absorbed into the publication flow).

args <- commandArgs(trailingOnly = TRUE)
save_rds <- "--save-rds" %in% args
requested <- setdiff(args, "--save-rds")

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else getwd()
if (!file.exists(file.path(script_dir, "R", "pipeline_helpers.R"))) {
  script_dir <- file.path(getwd(), "snRNAseq_publication_pipeline")
}
pipeline_root <- normalizePath(script_dir)
if (!file.exists(file.path(pipeline_root, "R", "pipeline_helpers.R"))) {
  stop("Cannot locate the pipeline root containing R/pipeline_helpers.R", call. = FALSE)
}

source(file.path(pipeline_root, "R", "pipeline_helpers.R"))
source(file.path(pipeline_root, "R", "gate_cluster_helpers.R"))
source(file.path(pipeline_root, "config", "pipeline_params.R"))

# Same package set and load order as the 01/02 Rmds.
check_packages(c(
  "Seurat", "SeuratObject", "SingleR", "celldex", "scuttle",
  "BiocParallel", "harmony", "dplyr", "ggplot2",
  "patchwork", "future", "glmGamPoi"
))

library(Seurat)
library(SingleR)
library(celldex)
library(scuttle)
library(BiocParallel)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

# Same anchors as the Rmds: thread caps → set.seed(seed) (per tissue) → plan/options.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

plan("multicore", workers = n_workers)
options(future.fork.enable = TRUE)
options(future.globals.maxSize = future_globals_maxsize)

paths <- get_pipeline_paths(project_root = dirname(pipeline_root))
defs <- gate_definitions()

# With no arguments, process all pending gates (those with a status file).
if (length(requested) == 0) {
  status_files <- list.files(paths$gate_dir, pattern = "^gate_status_.*\\.csv$", full.names = FALSE)
  requested <- intersect(names(defs), sub("^gate_status_(.*)\\.csv$", "\\1", status_files))
  if (length(requested) == 0) {
    message("No pending gates (no gate_status_*.csv under ", paths$gate_dir, "). Nothing to do.")
    if (!interactive()) quit(save = "no", status = 0)
  }
}
unknown_gates <- setdiff(requested, names(defs))
if (length(unknown_gates) > 0) {
  stop(
    "Unknown gate id(s): ", paste(unknown_gates, collapse = ", "),
    ". Known gates: ", paste(names(defs), collapse = ", "),
    call. = FALSE
  )
}

stamp <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

tissues <- unique(vapply(defs[requested], `[[`, "", "tissue"))
for (ti in tissues) {
  # Gate order for this tissue in the registry; run up to the last requested gate (keeps the RNG sequence identical).
  tissue_gates <- names(defs)[vapply(defs, `[[`, "", "tissue") == ti]
  requested_here <- intersect(tissue_gates, requested)
  run_seq <- tissue_gates[seq_len(max(match(requested_here, tissue_gates)))]

  stamp("tissue ", ti, ": loading ", defs[[tissue_gates[1]]]$qc_rds)
  set.seed(seed)
  obj <- readRDS(file.path(paths$rds_dir, defs[[tissue_gates[1]]]$qc_rds))
  obj$Group <- normalize_group_levels(obj$Group)
  stamp(ti, " cells: ", ncol(obj))

  stamp(ti, ": main clustering + SingleR (shared code path)")
  ref_immgen <- celldex::ImmGenData()
  obj <- run_main_clustering_annotation(
    obj,
    ref = ref_immgen,
    dims = dims_main,
    resolution = cluster_res_main,
    harmony_batch_var = harmony_batch_var,
    workers = n_workers
  )

  tissue_out <- ensure_dir(file.path(paths$gate_dir, ti))
  main_comp <- as.data.frame.matrix(table(obj$seurat_clusters, obj$CellType_Main))
  write.csv(
    data.frame(Cluster = rownames(main_comp), main_comp, check.names = FALSE),
    file.path(tissue_out, paste0(ti, "_main_clusters_x_CellTypeMain.csv")),
    row.names = FALSE
  )
  if (save_rds) saveRDS(obj, file.path(tissue_out, paste0(ti, "_main_clustered.rds")))

  for (gid in run_seq) {
    def <- defs[[gid]]
    stamp(gid, ": prepare compartment (", def$lineage, ")")
    compartment <- prepare_lineage_compartment(
      obj,
      label_pattern = def$label_pattern,
      dims = dims_subset,
      resolution = get(def$resolution),
      harmony_batch_var = harmony_batch_var,
      min_cells = min_cells_for_refinement,
      lineage_label = def$lineage
    )
    if (gid %in% requested_here) {
      out_dir <- ensure_dir(file.path(paths$gate_dir, gid))
      stamp(gid, ": dumping evidence -> ", out_dir)
      compartment <- dump_gate_evidence(
        compartment, gid,
        panel = get(def$panel),
        out_dir = out_dir,
        save_rds = save_rds
      )
      # Refresh the annotation template under the same rules as the gate trigger at render time, for side-by-side comparison with the evidence in the same directory.
      write_annotation_template(
        compartment,
        map_path = file.path(paths$map_dir, def$map_file),
        template_path = file.path(out_dir, paste0(gid, "_annotation_template.csv"))
      )
      stamp(gid, ": done (", length(unique(compartment$seurat_clusters)), " clusters, ",
            ncol(compartment), " cells)")
    }
  }
}

future::plan(future::sequential)
stamp("DONE")
