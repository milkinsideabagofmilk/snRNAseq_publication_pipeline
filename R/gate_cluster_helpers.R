# gate_cluster_helpers.R — code path shared by the 01/02 clustering annotation and gate evidence.
#
# Purpose: eliminate cluster-number drift caused by "evidence scripts mirroring the Rmd".
# The 01/02 Rmd and run_gate_evidence.R must call the same set of functions defined here
# (same call sequence + same set.seed(seed) anchors) for cluster numbering to match exactly.
# Annotation-map validation itself (validate_annotation_map / check_annotation_gate) lives
# in pipeline_helpers.R; this file holds only the clustering path and evidence export.

# Registry of manual-annotation gates: tissue of origin, lineage filter regex, and map file
# for each gate. The regexes and lineage names are "mirror-sensitive" information (change one
# place and it must change everywhere); their single definition lives here. resolution and
# marker panel remain solely defined in config/pipeline_params.R and are not duplicated here.
# config/pipeline_params.R must already be sourced before calling (the table itself references
# no config variables; this is a convention only).
gate_definitions <- function() {
  list(
    "01_bcell" = list(
      tissue = "spleen",
      qc_rds = "00_spleen_qc.rds",
      label_pattern = "B cells|B cell|B lymph",
      lineage = "spleen B-lineage",
      map_file = "spleen_bcell_manual_map.csv",
      panel = "spleen_b_markers",
      resolution = "cluster_res_spleen_subset"
    ),
    "01_tcell" = list(
      tissue = "spleen",
      qc_rds = "00_spleen_qc.rds",
      label_pattern = "T cells|T cell|T lymph",
      lineage = "spleen T-lineage",
      map_file = "spleen_tcell_manual_map.csv",
      panel = "spleen_t_markers",
      resolution = "cluster_res_spleen_subset"
    ),
    "02_blineage" = list(
      tissue = "bonemarrow",
      qc_rds = "00_bonemarrow_qc.rds",
      label_pattern = "B cells|B cell|B lymph|plasma",
      lineage = "bone marrow B-lineage",
      map_file = "bm_blineage_manual_map.csv",
      panel = "bm_blineage_markers",
      resolution = "cluster_res_bm_subset"
    )
  )
}

# Main clustering + SingleR main annotation: the first stage shared by the 01/02 main
# pipeline and the evidence script.
# Runs SCTransform embedding, SingleR (ImmGen), then initializes CellType_Fine from CellType_Main.
run_main_clustering_annotation <- function(obj, ref, dims, resolution, harmony_batch_var, workers) {
  obj <- run_sctransform_embedding(
    obj,
    dims = dims,
    resolution = resolution,
    harmony_batch_var = harmony_batch_var
  )
  obj <- run_singler_main(obj, ref, workers = workers)
  obj$CellType_Fine <- obj$CellType_Main
  obj
}

# Subset the compartment by lineage-label regex and re-cluster it: shared by 01/02 and the
# evidence script.
# Filtering considers both the raw SingleR labels and the majority-vote main labels; stops
# when the cell count falls below min_cells.
prepare_lineage_compartment <- function(
    obj,
    label_pattern,
    dims,
    resolution,
    harmony_batch_var,
    min_cells,
    lineage_label) {
  labels <- grep(
    label_pattern,
    unique(c(obj$CellType_SingleR_Raw, obj$CellType_Main)),
    value = TRUE, ignore.case = TRUE
  )
  keep <- rownames(obj@meta.data)[
    obj$CellType_SingleR_Raw %in% labels | obj$CellType_Main %in% labels
  ]
  if (length(keep) < min_cells) {
    stop("Too few ", lineage_label, " cells for refinement. Check SingleR labels.", call. = FALSE)
  }

  compartment <- subset(obj, cells = keep)
  DefaultAssay(compartment) <- "RNA"
  compartment <- run_sctransform_embedding(
    compartment,
    dims = dims,
    resolution = resolution,
    harmony_batch_var = harmony_batch_var
  )
  compartment
}

# Export annotation evidence for a single gate: per-cluster cell counts, SingleR/main-annotation
# composition, config marker-panel expression, FindAllMarkers top genes, and a marker dotplot.
# Outputs land in outputs/gate_evidence/<gate_id>/ (step 07 only collects outputs/tables and
# outputs/figures; this directory is not folded into the publication pipeline).
dump_gate_evidence <- function(obj, gate_id, panel, out_dir, top_n = 8L, save_rds = FALSE) {
  ensure_dir(out_dir)
  cl <- as.character(obj$seurat_clusters)

  sizes <- as.data.frame.table(table(Cluster = cl, Group = obj$Group))
  write.csv(sizes, file.path(out_dir, paste0(gate_id, "_cluster_sizes_by_group.csv")), row.names = FALSE)

  sr <- as.data.frame.matrix(table(cl, obj$CellType_SingleR_Raw))
  write.csv(
    data.frame(Cluster = rownames(sr), sr, check.names = FALSE),
    file.path(out_dir, paste0(gate_id, "_clusters_x_SingleR_raw.csv")),
    row.names = FALSE
  )
  sm <- as.data.frame.matrix(table(cl, obj$CellType_Main))
  write.csv(
    data.frame(Cluster = rownames(sm), sm, check.names = FALSE),
    file.path(out_dir, paste0(gate_id, "_clusters_x_CellTypeMain.csv")),
    row.names = FALSE
  )

  mvt <- marker_validation_table(obj, panel, group_col = "seurat_clusters")
  write.csv(mvt, file.path(out_dir, paste0(gate_id, "_marker_panel_by_cluster.csv")), row.names = FALSE)

  obj <- Seurat::PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
  fim <- Seurat::FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, verbose = FALSE)
  write.csv(fim, file.path(out_dir, paste0(gate_id, "_findallmarkers_full.csv")), row.names = FALSE)
  top <- do.call(rbind, lapply(split(fim, fim$cluster), function(d) {
    d <- d[order(-d$avg_log2FC), , drop = FALSE]
    d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
  }))
  write.csv(top, file.path(out_dir, paste0(gate_id, "_findallmarkers_top", top_n, ".csv")), row.names = FALSE)

  p <- make_marker_dotplot(
    obj, panel,
    group.by = "seurat_clusters",
    title = paste0(gate_id, " marker panel by cluster")
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0(gate_id, "_dotplot_by_cluster.png")),
    p, width = 11, height = 6, dpi = 150
  )

  if (isTRUE(save_rds)) {
    saveRDS(obj, file.path(out_dir, paste0(gate_id, "_compartment.rds")))
  }
  invisible(obj)
}
