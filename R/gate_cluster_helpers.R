# gate_cluster_helpers.R — 01/02 聚类注释与 gate 证据共用的代码路径。
#
# 目的：消灭"证据脚本镜像 Rmd"造成的 cluster 编号漂移。01/02 Rmd 与
# run_gate_evidence.R 必须调用这里定义的同一组函数（相同的调用序列 +
# 相同的 set.seed(seed) 锚点），聚类编号才能严格一致。
# 注释表校验本身（validate_annotation_map / check_annotation_gate）在
# pipeline_helpers.R；本文件只放聚类路径与证据导出。

# 手动注释 gate 登记表：每个 gate 的组织来源、lineage 筛选正则、map 文件。
# 正则与 lineage 名称是"镜像敏感"信息（改一处必须处处一致），唯一定义在这里；
# resolution / marker panel 仍由 config/pipeline_params.R 唯一定义，本表不重复。
# 调用前须已 source config/pipeline_params.R（表内不引用 config 变量，仅约定）。
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

# 主聚类 + SingleR 主注释：01/02 主流程与证据脚本共用的第一段。
# 依次执行 SCTransform 嵌入、SingleR（ImmGen）、并以 CellType_Main 初始化 CellType_Fine。
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

# 按 lineage 标签正则切出 compartment 并重聚类：01/02 与证据脚本共用。
# 筛选同时参考 SingleR 原始标签与多数投票主标签；细胞数低于 min_cells 时停止。
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

# 导出单个 gate 的注释证据：逐 cluster 细胞数、SingleR/主注释组成、
# config marker panel 表达、FindAllMarkers top 基因和 marker dotplot。
# 产物落在 outputs/gate_evidence/<gate_id>/（07 只汇总 outputs/tables 与
# outputs/figures，本目录不会被 publication 流程收编）。
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
