# 常用空值兜底：x 为空或全 NA 时返回 y。
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# 确保输出目录存在，并返回目录路径。
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# 检查脚本运行所需 R 包是否已安装。
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Install required packages before running this script: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# 将细胞类型等标签转换为适合文件名使用的安全字符串。
sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

# 统一管理项目中数据、注释表、RDS、图和表的输出路径。
get_pipeline_paths <- function(project_root = here::here()) {
  root <- file.path(project_root, "snRNAseq_publication_pipeline")
  list(
    project_root = project_root,
    pipeline_root = root,
    data_dir = file.path(project_root, "dnbc4tools_results"),
    map_dir = file.path(root, "annotation_maps"),
    rds_dir = ensure_dir(file.path(root, "outputs", "rds")),
    fig_dir = ensure_dir(file.path(root, "outputs", "figures")),
    table_dir = ensure_dir(file.path(root, "outputs", "tables")),
    report_dir = ensure_dir(file.path(root, "outputs", "reports"))
  )
}

# 根据样本名生成样本信息表，包括分组、组织和 MouseID。
make_sample_manifest <- function(samples = c(
  "CTM1", "CTM2", "CTM3",
  "CTS1", "CTS2", "CTS3",
  "MIXM1", "MIXM2", "MIXM3",
  "MIXS1", "MIXS2", "MIXS3"
)) {
  stopifnot(length(samples) > 0)
  data.frame(
    SampleID = samples,
    Group = ifelse(grepl("^CT", samples), "Control", "Experimental"),
    Tissue = ifelse(grepl("M[0-9]+$", samples), "BoneMarrow", "Spleen"),
    MouseNumber = sub("^.*?([0-9]+)$", "\\1", samples),
    stringsAsFactors = FALSE
  ) |>
    transform(MouseID = paste(Group, MouseNumber, sep = "_"))
}

# 读取 10X 矩阵；如果包含多个层，优先选择 Gene Expression/RNA 层。
read_10x_counts <- function(data_dir) {
  x <- Seurat::Read10X(data.dir = data_dir)
  if (is.list(x)) {
    gene_layer <- grep("Gene Expression|RNA", names(x), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(gene_layer)) gene_layer <- names(x)[1]
    x <- x[[gene_layer]]
  }
  x
}

# 给矩阵细胞条码加样本前缀，避免合并样本时 barcode 重名。
prefix_cells <- function(mat, prefix) {
  colnames(mat) <- paste(prefix, colnames(mat), sep = "_")
  mat
}

# 多个矩阵只保留共同基因后按列合并。
cbind_common_genes <- function(mats) {
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  if (length(common_genes) == 0) stop("No shared genes across matrices.", call. = FALSE)
  mats <- lapply(mats, function(x) x[common_genes, , drop = FALSE])
  Reduce(Matrix::cbind2, mats)
}

# 把 manifest 中的样本元数据写入 Seurat 对象。
add_basic_metadata <- function(obj, sample_row) {
  obj$SampleID <- sample_row$SampleID
  obj$MouseID <- sample_row$MouseID
  obj$Group <- sample_row$Group
  obj$Tissue <- sample_row$Tissue
  obj$MouseNumber <- sample_row$MouseNumber
  obj
}

# 计算线粒体和核糖体基因比例，用于 snRNA-seq QC。
add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^MT-")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^Rpl|^Rps|^RPL|^RPS")
  obj
}

# 对单个 Seurat 对象执行 snRNA-seq QC，并返回过滤后对象和统计摘要。
sn_qc_filter <- function(obj, thresholds) {
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    stop("sn_qc_filter expected a Seurat object, but received NULL or a non-Seurat object.", call. = FALSE)
  }

  sample_id <- unique(obj$SampleID)
  if (length(sample_id) != 1 || is.na(sample_id)) sample_id <- obj@project.name %||% "unknown"

  required_thresholds <- c("min_features", "max_features", "min_counts", "max_counts", "max_percent_mt")
  missing_thresholds <- setdiff(required_thresholds, names(thresholds))
  if (is.null(thresholds) || length(missing_thresholds) > 0) {
    stop(
      "QC thresholds are missing for sample ",
      sample_id,
      if (length(missing_thresholds) > 0) paste0(": ", paste(missing_thresholds, collapse = ", ")) else ".",
      call. = FALSE
    )
  }

  meta_before <- obj@meta.data
# 根据 feature、UMI 和线粒体比例阈值筛选细胞。
  keep <- with(
    meta_before,
    nFeature_RNA >= thresholds$min_features &
      nFeature_RNA <= thresholds$max_features &
      nCount_RNA >= thresholds$min_counts &
      nCount_RNA <= thresholds$max_counts &
      percent.mt <= thresholds$max_percent_mt
  )
  keep[is.na(keep)] <- FALSE

  if (!any(keep)) {
    stop(
      "No cells passed snRNA-seq QC for sample ",
      sample_id,
      ". Thresholds: ",
      paste(paste(required_thresholds, unlist(thresholds[required_thresholds]), sep = "="), collapse = ", "),
      call. = FALSE
    )
  }

# 在原对象中记录 QC 结果，再按通过细胞子集化。
  obj$pass_sn_qc <- keep
  filtered <- subset(obj, cells = rownames(meta_before)[keep])
  if (is.null(filtered) || !inherits(filtered, "Seurat")) {
    stop("QC filtering returned NULL for sample ", sample_id, ".", call. = FALSE)
  }

  summary <- rbind(
    data.frame(Stage = "Before_QC", Cells = nrow(meta_before)),
    data.frame(Stage = "After_QC", Cells = ncol(filtered))
  )
  summary$Removed <- c(NA_integer_, nrow(meta_before) - ncol(filtered))
  summary$RemovedPercent <- c(NA_real_, 100 * (nrow(meta_before) - ncol(filtered)) / nrow(meta_before))

  list(obj = filtered, summary = summary)
}

# 统一的出版风格 ggplot 主题。
theme_publication <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

# 同时保存 PDF 和 300 dpi PNG 图。
save_pub_plot <- function(plot, filename, fig_dir, width = 7, height = 5, dpi = 300) {
  ensure_dir(fig_dir)
  base <- file.path(fig_dir, filename)
  ggplot2::ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height, useDingbats = FALSE)
  ggplot2::ggsave(paste0(base, ".png"), plot = plot, width = width, height = height, dpi = dpi)
  invisible(base)
}

# 绘制每个样本的 QC 指标小提琴图。
plot_qc_violin <- function(meta, title = "snRNA-seq QC metrics") {
  meta |>
    dplyr::select(SampleID, Tissue, nFeature_RNA, nCount_RNA, percent.mt, percent.ribo) |>
    tidyr::pivot_longer(
      cols = c(nFeature_RNA, nCount_RNA, percent.mt, percent.ribo),
      names_to = "Metric",
      values_to = "Value"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = SampleID, y = Value, fill = Tissue)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# 运行 SCTransform、PCA、可选 Harmony、邻居图、聚类和 UMAP。
run_sctransform_embedding <- function(
    obj,
    dims = 1:30,
    resolution = 0.6,
    vars_to_regress = "percent.mt",
    harmony_batch_var = NULL,
    n_neighbors = 30,
    min_dist = 0.3,
    sct_method = "glmGamPoi",
    sct_vst_flavor = "v2") {
# 禁止把样本或 mouse 这类生物重复变量当作 Harmony batch。
  forbidden <- c("sample_id", "sampleid", "sample", "mouse_id", "mouseid", "mouse", "SampleID", "MouseID")
  if (!is.null(harmony_batch_var) && harmony_batch_var %in% forbidden) {
    stop("Do not use SampleID, sample_id, MouseID, or mouse_id as the Harmony batch variable.", call. = FALSE)
  }

# 默认使用 glmGamPoi 加速 SCTransform，并回归线粒体比例。
  obj <- Seurat::SCTransform(
    obj,
    method = sct_method,
    vst.flavor = sct_vst_flavor,
    vars.to.regress = vars_to_regress,
    verbose = FALSE
  )
  obj <- Seurat::RunPCA(obj, npcs = max(dims), verbose = FALSE)

# 如果提供了真正的技术 batch 变量，则基于 PCA 运行 Harmony。
  reduction_for_graph <- "pca"
  if (!is.null(harmony_batch_var) && harmony_batch_var %in% colnames(obj@meta.data)) {
    obj <- harmony::RunHarmony(
      obj,
      group.by.vars = harmony_batch_var,
      assay.use = "SCT",
      reduction.use = "pca",
      verbose = FALSE
    )
    reduction_for_graph <- "harmony"
  }

# 使用 PCA 或 Harmony 低维空间构建图、聚类并计算 UMAP。
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_for_graph, dims = dims)
  obj <- Seurat::FindClusters(obj, resolution = resolution)
  obj <- Seurat::RunUMAP(
    obj,
    reduction = reduction_for_graph,
    dims = dims,
    n.neighbors = n_neighbors,
    min.dist = min_dist,
    spread = 1
  )
  obj@misc$analysis_reduction <- reduction_for_graph
  obj
}

# 使用 SingleR 进行单细胞注释，并按 Seurat cluster 做多数投票生成主标签。
run_singler_main <- function(obj, ref, labels = ref$label.main, assay = "SCT", workers = 6) {
  sce <- as.SingleCellExperiment(obj, assay = assay)
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = labels,
    de.method = "wilcox",
    BPPARAM = BiocParallel::MulticoreParam(workers = workers)
  )
  obj$CellType_SingleR_Raw <- pred$labels

# 每个 cluster 取 SingleR 原始标签中的多数标签作为 CellType_Main。
  majority_vote <- tapply(obj$CellType_SingleR_Raw, Seurat::Idents(obj), function(x) {
    tbl <- table(x)
    names(tbl)[which.max(tbl)]
  })
  obj$CellType_Main <- unname(majority_vote[as.character(Seurat::Idents(obj))])
  obj
}

# 读取手动注释表，并返回 cluster 到 CellType 的命名向量。
read_annotation_map <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Cluster", "CellType")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) stop("Annotation map is missing columns: ", paste(missing, collapse = ", "))
  x$Cluster <- as.character(x$Cluster)
  stats::setNames(x$CellType, as.character(x$Cluster))
}

# 根据当前对象实际 cluster 生成手动注释模板，并合并旧 map 中已有标签。
write_annotation_template <- function(obj, map_path, template_path, cluster_col = "seurat_clusters") {
  old_map <- if (file.exists(map_path)) {
    read.csv(map_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(Cluster = character(), CellType = character(), MarkerEvidence = character())
  }
  old_map$Cluster <- as.character(old_map$Cluster)

# observed 包含当前聚类实际存在的 cluster 和每个 cluster 的细胞数。
  observed <- data.frame(
    Cluster = sort(unique(as.character(obj@meta.data[[cluster_col]]))),
    stringsAsFactors = FALSE
  )
  observed <- dplyr::left_join(observed, old_map, by = "Cluster")
  observed$CellType[is.na(observed$CellType)] <- ""
  observed$MarkerEvidence[is.na(observed$MarkerEvidence)] <- ""
  observed$Cells <- as.integer(table(as.character(obj@meta.data[[cluster_col]]))[observed$Cluster])
  write.csv(observed, template_path, row.names = FALSE)
  invisible(observed)
}

# 检查手动注释表是否覆盖当前所有 cluster，且 CellType 不为空。
validate_annotation_map <- function(obj, map_path, cluster_col = "seurat_clusters") {
  map <- read.csv(map_path, stringsAsFactors = FALSE, check.names = FALSE)
  map$Cluster <- as.character(map$Cluster)
  observed_clusters <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
  mapped_clusters <- sort(unique(as.character(map$Cluster)))
# missing 表示当前聚类有但 map 没有；extra 表示 map 有但当前聚类没有。
  missing_clusters <- setdiff(observed_clusters, mapped_clusters)
  extra_clusters <- setdiff(mapped_clusters, observed_clusters)
  blank_labels <- map$Cluster[is.na(map$CellType) | trimws(map$CellType) == ""]

  list(
    ok = length(missing_clusters) == 0 && length(blank_labels) == 0,
    missing_clusters = missing_clusters,
    extra_clusters = extra_clusters,
    blank_labels = blank_labels,
    map = stats::setNames(map$CellType, as.character(map$Cluster))
  )
}

# 按 cluster 到 CellType 的映射写入最终 CellType_Fine 标签。
apply_manual_annotation <- function(obj, map, cluster_col = "seurat_clusters", unknown_label = "Unknown") {
  clusters <- as.character(obj@meta.data[[cluster_col]])
  labels <- unname(map[clusters])
  labels[is.na(labels)] <- paste(unknown_label, clusters[is.na(labels)])
  obj$CellType_Fine <- labels
  obj
}

# 统计 marker 在每个分组中的平均表达和表达细胞比例。
marker_validation_table <- function(obj, markers, group_col = "CellType_Fine", assay = "RNA") {
  obj <- safe_join_layers(obj, assay = assay)
  data_mat <- tryCatch(get_assay_matrix(obj, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(data_mat) || ncol(data_mat) == 0) {
    obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
    data_mat <- get_assay_matrix(obj, assay = assay, layer = "data")
  }
  markers <- intersect(unique(markers), rownames(data_mat))
  groups <- obj@meta.data[[group_col]]
# 按细胞类型或 cluster 分组后逐组计算 marker 统计量。
  split_cells <- split(seq_along(groups), groups)

  out <- lapply(names(split_cells), function(g) {
    idx <- split_cells[[g]]
    mat <- data_mat[markers, idx, drop = FALSE]
    data.frame(
      Group = g,
      Gene = markers,
      AvgExpression = Matrix::rowMeans(mat),
      PctExpressing = Matrix::rowMeans(mat > 0) * 100,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

# Seurat v5 多 layer assay 在取矩阵前先合并 layer。
safe_join_layers <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5") && length(SeuratObject::Layers(obj[[assay]])) > 1) {
    obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
  }
  obj
}

# 兼容 Seurat v5 layer 和旧版 slot 的矩阵读取函数。
get_assay_matrix <- function(obj, assay = "RNA", layer = "counts") {
  if (inherits(obj[[assay]], "Assay5")) {
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer)
  } else {
    SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
  }
}

# 绘制 marker dotplot，并自动忽略对象中不存在的基因。
make_marker_dotplot <- function(obj, features, group.by = "CellType_Fine", title = NULL) {
  features <- intersect(features, rownames(obj))
  Seurat::DotPlot(obj, features = features, group.by = group.by, cols = c("grey90", "#B2182B")) +
    Seurat::RotatedAxis() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# 统一 Control/Experimental 分组水平，兼容 CT/MIX 等样本名前缀。
normalize_group_levels <- function(x) {
  x <- as.character(x)
  x[grepl("^control|^ct", x, ignore.case = TRUE)] <- "Control"
  x[grepl("experimental|^mix", x, ignore.case = TRUE)] <- "Experimental"
  factor(x, levels = c("Control", "Experimental"))
}

# 计算每只 mouse 中各细胞类型的细胞数和比例。
calc_mouse_proportions <- function(obj, celltype_col = "CellType_Fine") {
  meta <- obj@meta.data |>
    dplyr::mutate(Group = normalize_group_levels(Group))

  totals <- meta |>
    dplyr::count(Tissue, Group, MouseID, name = "TotalCells")

  counts <- meta |>
    dplyr::count(Tissue, Group, MouseID, CellType = .data[[celltype_col]], name = "CellCount")

  all_celltypes <- sort(unique(meta[[celltype_col]]))

# complete 补齐缺失组合，使某 mouse 中不存在的细胞类型比例记为 0。
  counts |>
    tidyr::complete(
      Tissue,
      tidyr::nesting(Group, MouseID),
      CellType = all_celltypes,
      fill = list(CellCount = 0)
    ) |>
    dplyr::left_join(totals, by = c("Tissue", "Group", "MouseID")) |>
    dplyr::mutate(Proportion = CellCount / TotalCells, Percentage = 100 * Proportion)
}

# 对每个组织和细胞类型做 mouse-level 组间比例检验。
test_mouse_proportions <- function(prop_df) {
  prop_df |>
    dplyr::group_by(Tissue, CellType) |>
    dplyr::summarise(
      N_Control = sum(Group == "Control"),
      N_Experimental = sum(Group == "Experimental"),
      Mean_Control = mean(Percentage[Group == "Control"], na.rm = TRUE),
      Mean_Experimental = mean(Percentage[Group == "Experimental"], na.rm = TRUE),
      Delta_Percentage = Mean_Experimental - Mean_Control,
# 每组至少 2 只 mouse 时才进行 Wilcoxon 检验。
      P_Value = ifelse(
        N_Control >= 2 && N_Experimental >= 2,
        stats::wilcox.test(Percentage ~ Group, exact = FALSE)$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(FDR = stats::p.adjust(P_Value, method = "BH"))
}

# 对指定细胞类型按 MouseID 聚合 raw counts，生成 pseudobulk 矩阵。
make_pseudobulk <- function(
    obj,
    celltype,
    celltype_col = "CellType_Fine",
    sample_col = "MouseID",
    group_col = "Group",
    assay = "RNA",
    min_cells_per_mouse = 20) {
  obj <- safe_join_layers(obj, assay = assay)
  counts <- get_assay_matrix(obj, assay = assay, layer = "counts")
  meta <- obj@meta.data
# 只保留目标细胞类型的细胞。
  cells <- rownames(meta)[meta[[celltype_col]] == celltype]
  if (length(cells) == 0) return(NULL)

  meta_sub <- meta[cells, , drop = FALSE] |>
    dplyr::mutate(Group = normalize_group_levels(.data[[group_col]]))

  cell_counts <- meta_sub |>
    dplyr::count(.data[[sample_col]], name = "N_Cells") |>
    dplyr::filter(N_Cells >= min_cells_per_mouse)

# 过滤细胞数过少的 mouse，并要求至少 4 个 pseudobulk 样本且包含两组。
  keep_samples <- cell_counts[[sample_col]]
  meta_sub <- meta_sub[meta_sub[[sample_col]] %in% keep_samples, , drop = FALSE]
  if (length(unique(meta_sub[[sample_col]])) < 4) return(NULL)
  if (length(unique(meta_sub$Group)) < 2) return(NULL)

# 对每只 mouse 内所有目标细胞的 raw counts 求和。
  split_cells <- split(rownames(meta_sub), meta_sub[[sample_col]])
  pb_counts <- do.call(cbind, lapply(split_cells, function(cell_ids) {
    Matrix::rowSums(counts[, cell_ids, drop = FALSE])
  }))
  colnames(pb_counts) <- names(split_cells)
  pb_counts <- round(as.matrix(pb_counts))

# 构建与 pseudobulk 矩阵列顺序一致的样本信息表。
  sample_info <- meta_sub |>
    dplyr::distinct(MouseID = .data[[sample_col]], Group, Tissue) |>
    dplyr::left_join(cell_counts, by = stats::setNames(sample_col, "MouseID")) |>
    dplyr::arrange(match(MouseID, colnames(pb_counts)))
  rownames(sample_info) <- sample_info$MouseID
  sample_info <- sample_info[colnames(pb_counts), , drop = FALSE]

  list(counts = pb_counts, sample_info = sample_info)
}

# 使用 edgeR quasi-likelihood 框架进行 pseudobulk 差异分析。
run_edgeR_pseudobulk <- function(pb, ref_group = "Control", contrast_group = "Experimental") {
  group <- stats::relevel(factor(pb$sample_info$Group), ref = ref_group)
  dge <- edgeR::DGEList(counts = pb$counts, group = group)
  keep <- edgeR::filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~ group)
  colnames(design) <- make.names(colnames(design))
# 估计离散度后拟合 GLM，并检验 Experimental 相对 Control 的系数。
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  coef_name <- grep(paste0("group", contrast_group), colnames(design), value = TRUE)[1]
  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  res <- edgeR::topTags(qlf, n = Inf)$table
  res$Gene <- rownames(res)
  res |>
    dplyr::select(Gene, logFC, logCPM, F, PValue, FDR) |>
    dplyr::arrange(FDR, dplyr::desc(abs(logFC)))
}

# 可选 DESeq2 pseudobulk 差异分析，作为 edgeR 结果的补充。
run_deseq2_pseudobulk <- function(pb, ref_group = "Control", contrast_group = "Experimental") {
  coldata <- pb$sample_info
  coldata$Group <- stats::relevel(factor(coldata$Group), ref = ref_group)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(pb$counts),
    colData = coldata,
    design = ~Group
  )
# 过滤低表达基因，减少不稳定检验。
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("Group", contrast_group, ref_group))
  res <- as.data.frame(res)
  res$Gene <- rownames(res)
  res |>
    dplyr::rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj) |>
    dplyr::select(Gene, baseMean, logFC, lfcSE, stat, PValue, FDR) |>
    dplyr::arrange(FDR, dplyr::desc(abs(logFC)))
}

# 绘制 DEG 火山图，按 FDR 和 logFC 标记上调、下调和不显著基因。
plot_deg_volcano <- function(res, title, fdr_cutoff = 0.05, logfc_cutoff = 0.5) {
  res |>
    dplyr::mutate(
      Direction = dplyr::case_when(
        FDR < fdr_cutoff & logFC >= logfc_cutoff ~ "Up",
        FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      NegLog10FDR = -log10(pmax(FDR, .Machine$double.xmin))
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = logFC, y = NegLog10FDR, color = Direction)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(Down = "#2166AC", `Not significant` = "grey70", Up = "#B2182B")) +
    ggplot2::geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = 2, linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, linewidth = 0.25) +
    ggplot2::labs(title = title, x = "log2 fold-change", y = "-log10 FDR") +
    theme_publication()
}
