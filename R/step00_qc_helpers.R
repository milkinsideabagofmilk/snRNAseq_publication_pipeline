# step00_qc_helpers.R — 00_data_loading_soupx_qc.Rmd 的专用函数库。
#
# 本文件只被 00 脚本 source()，收纳其重型逻辑：输入校验、单样本处理
# （scDblFinder + SoupX + QC）、细胞级审计、QC 表/图构建、合并校验。
# 全流程共用的函数在 R/pipeline_helpers.R；关键参数在 config/pipeline_params.R。

# ============================================================================
# 输入校验
# ============================================================================

# 硬校验样本 manifest：实验设计（2 组 × 2 组织 × 3 重复）、必填字段、
# doublet/rho 覆盖值范围、输入矩阵目录与 10X 三件套完整性。任何不满足即停止。
validate_manifest <- function(manifest) {
  required_columns <- c(
    "SampleID", "MouseID", "MouseNumber", "Group", "Tissue",
    "FilteredMatrixDir", "RawMatrixDir"
  )
  missing_columns <- setdiff(required_columns, colnames(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Manifest is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(manifest) != 12L || anyDuplicated(manifest$SampleID)) {
    stop("The expected design contains 12 unique SampleID values.", call. = FALSE)
  }
  if (anyNA(manifest[, required_columns, drop = FALSE]) ||
      any(trimws(manifest$SampleID) == "") ||
      any(trimws(manifest$MouseID) == "")) {
    stop("Manifest required fields cannot contain NA or empty values.", call. = FALSE)
  }
  if (!setequal(unique(manifest$Group), c("Control", "Experimental"))) {
    stop("Group must contain Control and Experimental.", call. = FALSE)
  }
  if (!setequal(unique(manifest$Tissue), c("Spleen", "BoneMarrow"))) {
    stop("Tissue must contain Spleen and BoneMarrow.", call. = FALSE)
  }

  design_counts <- dplyr::count(manifest, Group, Tissue, name = "n")
  if (nrow(design_counts) != 4L || any(design_counts$n != 3L)) {
    stop("Each Group x Tissue combination must contain exactly three samples.", call. = FALSE)
  }
  mouse_counts <- table(manifest$MouseID)
  if (length(mouse_counts) != 6L || any(mouse_counts != 2L)) {
    stop("Each of the six mice must contribute exactly two tissues.", call. = FALSE)
  }
  mouse_tissue_counts <- dplyr::count(manifest, MouseID, Tissue, name = "n")
  if (any(mouse_tissue_counts$n != 1L)) {
    stop("Each MouseID must contribute one spleen and one bone marrow sample.", call. = FALSE)
  }
  mouse_group_counts <- dplyr::count(manifest, MouseID, Group, name = "n")
  if (nrow(mouse_group_counts) != 6L) {
    stop("Each MouseID must belong to exactly one Group.", call. = FALSE)
  }

  manual_dbr <- suppressWarnings(as.numeric(manifest$ExpectedDoubletRate))
  bad_dbr <- !is.na(manual_dbr) &
    (!is.finite(manual_dbr) | manual_dbr <= 0 | manual_dbr >= 1)
  manual_rho <- suppressWarnings(as.numeric(manifest$SoupXManualRho))
  bad_rho <- !is.na(manual_rho) &
    (!is.finite(manual_rho) | manual_rho < 0 | manual_rho > 1)
  if (any(bad_dbr) || any(bad_rho)) {
    stop("Manifest doublet/rho overrides are outside valid ranges.", call. = FALSE)
  }

  missing_dirs <- c(
    manifest$FilteredMatrixDir[!dir.exists(manifest$FilteredMatrixDir)],
    manifest$RawMatrixDir[!dir.exists(manifest$RawMatrixDir)]
  )
  if (length(missing_dirs) > 0L) {
    stop(
      "Missing input matrix directories:\n",
      paste(unique(missing_dirs), collapse = "\n"),
      call. = FALSE
    )
  }

  has_10x_files <- function(path) {
    matrix_ok <- any(file.exists(file.path(
      path,
      c("matrix.mtx", "matrix.mtx.gz")
    )))
    barcode_ok <- any(file.exists(file.path(
      path,
      c("barcodes.tsv", "barcodes.tsv.gz")
    )))
    feature_ok <- any(file.exists(file.path(
      path,
      c(
        "features.tsv", "features.tsv.gz",
        "genes.tsv", "genes.tsv.gz"
      )
    )))
    matrix_ok && barcode_ok && feature_ok
  }
  all_matrix_dirs <- c(
    manifest$FilteredMatrixDir,
    manifest$RawMatrixDir
  )
  incomplete_dirs <- all_matrix_dirs[
    !vapply(all_matrix_dirs, has_10x_files, logical(1))
  ]
  if (length(incomplete_dirs) > 0L) {
    stop(
      "Incomplete 10X matrix triplets:\n",
      paste(unique(incomplete_dirs), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ============================================================================
# 单样本处理（scDblFinder + SoupX + QC）
# ============================================================================

# 计算 filtered 矩阵的基础 QC 指标（UMI 数、基因数、线粒体比例），保留原始 barcode。
matrix_qc <- function(counts, raw_barcodes = colnames(counts)) {
  if (length(raw_barcodes) != ncol(counts)) {
    stop("raw_barcodes must contain one value per matrix column.", call. = FALSE)
  }
  mt_genes <- grepl("^mt-|^MT-", rownames(counts))
  n_count <- as.numeric(Matrix::colSums(counts))
  n_feature <- as.numeric(Matrix::colSums(counts > 0))
  mt_count <- if (any(mt_genes)) {
    as.numeric(Matrix::colSums(counts[mt_genes, , drop = FALSE]))
  } else {
    rep(0, ncol(counts))
  }
  percent_mt <- ifelse(n_count > 0, 100 * mt_count / n_count, NA_real_)

  out <- data.frame(
    RawBarcode = raw_barcodes,
    raw_nCount_RNA = n_count,
    raw_nFeature_RNA = n_feature,
    raw_percent.mt = percent_mt,
    stringsAsFactors = FALSE,
    row.names = colnames(counts)
  )
  out
}

# 把逐细胞 QC 失败标记折叠为分号连接的原因字符串；全部通过记为 "pass"。
collapse_qc_reasons <- function(flags) {
  apply(flags, 1L, function(x) {
    failed <- sub("^fail_", "", colnames(flags)[as.logical(x)])
    if (length(failed) == 0L) "pass" else paste(failed, collapse = ";")
  })
}

# 处理单个样本的完整流程：
#   读取 raw/filter 矩阵 → 完整性检查 → 轻预过滤 → scDblFinder（原始整数计数）
#   → SoupX（复用 scDblFinder cluster，autoEstCont 或手动 rho）→ 校验校正矩阵
#   → 建 Seurat 对象并写入全部审计元数据 → 组织特异 QC 打标 → 按 QC 过滤。
# 返回过滤后对象、样本级 summary、raw QC 表和最终过滤前的完整元数据。
process_one_sample <- function(
    sample_row,
    sample_seed,
    pre_min_counts,
    pre_min_features,
    dbr_per_1000,
    dbr_sd,
    dbr_warning_difference,
    dbr_warning_fold,
    remove_doublets,
    thresholds) {
  if (!is.finite(dbr_per_1000) || dbr_per_1000 <= 0 ||
      !is.finite(dbr_sd) || dbr_sd < 0 ||
      !is.finite(dbr_warning_difference) || dbr_warning_difference < 0 ||
      !is.finite(dbr_warning_fold) || dbr_warning_fold <= 1) {
    stop("Invalid doublet-rate parameters.", call. = FALSE)
  }
  if (!is.logical(remove_doublets) ||
      length(remove_doublets) != 1L ||
      is.na(remove_doublets)) {
    stop("remove_doublets must be TRUE or FALSE.", call. = FALSE)
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }

  sample_id <- as.character(sample_row$SampleID[[1]])
  message("\n[", sample_id, "] Reading raw and filtered matrices")
  set.seed(sample_seed)

  filtered_counts <- read_10x_counts(sample_row$FilteredMatrixDir[[1]])
  raw_counts <- read_10x_counts(sample_row$RawMatrixDir[[1]])
  filtered_counts <- methods::as(filtered_counts, "CsparseMatrix")
  raw_counts <- methods::as(raw_counts, "CsparseMatrix")

  if (!identical(rownames(raw_counts), rownames(filtered_counts))) {
    stop("[", sample_id, "] Raw and filtered feature names/order differ.", call. = FALSE)
  }
  if (!all(colnames(filtered_counts) %in% colnames(raw_counts))) {
    stop("[", sample_id, "] Filtered barcodes are not a subset of raw barcodes.", call. = FALSE)
  }
  if (anyDuplicated(colnames(filtered_counts)) ||
      anyDuplicated(colnames(raw_counts))) {
    stop("[", sample_id, "] Matrix barcodes are not unique.", call. = FALSE)
  }
  if (any(raw_counts@x < 0) || any(filtered_counts@x < 0)) {
    stop("[", sample_id, "] Counts contain negative values.", call. = FALSE)
  }
  if (any(raw_counts@x != floor(raw_counts@x)) ||
      any(filtered_counts@x != floor(filtered_counts@x))) {
    stop("[", sample_id, "] scDblFinder requires original integer counts.", call. = FALSE)
  }

  # From this point onward, every cell ID is globally unique. Keep the original
  # vendor barcode separately for provenance and raw/filtered cross-checking.
  filtered_barcodes <- colnames(filtered_counts)
  raw_barcodes <- colnames(raw_counts)
  colnames(filtered_counts) <- paste(
    sample_id, filtered_barcodes, sep = "_"
  )
  colnames(raw_counts) <- paste(sample_id, raw_barcodes, sep = "_")

  raw_qc <- matrix_qc(
    filtered_counts,
    raw_barcodes = filtered_barcodes
  )
  pre_keep <- with(
    raw_qc,
    raw_nCount_RNA >= pre_min_counts &
      raw_nFeature_RNA >= pre_min_features
  )
  pre_keep[is.na(pre_keep)] <- FALSE
  raw_qc$pass_pre_qc <- pre_keep
  if (!any(pre_keep)) {
    stop("[", sample_id, "] No cells passed the minimal pre-filter.", call. = FALSE)
  }
  toc <- filtered_counts[, pre_keep, drop = FALSE]

  message(
    "[", sample_id, "] scDblFinder on ", ncol(toc),
    " cells using original integer counts"
  )
  expressed_genes <- Matrix::rowSums(toc) > 0
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = toc[expressed_genes, , drop = FALSE])
  )
  bp <- BiocParallel::SerialParam(
    RNGseed = sample_seed,
    progressbar = FALSE
  )

  manual_dbr <- suppressWarnings(as.numeric(sample_row$ExpectedDoubletRate[[1]]))
  dbl_args <- list(
    sce = sce,
    clusters = TRUE,
    dbr.sd = dbr_sd,
    returnType = "sce",
    verbose = FALSE,
    BPPARAM = bp
  )
  if (is.finite(manual_dbr)) {
    if (manual_dbr <= 0 || manual_dbr >= 1) {
      stop(
        "[", sample_id,
        "] ExpectedDoubletRate must be in (0, 1).",
        call. = FALSE
      )
    }
    dbl_args$dbr <- manual_dbr
    expected_dbr <- manual_dbr
    doublet_rate_method <- "per-sample override"
  } else {
    dbl_args$dbr.per1k <- dbr_per_1000
    expected_dbr <- dbr_per_1000 * ncol(toc) / 1000
    doublet_rate_method <- "dbr.per1k"
  }
  if (!is.finite(expected_dbr) || expected_dbr <= 0 || expected_dbr >= 1) {
    stop("[", sample_id, "] Invalid expected doublet rate.", call. = FALSE)
  }
  sce <- do.call(scDblFinder::scDblFinder, dbl_args)
  dbl_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  required_dbl_cols <- c(
    "scDblFinder.score", "scDblFinder.class", "scDblFinder.cluster"
  )
  missing_dbl_cols <- setdiff(required_dbl_cols, colnames(dbl_meta))
  if (length(missing_dbl_cols) > 0L) {
    stop(
      "[", sample_id, "] scDblFinder output is missing: ",
      paste(missing_dbl_cols, collapse = ", "),
      call. = FALSE
    )
  }
  dbl_meta <- dbl_meta[colnames(toc), required_dbl_cols, drop = FALSE]
  if (any(!is.finite(dbl_meta$scDblFinder.score)) ||
      anyNA(dbl_meta$scDblFinder.class) ||
      anyNA(dbl_meta$scDblFinder.cluster) ||
      !all(dbl_meta$scDblFinder.class %in% c("singlet", "doublet"))) {
    stop("[", sample_id, "] Invalid scDblFinder score/class/cluster output.", call. = FALSE)
  }
  observed_dbr <- mean(dbl_meta$scDblFinder.class == "doublet")
  if (abs(observed_dbr - expected_dbr) > dbr_warning_difference ||
      observed_dbr > dbr_warning_fold * expected_dbr) {
    warning(
      "[", sample_id, "] scDblFinder called ",
      round(100 * observed_dbr, 2), "% doublets versus ",
      round(100 * expected_dbr, 2),
      "% expected. Review the score distribution and rate assumption.",
      call. = FALSE
    )
  }

  message("[", sample_id, "] SoupX contamination estimation")
  soup_channel <- SoupX::SoupChannel(tod = raw_counts, toc = toc)
  soup_clusters <- stats::setNames(
    as.character(dbl_meta$scDblFinder.cluster),
    rownames(dbl_meta)
  )
  soup_channel <- SoupX::setClusters(
    soup_channel,
    soup_clusters[colnames(toc)]
  )

  manual_rho <- suppressWarnings(as.numeric(sample_row$SoupXManualRho[[1]]))
  soupx_rho_fwhm <- c(NA_real_, NA_real_)
  soupx_independent_estimates <- NA_integer_
  soupx_genes_used <- NA_integer_
  if (is.finite(manual_rho)) {
    if (manual_rho < 0 || manual_rho > 1) {
      stop("[", sample_id, "] SoupXManualRho must be in [0, 1].", call. = FALSE)
    }
    soup_channel <- SoupX::setContaminationFraction(soup_channel, manual_rho)
    soupx_method <- "manual override"
    estimated_rho <- NA_real_
  } else {
    soup_channel <- tryCatch(
      SoupX::autoEstCont(
        soup_channel,
        doPlot = FALSE,
        # verbose 必须为 FALSE：multisession worker 的 message 需经 socket 中继回
        # 主进程，多 worker 交叉输出会造成中继死锁（2026-08-03 排查确认）。
        verbose = FALSE
      ),
      error = function(e) {
        stop(
          "[", sample_id, "] SoupX autoEstCont failed: ",
          conditionMessage(e),
          ". Inspect this sample; do not silently substitute a fixed rho.",
          call. = FALSE
        )
      }
    )
    soupx_method <- "autoEstCont"
    candidate_rho <- as.numeric(soup_channel$fit$rhoEst)
    estimated_rho <- if (
      length(candidate_rho) == 1L && is.finite(candidate_rho)
    ) {
      candidate_rho
    } else {
      NA_real_
    }
    soupx_fit <- soup_channel$fit
    if (is.list(soupx_fit)) {
      candidate_fwhm <- as.numeric(soupx_fit$rhoFWHM)
      if (length(candidate_fwhm) == 2L &&
          all(is.finite(candidate_fwhm))) {
        soupx_rho_fwhm <- candidate_fwhm
      }
      if (is.data.frame(soupx_fit$dd) &&
          all(c("gene", "useEst") %in% colnames(soupx_fit$dd))) {
        used_estimates <- soupx_fit$dd$useEst %in% TRUE
        soupx_independent_estimates <- sum(used_estimates)
        soupx_genes_used <- length(unique(
          soupx_fit$dd$gene[used_estimates]
        ))
      }
    }
  }
  used_rho <- unique(as.numeric(soup_channel$metaData$rho))
  used_rho <- used_rho[is.finite(used_rho)]
  if (length(used_rho) != 1L || used_rho < 0 || used_rho > 1) {
    stop("[", sample_id, "] SoupX did not produce one valid contamination rate.", call. = FALSE)
  }

  set.seed(sample_seed)
  adjusted_counts <- SoupX::adjustCounts(
    soup_channel,
    method = "subtraction",
    roundToInt = TRUE,
    verbose = 0
  )
  adjusted_counts <- methods::as(adjusted_counts, "CsparseMatrix")
  if (!identical(dim(adjusted_counts), dim(toc)) ||
      !identical(rownames(adjusted_counts), rownames(toc)) ||
      !identical(colnames(adjusted_counts), colnames(toc))) {
    stop("[", sample_id, "] SoupX changed matrix dimensions or names.", call. = FALSE)
  }
  if (any(!is.finite(adjusted_counts@x)) ||
      any(adjusted_counts@x < 0) ||
      any(adjusted_counts@x != floor(adjusted_counts@x))) {
    stop("[", sample_id, "] SoupX output is not a non-negative integer matrix.", call. = FALSE)
  }
  if (sum(adjusted_counts) > sum(toc) + 1e-6) {
    stop("[", sample_id, "] SoupX-adjusted UMI total exceeds the input total.", call. = FALSE)
  }
  if (any(Matrix::colSums(adjusted_counts) == 0)) {
    stop("[", sample_id, "] SoupX produced one or more zero-count cells.", call. = FALSE)
  }

  cell_ids <- colnames(adjusted_counts)
  original_barcodes <- raw_qc[cell_ids, "RawBarcode"]

  obj <- Seurat::CreateSeuratObject(
    counts = adjusted_counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  obj <- add_basic_metadata(obj, sample_row)

  cell_metadata <- data.frame(
    RawBarcode = original_barcodes,
    raw_nCount_RNA = raw_qc[cell_ids, "raw_nCount_RNA"],
    raw_nFeature_RNA = raw_qc[cell_ids, "raw_nFeature_RNA"],
    raw_percent.mt = raw_qc[cell_ids, "raw_percent.mt"],
    pass_pre_qc = raw_qc[cell_ids, "pass_pre_qc"],
    scDblFinder.score = dbl_meta[cell_ids, "scDblFinder.score"],
    scDblFinder.class = as.character(
      dbl_meta[cell_ids, "scDblFinder.class"]
    ),
    scDblFinder.cluster = as.character(
      dbl_meta[cell_ids, "scDblFinder.cluster"]
    ),
    DoubletRemovalApplied = remove_doublets,
    SoupXMethod = soupx_method,
    SoupXEstimatedRho = estimated_rho,
    SoupXUsedRho = used_rho[[1]],
    stringsAsFactors = FALSE,
    row.names = cell_ids
  )
  obj <- Seurat::AddMetaData(obj, metadata = cell_metadata)
  obj <- add_qc_metrics(obj)

  md <- obj[[]]
  is_doublet <- md$scDblFinder.class == "doublet"
  qc_flags <- data.frame(
    fail_low_features = md$nFeature_RNA < thresholds$min_features,
    fail_high_features = md$nFeature_RNA > thresholds$max_features,
    fail_low_counts = md$nCount_RNA < thresholds$min_counts,
    fail_high_counts = if (is.finite(thresholds$max_counts)) {
      md$nCount_RNA > thresholds$max_counts
    } else {
      rep(FALSE, nrow(md))
    },
    fail_high_mt = md$percent.mt > thresholds$max_percent_mt,
    fail_doublet = is_doublet & remove_doublets,
    row.names = rownames(md)
  )
  qc_flags[is.na(qc_flags)] <- TRUE
  pass_expression_qc <- !apply(
    qc_flags[, setdiff(colnames(qc_flags), "fail_doublet"), drop = FALSE],
    1L,
    any
  )
  pass_final_qc <- pass_expression_qc & !qc_flags$fail_doublet
  doublets_among_expression_pass <- sum(
    is_doublet & pass_expression_qc
  )
  doublet_expression_fail_overlap <- sum(
    is_doublet & !pass_expression_qc
  )
  qc_metadata <- cbind(
    qc_flags,
    pass_expression_qc = pass_expression_qc,
    pass_final_qc = pass_final_qc,
    qc_failure_reason = collapse_qc_reasons(qc_flags)
  )
  obj <- Seurat::AddMetaData(obj, metadata = qc_metadata)

  cells_after_final_qc <- sum(pass_final_qc)
  if (cells_after_final_qc == 0L) {
    stop("[", sample_id, "] No cells passed final QC.", call. = FALSE)
  }

  summary <- data.frame(
    SampleID = sample_id,
    MouseID = as.character(sample_row$MouseID[[1]]),
    Tissue = as.character(sample_row$Tissue[[1]]),
    Group = as.character(sample_row$Group[[1]]),
    RawDroplets = ncol(raw_counts),
    FilteredCells = ncol(filtered_counts),
    CellsAfterPreQC = ncol(toc),
    PreQCRemoved = ncol(filtered_counts) - ncol(toc),
    ExpectedDoubletRate = expected_dbr,
    DoubletRateSD = dbr_sd,
    DoubletRateMethod = doublet_rate_method,
    DoubletRemovalApplied = remove_doublets,
    scDblFinderVersion = as.character(
      utils::packageVersion("scDblFinder")
    ),
    DoubletsCalled = sum(dbl_meta$scDblFinder.class == "doublet"),
    DoubletPercent = 100 * observed_dbr,
    DoubletsAmongExpressionPass = doublets_among_expression_pass,
    DoubletExpressionFailOverlap = doublet_expression_fail_overlap,
    DoubletsRemoved = sum(qc_flags$fail_doublet & pass_expression_qc),
    SoupXMethod = soupx_method,
    SoupXEstimatedRho = estimated_rho,
    SoupXUsedRho = used_rho[[1]],
    SoupXRhoFWHMLower = soupx_rho_fwhm[[1]],
    SoupXRhoFWHMUpper = soupx_rho_fwhm[[2]],
    SoupXIndependentEstimates = soupx_independent_estimates,
    SoupXGenesUsed = soupx_genes_used,
    SoupXVersion = as.character(utils::packageVersion("SoupX")),
    InputUMIs = as.numeric(sum(toc)),
    AdjustedUMIs = as.numeric(sum(adjusted_counts)),
    CellsPassingExpressionQC = sum(pass_expression_qc),
    CellsAfterFinalQC = cells_after_final_qc,
    FinalRetentionPercent = 100 * cells_after_final_qc / ncol(filtered_counts),
    stringsAsFactors = FALSE
  )

  raw_qc$SampleID <- sample_id
  raw_qc$MouseID <- as.character(sample_row$MouseID[[1]])
  raw_qc$Tissue <- as.character(sample_row$Tissue[[1]])
  raw_qc$Group <- as.character(sample_row$Group[[1]])
  raw_qc$CellID <- rownames(raw_qc)

  qc_metadata_before_final <- obj[[]]
  qc_metadata_before_final$CellID <- rownames(qc_metadata_before_final)
  filtered_obj <- subset(obj, cells = rownames(md)[pass_final_qc])

  rm(
    raw_counts, filtered_counts, toc, sce, soup_channel,
    adjusted_counts, obj
  )
  gc()

  list(
    obj = filtered_obj,
    summary = summary,
    raw_qc = raw_qc,
    qc_before_final = qc_metadata_before_final
  )
}

# ============================================================================
# 汇总与审计
# ============================================================================

# 设置外层样本级并行，workers 取 n_workers 与可用核数的较小值。
# backend 由 config/pipeline_params.R 的 parallel_backend 决定：
#   "multicore"    —— fork（Linux/macOS 推荐；与 01/02/08 的后端一致）
#   "multisession" —— socket 集群（当前 WSL2 环境下对长任务收集结果会挂死，
#                    2026-08-03 以 base R parallel::makeCluster("PSOCK") 复现确认，勿用）
#   "sequential"   —— 串行
# fork 不可用（典型场景：RStudio Console/Knit；R Core 与 parallelly 均不建议在
# GUI 前端 fork）时回退 sequential 并告警——绝不回退 multisession。2026-08-03 核实：
# 旧实现在此静默落入 multisession（本主机已知挂死后端）；且 future.fork.enable
# 必须在 supportsMulticore() 检查之前设置才有效，与其依赖该选项在 RStudio 强行
# fork，不如回退串行并提示改用 headless Rscript。plan 启动失败同样回退串行。
# 返回实际 workers 数和生效的 backend 名称。
setup_sample_parallel <- function(n_workers, backend = "multicore") {
  available_workers <- max(1L, as.integer(future::availableCores()[[1]]))
  requested_workers <- min(as.integer(n_workers), available_workers)
  workers <- 1L
  active_backend <- "sequential"
  if (requested_workers > 1L && identical(backend, "multicore") &&
      !isTRUE(future::supportsMulticore())) {
    warning(
      paste0(
        "multicore (fork) is not supported in this R session ",
        "(typical cause: running inside RStudio Console/Knit). ",
        "Falling back to sequential; multisession is never used as a fallback ",
        "because PSOCK result collection hangs on this WSL2 host. ",
        "Re-run via headless Rscript to get parallel execution."
      ),
      call. = FALSE
    )
    backend <- "sequential"
  }
  if (requested_workers > 1L && !identical(backend, "sequential")) {
    tryCatch(
      {
        if (identical(backend, "multicore")) {
          future::plan(future::multicore, workers = requested_workers)
        } else {
          future::plan(future::multisession, workers = requested_workers)
        }
        workers <- requested_workers
        active_backend <- backend
      },
      error = function(e) {
        warning(
          "Could not start parallel workers; continuing sequentially: ",
          conditionMessage(e),
          call. = FALSE
        )
        future::plan(future::sequential)
      }
    )
  } else {
    future::plan(future::sequential)
  }
  list(workers = workers, backend = active_backend)
}

# 合并 raw 阶段与 SoupX 后阶段的元数据，得到逐细胞 QC 审计表：
# 未过预过滤的细胞标记为 pre_qc_low_coverage，缺失标记统一回填。
build_cell_qc_audit <- function(raw_qc_metadata, qc_metadata_before_final, remove_doublets) {
  audit_after_pre_qc <- dplyr::select(
    qc_metadata_before_final,
    CellID,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    percent.ribo,
    scDblFinder.score,
    scDblFinder.class,
    scDblFinder.cluster,
    DoubletRemovalApplied,
    SoupXMethod,
    SoupXEstimatedRho,
    SoupXUsedRho,
    dplyr::starts_with("fail_"),
    pass_expression_qc,
    pass_final_qc,
    qc_failure_reason
  )
  cell_qc_audit <- dplyr::left_join(
    raw_qc_metadata,
    audit_after_pre_qc,
    by = "CellID"
  )
  cell_qc_audit$pass_final_qc[
    is.na(cell_qc_audit$pass_final_qc)
  ] <- FALSE
  cell_qc_audit$DoubletRemovalApplied[
    is.na(cell_qc_audit$DoubletRemovalApplied)
  ] <- remove_doublets
  cell_qc_audit$qc_failure_reason[
    !cell_qc_audit$pass_pre_qc
  ] <- "pre_qc_low_coverage"
  cell_qc_audit
}

# 把 QC 阈值列表转成可写表的 data.frame。
thresholds_to_data_frame <- function(thresholds, tissue) {
  data.frame(
    Tissue = tissue,
    min_features = thresholds$min_features,
    max_features = thresholds$max_features,
    min_counts = thresholds$min_counts,
    max_counts = thresholds$max_counts,
    max_percent_mt = thresholds$max_percent_mt,
    stringsAsFactors = FALSE
  )
}

# 写出本步骤全部 QC 表：样本级 summary、逐细胞审计（csv.gz）、SoupX/doublet
# 专项表、分阶段细胞数、QC 阈值和失败原因统计。
write_qc_tables <- function(sample_summary, cell_qc_audit, qc_thresholds, table_dir) {
  write.csv(
    sample_summary,
    file.path(table_dir, "00_qc_summary_by_sample.csv"),
    row.names = FALSE
  )

  cell_audit_connection <- gzfile(
    file.path(table_dir, "00_cell_qc_audit.csv.gz"),
    open = "wt"
  )
  tryCatch(
    write.csv(
      cell_qc_audit,
      cell_audit_connection,
      row.names = FALSE
    ),
    finally = close(cell_audit_connection)
  )

  write.csv(
    dplyr::select(
      sample_summary,
      SampleID, MouseID, Tissue, Group, RawDroplets, FilteredCells,
      CellsAfterPreQC, SoupXMethod, SoupXEstimatedRho, SoupXUsedRho,
      SoupXRhoFWHMLower, SoupXRhoFWHMUpper,
      SoupXIndependentEstimates, SoupXGenesUsed, SoupXVersion,
      InputUMIs, AdjustedUMIs
    ),
    file.path(table_dir, "00_soupx_per_sample_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    dplyr::select(
      sample_summary,
      SampleID, MouseID, Tissue, Group, CellsAfterPreQC,
      ExpectedDoubletRate, DoubletRateSD, DoubletRateMethod,
      DoubletRemovalApplied, DoubletsCalled, DoubletPercent,
      DoubletsAmongExpressionPass, DoubletExpressionFailOverlap,
      DoubletsRemoved, scDblFinderVersion
    ),
    file.path(table_dir, "00_doublet_summary_by_sample.csv"),
    row.names = FALSE
  )

  qc_stage_summary <- sample_summary |>
    dplyr::select(
      SampleID, MouseID, Tissue, Group,
      FilteredCells, CellsAfterPreQC,
      CellsPassingExpressionQC, CellsAfterFinalQC
    ) |>
    tidyr::pivot_longer(
      cols = c(
        "FilteredCells", "CellsAfterPreQC",
        "CellsPassingExpressionQC", "CellsAfterFinalQC"
      ),
      names_to = "Stage",
      values_to = "Cells"
    ) |>
    dplyr::mutate(
      StageOrder = match(
        Stage,
        c(
          "FilteredCells",
          "CellsAfterPreQC",
          "CellsPassingExpressionQC",
          "CellsAfterFinalQC"
        )
      )
    ) |>
    dplyr::arrange(SampleID, StageOrder) |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(
      RemovedFromPreviousStage = dplyr::lag(Cells) - Cells,
      RetentionPercentOfFiltered = 100 * Cells / dplyr::first(Cells)
    ) |>
    dplyr::ungroup()
  write.csv(
    qc_stage_summary,
    file.path(table_dir, "00_qc_filter_summary_by_sample.csv"),
    row.names = FALSE
  )

  write.csv(
    thresholds_to_data_frame(qc_thresholds$Spleen, "Spleen"),
    file.path(table_dir, "00_qc_thresholds_spleen.csv"),
    row.names = FALSE
  )
  write.csv(
    thresholds_to_data_frame(qc_thresholds$BoneMarrow, "BoneMarrow"),
    file.path(table_dir, "00_qc_thresholds_bonemarrow.csv"),
    row.names = FALSE
  )

  qc_failure_counts <- cell_qc_audit |>
    dplyr::count(
      SampleID,
      Tissue,
      Group,
      qc_failure_reason,
      name = "Cells"
    )
  write.csv(
    qc_failure_counts,
    file.path(table_dir, "00_qc_failure_reasons_by_sample.csv"),
    row.names = FALSE
  )
  invisible(TRUE)
}

# ============================================================================
# QC 图
# ============================================================================

# 按样本绘制 QC 指标小提琴图（nCount / nFeature / percent.mt 三个 facet）。
plot_qc_metrics <- function(
    meta,
    count_col,
    feature_col,
    mt_col,
    title) {
  plot_data <- data.frame(
    SampleID = meta$SampleID,
    Tissue = meta$Tissue,
    nCount = meta[[count_col]],
    nFeature = meta[[feature_col]],
    percent.mt = meta[[mt_col]],
    stringsAsFactors = FALSE
  )
  plot_data |>
    tidyr::pivot_longer(
      cols = c("nCount", "nFeature", "percent.mt"),
      names_to = "Metric",
      values_to = "Value"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = SampleID, y = Value, fill = Tissue)
    ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.15
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(Metric),
      scales = "free_y",
      nrow = 1
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_publication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# 构建本步骤全部 QC 图，返回命名列表：
# raw_qc（校正前小提琴）、final_qc（校正后小提琴）、raw_scatter（UMI vs 基因数散点）、
# doublet_scores（scDblFinder 分数分布）、method_summary（rho 与 doublet 率柱状图）。
build_qc_figures <- function(raw_qc_metadata, final_metadata, qc_metadata_before_final, sample_summary) {
  p_raw_qc <- plot_qc_metrics(
    raw_qc_metadata,
    count_col = "raw_nCount_RNA",
    feature_col = "raw_nFeature_RNA",
    mt_col = "raw_percent.mt",
    title = "Raw filtered-matrix QC metrics before correction"
  )
  p_final_qc <- plot_qc_metrics(
    final_metadata,
    count_col = "nCount_RNA",
    feature_col = "nFeature_RNA",
    mt_col = "percent.mt",
    title = "SoupX-corrected metrics after configured final QC"
  )

  p_raw_scatter <- ggplot2::ggplot(
    raw_qc_metadata,
    ggplot2::aes(
      x = raw_nCount_RNA,
      y = raw_nFeature_RNA,
      color = Group
    )
  ) +
    ggplot2::geom_point(alpha = 0.15, size = 0.25) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(ggplot2::vars(SampleID), scales = "free") +
    ggplot2::labs(
      title = "Raw UMI count versus detected features",
      x = "Raw UMI count (log10)",
      y = "Raw detected features (log10)",
      color = "Group"
    ) +
    theme_publication()

  p_doublet_scores <- ggplot2::ggplot(
    qc_metadata_before_final,
    ggplot2::aes(
      x = scDblFinder.score,
      fill = scDblFinder.class
    )
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      position = "identity",
      alpha = 0.55
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(SampleID),
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "Per-sample scDblFinder score distributions",
      x = "scDblFinder score",
      y = "Cells",
      fill = "Call"
    ) +
    theme_publication()

  method_plot_data <- sample_summary |>
    dplyr::transmute(
      SampleID,
      Tissue,
      `SoupX rho (%)` = 100 * SoupXUsedRho,
      `Doublets called (%)` = DoubletPercent
    ) |>
    tidyr::pivot_longer(
      cols = c("SoupX rho (%)", "Doublets called (%)"),
      names_to = "Metric",
      values_to = "Percent"
    )
  p_method_summary <- ggplot2::ggplot(
    method_plot_data,
    ggplot2::aes(x = SampleID, y = Percent, fill = Tissue)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(
      ggplot2::vars(Metric),
      scales = "free_y",
      nrow = 1
    ) +
    ggplot2::labs(
      title = "Per-sample SoupX and scDblFinder estimates",
      x = NULL,
      y = "Percent"
    ) +
    theme_publication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  list(
    raw_qc = p_raw_qc,
    final_qc = p_final_qc,
    raw_scatter = p_raw_scatter,
    doublet_scores = p_doublet_scores,
    method_summary = p_method_summary
  )
}

# ============================================================================
# 合并与校验
# ============================================================================

# 合并同组织的样本对象；只有一个样本时直接返回。
merge_sample_objects <- function(objects, project) {
  if (length(objects) == 0L) {
    stop("No objects supplied for ", project, ".", call. = FALSE)
  }
  if (length(objects) == 1L) {
    objects[[1]]
  } else {
    merge(
      x = objects[[1]],
      y = objects[-1],
      project = project
    )
  }
}

# 校验合并后的组织对象：关键元数据齐全、细胞 ID 唯一且为 SampleID_RawBarcode 格式、
# 无缺失值、全部通过最终 QC、包含 RNA assay、样本集合和细胞总数与预期一致。
validate_merged_object <- function(obj, expected_samples, expected_cells, tissue) {
  required_metadata <- c(
    "SampleID", "MouseID", "Group", "Tissue", "RawBarcode",
    "scDblFinder.score", "scDblFinder.class", "DoubletRemovalApplied",
    "SoupXUsedRho", "pass_final_qc"
  )
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata) > 0L) {
    stop(
      tissue, " object is missing metadata: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(colnames(obj))) {
    stop(tissue, " object contains duplicate cell IDs.", call. = FALSE)
  }
  expected_cell_ids <- paste(obj$SampleID, obj$RawBarcode, sep = "_")
  if (!identical(colnames(obj), expected_cell_ids)) {
    stop(tissue, " cell IDs do not match SampleID_RawBarcode.", call. = FALSE)
  }
  if (anyNA(obj[[]][, required_metadata, drop = FALSE])) {
    stop(tissue, " object contains missing critical metadata.", call. = FALSE)
  }
  if (!all(obj$pass_final_qc)) {
    stop(tissue, " object contains cells that did not pass final QC.", call. = FALSE)
  }
  if (!"RNA" %in% SeuratObject::Assays(obj)) {
    stop(tissue, " object is missing the RNA assay.", call. = FALSE)
  }
  if (!setequal(unique(obj$SampleID), expected_samples)) {
    stop(tissue, " object does not contain the expected six samples.", call. = FALSE)
  }
  if (ncol(obj) != expected_cells) {
    stop(
      tissue, " merged cell count does not equal the sum of sample objects.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
