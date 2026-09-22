if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown before running render_all.R", call. = FALSE)
}

# render_all.R — 全量驱动器（00–07 顺序渲染）。
#
# 用法：
#   Rscript snRNAseq_publication_pipeline/render_all.R           # resume（默认）
#   Rscript snRNAseq_publication_pipeline/render_all.R --fresh   # 归档旧 outputs/ 后全量重跑
#
# 08 已拆出（2026-08-06）：CellChat 峰值内存巨大且不稳定，不再随全量流程跑，
# 单独用 snRNAseq_publication_pipeline/run_08_cellchat.sh（tmux + cgroup 笼子）启动。
#
# resume 模式：逐步判断"报告 HTML 是否新于全部输入"（Rmd 本身 + config +
# R/ 函数库 + annotation_maps 的 map 表），新则跳过；一旦某步重跑，下游全部
# 级联重跑（下游吃上游的 rds）。报告是渲染最后一刻才生成的，存在即代表该步
# 完整跑完。gate_decisions.csv 是决策日志、不是任何步骤的输入，不参与判定。
#
# 退出码契约：0 = 全部完成；10 = 手动注释 gate 触发（非错误，填 map 后重跑
# 本命令即可从断点续跑）；1 = 真报错。

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else "snRNAseq_publication_pipeline"
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

fresh <- "--fresh" %in% commandArgs(trailingOnly = TRUE)

# --fresh 时归档旧 outputs/，做版本隔离（同一文件系统下 rename 是瞬时操作，
# 不额外占用磁盘）。归档只发生在显式 --fresh；resume 与单个 Rmd 重跑都不触发。
# 运行结束后由人工根据新产物质量决定旧归档的去留（删除或手动移回）。
archive_outputs <- function(outputs_dir, archive_root) {
  if (!dir.exists(outputs_dir)) return(invisible(NULL))
  n_files <- length(list.files(outputs_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE))
  if (n_files == 0) return(invisible(NULL))
  dir.create(archive_root, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(archive_root, format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (dir.exists(dest) || !file.rename(outputs_dir, dest)) {
    stop("归档旧 outputs 失败：", outputs_dir, " -> ", dest, call. = FALSE)
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
  # 08_optional_cellchat.Rmd 不在此列：用 run_08_cellchat.sh 单独跑（见文件头注释）。
)

# resume 判定输入：公共输入 = 该 Rmd 本身 + config 参数 + R/ 函数库；注释 map
# 只作为 01/02 的输入（00 不读 map，03+ 经 rds 间接受 01/02 级联覆盖）。
# gate_decisions.csv 是决策日志、不是任何步骤的输入，不参与判定。
# 原始数据（dnbc4tools_results）约定为静态，数据更新请用 --fresh。
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
      message("\n==== 手动注释 gate 触发：", gate_id, "（非错误，等待注释决策） ====")
      message("状态文件（待标 cluster / 细胞数 / map 与模板路径）: ", status_file)
      message("标准流程：")
      message("  1) Rscript snRNAseq_publication_pipeline/run_gate_evidence.R ", gate_id)
      message("  2) 依据 outputs/gate_evidence/", gate_id, "/ 证据填写 map（模糊簇标 Unknown，标签复用既有词表）")
      message("  3) 在 annotation_maps/gate_decisions.csv 记录决策")
      message("  4) 重跑本命令——resume 会跳过已完成步骤，从本 gate 续跑")
      if (!interactive()) quit(save = "no", status = 10) else stop(msg, call. = FALSE)
    }
    message("Render failed: ", rmd, "\n", msg)
    if (!interactive()) quit(save = "no", status = 1) else stop(msg, call. = FALSE)
  }
  # 级联：本步重跑过，下游不能再按"报告仍新"跳过。
  rendered_upstream <- TRUE
}
message("render_all 完成：全部步骤已渲染。")
