# pipeline_params.R — snRNA-seq publication pipeline 集中参数文件
#
# 关键参数的唯一定义处：调整流程参数只改本文件，不要改各 Rmd。
# 00–08 各 Rmd 在 configuration chunk 开头 source() 本文件；render_all.R 同样读取。
# 修改本文件属于行为性修改，须同步 README.md 实时总结（见 AGENT.md 维护规则）。

## ===== 全局 =====
seed <- 20260717L                       # 随机种子（00 逐样本派生 scDblFinder/SoupX 种子）
n_workers <- 6L                         # 外层并行 workers 上限（00 样本循环、SingleR）
cellchat_workers <- 2L                  # 08 CellChat 专用 workers：峰值内存随 worker 数线性叠加，
                                        # 2026-08-06 用全局 n_workers=6 时 computeCommunProb 把
                                        # WSL 82GB 内存吃光触发 OOM（宿主机仅 94GB），降为 2
parallel_backend <- "multicore"         # 样本级并行后端：multicore(fork，推荐) / multisession / sequential。
                                        # 当前环境为 WSL2，PSOCK(socket) 集群对长任务收集结果会挂死
                                        # （2026-08-03 用 base R parallel::makeCluster("PSOCK") 复现确认），
                                        # 因此禁用 multisession；multicore 走 fork，不受影响且更快。
future_globals_maxsize <- 88 * 1024^3   # future.globals.maxSize（SCTransform/CellChat 大对象）

## ===== 00 数据读取 / scDblFinder / SoupX / QC =====
pre_qc_min_counts <- 200                # 轻预过滤：仅去除空滴/极低复杂度 barcode
pre_qc_min_features <- 20
doublet_rate_per_1000 <- 0.004          # DNBelab C4 经验双联体率：每 1000 细胞 0.4%
doublet_rate_sd <- 0                    # scDblFinder dbr.sd；0 = 强先验
doublet_rate_warning_difference <- 0.05 # 观察率 vs 期望率告警阈值（绝对差）
doublet_rate_warning_fold <- 2          # 观察率 vs 期望率告警阈值（倍数）
remove_scDblFinder_doublets <- FALSE    # first-pass 默认只标记不剔除；复核分数分布后再开

# max_percent_mt 两组织统一 0.5%（2026-08-10 用户拍板：所有样本 mt<0.5%；原脾 5%/髓 3%）
qc_thresholds <- list(
  Spleen = list(
    min_features = 100, max_features = 3000,
    min_counts = 500, max_counts = Inf,
    max_percent_mt = 0.5
  ),
  BoneMarrow = list(
    min_features = 100, max_features = 4500,
    min_counts = 500, max_counts = Inf,
    max_percent_mt = 0.5
  )
)

## ===== 01 / 02 聚类与注释 =====
harmony_batch_var <- NULL        # 真正的技术批次列名；NULL = 不跑 Harmony。禁止 SampleID/MouseID
dims_main <- 1:30                # 主聚类 PCA 维度
dims_subset <- 1:20              # compartment 重聚类 PCA 维度
cluster_res_main <- 0.6          # 主聚类分辨率
cluster_res_spleen_subset <- 0.3 # 脾 B/T compartment 重聚类分辨率
cluster_res_bm_subset <- 0.5     # 骨髓 B lineage compartment 重聚类分辨率
min_cells_for_refinement <- 100  # compartment 细胞数下限，低于则停止并提示检查 SingleR 标签

spleen_b_markers <- c(
  "Cd19", "Ms4a1", "Cd79a", "Ighd", "Ighm", "Fas", "Bcl6", "Aicda",
  "Mki67", "Jchain", "Irf4", "Prdm1", "Sdc1", "Ighg1", "Cd80"
)

spleen_t_markers <- c(
  "Cd3d", "Cd3e", "Cd4", "Cd8a", "Sell", "Ccr7", "Cd44",
  "Cxcr5", "Bcl6", "Icos", "Pdcd1", "Gzmb", "Prf1",
  "Ifng", "Tbx21", "Foxp3", "Il2ra", "Ctla4"
)

bm_blineage_markers <- c(
  "Vpreb1", "Dntt", "Cd79a", "Cd19", "Ms4a1", "Ighm", "Ighd",
  "Irf4", "Prdm1", "Xbp1", "Sdc1", "Tnfrsf17", "Jchain",
  "Mki67", "Mcl1", "Slc3a2"
)

## ===== 03 marker 复核 panel =====
spleen_validation_markers <- c(
  "Ptprc", "Cd19", "Ms4a1", "Cd79a", "Ighd", "Ighm",
  "Fas", "Bcl6", "Aicda", "Mki67", "Jchain", "Irf4", "Prdm1", "Sdc1",
  "Cd3d", "Cd3e", "Cd4", "Cd8a", "Sell", "Ccr7", "Cxcr5", "Icos", "Pdcd1",
  "Foxp3", "Il2ra", "Gzmb", "Prf1", "Nkg7", "Lyz2", "Itgam", "Cst3"
)

bm_validation_markers <- c(
  "Ptprc", "Vpreb1", "Dntt", "Cd79a", "Cd19", "Ms4a1", "Ighm", "Ighd",
  "Irf4", "Prdm1", "Xbp1", "Sdc1", "Tnfrsf17", "Jchain", "Mcl1", "Slc3a2",
  "Mki67", "Cd3d", "Cd3e", "Nkg7", "Lyz2", "S100a8", "S100a9", "Csf1r", "Cst3"
)

## ===== 04 细胞比例 =====
# focus 列表之外的类型在图中合并为 "Other cells"；CellType_Fine 统计不受影响。
spleen_focus <- c(
  "Naive B cells", "Activated or transitional B cells", "Memory B cell-like",
  "Early GC or pre-GC B cells", "GC B cells", "Cycling B cells", "Plasma cell-like",
  "Naive CD4 T cells", "Central memory CD4 T cells", "Tfh cells",
  "Regulatory T cells", "Naive CD8 T cells", "Effector memory CD8 T cells",
  "NK cells", "Monocytes", "Macrophages", "DC"
)

bm_focus <- c(
  "Pre or pro B cells", "Naive or mature B cells", "Activated B cells",
  "Plasmablast-like", "Plasma cell-like",
  "T cells", "NK cells", "Monocytes", "Macrophages", "Granulocytes",
  "Neutrophils", "Eosinophils", "Basophils", "DC"
)

## ===== 05 pseudobulk DEG =====
min_cells_per_mouse <- 20   # 每只 mouse 目标细胞数下限，低于则该 mouse 不参与 pseudobulk
fdr_cutoff <- 0.05          # 显著 DEG 判定阈值（05 与 06 共用）
logfc_cutoff <- 0.5

target_spleen <- c(
  "GC B cells",
  "Early GC or pre-GC B cells",
  "Plasma cell-like",
  "Memory B cell-like",
  "Tfh cells",
  "Effector memory CD8 T cells"
)

target_bm <- c(
  "Plasma cell-like",
  "Plasmablast-like",
  "Naive or mature B cells",
  "Pre or pro B cells",
  "Activated B cells"
)

## ===== 06 GO/KEGG 富集 =====
deg_method_for_enrichment <- "edgeR"  # 富集输入的 DEG 方法（05 结果中的方法名）
min_genes_for_enrichment <- 10        # 可映射 Entrez 基因数下限，低于则跳过
go_pvalue_cutoff <- 0.05
go_qvalue_cutoff <- 0.20
kegg_pvalue_cutoff <- 0.05

## ===== 07 publication 导出 =====
top_n_deg_per_celltype <- 25          # 合并 DEG 表中每个组织+细胞类型保留的 top DEG 数

## ===== 08 CellChat =====
cellchat_min_cells <- 10              # filterCommunication 的 min.cells
cellchat_spleen_sources <- c("Tfh cells", "Naive CD4 T cells", "Central memory CD4 T cells")
cellchat_spleen_targets <- c("GC B cells", "Plasma cell-like", "Naive B cells", "Memory B cell-like")
cellchat_bm_targets <- c("Plasma cell-like", "Plasmablast-like", "Naive or mature B cells")
