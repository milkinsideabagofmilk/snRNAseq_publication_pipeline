# pipeline_params.R — snRNA-seq publication pipeline central parameter file
#
# Single source of truth for key parameters: to adjust pipeline parameters, edit only this file, not the individual Rmds.
# Each Rmd in steps 00–08 calls source() on this file at the start of its configuration chunk; render_all.R reads it as well.
# Modifying this file is a behavioral change and must be accompanied by an up-to-date summary in README.md (see AGENT.md maintenance rules).

## ===== Global =====
seed <- 20260717L                       # random seed (00 derives per-sample scDblFinder/SoupX seeds from it)
n_workers <- 6L                         # upper limit on outer-level parallel workers (00 sample loop, SingleR)
cellchat_workers <- 2L                  # dedicated workers for 08 CellChat: peak memory scales linearly with worker count;
                                        # on 2026-08-06, with global n_workers=6, computeCommunProb exhausted the
                                        # WSL 82GB memory allocation and triggered OOM (host has only 94GB), so reduced to 2
parallel_backend <- "multicore"         # sample-level parallel backend: multicore (fork, recommended) / multisession / sequential.
                                        # The current environment is WSL2, where PSOCK (socket) clusters hang when
                                        # collecting results of long tasks (reproduced and confirmed on 2026-08-03
                                        # with base R parallel::makeCluster("PSOCK")), so multisession is disabled;
                                        # multicore uses fork and is unaffected and faster.
future_globals_maxsize <- 88 * 1024^3   # future.globals.maxSize (large objects in SCTransform/CellChat)

## ===== 00 Data loading / scDblFinder / SoupX / QC =====
pre_qc_min_counts <- 200                # light pre-filtering: removes only empty droplets / very-low-complexity barcodes
pre_qc_min_features <- 20
doublet_rate_per_1000 <- 0.004          # DNBelab C4 empirical doublet rate: 0.4% per 1000 cells
doublet_rate_sd <- 0                    # scDblFinder dbr.sd; 0 = strong prior
doublet_rate_warning_difference <- 0.05 # observed vs expected rate warning threshold (absolute difference)
doublet_rate_warning_fold <- 2          # observed vs expected rate warning threshold (fold change)
remove_scDblFinder_doublets <- FALSE    # first-pass default is to flag only, not remove; enable after reviewing the score distribution

# max_percent_mt unified at 0.5% for both tissues (decided by the user on 2026-08-10: mt<0.5% for all samples; previously 5% spleen / 3% marrow)
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

## ===== 01 / 02 Clustering and annotation =====
harmony_batch_var <- NULL        # true technical batch column name; NULL = do not run Harmony. SampleID/MouseID are forbidden
dims_main <- 1:30                # PCA dimensions for main clustering
dims_subset <- 1:20              # PCA dimensions for compartment re-clustering
cluster_res_main <- 0.6          # main clustering resolution
cluster_res_spleen_subset <- 0.3 # spleen B/T compartment re-clustering resolution
cluster_res_bm_subset <- 0.5     # bone marrow B lineage compartment re-clustering resolution
min_cells_for_refinement <- 100  # lower bound on compartment cell count; below this, stop and prompt to check SingleR labels

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

## ===== 03 marker validation panel =====
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

## ===== 04 Cell proportions =====
# Types outside the focus lists are merged into "Other cells" in plots; CellType_Fine statistics are unaffected.
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
min_cells_per_mouse <- 20   # lower bound on target cells per mouse; below this, the mouse is excluded from pseudobulk
fdr_cutoff <- 0.05          # significance threshold for DEG calls (shared by 05 and 06)
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

## ===== 06 GO/KEGG enrichment =====
deg_method_for_enrichment <- "edgeR"  # DEG method used as enrichment input (method name in the 05 results)
min_genes_for_enrichment <- 10        # lower bound on mappable Entrez genes; below this, skip
go_pvalue_cutoff <- 0.05
go_qvalue_cutoff <- 0.20
kegg_pvalue_cutoff <- 0.05

## ===== 07 publication export =====
top_n_deg_per_celltype <- 25          # number of top DEGs retained per tissue + cell type in the merged DEG table

## ===== 08 CellChat =====
cellchat_min_cells <- 10              # min.cells for filterCommunication
cellchat_spleen_sources <- c("Tfh cells", "Naive CD4 T cells", "Central memory CD4 T cells")
cellchat_spleen_targets <- c("GC B cells", "Plasma cell-like", "Naive B cells", "Memory B cell-like")
cellchat_bm_targets <- c("Plasma cell-like", "Plasmablast-like", "Naive or mature B cells")
