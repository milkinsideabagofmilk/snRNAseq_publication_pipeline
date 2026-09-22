# snRNA-seq publication pipeline

## Live script-by-script summary (00–08)

This section corresponds one-to-one with the nine Rmd scripts 00–08 and the function libraries under `R/`; its content reflects the current code. Update this section in sync whenever the code changes. Key parameters are defined centrally in `config/pipeline_params.R`; the parameter values listed in each subsection match that file's current values.

### 00_data_loading_soupx_qc.Rmd

- **Purpose**: Reads dnbc4tools matrices sample by sample, performs doublet labeling, SoupX correction, and QC, then merges and saves by tissue. This script is a slim pipeline shell; heavy logic — input validation, per-sample processing, audit tables, and QC figure construction — is encapsulated in `R/step00_qc_helpers.R`.
- **Inputs**: `raw_matrix` and `filter_matrix` under `dnbc4tools_results/<SampleID>/<SampleID>/outs/` (12 samples; the manifest hard-validates the experimental design: 2 groups × 2 tissues × 3 replicates).
- **Key steps and parameters**: Light pre-filtering (≥200 UMI, ≥20 genes) removes only empty droplets; scDblFinder runs on raw integer counts before SoupX (expected rate 0.4%/1000 cells, `dbr.sd = 0`), and by default only flags doublets without removing them (`remove_scDblFinder_doublets = FALSE`); SoupX `autoEstCont` estimates the contamination fraction (the manifest can set `SoupXManualRho` to override; `verbose = FALSE` to avoid deadlocks caused by parallel-worker messages being relayed over sockets), reuses the scDblFinder clusters as soup clusters, and `adjustCounts` applies subtractive correction with rounding; tissue-specific QC (spleen: 100–3000 features, ≥500 UMI; bone marrow: 100–4500 features, ≥500 UMI; mt threshold unified at <0.5% for both tissues since 2026-08-10), with per-cell failure reasons recorded; 6 workers in parallel, with the backend determined by config `parallel_backend` (default multicore/fork); the outer `future_lapply` schedules per sample dynamically with `future.chunk.size = 1` (avoiding the default chunking of 2 samples/worker, which causes load imbalance and correlated failures).
- **Outputs**: `outputs/rds/00_spleen_qc.rds`, `00_bonemarrow_qc.rds`; CSVs for QC summary/thresholds/failure reasons/SoupX/doublets, QC figures, `00_RUN_STATUS.txt`, sessionInfo.
- **Downstream handoff**: 01 reads the spleen RDS; 02 reads the bone marrow RDS.

### 01_spleen_clustering_annotation.Rmd

- **Purpose**: Spleen cell clustering, SingleR main annotation, and refined B/T compartment annotation.
- **Inputs**: `outputs/rds/00_spleen_qc.rds`.
- **Key steps and parameters**: `set.seed(seed)` first fixes the random seed (guaranteeing reproducible cluster numbering, which `run_gate_evidence.R` depends on); main clustering + SingleR main annotation are performed via the shared function `run_main_clustering_annotation()` (SCTransform glmGamPoi regressing percent.mt → PCA → clustering dims 1:30 res 0.6 → UMAP; Harmony off by default; SingleR ImmGen majority vote per cluster generates `CellType_Main`); B and T compartments are extracted according to the label regexes registered in `gate_definitions()` and reclustered (`prepare_lineage_compartment()`, dims 1:20, resolution 0.3); manual-annotation gates are handled uniformly by `check_annotation_gate()` (writes the annotation template + validates the map; if not covered, writes `outputs/gate_evidence/gate_status_01_bcell.csv` / `gate_status_01_tcell.csv` and stops with the `GATE[01_bcell]` / `GATE[01_tcell]` prefix); marker dotplots and validation tables; fine labels written back to the main object.
- **Outputs**: `outputs/rds/01_spleen_annotated.rds`, `01_spleen_bcells_refined.rds`, `01_spleen_tcells_refined.rds`; UMAP plots, annotation count tables, marker validation tables.
- **Downstream handoff**: 03/04/05/08 read `01_spleen_annotated.rds`.

### 02_bm_clustering_annotation.Rmd

- **Purpose**: Bone marrow cell clustering, SingleR main annotation, and refined B lineage/plasma compartment annotation.
- **Inputs**: `outputs/rds/00_bonemarrow_qc.rds`.
- **Key steps and parameters**: Main flow as in 01 (also `set.seed(seed)` first; main clustering + SingleR via `run_main_clustering_annotation()`, dims 1:30, resolution 0.6); compartment extraction includes B- and plasma-related labels (regexes per `gate_definitions()`), reclustered at resolution 0.5; uses `annotation_maps/bm_blineage_manual_map.csv`, with the same gate mechanism as 01 (`check_annotation_gate()`, gate id `02_blineage`); marker panel covers pro/pre B → mature B → plasmablast/plasma-like (Vpreb1, Dntt, Cd79a, Cd19, Ms4a1, Ighm, Ighd, Irf4, Prdm1, Xbp1, Sdc1, Tnfrsf17, Jchain, Mki67, Mcl1, Slc3a2).
- **Outputs**: `outputs/rds/02_bonemarrow_annotated.rds`, `02_bm_blineage_refined.rds`; UMAP plots, annotation count tables, marker validation tables.
- **Downstream handoff**: 03/04/05/08 read `02_bonemarrow_annotated.rds`.

### 03_manual_marker_validation.Rmd

- **Purpose**: Manual marker review of the final `CellType_Fine` (pure reporting step; does not modify objects).
- **Inputs**: `outputs/rds/01_spleen_annotated.rds`, `02_bonemarrow_annotated.rds`.
- **Key steps and parameters**: Draws dotplots by fine label with a broader marker panel and generates validation tables; outputs per-mouse cell counts for each label; merges the three manual maps into a review summary table.
- **Outputs**: Spleen/bone marrow dotplots, `03_*_marker_validation_by_fine_label.csv`, `03_*_fine_label_counts_by_mouse.csv`, `03_manual_annotation_review_sheet.csv`.
- **Downstream handoff**: None (for manual review; 07 collects its figures and tables).

### 04_mouse_level_cell_proportion.Rmd

- **Purpose**: Computes cell proportions with mouse as the statistical unit and tests between groups.
- **Inputs**: Annotated RDS from 01/02.
- **Key steps and parameters**: Computes `CellType_Fine` proportions per `MouseID` (missing combinations filled with 0); additionally generates `CellType_Plot` (labels outside the focus list merged into "Other cells"); Wilcoxon rank-sum test per cell type (only tested when each group has ≥2 mice) + BH correction; computed separately for the fine and plot label sets.
- **Outputs**: Proportion and statistics tables for both tissues (one each for fine/plot), stacked bars, boxplots+dots, heatmaps (PDF + PNG).
- **Downstream handoff**: 07 collects its figures and tables.

### 05_pseudobulk_deg_edgeR_DESeq2.Rmd

- **Purpose**: Mouse-level pseudobulk differential expression analysis of target cell types.
- **Inputs**: Annotated RDS from 01/02 (RNA assay counts).
- **Key steps and parameters**: Aggregates raw counts by cell type + `MouseID`; keeps only mice with ≥20 cells; analyzes only when ≥4 pseudobulk samples with both groups present; edgeR QLF (filterByExpr + TMM + robust, design `~Group`) is the default method, with DESeq2 run in parallel when installed; the target cell type list is defined in `config/pipeline_params.R` (currently: 6 spleen, 5 bone marrow); significance called at FDR<0.05 and |logFC|≥0.5.
- **Outputs**: Per-cell-type `05_*_edgeR_pseudobulk_DEG.csv` (+ optional DESeq2 tables), volcano plots, pseudobulk counts RDS, `05_pseudobulk_deg_results.rds`, `05_pseudobulk_deg_summary.csv`.
- **Downstream handoff**: 06 reads the results RDS; 07 consolidates the DEG tables.

### 06_pseudobulk_go_kegg.Rmd

- **Purpose**: GO/KEGG enrichment analysis based on pseudobulk DEGs.
- **Inputs**: `outputs/rds/05_pseudobulk_deg_results.rds` (uses edgeR results only).
- **Key steps and parameters**: Up- and down-regulated DEGs analyzed separately (FDR<0.05 and |logFC|≥0.5); the universe is the genes actually tested in that cell type; GO BP (BH, `qvalueCutoff = 0.2`) and KEGG (mmu); skipped when fewer than 10 genes map to Entrez IDs.
- **Outputs**: Per-direction GO_BP/KEGG CSVs and dotplots, `06_pseudobulk_enrichment_summary.csv`, `06_pseudobulk_enrichment_results.rds`.
- **Downstream handoff**: 07 consolidates the enrichment tables.

### 07_publication_figures_tables.Rmd

- **Purpose**: Consolidates all figures and tables into a publication-ready collection.
- **Inputs**: All products under `outputs/figures/` and `outputs/tables/`.
- **Key steps and parameters**: Copies all PDF/PNG/CSV into `outputs/publication_ready/` (subdirectory names embedded in file names to prevent collisions); merges the edgeR DEG master table and extracts the top 25 by tissue + cell type (by FDR); merges GO/KEGG enrichment tables; generates figure, table, and master manifests recording the provenance of every exported file.
- **Outputs**: `outputs/publication_ready/figures/`, `tables/`, three manifests, merged DEG/enrichment tables.
- **Downstream handoff**: None (terminal step).

### 08_optional_cellchat.Rmd

- **Purpose**: Cell–cell communication analysis (split out of `render_all.R` since 2026-08-06 — CellChat has huge and unstable peak memory, so it is run separately via `run_08_cellchat.sh` and is no longer part of the full render).
- **Inputs**: Annotated RDS from 01/02.
- **Key steps and parameters**: Builds CellChat objects per tissue by group (CellChatDB.mouse, triMean, `min.cells = 10`), merges them and compares interaction counts/strengths with diff network plots; spleen produces T→B bubble plots, bone marrow produces plasma-related target bubble plots; parallelism uses the dedicated `cellchat_workers` (default 2, independent of the global `n_workers` since the 2026-08-06 OOM incident, as peak memory scales with worker count).
- **Outputs**: `08_*_cellchat_merged.rds`, global comparison plots, diff network PDFs, bubble plots.
- **Downstream handoff**: None.

### R/pipeline_helpers.R

- **Purpose**: Shared function library loaded via `source()` by every Rmd 00–08; centralizes common logic for paths, QC, clustering, annotation, plotting, and statistics to keep behavior consistent across scripts.
- **Function inventory (grouped by purpose)**:
  - Paths and environment: `get_pipeline_paths()` (unified management of data, annotation table, RDS, figure, table, report, and gate-evidence directories), `ensure_dir()`, `check_packages()`, `make_sample_manifest()` (derives Group/Tissue/MouseID from sample names).
  - Data loading and metadata: `read_10x_counts()` (prefers the Gene Expression layer), `add_basic_metadata()`, `add_qc_metrics()` (percent.mt / percent.ribo).
  - Clustering and annotation: `run_sctransform_embedding()` (SCTransform + PCA + optional Harmony + clustering + UMAP; SampleID/MouseID forbidden as batch), `run_singler_main()` (SingleR + cluster majority vote), `write_annotation_template()` / `validate_annotation_map()` / `apply_manual_annotation()` / `check_annotation_gate()` (the manual-annotation gate quartet — `check_annotation_gate()` is the unified entry point: writes the template + validates; if not covered, writes `outputs/gate_evidence/gate_status_<gate_id>.csv` and stops with the `GATE[<gate_id>]` prefix; on pass, clears stale status files).
  - Plotting: `theme_publication()`, `save_pub_plot()` (saves PDF + 300 dpi PNG together), `make_marker_dotplot()`, `plot_deg_volcano()`.
  - Statistics and pseudobulk: `normalize_group_levels()`, `calc_mouse_proportions()`, `test_mouse_proportions()` (Wilcoxon + BH), `make_pseudobulk()` (aggregates counts by MouseID), `run_edgeR_pseudobulk()` (QLF, default DEG method), `run_deseq2_pseudobulk()` (optional).
  - Compatibility and utilities: `safe_join_layers()` / `get_assay_matrix()` (Seurat v5 layer compatibility), `marker_validation_table()`, `sanitize_label()`, `%||%`.
  - Legacy functions currently not called by 00–08: `sn_qc_filter()`, `plot_qc_violin()`, `prefix_cells()`, `cbind_common_genes()`, `read_annotation_map()` (00's QC and plotting use script-local inline implementations; modifying or deleting these functions does not affect the pipeline, but this section must be updated accordingly).
- **Call sites**: Every Rmd 00–08 sources this file in its setup chunk; specific call sites are listed under "Key steps and parameters" in each script subsection.

### R/step00_qc_helpers.R

- **Purpose**: Function library dedicated to script 00 (sourced only by 00), holding its heavy logic so that 00 remains a slim pipeline shell.
- **Function inventory (grouped by purpose)**:
  - Input validation: `validate_manifest()` (hard validation of experimental design and matrix integrity; stops on failure).
  - Per-sample processing: `matrix_qc()`, `collapse_qc_reasons()`, `process_one_sample()` (complete per-sample flow: read matrix → light pre-filter → scDblFinder → SoupX → tissue-specific QC).
  - Aggregation and audit: `setup_sample_parallel()` (sample-level parallelism: prefers multicore/fork; falls back to sequential with a warning when fork is unavailable, never to multisession — PSOCK result collection hangs on this WSL2 host), `build_cell_qc_audit()` (per-cell QC audit table), `thresholds_to_data_frame()`, `write_qc_tables()` (all QC tables for this step).
  - QC figures: `plot_qc_metrics()`, `build_qc_figures()` (5 figures: pre/post-correction violins, scatter, doublet-score distribution, rho/doublet bar plots).
  - Merge and validation: `merge_sample_objects()`, `validate_merged_object()`.
- **Call sites**: Sourced only by `00_data_loading_soupx_qc.Rmd` in its setup chunk.

### R/gate_cluster_helpers.R

- **Purpose**: Shared code path for 01/02 clustering/annotation and gate evidence. The 01/02 Rmds and `run_gate_evidence.R` call the same set of functions (with the same `set.seed(seed)` anchor), mechanically guaranteeing that cluster numbering in the evidence matches the official render (replacing the hand-written mirror scripts in agent_scratch from before 2026-08-03, eliminating a structural source of numbering drift).
- **Function inventory (grouped by purpose)**:
  - Gate registry: `gate_definitions()` — the tissue, QC rds, lineage-selection regex, map file, and config variable names for panel/resolution of the three gates (`01_bcell` / `01_tcell` / `02_blineage`); the single definition site for "mirror-sensitive" information such as the regexes.
  - Clustering path: `run_main_clustering_annotation()` (SCTransform embedding + SingleR + `CellType_Fine` initialization), `prepare_lineage_compartment()` (slices compartment by regex + reclusters; stops when the cell count is below `min_cells_for_refinement`).
  - Evidence export: `dump_gate_evidence()` — per-cluster cell counts (by Group), SingleR/main annotation composition, marker panel expression, FindAllMarkers full + top 8, dotplot PNG, written to `outputs/gate_evidence/<gate_id>/`, optionally saving the compartment RDS.
- **Call sites**: Sourced in the setup chunks of 01 and 02; also loaded by `run_gate_evidence.R`.

## Parallel backends and WSL2 performance notes

### Root cause of the hangs (confirmed 2026-08-03)

This environment is WSL2 (kernel `6.18.33.2-microsoft-standard-WSL2`). WSL2's localhost socket readiness-notification mechanism is defective: **PSOCK/socket clusters hang when collecting results of long-running tasks** — the workers have finished computing, but the master process never receives the results. Affected backends include `future::multisession`, `BiocParallel::SnowParam`, and `parallel::makeCluster("PSOCK")`. Tasks that return immediately are unaffected, so the problem is easily misdiagnosed as a package, data, or threading issue. Reproducible with base R alone (independent of pipeline code):

```r
cl <- parallel::makeCluster(2, type = "PSOCK")
parallel::clusterEvalQ(cl, { Sys.sleep(80); "done" })  # never returns
parallel::stopCluster(cl)
```

### Squeezing maximum performance out of the machine under WSL2

Core principle: **use fork-based backends to bypass the socket collection path**.

1. Choose fork-based parallel backends:
   - future: first `options(future.fork.enable = TRUE)`, then `future::plan(future::multicore, workers = N)`
   - BiocParallel: `BiocParallel::MulticoreParam(workers = N)` (SingleR, etc.)
   - base R: `parallel::mclapply(..., mc.cores = N)`
2. Cap threads: at each worker's entry point, run `RhpcBLASctl::blas_set_num_threads(1)` and `RhpcBLASctl::omp_set_num_threads(1)` to avoid the oversubscription livelock of N workers × 28 threads.
3. Exploit fork's copy-on-write: workers inherit the parent's already-loaded packages and large objects, saving memory and avoiding reloads compared with socket clusters, and run faster.
4. Avoid heavy message/print inside workers (socket relay carries a deadlock risk); write progress to log files instead.
5. The local BLAS is OpenBLAS-pthread (r0.3.26), which behaves correctly after fork; if mysterious hangs occur, set `OPENBLAS_NUM_THREADS=1` before starting R as a fallback.

This pipeline is already wrapped accordingly: `parallel_backend = "multicore"` and `n_workers` in `config/pipeline_params.R` control outer parallelism; `setup_sample_parallel()` handles backend selection and fallback (falls back to sequential with a warning when fork is unavailable — typically in RStudio Console/Knit — and never falls back to multisession; use headless Rscript when parallelism is needed under RStudio); 00's outer `future_lapply` schedules per sample dynamically with `future.chunk.size = 1`; `process_one_sample()` caps threads at entry.

### Known residual risks (verified 2026-08-03; assessed and intentionally left unaddressed for now)

1. **scDblFinder 1.24.10 `.xgbtrain`'s `nthreads` is a dead parameter**: `.scDblscore` passes `nthreads = BiocParallel::bpnworkers(BPPARAM)`, but the `.xgbtrain` function body never forwards it to `xgb.cv`/`xgboost` (the parameter name also differs from xgboost's `nthread`), leaving no injection path via the public API. Measured: even after capping both BLAS and OMP at 1 thread, a single `.xgbtrain` training process peaks at 29 threads (xgboost saturates 28 cores; its explicit `num_threads` clause overrides the OMP runtime cap); with 6 workers training concurrently there is hidden thread oversubscription (theoretical peak ~168 threads). Candidate mitigations: set `OMP_THREAD_LIMIT` before render, lower `n_workers`, or inject a patched `.xgbtrain` via `assignInNamespace` (requires a scDblFinder version guard).
2. **00 has no per-sample retry/checkpoint**: if any sample fails, the entire `future_lapply` errors and rendering aborts, discarding the results of completed samples (the RDS is written only after the final merge); the same applies when a fork worker is OOM-killed.
