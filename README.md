# snRNA-seq publication pipeline

## 脚本流程实时总结（00–08）

本章节与 00–08 九个 Rmd 及 `R/` 下的函数库一一对应，内容以当前实际代码为准；修改代码后须同步更新本节。关键参数集中定义于 `config/pipeline_params.R`，各小节列出的参数值对应该文件当前值。

### 00_data_loading_soupx_qc.Rmd

- **作用**：逐样本读取 dnbc4tools 矩阵，完成双联体标记、SoupX 校正和 QC，按组织合并保存。本脚本为精简流程壳；输入校验、单样本处理、审计表与 QC 图构建等重型逻辑封装在 `R/step00_qc_helpers.R`。
- **输入**：`dnbc4tools_results/<SampleID>/<SampleID>/outs/` 下的 `raw_matrix` 与 `filter_matrix`（12 个样本；manifest 硬校验实验设计：2 组 × 2 组织 × 3 重复）。
- **关键步骤与参数**：轻预过滤（≥200 UMI、≥20 genes）仅去空滴；scDblFinder 在 SoupX 之前的原始整数计数上运行（期望率 0.4%/1000 细胞，`dbr.sd = 0`），默认只标记不剔除（`remove_scDblFinder_doublets = FALSE`）；SoupX `autoEstCont` 估计污染率（manifest 可设 `SoupXManualRho` 覆盖；`verbose = FALSE`，避免并行 worker 的 message 经 socket 中继造成死锁），复用 scDblFinder 的 cluster 作为 soup cluster，`adjustCounts` 减法校正并取整；组织特异 QC（脾：100–3000 features、≥500 UMI；骨髓：100–4500 features、≥500 UMI；mt 阈值 2026-08-10 起两组织统一 <0.5%），逐细胞记录失败原因；6 个 workers 并行，后端由 config 的 `parallel_backend` 决定（默认 multicore/fork）；外层 `future_lapply` 以 `future.chunk.size = 1` 逐样本动态调度（避免默认分块 2 样本/worker 的负载不均与失败连坐）。
- **输出**：`outputs/rds/00_spleen_qc.rds`、`00_bonemarrow_qc.rds`；QC summary/阈值/失败原因/SoupX/doublet 等 CSV、QC 图、`00_RUN_STATUS.txt`、sessionInfo。
- **下游衔接**：01 读脾脏 RDS，02 读骨髓 RDS。

### 01_spleen_clustering_annotation.Rmd

- **作用**：脾脏细胞聚类、SingleR 主注释、B/T compartment 精细注释。
- **输入**：`outputs/rds/00_spleen_qc.rds`。
- **关键步骤与参数**：先 `set.seed(seed)` 固定随机种子（保证聚类编号可复现，`run_gate_evidence.R` 依赖此一致性）；主聚类 + SingleR 主注释经共享函数 `run_main_clustering_annotation()` 完成（SCTransform glmGamPoi 回归 percent.mt → PCA → 聚类 dims 1:30 res 0.6 → UMAP；Harmony 默认关闭；SingleR ImmGen 按 cluster 多数票生成 `CellType_Main`）；按 `gate_definitions()` 登记的标签正则提取 B、T compartment 重聚类（`prepare_lineage_compartment()`，dims 1:20，resolution 0.3）；手动注释 gate 由 `check_annotation_gate()` 统一处理（写注释模板 + 校验 map，未覆盖则写 `outputs/gate_evidence/gate_status_01_bcell.csv` / `gate_status_01_tcell.csv` 并以 `GATE[01_bcell]` / `GATE[01_tcell]` 前缀停止）；marker dotplot 与验证表；fine label 回写总对象。
- **输出**：`outputs/rds/01_spleen_annotated.rds`、`01_spleen_bcells_refined.rds`、`01_spleen_tcells_refined.rds`；UMAP 图、注释计数表、marker 验证表。
- **下游衔接**：03/04/05/08 读取 `01_spleen_annotated.rds`。

### 02_bm_clustering_annotation.Rmd

- **作用**：骨髓细胞聚类、SingleR 主注释、B lineage/plasma compartment 精细注释。
- **输入**：`outputs/rds/00_bonemarrow_qc.rds`。
- **关键步骤与参数**：主流程同 01（同样先 `set.seed(seed)`；主聚类 + SingleR 经 `run_main_clustering_annotation()`，dims 1:30，resolution 0.6）；compartment 提取纳入 B 与 plasma 相关标签（正则以 `gate_definitions()` 为准），重聚类 resolution 0.5；使用 `annotation_maps/bm_blineage_manual_map.csv`，gate 机制同 01（`check_annotation_gate()`，gate id `02_blineage`）；marker panel 覆盖 pro/pre B → 成熟 B → 浆母/浆细胞样（Vpreb1、Dntt、Cd79a、Cd19、Ms4a1、Ighm、Ighd、Irf4、Prdm1、Xbp1、Sdc1、Tnfrsf17、Jchain、Mki67、Mcl1、Slc3a2）。
- **输出**：`outputs/rds/02_bonemarrow_annotated.rds`、`02_bm_blineage_refined.rds`；UMAP 图、注释计数表、marker 验证表。
- **下游衔接**：03/04/05/08 读取 `02_bonemarrow_annotated.rds`。

### 03_manual_marker_validation.Rmd

- **作用**：对最终 `CellType_Fine` 做手动 marker 复核（纯报告步骤，不修改对象）。
- **输入**：`outputs/rds/01_spleen_annotated.rds`、`02_bonemarrow_annotated.rds`。
- **关键步骤与参数**：用更宽的 marker panel 按 fine label 绘制 dotplot 并生成验证表；输出每只 mouse 各标签细胞数；合并三张手动 map 生成复核汇总表。
- **输出**：脾/骨髓 dotplot 图、`03_*_marker_validation_by_fine_label.csv`、`03_*_fine_label_counts_by_mouse.csv`、`03_manual_annotation_review_sheet.csv`。
- **下游衔接**：无（供人工复核；07 收集其图表）。

### 04_mouse_level_cell_proportion.Rmd

- **作用**：以 mouse 为统计单位计算细胞比例并做组间检验。
- **输入**：01/02 的 annotated RDS。
- **关键步骤与参数**：按 `MouseID` 计算各 `CellType_Fine` 比例（缺失组合补 0）；另生成 `CellType_Plot`（focus 列表之外合并为 "Other cells"）；每个细胞类型做 Wilcoxon 秩和检验（每组 ≥2 只 mouse 才检验）+ BH 校正；fine 与 plot 两套标签各算一份。
- **输出**：两组织的比例表与统计表（fine/plot 各一份）、stacked bar、boxplot+dots、heatmap（PDF + PNG）。
- **下游衔接**：07 收集其图表。

### 05_pseudobulk_deg_edgeR_DESeq2.Rmd

- **作用**：目标细胞类型的 mouse-level pseudobulk 差异表达分析。
- **输入**：01/02 的 annotated RDS（RNA assay 计数）。
- **关键步骤与参数**：按 cell type + `MouseID` 聚合 raw counts；每只 mouse ≥20 细胞才保留，≥4 个 pseudobulk 样本且两组齐全才分析；edgeR QLF（filterByExpr + TMM + robust，design `~Group`）为默认方法，安装 DESeq2 时平行运行；目标细胞类型清单在 `config/pipeline_params.R` 中定义（当前：脾 6 种、骨髓 5 种）；FDR<0.05 且 |logFC|≥0.5 判定显著。
- **输出**：每细胞类型 `05_*_edgeR_pseudobulk_DEG.csv`（+可选 DESeq2 表）、火山图、pseudobulk counts RDS、`05_pseudobulk_deg_results.rds`、`05_pseudobulk_deg_summary.csv`。
- **下游衔接**：06 读取结果 RDS；07 汇总 DEG 表。

### 06_pseudobulk_go_kegg.Rmd

- **作用**：基于 pseudobulk DEG 的 GO/KEGG 富集分析。
- **输入**：`outputs/rds/05_pseudobulk_deg_results.rds`（只用 edgeR 结果）。
- **关键步骤与参数**：上/下调 DEG 分开做（FDR<0.05 且 |logFC|≥0.5）；universe 取该细胞类型实际被检验的基因；GO BP（BH，`qvalueCutoff = 0.2`）与 KEGG（mmu）；可映射 Entrez 基因 <10 个则跳过。
- **输出**：各方向 GO_BP/KEGG CSV 与 dotplot、`06_pseudobulk_enrichment_summary.csv`、`06_pseudobulk_enrichment_results.rds`。
- **下游衔接**：07 汇总富集表。

### 07_publication_figures_tables.Rmd

- **作用**：汇总全部图表为 publication-ready 集合。
- **输入**：`outputs/figures/`、`outputs/tables/` 下的全部产物。
- **关键步骤与参数**：复制所有 PDF/PNG/CSV 到 `outputs/publication_ready/`（子目录名编入文件名防重名）；合并 edgeR DEG 总表并按组织+细胞类型提取 top 25（按 FDR）；合并 GO/KEGG 富集表；生成图、表和总 manifest 记录每个导出文件的来源。
- **输出**：`outputs/publication_ready/figures/`、`tables/`、三份 manifest、合并 DEG/富集表。
- **下游衔接**：无（终端步骤）。

### 08_optional_cellchat.Rmd

- **作用**：细胞间通讯分析（2026-08-06 起从 `render_all.R` 拆出——CellChat 峰值内存巨大且不稳定，单独用 `run_08_cellchat.sh` 跑，不再是全量的一环）。
- **输入**：01/02 的 annotated RDS。
- **关键步骤与参数**：每组织按组分别构建 CellChat 对象（CellChatDB.mouse，triMean，`min.cells = 10`），merge 后比较互作数量/强度并画 diff network；脾出 T→B bubble plot，骨髓出 plasma 相关靶点 bubble plot；并行用专用的 `cellchat_workers`（默认 2，2026-08-06 OOM 事故后独立于全局 `n_workers`，峰值内存随 worker 数叠加）。
- **输出**：`08_*_cellchat_merged.rds`、全局比较图、diff network PDF、bubble 图。
- **下游衔接**：无。

### R/pipeline_helpers.R

- **作用**：公共函数库，00–08 各 Rmd 均通过 `source()` 加载，集中管理路径、QC、聚类、注释、绘图、统计等公共逻辑，保证各脚本行为一致。
- **函数清单（按用途分组）**：
  - 路径与环境：`get_pipeline_paths()`（统一管理数据、注释表、RDS、图、表、报告与 gate 证据目录）、`ensure_dir()`、`check_packages()`、`make_sample_manifest()`（由样本名推导 Group/Tissue/MouseID）。
  - 数据读取与元数据：`read_10x_counts()`（优先取 Gene Expression 层）、`add_basic_metadata()`、`add_qc_metrics()`（percent.mt / percent.ribo）。
  - 聚类与注释：`run_sctransform_embedding()`（SCTransform + PCA + 可选 Harmony + 聚类 + UMAP，禁止 SampleID/MouseID 作 batch）、`run_singler_main()`（SingleR + cluster 多数票）、`write_annotation_template()` / `validate_annotation_map()` / `apply_manual_annotation()` / `check_annotation_gate()`（手动注释 gate 四件套——`check_annotation_gate()` 为统一入口：写模板 + 校验，未覆盖则写 `outputs/gate_evidence/gate_status_<gate_id>.csv` 并以 `GATE[<gate_id>]` 前缀停止，通过时清除残留状态文件）。
  - 绘图：`theme_publication()`、`save_pub_plot()`（同时存 PDF + 300 dpi PNG）、`make_marker_dotplot()`、`plot_deg_volcano()`。
  - 统计与 pseudobulk：`normalize_group_levels()`、`calc_mouse_proportions()`、`test_mouse_proportions()`（Wilcoxon + BH）、`make_pseudobulk()`（按 MouseID 聚合 counts）、`run_edgeR_pseudobulk()`（QLF，默认 DEG 方法）、`run_deseq2_pseudobulk()`（可选）。
  - 兼容与工具：`safe_join_layers()` / `get_assay_matrix()`（Seurat v5 layer 兼容）、`marker_validation_table()`、`sanitize_label()`、`%||%`。
  - 当前未被 00–08 调用的遗留函数：`sn_qc_filter()`、`plot_qc_violin()`、`prefix_cells()`、`cbind_common_genes()`、`read_annotation_map()`（00 的 QC 与绘图使用脚本内联实现；修改或删除这些函数不影响流程，但须同步本节）。
- **被调用情况**：00–08 全部 Rmd 在 setup chunk 中 `source()` 本文件；具体调用点见各脚本小节的"关键步骤与参数"。

### R/step00_qc_helpers.R

- **作用**：00 脚本专用的函数库（只被 00 `source()`），收纳其重型逻辑，使 00 保持精简流程壳。
- **函数清单（按用途分组）**：
  - 输入校验：`validate_manifest()`（实验设计与矩阵完整性硬校验，不满足即停止）。
  - 单样本处理：`matrix_qc()`、`collapse_qc_reasons()`、`process_one_sample()`（读取矩阵 → 轻预过滤 → scDblFinder → SoupX → 组织特异 QC 的完整单样本流程）。
  - 汇总与审计：`setup_sample_parallel()`（样本级并行：优先 multicore/fork；fork 不可用时回退 sequential 并告警，绝不回退 multisession——PSOCK 结果收集在本 WSL2 主机挂死）、`build_cell_qc_audit()`（逐细胞 QC 审计表）、`thresholds_to_data_frame()`、`write_qc_tables()`（本步骤全部 QC 表）。
  - QC 图：`plot_qc_metrics()`、`build_qc_figures()`（校正前/后小提琴、散点、双联体分数分布、rho/doublet 柱状图共 5 张）。
  - 合并与校验：`merge_sample_objects()`、`validate_merged_object()`。
- **被调用情况**：仅 `00_data_loading_soupx_qc.Rmd` 在 setup chunk 中 `source()`。

### R/gate_cluster_helpers.R

- **作用**：01/02 聚类注释与 gate 证据共用的代码路径。01/02 Rmd 与 `run_gate_evidence.R` 调用同一组函数（加同一 `set.seed(seed)` 锚点），从机制上保证证据的 cluster 编号与正式 render 一致（替代 2026-08-03 之前 agent_scratch 里的手写镜像脚本，消除编号漂移的结构性来源）。
- **函数清单（按用途分组）**：
  - gate 登记：`gate_definitions()`——三个 gate（`01_bcell` / `01_tcell` / `02_blineage`）的组织、QC rds、lineage 筛选正则、map 文件及 panel/resolution 对应的 config 变量名；正则等"镜像敏感"信息的唯一定义处。
  - 聚类路径：`run_main_clustering_annotation()`（SCTransform 嵌入 + SingleR + `CellType_Fine` 初始化）、`prepare_lineage_compartment()`（按正则切 compartment + 重聚类，细胞数低于 `min_cells_for_refinement` 即停止）。
  - 证据导出：`dump_gate_evidence()`——逐 cluster 细胞数（按 Group）、SingleR/主注释组成、marker panel 表达、FindAllMarkers full + top 8、dotplot PNG，输出到 `outputs/gate_evidence/<gate_id>/`，可选保存 compartment RDS。
- **被调用情况**：01、02 的 setup chunk 均 `source()` 本文件；`run_gate_evidence.R` 同样加载。

## 并行后端与 WSL2 性能说明

### 挂死主因（2026-08-03 排查确认）

本环境为 WSL2（内核 `6.18.33.2-microsoft-standard-WSL2`）。WSL2 的 localhost socket 就绪通知机制存在缺陷：**PSOCK/socket 集群对长任务收集结果时会挂死**——worker 已完成计算，主进程却永远收不到结果。受影响的后端包括 `future::multisession`、`BiocParallel::SnowParam`、`parallel::makeCluster("PSOCK")`。瞬时返回的任务不受影响，因此很容易被误判为包、数据或线程问题。用 base R 即可复现（与 pipeline 代码无关）：

```r
cl <- parallel::makeCluster(2, type = "PSOCK")
parallel::clusterEvalQ(cl, { Sys.sleep(80); "done" })  # 永不返回
parallel::stopCluster(cl)
```

### WSL2 下榨干机器性能的做法

核心原则：**用 fork 系后端，绕开 socket 收集路径**。

1. 并行后端选 fork：
   - future：先 `options(future.fork.enable = TRUE)`，再 `future::plan(future::multicore, workers = N)`
   - BiocParallel：`BiocParallel::MulticoreParam(workers = N)`（SingleR 等）
   - base R：`parallel::mclapply(..., mc.cores = N)`
2. 线程封顶：每个 worker 入口执行 `RhpcBLASctl::blas_set_num_threads(1)` 和 `RhpcBLASctl::omp_set_num_threads(1)`，避免 N workers × 28 线程的超订阅 livelock。
3. 利用 fork 的按写复制：worker 继承父进程已加载的包和大对象，比 socket 集群省内存、免重载，速度更快。
4. worker 内避免大量 message/print（socket 中继有死锁风险）；要记录进度就写日志文件。
5. 本机 BLAS 为 OpenBLAS-pthread（r0.3.26），fork 后行为正常；若遇诡异挂死，可在启动 R 前设 `OPENBLAS_NUM_THREADS=1` 兜底。

本 pipeline 已按此封装：`config/pipeline_params.R` 的 `parallel_backend = "multicore"` 与 `n_workers` 控制外层并行，`setup_sample_parallel()` 负责后端选择与回退（fork 不可用——典型如 RStudio Console/Knit——时回退 sequential 并告警，绝不回退 multisession；RStudio 下需要并行请改用 headless Rscript），00 外层 `future_lapply` 以 `future.chunk.size = 1` 逐样本动态调度，`process_one_sample()` 入口做线程封顶。

### 已知遗留风险（2026-08-03 核实，评估后暂不处理）

1. **scDblFinder 1.24.10 `.xgbtrain` 的 `nthreads` 是死参数**：`.scDblscore` 传入 `nthreads = BiocParallel::bpnworkers(BPPARAM)`，但 `.xgbtrain` 函数体从未把它转发给 `xgb.cv`/`xgboost`（形参名亦与 xgboost 的 `nthread` 不一致），经公开 API 无注入路径。实测在 BLAS+OMP 均 cap 1 线程后，单次 `.xgbtrain` 训练进程峰值 29 线程（xgboost 用满 28 核；其显式 `num_threads` 子句会覆盖 OMP 运行时 cap）；6 workers 并发训练时存在隐藏线程超额订阅（理论峰值 ~168 线程）。候选缓解：render 前设 `OMP_THREAD_LIMIT`、降 `n_workers`、或 `assignInNamespace` 注入打补丁版 `.xgbtrain`（需加 scDblFinder 版本守卫）。
2. **00 无逐样本重试/checkpoint**：任一样本失败则整个 `future_lapply` 报错、渲染中止，已完成样本的结果一并丢弃（RDS 仅在末尾合并后写出）；fork worker 被 OOM kill 同理。
