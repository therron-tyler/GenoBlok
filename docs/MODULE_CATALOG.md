# Module catalog

Every block, what it consumes, and what it emits. Blocks in the same track share
the input contracts in [INPUT_SPECS.md](INPUT_SPECS.md), so any block downstream
can pick up where an upstream block left off.

Legend: 🟦 R · 🐍 Python

---

## Bulk RNA-seq track

### `bulk/01_qc/BulkSeq_QC.R` 🟦
Bulk RNA-seq quality-control visualizations.
- **In:** CPM file, raw-counts file (or `none`), groupfile, output dir.
- **Out:** filter grid, per-sample & per-group expression histograms, IQR-outlier
  flags, sample×sample Pearson correlation heatmap, PCA (one per metadata column),
  plus the underlying CSVs (PCA coordinates, correlation drivers, bin data).
- **CLI:** `Rscript BulkSeq_QC.R <cpm> <raw> <groupfile> <outdir>`

### `bulk/02_normalize_cluster/bulk_normalize_cluster_DE.R` 🟦
The bulk workhorse: normalization → gene clustering → differential expression → GSEA.
- **In:** raw counts, CPM, groupfile + a continuous/grouped phenotype column.
- **Does:** CPM filtering → VST → variance-based variable-gene selection →
  hierarchical sample clustering → k-means gene clustering (k = 3–6) with heatmaps →
  DESeq2 (`~ continuous phenotype` and group contrasts) → top-gene box/jitter plots →
  fgsea (Hallmark) dotplots.
- **Out:** gene→cluster maps, per-sample cluster mean z-scores (the input to the
  modeling block), DE tables, heatmaps, GSEA figures.

### `bulk/03_robustness/gene_cluster_robustness.R` 🟦
Are the gene clusters reproducible at a second timepoint / in a second cohort?
- **In:** the V1 gene→cluster map + CPM + groupfile, restricted to paired subjects.
- **Does:** (A) project genes onto fixed V1 centroids; (B) fresh recluster + Hungarian
  label alignment; reports Adjusted Rand Index; Sankey (V1→V2) flow plots; emulated
  silhouette / gene-PCA / per-cluster ORA / cluster-vs-metadata correlation.
- **Out:** robustness metrics, alluvial plots, validation figures.

### `bulk/04_modeling/` 🐍
Phenotype modeling from cluster scores (small-n, leakage-safe).
- `phenotype_regression_nestedCV.py` — nested-CV regression zoo (OLS/Ridge/Lasso/ElasticNet/SVR/kNN)
  predicting a continuous phenotype from per-sample cluster z-scores; VIF, permutation
  importance, LOO predicted-vs-actual, leave-one-subject-out sensitivity.
- `cluster_phenotype_partial_correlation.py` — per-cluster effect on the phenotype adjusted for age + sex
  (partial correlation, Bonferroni).
- `compare_target_transforms.py` — none vs log vs Yeo-Johnson target transforms, scored
  on the original scale.
- **In:** cluster mean-z-score CSV (from block 02) + groupfile. **Out:** results tables + plots.

---

## Single-cell track

### `singlecell/01_build_object/merge_hashed_replicates.R` 🟦
Build a standardized merged Seurat object from hashed per-pool replicates.
- **Does:** add mito/ribo metrics, set HTO order, merge pools, standardize metadata.
- **Out:** merged `.rds` ready for any downstream single-cell block.

### `singlecell/02_normalization_variants/` 🟦
Choose / compare a normalization + integration strategy.
- `compare_normalization_variants.R` — runs **merged / SCTv2 / RPCA / Harmony**
  on identical dims+resolution and emits per-variant UMAP/QC/ADT/clustree/marker figures
  plus a cross-variant comparison panel. `--pca_dims` accepts multiple dim specs
  (e.g. `"1:9,11:12;1:12"`).
- `sctransform_v2_pipeline.R` — focused SCTransform-v2-only path
  (glmGamPoi, regress mito + cell cycle) with PrepSCTFindMarkers DE.

### `singlecell/03_pca_investigation/` 🟦
Understand and trust the embedding before clustering.
- `pca_driver_investigation.R` — which genes drive each PC **before vs after**
  SCT; per-PC driver barplots, biplots, DimHeatmaps, per-sample PC-score violins (catches
  single-donor axes); then a dims×resolution cluster sweep with de novo marker heatmaps.
- `single_sample_embedding_deepdive.R` — single-sample/batch deep-dive: SCT, Harmony,
  integration-anchor heatmap, pseudobulk PCA + fgsea on PC loadings, pseudobulk DESeq2
  (one sample vs rest), feature/marker panels.

### `singlecell/04_subclustering/` 🟦
Subcluster a lineage and run the group comparison.
- `subcluster_group_comparison.R` — subcluster (CD8/Memory-CD4/Naive-CD4),
  compare groups via **propeller** + **MASC (lme4)** proportion testing + pseudobulk DESeq2.
- `subcluster_survey.R` — per-dims survey PDFs: clustree, DE heatmaps,
  canonical violins/feature plots, ADT panels.
- `dims_resolution_scan.R` — scan dims × resolution across an RDS list.
- `adt_umap_survey.R` — ADT/RNA UMAP survey, one page per ADT marker.

### `singlecell/05_figures_report/singlecell_figure_engine.R` 🟦
The single-cell figure engine (~3,300 lines, ~60 plotting functions).
- **Does:** standardized QC, UMAPs (celltype/pain/pool/sample), ADT feature + violin pages,
  stacked cell-composition + propeller/limma proportion tests, pseudobulk DESeq2 per cell
  type, volcano grids, DE heatmaps, gene-overlap Venn/UpSet, GSEA bars, highlight-gene
  box/dot/violin panels — then optionally drops every figure into a PowerPoint template.
- **Out:** `Figures/*.png`, `Tables/*.csv`, `figures_auto_filled.pptx`.

---

## Reporting blocks (track-agnostic)

### `report/decks/` 🐍
- `build_pca_sct_deck.py` — assemble a captioned PowerPoint from a figure folder
  in pipeline order (title → pre-SCT → post-SCT → cluster sweep), with reusable
  text-slide and image-slide helpers.
- `make_deck.py` — generic deck builder that turns labeled figure subfolders into slides.

### `report/apps/fate_mapping_explorer.R` 🟦
An interactive R Shiny app for browsing reporter levels + gene expression across
myeloid fate-mapping models (model selector + tabbed cell-compartment panels). Pattern
to copy when you want to ship an analysis as an interactive explorer rather than static figures.
