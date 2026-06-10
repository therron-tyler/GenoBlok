#!/usr/bin/env Rscript
# =============================================================================
# IMP12_SCTv2_Investigation_Pipeline.R
#
# Mirrors IMP12_Investigation.pptx on the merged WHL1-5 PBMC object.
# Output subfolder: <out_dir>/IMP12_SCTv2_Integration_Analysis/
#
# Stages (comma-separated in --stages, or "all"):
#   qc          - QC violin + scatter by IMP sample and WHL pool
#   sct         - Cell cycle scoring + SCTransform v2 (3000 HVGs) + PCA
#                 + per-sample HVG Jaccard-overlap boxplot (IMP12 vs others)
#   umap_uncorr - Uncorrected (SCT, no Harmony) UMAP + clustering + IMP12 isolation
#   harmony     - Harmony (hash_run) + re-cluster + all Harmony UMAP plots
#   compare     - Facet of raw vs SCT vs SCT+Harmony clustered UMAPs (same params)
#   markers     - PrepSCTFindMarkers + FindAllMarkers + de novo top-5 marker heatmap
#   anchors     - Pairwise normalized integration-anchor heatmap across samples
#   pb_pca      - Pseudobulk PCA all genes per donor + fgsea on PC loadings (Slide 21)
#   deseq2      - Pseudobulk DESeq2 IMP12 vs others pool-tagged replicates + fgsea (Slide 20)
#   features    - FeaturePlot/VlnPlot CX3CR1, EGR1 + canonical marker panels + heatmap (Slide 23)
#
# Example QUEST sbatch call:
#   Rscript IMP12_SCTv2_Investigation_Pipeline.R \
#     --input_rds /path/to/data/20260506_WHL1_SNPrecovery/Merged_nCountRibo150_nCountRNA2000_UMAP30_clust1-9_11-12_res03.rds \
#     --out_dir /path/to/data/20260506_WHL1_SNPrecovery \
#     --pca_dims 20 \
#     --final_res 0.5 \
#     --stages all
# =============================================================================

.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(sctransform)
  library(glmGamPoi)
  library(harmony)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(colorspace)
  library(clustree)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(ggrepel)
  library(Matrix)
  library(optparse)
})

# =============================================================================
## CLI
# =============================================================================
option_list <- list(
  make_option("--input_rds",    type = "character", default = NULL),
  make_option("--out_dir",      type = "character", default = NULL),
  make_option("--focus_sample", type = "character", default = "IMP12",
              help = "IMP sample ID to highlight [default: IMP12]"),
  make_option("--sample_col",   type = "character", default = "HTO_maxID",
              help = "Metadata column for patient IMP IDs [default: HTO_maxID]"),
  make_option("--pool_col",     type = "character", default = "hash_run",
              help = "Metadata column for WHL pool [default: hash_run]"),
  make_option("--n_hvg",        type = "integer",   default = 3000L),
  make_option("--n_pcs",        type = "integer",   default = 30L,
              help = "Number of PCs to compute [default: 30]"),
  make_option("--pca_dims",     type = "integer",   default = 20L,
              help = "PCA dims to use downstream (e.g. 20 → 1:20) [default: 20]"),
  make_option("--final_res",    type = "double",    default = 0.5,
              help = "Clustering resolution to use as final [default: 0.5]"),
  make_option("--gsea_nperm",   type = "integer",   default = 1000L),
  make_option("--gsea_top_n",   type = "integer",   default = 30L),
  make_option("--stages",       type = "character", default = "all",
              help = "Comma-separated stages or 'all' [default: all]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds) || is.null(opt$out_dir))
  stop("--input_rds and --out_dir are required.")

# =============================================================================
## CONFIG
# =============================================================================
cfg <- list(
  input_rds    = opt$input_rds,
  focus_sample = opt$focus_sample,
  sample_col   = opt$sample_col,
  pool_col     = opt$pool_col,

  n_hvg        = opt$n_hvg,
  n_pcs        = opt$n_pcs,
  pca_dims     = seq_len(opt$pca_dims),
  final_res    = opt$final_res,

  regress_vars     = c("percent.mt", "S.Score", "G2M.Score"),
  exclude_patterns = c("^MT-", "^RPS", "^RPL",
                       "^IGHV", "^IGLV", "^IGKV",
                       "^TRAV", "^TRBV", "^TRDV", "^TRGV"),
  resolutions  = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),

  gsea_nperm   = opt$gsea_nperm,
  gsea_top_n   = opt$gsea_top_n,

  umap_w = 10, umap_h = 8,
  split_w = 16
)

ALL_STAGES <- c("qc","sct","umap_uncorr","harmony","compare","markers","anchors",
                "pb_pca","deseq2","features")
active_stages <- if (grepl("^all$", trimws(opt$stages), ignore.case = TRUE)) ALL_STAGES else
  trimws(strsplit(opt$stages, ",")[[1]])
run_stage <- function(s) s %in% active_stages

out_root <- file.path(opt$out_dir, "IMP12_SCTv2_Integration_Analysis")
plot_dir <- file.path(out_root, "plots")
data_dir <- file.path(out_root, "data")
for (d in c(out_root, plot_dir, data_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

ckpt_sct     <- file.path(out_root, "ckpt_post_sct.rds")
ckpt_uncorr  <- file.path(out_root, "ckpt_post_umap_uncorr.rds")
ckpt_harmony <- file.path(out_root, "ckpt_post_harmony.rds")

message("== IMP12 SCTv2 Investigation Pipeline ==")
message("  input_rds    : ", cfg$input_rds)
message("  out_dir      : ", out_root)
message("  focus_sample : ", cfg$focus_sample)
message("  stages       : ", paste(active_stages, collapse = ", "))
message("  n_hvg        : ", cfg$n_hvg, "  pca_dims: 1-", max(cfg$pca_dims),
        "  final_res: ", cfg$final_res)

set.seed(42)

# =============================================================================
## HELPERS
# =============================================================================

save_png <- function(p, name, dir, width, height, dpi = 300) {
  path <- file.path(dir, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, units = "in", dpi = dpi,
         limitsize = FALSE)
  message("  Saved: ", basename(path))
  invisible(path)
}

# DoHeatmap + ggsave frequently yields a blank PNG; render through a real device.
save_png_device <- function(p, name, dir, width, height, dpi = 300) {
  path <- file.path(dir, paste0(name, ".png"))
  png(path, width = width, height = height, units = "in", res = dpi)
  print(p)
  dev.off()
  message("  Saved: ", basename(path))
  invisible(path)
}

# Standard subtitle for every UMAP panel: clustering resolution + PCA dims used.
umap_subtitle <- function(cfg, show_res = TRUE) {
  if (show_res)
    paste0("res = ", cfg$final_res, "  |  dims = 1:", max(cfg$pca_dims))
  else
    paste0("dims = 1:", max(cfg$pca_dims))
}

# Draw a dashed line on a clustree plot marking the final clustering resolution.
# clustree stacks resolutions bottom-up (0-indexed), lowest at the bottom.
add_clustree_res_marker <- function(p, resolutions, final_res) {
  res_sorted <- sort(unique(resolutions))
  idx <- match(final_res, res_sorted)
  if (is.na(idx)) return(p)
  p +
    geom_hline(yintercept = idx - 1, linetype = "dashed",
               colour = "red", linewidth = 0.9) +
    labs(caption = paste0("red dashed line = final resolution (", final_res, ")"))
}

filter_hvg <- function(genes, patterns) {
  keep <- rep(TRUE, length(genes))
  for (pat in patterns) keep <- keep & !grepl(pat, genes)
  genes[keep]
}

safe_counts <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
}

ensure_mito_ribo <- function(obj) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt <- grep("^MT-", rownames(obj), value = TRUE)
    obj[["percent.mt"]] <- if (length(mt) > 0)
      PercentageFeatureSet(obj, pattern = "^MT-") else 0
  }
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    rb <- grep("^RP[SL][[:digit:]]|^RPSA", rownames(obj), value = TRUE)
    obj[["percent.ribo"]] <- if (length(rb) > 0)
      PercentageFeatureSet(obj, features = rb) else 0
  }
  obj
}

ensure_is_focus <- function(obj, cfg) {
  if (!"is_focus" %in% colnames(obj@meta.data))
    obj$is_focus <- ifelse(
      as.character(obj@meta.data[[cfg$sample_col]]) == cfg$focus_sample,
      cfg$focus_sample, "Other")
  obj
}

make_sample_pal <- function(samples, focus) {
  samples <- sort(unique(samples))
  others  <- setdiff(samples, focus)
  # Drop Set1's leading red so no non-focus sample collides with the focus red
  # (previously IMP11 got #E41A1C while IMP12 got #E31A1C — visually identical).
  base_cols <- c(brewer.pal(8, "Set1")[-1], brewer.pal(8, "Set2"))
  pal <- setNames(colorRampPalette(base_cols)(length(others)), others)
  pal <- c(pal, setNames("#E31A1C", focus))   # focus stays pure red
  pal[samples]
}

make_pool_pal <- function(pools) {
  pools <- sort(unique(pools))
  setNames(
    colorRampPalette(brewer.pal(min(length(pools), 8), "Dark2"))(length(pools)),
    pools
  )
}

load_hallmark <- function() {
  msig <- msigdbr(species = "Homo sapiens", category = "H")
  split(msig$gene_symbol, msig$gs_name)
}

gsea_dotplot <- function(gsea_res, title, top_n) {
  n <- min(top_n %/% 2, 50)
  top_p <- head(gsea_res[!is.na(gsea_res$NES) & gsea_res$NES > 0, ], n)
  top_n2 <- tail(gsea_res[!is.na(gsea_res$NES) & gsea_res$NES < 0, ], n)
  df <- rbind(top_p, top_n2)
  if (nrow(df) == 0) return(NULL)
  df$pathway <- sub("HALLMARK_", "", df$pathway)
  df$pathway <- factor(df$pathway, levels = df$pathway[order(df$NES)])
  ggplot(df, aes(x = NES, y = pathway, color = padj, size = size)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_gradient(low = "#B2182B", high = "grey80", name = "padj",
                         na.value = "grey80") +
    theme_cowplot(11) +
    labs(title = title, x = "NES", y = NULL)
}

# =============================================================================
## STAGE: qc
# =============================================================================
stage_qc <- function(obj, cfg, plot_dir) {
  message("\n[qc] QC plots by IMP sample and WHL pool...")

  qc_feats <- intersect(
    c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
    colnames(obj@meta.data)
  )
  n_feat <- length(qc_feats)
  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  pools    <- sort(unique(as.character(obj@meta.data[[cfg$pool_col]])))
  samp_pal <- make_sample_pal(samples, cfg$focus_sample)
  pool_pal <- make_pool_pal(pools)

  # Violin by IMP sample
  obj@meta.data[[cfg$sample_col]] <- factor(obj@meta.data[[cfg$sample_col]],
                                             levels = samples)
  plist <- VlnPlot(obj, features = qc_feats, group.by = cfg$sample_col,
                   pt.size = 0, cols = samp_pal, combine = FALSE)
  for (i in seq_along(plist))
    plist[[i]] <- plist[[i]] + NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.text.x = if (i < n_feat) element_blank() else
              element_text(angle = 45, hjust = 1, size = 9))
  p <- wrap_plots(plist, ncol = 1) + plot_annotation(title = "QC by IMP Sample")
  save_png(p, "00a_qc_violin_by_IMP_sample", plot_dir,
           width = max(14, length(samples) * 0.5), height = 3.5 * n_feat)

  # Violin by WHL pool
  obj@meta.data[[cfg$pool_col]] <- factor(obj@meta.data[[cfg$pool_col]], levels = pools)
  plist2 <- VlnPlot(obj, features = qc_feats, group.by = cfg$pool_col,
                    pt.size = 0, cols = pool_pal, combine = FALSE)
  for (i in seq_along(plist2))
    plist2[[i]] <- plist2[[i]] + NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.text.x = if (i < n_feat) element_blank() else
              element_text(angle = 45, hjust = 1))
  p2 <- wrap_plots(plist2, ncol = 1) + plot_annotation(title = "QC by WHL Pool")
  save_png(p2, "00b_qc_violin_by_pool", plot_dir, width = 9, height = 3.5 * n_feat)

  # Scatter nCount_RNA vs nCount_Ribo
  if (all(c("nCount_RNA","nCount_Ribo") %in% colnames(obj@meta.data))) {
    md <- obj@meta.data[, c("nCount_RNA","nCount_Ribo",
                             cfg$pool_col, cfg$sample_col), drop = FALSE]
    md$is_focus <- ifelse(as.character(md[[cfg$sample_col]]) == cfg$focus_sample,
                           cfg$focus_sample, "Other")
    p3 <- ggplot(md, aes(x = nCount_RNA, y = nCount_Ribo,
                          color = .data[[cfg$pool_col]])) +
      geom_point(size = 0.15, alpha = 0.2) +
      scale_color_manual(values = pool_pal) +
      theme_cowplot(12) +
      labs(title = "nCount_RNA vs nCount_Ribo by WHL Pool", color = "Pool")
    save_png(p3, "00c_scatter_nCountRNA_vs_nCountRibo", plot_dir, 8, 6)
  }

  invisible(NULL)
}

# Per-sample HVGs (pools merged) → pairwise Jaccard overlap, split by whether the
# pair contains the focus sample. Reproduces the "HVG overlap by pair type" boxplot.
plot_hvg_overlap <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[hvg] Per-sample HVG Jaccard overlap (", cfg$focus_sample, " vs others)...")
  DefaultAssay(obj) <- "RNA"
  samples <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))

  hvg_list <- list()
  for (s in samples) {
    cells <- rownames(obj@meta.data)[
      as.character(obj@meta.data[[cfg$sample_col]]) == s]
    if (length(cells) < 20) { message("  ", s, ": <20 cells, skipped"); next }
    sub <- subset(obj, cells = cells)
    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- FindVariableFeatures(sub, selection.method = "vst",
                                nfeatures = cfg$n_hvg, verbose = FALSE)
    hvg_list[[s]] <- filter_hvg(VariableFeatures(sub), cfg$exclude_patterns)
  }

  samples <- names(hvg_list)
  if (length(samples) < 3 || !(cfg$focus_sample %in% samples)) {
    message("  Too few samples (or focus missing) — skipping HVG-overlap boxplot.")
    return(invisible(NULL))
  }

  pairs <- combn(samples, 2, simplify = FALSE)
  jac   <- vapply(pairs, function(p) {
    a <- hvg_list[[p[1]]]; b <- hvg_list[[p[2]]]
    length(intersect(a, b)) / length(union(a, b))
  }, numeric(1))
  has_focus <- vapply(pairs, function(p) cfg$focus_sample %in% p, logical(1))

  lvl_other <- paste0("non-", cfg$focus_sample, " pairs")
  lvl_focus <- paste0(cfg$focus_sample, " pairs")
  df <- data.frame(
    sample1   = vapply(pairs, `[`, character(1), 1),
    sample2   = vapply(pairs, `[`, character(1), 2),
    jaccard   = jac,
    pair_type = factor(ifelse(has_focus, lvl_focus, lvl_other),
                       levels = c(lvl_other, lvl_focus)),
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(data_dir, "hvg_overlap_by_pair_type.csv"), row.names = FALSE)

  pal <- setNames(c("grey75", "#E07B6B"), c(lvl_other, lvl_focus))
  p <- ggplot(df, aes(x = pair_type, y = jaccard, fill = pair_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "HVG overlap by pair type", x = NULL, y = "Jaccard overlap") +
    theme_cowplot(14)
  save_png(p, "01b_hvg_overlap_by_pair_type", plot_dir, 7, 7)
  invisible(df)
}

# =============================================================================
## STAGE: sct — cell cycle + SCTransform v2 + PCA
# =============================================================================
stage_sct <- function(obj, cfg, plot_dir, data_dir, ckpt) {
  message("\n[sct] Cell cycle scoring + SCTransform v2 + PCA...")
  DefaultAssay(obj) <- "RNA"

  # Join Seurat v5 layers if needed
  tryCatch(obj <- JoinLayers(obj, assay = "RNA"), error = function(e) NULL)

  obj <- ensure_mito_ribo(obj)

  # Quick log-norm for cell cycle scoring only; SCT re-normalizes from raw counts
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(
    obj,
    s.features   = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident    = FALSE, verbose = FALSE
  )
  message("  Cell cycle: ", paste(names(table(obj$Phase)), table(obj$Phase),
                                   sep = "=", collapse = " | "))

  # SCTransform v2
  message("  SCTransform v2 (n_hvg=", cfg$n_hvg, ", regress: ",
          paste(cfg$regress_vars, collapse = ", "), ")...")
  obj <- SCTransform(
    obj,
    vst.flavor          = "v2",
    method              = "glmGamPoi",
    vars.to.regress     = cfg$regress_vars,
    variable.features.n = cfg$n_hvg,
    verbose             = TRUE
  )

  all_hvg     <- VariableFeatures(obj)
  cleaned_hvg <- filter_hvg(all_hvg, cfg$exclude_patterns)
  VariableFeatures(obj) <- cleaned_hvg
  message("  HVGs: ", length(all_hvg), " → ", length(cleaned_hvg),
          " after stripping MT/ribo/TCR/Ig")
  write.csv(data.frame(gene = cleaned_hvg),
            file.path(data_dir, "sct_hvg_list.csv"), row.names = FALSE)

  # PCA
  message("  RunPCA (n=", cfg$n_pcs, ")...")
  obj <- RunPCA(obj, features = cleaned_hvg, npcs = cfg$n_pcs, verbose = FALSE)

  p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs) +
    ggtitle("PCA Elbow Plot (SCTransform v2)") + theme_cowplot(14)
  save_png(p_elbow, "01_elbow_pca", plot_dir, 8, 5)

  # Per-sample HVG Jaccard overlap boxplot (IMP12 vs others)
  tryCatch(plot_hvg_overlap(obj, cfg, plot_dir, data_dir),
           error = function(e) message("  HVG-overlap boxplot failed: ",
                                        conditionMessage(e)))

  saveRDS(obj, ckpt); message("  Checkpoint: ", basename(ckpt))
  obj
}

# =============================================================================
## STAGE: umap_uncorr
# =============================================================================
stage_umap_uncorr <- function(obj, cfg, plot_dir, data_dir, ckpt) {
  message("\n[umap_uncorr] Uncorrected UMAP + clustering + IMP12 split plots...")

  obj <- ensure_is_focus(obj, cfg)
  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  samp_pal <- make_sample_pal(samples, cfg$focus_sample)
  pool_pal <- make_pool_pal(unique(as.character(obj@meta.data[[cfg$pool_col]])))

  obj <- RunUMAP(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)

  for (res in cfg$resolutions) {
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    message("  res=", res, " → ", length(unique(Idents(obj))), " clusters")
  }

  p_ct <- clustree(obj@meta.data, prefix = "SCT_snn_res.") +
    ggtitle("Clustree — uncorrected SCTv2") + theme(legend.position = "right")
  p_ct <- add_clustree_res_marker(p_ct, cfg$resolutions, cfg$final_res)
  save_png(p_ct, "02_clustree_uncorr", plot_dir, 12, 14)

  Idents(obj) <- paste0("SCT_snn_res.", cfg$final_res)
  obj$seurat_clusters <- Idents(obj)
  message("  Active res=", cfg$final_res, " → ", length(unique(Idents(obj))), " clusters")

  # Slide 6: clusters
  p_clust <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
                     label.size = 5) + NoLegend() +
    labs(title = "Uncorrected clusters (SCT, no Harmony)",
         subtitle = umap_subtitle(cfg)) + theme_cowplot(13)
  save_png(p_clust, "03a_umap_uncorr_clusters", plot_dir, cfg$umap_w, cfg$umap_h)

  # Slide 5: by WHL pool
  p_pool <- DimPlot(obj, reduction = "umap", group.by = cfg$pool_col,
                    cols = pool_pal, pt.size = 0.3) +
    labs(title = "Samples by WHL Pool (uncorrected)",
         subtitle = umap_subtitle(cfg, show_res = FALSE)) + theme_cowplot(13)
  save_png(p_pool, "03b_umap_uncorr_by_pool", plot_dir, cfg$umap_w + 2, cfg$umap_h)

  # by IMP sample
  p_imp <- DimPlot(obj, reduction = "umap", group.by = cfg$sample_col,
                   cols = samp_pal, pt.size = 0.3) +
    labs(title = "Samples by IMP ID (uncorrected)",
         subtitle = umap_subtitle(cfg, show_res = FALSE)) + theme_cowplot(13)
  save_png(p_imp, "03c_umap_uncorr_by_IMP_sample", plot_dir, cfg$umap_w + 4, cfg$umap_h)

  # Slide 7: IMP12 (right) vs others (left)
  focus_col <- setNames("#E31A1C", cfg$focus_sample)
  p_split <- DimPlot(obj, reduction = "umap", group.by = cfg$sample_col,
                     cells.highlight = rownames(obj@meta.data)[
                       obj@meta.data[[cfg$sample_col]] == cfg$focus_sample],
                     cols.highlight = "#E31A1C", cols = "grey80", pt.size = 0.3,
                     split.by = "is_focus") +
    labs(title = paste0(cfg$focus_sample, " vs All Others — uncorrected UMAP"),
         subtitle = umap_subtitle(cfg, show_res = FALSE)) +
    theme_cowplot(12)
  save_png(p_split, "04a_umap_IMP12_vs_others_split", plot_dir, cfg$split_w, cfg$umap_h)

  # Slide 8: identify clusters >80% IMP12
  clust_vec <- as.character(obj$seurat_clusters)
  samp_vec  <- as.character(obj@meta.data[[cfg$sample_col]])
  keep      <- !is.na(clust_vec) & !is.na(samp_vec)
  clust_tab <- as.data.frame(
    table(seurat_clusters = clust_vec[keep], sample_id = samp_vec[keep], useNA = "no"),
    stringsAsFactors = FALSE
  )
  colnames(clust_tab) <- c("seurat_clusters", cfg$sample_col, "n")
  clust_comp <- clust_tab %>%
    group_by(seurat_clusters) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()
  focus_comp <- clust_comp %>%
    filter(.data[[cfg$sample_col]] == cfg$focus_sample) %>%
    arrange(desc(pct))
  write.csv(focus_comp,
            file.path(data_dir, "cluster_IMP12_pct_uncorr.csv"), row.names = FALSE)

  unique_clusts <- focus_comp$seurat_clusters[focus_comp$pct > 0.8]
  message("  Clusters >80% ", cfg$focus_sample, ": ",
          if (length(unique_clusts) == 0) "none" else paste(unique_clusts, collapse = ", "))

  if (length(unique_clusts) > 0) {
    obj$cluster_type <- ifelse(obj$seurat_clusters %in% unique_clusts,
                                paste0(cfg$focus_sample, "-unique"), "Shared")
    p_uniq <- DimPlot(obj, reduction = "umap", group.by = "cluster_type",
                      cols = c(setNames("#E31A1C", paste0(cfg$focus_sample, "-unique")),
                               "Shared" = "grey80"), pt.size = 0.3) +
      labs(title = paste0("Clusters >80% ", cfg$focus_sample, " cells"),
           subtitle = umap_subtitle(cfg)) + theme_cowplot(13)
    save_png(p_uniq, "04b_umap_IMP12_unique_clusters", plot_dir, cfg$umap_w, cfg$umap_h)
  }

  # Slide 9: IMP12 subset only — its own UMAP
  focus_cells <- rownames(obj@meta.data)[
    as.character(obj@meta.data[[cfg$sample_col]]) == cfg$focus_sample]
  obj_focus   <- subset(obj, cells = focus_cells)
  obj_focus   <- RunUMAP(obj_focus, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)
  p_foc <- DimPlot(obj_focus, reduction = "umap", label = TRUE, repel = TRUE,
                   label.size = 5) + NoLegend() +
    labs(title = paste0(cfg$focus_sample, " alone"),
         subtitle = umap_subtitle(cfg)) + theme_cowplot(13)
  save_png(p_foc, "04c_umap_IMP12_subset_only", plot_dir, cfg$umap_w, cfg$umap_h)
  rm(obj_focus)

  saveRDS(obj, ckpt); message("  Checkpoint: ", basename(ckpt))
  obj
}

# =============================================================================
## STAGE: harmony
# =============================================================================
stage_harmony <- function(obj, cfg, plot_dir, data_dir, ckpt) {
  message("\n[harmony] RunHarmony (", cfg$pool_col, ") + re-cluster + UMAP plots...")

  obj <- ensure_is_focus(obj, cfg)
  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  samp_pal <- make_sample_pal(samples, cfg$focus_sample)
  pool_pal <- make_pool_pal(unique(as.character(obj@meta.data[[cfg$pool_col]])))

  obj$clusters_uncorr <- obj$seurat_clusters

  # Harmony on PCA — correct only for WHL pool (technical batch)
  obj <- RunHarmony(
    obj,
    group.by.vars  = cfg$pool_col,
    reduction      = "pca",
    dims.use       = cfg$pca_dims,
    verbose        = TRUE,
    reduction.save = "harmony"
  )

  obj <- RunUMAP(obj, reduction = "harmony", dims = cfg$pca_dims,
                 reduction.name = "umap_harmony", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = cfg$pca_dims,
                       graph.name = "SCT_harmony_snn", verbose = FALSE)

  for (res in cfg$resolutions) {
    obj <- FindClusters(obj, graph.name = "SCT_harmony_snn",
                         resolution = res, verbose = FALSE)
    message("  res=", res, " → ", length(unique(Idents(obj))), " clusters")
  }

  p_ct <- clustree(obj@meta.data, prefix = "SCT_harmony_snn_res.") +
    ggtitle("Clustree — Harmony-corrected") + theme(legend.position = "right")
  p_ct <- add_clustree_res_marker(p_ct, cfg$resolutions, cfg$final_res)
  save_png(p_ct, "05_clustree_harmony", plot_dir, 12, 14)

  Idents(obj) <- paste0("SCT_harmony_snn_res.", cfg$final_res)
  obj$seurat_clusters_harmony <- Idents(obj)
  n_h <- length(unique(Idents(obj)))
  message("  Harmony clusters: ", n_h, " (res=", cfg$final_res, ")")

  # ── Harmony UMAP visualizations ──────────────────────────────────────────────

  # Slide 11: corrected for pool — color by pool
  p_pool <- DimPlot(obj, reduction = "umap_harmony", group.by = cfg$pool_col,
                    cols = pool_pal, pt.size = 0.3) +
    labs(title = "Harmony UMAP — by WHL Pool (corrected)",
         subtitle = umap_subtitle(cfg, show_res = FALSE)) + theme_cowplot(13)
  save_png(p_pool, "06a_harmony_umap_by_pool", plot_dir, cfg$umap_w + 2, cfg$umap_h)

  # Slide 12: new Harmony clusters
  p_nh <- DimPlot(obj, reduction = "umap_harmony", label = TRUE,
                   repel = TRUE, label.size = 5) + NoLegend() +
    labs(title = "Harmony UMAP — new clusters",
         subtitle = umap_subtitle(cfg)) +
    theme_cowplot(13)
  save_png(p_nh, "06b_harmony_umap_new_clusters", plot_dir, cfg$umap_w, cfg$umap_h)

  # Slide 13: original (uncorrected) clusters on Harmony UMAP
  p_oh <- DimPlot(obj, reduction = "umap_harmony", group.by = "clusters_uncorr",
                   label = TRUE, repel = TRUE, label.size = 5) + NoLegend() +
    labs(title = "Harmony UMAP — original uncorrected clusters",
         subtitle = umap_subtitle(cfg)) + theme_cowplot(13)
  save_png(p_oh, "06c_harmony_umap_orig_clusters", plot_dir, cfg$umap_w, cfg$umap_h)

  # by IMP sample
  p_imp <- DimPlot(obj, reduction = "umap_harmony", group.by = cfg$sample_col,
                   cols = samp_pal, pt.size = 0.3) +
    labs(title = "Harmony UMAP — by IMP Sample",
         subtitle = umap_subtitle(cfg, show_res = FALSE)) + theme_cowplot(13)
  save_png(p_imp, "06d_harmony_umap_by_IMP_sample", plot_dir, cfg$umap_w + 4, cfg$umap_h)

  # Slide 14: IMP12 vs others, old clusters
  p_old_split <- DimPlot(obj, reduction = "umap_harmony",
                          group.by = "clusters_uncorr",
                          split.by = "is_focus",
                          label = TRUE, repel = TRUE, label.size = 4) + NoLegend() +
    labs(title = paste0("Harmony UMAP — original clusters | ",
                        cfg$focus_sample, " vs Others"),
         subtitle = umap_subtitle(cfg)) + theme_cowplot(12)
  save_png(p_old_split, "07a_harmony_IMP12_vs_others_old_clusters",
           plot_dir, cfg$split_w, cfg$umap_h)

  # Slide 15: IMP12 vs others, new clusters
  p_new_split <- DimPlot(obj, reduction = "umap_harmony",
                          group.by = "seurat_clusters_harmony",
                          split.by = "is_focus",
                          label = TRUE, repel = TRUE, label.size = 4) + NoLegend() +
    labs(title = paste0("Harmony UMAP — new clusters | ",
                        cfg$focus_sample, " vs Others"),
         subtitle = umap_subtitle(cfg)) + theme_cowplot(12)
  save_png(p_new_split, "07b_harmony_IMP12_vs_others_new_clusters",
           plot_dir, cfg$split_w, cfg$umap_h)

  # ── Slide 16: cluster composition stacked bar + identify missing clusters ─────
  clust_vec_h <- as.character(obj$seurat_clusters_harmony)
  samp_vec_h  <- as.character(obj@meta.data[[cfg$sample_col]])
  keep_h      <- !is.na(clust_vec_h) & !is.na(samp_vec_h)
  prop_tab    <- as.data.frame(
    table(seurat_clusters_harmony = clust_vec_h[keep_h],
          sample_id = samp_vec_h[keep_h], useNA = "no"),
    stringsAsFactors = FALSE
  )
  colnames(prop_tab) <- c("seurat_clusters_harmony", cfg$sample_col, "n")
  prop_tab$n <- as.integer(prop_tab$n)
  prop_df <- prop_tab %>%
    group_by(.data[[cfg$sample_col]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  clust_lvls <- levels(factor(as.integer(as.character(
    prop_df$seurat_clusters_harmony))))
  prop_df$seurat_clusters_harmony <- factor(prop_df$seurat_clusters_harmony,
                                             levels = clust_lvls)
  n_clust <- length(clust_lvls)
  clust_pal <- setNames(
    colorRampPalette(c(brewer.pal(8,"Set2"), brewer.pal(8,"Set1")))(n_clust),
    clust_lvls
  )
  samp_order <- c(cfg$focus_sample,
                  sort(setdiff(unique(as.character(prop_df[[cfg$sample_col]])),
                               cfg$focus_sample)))
  prop_df[[cfg$sample_col]] <- factor(prop_df[[cfg$sample_col]], levels = samp_order)

  p_bar <- ggplot(prop_df,
                  aes(x = .data[[cfg$sample_col]], y = prop,
                      fill = seurat_clusters_harmony)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = clust_pal, name = "Cluster") +
    labs(title = "Cluster Composition by IMP Sample (Harmony)",
         x = NULL, y = "Proportion") +
    theme_cowplot(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  save_png(p_bar, "08_cluster_composition_by_IMP_sample",
           plot_dir, width = max(18, length(samp_order) * 0.6), height = 7)

  write.csv(prop_df, file.path(data_dir, "cluster_composition_harmony.csv"),
            row.names = FALSE)

  saveRDS(obj, ckpt)
  message("  Checkpoint: ", basename(ckpt))
  obj
}

# =============================================================================
## STAGE: markers
# =============================================================================
stage_markers <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[markers] PrepSCTFindMarkers + FindAllMarkers + de novo top-5 heatmap...")
  DefaultAssay(obj) <- "SCT"
  Idents(obj) <- "seurat_clusters_harmony"

  obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)

  markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE,
                             min.pct = 0.05, logfc.threshold = 0.2,
                             test.use = "wilcox", verbose = FALSE)
  markers_sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, ]
  write.csv(markers_sig, file.path(data_dir, "harmony_cluster_markers_sig.csv"),
            row.names = FALSE)
  write.csv(markers,     file.path(data_dir, "harmony_cluster_markers_all.csv"),
            row.names = FALSE)
  message("  Significant markers: ", nrow(markers_sig))

  # Slide 19: cluster IMP-sample composition table (who contributes to each cluster)
  clust_vec_m <- as.character(obj$seurat_clusters_harmony)
  samp_vec_m  <- as.character(obj@meta.data[[cfg$sample_col]])
  keep_m      <- !is.na(clust_vec_m) & !is.na(samp_vec_m)
  imp_tab     <- as.data.frame(
    table(seurat_clusters_harmony = clust_vec_m[keep_m],
          sample_id = samp_vec_m[keep_m], useNA = "no"),
    stringsAsFactors = FALSE
  )
  colnames(imp_tab) <- c("seurat_clusters_harmony", cfg$sample_col, "n")
  imp_tab$n <- as.integer(imp_tab$n)
  clust_imp <- imp_tab %>%
    group_by(seurat_clusters_harmony) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    arrange(seurat_clusters_harmony, desc(pct))
  write.csv(clust_imp, file.path(data_dir, "cluster_per_IMP_sample_pct.csv"),
            row.names = FALSE)

  if (nrow(markers_sig) == 0) {
    message("  No significant markers found — skipping de novo heatmap.")
    return(invisible(NULL))
  }

  # ── De novo top-5 markers per cluster heatmap (replaces the small dotplots) ───
  # Large, italic gene-name rows; rendered via png device (DoHeatmap + ggsave
  # frequently produces a blank PNG, and SCT scale.data must hold these genes).
  top5 <- markers_sig %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>% pull(gene) %>% unique()
  if (length(top5) > 0) {
    cl_levels <- levels(factor(suppressWarnings(
      as.integer(as.character(obj$seurat_clusters_harmony)))))
    if (anyNA(cl_levels) || length(cl_levels) == 0)
      cl_levels <- sort(unique(as.character(obj$seurat_clusters_harmony)))
    Idents(obj) <- factor(as.character(obj$seurat_clusters_harmony), levels = cl_levels)

    obj <- ScaleData(obj, features = top5, assay = "SCT", verbose = FALSE)
    n_cl <- length(levels(Idents(obj)))

    p_heat <- DoHeatmap(obj, features = top5, assay = "SCT", slot = "scale.data",
                        size = 5, angle = 45) +
      scale_fill_gradientn(colors = c("#2166AC", "white", "#B2182B"),
                           name = "Scaled\nexpression") +
      theme(axis.text.y = element_text(size = 13, face = "italic")) +
      ggtitle("De novo top 5 markers per Harmony cluster (avg_log2FC)")
    save_png_device(p_heat, "10_denovo_top5_marker_heatmap", plot_dir,
                    width  = max(12, n_cl * 0.9),
                    height = min(40, max(10, length(top5) * 0.34)))
  }

  invisible(markers_sig)
}

# =============================================================================
## STAGE: pb_pca — Pseudobulk PCA all genes per donor + fgsea (Slide 21)
# =============================================================================
stage_pb_pca <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[pb_pca] Pseudobulk PCA (all genes, per-donor) + GSEA on PC loadings...")

  counts <- safe_counts(obj, "RNA")
  md     <- obj@meta.data

  # Aggregate per HTO_maxID (all pools merged → one point per patient)
  donors <- sort(unique(as.character(md[[cfg$sample_col]])))
  pb_mat <- vapply(donors, function(d) {
    cells <- rownames(md)[as.character(md[[cfg$sample_col]]) == d]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(pb_mat) <- rownames(counts)
  colnames(pb_mat) <- donors

  # Filter: gene expressed (count > 0) in ≥20% of donors
  min_d  <- max(2, round(0.2 * length(donors)))
  keep   <- rowSums(pb_mat > 0) >= min_d
  pb_mat <- pb_mat[keep, ]
  message("  Genes after 20%-donor filter: ", nrow(pb_mat))

  # CPM + log2
  lib_sz  <- colSums(pb_mat)
  cpm_mat <- sweep(pb_mat, 2, lib_sz / 1e6, "/")
  log_cpm <- log2(cpm_mat + 1)

  # Export for PROMOTER pipeline compatibility (tab-delimited, genes as rows)
  write.table(pb_mat,  file.path(data_dir, "pseudobulk_perdonor_raw_counts.txt"),
              sep = "\t", quote = FALSE)
  write.table(log_cpm, file.path(data_dir, "pseudobulk_perdonor_logCPM.txt"),
              sep = "\t", quote = FALSE)

  # PCA
  pca_res <- prcomp(t(log_cpm), center = TRUE, scale. = TRUE)
  scores  <- as.data.frame(pca_res$x)
  var_pct <- round(summary(pca_res)$importance[2, 1:min(10, ncol(scores))] * 100, 1)

  # Biplot dataframe
  scores$sample   <- rownames(scores)
  scores$is_focus <- ifelse(scores$sample == cfg$focus_sample,
                             cfg$focus_sample, "Other")
  focus_col_pal   <- c(setNames("#E31A1C", cfg$focus_sample), "Other" = "#6BAED6")

  p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = is_focus, label = sample)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3.5, max.overlaps = 40) +
    scale_color_manual(values = focus_col_pal, name = NULL) +
    labs(title = "Pseudobulk PCA — All Genes (per donor)",
         x = paste0("PC1 (", var_pct[1], "%)"),
         y = paste0("PC2 (", var_pct[2], "%)")) +
    theme_cowplot(13)
  save_png(p_pca, "12a_pseudobulk_PCA_biplot", plot_dir, 9, 7)

  if (ncol(scores) >= 4) {
    p_pca13 <- ggplot(scores, aes(x = PC1, y = PC3, color = is_focus, label = sample)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3.5, max.overlaps = 40) +
      scale_color_manual(values = focus_col_pal, name = NULL) +
      labs(title = "Pseudobulk PCA — PC1 vs PC3",
           x = paste0("PC1 (", var_pct[1], "%)"),
           y = paste0("PC3 (", var_pct[3], "%)")) +
      theme_cowplot(13)
    save_png(p_pca13, "12b_pseudobulk_PCA_PC1_vs_PC3", plot_dir, 9, 7)
  }

  # Scree
  scree_df <- data.frame(PC = seq_along(var_pct), var = var_pct)
  p_scree <- ggplot(scree_df, aes(x = PC, y = var)) +
    geom_line() + geom_point(size = 2) + theme_cowplot(13) +
    labs(title = "Pseudobulk PCA Scree", x = "PC", y = "% Variance Explained")
  save_png(p_scree, "12c_pseudobulk_PCA_scree", plot_dir, 7, 5)

  # ── GSEA on combined PC1+PC2 loadings in direction of IMP12 separation ────────
  focus_row <- rownames(pca_res$x) == cfg$focus_sample
  if (!any(focus_row)) {
    message("  WARNING: ", cfg$focus_sample,
            " not in pseudobulk donors — skipping GSEA.")
    return(invisible(NULL))
  }

  other_mean <- colMeans(pca_res$x[!focus_row, 1:2, drop = FALSE])
  focus_score <- pca_res$x[focus_row, 1:2]

  # Sign convention (matches PPTX "Ranking relative to IMP12 (-PC1, +PC2)")
  sign1 <- sign(as.numeric(focus_score[1]) - other_mean[1])
  sign2 <- sign(as.numeric(focus_score[2]) - other_mean[2])
  message("  ", cfg$focus_sample, " PC1=", round(focus_score[1], 2),
          " sign=", sign1, "  PC2=", round(focus_score[2], 2), " sign=", sign2)

  loadings     <- pca_res$rotation[, 1:2]
  combined     <- sign1 * loadings[, 1] + sign2 * loadings[, 2]
  gene_ranks   <- sort(combined, decreasing = TRUE)

  write.csv(data.frame(gene = names(gene_ranks), loading_score = gene_ranks),
            file.path(data_dir, "pseudobulk_PCA_combined_loading_ranks.csv"),
            row.names = FALSE)

  pathways_h <- load_hallmark()
  gsea_pca   <- fgsea(pathways = pathways_h, stats = gene_ranks,
                       nPermSimple = cfg$gsea_nperm, minSize = 15, maxSize = 500)
  gsea_pca   <- gsea_pca[order(gsea_pca$NES, decreasing = TRUE), ]
  write.csv(gsea_pca[, c("pathway","pval","padj","NES","size")],
            file.path(data_dir, "pseudobulk_PCA_GSEA_hallmark.csv"), row.names = FALSE)

  p_g <- gsea_dotplot(gsea_pca,
                       paste0("GSEA on Pseudobulk PC Loadings\n(",
                              cfg$focus_sample, " vs others direction)"),
                       cfg$gsea_top_n)
  if (!is.null(p_g)) {
    # Make the IMP12-vs-others direction explicit: ranks point toward IMP12, so
    # positive NES = enriched in the IMP12 direction, negative = toward others.
    p_g <- p_g +
      labs(x = paste0("NES   (positive → enriched in ", cfg$focus_sample,
                      ";  negative → enriched in others)"),
           subtitle = paste0("Ranking points toward ", cfg$focus_sample,
                             " separation on PC1/PC2")) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.02, vjust = -0.6,
               label = paste0("→ ", cfg$focus_sample, "-like"),
               colour = "#E31A1C", fontface = "bold", size = 4.2) +
      annotate("text", x = -Inf, y = -Inf, hjust = -0.02, vjust = -0.6,
               label = "Other-like ←",
               colour = "#2166AC", fontface = "bold", size = 4.2) +
      coord_cartesian(clip = "off")
    save_png(p_g, "13_pseudobulk_PCA_GSEA_hallmark", plot_dir, 13, 11)
  }

  message("  Pseudobulk PCA GSEA: ",
          sum(!is.na(gsea_pca$padj) & gsea_pca$padj < 0.05), " sig pathways")
  invisible(list(pca = pca_res, gsea = gsea_pca))
}

# =============================================================================
## STAGE: deseq2 — Pseudobulk DESeq2 IMP12 vs others + fgsea (Slide 20)
# =============================================================================
stage_deseq2 <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[deseq2] Pseudobulk DESeq2 (", cfg$focus_sample,
          " vs others, pool-tagged replicates)...")

  counts <- safe_counts(obj, "RNA")
  md     <- obj@meta.data

  # Pool-tagged pseudobulk ID: IMP12_WHL2, IMP12_WHL3, IMP11_WHL2, etc.
  md$pb_id <- paste0(as.character(md[[cfg$sample_col]]), "_",
                     toupper(as.character(md[[cfg$pool_col]])))

  pb_ids <- sort(unique(md$pb_id))
  pb_mat <- vapply(pb_ids, function(id) {
    cells <- rownames(md)[md$pb_id == id]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(pb_mat) <- rownames(counts)
  colnames(pb_mat) <- pb_ids

  # Filter: expressed in ≥20% of pseudosamples
  min_n  <- max(2, round(0.2 * ncol(pb_mat)))
  keep   <- rowSums(pb_mat > 0) >= min_n
  pb_mat <- pb_mat[keep, ]
  message("  Pseudosamples: ", ncol(pb_mat), "  Genes: ", nrow(pb_mat))

  # Export for PROMOTER pipeline
  write.table(pb_mat, file.path(data_dir, "pseudobulk_pooltagged_raw_counts.txt"),
              sep = "\t", quote = FALSE)

  # colData
  col_df <- data.frame(
    pb_id  = pb_ids,
    sample = sub("_WHL[0-9]+$", "", pb_ids),
    pool   = sub(".*_(WHL[0-9]+)$", "\\1", pb_ids),
    stringsAsFactors = FALSE
  )
  col_df$group <- factor(
    ifelse(col_df$sample == cfg$focus_sample, cfg$focus_sample, "Other"),
    levels = c("Other", cfg$focus_sample)
  )
  rownames(col_df) <- col_df$pb_id
  col_df <- col_df[colnames(pb_mat), , drop = FALSE]

  n_focus <- sum(col_df$group == cfg$focus_sample)
  n_other <- sum(col_df$group == "Other")
  message("  ", cfg$focus_sample, ": n=", n_focus, "  Other: n=", n_other)

  if (n_focus < 2) {
    warning("  Only ", n_focus, " replicate(s) for ", cfg$focus_sample,
            " — DESeq2 requires ≥2. Skipping deseq2 stage.")
    return(invisible(NULL))
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(pb_mat)),
    colData   = col_df,
    design    = ~ group
  )
  dds <- DESeq(dds, quiet = FALSE)
  res <- results(dds, contrast = c("group", cfg$focus_sample, "Other"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

  write.csv(res_df,
            file.path(data_dir, paste0("deseq2_", cfg$focus_sample,
                                        "_vs_Others_results.csv")),
            row.names = FALSE)

  n_up   <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0)
  n_down <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < 0)
  message("  DESeq2: ", n_up, " up, ", n_down, " down in ", cfg$focus_sample,
          " (padj<0.05)")

  # Volcano
  vdf <- res_df[!is.na(res_df$pvalue), ]
  vdf$neglogp <- -log10(vdf$pvalue + 1e-300)
  vdf$sig <- "NS"
  vdf$sig[vdf$log2FoldChange >  0.585 & !is.na(vdf$padj) & vdf$padj < 0.05] <-
    paste0("Up in ", cfg$focus_sample)
  vdf$sig[vdf$log2FoldChange < -0.585 & !is.na(vdf$padj) & vdf$padj < 0.05] <-
    paste0("Down in ", cfg$focus_sample)

  top_lab <- head(vdf[order(vdf$neglogp, decreasing = TRUE), ], 20)

  p_v <- ggplot(vdf, aes(x = log2FoldChange, y = neglogp, color = sig)) +
    geom_point(size = 0.7, alpha = 0.5) +
    geom_text_repel(data = top_lab, aes(label = gene), size = 3,
                    max.overlaps = 20, color = "black") +
    scale_color_manual(values = c(
      "NS" = "grey80",
      setNames("#E31A1C", paste0("Up in ",   cfg$focus_sample)),
      setNames("#2166AC", paste0("Down in ", cfg$focus_sample))
    ), name = NULL) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    theme_cowplot(13) +
    labs(title = paste0("Pseudobulk DESeq2: ", cfg$focus_sample, " vs Others"),
         x = "log2 Fold Change", y = "-log10(p-value)")
  save_png(p_v, "14a_deseq2_IMP12_vs_Others_volcano", plot_dir, 10, 8)

  # GSEA on DESeq2 Wald statistic
  stat_ranks <- setNames(res_df$stat, res_df$gene)
  stat_ranks <- sort(stat_ranks[!is.na(stat_ranks)], decreasing = TRUE)

  pathways_h  <- load_hallmark()
  gsea_deseq  <- fgsea(pathways = pathways_h, stats = stat_ranks,
                        nPermSimple = cfg$gsea_nperm, minSize = 15, maxSize = 500)
  gsea_deseq  <- gsea_deseq[order(gsea_deseq$NES, decreasing = TRUE), ]
  write.csv(gsea_deseq[, c("pathway","pval","padj","NES","size")],
            file.path(data_dir, paste0("deseq2_", cfg$focus_sample,
                                        "_vs_Others_GSEA_hallmark.csv")),
            row.names = FALSE)

  p_g <- gsea_dotplot(gsea_deseq,
                       paste0("GSEA — DESeq2 ", cfg$focus_sample,
                              " vs Others (Hallmark)"),
                       cfg$gsea_top_n)
  if (!is.null(p_g))
    save_png(p_g, "14b_deseq2_IMP12_vs_Others_GSEA_hallmark", plot_dir, 13, 11)

  message("  DESeq2 GSEA: ",
          sum(!is.na(gsea_deseq$padj) & gsea_deseq$padj < 0.05), " sig pathways")
  invisible(list(deseq2 = res_df, gsea = gsea_deseq))
}

# =============================================================================
## STAGE: features — CX3CR1, EGR1, PBMC canonical panel (Slide 23)
# =============================================================================
PBMC_MARKERS <- list(
  T_cell   = c("CD3D","CD3E","TRAC"),
  CD4_T    = c("CD4","IL7R","CCR7","SELL"),
  CD8_T    = c("CD8A","CD8B","GZMK","GZMB"),
  NK       = c("GNLY","NKG7","NCAM1","KLRD1"),
  NKT      = c("NKG7","CD3D","FCGR3A"),
  B_cell   = c("CD19","MS4A1","CD79A"),
  Mono_cl  = c("CD14","LYZ","CST3","S100A9"),
  Mono_nc  = c("FCGR3A","MS4A7","CX3CR1"),
  DC       = c("FCER1A","CD1C","CLEC4C"),
  IFN_stim = c("ISG15","MX1","IFIT1","RSAD2"),
  Cycling  = c("MKI67","TOP2A"),
  Key      = c("CX3CR1","EGR1")
)

# Canonical PBMC markers shown as a multi-reduction FeaturePlot panel (Image #3)
KEY_MARKER_PANEL <- c("CD3E","CD4","CD8A","NKG7","LYZ","CD68","CD14",
                      "FCGR3A","FCER1A","CD79A","MZB1","IL3RA")

stage_features <- function(obj, cfg, plot_dir) {
  message("\n[features] Feature + violin plots: CX3CR1, EGR1, PBMC canonical panel...")
  DefaultAssay(obj) <- "SCT"
  obj <- ensure_is_focus(obj, cfg)

  # Use umap_harmony if present, otherwise fall back to umap
  reduc <- if ("umap_harmony" %in% names(obj@reductions)) "umap_harmony" else "umap"
  message("  Using reduction: ", reduc)

  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  samp_pal <- make_sample_pal(samples, cfg$focus_sample)

  # Slide 23: CX3CR1 and EGR1
  for (gene in c("CX3CR1", "EGR1")) {
    if (!gene %in% rownames(obj)) {
      message("  ", gene, " not found — skipping."); next
    }
    # FeaturePlot split by is_focus
    p_fp <- FeaturePlot(obj, features = gene, reduction = reduc,
                         split.by = "is_focus", order = TRUE, pt.size = 0.4) +
      plot_annotation(title = paste0(gene, ": ", cfg$focus_sample, " vs Others"))
    save_png(p_fp, paste0("15a_featureplot_", gene), plot_dir, cfg$split_w, 7)

    # VlnPlot by IMP sample
    obj@meta.data[[cfg$sample_col]] <- factor(obj@meta.data[[cfg$sample_col]],
                                               levels = samples)
    p_vln <- VlnPlot(obj, features = gene, group.by = cfg$sample_col,
                      cols = samp_pal, pt.size = 0, assay = "SCT") +
      ggtitle(paste0(gene, " by IMP Sample")) + NoLegend() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.title.x = element_blank())
    save_png(p_vln, paste0("15b_violin_", gene, "_by_IMP_sample"),
             plot_dir, max(14, length(samples) * 0.5), 5)
  }

  # PBMC canonical panel — FeaturePlot
  all_genes <- unique(unlist(PBMC_MARKERS))
  present   <- intersect(all_genes, rownames(obj))
  message("  PBMC panel genes present: ", length(present), "/", length(all_genes))

  if (length(present) > 0) {
    ncol_fp <- 5
    nrow_fp <- ceiling(length(present) / ncol_fp)
    p_panel <- (FeaturePlot(obj, features = present, reduction = reduc,
                            ncol = ncol_fp, order = TRUE, pt.size = 0.2,
                            raster = TRUE) & theme_cowplot(9)) +
      plot_annotation(subtitle = umap_subtitle(cfg))
    save_png(p_panel, "16_featureplot_PBMC_canonical_panel",
             plot_dir, 25, 4 * nrow_fp)

    # Slide 17 (was a dotplot): canonical markers as a scaled mean-expression
    # heatmap per Harmony cluster — large italic gene rows for cell-type calling.
    avg <- AverageExpression(obj, assays = "SCT", features = present,
                             group.by = "seurat_clusters_harmony",
                             layer = "data")$SCT
    avg <- as.matrix(avg)
    z <- t(scale(t(avg)))           # z-score each gene across clusters
    z[is.na(z)] <- 0
    cl_ord <- order(suppressWarnings(as.integer(colnames(z))))
    if (anyNA(suppressWarnings(as.integer(colnames(z))))) cl_ord <- order(colnames(z))
    z <- z[, cl_ord, drop = FALSE]

    zdf <- data.frame(
      gene    = factor(rep(rownames(z), times = ncol(z)), levels = rev(present)),
      cluster = factor(rep(colnames(z), each = nrow(z)), levels = colnames(z)),
      z       = as.vector(z),
      stringsAsFactors = FALSE
    )
    p_hm <- ggplot(zdf, aes(x = cluster, y = gene, fill = z)) +
      geom_tile(colour = "grey92") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Scaled\nmean expr") +
      labs(title = "Canonical PBMC markers by Harmony cluster",
           subtitle = umap_subtitle(cfg), x = "Harmony cluster", y = NULL) +
      theme_cowplot(13) +
      theme(axis.text.y = element_text(size = 12, face = "italic"))
    save_png(p_hm, "17_canonical_marker_heatmap_by_cluster", plot_dir,
             width  = max(8, ncol(z) * 0.6 + 2),
             height = max(8, length(present) * 0.3 + 1))
  }

  # ── Key marker panel (Image #3) across three reductions ──────────────────────
  key_present <- intersect(KEY_MARKER_PANEL, rownames(obj))
  message("  Key panel genes present: ", length(key_present), "/",
          length(KEY_MARKER_PANEL))

  if (length(key_present) > 0) {
    save_key_panel <- function(o, red, tag, ttl, sub) {
      nc <- 4; nr <- ceiling(length(key_present) / nc)
      p <- (FeaturePlot(o, features = key_present, reduction = red, ncol = nc,
                        order = TRUE, pt.size = 0.2, raster = TRUE,
                        cols = c("lightgrey", "#08306B")) & theme_cowplot(9)) +
        plot_annotation(title = ttl, subtitle = sub)
      save_png(p, paste0("19_keypanel_", tag), plot_dir,
               width = 4 * nc, height = 3.2 * nr + 0.6)
    }

    if ("umap" %in% names(obj@reductions))
      save_key_panel(obj, "umap", "uncorr",
                     "Key markers — uncorrected (SCT) UMAP", umap_subtitle(cfg))
    if ("umap_harmony" %in% names(obj@reductions))
      save_key_panel(obj, "umap_harmony", "harmony",
                     "Key markers — Harmony-corrected UMAP", umap_subtitle(cfg))

    # IMP12-only UMAP (recompute on the focus subset using existing PCA)
    focus_cells <- rownames(obj@meta.data)[
      as.character(obj@meta.data[[cfg$sample_col]]) == cfg$focus_sample]
    if (length(focus_cells) > 20 && "pca" %in% names(obj@reductions)) {
      of <- subset(obj, cells = focus_cells)
      of <- RunUMAP(of, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)
      save_key_panel(of, "umap", "IMP12only",
                     paste0("Key markers — ", cfg$focus_sample, " only UMAP"),
                     umap_subtitle(cfg, show_res = FALSE))
      rm(of)
    }
  }

  invisible(NULL)
}

# =============================================================================
## STAGE: compare — raw vs SCT vs SCT+Harmony clustered UMAPs (same params)
# =============================================================================
stage_compare_umaps <- function(obj, cfg, plot_dir) {
  message("\n[compare] Raw vs SCT vs SCT+Harmony UMAPs grouped by ",
          cfg$sample_col, "...")

  have_sct  <- "umap" %in% names(obj@reductions)
  have_harm <- "umap_harmony" %in% names(obj@reductions)
  if (!have_sct || !have_harm) {
    message("  Need uncorrected + Harmony UMAPs first — skipping comparison.")
    return(invisible(NULL))
  }

  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  samp_pal <- make_sample_pal(samples, cfg$focus_sample)

  # ── Raw branch: standard log-normalize → PCA → UMAP (no SCT/Harmony),
  #    using the SAME dims as the SCT/Harmony branches. ───────────────────────
  raw <- obj
  DefaultAssay(raw) <- "RNA"
  raw <- NormalizeData(raw, verbose = FALSE)
  raw <- FindVariableFeatures(raw, selection.method = "vst",
                              nfeatures = cfg$n_hvg, verbose = FALSE)
  VariableFeatures(raw) <- filter_hvg(VariableFeatures(raw), cfg$exclude_patterns)
  raw <- ScaleData(raw, verbose = FALSE)
  raw <- RunPCA(raw, npcs = cfg$n_pcs, reduction.name = "pca_raw", verbose = FALSE)
  raw <- RunUMAP(raw, reduction = "pca_raw", dims = cfg$pca_dims,
                 reduction.name = "umap_raw", verbose = FALSE)

  base_thm <- theme_cowplot(12)
  p_raw <- DimPlot(raw, reduction = "umap_raw", group.by = cfg$sample_col,
                   cols = samp_pal, pt.size = 0.3) +
    labs(title = "Raw (log-normalized)") + base_thm + NoLegend()
  p_sct <- DimPlot(obj, reduction = "umap", group.by = cfg$sample_col,
                   cols = samp_pal, pt.size = 0.3) +
    labs(title = "SCTransform (uncorrected)") + base_thm + NoLegend()
  p_harm <- DimPlot(obj, reduction = "umap_harmony", group.by = cfg$sample_col,
                    cols = samp_pal, pt.size = 0.3) +
    labs(title = "SCTransform + Harmony") + base_thm

  p <- (p_raw | p_sct | p_harm) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste0("UMAP comparison — colored by ", cfg$sample_col),
      subtitle = umap_subtitle(cfg, show_res = FALSE))
  save_png(p, "20_umap_comparison_raw_sct_harmony", plot_dir,
           width = 3 * cfg$umap_w + 2, height = cfg$umap_h)
  rm(raw)
  invisible(NULL)
}

# =============================================================================
## STAGE: anchors — pairwise normalized integration-anchor heatmap
# =============================================================================
stage_anchor_heatmap <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[anchors] Pairwise normalized integration-anchor heatmap...")
  DefaultAssay(obj) <- "RNA"

  obj_list <- SplitObject(obj, split.by = cfg$sample_col)
  keep     <- vapply(obj_list, ncol, integer(1)) >= 50
  obj_list <- obj_list[keep]
  if (length(obj_list) < 4) {
    message("  <4 samples with >=50 cells — skipping anchor heatmap.")
    return(invisible(NULL))
  }
  obj_list <- lapply(obj_list, function(o) {
    DefaultAssay(o) <- "RNA"
    o <- NormalizeData(o, verbose = FALSE)
    FindVariableFeatures(o, selection.method = "vst",
                         nfeatures = cfg$n_hvg, verbose = FALSE)
  })

  features <- SelectIntegrationFeatures(obj_list, nfeatures = cfg$n_hvg,
                                        verbose = FALSE)
  min_cells <- min(vapply(obj_list, ncol, integer(1)))
  anchors <- tryCatch(
    FindIntegrationAnchors(
      obj_list, anchor.features = features,
      normalization.method = "LogNormalize", dims = cfg$pca_dims,
      k.filter = min(200, max(5, min_cells - 1)), verbose = TRUE),
    error = function(e) { message("  FindIntegrationAnchors failed: ",
                                  conditionMessage(e)); NULL })
  if (is.null(anchors)) return(invisible(NULL))

  ad       <- as.data.frame(slot(anchors, "anchors"))
  ds_names <- names(obj_list)
  n        <- length(ds_names)
  cells_n  <- vapply(obj_list, ncol, integer(1))

  # Count anchors per (unordered) sample pair → symmetric matrix
  M <- matrix(0, n, n, dimnames = list(ds_names, ds_names))
  for (r in seq_len(nrow(ad))) {
    i <- ad$dataset1[r]; j <- ad$dataset2[r]
    M[i, j] <- M[i, j] + 1
    M[j, i] <- M[j, i] + 1
  }

  # Size-normalize (anchor counts scale with dataset size), then center on the
  # median off-diagonal value so typical pairs sit near 1.0 (as in the example).
  Msz     <- M / outer(sqrt(cells_n), sqrt(cells_n))
  offdiag <- Msz[upper.tri(Msz)]
  med     <- median(offdiag[offdiag > 0])
  norm    <- if (is.finite(med) && med > 0) Msz / med else Msz
  diag(norm) <- max(norm[upper.tri(norm)], na.rm = TRUE)

  write.csv(data.frame(norm, check.names = FALSE),
            file.path(data_dir, "pairwise_normalized_anchor_counts.csv"))

  pal  <- colorRampPalette(c("#FFFFCC", "#FD8D3C", "#E31A1C"))(100)
  path <- file.path(plot_dir, "21_pairwise_anchor_heatmap.png")
  png(path, width = 11, height = 10, units = "in", res = 300)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(norm, color = pal,
                       main = "Pairwise normalized anchor counts",
                       treeheight_row = 50, treeheight_col = 50,
                       border_color = NA, fontsize = 11)
  } else {
    heatmap(norm, col = pal, symm = TRUE,
            main = "Pairwise normalized anchor counts")
  }
  dev.off()
  message("  Saved: 21_pairwise_anchor_heatmap.png")
  invisible(norm)
}

# =============================================================================
## MAIN
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Starting pipeline: ", paste(active_stages, collapse = ", "))
message(paste(rep("=", 60), collapse = ""))

# Determine load point
obj <- NULL

if (run_stage("sct")) {
  message("[load] Reading: ", cfg$input_rds)
  obj <- readRDS(cfg$input_rds)
  DefaultAssay(obj) <- "RNA"
  message("  Cells: ", ncol(obj), "  Genes: ", nrow(obj))
  print(table(obj@meta.data[[cfg$sample_col]]))

  # Preserve original cluster labels before FindClusters overwrites seurat_clusters
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    obj$clusters_premerge <- obj$seurat_clusters
    message("  Original clusters saved to: clusters_premerge (",
            length(unique(obj$clusters_premerge)), " clusters)")
  }
} else if (run_stage("umap_uncorr")) {
  message("[load] Resuming from SCT checkpoint: ", ckpt_sct)
  obj <- readRDS(ckpt_sct)
} else if (run_stage("harmony")) {
  message("[load] Resuming from uncorrected UMAP checkpoint: ", ckpt_uncorr)
  obj <- readRDS(ckpt_uncorr)
} else {
  message("[load] Resuming from Harmony checkpoint: ", ckpt_harmony)
  harmony_ckpt_data <- readRDS(ckpt_harmony)
  # Checkpoint now stores the Seurat object directly; tolerate the old
  # list(obj=, missing_clusters=) format for backward compatibility.
  if (is.list(harmony_ckpt_data) && "obj" %in% names(harmony_ckpt_data)) {
    obj <- harmony_ckpt_data$obj
  } else {
    obj <- harmony_ckpt_data
  }
}

# For qc, we can always run it as long as obj is loaded
if (run_stage("qc"))          stage_qc(obj, cfg, plot_dir)
if (run_stage("sct"))         obj <- stage_sct(obj, cfg, plot_dir, data_dir, ckpt_sct)
if (run_stage("umap_uncorr")) obj <- stage_umap_uncorr(obj, cfg, plot_dir, data_dir, ckpt_uncorr)
if (run_stage("harmony"))     obj <- stage_harmony(obj, cfg, plot_dir, data_dir, ckpt_harmony)
if (run_stage("compare"))
  tryCatch(stage_compare_umaps(obj, cfg, plot_dir),
           error = function(e) message("  [compare] failed: ", conditionMessage(e)))
if (run_stage("markers"))  stage_markers(obj, cfg, plot_dir, data_dir)
if (run_stage("anchors"))
  tryCatch(stage_anchor_heatmap(obj, cfg, plot_dir, data_dir),
           error = function(e) message("  [anchors] failed: ", conditionMessage(e)))
if (run_stage("pb_pca"))   stage_pb_pca(obj, cfg, plot_dir, data_dir)
if (run_stage("deseq2"))   stage_deseq2(obj, cfg, plot_dir, data_dir)
if (run_stage("features")) stage_features(obj, cfg, plot_dir)

message("\n", paste(rep("=", 60), collapse = ""))
message("Pipeline complete. Outputs: ", out_root)
message(paste(rep("=", 60), collapse = ""))
