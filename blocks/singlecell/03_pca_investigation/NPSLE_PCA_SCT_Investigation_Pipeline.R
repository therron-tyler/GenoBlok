#!/usr/bin/env Rscript
# =============================================================================
# NPSLE_PCA_SCT_Investigation_Pipeline.R
#
# Investigate which genes drive the principal components BEFORE and AFTER
# SCTransform v2, with every plot resolved by individual sample (HTO_maxID),
# then sweep clustering across user-supplied dims x resolutions and emit a
# de novo marker heatmap for each combination.
#
# Output subfolder: <out_dir>/PCA_SCT_Investigation/{plots,data}/
#
# Stages (comma-separated in --stages, or "all"):
#   presct   - LogNormalize + cell-cycle scoring + ScaleData + PCA  (reduction "pca_lognorm")
#   sct      - SCTransform v2 (n_hvg HVGs) + PCA                    (reduction "pca_sct")
#   pcviz    - For each available reduction: per-PC-pair sample-coloured PCA,
#              diverging top-loadings barplots, DimHeatmap, biplot with loading
#              arrows, and per-sample PC-score violins (genes -> PCs -> samples)
#   cluster  - For each available reduction (pca_lognorm and/or pca_sct) and each
#              (dim, resolution): FindNeighbors/FindClusters/UMAP, cluster +
#              sample UMAPs, cluster-by-sample composition, and a de novo top-N
#              marker heatmap (Seurat purple/yellow palette). Outputs are tagged
#              with the reduction (e.g. cluster_presct_d10_res0.4_umap.png vs
#              cluster_sct_d10_res0.4_umap.png) so the two embeddings can be
#              compared side-by-side.
#
# Example QUEST sbatch call:
#   Rscript NPSLE_PCA_SCT_Investigation_Pipeline.R \
#     --input_rds  /path/to/.../20260213_GroupComps_NPSLEsamples/NPSLE_Tcell_Subcluster.rds \
#     --out_dir    /path/to/.../NPSLE_PCA_SCT \
#     --n_pcs      30 \
#     --n_viz_pcs  10 \
#     --dims       8,10,15 \
#     --resolutions 0.3,0.5,0.8 \
#     --stages     all
# =============================================================================

.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(sctransform)
  library(glmGamPoi)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  library(Matrix)
  library(optparse)
})

# =============================================================================
## CLI
# =============================================================================
option_list <- list(
  make_option("--input_rds",    type = "character", default = NULL,
              help = "Input Seurat .rds (REQUIRED)"),
  make_option("--out_dir",      type = "character", default = NULL,
              help = "Output directory (REQUIRED)"),
  make_option("--sample_col",   type = "character", default = "HTO_maxID",
              help = "Metadata column with per-sample IDs [default: HTO_maxID]"),
  make_option("--n_hvg",        type = "integer",   default = 3000L,
              help = "Highly variable genes for SCT / lognorm [default: 3000]"),
  make_option("--n_pcs",        type = "integer",   default = 30L,
              help = "PCs to compute [default: 30]"),
  make_option("--n_viz_pcs",    type = "integer",   default = 10L,
              help = "Leading PCs to visualise as pairs (1,2)(3,4)... [default: 10]"),
  make_option("--top_loadings", type = "integer",   default = 15L,
              help = "Top +/- genes per PC in loadings barplots [default: 15]"),
  make_option("--biplot_genes", type = "integer",   default = 10L,
              help = "Loading-vector arrows drawn on each biplot [default: 10]"),
  make_option("--regress_cc",   type = "logical",   default = TRUE,
              help = "Regress percent.mt + S/G2M scores in ScaleData/SCT [default: TRUE]"),
  make_option("--dims",         type = "character", default = "10",
              help = "Comma list of PCA dims to sweep, each used as 1:d [default: 10]"),
  make_option("--resolutions",  type = "character", default = "0.5",
              help = "Comma list of clustering resolutions to sweep [default: 0.5]"),
  make_option("--heatmap_top_n",type = "integer",   default = 5L,
              help = "Top markers per cluster in de novo heatmap [default: 5]"),
  make_option("--from_checkpoint", type = "character", default = NULL,
              help = "Post-SCT checkpoint .rds to resume from (skips presct/sct)"),
  make_option("--gene_violin_pcs", type = "character", default = "8",
              help = "Comma list of PCs to emit per-sample driver-gene violins for; '' = none [default: 8]"),
  make_option("--inject_rds",   type = "character", default = NULL,
              help = "Optional .rds whose cells matching --inject_samples are merged into the input object before processing"),
  make_option("--inject_samples", type = "character", default = NULL,
              help = "Comma list of sample IDs (matched on --sample_col) to extract from --inject_rds and merge in"),
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
  sample_col   = opt$sample_col,
  n_hvg        = opt$n_hvg,
  n_pcs        = opt$n_pcs,
  n_viz_pcs    = opt$n_viz_pcs,
  top_loadings = opt$top_loadings,
  biplot_genes = opt$biplot_genes,
  regress_cc   = opt$regress_cc,
  heatmap_top_n= opt$heatmap_top_n,

  dims         = as.integer(trimws(strsplit(opt$dims, ",")[[1]])),
  resolutions  = as.numeric(trimws(strsplit(opt$resolutions, ",")[[1]])),
  gene_violin_pcs = if (nzchar(opt$gene_violin_pcs))
                      as.integer(trimws(strsplit(opt$gene_violin_pcs, ",")[[1]]))
                    else integer(0),
  inject_rds      = opt$inject_rds,
  inject_samples  = if (!is.null(opt$inject_samples) && nzchar(opt$inject_samples))
                      trimws(strsplit(opt$inject_samples, ",")[[1]])
                    else character(0),

  # Strip technical / receptor genes so PCs reflect biology, not chemistry.
  exclude_patterns = c("^MT-", "^MTRNR", "^RPS", "^RPL",
                       "^IGHV", "^IGLV", "^IGKV",
                       "^TRAV", "^TRBV", "^TRDV", "^TRGV")
)
cfg$regress_vars <- if (cfg$regress_cc)
  c("percent.mt", "S.Score", "G2M.Score") else character(0)

ALL_STAGES <- c("presct", "sct", "pcviz", "cluster", "canon")
active_stages <- if (grepl("^all$", trimws(opt$stages), ignore.case = TRUE)) ALL_STAGES else
  trimws(strsplit(opt$stages, ",")[[1]])
run_stage <- function(s) s %in% active_stages

out_root <- file.path(opt$out_dir, "PCA_SCT_Investigation")
plot_dir <- file.path(out_root, "plots")
data_dir <- file.path(out_root, "data")
for (d in c(out_root, plot_dir, data_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

ckpt_sct <- file.path(out_root, "ckpt_post_sct.rds")

message("== NPSLE PCA / SCTv2 Investigation Pipeline ==")
message("  input_rds   : ", cfg$input_rds)
message("  out_dir     : ", out_root)
message("  sample_col  : ", cfg$sample_col)
message("  stages      : ", paste(active_stages, collapse = ", "))
message("  n_hvg       : ", cfg$n_hvg, "  n_pcs: ", cfg$n_pcs,
        "  viz PCs: 1-", cfg$n_viz_pcs)
message("  dims sweep  : ", paste(cfg$dims, collapse = ", "),
        "   res sweep: ", paste(cfg$resolutions, collapse = ", "))
message("  regress     : ",
        if (length(cfg$regress_vars)) paste(cfg$regress_vars, collapse = ", ") else "(none)")
message("  gene violin PCs: ",
        if (length(cfg$gene_violin_pcs)) paste(cfg$gene_violin_pcs, collapse = ", ") else "(none)")
message("  inject_rds  : ",
        if (!is.null(cfg$inject_rds)) cfg$inject_rds else "(none)",
        "  samples: ",
        if (length(cfg$inject_samples)) paste(cfg$inject_samples, collapse = ", ") else "(none)")

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

# Align the surface-protein (ADT) assay across two objects so a Seurat v5
# merge() preserves it — merge keeps only assays common to BOTH objects, so an
# ADT assay living on just one side is silently dropped. When both carry ADT we
# restrict to shared antibodies; when only one does, the other is zero-filled
# (its cells genuinely have no surface measurement). Returns the ADT-compatible
# pair plus the shared antibody count. Each rebuilt assay keeps only a counts
# layer; the caller re-normalizes (CLR) after the merge.
align_adt_for_merge <- function(a, b) {
  ha <- "ADT" %in% Assays(a); hb <- "ADT" %in% Assays(b)
  if (!ha && !hb) return(list(a = a, b = b, n_feat = 0L, both = FALSE))
  feats <- if (ha && hb) intersect(rownames(a[["ADT"]]), rownames(b[["ADT"]])) else
           if (ha) rownames(a[["ADT"]]) else rownames(b[["ADT"]])
  if (length(feats) == 0) return(list(a = a, b = b, n_feat = 0L, both = FALSE))

  adt_counts <- function(o) {
    asy  <- o[["ADT"]]
    lyrs <- tryCatch(SeuratObject::Layers(asy), error = function(e) character(0))
    lyr  <- if ("counts" %in% lyrs) "counts" else if (length(lyrs)) lyrs[[1]] else "data"
    as.matrix(SeuratObject::LayerData(asy, layer = lyr))[feats, , drop = FALSE]
  }
  rebuild <- function(o, has) {
    cnt <- if (has) adt_counts(o) else
           matrix(0, length(feats), ncol(o), dimnames = list(feats, colnames(o)))
    o[["ADT"]] <- CreateAssay5Object(counts = as(cnt, "CsparseMatrix"))
    o
  }
  list(a = rebuild(a, ha), b = rebuild(b, hb), n_feat = length(feats), both = ha && hb)
}

filter_hvg <- function(genes, patterns) {
  keep <- rep(TRUE, length(genes))
  for (pat in patterns) keep <- keep & !grepl(pat, genes)
  genes[keep]
}

ensure_mito_ribo <- function(obj) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt <- grep("^MT-", rownames(obj), value = TRUE)
    obj[["percent.mt"]] <- if (length(mt) > 0)
      PercentageFeatureSet(obj, pattern = "^MT-") else 0
  }
  obj
}

# Stable sample palette (no focus sample; one colour per HTO_maxID).
make_sample_pal <- function(samples) {
  samples <- sort(unique(as.character(samples)))
  base_cols <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"),
                 brewer.pal(9, "Set1"))
  setNames(colorRampPalette(base_cols)(length(samples)), samples)
}

# Pairs (1,2)(3,4)... up to n_viz_pcs.
pc_pairs <- function(n) {
  starts <- seq(1, n - 1, by = 2)
  lapply(starts, function(i) c(i, i + 1))
}

# Diverging horizontal barplot of the top +/- loading genes for one PC.
pc_loading_barplot <- function(obj, reduction, pc, top_n, key) {
  load <- Loadings(obj, reduction = reduction)[, pc]
  ord  <- order(load)
  genes <- c(head(names(load)[ord], top_n), tail(names(load)[ord], top_n))
  df <- data.frame(gene = genes, loading = as.numeric(load[genes]))
  df$dir  <- ifelse(df$loading >= 0, "positive", "negative")
  df$gene <- factor(df$gene, levels = df$gene[order(df$loading)])
  ggplot(df, aes(x = loading, y = gene, fill = dir)) +
    geom_col() +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = c(positive = "#B2182B", negative = "#2166AC"),
                      guide = "none") +
    labs(title = paste0(key, pc, " loadings"), x = "loading", y = NULL) +
    theme_cowplot(10) +
    theme(axis.text.y = element_text(size = 8, face = "italic"))
}

# Sample-coloured scatter for one PC pair.
pca_pair_scatter <- function(obj, reduction, d1, d2, sample_col, pal, key) {
  DimPlot(obj, reduction = reduction, dims = c(d1, d2),
          group.by = sample_col, cols = pal, pt.size = 0.4) +
    labs(title = paste0("PCA ", key, d1, " vs ", key, d2, " (by ", sample_col, ")")) +
    theme_cowplot(12)
}

# Biplot: sample-coloured cells + top loading-vector arrows (genes -> PC axes).
pca_biplot <- function(obj, reduction, d1, d2, sample_col, pal, top_n, key) {
  emb  <- Embeddings(obj, reduction)[, c(d1, d2)]
  df   <- data.frame(x = emb[, 1], y = emb[, 2],
                     sample = obj@meta.data[[sample_col]])
  load <- Loadings(obj, reduction = reduction)[, c(d1, d2)]
  mag  <- sqrt(load[, 1]^2 + load[, 2]^2)
  top  <- names(sort(mag, decreasing = TRUE))[seq_len(min(top_n, length(mag)))]
  ld   <- data.frame(x = load[top, 1], y = load[top, 2], gene = top)
  scale_f <- 0.7 * max(abs(c(df$x, df$y))) / max(abs(as.matrix(ld[, 1:2])))
  ggplot(df, aes(x, y)) +
    geom_point(aes(colour = sample), size = 0.35, alpha = 0.55) +
    scale_colour_manual(values = pal, name = sample_col) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = x * scale_f, yend = y * scale_f),
                 arrow = arrow(length = unit(0.15, "cm")), colour = "black",
                 linewidth = 0.4) +
    geom_text_repel(data = ld, aes(x = x * scale_f, y = y * scale_f, label = gene),
                    size = 3, fontface = "italic", colour = "black",
                    max.overlaps = Inf) +
    labs(title = paste0("Biplot ", key, d1, " vs ", key, d2),
         x = paste0(key, d1), y = paste0(key, d2)) +
    theme_cowplot(12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# Per-sample driver-gene violins for one PC: top +/- loading genes (same
# selection as pc_loading_barplot) shown as per-sample expression violins.
# layer = "data" gives log-normalised expression for RNA, or log of Pearson
# residuals for SCT — appropriate for cross-sample comparison.
pc_driver_gene_violins <- function(obj, reduction, pc, top_n, sample_col, pal, key, assay) {
  load <- Loadings(obj, reduction = reduction)[, pc]
  ord  <- order(load)
  genes <- c(head(names(load)[ord], top_n), tail(names(load)[ord], top_n))
  genes <- intersect(genes, rownames(obj[[assay]]))
  if (length(genes) == 0) return(NULL)

  expr <- FetchData(obj, vars = genes, layer = "data", assay = assay)
  expr$sample <- obj@meta.data[[sample_col]]
  dfl <- pivot_longer(expr, cols = all_of(genes),
                      names_to = "gene", values_to = "expr")
  dfl$gene <- factor(dfl$gene, levels = genes)   # preserves loading order (- -> +)

  ncols <- min(6L, length(genes))
  ggplot(dfl, aes(x = sample, y = expr, fill = sample)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white", linewidth = 0.2) +
    facet_wrap(~ gene, ncol = ncols, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title    = paste0("Per-sample expression of ", key, pc,
                           " top ", top_n, " +/- loading genes (", assay, ")"),
         subtitle = "Gene order: most negative loading -> most positive",
         x = NULL, y = "log expression") +
    theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          strip.text  = element_text(face = "italic", size = 10))
}

# Per-sample PC-score violins: which samples separate along each PC.
pc_score_violins <- function(obj, reduction, pcs, sample_col, pal, key) {
  emb <- Embeddings(obj, reduction)[, pcs, drop = FALSE]
  colnames(emb) <- paste0(key, pcs)
  df <- data.frame(emb, sample = obj@meta.data[[sample_col]], check.names = FALSE)
  dfl <- pivot_longer(df, cols = all_of(colnames(emb)),
                      names_to = "PC", values_to = "score")
  dfl$PC <- factor(dfl$PC, levels = paste0(key, pcs))
  ggplot(dfl, aes(x = sample, y = score, fill = sample)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white", linewidth = 0.2) +
    facet_wrap(~PC, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "Per-sample PC-score distributions", x = NULL, y = "PC score") +
    theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
}

# =============================================================================
## STAGE: presct — LogNormalize + cell cycle + ScaleData + PCA
# =============================================================================
stage_presct <- function(obj, cfg, plot_dir) {
  message("\n[presct] LogNormalize + cell-cycle scoring + ScaleData + PCA...")
  DefaultAssay(obj) <- "RNA"
  tryCatch(obj <- JoinLayers(obj, assay = "RNA"), error = function(e) NULL)
  obj <- ensure_mito_ribo(obj)

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(
    obj,
    s.features   = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident    = FALSE, verbose = FALSE
  )
  message("  Cell cycle: ", paste(names(table(obj$Phase)), table(obj$Phase),
                                   sep = "=", collapse = " | "))

  obj <- FindVariableFeatures(obj, nfeatures = cfg$n_hvg, verbose = FALSE)
  VariableFeatures(obj) <- filter_hvg(VariableFeatures(obj), cfg$exclude_patterns)
  obj <- ScaleData(obj, features = VariableFeatures(obj),
                   vars.to.regress = cfg$regress_vars, verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = cfg$n_pcs,
                reduction.name = "pca_lognorm", reduction.key = "PClog_",
                verbose = FALSE)

  p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs, reduction = "pca_lognorm") +
    ggtitle("Elbow (LogNormalize PCA, pre-SCT)") + theme_cowplot(14)
  save_png(p_elbow, "00_elbow_pca_lognorm", plot_dir, 8, 5)
  obj
}

# =============================================================================
## STAGE: sct — SCTransform v2 + PCA
# =============================================================================
stage_sct <- function(obj, cfg, plot_dir, data_dir, ckpt) {
  message("\n[sct] SCTransform v2 + PCA...")
  DefaultAssay(obj) <- "RNA"
  tryCatch(obj <- JoinLayers(obj, assay = "RNA"), error = function(e) NULL)
  obj <- ensure_mito_ribo(obj)

  # Cell-cycle scores are needed for regression even if presct was skipped.
  if (length(cfg$regress_vars) && !"S.Score" %in% colnames(obj@meta.data)) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- CellCycleScoring(obj,
                            s.features   = Seurat::cc.genes.updated.2019$s.genes,
                            g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
                            set.ident = FALSE, verbose = FALSE)
  }

  message("  SCTransform v2 (n_hvg=", cfg$n_hvg, ", regress: ",
          if (length(cfg$regress_vars)) paste(cfg$regress_vars, collapse = ", ") else "none",
          ")...")
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
  message("  HVGs: ", length(all_hvg), " -> ", length(cleaned_hvg),
          " after stripping MT/MTRNR/ribo/TCR/Ig")
  write.csv(data.frame(gene = cleaned_hvg),
            file.path(data_dir, "sct_hvg_list.csv"), row.names = FALSE)

  obj <- RunPCA(obj, features = cleaned_hvg, npcs = cfg$n_pcs,
                reduction.name = "pca_sct", reduction.key = "PCsct_",
                verbose = FALSE)

  p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs, reduction = "pca_sct") +
    ggtitle("Elbow (SCTransform v2 PCA)") + theme_cowplot(14)
  save_png(p_elbow, "01_elbow_pca_sct", plot_dir, 8, 5)

  saveRDS(obj, ckpt); message("  Checkpoint: ", basename(ckpt))
  obj
}

# =============================================================================
## STAGE: pcviz — PC-driver visualisations for each available reduction
# =============================================================================
viz_one_reduction <- function(obj, reduction, key, tag, assay, cfg, plot_dir, data_dir, pal) {
  message("  [pcviz] ", reduction, " ...")
  npc   <- min(cfg$n_viz_pcs, ncol(Embeddings(obj, reduction)))
  pairs <- pc_pairs(npc)

  # Export the loading matrix so PC drivers can be inspected outside R.
  write.csv(round(Loadings(obj, reduction = reduction)[, seq_len(npc)], 5),
            file.path(data_dir, paste0(tag, "_pc_loadings.csv")))

  for (pr in pairs) {
    d1 <- pr[1]; d2 <- pr[2]
    scat  <- pca_pair_scatter(obj, reduction, d1, d2, cfg$sample_col, pal, key)
    bar_x <- pc_loading_barplot(obj, reduction, d1, cfg$top_loadings, key)
    bar_y <- pc_loading_barplot(obj, reduction, d2, cfg$top_loadings, key)
    combo <- scat / (bar_x | bar_y) +
      plot_layout(heights = c(1.3, 1)) +
      plot_annotation(title = paste0(tag, ": ", key, d1, "/", key, d2,
                                     " scores + driver genes"))
    save_png(combo, paste0(tag, "_", key, d1, "_", key, d2, "_drivers"),
             plot_dir, width = 11, height = 12)

    bip <- pca_biplot(obj, reduction, d1, d2, cfg$sample_col, pal,
                      cfg$biplot_genes, key)
    save_png(bip, paste0(tag, "_", key, d1, "_", key, d2, "_biplot"),
             plot_dir, width = 9, height = 8)
  }

  # Seurat DimHeatmap: cells (ranked by PC score) x top loading genes.
  # Use the reduction's own assay — pca_lognorm genes live in RNA scale.data,
  # pca_sct genes in SCT scale.data; mixing them gives "undefined columns".
  # Chunked into groups of up to 4 PCs/panel so gene labels stay readable.
  tryCatch({
    n_cells_hm <- min(500L, ncol(obj))
    chunks     <- split(seq_len(npc), ceiling(seq_len(npc) / 4L))
    for (pcs_c in chunks) {
      lo <- min(pcs_c); hi <- max(pcs_c)
      dh_list <- DimHeatmap(obj, reduction = reduction, dims = pcs_c,
                            cells = n_cells_hm, balanced = TRUE, assays = assay,
                            fast = FALSE, combine = FALSE)
      dh_list <- lapply(dh_list, function(p) p +
        theme(axis.text.y = element_text(size = 14, face = "italic"),
              plot.title  = element_text(size = 14)))
      dh <- patchwork::wrap_plots(dh_list, ncol = 1) +
        patchwork::plot_annotation(
          title = paste0(tag, ": DimHeatmap PC", lo, "-", hi, " (", assay, ")"))
      save_png_device(dh,
                      paste0(tag, "_dimheatmap_PC", lo, "_", hi),
                      plot_dir,
                      width  = 14,
                      height = max(6, length(pcs_c) * 3.5))
    }
  }, error = function(e)
     message("  DimHeatmap (", reduction, ") skipped: ", conditionMessage(e)))

  # Per-sample PC-score violins: connects PCs to samples directly.
  save_png(pc_score_violins(obj, reduction, seq_len(npc), cfg$sample_col, pal, key),
           paste0(tag, "_persample_pc_score_violins"), plot_dir,
           width = max(10, length(pal) * 0.8), height = 2.2 * ceiling(npc / 3))

  # Per-sample driver-gene violins for each PC listed in --gene_violin_pcs.
  # One faceted figure per (reduction, PC): top +/- loading genes shown as
  # per-sample expression violins to surface batch-driven loading effects.
  for (pc in cfg$gene_violin_pcs) {
    if (pc < 1 || pc > npc) {
      message("  Skipping driver-gene violins for ", key, pc,
              " (out of range 1..", npc, ")")
      next
    }
    vplot <- pc_driver_gene_violins(obj, reduction, pc, cfg$top_loadings,
                                    cfg$sample_col, pal, key, assay)
    if (is.null(vplot)) {
      message("  No matching driver genes in assay ", assay, " for ", key, pc,
              " — skipping violins.")
      next
    }
    n_g   <- 2L * cfg$top_loadings
    ncols <- min(6L, n_g)
    save_png(vplot,
             paste0(tag, "_", key, pc, "_driver_gene_persample_violins"),
             plot_dir,
             width  = max(12, ncols * (1.5 + length(pal) * 0.18)),
             height = max(5,  2.0 * ceiling(n_g / ncols)))
  }
}

stage_pcviz <- function(obj, cfg, plot_dir, data_dir, pal) {
  message("\n[pcviz] PC-driver visualisations (before/after SCT)...")
  reds <- list()
  if ("pca_lognorm" %in% Reductions(obj))
    reds[["presct"]] <- list(reduction = "pca_lognorm", key = "PClog_", assay = "RNA")
  if ("pca_sct" %in% Reductions(obj))
    reds[["sct"]]    <- list(reduction = "pca_sct",     key = "PCsct_", assay = "SCT")
  if (length(reds) == 0) { message("  No PCA reductions present; skipping."); return(invisible(obj)) }
  for (tag in names(reds))
    viz_one_reduction(obj, reds[[tag]]$reduction, reds[[tag]]$key, tag,
                      reds[[tag]]$assay, cfg, plot_dir, data_dir, pal)
  invisible(obj)
}

# =============================================================================
## STAGE: cluster — sweep dims x resolutions, de novo heatmap per combination
# =============================================================================
cluster_one_reduction <- function(obj, reduction, assay, tag, cfg, plot_dir, data_dir, pal) {
  message("\n[cluster:", tag, "] sweeping dims x resolutions on ", reduction,
          " (assay = ", assay, ")...")
  DefaultAssay(obj) <- assay
  # SCT sweep follows the Seurat SCT DE tutorial: recorrect counts ONCE here,
  # then FindAllMarkers runs with recorrect_umi = FALSE so it doesn't re-correct.
  if (assay == "SCT")
    obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)

  umap_name <- paste0("umap_sweep_", tag)
  summary_rows <- list()
  for (d in cfg$dims) {
    obj <- FindNeighbors(obj, reduction = reduction, dims = seq_len(d), verbose = FALSE)
    obj <- RunUMAP(obj, reduction = reduction, dims = seq_len(d),
                   reduction.name = umap_name, verbose = FALSE)
    for (res in cfg$resolutions) {
      sub_tag <- paste0(tag, "_d", d, "_res", res)
      message("\n  -- ", sub_tag, " --")
      obj <- FindClusters(obj, resolution = res, verbose = FALSE)
      clcol <- paste0("clust_", sub_tag)
      obj[[clcol]] <- Idents(obj)
      n_cl <- length(levels(Idents(obj)))
      message("    clusters: ", n_cl)
      summary_rows[[sub_tag]] <- data.frame(reduction = tag, dims = d,
                                            resolution = res, n_clusters = n_cl)

      # Cluster + sample UMAPs
      p_cl <- DimPlot(obj, reduction = umap_name, label = TRUE, repel = TRUE) +
        labs(title = paste0("Clusters (", sub_tag, ")"),
             subtitle = paste0("reduction = ", reduction,
                               "  |  dims = 1:", d, "  |  res = ", res)) +
        theme_cowplot(12) + NoLegend()
      p_sm <- DimPlot(obj, reduction = umap_name, group.by = cfg$sample_col,
                      cols = pal) +
        labs(title = paste0("Sample (", cfg$sample_col, ")")) + theme_cowplot(12)
      save_png(p_cl | p_sm, paste0("cluster_", sub_tag, "_umap"), plot_dir, 15, 7)

      # Cluster-by-sample composition (tracks samples into clusters)
      comp <- as.data.frame(table(cluster = as.character(Idents(obj)),
                                   sample  = as.character(obj@meta.data[[cfg$sample_col]])),
                            stringsAsFactors = FALSE)
      comp$Freq <- as.integer(comp$Freq)
      comp <- comp %>% group_by(cluster) %>%
        mutate(pct = round(Freq / sum(Freq) * 100, 1)) %>% ungroup()
      write.csv(comp, file.path(data_dir, paste0("cluster_", sub_tag, "_sample_composition.csv")),
                row.names = FALSE)
      p_comp <- ggplot(comp, aes(x = cluster, y = pct, fill = sample)) +
        geom_col() + scale_fill_manual(values = pal, name = cfg$sample_col) +
        labs(title = paste0("Sample composition per cluster (", sub_tag, ")"),
             x = "cluster", y = "% of cluster") + theme_cowplot(12)
      save_png(p_comp, paste0("cluster_", sub_tag, "_sample_composition"), plot_dir, 11, 6)

      # De novo markers + heatmap (Seurat purple/yellow palette).
      # Pre-SCT: markers from RNA 'data' (LogNormalize). SCT: markers from the
      # recorrected SCT assay (recorrect_umi = FALSE — PrepSCTFindMarkers above
      # already corrected once). Each sweep's heatmap then uses ITS OWN assay so
      # the LogNormalize and SCTransform approaches can be compared head-to-head.
      markers <- if (assay == "SCT")
        FindAllMarkers(obj, assay = "SCT", recorrect_umi = FALSE, only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25,
                       test.use = "wilcox", verbose = FALSE)
      else
        FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                       min.pct = 0.1, logfc.threshold = 0.25,
                       test.use = "wilcox", verbose = FALSE)
      markers_sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, ]
      write.csv(markers_sig,
                file.path(data_dir, paste0("cluster_", sub_tag, "_markers_sig.csv")),
                row.names = FALSE)
      if (nrow(markers_sig) == 0) {
        message("    No significant markers — skipping heatmap."); next
      }
      topN <- markers_sig %>% group_by(cluster) %>%
        slice_max(avg_log2FC, n = cfg$heatmap_top_n, with_ties = FALSE) %>%
        pull(gene) %>% unique()
      p_heat <- DoHeatmap(obj, features = topN, assay = assay,
                          size = 4.5, angle = 45) +
        theme(axis.text.y = element_text(size = 11, face = "italic")) +
        ggtitle(paste0("De novo top ", cfg$heatmap_top_n,
                       " markers per cluster (", sub_tag, ")"))
      save_png_device(p_heat, paste0("cluster_", sub_tag, "_denovo_heatmap"), plot_dir,
                      width  = max(12, n_cl * 0.9),
                      height = min(40, max(8, length(topN) * 0.34)))
    }
  }
  if (length(summary_rows))
    write.csv(do.call(rbind, summary_rows),
              file.path(data_dir, paste0("cluster_", tag, "_sweep_summary.csv")),
              row.names = FALSE)
  invisible(obj)
}

stage_cluster <- function(obj, cfg, plot_dir, data_dir, pal) {
  reds <- list()
  if ("pca_lognorm" %in% Reductions(obj))
    reds[["presct"]] <- list(reduction = "pca_lognorm", assay = "RNA")
  if ("pca_sct" %in% Reductions(obj))
    reds[["sct"]]    <- list(reduction = "pca_sct",     assay = "SCT")
  if (length(reds) == 0) {
    message("\n[cluster] no PCA reductions present (need pca_lognorm or pca_sct); skipping.")
    return(invisible(obj))
  }
  message("\n[cluster] reductions to sweep: ", paste(names(reds), collapse = ", "))
  for (tag in names(reds))
    obj <- cluster_one_reduction(obj, reds[[tag]]$reduction, reds[[tag]]$assay,
                                  tag, cfg, plot_dir, data_dir, pal)
  invisible(obj)
}

# =============================================================================
## STAGE: canon — canonical PBMC RNA + ADT markers for a chosen pre-SCT cluster
# =============================================================================
# Characterises one cluster from the pre-SCT (LogNormalize / pca_lognorm) sweep
# with well-known PBMC lineage markers, from BOTH the RNA assay and the surface-
# protein ADT assay (preserved through the inject merge above). Shown on the
# pre-SCT UMAP and as per-cluster violins so the target cluster — where the
# injected sample is sparse — is read against every other cluster.
# Canonical RNA markers for cell-type calling (shown 1-column, violin per row).
CANONICAL_RNA_MARKERS <- c(
  "IL7R",    # Naive/Memory CD4+ T
  "CCR7",    # Naive T
  "ICOS",    # CD8 T
  "CXCR3",   # Memory CD4+ T
  "CD8A",    # CD8+ T
  "ROR1",    # RTK (Wnt5a receptor)
  "ROR2",    # RTK (Wnt5a receptor)
  "RORC",    # RORgammat (Th17 / type-3)
  "IFIT2",   # Interferon response
  "NKG7",    # NK
  "MS4A1",   # B cells
  "CD14",    # CD14+ Monocyte
  "FCGR3A",  # FCGR3A+ Monocyte
  "CD1C",    # cDC
  "MZB1"     # Plasma B (if present)
)

# Canonical surface (ADT) markers: ADT-assay rowname -> display name. Intersected
# with rownames(obj[["ADT"]]) at run time; unmatched names are reported/dropped.
ADT_MARKERS <- c(
  "Hu.CD3-UCHT1"  = "CD3",
  "Hu.CD4-RPA.T4" = "CD4",
  "Hu.CD8"        = "CD8",
  "Hu.CD45RA"     = "CD45RA",
  "Hu.CD45RO"     = "CD45RO",
  "Hu.CD335"      = "NKp46",
  "Hu.CD20-2H7"   = "CD20",
  "Hu.CD14-M5E2"  = "CD14",
  "Hu.CD16"       = "CD16",
  "Hu.CD1c"       = "CD1c",
  "Hu.CD303"      = "CD303"
)

stage_canon <- function(obj, cfg, plot_dir, target = "4") {
  message("\n[canon] Canonical RNA + ADT markers for pre-SCT cluster ", target, "...")

  # The pre-SCT sweep leaves one cluster column per (dims, res) combo and a
  # single umap_sweep_presct holding the last dims' embedding — target both.
  d   <- cfg$dims[length(cfg$dims)]
  res <- cfg$resolutions[length(cfg$resolutions)]
  clcol <- paste0("clust_presct_d", d, "_res", res)
  reduc <- "umap_sweep_presct"

  if (!clcol %in% colnames(obj@meta.data)) {
    message("  cluster column '", clcol, "' not found — run the 'cluster' stage ",
            "on the pre-SCT reduction in the same invocation; skipping.")
    return(invisible(obj))
  }
  if (!reduc %in% Reductions(obj)) {
    message("  reduction '", reduc, "' absent — skipping."); return(invisible(obj))
  }

  cl        <- as.character(obj@meta.data[[clcol]])
  num_ok    <- !anyNA(suppressWarnings(as.integer(unique(cl))))
  cl_levels <- if (num_ok) as.character(sort(unique(as.integer(cl)))) else sort(unique(cl))
  obj@meta.data[[clcol]] <- factor(cl, levels = cl_levels)
  Idents(obj) <- obj@meta.data[[clcol]]

  if (!target %in% cl_levels) {
    message("  cluster ", target, " not present in ", clcol, " (clusters: ",
            paste(cl_levels, collapse = ", "), ") — skipping.")
    return(invisible(obj))
  }

  inj <- if (length(cfg$inject_samples))
    as.character(obj@meta.data[[cfg$sample_col]]) %in% cfg$inject_samples else
    rep(FALSE, ncol(obj))
  n4     <- sum(cl == target)
  n4_inj <- sum(cl == target & inj)
  message("  cluster ", target, ": ", n4, " cells (injected ",
          paste(cfg$inject_samples, collapse = "/"), ": ", n4_inj, ")")

  # ── Where cluster `target` (+ injected cells) sit on the pre-SCT UMAP ─────────
  hi <- list(); hi[[paste0("cluster ", target)]] <- colnames(obj)[cl == target]
  if (any(inj)) hi[["injected"]] <- colnames(obj)[inj]
  p_loc <- DimPlot(obj, reduction = reduc, cells.highlight = hi,
                   cols.highlight = rev(c("#E31A1C", "#1F78B4")[seq_along(hi)]),
                   sizes.highlight = 0.6, pt.size = 0.3) +
    labs(title = paste0("pre-SCT cluster ", target, " (n=", n4,
                        ") + injected cells on ", reduc)) + theme_cowplot(12)
  save_png(p_loc, paste0("canon_cluster", target, "_location_umap"), plot_dir, 9, 7)

  # ── RNA canonical markers (LogNormalize / RNA assay) ─────────────────────────
  DefaultAssay(obj) <- "RNA"
  rna_present <- intersect(CANONICAL_RNA_MARKERS, rownames(obj))
  message("  RNA canonical present: ", length(rna_present), "/",
          length(CANONICAL_RNA_MARKERS))
  if (length(rna_present) > 0) {
    # FeaturePlots: max 4 panels (2x2) per page, larger dots.
    rna_fp_list <- FeaturePlot(obj, features = rna_present, reduction = reduc,
                               order = TRUE, pt.size = 1.0, raster = TRUE,
                               combine = FALSE)
    rna_fp_list <- lapply(rna_fp_list, function(p) p + theme_cowplot(11))
    rna_pages <- split(seq_along(rna_fp_list),
                       ceiling(seq_along(rna_fp_list) / 4))
    for (pg in seq_along(rna_pages)) {
      p_rna_fp <- wrap_plots(rna_fp_list[rna_pages[[pg]]], ncol = 2) +
        plot_annotation(title = paste0("Canonical PBMC RNA markers (", reduc,
                                       ") — page ", pg, "/", length(rna_pages)))
      save_png(p_rna_fp,
               paste0("canon_cluster", target, "_RNA_featureplot_p", pg),
               plot_dir, 12, 11)
    }

    p_rna_vln <- VlnPlot(obj, features = rna_present, group.by = clcol,
                         assay = "RNA", pt.size = 0, stack = TRUE, flip = TRUE) +
      NoLegend() +
      ggtitle(paste0("Canonical RNA markers by pre-SCT cluster (cluster ",
                     target, " of interest)")) +
      theme(axis.title.x = element_blank())
    save_png(p_rna_vln, paste0("canon_cluster", target, "_RNA_violin"),
             plot_dir, max(10, length(cl_levels) * 0.7),
             max(8, length(rna_present) * 0.35))
  }

  # ── ADT surface markers (carried through the inject merge) ───────────────────
  if (!"ADT" %in% Assays(obj)) {
    message("  No ADT assay in merged object — skipping surface-protein plots.")
    return(invisible(obj))
  }
  DefaultAssay(obj) <- "ADT"
  adt_names   <- names(ADT_MARKERS)
  adt_present <- intersect(adt_names, rownames(obj[["ADT"]]))
  adt_missing <- setdiff(adt_names, rownames(obj[["ADT"]]))
  message("  ADT canonical present: ", length(adt_present), "/", length(adt_names),
          if (length(adt_missing))
            paste0("  (absent: ", paste(adt_missing, collapse = ", "), ")") else "")
  if (length(adt_present) == 0) {
    message("  No canonical ADT names matched. Available ADT features: ",
            paste(head(rownames(obj[["ADT"]]), 40), collapse = ", "))
    return(invisible(obj))
  }
  adt_disp <- ADT_MARKERS[adt_present]   # display names, aligned to adt_present

  # FeaturePlot — one panel per antibody, CLR colour scale (margin=2, limits 0–3
  # squished) as in NPSLE_Tcell_Subcluster_Survey.R, titled with display names.
  # max 4 panels (2x2) per page, larger dots.
  fp_list <- FeaturePlot(obj, features = adt_present, reduction = reduc,
                         order = TRUE, pt.size = 1.0, raster = TRUE,
                         combine = FALSE)
  fp_list <- Map(function(p, disp) p +
    scale_color_gradientn(colors = c("lightgrey", "steelblue2", "navy"),
                          limits = c(0, 3), oob = scales::squish,
                          name = "CLR\n(m=2)") +
    ggtitle(paste0("ADT: ", disp)) + theme_cowplot(11),
    fp_list, adt_disp)
  adt_pages <- split(seq_along(fp_list), ceiling(seq_along(fp_list) / 4))
  for (pg in seq_along(adt_pages)) {
    p_adt_fp <- wrap_plots(fp_list[adt_pages[[pg]]], ncol = 2) +
      plot_annotation(title = paste0("Canonical surface markers (ADT, ", reduc,
                                     ") — page ", pg, "/", length(adt_pages)))
    save_png(p_adt_fp,
             paste0("canon_cluster", target, "_ADT_featureplot_p", pg),
             plot_dir, 12, 11)
  }

  # Violin by pre-SCT cluster — one row per antibody (1-column), display titles.
  vln_list <- VlnPlot(obj, features = adt_present, group.by = clcol,
                      assay = "ADT", pt.size = 0, combine = FALSE)
  vln_list <- Map(function(p, disp) p + NoLegend() +
    ggtitle(paste0("ADT: ", disp)) +
    theme(axis.title.x = element_blank(),
          plot.title  = element_text(size = 10, face = "bold")),
    vln_list, adt_disp)
  p_adt_vln <- wrap_plots(vln_list, ncol = 1) +
    plot_annotation(title = paste0("Canonical surface (ADT) markers by pre-SCT cluster",
                                   " (cluster ", target, " of interest)"))
  save_png(p_adt_vln, paste0("canon_cluster", target, "_ADT_violin"),
           plot_dir, max(10, length(cl_levels) * 0.7),
           max(8, length(adt_present) * 2.0))

  invisible(obj)
}

# =============================================================================
## DRIVER
# =============================================================================
if (!is.null(opt$from_checkpoint)) {
  message("\nLoading post-SCT checkpoint: ", opt$from_checkpoint)
  obj <- readRDS(opt$from_checkpoint)
} else {
  message("\nLoading input object: ", cfg$input_rds)
  obj <- readRDS(cfg$input_rds)
}

if (!cfg$sample_col %in% colnames(obj@meta.data))
  stop("sample_col '", cfg$sample_col, "' not found in meta.data. Available: ",
       paste(head(colnames(obj@meta.data), 40), collapse = ", "))

# Optional cell injection from a second .rds (batch-effect spot check).
# Pulls cells whose --sample_col value is in --inject_samples and merges them
# into obj before any normalisation/PCA, so downstream stages treat them as
# part of the same object. No anchor-based integration — straight Seurat merge.
if (!is.null(cfg$inject_rds) && length(cfg$inject_samples) > 0) {
  message("\n[inject] Loading ", cfg$inject_rds,
          " to extract sample(s): ", paste(cfg$inject_samples, collapse = ", "))
  inj_full <- readRDS(cfg$inject_rds)
  if (!cfg$sample_col %in% colnames(inj_full@meta.data))
    stop("[inject] sample_col '", cfg$sample_col, "' missing in inject_rds metadata.")
  inj_mask  <- as.character(inj_full@meta.data[[cfg$sample_col]]) %in% cfg$inject_samples
  inj_cells <- rownames(inj_full@meta.data)[inj_mask]
  if (length(inj_cells) == 0)
    stop("[inject] No cells matched ", paste(cfg$inject_samples, collapse = ","),
         " on column '", cfg$sample_col, "' in inject_rds.")
  inj_sub <- subset(inj_full, cells = inj_cells)
  message("  Extracted ", length(inj_cells), " cells from inject_rds (",
          paste(names(table(inj_sub@meta.data[[cfg$sample_col]])),
                table(inj_sub@meta.data[[cfg$sample_col]]),
                sep = "=", collapse = " | "), ")")

  # Align to the RNA assay on both objects and join layers (Seurat v5) so the
  # merge produces a clean combined counts matrix without split layers.
  DefaultAssay(obj)     <- "RNA"
  DefaultAssay(inj_sub) <- "RNA"
  tryCatch(obj     <- JoinLayers(obj,     assay = "RNA"), error = function(e) NULL)
  tryCatch(inj_sub <- JoinLayers(inj_sub, assay = "RNA"), error = function(e) NULL)

  # Carry the surface-protein (ADT) assay through the merge (merge() would
  # otherwise drop an assay not common to both objects). Align first, then
  # JoinLayers ADT on each side so the post-merge assay is single-layer.
  adt_al  <- align_adt_for_merge(obj, inj_sub)
  obj     <- adt_al$a; inj_sub <- adt_al$b
  if (adt_al$n_feat > 0) {
    tryCatch(obj     <- JoinLayers(obj,     assay = "ADT"), error = function(e) NULL)
    tryCatch(inj_sub <- JoinLayers(inj_sub, assay = "ADT"), error = function(e) NULL)
    if (adt_al$both)
      message("  ADT preserved across merge: ", adt_al$n_feat, " shared antibodies")
    else
      message("  ADT on one object only; other side zero-filled (", adt_al$n_feat,
              " antibodies) — cells without a real measurement read as 0.")
  } else {
    message("  No ADT assay to preserve across merge.")
  }

  n_before <- ncol(obj)
  obj <- merge(obj, y = inj_sub)
  tryCatch(obj <- JoinLayers(obj, assay = "RNA"), error = function(e) NULL)
  if ("ADT" %in% Assays(obj)) {
    tryCatch(obj <- JoinLayers(obj, assay = "ADT"), error = function(e) NULL)
    # CLR (Seurat default for ADT, margin = 2) so a 'data' layer exists to plot.
    obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                         margin = 2, verbose = FALSE)
    DefaultAssay(obj) <- "RNA"
  }
  rm(inj_full, inj_sub); invisible(gc(verbose = FALSE))
  message("  Merged: ", n_before, " -> ", ncol(obj), " cells total")
}

pal <- make_sample_pal(obj@meta.data[[cfg$sample_col]])
message("  samples (", cfg$sample_col, "): ",
        paste(names(pal), collapse = ", "))

if (run_stage("presct")) obj <- stage_presct(obj, cfg, plot_dir)
if (run_stage("sct"))    obj <- stage_sct(obj, cfg, plot_dir, data_dir, ckpt_sct)
if (run_stage("pcviz"))  obj <- stage_pcviz(obj, cfg, plot_dir, data_dir, pal)
if (run_stage("cluster")) obj <- stage_cluster(obj, cfg, plot_dir, data_dir, pal)
if (run_stage("canon"))   obj <- stage_canon(obj, cfg, plot_dir)

message("\n== Done. Outputs in: ", out_root, " ==")
