#!/usr/bin/env Rscript
# =============================================================================
#  NPSLE T Cell — SCTransform v2 Normalization Pipeline
#
#  Input : NPSLE_Tcell_Subcluster.rds  (T cells already subset + labelled)
#  Output: SCT_v2_Output/
#
#  Normalization: SCTransform v2 (vst.flavor = "v2", glmGamPoi backend)
#    • Replaces NormalizeData → FindVariableFeatures → ScaleData
#    • Regresses: percent.mt, S.Score, G2M.Score
#    • nCount_RNA is handled implicitly by the GLM (NOT in vars.to.regress)
#  DE testing  : PrepSCTFindMarkers() → FindAllMarkers(assay="SCT")
#  Pseudobulk  : DESeq2 on raw RNA counts (RNA assay, not SCT)
# =============================================================================

.libPaths(c("/path/to/your/R_library", .libPaths()))

#BiocManager::install("glmGamPoi")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(sctransform)      # >= 0.3.3 for vst.flavor = "v2"
  library(glmGamPoi)        # fast GLM backend for SCTransform v2
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(DESeq2)
  library(speckle)
  library(limma)
  library(clustree)
  library(colorspace)
  library(lme4)
  library(optparse)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input_rds",  type = "character", default = NULL),
  make_option("--out_dir",    type = "character", default = NULL),
  make_option("--final_res",  type = "double",    default = NULL),
  make_option("--pca_dims",   type = "integer",   default = NULL),
  make_option("--n_hvg",      type = "integer",   default = NULL),
  make_option("--scan_dims",  type = "character", default = NULL,
              help = "Comma-separated dims to scan, e.g. '10,15,20,25,30'"),
  make_option("--scan_res",   type = "character", default = NULL,
              help = "Comma-separated resolutions to scan, e.g. '0.2,0.4,0.6,0.8'"),
  make_option("--stages",     type = "character", default = NULL,
              help = paste("Comma-separated stages to run (or 'all'):",
                           "qc, sct, pca, scan, cluster, umap, markers, de"))
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Configuration ─────────────────────────────────────────────────────────────
cfg <- list(

  input_rds = file.path(
    "./NPSLE_Tcell_Subcluster.rds"),

  out_dir = file.path(
    "Tcell_SCTv2_Output"),

  # Metadata columns
  celltype_col = "celltype",
  donor_col    = "HTO_maxID",
  group_col    = "Status",
  sample_col   = "sample_origin",
  hi_group     = "Hi_Cog",
  lo_group     = "Lo_Cog",

  pool_cols = c("A" = "#4DAFE0", "B" = "#E0874D"),

  # SCTransform v2 settings
  n_hvg            = 3000,          # SCTransform v2 default; more stable than 2000
  exclude_patterns = c("^MT-", "^RPS", "^RPL", "^IGHV", "^IGLV", "^IGKV",
                       "^TRAV", "^TRBV", "^TRDV", "^TRGV"),
  # nCount_RNA is NOT included — SCT models it internally via GLM offset
  regress_vars = c("percent.mt", "S.Score", "G2M.Score"),

  # PCA / UMAP
  n_pcs    = 30,
  pca_dims = 1:14,
  final_res = 0.8,

  # DimRezScanner sweep (adapted from SingleCell_RDS_list_DimRezScanner_UMAP.R)
  scan_dims = c(14, 16, 18, 20),
  scan_res  = c(0.4, 0.6, 0.8),

  # Resolution sweep for clustree
  resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),

  # Figure dimensions
  umap_width  = 10, umap_height = 8,
  dot_width   = 14, dot_height  = 9,

  # ── Stage control ────────────────────────────────────────────────────────────
  # Which pipeline stages to run. Set to "all" or a character vector of names:
  #   "qc"      — QC violin plots (by group + by sample)
  #   "sct"     — Cell cycle scoring + SCTransform v2
  #   "pca"     — RunPCA + elbow plot
  #   "scan"    — DimRezScanner UMAP sweep (dims × resolution grid)
  #   "cluster" — Final RunUMAP + FindNeighbors + FindClusters + clustree
  #   "umap"    — UMAP visualization panels
  #   "markers" — PrepSCTFindMarkers + FindAllMarkers + dotplot + heatmap
  #   "de"      — Propeller + MASC + pseudobulk DESeq2
  #
  # For steps 1-7 only:  stages = c("qc","sct","pca","scan","cluster")
  # For downstream only: stages = c("umap","markers","de")   (loads cluster ckpt)
  stages = "all"
)

# CLI overrides
if (!is.null(opt$input_rds)) cfg$input_rds  <- opt$input_rds
if (!is.null(opt$out_dir))   cfg$out_dir    <- opt$out_dir
if (!is.null(opt$final_res)) cfg$final_res  <- opt$final_res
if (!is.null(opt$pca_dims))  cfg$pca_dims   <- seq_len(opt$pca_dims)
if (!is.null(opt$n_hvg))     cfg$n_hvg      <- opt$n_hvg
if (!is.null(opt$scan_dims)) {
  cfg$scan_dims <- as.integer(trimws(strsplit(opt$scan_dims, ",")[[1]]))
}
if (!is.null(opt$scan_res)) {
  cfg$scan_res <- as.numeric(trimws(strsplit(opt$scan_res, ",")[[1]]))
}
if (!is.null(opt$stages)) {
  cfg$stages <- trimws(strsplit(opt$stages, ",")[[1]])
}

# Checkpoint paths — derived from out_dir so they survive across sessions
cfg$ckpt_sct     <- file.path(cfg$out_dir, "NPSLE_Tcell_SCTv2_post_sct.rds")
cfg$ckpt_pca     <- file.path(cfg$out_dir, "NPSLE_Tcell_SCTv2_post_pca.rds")
cfg$ckpt_cluster <- file.path(cfg$out_dir, "NPSLE_Tcell_SCTv2_clustered.rds")

active_stages <- if ("all" %in% cfg$stages) {
  c("qc","sct","pca","scan","cluster","umap","markers","de")
} else {
  cfg$stages
}

message("Pipeline config:")
message("  input_rds : ", cfg$input_rds)
message("  out_dir   : ", cfg$out_dir)
message("  stages    : ", paste(active_stages, collapse = ", "))
message("  n_hvg     : ", cfg$n_hvg, "  (SCTransform v2)")
message("  pca_dims  : 1–", max(cfg$pca_dims))
message("  final_res : ", cfg$final_res)
message("  scan_dims : ", paste(cfg$scan_dims, collapse = ", "))
message("  scan_res  : ", paste(cfg$scan_res,  collapse = ", "))

# ── Canonical T cell marker panels ────────────────────────────────────────────
TCELL_MARKERS <- list(
  T_lineage      = c("CD3D", "CD3E", "TRAC"),
  CD4            = c("CD4", "IL7R"),
  CD8            = c("CD8A", "CD8B"),
  Naive          = c("CCR7", "SELL", "TCF7", "LEF1", "KLF2"),
  Tcm            = c("IL7R", "FN1", "ANXA1"),
  Tem            = c("GZMK", "EOMES", "PRDM1"),
  Temra          = c("GZMB", "GNLY", "PRF1", "CX3CR1", "FGFBP2", "KLRG1"),
  Treg           = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  Tfh            = c("CXCR5", "ICOS", "BCL6", "SH2D1A"),
  Th1            = c("TBX21", "CXCR3", "IFNG"),
  Th17           = c("RORC", "CCR6", "IL23R"),
  Exhausted      = c("LAG3", "PDCD1", "HAVCR2", "TIGIT", "TOX", "ENTPD1"),
  IFN_stimulated = c("ISG15", "MX1", "IFIT1", "IFIT3", "RSAD2"),
  Cycling        = c("MKI67", "TOP2A", "PCNA"),
  NKT_like       = c("NKG7", "GNLY", "NCAM1"),
  MAIT           = c("SLC4A10", "KLRB1", "NCR3"),
  GammaDelta     = c("TRDC", "TRGC1", "TRGC2")
)

# T cell ADT panel (from SingleCell_RDS_list_DimRezScanner_UMAP.R)
TCELL_ADT_MARKERS <- c(
  "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8",
  "Hu.CD2",       "Hu.CD5",         "Hu.CD7",
  "Hu.TCR.AB",    "Hu.TCR.Vd2",     "Hu.TCR.Va7.2",
  "Hu.CD45RA",  "Hu.CD45RO",  "Hu.CD62L",
  "Hu.CD27",    "Hu.CD28",    "Hu.CD95",
  "Hu.CD127",   "Hu.CD122",   "HuMs.CD44",
  "Hu.CD25",
  "Hu.CD57",    "Hu.CD56",    "Hu.CD94",
  "Hu.CX3CR1",  "Hu.KLRG1",  "Hu.CD161",
  "Hu.CD69",   "Hu.CD38-HIT2", "Hu.HLA.DR",
  "Hu.CD26",   "Hu.CD137",     "Hu.CD154",
  "Hu.CD223",  "Hu.CD279",  "Hu.TIGIT",
  "Hu.CD152",  "Hu.CD244",
  "Hu.CD183",  "Hu.CD185",  "Hu.CD194",
  "Hu.CD195",  "Hu.CD196",
  "Hu.CD49d",  "Hu.CD58",   "Hu.CD52"
)

# =============================================================================
##  HELPER FUNCTIONS
# =============================================================================

save_figure <- function(p, name, dir, width, height, units = "in", dpi = 300) {
  png_path <- file.path(dir, paste0(name, ".png"))
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = width, height = height, units = units, dpi = dpi)
  ggsave(pdf_path, p, width = width, height = height, units = units)
  message("  Saved: ", basename(png_path))
  invisible(list(png = png_path, pdf = pdf_path))
}

filter_hvg <- function(var_features, exclude_patterns) {
  keep <- rep(TRUE, length(var_features))
  for (pat in exclude_patterns) keep <- keep & !grepl(pat, var_features)
  var_features[keep]
}

add_mito_pct <- function(obj, assay = "RNA") {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt_genes <- rownames(obj)[grepl("^MT-", rownames(obj))]
    if (length(mt_genes) == 0) {
      obj$percent.mt <- 0
    } else {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    }
  }
  obj
}

# ── QC violin plots: by group (Status) ───────────────────────────────────────
plot_qc_by_group <- function(obj, group_col, title = "QC", ncol = 2) {
  qc_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_features) == 0) {
    warning("No QC fields found.")
    return(NULL)
  }
  plist <- VlnPlot(obj, features = qc_features, group.by = group_col,
                   pt.size = 0, combine = FALSE)
  n <- length(plist)
  for (i in seq_along(plist)) {
    plist[[i]] <- plist[[i]] +
      ggtitle(qc_features[i]) +
      theme(
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 22),
        axis.title.x     = element_blank(),
        axis.text.x      = if (i < n) element_blank()
                           else element_text(angle = 45, hjust = 1, size = 19),
        axis.ticks.x     = if (i < n) element_blank() else element_line()
      )
  }
  wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))
}

# ── QC violin plots: by sample, coloured by pool ─────────────────────────────
plot_qc_by_sample_pooled <- function(obj, cfg, title = "Single-cell QC by Sample") {
  qc_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_features) == 0) return(NULL)
  sample_col <- cfg$sample_col

  md <- data.frame(Sample = as.character(obj@meta.data[[sample_col]]),
                   stringsAsFactors = FALSE)
  md$PoolGroup <- sub(".*_([^_]+)$", "\\1", md$Sample)
  md$PoolGroup[!md$PoolGroup %in% names(cfg$pool_cols)] <- "Unknown"
  md$PoolGroup <- factor(md$PoolGroup, levels = c(names(cfg$pool_cols), "Unknown"))

  pool_df <- md %>% distinct(Sample, PoolGroup) %>% arrange(PoolGroup, Sample)
  sample_order <- pool_df$Sample
  sample_cols  <- setNames(
    ifelse(as.character(pool_df$PoolGroup) %in% names(cfg$pool_cols),
           cfg$pool_cols[as.character(pool_df$PoolGroup)], "grey60"),
    pool_df$Sample
  )
  obj@meta.data[[sample_col]] <- factor(
    as.character(obj@meta.data[[sample_col]]), levels = sample_order
  )

  plist <- VlnPlot(obj, features = qc_features, group.by = sample_col,
                   pt.size = 0, combine = FALSE, cols = sample_cols)
  n <- length(plist)
  for (i in seq_along(plist)) {
    plist[[i]] <- plist[[i]] +
      ggtitle(qc_features[i]) +
      theme(
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 22),
        axis.title.x     = element_blank(),
        axis.text.x      = if (i < n) element_blank()
                           else element_text(angle = 45, hjust = 1, size = 19),
        axis.ticks.x     = if (i < n) element_blank() else element_line()
      )
  }

  legend_df <- data.frame(x = 1, y = 1,
                           PoolGroup = factor(names(cfg$pool_cols),
                                              levels = names(cfg$pool_cols)))
  p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = PoolGroup)) +
    geom_col() +
    scale_fill_manual(values = cfg$pool_cols, name = "Pool Group") +
    theme_void() + theme(legend.position = "right")
  legend_grob <- cowplot::get_legend(p_legend)

  combined <- wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))

  combined + inset_element(legend_grob, left = 0.85, bottom = 0.02,
                            right = 1.0, top = 0.12)
}

# ── SCT residuals QC: nCount_SCT and nFeature_SCT on UMAP ────────────────────
# After SCTransform, the SCT assay exposes corrected counts and residuals.
# Plotting nCount_SCT / nFeature_SCT lets you verify that sequencing-depth
# variation is no longer a major driver of cell separation in the UMAP.
plot_sct_qc_umaps <- function(obj, plot_dir) {
  qc_scale <- scale_color_gradientn(
    colours = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C",
                "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"),
    na.value = "grey90"
  )

  sct_qc_feats <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
      "nCount_SCT", "nFeature_SCT"),
    colnames(obj@meta.data)
  )
  qc_titles <- c(
    nCount_RNA   = "UMI Count (RNA)",
    nFeature_RNA = "Genes Detected (RNA)",
    percent.mt   = "% Mitochondrial",
    percent.ribo = "% Ribosomal",
    nCount_SCT   = "Corrected Counts (SCT)",
    nFeature_SCT = "Genes Detected (SCT)"
  )
  if (length(sct_qc_feats) == 0) return(invisible(NULL))

  qc_plts <- lapply(sct_qc_feats, function(feat) {
    FeaturePlot(obj, features = feat, reduction = "umap",
                order = TRUE, pt.size = 0.3) +
      qc_scale +
      ggtitle(qc_titles[feat]) +
      theme_cowplot(12) +
      theme(
        legend.position   = "right",
        plot.title        = element_text(size = 13, face = "bold"),
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.35, "cm")
      )
  })

  ncols_qc <- min(3L, length(qc_plts))
  nrows_qc <- ceiling(length(qc_plts) / ncols_qc)
  p <- wrap_plots(qc_plts, ncol = ncols_qc) +
    plot_annotation(
      title = "QC Metrics on UMAP (post-SCTransform v2)",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  save_figure(p, "08b_umap_qc_metrics_SCT", plot_dir,
              width  = ncols_qc * 6,
              height = nrows_qc * 5.5)
}

# ── Paginated panel saver (from SingleCell_RDS_list_DimRezScanner_UMAP.R) ────
save_panel_pages <- function(panel_list, base_file, ncol,
                              panel_w, panel_h, dpi,
                              title = NULL, subtitle = NULL,
                              max_per_page = 12) {
  chunks    <- split(panel_list, ceiling(seq_along(panel_list) / max_per_page))
  out_files <- character(length(chunks))
  for (pg in seq_along(chunks)) {
    ch     <- chunks[[pg]]
    n_cols <- min(ncol, length(ch))
    n_rows <- ceiling(length(ch) / n_cols)
    pg_ttl <- if (!is.null(title) && length(chunks) > 1)
      paste0(title, "  (page ", pg, "/", length(chunks), ")") else title
    p <- wrap_plots(ch, ncol = n_cols) +
      plot_annotation(title = pg_ttl, subtitle = subtitle)
    suffix   <- if (length(chunks) > 1) paste0("_page", pg) else ""
    out_file <- sub("\\.jpeg$", paste0(suffix, ".jpeg"), base_file)
    ggsave(out_file, p,
           width  = n_cols * panel_w,
           height = n_rows * panel_h,
           dpi    = dpi)
    out_files[pg] <- out_file
  }
  out_files
}

# ── DimRezScanner: UMAP grid × dims × resolution (adapted from DimRezScanner) ─
# Sweeps all combinations of (dims, resolution) and saves:
#   - Individual DimPlot per combo (UMAP_dims{d}_res{r}.jpeg)
#   - Paginated DimPlot grid summary
#   - RNA FeaturePlot panels per combo
#   - ADT FeaturePlot panels per combo (if ADT assay present)
#   - VlnPlot by cluster per combo
run_dimrez_scanner <- function(
  obj,
  scan_dir,
  dims_vec    = c(10, 15, 20, 25, 30),
  resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0),
  features    = c("CD3D", "CD4", "CD8A", "CCR7", "SELL", "IL7R",
                  "GZMB", "GZMK", "FOXP3", "MKI67", "LAG3", "PDCD1",
                  "ISG15", "CX3CR1", "NKG7"),
  adt_markers = TCELL_ADT_MARKERS,
  dpi         = 150    # lower DPI for scanner grids to keep file sizes manageable
) {
  dir.create(scan_dir, showWarnings = FALSE, recursive = TRUE)

  # Clamp dims to available PCs
  max_pcs  <- ncol(obj@reductions$pca)
  bad_dims <- dims_vec[dims_vec > max_pcs]
  if (length(bad_dims)) {
    warning("Clamping dims > max PCs (", max_pcs, "): ",
            paste(bad_dims, collapse = ", "))
    dims_vec <- dims_vec[dims_vec <= max_pcs]
  }

  has_adt     <- "ADT" %in% names(obj@assays)
  adt_present <- if (has_adt) intersect(adt_markers, rownames(obj@assays[["ADT"]]))
                 else character(0)

  grid_df <- expand.grid(dims = dims_vec, resolution = resolutions,
                         stringsAsFactors = FALSE)
  panel_size <- 6   # inches per DimPlot panel

  jobs <- lapply(seq_len(nrow(grid_df)), function(i) {
    dims_i    <- grid_df$dims[i]
    res_i     <- grid_df$resolution[i]
    dims_use  <- seq_len(dims_i)
    res_tag   <- gsub("\\.", "p", as.character(res_i))
    tag       <- paste0("d", dims_i, "_r", res_tag)
    umap_name <- paste0("umap_", tag)
    snn_name  <- paste0("snn_",  tag)

    message("  [DimRezScanner] dims=1:", dims_i, "  res=", res_i)

    obj2 <- FindNeighbors(obj, reduction = "pca", dims = dims_use,
                          graph.name = snn_name, verbose = FALSE)
    obj2 <- FindClusters(obj2, graph.name = snn_name,
                          resolution = res_i, verbose = FALSE)
    obj2 <- RunUMAP(obj2, reduction = "pca", dims = dims_use,
                    reduction.name = umap_name, verbose = FALSE)

    n_clust   <- length(unique(Idents(obj2)))
    combo_lbl <- paste0("dims 1:", dims_i, "  |  res ", res_i,
                        "  |  ", n_clust, " clusters")

    # DimPlot
    p_dim <- DimPlot(obj2, reduction = umap_name, label = TRUE, repel = TRUE) +
      ggtitle(combo_lbl) +
      theme(plot.title = element_text(size = 11))
    ggsave(
      file.path(scan_dir, paste0("UMAP_dims", dims_i, "_res", res_i, ".jpeg")),
      p_dim, width = panel_size, height = panel_size, dpi = dpi
    )

    # RNA FeaturePlot
    genes_present <- intersect(features, rownames(obj2))
    if (length(genes_present) > 0) {
      DefaultAssay(obj2) <- "SCT"
      fp_list <- FeaturePlot(obj2, features = genes_present,
                             reduction = umap_name, pt.size = 2,
                             order = TRUE, combine = FALSE, raster = TRUE)
      fp_list <- lapply(seq_along(fp_list), function(k)
        fp_list[[k]] +
          ggtitle(genes_present[[k]]) +
          theme(plot.title = element_text(size = 11, face = "bold")))
      save_panel_pages(
        fp_list,
        base_file = file.path(scan_dir,
          paste0("FeaturePlots_RNA_dims", dims_i, "_res", res_i, ".jpeg")),
        ncol = 3, panel_w = 5.5, panel_h = 5.0, dpi = dpi,
        title = paste0("RNA Feature Plots  |  ", combo_lbl)
      )
    }

    # ADT FeaturePlot
    if (has_adt && length(adt_present) > 0) {
      DefaultAssay(obj2) <- "ADT"
      adt_fp_list <- FeaturePlot(obj2, features = adt_present,
                                  reduction = umap_name, pt.size = 2,
                                  order = TRUE, combine = FALSE, raster = TRUE)
      adt_display <- sub("^Hu[A-Za-z]*\\.", "", adt_present)
      adt_fp_list <- lapply(seq_along(adt_fp_list), function(k)
        adt_fp_list[[k]] +
          ggtitle(adt_display[[k]]) +
          theme(plot.title = element_text(size = 10, face = "bold")))
      save_panel_pages(
        adt_fp_list,
        base_file = file.path(scan_dir,
          paste0("FeaturePlots_ADT_dims", dims_i, "_res", res_i, ".jpeg")),
        ncol = 4, panel_w = 5.0, panel_h = 4.5, dpi = dpi,
        title = paste0("ADT Feature Plots  |  ", combo_lbl)
      )
      DefaultAssay(obj2) <- "SCT"
    }

    # VlnPlot by cluster
    if (length(genes_present) > 0) {
      DefaultAssay(obj2) <- "SCT"
      vln_list <- VlnPlot(obj2, features = genes_present,
                           group.by = "seurat_clusters",
                           pt.size = 0, combine = FALSE)
      vln_list <- lapply(seq_along(vln_list), function(k)
        vln_list[[k]] +
          ggtitle(genes_present[[k]]) +
          theme(
            legend.position = "none",
            plot.title      = element_text(size = 11, face = "bold"),
            axis.title.x    = element_blank(),
            axis.text.x     = element_text(angle = 45, hjust = 1, size = 8)
          ))
      save_panel_pages(
        vln_list,
        base_file = file.path(scan_dir,
          paste0("VlnPlot_clusters_dims", dims_i, "_res", res_i, ".jpeg")),
        ncol = 3, panel_w = 5.0, panel_h = 3.5, dpi = dpi,
        title = paste0("RNA Violin by Cluster  |  ", combo_lbl)
      )
    }

    list(dims = dims_i, resolution = res_i, tag = tag, plot = p_dim)
  })

  # Paginated DimPlot summary grid
  plot_lookup <- setNames(lapply(jobs, `[[`, "plot"),
                          vapply(jobs, `[[`, character(1), "tag"))
  ordered <- list()
  for (d in dims_vec)
    for (r in resolutions) {
      key <- paste0("d", d, "_r", gsub("\\.", "p", as.character(r)))
      if (!is.null(plot_lookup[[key]])) ordered[[key]] <- plot_lookup[[key]]
    }

  page_ncol <- min(3L, length(resolutions))
  chunks    <- split(ordered, ceiling(seq_along(ordered) / 9))
  for (pg in seq_along(chunks)) {
    ch     <- chunks[[pg]]
    n_cols <- min(page_ncol, length(ch))
    n_rows <- ceiling(length(ch) / n_cols)
    suffix <- if (length(chunks) > 1) paste0("_page", pg) else ""
    pg_lbl <- if (length(chunks) > 1)
      paste0("  (page ", pg, "/", length(chunks), ")") else ""
    pg_plot <- wrap_plots(ch, ncol = n_cols) +
      plot_annotation(
        title    = paste0("SCT v2 UMAP Parameter Scan", pg_lbl),
        subtitle = paste0("Rows = dims (", paste(dims_vec, collapse = ", "),
                          ")  |  Cols = res (", paste(resolutions, collapse = ", "), ")")
      )
    ggsave(file.path(scan_dir, paste0("UMAP_grid", suffix, ".jpeg")),
           pg_plot, width = n_cols * panel_size, height = n_rows * panel_size,
           dpi = dpi)
    message("  DimRezScanner grid page ", pg, " saved")
  }

  invisible(jobs)
}

# ── Pseudobulk DESeq2 ─────────────────────────────────────────────────────────
# Uses raw RNA counts — NOT the SCT corrected counts.
# SCT residuals are designed for dimensionality reduction and clustering,
# not for differential expression at the pseudobulk level. DESeq2 requires
# raw integer count data so its negative binomial dispersion model is valid.
pseudobulk_tcell_deseq2 <- function(
  obj,
  cluster_col = "seurat_clusters",
  donor_col   = "HTO_maxID",
  group_col   = "Status",
  hi_group    = "Hi_Cog",
  lo_group    = "Lo_Cog",
  min_cells   = 5,
  min_counts  = 10
) {
  stopifnot(inherits(obj, "Seurat"))
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot  = "counts")
  )
  counts <- as(counts, "dgCMatrix")

  md <- obj@meta.data[colnames(counts), , drop = FALSE]
  md$cluster <- as.character(md[[cluster_col]])
  md$donor   <- as.character(md[[donor_col]])
  md$group   <- as.character(md[[group_col]])
  md <- md[md$group %in% c(hi_group, lo_group), , drop = FALSE]
  counts <- counts[, rownames(md), drop = FALSE]

  results <- list()
  for (cl in sort(unique(md$cluster))) {
    cl_md <- md[md$cluster == cl, , drop = FALSE]
    cl_ct <- counts[, rownames(cl_md), drop = FALSE]

    pb_list <- lapply(unique(cl_md$donor), function(d) {
      idx <- rownames(cl_md)[cl_md$donor == d]
      if (length(idx) < min_cells) return(NULL)
      Matrix::rowSums(cl_ct[, idx, drop = FALSE])
    })
    names(pb_list) <- unique(cl_md$donor)
    pb_list <- Filter(Negate(is.null), pb_list)

    donor_groups <- cl_md %>%
      distinct(donor, group) %>%
      filter(donor %in% names(pb_list))
    has_hi <- any(donor_groups$group == hi_group)
    has_lo <- any(donor_groups$group == lo_group)
    if (!has_hi || !has_lo || nrow(donor_groups) < 3) {
      message("  Cluster ", cl, ": skipped (need ≥3 donors with both groups)")
      next
    }

    pb_mat <- do.call(cbind, pb_list)
    pb_mat <- pb_mat[Matrix::rowSums(pb_mat) >= min_counts, , drop = FALSE]

    col_d <- donor_groups[match(colnames(pb_mat), donor_groups$donor), ]
    col_d$group <- factor(col_d$group, levels = c(lo_group, hi_group))
    rownames(col_d) <- col_d$donor

    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(pb_mat)),
        colData   = col_d,
        design    = ~ group
      )
      dds <- DESeq(dds, quiet = TRUE)
      res <- results(dds, contrast = c("group", hi_group, lo_group))
      res_df <- as.data.frame(res)
      res_df$gene    <- rownames(res_df)
      res_df$cluster <- cl
      results[[cl]]  <- res_df[order(res_df$padj), ]
      message("  DESeq2 [", cl, "]: ",
              sum(!is.na(res_df$padj) & res_df$padj < 0.05),
              " sig genes (padj<0.05)")
    }, error = function(e)
      message("  DESeq2 failed cluster ", cl, ": ", conditionMessage(e)))
  }
  results
}

# ── Propeller ─────────────────────────────────────────────────────────────────
run_tcell_propeller <- function(
  obj,
  cluster_col = "seurat_clusters",
  donor_col   = "HTO_maxID",
  group_col   = "Status",
  hi_group    = "Hi_Cog",
  lo_group    = "Lo_Cog"
) {
  md   <- obj@meta.data
  keep <- md[[group_col]] %in% c(hi_group, lo_group)
  props <- speckle::getTransformedProps(
    clusters  = md[[cluster_col]][keep],
    sample    = md[[donor_col]][keep],
    transform = "logit"
  )
  P <- props$Proportions
  donor_grp <- md[keep, ] %>%
    distinct(.data[[donor_col]], .data[[group_col]]) %>%
    rename(donor = 1, group = 2) %>%
    filter(donor %in% colnames(P))
  P <- P[, donor_grp$donor, drop = FALSE]
  donor_grp$group <- factor(donor_grp$group, levels = c(lo_group, hi_group))
  design   <- model.matrix(~ 0 + group, data = donor_grp)
  colnames(design) <- c(lo_group, hi_group)
  contrast <- limma::makeContrasts(
    contrasts = paste0(hi_group, "-", lo_group),
    levels    = design
  )
  speckle::propeller.ttest(props, design = design,
                            contrasts = contrast, robust = TRUE,
                            trend = FALSE, sort = TRUE)
}

# ── MASC ──────────────────────────────────────────────────────────────────────
run_masc <- function(
  obj,
  cluster_col = "seurat_clusters",
  donor_col   = "HTO_maxID",
  group_col   = "Status",
  hi_group    = "Hi_Cog",
  lo_group    = "Lo_Cog"
) {
  suppressPackageStartupMessages(library(lme4))
  md <- as.data.frame(obj@meta.data)
  md$donor   <- as.character(md[[donor_col]])
  md$group   <- factor(md[[group_col]], levels = c(lo_group, hi_group))
  md$cluster <- as.character(md[[cluster_col]])
  md <- md[md$group %in% c(hi_group, lo_group), , drop = FALSE]

  clusters <- sort(unique(md$cluster))
  results  <- lapply(clusters, function(cl) {
    md$is_cl <- as.integer(md$cluster == cl)
    fit <- tryCatch(
      glmer(is_cl ~ group + (1 | donor), data = md, family = binomial,
            control = glmerControl(optimizer = "bobyqa")),
      error = function(e) NULL
    )
    if (is.null(fit))
      return(data.frame(cluster = cl, estimate = NA, se = NA, z = NA, p_value = NA))
    cf  <- summary(fit)$coefficients
    row <- cf[grep(hi_group, rownames(cf)), , drop = FALSE]
    if (nrow(row) == 0)
      return(data.frame(cluster = cl, estimate = NA, se = NA, z = NA, p_value = NA))
    data.frame(cluster  = cl,
               estimate = row[1, "Estimate"],
               se       = row[1, "Std. Error"],
               z        = row[1, "z value"],
               p_value  = row[1, "Pr(>|z|)"])
  })
  out <- do.call(rbind, results)
  out$p_adj     <- p.adjust(out$p_value, method = "BH")
  out$direction <- ifelse(is.na(out$estimate), "NA",
                          ifelse(out$estimate > 0, hi_group, lo_group))
  out[order(out$p_value), ]
}

# =============================================================================
##  MAIN PIPELINE
# =============================================================================

run_sct_pipeline <- function(cfg) {

  # ── Stage helpers ──────────────────────────────────────────────────────────
  ALL_STAGES <- c("qc","sct","pca","scan","cluster","umap","markers","de")
  active     <- if ("all" %in% cfg$stages) ALL_STAGES else cfg$stages
  run_stage  <- function(s) s %in% active

  message("\n== Stages active: ", paste(active, collapse = ", "), " ==")

  # ── Directories ───────────────────────────────────────────────────────────
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(cfg$out_dir, "plots")
  data_dir <- file.path(cfg$out_dir, "data")
  scan_dir <- file.path(cfg$out_dir, "dimrez_scan")
  for (d in c(plot_dir, data_dir, scan_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # ── Load: pick starting checkpoint based on earliest active stage ──────────
  # Dependency chain: qc/sct → raw | pca → sct_ckpt | scan/cluster → pca_ckpt
  #                   umap/markers/de → cluster_ckpt
  start_from <- if (run_stage("sct")) {
    "raw"
  } else if (run_stage("pca")) {
    "sct_ckpt"
  } else if (any(c("scan","cluster") %in% active)) {
    "pca_ckpt"
  } else if (any(c("umap","markers","de") %in% active)) {
    "cluster_ckpt"
  } else {
    "raw"
  }

  obj <- switch(start_from,
    raw = {
      message("\n[load] Reading input: ", cfg$input_rds)
      tmp <- readRDS(cfg$input_rds)
      tmp <- add_mito_pct(tmp, assay = "RNA")
      DefaultAssay(tmp) <- "RNA"
      tmp
    },
    sct_ckpt = {
      message("\n[load] Resuming from SCT checkpoint: ", cfg$ckpt_sct)
      readRDS(cfg$ckpt_sct)
    },
    pca_ckpt = {
      message("\n[load] Resuming from PCA checkpoint: ", cfg$ckpt_pca)
      readRDS(cfg$ckpt_pca)
    },
    cluster_ckpt = {
      message("\n[load] Resuming from cluster checkpoint: ", cfg$ckpt_cluster)
      readRDS(cfg$ckpt_cluster)
    }
  )

  message("  Cells: ", ncol(obj))
  message("  Groups: ",
          paste(names(table(obj@meta.data[[cfg$group_col]])),
                table(obj@meta.data[[cfg$group_col]]), sep = "=", collapse = " | "))
  if (cfg$celltype_col %in% colnames(obj@meta.data))
    print(table(obj@meta.data[[cfg$celltype_col]]))

  # Initialize return values (populated only if respective stage runs)
  markers_sig <- NULL
  de_results  <- NULL
  masc_res    <- NULL
  prop_res    <- NULL

  # ── Stage: qc ─────────────────────────────────────────────────────────────
  if (run_stage("qc")) {
    message("\n[qc] QC plots...")

    p_qc_group <- plot_qc_by_group(obj, group_col = cfg$group_col,
                                    title = "T Cell QC by Cognitive Status")
    if (!is.null(p_qc_group)) {
      n_qc <- length(intersect(
        c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
        colnames(obj@meta.data)))
      save_figure(p_qc_group, "00a_qc_by_status", plot_dir,
                  width = 10, height = 4 * n_qc)
    }

    p_qc_sample <- tryCatch(
      plot_qc_by_sample_pooled(obj, cfg, title = "T Cell QC by Sample (Pool A vs B)"),
      error = function(e) { message("  QC by sample failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(p_qc_sample)) {
      n_qc <- length(intersect(
        c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
        colnames(obj@meta.data)))
      save_figure(p_qc_sample, "00b_qc_by_sample", plot_dir,
                  width = 18, height = 4 * n_qc)
    }
  }

  # ── Stage: sct ────────────────────────────────────────────────────────────
  # Cell cycle scoring must happen first so S.Score/G2M.Score are available
  # for vars.to.regress. A quick log-norm is run solely for CellCycleScoring;
  # SCTransform then re-normalizes from raw counts independently.
  if (run_stage("sct")) {
    message("\n[sct] Cell cycle scoring + SCTransform v2...")
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- CellCycleScoring(obj,
      s.features   = Seurat::cc.genes.updated.2019$s.genes,
      g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
      set.ident    = FALSE)
    print(table(obj$Phase))

    p_cc_violin <- VlnPlot(obj, features = c("S.Score","G2M.Score"),
                            group.by = cfg$group_col, pt.size = 0, combine = TRUE) +
      plot_annotation(title = "Cell Cycle Scores by Group")
    save_figure(p_cc_violin, "00c_cc_scores_by_group", plot_dir, 10, 5)

    message("  SCTransform v2 (vars.to.regress: ",
            paste(cfg$regress_vars, collapse = ", "), ")...")
    obj <- SCTransform(
      obj,
      vst.flavor          = "v2",
      method              = "glmGamPoi",
      vars.to.regress     = cfg$regress_vars,
      variable.features.n = cfg$n_hvg,
      verbose             = TRUE
    )
    message("  DefaultAssay: ", DefaultAssay(obj))

    all_hvg     <- VariableFeatures(obj)
    cleaned_hvg <- filter_hvg(all_hvg, cfg$exclude_patterns)
    VariableFeatures(obj) <- cleaned_hvg
    message("  HVGs: ", length(all_hvg), " → ", length(cleaned_hvg),
            " after removing MT/ribo/TCR/Ig")
    write.csv(data.frame(gene = cleaned_hvg),
              file.path(data_dir, "sct_hvg_list.csv"), row.names = FALSE)

    if ("ADT" %in% names(obj@assays)) {
      message("  Re-normalizing ADT (CLR) on T cell subset...")
      obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR",
                           margin = 2, verbose = FALSE)
    }

    saveRDS(obj, cfg$ckpt_sct)
    message("  Checkpoint saved: ", basename(cfg$ckpt_sct))
  }

  # ── Stage: pca ────────────────────────────────────────────────────────────
  if (run_stage("pca")) {
    cleaned_hvg <- VariableFeatures(obj)
    message("\n[pca] RunPCA (n=", cfg$n_pcs, " dims)...")
    obj <- RunPCA(obj, features = cleaned_hvg, npcs = cfg$n_pcs, verbose = FALSE)

    p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs) +
      ggtitle("PCA Elbow Plot — T Cell Subset (SCTransform v2)") +
      theme_cowplot(14)
    save_figure(p_elbow, "01_elbow_plot_SCT", plot_dir, 8, 5)

    saveRDS(obj, cfg$ckpt_pca)
    message("  Checkpoint saved: ", basename(cfg$ckpt_pca))
  }

  # ── Stage: scan ───────────────────────────────────────────────────────────
  if (run_stage("scan")) {
    message("\n[scan] DimRezScanner (dims: ",
            paste(cfg$scan_dims, collapse = ", "),
            "  |  res: ", paste(cfg$scan_res, collapse = ", "), ")...")
    key_features <- c("CD3D","CD4","CD8A","CCR7","SELL","IL7R",
                      "GZMB","GZMK","FOXP3","MKI67","LAG3","PDCD1",
                      "ISG15","CX3CR1","NKG7")
    run_dimrez_scanner(
      obj         = obj,
      scan_dir    = scan_dir,
      dims_vec    = cfg$scan_dims,
      resolutions = cfg$scan_res,
      features    = key_features,
      adt_markers = TCELL_ADT_MARKERS
    )
  }

  # ── Stage: cluster ────────────────────────────────────────────────────────
  if (run_stage("cluster")) {
    message("\n[cluster] RunUMAP + FindNeighbors + FindClusters + clustree",
            " (dims=1–", max(cfg$pca_dims), ")...")
    obj <- RunUMAP(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)

    for (res in cfg$resolutions) {
      obj <- FindClusters(obj, resolution = res, verbose = FALSE)
      message("  res=", res, " → ",
              length(unique(obj@meta.data[[paste0("SCT_snn_res.", res)]])),
              " clusters")
    }

    p_clustree <- clustree::clustree(obj@meta.data, prefix = "SCT_snn_res.") +
      ggtitle("Clustree — SCT v2 Resolution Sweep") +
      theme(legend.position = "right")
    save_figure(p_clustree, "02_clustree_SCT", plot_dir, 12, 14)

    Idents(obj) <- paste0("SCT_snn_res.", cfg$final_res)
    obj$seurat_clusters <- Idents(obj)
    message("  Active resolution: ", cfg$final_res, " → ",
            length(unique(Idents(obj))), " clusters")

    saveRDS(obj, cfg$ckpt_cluster)
    message("  Checkpoint saved: ", basename(cfg$ckpt_cluster))
  }

  # ── Stage: umap ───────────────────────────────────────────────────────────
  if (run_stage("umap")) {
    message("\n[umap] UMAP visualization plots...")
    key_features <- c("CD3D","CD4","CD8A","CCR7","SELL","IL7R",
                      "GZMB","GZMK","FOXP3","MKI67","LAG3","PDCD1",
                      "ISG15","CX3CR1","NKG7")

    p_clust <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                       label = TRUE, label.size = 5, repel = TRUE) +
      ggtitle(paste0("T Cell Clusters — SCT v2 (res=", cfg$final_res, ")")) +
      theme_cowplot(14) + NoLegend()
    save_figure(p_clust, "03_umap_clusters_SCT", plot_dir,
                cfg$umap_width, cfg$umap_height)

    parent_cols <- c("CD8+ T" = "#2196A3", "Memory CD4+ T" = "#3F51B5",
                     "Naive CD4+ T" = "#7B2D8B")
    if (cfg$celltype_col %in% colnames(obj@meta.data)) {
      p_parent <- DimPlot(obj, reduction = "umap", group.by = cfg$celltype_col,
                          cols = parent_cols, pt.size = 0.5) +
        ggtitle("Parent Cell Type") + theme_cowplot(14)
      save_figure(p_parent, "04_umap_parent_celltype", plot_dir,
                  cfg$umap_width + 2, cfg$umap_height)
    }

    status_cols <- c("Hi_Cog" = "#D95F5F", "Lo_Cog" = "#6BAED6")
    p_status <- DimPlot(obj, reduction = "umap", group.by = cfg$group_col,
                        cols = status_cols, pt.size = 0.5) +
      ggtitle("Cognitive Status (Hi_Cog vs Lo_Cog)") + theme_cowplot(14)
    save_figure(p_status, "05_umap_status", plot_dir,
                cfg$umap_width + 2, cfg$umap_height)

    donor_pal <- setNames(
      colorRampPalette(brewer.pal(8, "Set1"))(
        length(unique(obj@meta.data[[cfg$donor_col]]))),
      sort(unique(obj@meta.data[[cfg$donor_col]]))
    )
    p_donor <- DimPlot(obj, reduction = "umap", group.by = cfg$donor_col,
                       cols = donor_pal, pt.size = 0.4) +
      ggtitle("Donor") + theme_cowplot(14)
    save_figure(p_donor, "06_umap_donor", plot_dir,
                cfg$umap_width + 2, cfg$umap_height)

    p_phase <- DimPlot(obj, reduction = "umap", group.by = "Phase",
                       cols = c(G1 = "#A6A6A6", S = "#F4A460", G2M = "#DC143C"),
                       pt.size = 0.4) +
      ggtitle("Cell Cycle Phase (regressed in SCT)") + theme_cowplot(14)
    save_figure(p_phase, "07_umap_cell_cycle", plot_dir,
                cfg$umap_width, cfg$umap_height)

    DefaultAssay(obj) <- "SCT"
    feat_present <- intersect(key_features, rownames(obj))
    if (length(feat_present) > 0) {
      p_feat <- FeaturePlot(obj, features = feat_present, reduction = "umap",
                            ncol = 4, order = TRUE, pt.size = 0.3) &
        theme_cowplot(10)
      save_figure(p_feat, "08_feature_key_markers_SCT", plot_dir, 18, 10)
    }

    plot_sct_qc_umaps(obj, plot_dir)
  }

  # ── Stage: markers ────────────────────────────────────────────────────────
  if (run_stage("markers")) {
    message("\n[markers] PrepSCTFindMarkers + FindAllMarkers + dotplot...")
    obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)

    plot_genes <- intersect(unique(unlist(TCELL_MARKERS)), rownames(obj))
    p_dot <- DotPlot(obj, features = plot_genes, group.by = "seurat_clusters",
                     assay = "SCT", dot.scale = 6, col.min = -2, col.max = 2) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0) +
      coord_flip() + theme_cowplot(11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)) +
      labs(title = "Canonical T Cell Markers by Cluster (SCT v2)",
           x = NULL, y = "Cluster")
    save_figure(p_dot, "09_dotplot_canonical_markers_SCT", plot_dir,
                cfg$dot_width, cfg$dot_height + 4)

    message("  FindAllMarkers (SCT assay, Wilcoxon)...")
    markers     <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE,
                                   min.pct = 0.1, logfc.threshold = 0.25,
                                   test.use = "wilcox", verbose = FALSE)
    markers_sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, ]
    write.csv(markers_sig,
              file.path(data_dir, "cluster_markers_significant_SCT.csv"),
              row.names = FALSE)
    write.csv(markers,
              file.path(data_dir, "cluster_markers_all_SCT.csv"),
              row.names = FALSE)
    message("  Significant markers: ", nrow(markers_sig))

    top5 <- markers_sig %>% group_by(cluster) %>%
      slice_max(avg_log2FC, n = 5) %>% pull(gene) %>% unique()
    if (length(top5) > 0 && length(top5) <= 200) {
      p_heat <- DoHeatmap(obj, features = top5, group.by = "seurat_clusters",
                          assay = "SCT", size = 3, angle = 45) +
        scale_fill_gradientn(colors = c("#2166AC","white","#B2182B")) +
        theme(axis.text.y = element_text(size = 6)) +
        ggtitle("Top 5 Markers per Cluster — SCT v2 Pearson Residuals")
      save_figure(p_heat, "10_heatmap_top5_markers_SCT", plot_dir, 18, 12)
    }
  }

  # ── Stage: de ─────────────────────────────────────────────────────────────
  if (run_stage("de")) {
    message("\n[de] Proportions (Propeller + MASC) + pseudobulk DESeq2...")

    # Proportion bar by donor
    .donor   <- as.vector(as.character(obj@meta.data[[cfg$donor_col]]))
    .group   <- as.vector(as.character(obj@meta.data[[cfg$group_col]]))
    .cluster <- as.vector(as.character(obj$seurat_clusters))
    prop_df  <- as.data.frame(
      table(setNames(list(.donor, .group, .cluster),
                     c(cfg$donor_col, cfg$group_col, "seurat_clusters"))),
      stringsAsFactors = FALSE
    )
    prop_df  <- prop_df[prop_df$Freq > 0, ]
    names(prop_df)[names(prop_df) == "Freq"] <- "n"
    donor_totals <- tapply(prop_df$n, prop_df[[cfg$donor_col]], sum)
    prop_df$prop <- prop_df$n / donor_totals[prop_df[[cfg$donor_col]]]
    prop_df$seurat_clusters <- factor(prop_df$seurat_clusters)

    p_prop <- ggplot(prop_df,
                     aes(x = .data[[cfg$donor_col]], y = prop,
                         fill = seurat_clusters)) +
      geom_bar(stat = "identity", position = "stack", width = 0.8) +
      facet_grid(~ .data[[cfg$group_col]], scales = "free_x", space = "free_x") +
      labs(title = "T Cell Subcluster Proportions by Donor (SCT v2)",
           x = "Donor", y = "Proportion", fill = "Cluster") +
      theme_cowplot(13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_figure(p_prop, "11_proportions_by_donor_SCT", plot_dir, 14, 7)

    prop_res <- tryCatch(
      run_tcell_propeller(obj, cluster_col = "seurat_clusters",
        donor_col = cfg$donor_col, group_col = cfg$group_col,
        hi_group  = cfg$hi_group, lo_group   = cfg$lo_group),
      error = function(e) { message("  Propeller failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(prop_res)) {
      write.csv(prop_res,
                file.path(data_dir, "propeller_HiCog_vs_LoCog_SCT.csv"))
      message("  Propeller complete:"); print(head(prop_res))
    }

    masc_res <- tryCatch(
      run_masc(obj, cluster_col = "seurat_clusters",
        donor_col = cfg$donor_col, group_col = cfg$group_col,
        hi_group  = cfg$hi_group, lo_group   = cfg$lo_group),
      error = function(e) { message("  MASC failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(masc_res)) {
      write.csv(masc_res,
                file.path(data_dir, "MASC_HiCog_vs_LoCog_SCT.csv"),
                row.names = FALSE)
      p_masc <- ggplot(masc_res[!is.na(masc_res$p_value), ],
                       aes(x = reorder(cluster, estimate), y = estimate,
                           fill = direction, alpha = p_adj < 0.05)) +
        geom_col() +
        geom_errorbar(aes(ymin = estimate - 1.96 * se,
                          ymax = estimate + 1.96 * se), width = 0.3) +
        scale_fill_manual(values = c("Hi_Cog"="#D95F5F","Lo_Cog"="#6BAED6")) +
        scale_alpha_manual(values = c("TRUE"=1,"FALSE"=0.4), name="padj<0.05") +
        coord_flip() +
        labs(title = "MASC: T Cell Subcluster Abundance\n(Hi_Cog vs Lo_Cog — SCT v2)",
             x = "Cluster", y = "Log-odds (Hi_Cog vs Lo_Cog)",
             fill = "Enriched in") +
        theme_cowplot(13)
      save_figure(p_masc, "12_MASC_cluster_abundance_SCT", plot_dir, 10, 7)
    }

    message("  Pseudobulk DESeq2 (raw RNA counts)...")
    de_results <- pseudobulk_tcell_deseq2(obj,
      cluster_col = "seurat_clusters", donor_col = cfg$donor_col,
      group_col   = cfg$group_col,     hi_group  = cfg$hi_group,
      lo_group    = cfg$lo_group)
    if (length(de_results) > 0) {
      de_all <- do.call(rbind, de_results)
      write.csv(de_all,
                file.path(data_dir, "pseudobulk_DESeq2_HiCog_vs_LoCog_SCT.csv"),
                row.names = FALSE)
      message("  DESeq2 results saved (", length(de_results), " clusters)")
    }
  }

  # ── Done ──────────────────────────────────────────────────────────────────
  message("\n== Stages complete: ", paste(active, collapse = ", "), " ==")
  message("Outputs in: ", cfg$out_dir)

  invisible(list(obj        = obj,
                 markers    = markers_sig,
                 de_results = de_results,
                 masc       = masc_res,
                 propeller  = prop_res))
}

# =============================================================================
##  RUN
# =============================================================================
res <- run_sct_pipeline(cfg)

# =============================================================================
##  HOW SCTRANSFORM V2 WORKS — IN DEPTH
# =============================================================================
#
#  PROBLEM SCTRANSFORM SOLVES
#  ─────────────────────────
#  Standard log-normalization (log(CPM+1)) assumes every cell had the same
#  amount of mRNA captured. In practice, total UMI count varies 10–100-fold
#  across cells (due to cell size, RNA content, capture efficiency). When you
#  log-normalize, highly sequenced cells still appear artificially high for
#  thousands of genes — this technical variation drives the first PCs and
#  dominates clustering. You can't simply "regress out nCount" afterwards
#  because the relationship is non-linear and gene-specific.
#
#  THE GENERALIZED LINEAR MODEL (GLM) FRAMEWORK
#  ─────────────────────────────────────────────
#  For each gene g and cell i, SCTransform fits:
#
#    E[Y_gi] = mu_gi
#    log(mu_gi) = beta_0g + beta_1g * log10(N_i) + gamma_g * X_i
#
#  where:
#    Y_gi   = raw UMI count for gene g in cell i
#    N_i    = total UMI count in cell i (sequencing depth)
#    X_i    = additional covariates to regress (e.g. percent.mt, cc scores)
#    beta_0g = gene-specific intercept (average expression)
#    beta_1g = gene-specific slope for sequencing depth
#    gamma_g = coefficients for covariates in vars.to.regress
#
#  The count distribution follows a NEGATIVE BINOMIAL (overdispersed Poisson),
#  which is appropriate for scRNA-seq because counts have extra variability
#  beyond what a Poisson would predict (biological + technical noise).
#
#  V2 IMPROVEMENT — FIXED SLOPE
#  ─────────────────────────────
#  In v1, beta_1g (the depth slope) was estimated separately per gene.
#  This was unstable for lowly expressed genes (few observations, noisy fits).
#  Lause et al. 2021 showed theoretically that the "natural" slope for a
#  sequencing depth predictor on a log10 scale is ln(10) ≈ 2.303.
#
#  In v2: beta_1g is FIXED to ln(10) for all genes.
#  Only beta_0g and gamma_g are estimated. This makes fitting much more
#  stable, especially for sparse genes, and removes a major source of
#  estimation bias that caused v1 to artificially inflate residuals for
#  lowly expressed genes.
#
#  PEARSON RESIDUALS AS NORMALIZED EXPRESSION
#  ───────────────────────────────────────────
#  After fitting the GLM, for each gene/cell the Pearson residual is:
#
#    r_gi = (Y_gi - mu_gi) / sqrt(mu_gi + mu_gi^2 / theta_g)
#
#  where theta_g is the per-gene overdispersion parameter from the NB model.
#
#  Interpretation: r_gi is roughly "how many standard deviations above/below
#  the expected count is this observation, given the cell's depth and the
#  gene's mean expression?" A residual of 0 means the gene is expressed
#  exactly as much as the model predicts given depth. A residual of +3 means
#  the gene is dramatically overexpressed relative to expectations.
#
#  Because depth is already modeled (and its effect subtracted), these
#  residuals are NOT driven by total UMI — which is the core goal.
#
#  V2 IMPROVEMENT — LOWER BOUND ON SD
#  ─────────────────────────────────────
#  In v1, genes with near-zero expression could have extremely small predicted
#  mu_gi and very small denominators in the Pearson residual, producing
#  enormous inflated residuals from even a single count. v2 places a lower
#  bound on the gene-level standard deviation used to scale residuals,
#  preventing these pathological cases from distorting HVG selection and PCA.
#
#  HVG SELECTION FROM RESIDUAL VARIANCE
#  ──────────────────────────────────────
#  After computing Pearson residuals for all cells, SCTransform ranks genes
#  by the variance of their residuals across all cells. Genes with high
#  residual variance are genuinely variable in expression AFTER accounting
#  for depth — these are the biologically meaningful HVGs. The top n_hvg
#  (here 3000) are selected as variable features for PCA.
#
#  This is fundamentally different from the "mean-variance trend" approach
#  (FindVariableFeatures with method="vst"), which corrects only empirically
#  for the mean-variance relationship in log-normalized data and is less
#  theoretically grounded.
#
#  glmGamPoi BACKEND
#  ─────────────────
#  Fitting a GLM per gene for 10,000–50,000 cells × 20,000+ genes would be
#  prohibitively slow in standard R. glmGamPoi uses a specialized algorithm
#  for gamma-Poisson (= negative binomial) GLMs that leverages sparse matrix
#  structure and vectorized C++ code, achieving 5–20× speedup over the
#  original IRLS fitting in v1.
#
#  WHAT IS STORED IN THE SCT ASSAY
#  ──────────────────────────────────
#  After SCTransform, the Seurat object gains an "SCT" assay with three layers:
#    counts  : "corrected counts" — reverse-transformed from the fitted model
#              at a fixed sequencing depth (geometric mean across cells). These
#              are useful for visualization (FeaturePlot, DotPlot, VlnPlot)
#              and for PrepSCTFindMarkers / FindMarkers.
#    data    : log1p(corrected counts) — log-transformed version of above.
#    scale.data : Pearson residuals — used for PCA, UMAP, and clustering.
#
#  WHY NOT REGRESS nCount_RNA IN vars.to.regress
#  ───────────────────────────────────────────────
#  This is a common mistake. SCTransform already handles sequencing depth via
#  the GLM offset (the log10(N_i) term). Adding nCount_RNA to vars.to.regress
#  would double-subtract it — removing true biological signal correlated with
#  cell size. Only regress NUISANCE variables that are NOT modeled by the GLM:
#  percent.mt (apoptosis / low quality), S.Score, G2M.Score (cell cycle).
#
#  PrepSCTFindMarkers BEFORE FindAllMarkers
#  ─────────────────────────────────────────
#  SCTransform fits the NB model per gene using all cells in the object at
#  the time of fitting. After subsetting (e.g., to T cells only), those
#  corrected counts are no longer consistent — the model was fit on a larger
#  gene × cell space. PrepSCTFindMarkers re-runs the reverse transformation
#  ("recorrection") using raw RNA counts within the current object, producing
#  a corrected counts matrix that is valid for the current cell set and
#  appropriate for Wilcoxon / LR tests in FindAllMarkers.
#
#  PSEUDOBULK DESeq2 USES RAW RNA COUNTS, NOT SCT
#  ─────────────────────────────────────────────────
#  DESeq2's negative binomial dispersion estimation requires raw integer
#  counts. SCT corrected counts are on a continuous, floating-point scale
#  anchored to a reference depth, not to the actual per-donor library size.
#  Aggregating SCT counts into pseudobulks would produce invalid size factors
#  and miscalibrated dispersions. Always aggregate raw RNA counts for DESeq2.
