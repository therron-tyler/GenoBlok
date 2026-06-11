#!/usr/bin/env Rscript

## Optional: set custom lib path first, then append existing ones
.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(cowplot)
})

# If you want explicit HTO order per object, pass a named list:
# hto_names_map <- list(
#   Obj2 = c("HTO1","HTO2","HTO3"),
#   Obj3 = c("HTO1","HTO2","HTO3"),
#   Obj4 = c("HTO7","HTO8","HTO9"),
#   Obj5 = c("HTO7","HTO8","HTO9")
# )

add_ribo_metrics <- function(obj, ribo_features = NULL, assay = "RNA") {
  DefaultAssay(obj) <- assay

  if (is.null(ribo_features)) {
    ribo_features <- grep("^RP[SL][[:digit:]]|^RPSA", rownames(obj), value = TRUE)
  }
  ribo_features <- intersect(ribo_features, rownames(obj))
  if (length(ribo_features) == 0) {
    warning("No ribosomal genes matched; setting percent.ribo and nCount_Ribo to NA.")
    obj[["percent.ribo"]] <- NA_real_
    obj[["nCount_Ribo"]]  <- NA_real_
    return(obj)
  }

  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_features)

  ribo_counts <- GetAssayData(obj, assay = assay, layer = "counts")[ribo_features, , drop = FALSE]
  obj[["nCount_Ribo"]] <- Matrix::colSums(ribo_counts)

  obj
}

get_hto_order_for_object <- function(obj, nm, hto_names = NULL, hto_names_map = NULL) {
  if (!is.null(hto_names_map) && nm %in% names(hto_names_map)) {
    return(hto_names_map[[nm]])
  }
  if (!is.null(hto_names)) {
    return(hto_names)
  }

  lev <- levels(Idents(obj))
  lev <- setdiff(lev, c("Negative", "Doublet"))
  return(lev)
}

infer_techrep_tag <- function(nm, techrep_map = NULL) {
  if (!is.null(techrep_map) && nm %in% names(techrep_map)) {
    return(techrep_map[[nm]])
  }

  m1 <- regexpr("([A-D])$", nm, perl = TRUE)
  if (m1[1] != -1) return(substr(nm, m1[1], m1[1]))

  m2 <- regexpr("Hash[0-9]+([A-D])$", nm, perl = TRUE)
  if (m2[1] != -1) return(sub("^.*Hash[0-9]+", "", nm))

  m3 <- regexpr("(run[0-9]+)$", tolower(nm), perl = TRUE)
  if (m3[1] != -1) return(substr(tolower(nm), m3[1], m3[1] + attr(m3, "match.length") - 1))

  return(nm)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  get_flag <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx == length(args)) stop("Missing value after ", flag)
    args[[idx + 1]]
  }

  get_multi_flag <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(character(0))
    vals <- vapply(idx, function(i) {
      if (i == length(args)) stop("Missing value after ", flag)
      args[[i + 1]]
    }, character(1))
    vals
  }

  a <- list(
    rds_paths   = get_multi_flag("--rds"),
    names       = get_flag("--names", ""),
    outdir      = get_flag("--outdir"),

    nf_low      = as.numeric(get_flag("--nf_low",      "0")),
    nf_high     = as.numeric(get_flag("--nf_high",     "6000")),
    nc_low      = as.numeric(get_flag("--nc_low",      "2000")),
    nc_high     = as.numeric(get_flag("--nc_high",     "30000")),
    nc_ribo_min = as.numeric(get_flag("--nc_ribo_min", "150")),
    mt_max      = as.numeric(get_flag("--mt_max",      "18")),
    ribo_max    = as.numeric(get_flag("--ribo_max",    "50")),
    nfeatures   = as.integer(get_flag("--nfeatures",   "2000")),
    dpi         = as.integer(get_flag("--dpi",         "300")),
    vln_width   = as.numeric(get_flag("--vln_width",   "14")),
    vln_height  = as.numeric(get_flag("--vln_height",  "8")),
    mt_pattern  = get_flag("--mt_pattern", "^MT-"),
    hto_names   = get_flag("--hto_names",  "")
  )

  if (is.null(a$outdir) || length(a$rds_paths) < 2) {
    stop(
      "Usage:\n",
      "  Rscript merge_hashed_replicates.R \\\n",
      "    --rds <obj1.rds> --rds <obj2.rds> [--rds <obj3.rds> ...] \\\n",
      "    --outdir <dir> [--names name1,name2,...] [--hto_names HTO1,HTO2,...]\n\n",
      "Notes:\n",
      "  - Provide at least 2 --rds inputs.\n",
      "  - If --names is omitted, names are auto-generated from filenames.\n"
    )
  }

  if (!is.null(a$names) && nchar(a$names) > 0) {
    a$names <- trimws(strsplit(a$names, ",")[[1]])
    if (length(a$names) != length(a$rds_paths)) {
      stop("--names must have the same number of entries as --rds flags.")
    }
  } else {
    a$names <- sub("\\.rds$", "", basename(a$rds_paths), ignore.case = TRUE)
  }

  if (!is.null(a$hto_names) && nchar(a$hto_names) > 0) {
    a$hto_names <- trimws(strsplit(a$hto_names, ",")[[1]])
  } else {
    a$hto_names <- NULL
  }

  return(a)
}

order_idents_by_hto <- function(obj, hto_names) {
  if (is.null(hto_names)) return(obj)

  cur <- as.character(Idents(obj))
  ordered_levels <- c(intersect(hto_names, unique(cur)),
                      setdiff(unique(cur), hto_names))
  Idents(obj) <- factor(cur, levels = ordered_levels)
  obj
}

# ---- QC plotting helpers (style matches singlecell_figure_engine) ----

plot_qc_by_group <- function(obj, group_col, title = "QC") {
  qc_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_features) == 0) { warning("No QC fields found."); return(NULL) }

  plist <- VlnPlot(obj, features = qc_features, group.by = group_col,
                   pt.size = 0, combine = FALSE)
  n <- length(plist)
  for (i in seq_along(plist)) {
    plist[[i]] <- plist[[i]] +
      ggtitle(qc_features[i]) +
      theme(
        legend.position = "none",
        plot.title      = element_text(face = "bold", size = 22),
        axis.title.x    = element_blank(),
        axis.text.x     = if (i < n) element_blank() else element_text(angle = 45, hjust = 1, size = 19),
        axis.ticks.x    = if (i < n) element_blank() else element_line()
      )
  }
  wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_qc_by_sample <- function(obj, sample_col, pool_col, pool_cols,
                               title = "QC by Sample") {
  qc_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_features) == 0) return(NULL)

  md <- as.data.frame(obj@meta.data[, c(sample_col, pool_col), drop = FALSE])
  md[] <- lapply(md, as.character)
  colnames(md) <- c("Sample", "Pool")

  pool_df <- md %>%
    dplyr::count(Sample, Pool) %>%
    dplyr::group_by(Sample) %>%
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Sample, Pool)

  pool_df <- pool_df[order(pool_df$Pool, pool_df$Sample), ]
  sample_order <- pool_df$Sample

  sample_col_vals <- setNames(
    ifelse(pool_df$Pool %in% names(pool_cols), pool_cols[pool_df$Pool], "grey60"),
    pool_df$Sample
  )

  obj@meta.data[[sample_col]] <- factor(
    as.character(obj@meta.data[[sample_col]]),
    levels = sample_order
  )

  plist <- VlnPlot(obj, features = qc_features, group.by = sample_col,
                   pt.size = 0, combine = FALSE, cols = sample_col_vals)
  n <- length(plist)
  for (i in seq_along(plist)) {
    plist[[i]] <- plist[[i]] +
      ggtitle(qc_features[i]) +
      theme(
        legend.position = "none",
        plot.title      = element_text(face = "bold", size = 22),
        axis.title.x    = element_blank(),
        axis.text.x     = if (i < n) element_blank() else element_text(angle = 45, hjust = 1, size = 19),
        axis.ticks.x    = if (i < n) element_blank() else element_line()
      )
  }

  wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_qc_scatter_cutoff <- function(obj, x_feat = "nCount_RNA", y_feat = "nCount_Ribo",
                                    color_by = "hash_run",
                                    x_cutoff = NULL, y_cutoff = NULL,
                                    x_max = NULL, y_max = NULL,
                                    title = "") {
  feats <- intersect(c(x_feat, y_feat, color_by), colnames(obj@meta.data))
  md    <- as.data.frame(obj@meta.data[, feats, drop = FALSE])

  p <- ggplot(md, aes(x = .data[[x_feat]], y = .data[[y_feat]],
                       color = .data[[color_by]])) +
    geom_point(size = 0.2, alpha = 0.3) +
    labs(x = x_feat, y = y_feat, title = title, color = color_by) +
    theme_classic(base_size = 14) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  if (!is.null(x_cutoff))
    p <- p + geom_vline(xintercept = x_cutoff, linetype = "dashed",
                         color = "firebrick", linewidth = 0.9)
  if (!is.null(y_cutoff))
    p <- p + geom_hline(yintercept = y_cutoff, linetype = "dashed",
                         color = "firebrick", linewidth = 0.9)
  if (!is.null(x_max) || !is.null(y_max))
    p <- p + coord_cartesian(
      xlim = if (!is.null(x_max)) c(0, x_max) else NULL,
      ylim = if (!is.null(y_max)) c(0, y_max) else NULL
    )
  p
}

# ---- Main pipeline function ----

process_seurat_list_compare_merge <- function(
    rds_list,
    outdir,
    hto_names    = NULL,
    hto_names_map = NULL,
    techrep_map  = NULL,
    mt_pattern   = "^MT-",
    nf_low       = 0,
    nf_high      = 12000,
    mt_max       = 50,
    nfeatures    = 2000,
    dpi          = 300,
    vln_width    = 16,
    vln_height   = 9,
    nc_low       = 2000,
    nc_high      = 30000,
    ribo_max     = 50,
    nc_ribo_min  = 150,
    pool_cols    = c("WHL1" = "#E41A1C", "WHL2" = "#377EB8", "WHL3" = "#4DAF4A",
                     "WHL4" = "#FF7F00", "WHL5" = "#984EA3")
) {
  if (is.null(names(rds_list)) || any(names(rds_list) == "")) {
    stop("rds_list must be a *named* list (e.g., list(S1=obj1, S2=obj2)).")
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- Per-object prep: remove doublets/negatives, add metrics, tag metadata ----
  pre_objects <- lapply(names(rds_list), function(nm) {
    obj <- rds_list[[nm]]
    DefaultAssay(obj) <- "RNA"

    obj$hash_run <- nm

    if ("Negative" %in% levels(Idents(obj))) obj <- subset(obj, idents = "Negative", invert = TRUE)
    if ("Doublet"  %in% levels(Idents(obj))) obj <- subset(obj, idents = "Doublet",  invert = TRUE)

    obj <- add_ribo_metrics(obj, assay = "RNA")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

    hto_order <- get_hto_order_for_object(obj, nm, hto_names = hto_names,
                                           hto_names_map = hto_names_map)
    obj <- order_idents_by_hto(obj, hto_order)

    # Derive clean sample IDs from HTO metadata if present (WHL2-5 may have
    # cluster-number Idents after prior processing; HTO_maxID/hash.ID preserve
    # the original demux identity). Rename IMPACT* -> IMP* for WHL1.
    if ("HTO_maxID" %in% colnames(obj@meta.data)) {
      raw_ids <- as.character(obj@meta.data$HTO_maxID)
    } else if ("hash.ID" %in% colnames(obj@meta.data)) {
      raw_ids <- as.character(obj@meta.data$hash.ID)
    } else {
      raw_ids <- as.character(Idents(obj))
    }
    obj$HTO_maxID <- gsub("^IMPACT", "IMP", raw_ids)

    # Append techrep suffix to Idents for merge disambiguation only
    tag    <- infer_techrep_tag(nm, techrep_map = techrep_map)
    suffix <- paste0("_", tag)
    obj$techrep  <- tag

    Idents(obj) <- factor(paste0(as.character(Idents(obj)), suffix))
    desired_levels <- paste0(hto_order, suffix)
    present_levels <- intersect(desired_levels, levels(Idents(obj)))
    Idents(obj) <- factor(
      as.character(Idents(obj)),
      levels = c(present_levels, setdiff(levels(Idents(obj)), present_levels))
    )

    obj
  })
  names(pre_objects) <- names(rds_list)

  # ---- Merge ----
  merged <- merge(x = pre_objects[[1]], y = pre_objects[-1],
                  add.cell.ids = names(pre_objects))
  DefaultAssay(merged) <- "RNA"

  # Seurat v5: join layered assays before filtering / HVG
  if (exists("JoinLayers")) merged <- JoinLayers(merged, assay = "RNA")

  # Combine HTO identity + pool to distinguish tech replicates in per-sample plots
  merged$sample_id <- paste0(merged$HTO_maxID, "_", merged$hash_run)

  merged <- add_ribo_metrics(merged, assay = "RNA")
  merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = mt_pattern)

  # ---- Pre-filter QC plots ----
  p_pre_run <- plot_qc_by_group(
    merged, group_col = "hash_run",
    title = "Pre-filter QC by WHL run"
  )
  ggsave(file.path(outdir, "QC_pre_filter_violin_by_run.jpeg"),
         p_pre_run, width = vln_width, height = vln_height + 4, dpi = dpi)

  p_pre_sample <- plot_qc_by_sample(
    merged, sample_col = "sample_id", pool_col = "hash_run",
    pool_cols = pool_cols, title = "Pre-filter QC by Sample"
  )
  ggsave(file.path(outdir, "QC_pre_filter_violin_by_sample.jpeg"),
         p_pre_sample, width = vln_width + 4, height = vln_height + 4, dpi = dpi)

  p_pre_scatter <- plot_qc_scatter_cutoff(
    merged, x_feat = "nCount_RNA", y_feat = "nCount_Ribo",
    color_by = "hash_run",
    x_cutoff = nc_low, y_cutoff = nc_ribo_min,
    title = paste0("Pre-filter: nCount_RNA (cutoff=", nc_low,
                   ") vs nCount_Ribo (cutoff=", nc_ribo_min, ")")
  )
  ggsave(file.path(outdir, "QC_pre_filter_scatter_nCountRNA_vs_nCountRibo.png"),
         p_pre_scatter, width = 8, height = 6, dpi = dpi)

  p_pre_scatter_zoom <- plot_qc_scatter_cutoff(
    merged, x_feat = "nCount_RNA", y_feat = "nCount_Ribo",
    color_by = "hash_run",
    x_cutoff = nc_low, y_cutoff = nc_ribo_min,
    x_max = 20000, y_max = 5000,
    title = paste0("Pre-filter: nCount_RNA vs nCount_Ribo (zoomed)")
  )
  ggsave(file.path(outdir, "QC_pre_filter_scatter_nCountRNA_vs_nCountRibo_zoom.png"),
         p_pre_scatter_zoom, width = 8, height = 6, dpi = dpi)

  # ---- Single unified filter on merged object ----
  merged <- subset(
    merged,
    subset = nFeature_RNA > nf_low  & nFeature_RNA < nf_high  &
             percent.mt  < mt_max   &
             nCount_RNA  >= nc_low  & nCount_RNA   < nc_high  &
             percent.ribo < ribo_max &
             nCount_Ribo >= nc_ribo_min
  )

  merged <- add_ribo_metrics(merged, assay = "RNA")
  merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = mt_pattern)

  # ---- Post-filter QC plots ----
  p_post_run <- plot_qc_by_group(
    merged, group_col = "hash_run",
    title = "Post-filter QC by WHL run"
  )
  ggsave(file.path(outdir, "QC_post_filter_violin_by_run.jpeg"),
         p_post_run, width = vln_width, height = vln_height + 4, dpi = dpi)

  p_post_sample <- plot_qc_by_sample(
    merged, sample_col = "sample_id", pool_col = "hash_run",
    pool_cols = pool_cols, title = "Post-filter QC by Sample"
  )
  ggsave(file.path(outdir, "QC_post_filter_violin_by_sample.jpeg"),
         p_post_sample, width = vln_width + 4, height = vln_height + 4, dpi = dpi)

  p_post_scatter <- plot_qc_scatter_cutoff(
    merged, x_feat = "nCount_RNA", y_feat = "nCount_Ribo",
    color_by = "hash_run",
    title = "Post-filter: nCount_RNA vs nCount_Ribo"
  )
  ggsave(file.path(outdir, "QC_post_filter_scatter_nCountRNA_vs_nCountRibo.png"),
         p_post_scatter, width = 8, height = 6, dpi = dpi)

  p_post_scatter_zoom <- plot_qc_scatter_cutoff(
    merged, x_feat = "nCount_RNA", y_feat = "nCount_Ribo",
    color_by = "hash_run",
    x_max = 20000, y_max = 5000,
    title = "Post-filter: nCount_RNA vs nCount_Ribo (zoomed)"
  )
  ggsave(file.path(outdir, "QC_post_filter_scatter_nCountRNA_vs_nCountRibo_zoom.png"),
         p_post_scatter_zoom, width = 8, height = 6, dpi = dpi)

  # ---- HVG + Scaling + PCA ----
  # Rank all genes by VST variability, then exclude mito/ribo/chrY before taking top N
  all_genes  <- rownames(merged)
  mito_hvg   <- grep(mt_pattern, all_genes, value = TRUE)
  ribo_hvg   <- grep("^RP[SL][[:digit:]]|^RPSA", all_genes, value = TRUE)
  chrY_hvg   <- grep(
    "^DDX3Y|^EIF1AY|^KDM5D|^NLGN4Y|^RPS4Y|^TMSB4Y|^USP9Y|^UTY|^ZFY|^PCDH11Y|^TTTY",
    all_genes, value = TRUE
  )
  exclude_hvg <- unique(c(mito_hvg, ribo_hvg, chrY_hvg))

  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = nrow(merged))
  hvg_ranked  <- VariableFeatures(merged)
  VariableFeatures(merged) <- head(setdiff(hvg_ranked, exclude_hvg), nfeatures)

  top10  <- head(VariableFeatures(merged), 10)
  vf1    <- VariableFeaturePlot(merged)
  varplt <- LabelPoints(plot = vf1, points = top10, repel = TRUE) + ggtitle("Merged HVG")

  merged <- ScaleData(merged, features = rownames(merged), verbose = FALSE)
  merged <- RunPCA(merged, features = VariableFeatures(merged), verbose = FALSE)

  pca   <- DimPlot(merged, reduction = "pca", label = FALSE, repel = TRUE,
                   split.by = "hash_run") +
    NoLegend() + ggtitle("Merged PCA (split by WHL run)")
  elbow <- ElbowPlot(merged) + ggtitle("Merged Elbow")

  ggsave(file.path(outdir, "Merged_VariableFeatures.pdf"),    varplt, width = 7,  height = 5,  dpi = dpi)
  ggsave(file.path(outdir, "Merged_PCA_splitBy_run.jpeg"),    pca,    width = 14, height = 7,  dpi = dpi)
  ggsave(file.path(outdir, "Merged_Elbow.pdf"),               elbow,  width = 7,  height = 5,  dpi = dpi)

  # ---- UMAP (dims 1-30) ----
  merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

  # ---- Clustering (dims 1-9, 11-12 | res 0.3) ----
  merged <- FindNeighbors(merged, dims = c(1:9, 11:12), verbose = FALSE)
  merged <- FindClusters(merged,  resolution = 0.3,     verbose = FALSE)

  # ---- UMAP: clusters + WHL run ----
  umap_clust <- DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle("UMAP — initial clusters (dims 1-9,11-12, res 0.3)")
  ggsave(file.path(outdir, "UMAP_initial_clusters.jpeg"),
         umap_clust, width = 9, height = 7, dpi = dpi)

  umap_run <- DimPlot(merged, reduction = "umap", group.by = "hash_run",
                      raster = FALSE) +
    ggtitle("UMAP by WHL run") +
    theme_classic(base_size = 14)
  ggsave(file.path(outdir, "UMAP_by_hash_run.jpeg"),
         umap_run, width = 9, height = 7, dpi = dpi)

  umap_imp <- DimPlot(merged, reduction = "umap", group.by = "HTO_maxID",
                      raster = FALSE, label = FALSE) +
    ggtitle("UMAP by IMP sample") +
    theme_classic(base_size = 14)
  ggsave(file.path(outdir, "UMAP_by_IMP_sample.jpeg"),
         umap_imp, width = 10, height = 7, dpi = dpi)

  # ---- UMAP: QC metric overlays ----
  qc_feats_umap <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(merged@meta.data)
  )
  fp_list <- FeaturePlot(merged, features = qc_feats_umap, reduction = "umap",
                          raster = FALSE, combine = FALSE)
  fp_combined <- wrap_plots(fp_list, ncol = 3) +
    plot_annotation(title = "UMAP — QC metrics")
  ggsave(file.path(outdir, "UMAP_QC_metrics.jpeg"),
         fp_combined, width = 15, height = 10, dpi = dpi)

  # ---- Single output RDS ----
  saveRDS(merged, file = file.path(outdir,
    "Merged_nCountRibo150_nCountRNA2000_UMAP30_clust1-9_11-12_res03.rds"))

  invisible(merged)
}

main <- function() {
  a <- parse_args()

  objs     <- lapply(a$rds_paths, readRDS)
  rds_list <- objs
  names(rds_list) <- a$names

  process_seurat_list_compare_merge(
    rds_list    = rds_list,
    outdir      = a$outdir,
    mt_pattern  = a$mt_pattern,
    nf_low      = a$nf_low,
    nf_high     = a$nf_high,
    nc_low      = a$nc_low,
    nc_high     = a$nc_high,
    mt_max      = a$mt_max,
    ribo_max    = a$ribo_max,
    nc_ribo_min = a$nc_ribo_min,
    nfeatures   = a$nfeatures,
    dpi         = a$dpi,
    vln_width   = a$vln_width,
    vln_height  = a$vln_height,
    hto_names   = a$hto_names
  )
}

#main()


hto_map <- list(
  whl1 = c("IMPACT38","IMPACT15","IMPACT29","IMPACT25","IMPACT23","IMPACT35"),
  whl2 = c("IMP11","IMP12","IMP19","IMP21","IMP30","IMP33","IMP37","IMP39"),
  whl3 = c("IMP11","IMP12","IMP19","IMP21","IMP30","IMP33","IMP37","IMP39"),
  whl4 = c("IMP13","IMP16","IMP18","IMP24","IMP34","IMP36","IMP42","IMP46"),
  whl5 = c("IMP13","IMP16","IMP18","IMP24","IMP34","IMP36","IMP42","IMP46")
)

tech_map <- list(
  whl1 = "1",
  whl2 = "2",
  whl3 = "3",
  whl4 = "4",
  whl5 = "5"
)

WHL1 <- readRDS("/path/to/data/20260506_WHL1_SNPrecovery/WHL1_recovered_SNPconsensus.rds")
WHL2 <- readRDS("/path/to/data/20260506_WHL1_SNPrecovery/WHL2_20_singlets_RNA.rds")
WHL3 <- readRDS("/path/to/data/20260506_WHL1_SNPrecovery/WHL3_20_singlets_RNA.rds")
WHL4 <- readRDS("/path/to/data/20260506_WHL1_SNPrecovery/WHL4_20_singlets_RNA.rds")
WHL5 <- readRDS("/path/to/data/20260506_WHL1_SNPrecovery/WHL5_20_singlets_RNA.rds")

rds_list <- list(
  WHL1 = WHL1,
  WHL2 = WHL2,
  WHL3 = WHL3,
  WHL4 = WHL4,
  WHL5 = WHL5
)
outdir <- "/path/to/data/20260506_WHL1_SNPrecovery/QC_Analysis"

process_seurat_list_compare_merge(
  rds_list      = rds_list,
  outdir        = outdir,
  hto_names_map = hto_map,
  techrep_map   = tech_map,
  nc_low        = 2000,
  nc_ribo_min   = 150
)
