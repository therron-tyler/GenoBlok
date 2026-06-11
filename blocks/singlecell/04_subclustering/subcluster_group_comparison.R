#!/usr/bin/env Rscript
.libPaths(c("/path/to/your/R_library", .libPaths()))

## ============================================================================
##  NPSLE T Cell Subclustering Pipeline
##  Populations: CD8+ T, Memory CD4+ T, Naive CD4+ T
##  Group comparison: Hi_Cog vs Lo_Cog (Status column)

## ============================================================================
#BiocManager::install("speckle")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(DESeq2)
  library(speckle)       # propeller cell-type proportion testing
  library(limma)
  library(clustree)      # resolution selection via cluster tree
  library(colorspace)
  library(lme4)          # MASC mixed-effects proportion testing
  library(optparse)      # CLI argument parsing
})

# ── CLI arguments (override cfg defaults when running on HPC) ────────────────
# Usage:
#   Rscript subcluster_group_comparison.R \
#     --input_rds  /path/to/object.rds \
#     --out_dir    /path/to/output \
#     --final_res  0.4 \
#     --pca_dims   20 \
#     --n_hvg      2000
#
# All args are optional — omit to use the hardcoded cfg defaults below.

option_list <- list(
  make_option("--input_rds",  type="character", default=NULL,
              help="Path to input Seurat RDS file"),
  make_option("--out_dir",    type="character", default=NULL,
              help="Output directory"),
  make_option("--final_res",  type="double",    default=NULL,
              help="Final clustering resolution (e.g. 0.4)"),
  make_option("--pca_dims",   type="integer",   default=NULL,
              help="Number of PCA dims for UMAP/FindNeighbors (e.g. 20)"),
  make_option("--n_hvg",      type="integer",   default=NULL,
              help="Number of highly variable genes (default 2000)"),
  make_option("--regress_cc", type="logical",   default=NULL,
              help="Regress cell cycle scores TRUE/FALSE")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ── 0. Configuration ─────────────────────────────────────────────────────────
cfg <- list(

  # Paths
  input_rds = file.path(
    "/path/to/data",
    "20260213_GroupComps_NPSLEsamples",
    "NPSLE_CellTypeCog_Hash_dims8_rez3.rds"
  ),
  out_dir = file.path(
    "/path/to/data",
    "20260213_GroupComps_NPSLEsamples",
    "Tcell_Subclustering_Output"
  ),

  # Metadata columns
  celltype_col = "celltype",        # coarse cell-type labels
  donor_col    = "HTO_maxID",       # per-donor ID — pseudobulk unit
  group_col    = "Status",          # biological group: Hi_Cog / Lo_Cog
  sample_col   = "sample_origin",   # per-sample ID for QC plots (NPSLE##_A/B)
  hi_group     = "Hi_Cog",
  lo_group     = "Lo_Cog",
  assay        = "RNA",

  # Pool group colors — _A and _B pools (derived from sample_origin suffix)
  pool_cols = c("A" = "#4DAFE0", "B" = "#E0874D"),

  # T cell populations to subset
  tcell_types = c("CD8+ T", "Memory CD4+ T", "Naive CD4+ T"),

  # HVG / scaling
  n_hvg            = 2000,
  exclude_patterns = c("^MT-", "^RPS", "^RPL", "^IGHV", "^IGLV", "^IGKV",
                       "^TRAV", "^TRBV", "^TRDV", "^TRGV"),

  # Cell cycle regression
  regress_cc   = TRUE,
  regress_vars = c("nCount_RNA", "percent.mt"),

  # PCA / UMAP
  n_pcs    = 50,
  pca_dims = 1:20,   # dims for UMAP / FindNeighbors — adjust after elbow plot

  # Clustering resolutions to sweep
  resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),
  final_res   = 0.4,   # update after inspecting clustree + markers

  # Figure dimensions
  umap_width  = 10, umap_height = 8,
  dot_width   = 14, dot_height  = 9
)

# Apply any CLI overrides on top of the defaults above
if (!is.null(opt$input_rds))  cfg$input_rds  <- opt$input_rds
if (!is.null(opt$out_dir))    cfg$out_dir    <- opt$out_dir
if (!is.null(opt$final_res))  cfg$final_res  <- opt$final_res
if (!is.null(opt$pca_dims))   cfg$pca_dims   <- seq_len(opt$pca_dims)
if (!is.null(opt$n_hvg))      cfg$n_hvg      <- opt$n_hvg
if (!is.null(opt$regress_cc)) cfg$regress_cc <- opt$regress_cc

message("Pipeline config:")
message("  input_rds : ", cfg$input_rds)
message("  out_dir   : ", cfg$out_dir)
message("  final_res : ", cfg$final_res)
message("  pca_dims  : 1–", max(cfg$pca_dims))
message("  n_hvg     : ", cfg$n_hvg)

# ── Canonical T cell marker panels ───────────────────────────────────────────
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
    if (length(mt_genes) == 0) { obj$percent.mt <- 0
    } else obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
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
    warning("No QC fields found among: nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, nCount_Ribo")
    return(NULL)
  }
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

# ── QC violin plots: by sample, coloured by pool (_A vs _B) ──────────────────
# Pool group is derived from the suffix of sample_origin (NPSLE##_A → "A").
# cfg$pool_cols must be a named colour vector, e.g. c("A"="#4DAFE0","B"="#E0874D").
plot_qc_by_sample_pooled <- function(obj, cfg, title = "Single-cell QC by Sample") {
  qc_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_features) == 0) return(NULL)
  sample_col <- cfg$sample_col

  # Derive pool group from the trailing _A / _B in sample_origin
  md <- data.frame(
    Sample     = as.character(obj@meta.data[[sample_col]]),
    stringsAsFactors = FALSE
  )
  md$PoolGroup <- sub(".*_([^_]+)$", "\\1", md$Sample)   # extract suffix after last _
  md$PoolGroup[!md$PoolGroup %in% names(cfg$pool_cols)] <- "Unknown"
  md$PoolGroup <- factor(md$PoolGroup, levels = c(names(cfg$pool_cols), "Unknown"))

  pool_df <- md %>%
    distinct(Sample, PoolGroup) %>%
    arrange(PoolGroup, Sample)

  sample_order <- pool_df$Sample
  sample_cols  <- setNames(
    ifelse(as.character(pool_df$PoolGroup) %in% names(cfg$pool_cols),
           cfg$pool_cols[as.character(pool_df$PoolGroup)],
           "grey60"),
    pool_df$Sample
  )

  # Re-level so VlnPlot respects pool-sorted order
  obj@meta.data[[sample_col]] <- factor(
    as.character(obj@meta.data[[sample_col]]),
    levels = sample_order
  )

  plist <- VlnPlot(obj, features = qc_features, group.by = sample_col,
                   pt.size = 0, combine = FALSE, cols = sample_cols)
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

  # Pool group legend
  legend_df <- data.frame(
    x = 1, y = 1,
    PoolGroup = factor(names(cfg$pool_cols), levels = names(cfg$pool_cols))
  )
  p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = PoolGroup)) +
    geom_col() +
    scale_fill_manual(values = cfg$pool_cols, name = "Pool Group") +
    theme_void() +
    theme(legend.position = "right")
  legend_grob <- cowplot::get_legend(p_legend)

  combined <- wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))

  combined + inset_element(legend_grob, left = 0.85, bottom = 0.02, right = 1, top = 0.12)
}

# ── Pseudobulk DESeq2: Hi_Cog vs Lo_Cog within each T cell subcluster ────────
pseudobulk_tcell_deseq2 <- function(
  obj,
  cluster_col = "seurat_clusters",
  donor_col   = "HTO_maxID",
  group_col   = "Status",
  hi_group    = "Hi_Cog",
  lo_group    = "Lo_Cog",
  assay       = "RNA",
  min_cells   = 5,
  min_counts  = 10
) {
  stopifnot(inherits(obj, "Seurat"))
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, layer  = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
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

    # Pseudobulk: sum counts per donor
    pb_list <- lapply(unique(cl_md$donor), function(d) {
      idx <- rownames(cl_md)[cl_md$donor == d]
      if (length(idx) < min_cells) return(NULL)
      Matrix::rowSums(cl_ct[, idx, drop = FALSE])
    })
    names(pb_list) <- unique(cl_md$donor)
    pb_list <- Filter(Negate(is.null), pb_list)

    # Need donors from both groups
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
              " sig genes (padj<0.05), Hi_Cog vs Lo_Cog")
    }, error = function(e) message("  DESeq2 failed cluster ", cl, ": ", conditionMessage(e)))
  }
  results
}

# ── Propeller: Hi_Cog vs Lo_Cog subcluster proportions ───────────────────────
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

# ── MASC: Mixed-effects Association of Single Cells ──────────────────────────
# Tests whether each cluster's abundance changes with Status,
# using donor as a random effect to account for inter-individual variation.
# Reference: Fonseka et al. 2018 Journal of Experimental Medicine
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
    if (is.null(fit)) return(data.frame(cluster=cl, estimate=NA, se=NA, z=NA, p_value=NA))
    cf <- summary(fit)$coefficients
    row <- cf[grep(hi_group, rownames(cf)), , drop = FALSE]
    if (nrow(row) == 0) return(data.frame(cluster=cl, estimate=NA, se=NA, z=NA, p_value=NA))
    data.frame(
      cluster  = cl,
      estimate = row[1, "Estimate"],
      se       = row[1, "Std. Error"],
      z        = row[1, "z value"],
      p_value  = row[1, "Pr(>|z|)"]
    )
  })
  out <- do.call(rbind, results)
  out$p_adj   <- p.adjust(out$p_value, method = "BH")
  out$direction <- ifelse(is.na(out$estimate), "NA",
                          ifelse(out$estimate > 0, hi_group, lo_group))
  out[order(out$p_value), ]
}

# =============================================================================
##  MAIN PIPELINE
# =============================================================================

run_tcell_pipeline <- function(cfg) {

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(cfg$out_dir, "plots")
  data_dir <- file.path(cfg$out_dir, "data")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # ── 1. Load ──────────────────────────────────────────────────────────────
  message("\n[1/9] Loading Seurat object...")
  obj_full <- readRDS(cfg$input_rds)
  message("  Loaded: ", ncol(obj_full), " cells")
  message("  Groups: ", paste(table(obj_full@meta.data[[cfg$group_col]]), collapse = " | "))

  # ── 2. Subset T cells ────────────────────────────────────────────────────
  message("\n[2/9] Subsetting: ", paste(cfg$tcell_types, collapse = " | "))
  obj <- subset(obj_full, cells = colnames(obj_full)[
    obj_full@meta.data[[cfg$celltype_col]] %in% cfg$tcell_types
  ])
  rm(obj_full); gc()
  message("  T cells retained: ", ncol(obj))
  print(table(obj@meta.data[[cfg$celltype_col]]))
  print(table(obj@meta.data[[cfg$group_col]]))

  # Save raw subset immediately — useful checkpoint before long processing
  subset_rds <- file.path(cfg$out_dir, "subcluster_raw_subset.rds")
  saveRDS(obj, subset_rds)
  message("  Raw subset saved: ", subset_rds)

  obj <- add_mito_pct(obj, assay = cfg$assay)
  DefaultAssay(obj) <- cfg$assay

  # ── 2b. QC plots (T cell subset) ─────────────────────────────────────────
  message("\n[2b/9] QC plots (T cell subset)...")

  p_qc_group <- plot_qc_by_group(obj,
    group_col = cfg$group_col,
    title     = "T Cell QC by Cognitive Status")
  if (!is.null(p_qc_group))
    save_figure(p_qc_group, "00a_qc_by_status", plot_dir,
                width = 10, height = 4 * length(intersect(
                  c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
                  colnames(obj@meta.data))))

  p_qc_sample <- tryCatch(
    plot_qc_by_sample_pooled(obj, cfg,
      title = "T Cell QC by Sample (Pool A vs B)"),
    error = function(e) {
      message("  plot_qc_by_sample_pooled failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p_qc_sample))
    save_figure(p_qc_sample, "00b_qc_by_sample", plot_dir,
                width = 18, height = 4 * length(intersect(
                  c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
                  colnames(obj@meta.data))))

  # ── 3. Normalize ─────────────────────────────────────────────────────────
  message("\n[3/9] NormalizeData...")
  obj <- NormalizeData(obj, verbose = FALSE)   # RNA log-normalization

  # Re-run CLR on the T cell subset.
  # REQUIRED: CLR is compositional — each marker's value is computed relative
  # to the geometric mean of ALL markers in that cell. Values from the full
  # object reflect a whole-blood composition (myeloid, B cell, NK markers
  # present) and are no longer valid after subsetting to T cells only.
  if ("ADT" %in% names(obj@assays)) {
    message("  Re-normalizing ADT (CLR) on T cell subset...")
    obj <- NormalizeData(obj, assay = "ADT",
                         normalization.method = "CLR",
                         margin = 2,          # way to think about it - is this cell high for CD3 compared to other cells in my sample?
                         verbose = FALSE)
    message("  ADT CLR normalization complete")
  }

  # ── 4. HVGs ──────────────────────────────────────────────────────────────
  message("\n[4/9] FindVariableFeatures (n=", cfg$n_hvg, ")...")
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = cfg$n_hvg, verbose = FALSE)
  all_hvg     <- VariableFeatures(obj)
  cleaned_hvg <- filter_hvg(all_hvg, cfg$exclude_patterns)
  VariableFeatures(obj) <- cleaned_hvg
  message("  HVGs: ", length(all_hvg), " → ", length(cleaned_hvg),
          " after removing MT/ribo/TCR/Ig")
  write.csv(data.frame(gene = cleaned_hvg),
            file.path(data_dir, "hvg_list.csv"), row.names = FALSE)

  # ── 5. Cell cycle scoring ────────────────────────────────────────────────
  message("\n[5/9] Cell cycle scoring...")
  obj <- CellCycleScoring(obj,
    s.features   = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE)
  print(table(obj$Phase))

  regress_vars <- cfg$regress_vars
  if (cfg$regress_cc) regress_vars <- c(regress_vars, "S.Score", "G2M.Score")

  # ── 6. Scale + PCA ───────────────────────────────────────────────────────
  message("\n[6/9] ScaleData + RunPCA (regressing: ",
          paste(regress_vars, collapse = ", "), ")...")
  obj <- ScaleData(obj, features = cleaned_hvg,
                   vars.to.regress = regress_vars, verbose = FALSE)
  obj <- RunPCA(obj, features = cleaned_hvg, npcs = cfg$n_pcs, verbose = FALSE)

  p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs) +
    ggtitle("PCA Elbow Plot — T Cell Subset") + theme_cowplot(14)
  save_figure(p_elbow, "01_elbow_plot", plot_dir, 8, 5)

  # ── 7. UMAP + Clustering ─────────────────────────────────────────────────
  message("\n[7/9] UMAP + multi-resolution clustering (using PCA dims ",
          min(cfg$pca_dims), "–", max(cfg$pca_dims), ")...")
  obj <- RunUMAP(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = cfg$pca_dims, verbose = FALSE)

  for (res in cfg$resolutions) {
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    message("  res=", res, " → ",
            length(unique(obj@meta.data[[paste0(cfg$assay, "_snn_res.", res)]])),
            " clusters")
  }

  p_clustree <- clustree::clustree(obj@meta.data,
                                   prefix = paste0(cfg$assay, "_snn_res.")) +
    ggtitle("Clustree — Resolution Sweep") + theme(legend.position = "right")
  save_figure(p_clustree, "02_clustree", plot_dir, 12, 14)

  Idents(obj) <- paste0(cfg$assay, "_snn_res.", cfg$final_res)
  obj$seurat_clusters <- Idents(obj)
  message("  Active resolution: ", cfg$final_res, " → ",
          length(unique(Idents(obj))), " clusters")

  # ── Save clustered object checkpoint ─────────────────────────────────────
  clustered_rds <- file.path(cfg$out_dir, "subcluster_clustered.rds")
  saveRDS(obj, clustered_rds)
  message("  Clustered object saved: ", clustered_rds)

  # ── 8. UMAP visualizations ────────────────────────────────────────────────
  message("\n[8/9] UMAP plots...")

  p_clust <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle(paste0("T Cell Clusters (res=", cfg$final_res, ")")) +
    theme_cowplot(14) + NoLegend()
  save_figure(p_clust, "03_umap_clusters", plot_dir, cfg$umap_width, cfg$umap_height)

  parent_cols <- c("CD8+ T"="#2196A3", "Memory CD4+ T"="#3F51B5", "Naive CD4+ T"="#7B2D8B")
  p_parent <- DimPlot(obj, reduction = "umap", group.by = cfg$celltype_col,
                      cols = parent_cols, pt.size = 0.5) +
    ggtitle("Parent Cell Type") + theme_cowplot(14)
  save_figure(p_parent, "04_umap_parent_celltype", plot_dir,
              cfg$umap_width + 2, cfg$umap_height)

  status_cols <- c("Hi_Cog" = "#D95F5F", "Lo_Cog" = "#6BAED6")
  p_status <- DimPlot(obj, reduction = "umap", group.by = cfg$group_col,
                      cols = status_cols, pt.size = 0.5) +
    ggtitle("Cognitive Status (Hi_Cog vs Lo_Cog)") + theme_cowplot(14)
  save_figure(p_status, "05_umap_status", plot_dir, cfg$umap_width + 2, cfg$umap_height)

  donor_pal <- setNames(
    colorRampPalette(brewer.pal(8, "Set1"))(length(unique(obj@meta.data[[cfg$donor_col]]))),
    sort(unique(obj@meta.data[[cfg$donor_col]]))
  )
  p_donor <- DimPlot(obj, reduction = "umap", group.by = cfg$donor_col,
                     cols = donor_pal, pt.size = 0.4) +
    ggtitle("Donor") + theme_cowplot(14)
  save_figure(p_donor, "06_umap_donor", plot_dir, cfg$umap_width + 2, cfg$umap_height)

  p_phase <- DimPlot(obj, reduction = "umap", group.by = "Phase",
                     cols = c(G1="#A6A6A6", S="#F4A460", G2M="#DC143C"), pt.size = 0.4) +
    ggtitle("Cell Cycle Phase") + theme_cowplot(14)
  save_figure(p_phase, "07_umap_cell_cycle", plot_dir, cfg$umap_width, cfg$umap_height)

  key_features <- intersect(c("CD3D","CD4","CD8A","CCR7","SELL","GZMB",
                               "FOXP3","MKI67","LAG3","PDCD1","ISG15","CX3CR1"),
                             rownames(obj))
  p_feat <- FeaturePlot(obj, features = key_features, reduction = "umap",
                        ncol = 4, order = TRUE, pt.size = 0.3) & theme_cowplot(10)
  save_figure(p_feat, "08_feature_key_markers", plot_dir, 18, 10)

  # ── 8b. QC metric UMAPs (pale yellow → orange → deep red) ────────────────
  message("  QC metric UMAPs...")

  qc_scale <- scale_color_gradientn(
    colours = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C",
                "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"),
    na.value = "grey90"
  )

  qc_umap_features <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
    colnames(obj@meta.data)
  )

  qc_titles <- c(
    nCount_RNA   = "UMI Count (nCount_RNA)",
    nFeature_RNA = "Genes Detected (nFeature_RNA)",
    percent.mt   = "% Mitochondrial",
    percent.ribo = "% Ribosomal",
    nCount_Ribo  = "Ribosomal UMI Count"
  )

  if (length(qc_umap_features) > 0) {
    qc_plts <- lapply(qc_umap_features, function(feat) {
      FeaturePlot(obj, features = feat, reduction = "umap",
                  order = TRUE, pt.size = 0.3) +
        qc_scale +
        ggtitle(qc_titles[feat]) +
        theme_cowplot(12) +
        theme(
          legend.position  = "right",
          plot.title       = element_text(size = 13, face = "bold"),
          legend.key.height = unit(0.9, "cm"),
          legend.key.width  = unit(0.35, "cm")
        )
    })

    ncols_qc <- min(3L, length(qc_plts))
    nrows_qc <- ceiling(length(qc_plts) / ncols_qc)
    p_qc_umap <- wrap_plots(qc_plts, ncol = ncols_qc) +
      plot_annotation(
        title   = "QC Metrics on UMAP",
        theme   = theme(plot.title = element_text(face = "bold", size = 16))
      )
    save_figure(p_qc_umap, "08b_umap_qc_metrics", plot_dir,
                width  = ncols_qc * 6,
                height = nrows_qc * 5.5)
  }

  # ── 9. Marker visualizations + DE ─────────────────────────────────────────
  message("\n[9/9] Markers, proportions, DESeq2...")

  plot_genes <- intersect(unique(unlist(TCELL_MARKERS)), rownames(obj))
  p_dot <- DotPlot(obj, features = plot_genes, group.by = "seurat_clusters",
                   dot.scale = 6, col.min = -2, col.max = 2) +
    scale_color_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
    coord_flip() + theme_cowplot(11) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(size=8)) +
    labs(title="Canonical T Cell Markers by Cluster", x=NULL, y="Cluster")
  save_figure(p_dot, "09_dotplot_canonical_markers", plot_dir,
              cfg$dot_width, cfg$dot_height + 4)

  # Cluster marker genes
  message("  FindAllMarkers...")
  markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.1,
                            logfc.threshold=0.25, test.use="wilcox", verbose=FALSE)
  markers_sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, ]
  write.csv(markers_sig, file.path(data_dir, "cluster_markers_significant.csv"), row.names=FALSE)
  write.csv(markers,     file.path(data_dir, "cluster_markers_all.csv"),         row.names=FALSE)
  message("  Significant markers: ", nrow(markers_sig))

  top5 <- markers_sig %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n=5) %>% pull(gene) %>% unique()
  if (length(top5) > 0 && length(top5) <= 200) {
    p_heat <- DoHeatmap(obj, features=top5, group.by="seurat_clusters",
                        size=3, angle=45) +
      scale_fill_gradientn(colors=c("#2166AC","white","#B2182B")) +
      theme(axis.text.y=element_text(size=6))
    save_figure(p_heat, "10_heatmap_top5_markers", plot_dir, 18, 12)
  }

  # Proportion bar by donor + Status
  # Pure base-R: extract to plain atomic vectors BEFORE any dplyr/S4Vectors
  # method dispatch can intercept them. table() + tapply() are immune to the
  # Rle-column issue that kills dplyr::count() in a Seurat v5 context.
  .donor   <- as.vector(as.character(obj@meta.data[[cfg$donor_col]]))
  .group   <- as.vector(as.character(obj@meta.data[[cfg$group_col]]))
  .cluster <- as.vector(as.character(obj$seurat_clusters))

  prop_df <- as.data.frame(
    table(setNames(list(.donor, .group, .cluster),
                   c(cfg$donor_col, cfg$group_col, "seurat_clusters"))),
    stringsAsFactors = FALSE
  )
  prop_df <- prop_df[prop_df$Freq > 0, ]
  names(prop_df)[names(prop_df) == "Freq"] <- "n"
  donor_totals  <- tapply(prop_df$n, prop_df[[cfg$donor_col]], sum)
  prop_df$prop  <- prop_df$n / donor_totals[prop_df[[cfg$donor_col]]]
  prop_df$seurat_clusters <- factor(prop_df$seurat_clusters)

  p_prop <- ggplot(prop_df,
                   aes(x=.data[[cfg$donor_col]], y=prop, fill=seurat_clusters)) +
    geom_bar(stat="identity", position="stack", width=0.8) +
    facet_grid(~ .data[[cfg$group_col]], scales="free_x", space="free_x") +
    labs(title="T Cell Subcluster Proportions by Donor",
         x="Donor", y="Proportion", fill="Cluster") +
    theme_cowplot(13) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  save_figure(p_prop, "11_proportions_by_donor", plot_dir, 14, 7)

  # ── Propeller ─────────────────────────────────────────────────────────────
  message("  Propeller proportion test (Hi_Cog vs Lo_Cog)...")
  prop_res <- tryCatch(
    run_tcell_propeller(obj,
      cluster_col = "seurat_clusters",
      donor_col   = cfg$donor_col,
      group_col   = cfg$group_col,
      hi_group    = cfg$hi_group,
      lo_group    = cfg$lo_group),
    error = function(e) { message("  Propeller failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(prop_res)) {
    write.csv(prop_res, file.path(data_dir, "propeller_HiCog_vs_LoCog.csv"))
    message("  Propeller complete:")
    print(head(prop_res))
  }

  # ── MASC ──────────────────────────────────────────────────────────────────
  message("  MASC proportion test (mixed-effects, Hi_Cog vs Lo_Cog)...")
  masc_res <- tryCatch(
    run_masc(obj,
      cluster_col = "seurat_clusters",
      donor_col   = cfg$donor_col,
      group_col   = cfg$group_col,
      hi_group    = cfg$hi_group,
      lo_group    = cfg$lo_group),
    error = function(e) { message("  MASC failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(masc_res)) {
    write.csv(masc_res, file.path(data_dir, "MASC_HiCog_vs_LoCog.csv"), row.names=FALSE)
    message("  MASC complete:")
    print(masc_res)

    p_masc <- ggplot(masc_res[!is.na(masc_res$p_value), ],
                     aes(x=reorder(cluster, estimate), y=estimate,
                         fill=direction, alpha=p_adj < 0.05)) +
      geom_col() +
      geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se), width=0.3) +
      scale_fill_manual(values=c("Hi_Cog"="#D95F5F","Lo_Cog"="#6BAED6")) +
      scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.4), name="padj<0.05") +
      coord_flip() +
      labs(title="MASC: T Cell Subcluster Abundance\n(Hi_Cog vs Lo_Cog)",
           x="Cluster", y="Log-odds (Hi_Cog vs Lo_Cog)", fill="Enriched in") +
      theme_cowplot(13)
    save_figure(p_masc, "12_MASC_cluster_abundance", plot_dir, 10, 7)
  }

  # ── Pseudobulk DESeq2 ─────────────────────────────────────────────────────
  message("  Pseudobulk DESeq2 (Hi_Cog vs Lo_Cog) per cluster...")
  de_results <- pseudobulk_tcell_deseq2(obj,
    cluster_col = "seurat_clusters",
    donor_col   = cfg$donor_col,
    group_col   = cfg$group_col,
    hi_group    = cfg$hi_group,
    lo_group    = cfg$lo_group)
  if (length(de_results) > 0) {
    de_all <- do.call(rbind, de_results)
    write.csv(de_all,
              file.path(data_dir, "pseudobulk_DESeq2_HiCog_vs_LoCog.csv"),
              row.names=FALSE)
    message("  DESeq2 results saved (", length(de_results), " clusters)")
  }

  # ── Save object ───────────────────────────────────────────────────────────
  out_rds <- file.path(cfg$out_dir, "subcluster_final.rds")
  saveRDS(obj, out_rds)
  message("\n  Object saved: ", out_rds)
  message("\nPipeline complete. Outputs in:\n  ", cfg$out_dir)
  message(
    "\nNEXT STEPS:\n",
    "  1. 01_elbow_plot.png — choose cfg$pca_dims cutoff\n",
    "  2. 02_clustree.png   — choose cfg$final_res\n",
    "  3. 09_dotplot_canonical_markers.png — annotate clusters\n",
    "  4. 12_MASC_cluster_abundance.png — which clusters differ by cog status"
  )

  invisible(list(obj=obj, markers=markers_sig,
                 de_results=de_results, masc=masc_res, propeller=prop_res))
}

# =============================================================================
##  RUN
# =============================================================================
res <- run_tcell_pipeline(cfg)
