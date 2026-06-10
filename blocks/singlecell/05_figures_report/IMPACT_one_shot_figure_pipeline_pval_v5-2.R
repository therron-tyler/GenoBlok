
#!/usr/bin/env Rscript
# .libPaths(c("/path/to/your/R_library", .libPaths()))
# ============================================================
# IMPACT one-shot figure + PowerPoint generation pipeline
# Author: Tyler Therron
# Purpose:
#   Reproducibly generate IMPACT figures from a Seurat object and
#   optionally place them into the uploaded PowerPoint template.
#
# Main outputs:
#   - Figures/*.png
#   - Tables/*.csv
#   - IMPACT_auto_filled.pptx (optional)
#   - No analysis RDS is written by default
#
# Notes:
#   - This script is designed to be adapted once your final RDS arrives.
#   - It assumes sample IDs live in metadata column: HTO_maxID
#   - Pain group lives in metadata column: PainGroup
#   - Pool/hash run lives in metadata column: hash_run
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(colorspace)
  library(scales)
  library(forcats)
  # library(officer)
  # library(rvg)
  library(glue)
  library(limma)
})

# Required but loaded on demand: DESeq2, fgsea, msigdbr, UpSetR, ggVennDiagram

# ---------------------------- #
# 1) USER CONFIGURATION
# ---------------------------- #

cfg <- list(
  project_name = "IMPACT",
  input_rds    = "/path/to/data/WHLmerge_noIMP12_ManuscriptLabel_nCountRibo150_nCountRNA2000_dims1-9_11-12_res03.rds",
  output_dir   = "IMPACT_Figure_Output",
  covariate_csv = "IMPACT_Covariate_DEseq2_Table.csv",   # Age + Sex per sample
  template_pptx = "IMPACT figure template.pptx",   # optional
  make_pptx     = FALSE,

  # Standalone deck built from a directory of pre-rendered PNGs
  # (one image per slide, plain white background, no template required)
  plots_input_dir = "plots",
  plots_pptx_path = "analysis_plots_slides.pptx",
  make_plots_pptx = TRUE,

  # metadata columns
  meta = list(
    sample_col   = "HTO_maxID",
    pool_col     = "hash_run",
    pain_col     = "PainGroup",
    celltype_col = "CellType"
  ),

  # embeddings / assays
  reduction_umap = "umap",
  rna_assay      = "RNA",
  adt_assay      = "ADT",
  hto_assay      = "HTO",
  count_layer    = "counts",
  low_pain = c(
    "IMP38","IMP35","IMP29",
    "IMP11","IMP13","IMP21","IMP24","IMP37","IMP42","IMP30","IMP34"
  ),
  high_pain = c(
    "IMP16","IMP18","IMP19","IMP33","IMP36","IMP39","IMP46",
    "IMP25","IMP23","IMP15"
  ),

  # ADT markers: assay_name → display_name
  adt_markers = c(
    "Hu.CD3-UCHT1" = "CD3",
    "Hu.CD4-RPA.T4"   = "CD4",
    "Hu.CD8"          = "CD8",
    "Hu.CD335"        = "NKp46",
    "Hu.CD20-2H7"     = "CD20",
    "Hu.CD14-M5E2"    = "CD14",
    "Hu.CD16"         = "CD16",
    "Hu.CD1c"         = "CD1c",
    "Hu.CD303"        = "CD303"
  ),

  # canonical RNA markers for cell types (1-column violin layout)
  canonical_rna_markers = c(
    "IL7R",      # Naive/Memory CD4+ T
    "CCR7",      # Naive T
    "ICOS",      # CD8 T
    "CXCR3",      # Memory CD4+ T
    "IFIT2",     # Interferon response
    "NKG7",      # NK
    "MS4A1",     # B cells
    "CD14",       # CD14+ Monocyte
    "FCGR3A",    # FCGR3A+ Monocyte
    "CD1C",      # cDC
    "MZB1"       # Plasma B (if present)
  ),

  # highlight genes for pain group comparison across cell types
  # highlight_genes = c(
  #   "IFIT2","IFIT3","OASL","HERC5","PARP9","IFITM1","CMPK2",
  #   "PDCD1","SLC7A5","DUSP8",
  #   "KIR3DL1","PWP2",
  #   "TLE1","FCER1A","CD83",
  #   "IGHV1-69D","WARS","IGKC","BLK",
  #   "LMO4","A2M",
  #   "AREG","RASGRP2","CEBPD"
  # ),
  # highlight_genes = c(
  #   
  #   # IFN suppression in High Pain
  #   "IFIT2", "OASL", "CMPK2",
  #   
  #   # CD8+ T activation axis — TNF-α/NF-κB leading edge
  #   "CEBPD",    # Top of CD8+ T TNF-α NF-κB leading edge
  #   "TNFAIP3",  # CD8+ T TNF-α, Hypoxia, KRAS — NF-κB negative feedback; pain-relevant
  #   "DUSP4",    # CD8+ T TNF-α leading edge 
  #   
  #   # NK divergence
  #   "SOCS3",    # In IFN-γ leading edges of 5 cell types AND NK TNF-α — direction flip
  #   "JAG1",     # NK WNT/β-catenin leading edge; Notch-WNT crosstalk in tolerogenic NK
  #   "AXIN1"     # NK WNT/β-catenin leading edge; Notch-WNT crosstalk in tolerogenic NK
  # ),
  # highlight_genes = c(
  #   "PF4",
  #   "VEGFA",
  #   "VCAN",
  #   "TIMP1",
  #   "VAV1",
  #   "TNFAIP3",
  #   "SMAD3",
  #   "IL7R",
  #   "IRS2",
  #   "RELB",
  #   "JAG1",
  #   "PPARD",
  #   "DLL1",
  #   "AXIN1",
  #   "IFIT2",
  #   "IFIT3",
  #   "OASL",
  #   "OAS2"
  # ),
  
  
  
  
  highlight_genes = c(
    "IFIT3", # IFN
    "PF4", # angiogenesis
    "TNFAIP3", # TNF
    "JAG1", #WNT
    "OASL", # IFN
    "OLR1", # angiogenesis
    "SMAD3", # TNF
    "DLL1", #WNT
    "OAS2", # IFN
    "MEF2D", # angiogenesis
    "CEBPD", # TNF
    "AXIN1" #WNT
  ),

  # cell type colors — must match obj$CellType exactly
  # Order here defines display order for all figures
  cell_cols = c(
    "Naive CD4+ T"     = "#40B9B1",
    "Memory CD4+ T"    = "#00E3C6",
    "Interferon CD4+"  = "#47D5FE",
    "CD8+ T"           = "#0000FF",
    "Interferon CD8+"  = "#719AFF",
    "NK"               = "#A020F0",
    "B"                = "forestgreen",
    "CM"   = "#DC0000",
    "NCM" = "red4",
    "cDC"               = "coral1",
    "pDC"              = "#FFC700"
  ),

  # 6-celltype merged grouping for pseudobulk DE / GSEA / Venn
  merged_celltypes = list(
    "CD4+ T"    = c("Naive CD4+ T", "Memory CD4+ T", "Interferon CD4+"),
    "CD8+ T"    = c("CD8+ T", "Interferon CD8+"),
    "NK"        = "NK",
    "B"         = "B",
    "Monocytes" = c("CM", "NCM"),
    "DC"        = c("cDC", "pDC")
  ),

  # Colors for the 6 merged cell types
  merged_cell_cols = c(
    "CD4+ T"    = "#40B9B1",
    "CD8+ T"    = "#0000FF",
    "NK"        = "#A020F0",
    "B"         = "forestgreen",
    "Monocytes" = "#DC0000",
    "DC"        = "coral1"
  ),

  # Pool grouping: technical replicate pools collapsed
  pool_groups = list(
    "WHL1"   = c("IMP38","IMP35","IMP29","IMP25","IMP23","IMP15"),
    "WHL2/3" = c("IMP11","IMP12","IMP19","IMP21","IMP30","IMP33","IMP37","IMP39"),
    "WHL4/5" = c("IMP13","IMP16","IMP18","IMP24","IMP34","IMP36","IMP42","IMP46")
  ),
  pool_cols = c("WHL1" = "#E41A1C", "WHL2/3" = "#377EB8", "WHL4/5" = "#4DAF4A"),

  # DE / pseudobulk settings
  de = list(
    min_cells_per_sample = 20,
    pval_cutoff = 0.05,
    lfc_cutoff  = 0.585,  # log2(1.5) — 1.5-fold change
    top_n_labels = 12
  ),

  # powerpoint placement:
  # Adjust once after first run if needed.
  ppt = list(
    width = 13.333,
    height = 7.5,
    placements = list(
      slide1 = list(
        fig3A = c(left=0.40, top=1.20, width=4.00, height=2.35),
        fig3B = c(left=4.70, top=1.20, width=4.00, height=2.35),
        fig3C = c(left=8.95, top=1.20, width=3.95, height=2.35),
        fig3D = c(left=0.40, top=4.05, width=4.00, height=2.55),
        fig3E = c(left=4.70, top=4.05, width=8.20, height=2.55)
      ),
      slide2 = list(
        sfig3A = c(left=0.40, top=1.20, width=4.00, height=2.35),
        sfig3B = c(left=4.70, top=1.20, width=4.00, height=2.35),
        sfig3C = c(left=8.95, top=1.20, width=3.95, height=2.35),
        sfig3D = c(left=0.40, top=4.05, width=4.00, height=2.55),
        sfig3E = c(left=4.70, top=4.05, width=8.20, height=2.55)
      ),
      slide3 = list(
        fig4A = c(left=0.40, top=1.20, width=4.00, height=5.25),
        fig4B = c(left=4.70, top=1.20, width=4.00, height=5.25),
        fig4C = c(left=8.95, top=1.20, width=3.95, height=5.25)
      ),
      slide4 = list(
        sfig4A = c(left=0.40, top=1.20, width=6.10, height=5.25),
        sfig4B = c(left=6.80, top=1.20, width=6.10, height=5.25)
      )
    )
  )
)

# ---------------------------- #
# 2) UTILITIES
# ---------------------------- #

dir_create2 <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

save_plot_dual <- function(plot, filename_base, out_dir, width = 8, height = 6, dpi = 300) {
  dir_create2(out_dir)
  png_file <- file.path(out_dir, paste0(filename_base, ".png"))
  ggsave(png_file, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(list(png = png_file))
}

message2 <- function(...) cat(glue(..., .envir = parent.frame()), "\n")

# Helper: italic gene title for ggplot
gene_title <- function(g) bquote(italic(.(g)))

# Helper: assign pool group from sample ID
assign_pool_group <- function(sample_ids, cfg) {
  pool <- rep(NA_character_, length(sample_ids))
  for (grp in names(cfg$pool_groups)) {
    pool[sample_ids %in% cfg$pool_groups[[grp]]] <- grp
  }
  pool
}

validate_inputs <- function(obj, cfg) {
  # CellType is allowed to be missing — it will need to be added before running
  required_meta <- c(cfg$meta$sample_col, cfg$meta$pool_col)
  missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
  if (length(missing_meta) > 0) {
    stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
  }

  # CellType and PainGroup may need to be added manually before pipeline
  ct_col <- cfg$meta$celltype_col
  if (!(ct_col %in% colnames(obj@meta.data))) {
    stop("CellType column '", ct_col, "' not found in metadata. ",
         "You must assign cell type annotations before running this pipeline.\n",
         "  Example: obj$CellType <- Idents(obj)  # or your own mapping")
  }

  pain_col <- cfg$meta$pain_col
  if (!(pain_col %in% colnames(obj@meta.data))) {
    message("PainGroup column '", pain_col, "' not in metadata — will be created from low/high sample lists.")
  }

  if (!(cfg$reduction_umap %in% names(obj@reductions))) {
    stop("Reduction not found: ", cfg$reduction_umap)
  }

  if (!(cfg$rna_assay %in% names(obj@assays))) {
    stop("RNA assay not found: ", cfg$rna_assay)
  }

  if (!(cfg$adt_assay %in% names(obj@assays))) {
    warning("ADT assay '", cfg$adt_assay, "' not found. ADT FeaturePlots will be skipped.")
  }
}

standardize_metadata <- function(obj, cfg) {
  md <- obj@meta.data
  sample_col   <- cfg$meta$sample_col
  pain_col     <- cfg$meta$pain_col
  celltype_col <- cfg$meta$celltype_col

  md[[sample_col]] <- as.character(md[[sample_col]])
  md[[celltype_col]] <- as.character(md[[celltype_col]])

  # Create / overwrite PainGroup from provided sample lists
  md[[pain_col]] <- case_when(
    md[[sample_col]] %in% cfg$low_pain  ~ "Low Pain",
    md[[sample_col]] %in% cfg$high_pain ~ "High Pain",
    TRUE ~ NA_character_
  )
  md[[pain_col]] <- factor(md[[pain_col]], levels = c("Low Pain", "High Pain"))

  n_assigned <- sum(!is.na(md[[pain_col]]))
  n_total    <- nrow(md)
  message(glue("PainGroup assigned: {n_assigned}/{n_total} cells ",
               "({sum(md[[pain_col]] == 'Low Pain', na.rm=TRUE)} Low, ",
               "{sum(md[[pain_col]] == 'High Pain', na.rm=TRUE)} High)"))
  if (n_assigned == 0) {
    warning("No cells matched the low/high pain sample lists. ",
            "Check that HTO_maxID values match cfg$low_pain / cfg$high_pain.")
  }

  # enforce cell type order for all downstream plots
  present_levels <- intersect(names(cfg$cell_cols), unique(md[[celltype_col]]))
  extra_levels <- setdiff(unique(md[[celltype_col]]), names(cfg$cell_cols))
  md[[celltype_col]] <- factor(md[[celltype_col]], levels = c(present_levels, extra_levels))

  obj@meta.data <- md
  obj
}

plot_umap_celltypes <- function(obj, cfg) {
  ct_col <- cfg$meta$celltype_col
  umap_cols <- colnames(Embeddings(obj, cfg$reduction_umap))
  emb <- Embeddings(obj, cfg$reduction_umap) %>%
    as.data.frame() %>%
    mutate(CellType = obj@meta.data[[ct_col]]) %>%
    group_by(CellType) %>%
    summarize(across(all_of(umap_cols), median, na.rm = TRUE), .groups = "drop")

  p <- DimPlot(
    obj,
    reduction = cfg$reduction_umap,
    group.by = ct_col,
    cols = cfg$cell_cols,
    raster = FALSE
  ) +
#    geom_text_repel(
#      data = emb,
#      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], label = CellType),
#      size = 4,
#      segment.color = "grey30",
#      box.padding = 0.4,
#      point.padding = 0.25,
#      seed = 1
#    ) +
    labs(title = "Cell Type", x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 22) +
    theme(
#      legend.position = "none",
      plot.title = element_text(face = "bold", size = 23)
    )
  p
}

plot_umap_pain <- function(obj, cfg) {
  pain_col <- cfg$meta$pain_col
  pain_cols <- c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")
  DimPlot(
    obj,
    reduction = cfg$reduction_umap,
    group.by = pain_col,
    cols = pain_cols,
    raster = FALSE
  ) +
    labs(title = "Pain Group", x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 22) +
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_umap_pool <- function(obj, cfg) {
  pool_col <- cfg$meta$pool_col
  DimPlot(
    obj,
    reduction = cfg$reduction_umap,
    group.by = pool_col,
    raster = FALSE
  ) +
    labs(title = "UMAP by Pool", x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 22) +
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_umap_sample <- function(obj, cfg) {
  sample_col <- cfg$meta$sample_col
  DimPlot(
    obj,
    reduction = cfg$reduction_umap,
    group.by = sample_col,
    raster = FALSE,
    label = FALSE
  ) +
    labs(title = "UMAP by Sample", x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 22) +
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_qc_by_group <- function(obj, group_col, title = "QC", ncol = 2) {
  qc_features <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"), colnames(obj@meta.data))
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
        plot.title = element_text(face = "bold", size = 22),
        axis.title.x = element_blank(),
        axis.text.x = if (i < n) element_blank() else element_text(angle = 45, hjust = 1, size = 19),
        axis.ticks.x = if (i < n) element_blank() else element_line()
      )
  }

  wrap_plots(plist, ncol = 1) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 23))
}

plot_qc_by_sample_pooled <- function(obj, cfg, title = "Single-cell QC by Sample") {
  qc_features <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_Ribo"),
                            colnames(obj@meta.data))
  if (length(qc_features) == 0) return(NULL)

  sample_col <- cfg$meta$sample_col
  pool_col   <- cfg$meta$pool_col   # hash_run — used directly so no sample is lost

  # Pull sample -> pool mapping directly from metadata (no static lookup table).
  # Use the most-frequent hash_run per sample to guard against rare cells where
  # a sample ID appears with more than one pool tag (Seurat v5 Rle artefact or
  # genuine multi-pool overlap), which would produce duplicate factor levels.
  md <- as.data.frame(obj@meta.data[, c(sample_col, pool_col), drop = FALSE])
  md[] <- lapply(md, as.character)
  colnames(md) <- c("Sample", "PoolGroup")

  pool_df <- md %>%
    count(Sample, PoolGroup) %>%                          # tally each Sample × PoolGroup combo
    group_by(Sample) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>% # keep majority pool per sample
    ungroup() %>%
    select(Sample, PoolGroup) %>%
    mutate(PoolGroup = factor(PoolGroup, levels = names(cfg$pool_cols)))

  # Any pool not in cfg$pool_cols gets labelled Unknown rather than silently dropped
  pool_df$PoolGroup[is.na(pool_df$PoolGroup)] <- "Unknown"

  # Sort: pool group first, then sample ID alphabetically within pool
  pool_df <- pool_df[order(pool_df$PoolGroup, pool_df$Sample), ]
  sample_order <- pool_df$Sample

  # Sample-to-colour mapping derived from pool membership
  sample_cols <- setNames(
    ifelse(as.character(pool_df$PoolGroup) %in% names(cfg$pool_cols),
           cfg$pool_cols[as.character(pool_df$PoolGroup)],
           "grey60"),
    pool_df$Sample
  )

  # Apply ordering to the metadata factor before VlnPlot
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

plot_adt_feature_pages <- function(obj, cfg, per_page = 9, ncol = 3) {
  if (!(cfg$adt_assay %in% names(obj@assays))) return(NULL)

  DefaultAssay(obj) <- cfg$adt_assay
  assay_names <- names(cfg$adt_markers)
  feats <- assay_names[assay_names %in% rownames(obj[[cfg$adt_assay]])]
  if (length(feats) == 0) {
    warning("No requested ADT markers found in assay: ", cfg$adt_assay)
    return(NULL)
  }

  chunks <- split(feats, ceiling(seq_along(feats) / per_page))

  pages <- lapply(seq_along(chunks), function(i) {
    chunk     <- chunks[[i]]
    n_panels  <- length(chunk)
    nrow_page <- ceiling(n_panels / ncol)

    # Bottom row: last ncol positions; left column: every ncol-th starting at 1
    bottom_idx <- seq(from = (nrow_page - 1) * ncol + 1, to = n_panels)
    left_idx   <- seq(from = 1, to = n_panels, by = ncol)

    # Extract shared colorbar from the first feature using Seurat's actual color scale
    # so the legend exactly matches the plot points
    first_f       <- chunk[[1]]
    first_display <- cfg$adt_markers[[first_f]]
    adt_limits <- c(0, 3)   # cap at 3 to keep scale consistent across markers

    legend_source <- FeaturePlot(obj, features = first_f,
                                 reduction = cfg$reduction_umap, raster = FALSE) +
      scale_color_gradientn(
        colours  = c("lightgrey", "#4B0082"),
        limits   = adt_limits,
        oob      = scales::squish,
        name     = "Surface\nMarker\nIntensity",
        breaks   = adt_limits,
        labels   = as.character(adt_limits)
      ) +
      guides(color = guide_colorbar(
        title.position = "top",
        barwidth       = unit(0.5, "cm"),
        barheight      = unit(5,   "cm"),
        ticks          = FALSE
      )) +
      theme(
        legend.title         = element_text(size = 16, face = "bold"),
        legend.text          = element_text(size = 14),
        legend.position      = "right",
        legend.justification = "center"
      )
    shared_cbar <- cowplot::get_legend(legend_source)

    # Build all panels with per-panel legends suppressed
    plist <- lapply(seq_along(chunk), function(j) {
      f         <- chunk[[j]]
      display   <- cfg$adt_markers[[f]]
      is_bottom <- j %in% bottom_idx
      is_left   <- j %in% left_idx

      FeaturePlot(obj, features = f, reduction = cfg$reduction_umap, raster = FALSE) +
        scale_color_gradientn(
          colours = c("lightgrey", "#4B0082"),
          limits  = adt_limits,
          oob     = scales::squish
        ) +
        ggtitle(display) +
        theme(
          plot.title   = element_text(face = "bold", size = 23),
          axis.title.x = if (is_bottom) element_text(size = 23) else element_blank(),
          axis.text.x  = if (is_bottom) element_text(size = 21) else element_blank(),
          axis.ticks.x = if (is_bottom) element_line()          else element_blank(),
          axis.title.y = if (is_left)   element_text(size = 23) else element_blank(),
          axis.text.y  = if (is_left)   element_text(size = 21) else element_blank(),
          axis.ticks.y = if (is_left)   element_line()          else element_blank(),
          legend.position = "none"
        )
    })

    feature_grid <- wrap_plots(plist, ncol = ncol) +
      plot_annotation(title = paste0("ADT Surface Protein Markers (", i, "/", length(chunks), ")"))

    # Attach shared colorbar; rel_widths gives enough room for the full legend title
    cowplot::plot_grid(feature_grid, shared_cbar,
                       ncol = 2, rel_widths = c(1, 0.12))
  })

  DefaultAssay(obj) <- cfg$rna_assay
  pages
}

plot_adt_violin_pages <- function(obj, cfg, per_page = 9, ncol = 3) {
  if (!(cfg$adt_assay %in% names(obj@assays))) return(NULL)

  DefaultAssay(obj) <- cfg$adt_assay
  ct_col <- cfg$meta$celltype_col
  assay_names <- names(cfg$adt_markers)
  feats <- assay_names[assay_names %in% rownames(obj[[cfg$adt_assay]])]
  if (length(feats) == 0) {
    warning("No requested ADT markers found in assay: ", cfg$adt_assay)
    return(NULL)
  }

  # Enforce cell type order
  Idents(obj) <- factor(obj@meta.data[[ct_col]], levels = names(cfg$cell_cols))

  chunks <- split(feats, ceiling(seq_along(feats) / per_page))

  pages <- lapply(seq_along(chunks), function(i) {
    plist <- lapply(chunks[[i]], function(f) {
      display <- cfg$adt_markers[[f]]
      VlnPlot(obj, features = f, group.by = ct_col, pt.size = 0,
              cols = cfg$cell_cols) +
        ggtitle(display) +
        theme(
          plot.title = element_text(face = "bold", size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
          legend.position = "none"
        )
    })
    wrap_plots(plist, ncol = ncol) +
      plot_annotation(title = paste0("ADT Violin Plots (", i, "/", length(chunks), ")"))
  })

  DefaultAssay(obj) <- cfg$rna_assay
  pages
}

plot_canonical_rna_violins <- function(obj, cfg, ncol_layout = 1) {
  DefaultAssay(obj) <- cfg$rna_assay
  feats <- cfg$canonical_rna_markers[cfg$canonical_rna_markers %in% rownames(obj[[cfg$rna_assay]])]
  ct_col <- cfg$meta$celltype_col

  # Enforce cell type ordering
  Idents(obj) <- factor(obj@meta.data[[ct_col]], levels = names(cfg$cell_cols))

  plist <- VlnPlot(
    obj,
    features = feats,
    group.by = ct_col,
    pt.size = 0,
    combine = FALSE,
    cols = cfg$cell_cols
  )

  n <- length(plist)
  # Bottom row: last ncol_layout plots get x-axis labels
  bottom_row_start <- n - (n %% ncol_layout)
  if (bottom_row_start == n) bottom_row_start <- n - ncol_layout
  bottom_indices <- seq(bottom_row_start + 1, n)

  for (i in seq_along(plist)) {
    show_x <- i %in% bottom_indices
    plist[[i]] <- plist[[i]] +
      labs(title = bquote(bold(italic(.(feats[i]))))) +
      theme_classic(base_size = 23) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = if (show_x) element_text(angle = 45, hjust = 1, size = 23)
                       else element_blank(),
        axis.ticks.x = if (show_x) element_line() else element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 23),
        plot.margin = margin(3, 3, 3, 45)   # left margin keeps x-labels from clipping
      )
  }

  # Single shared y-axis label to the left of the entire stacked grid
  violin_grid <- wrap_plots(plist, ncol = ncol_layout)
  y_label     <- cowplot::ggdraw() +
    cowplot::draw_label("RNA Expression", angle = 90, size = 26, fontface = "bold",
                        x = 0.5, y = 0.65)
  cowplot::plot_grid(y_label, violin_grid, ncol = 2, rel_widths = c(0.04, 1))
}

celltype_sample_count_table <- function(obj, cfg) {
  ct_col <- cfg$meta$celltype_col
  sample_col <- cfg$meta$sample_col
  pain_col <- cfg$meta$pain_col

  obj@meta.data %>%
    count(.data[[sample_col]], .data[[pain_col]], .data[[ct_col]], name = "n") %>%
    group_by(.data[[sample_col]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    rename(Sample = 1, PainGroup = 2, CellType = 3)
}

plot_stacked_cell_composition_by_sample <- function(obj, cfg) {
  df <- celltype_sample_count_table(obj, cfg)
  ggplot(df, aes(x = Sample, y = prop, fill = CellType)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = cfg$cell_cols, drop = FALSE) +
    facet_grid(~ PainGroup, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "Cell Type Composition by Sample", x = NULL, y = "Proportion") +
    theme_classic(base_size = 23) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y     = element_text(size = 21),
      axis.title.y    = element_text(size = 22),
      plot.title      = element_text(face = "bold", size = 23),
      strip.text      = element_text(size = 16, face = "bold"),
      legend.text     = element_text(size = 20),
      legend.title    = element_text(size = 21),
      legend.position = "right"
    )
}

run_propeller_limma <- function(obj, cfg, covariates = NULL,
                                transform = c("asin", "logit"),
                                adjust.method = "fdr", eps = 1e-3) {
  transform <- match.arg(transform)
  ct_col <- cfg$meta$celltype_col
  sample_col <- cfg$meta$sample_col
  pain_col <- cfg$meta$pain_col

  # Convert to plain data.frame to avoid Seurat v5 Rle issues
  md <- as.data.frame(obj@meta.data[, c(sample_col, pain_col, ct_col), drop = FALSE])
  md[] <- lapply(md, as.character)
  colnames(md) <- c("Sample", "PainGroup", "CellType")
  md <- md[!is.na(md$PainGroup), ]

  counts <- table(md$CellType, md$Sample)
  rawP <- prop.table(counts, margin = 2)

  sample_order <- colnames(rawP)

  # Build sample-level metadata
  sample_meta <- md %>%
    distinct(Sample, PainGroup) %>%
    filter(Sample %in% sample_order) %>%
    arrange(match(Sample, sample_order))

  use_covariates <- !is.null(covariates) && all(c("Sex", "Age") %in% colnames(covariates))

  # Merge covariates if provided, then explicitly align dimensions
  if (use_covariates) {
    cov_slim <- covariates[, c("Sample", "Sex", "Age"), drop = FALSE]
    cov_slim <- cov_slim[!duplicated(cov_slim$Sample), ]   # guard against dup rows
    sample_meta <- merge(sample_meta, cov_slim, by = "Sample", all.x = TRUE)

    # Reorder to match sample_order, then drop samples missing Sex/Age
    sample_meta <- sample_meta[match(sample_order, sample_meta$Sample), ]
    complete_idx <- which(!is.na(sample_meta$Sex) & !is.na(sample_meta$Age) &
                            !is.na(sample_meta$Sample))
    if (length(complete_idx) < nrow(sample_meta)) {
      dropped <- sample_meta$Sample[setdiff(seq_len(nrow(sample_meta)), complete_idx)]
      message("  Propeller: dropping ", length(dropped),
              " sample(s) missing covariate data: ", paste(dropped, collapse = ", "))
      sample_meta  <- sample_meta[complete_idx, , drop = FALSE]
      sample_order <- sample_meta$Sample          # update order to match trimmed meta
      rawP         <- rawP[, sample_order, drop = FALSE]
    }
    sample_meta$Sex <- factor(sample_meta$Sex)
    sample_meta$Age <- as.numeric(sample_meta$Age)
  } else {
    sample_meta <- sample_meta[match(sample_order, sample_meta$Sample), ]
  }

  # Verify alignment before proceeding
  stopifnot(nrow(sample_meta) == length(sample_order),
            all(sample_meta$Sample == sample_order))

  grp <- factor(sample_meta$PainGroup, levels = c("Low Pain", "High Pain"))

  if (transform == "asin") {
    props <- asin(sqrt(rawP))
  } else {
    p2 <- (rawP + eps) / (1 + 2 * eps)
    props <- log(p2 / (1 - p2))
  }

  # Design matrix: include covariates if available
  if (use_covariates) {
    design <- model.matrix(~ Sex + Age + grp, data = sample_meta)
    message("  Propeller design: ~ Sex + Age + PainGroup")
    # PainGroup coefficient is the last column: "grpHigh Pain"
    coef_idx <- ncol(design)
  } else {
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- c("LowPain", "HighPain")
    message("  Propeller design: ~ 0 + PainGroup")
  }

  fit <- lmFit(props, design)

  if (!is.null(covariates) && all(c("Sex", "Age") %in% colnames(sample_meta))) {
    # With covariates, the PainGroup effect is already a coefficient
    fit2 <- eBayes(fit)
    tt <- topTable(fit2, coef = coef_idx, number = nrow(props),
                   adjust.method = adjust.method, sort.by = "none") %>%
      rownames_to_column("CellType")
  } else {
    # Without covariates, use contrast
    contrast <- makeContrasts(HighPain - LowPain, levels = design)
    fit2 <- contrasts.fit(fit, contrast) |> eBayes()
    tt <- topTable(fit2, coef = 1, number = nrow(props),
                   adjust.method = adjust.method, sort.by = "none") %>%
      rownames_to_column("CellType")
  }

  means <- data.frame(
    CellType = rownames(rawP),
    LowMean  = rowMeans(rawP[, grp == "Low Pain", drop = FALSE]),
    HighMean = rowMeans(rawP[, grp == "High Pain", drop = FALSE])
  )

  left_join(tt, means, by = "CellType")
}

plot_cell_composition_boxplots <- function(obj, cfg, prop_res = NULL, lighten_amount = 0.55) {
  df <- celltype_sample_count_table(obj, cfg)
  df$PainGroup <- factor(df$PainGroup, levels = c("Low Pain", "High Pain"))
  df$CellType  <- factor(df$CellType,  levels = names(cfg$cell_cols))

  dark_cols  <- cfg$cell_cols
  light_cols <- lighten(cfg$cell_cols, amount = lighten_amount)
  fill_map <- c(
    setNames(light_cols, paste0(names(light_cols), "::Low Pain")),
    setNames(dark_cols,  paste0(names(dark_cols), "::High Pain"))
  )

  df <- df %>%
    mutate(fill_key = factor(paste0(CellType, "::", PainGroup), levels = names(fill_map)))

  p <- ggplot(df, aes(x = CellType, y = prop, fill = fill_key)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
    geom_jitter(aes(group = PainGroup), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8), size = 1.8, alpha = 0.85) +
    scale_fill_manual(values = fill_map, guide = "none") +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "Cell Composition by Pain Group", x = NULL, y = "Per-sample fraction") +
    theme_classic(base_size = 21) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 23)
    )

  if (!is.null(prop_res) && nrow(prop_res) > 0) {
    # Show p-value labels for all cell types
    ann <- prop_res %>%
      mutate(
        CellType = factor(CellType, levels = names(cfg$cell_cols)),
        pval_label = case_when(
          P.Value < 0.001 ~ formatC(P.Value, format = "e", digits = 1),
          P.Value < 0.01  ~ paste0("p=", formatC(P.Value, format = "f", digits = 3)),
          P.Value < 0.1   ~ paste0("p=", formatC(P.Value, format = "f", digits = 2)),
          TRUE             ~ paste0("p=", formatC(P.Value, format = "f", digits = 2))
        ),
        star = case_when(
          P.Value < 0.001 ~ " ***",
          P.Value < 0.01  ~ " **",
          P.Value < 0.05  ~ " *",
          P.Value < 0.1   ~ " .",
          TRUE ~ ""
        ),
        label = paste0(pval_label, star)
      ) %>%
      filter(!is.na(CellType)) %>%
      # stagger y positions per cell type to avoid overlap
      mutate(y = 1.03 * max(df$prop, na.rm = TRUE))

    if (nrow(ann) > 0) {
      p <- p + geom_text(data = ann, aes(x = CellType, y = y, label = label),
                         inherit.aes = FALSE, size = 2.8, angle = 30, hjust = 0)
    }

    # Add subtitle with FDR note
    p <- p + labs(subtitle = "P-values: limma moderated t-test on arcsin(sqrt) transformed proportions\n(High Pain − Low Pain); '.' p<0.1, '*' p<0.05, '**' p<0.01, '***' p<0.001")
  }

  p
}

find_denovo_markers <- function(obj, cfg, top_n = 10) {
  DefaultAssay(obj) <- cfg$rna_assay
  ct_col <- cfg$meta$celltype_col

  old_id <- Idents(obj)
  Idents(obj) <- obj@meta.data[[ct_col]]

  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.2,
    verbose = FALSE
  )

  Idents(obj) <- old_id

  markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()
}

plot_denovo_marker_heatmap <- function(obj, cfg, markers_tbl, top_n = 8, max_cells_per_ct = 200) {
  DefaultAssay(obj) <- cfg$rna_assay

  top_features <- markers_tbl %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    pull(gene) %>%
    unique()

  # Set idents so DoHeatmap groups by cell type
  ct_col <- cfg$meta$celltype_col
  old_id <- Idents(obj)
  Idents(obj) <- factor(obj@meta.data[[ct_col]], levels = names(cfg$cell_cols))

  # Downsample to max_cells_per_ct per group to prevent PDF smearing
  set.seed(42)
  cells_keep <- obj@meta.data %>%
    as.data.frame() %>%
    mutate(.cell = rownames(obj@meta.data)) %>%
    group_by(.data[[ct_col]]) %>%
    slice_sample(prop = 1) %>%
    slice_head(n = max_cells_per_ct) %>%
    pull(.cell)
  sub <- subset(obj, cells = cells_keep)
  DefaultAssay(sub) <- cfg$rna_assay
  sub <- ScaleData(sub, features = top_features, verbose = FALSE)
  Idents(sub) <- factor(sub@meta.data[[ct_col]], levels = names(cfg$cell_cols))

  p <- DoHeatmap(
    sub,
    features = top_features,
    group.by = ct_col,
    group.colors = cfg$cell_cols,
    size = 5,
    raster = TRUE
  ) +
    scale_fill_gradientn(colors = c("purple4", "grey10", "yellow")) +
    labs(title = "De novo Marker Heatmap") +
    guides(color = "none") +   # remove cell type color bar; column annotation provides that info
    theme(
      plot.title        = element_text(face = "bold", size = 23),
      axis.text.y       = element_text(size = 21),
      legend.text       = element_text(size = 21),
      legend.title      = element_text(size = 22),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.5, "cm")
    )

  Idents(obj) <- old_id
  p
}

# DoHeatmap + ggsave = blank PNG; save with explicit device calls
save_heatmap_direct <- function(plot, filename_base, out_dir, width = 12, height = 14, dpi = 300) {
  dir_create2(out_dir)
  png_file <- file.path(out_dir, paste0(filename_base, ".png"))

  png(png_file, width = width, height = height, units = "in", res = dpi, type = "cairo")
  print(plot)
  invisible(dev.off())

  invisible(list(png = png_file))
}

pseudo_bulk_counts <- function(obj, cfg, celltype, covariates = NULL) {
  DefaultAssay(obj) <- cfg$rna_assay
  ct_col     <- cfg$meta$celltype_col
  sample_col <- cfg$meta$sample_col
  pain_col   <- cfg$meta$pain_col

  keep_cells <- rownames(obj@meta.data)[
    obj@meta.data[[ct_col]] == celltype &
      !is.na(obj@meta.data[[pain_col]])
  ]
  if (length(keep_cells) == 0) stop("No cells found for cell type: ", celltype)

  sub <- subset(obj, cells = keep_cells)
  md <- sub@meta.data

  sample_sizes <- table(md[[sample_col]])
  valid_samples <- names(sample_sizes)[sample_sizes >= cfg$de$min_cells_per_sample]
  if (length(valid_samples) < 4) {
    stop("Too few valid samples (", length(valid_samples), "/", length(sample_sizes),
         " with >=", cfg$de$min_cells_per_sample, " cells) for ", celltype)
  }

  sub <- subset(sub, cells = rownames(md)[md[[sample_col]] %in% valid_samples])
  md <- sub@meta.data

  counts <- GetAssayData(sub, assay = cfg$rna_assay, layer = cfg$count_layer)
  sample_ids <- md[[sample_col]]

  pb <- sapply(split(seq_along(sample_ids), sample_ids), function(idx) {
    Matrix::rowSums(counts[, idx, drop = FALSE])
  })

  pb <- as.matrix(pb)

  # One row per sample — convert to plain data.frame first (Seurat v5 Rle fix)
  pool_col <- cfg$meta$pool_col
  coldata <- as.data.frame(md[, c(sample_col, pain_col, pool_col), drop = FALSE])
  coldata[] <- lapply(coldata, as.character)  # strip Rle
  coldata <- coldata[!duplicated(coldata[[sample_col]]), , drop = FALSE]
  coldata <- coldata[coldata[[sample_col]] %in% colnames(pb), , drop = FALSE]
  coldata <- coldata[match(colnames(pb), coldata[[sample_col]]), , drop = FALSE]
  rownames(coldata) <- coldata[[sample_col]]

  # Merge covariates (Age, Sex) if provided — only bring in Sex and Age
  if (!is.null(covariates) && all(c("Sex", "Age") %in% colnames(covariates))) {
    cov_slim <- covariates[, c("Sample", "Sex", "Age"), drop = FALSE]
    coldata <- merge(coldata, cov_slim, by.x = sample_col, by.y = "Sample", all.x = TRUE)
    rownames(coldata) <- coldata[[sample_col]]
    coldata <- coldata[match(colnames(pb), coldata[[sample_col]]), , drop = FALSE]
  }

  list(counts = pb, coldata = coldata)
}

run_deseq_one_celltype <- function(obj, cfg, celltype, covariates = NULL) {
  suppressPackageStartupMessages(library(DESeq2))
  pb <- pseudo_bulk_counts(obj, cfg, celltype, covariates = covariates)

  pain_col <- cfg$meta$pain_col

  # ── Filter out uninformative genes before DESeq2 ──
  gene_names <- rownames(pb$counts)

  # Mitochondrial genes
  mt_genes <- grepl("^MT-", gene_names)
  # Ribosomal protein genes
  ribo_genes <- grepl("^RP[SL][0-9]", gene_names)
  # Y chromosome genes (curated GRCh38 list)
  y_chrom <- c(
    "AMELX", "AMELY", "BPY2", "BPY2B", "BPY2C", "CDY1", "CDY1B", "CDY2A",
    "CDY2B", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DDX3Y", "EIF1AY", "HSFY1",
    "HSFY2", "KDM5D", "LINC00278", "NLGN4Y", "PCDH11Y", "PRY", "PRY2",
    "RBMY1A1", "RBMY1B", "RBMY1D", "RBMY1E", "RBMY1F", "RBMY1J",
    "RPS4Y1", "RPS4Y2", "SRY", "TBL1Y", "TMSB4Y", "TTTY1", "TTTY2",
    "TTTY10", "TTTY14", "TTTY15", "TTTY21", "TTTY22", "TXLNGY",
    "USP9Y", "UTY", "VCY", "VCY1B", "XKRY", "XKRY2", "ZFY", "PRKY",
    "TSPY1", "TSPY2", "TSPY3", "TSPY4", "TSPY8", "TSPY9P", "TSPY10"
  )
  y_genes <- gene_names %in% y_chrom | grepl("^TTTY", gene_names)

  exclude <- mt_genes | ribo_genes | y_genes
  n_excluded <- sum(exclude)
  pb$counts <- pb$counts[!exclude, , drop = FALSE]
  message("    Excluded ", n_excluded, " genes (MT/ribo/chrY); ",
          nrow(pb$counts), " genes remaining")

  # ── Build design formula ──
  if (!is.null(covariates) && all(c("Sex", "Age") %in% colnames(pb$coldata))) {
    pb$coldata$Sex <- factor(pb$coldata$Sex)
    pb$coldata$Age <- as.numeric(pb$coldata$Age)
    design_formula <- reformulate(c("Sex", "Age", pain_col))
    message("    Design: ~ Sex + Age + ", pain_col)
  } else {
    design_formula <- reformulate(pain_col)
    message("    Design: ~ ", pain_col)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(pb$counts),
    colData   = pb$coldata,
    design    = design_formula
  )

  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  print(paste0("Hey man, theres like this many rows going into desheet2 analysis",nrow(dds)))
  dds[[pain_col]] <- relevel(factor(dds[[pain_col]]), ref = "Low Pain")
  dds <- DESeq(dds, quiet = TRUE)

  # positive log2FC = up in High Pain relative to Low Pain
  res <- results(dds, contrast = c(pain_col, "High Pain", "Low Pain")) %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    arrange(pvalue)

  list(dds = dds, res = res, pb = pb, celltype = celltype)
}

plot_volcano <- function(res_df, title = "Volcano", pval_cutoff = 0.05, lfc_cutoff = 1, top_n = 12) {
  df <- res_df %>%
    mutate(
      neglog10 = -log10(pvalue),
      `Expression Significance` = case_when(
        is.na(pvalue) ~ "NS",
        pvalue < pval_cutoff & log2FoldChange >=  lfc_cutoff ~ "High Pain",
        pvalue < pval_cutoff & log2FoldChange <= -lfc_cutoff ~ "Low Pain",
        TRUE ~ "NS"
      )
    )

  top_df <- df %>%
    filter(`Expression Significance` != "NS") %>%
    arrange(pvalue) %>%
    slice_head(n = top_n)

  ggplot(df, aes(x = log2FoldChange, y = neglog10, color = `Expression Significance`)) +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = 2) +
    geom_text_repel(
      data = top_df,
      aes(label = Gene),
      size = 4.5,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c("High Pain" = "#D95F5F", "Low Pain" = "#6BAED6", "NS" = "grey75"),
                       name = "Expression\nSignificance") +
    guides(color = guide_legend(override.aes = list(size = 8, label = ""))) +
    labs(title = title, x = "log2 fold-change (High/Low Pain)", y = expression(-log[10](p))) +
    theme_classic(base_size = 21) +
    theme(plot.title = element_text(face = "bold", size = 23),
          legend.text     = element_text(size = 23),
          legend.title    = element_text(size = 23))
}

# ── Combined volcano grid: all merged cell types, shared legend, 2:1 aspect ──

plot_volcano_grid <- function(de_results, cfg,
                               ct_order = c("CD4+ T", "CD8+ T", "NK",
                                            "B", "Monocytes", "DC"),
                               ncol = 3) {
  cts_use <- intersect(ct_order, names(de_results))
  if (length(cts_use) == 0) return(NULL)

  n_panels  <- length(cts_use)
  nrow_grid <- ceiling(n_panels / ncol)

  # Identify bottom row and left column positions
  bottom_idx <- seq(from = (nrow_grid - 1) * ncol + 1, to = n_panels)
  left_idx   <- seq(from = 1, to = n_panels, by = ncol)

  # Build one volcano per cell type with axis labels suppressed/shown per position
  plist <- lapply(seq_along(cts_use), function(j) {
    ct         <- cts_use[j]
    is_bottom  <- j %in% bottom_idx
    is_left    <- j %in% left_idx
    ct_color   <- if (!is.null(cfg$merged_cell_cols) && ct %in% names(cfg$merged_cell_cols))
                    cfg$merged_cell_cols[[ct]] else "black"

    plot_volcano(
      de_results[[ct]]$res,
      title       = ct,
      pval_cutoff = cfg$de$pval_cutoff,
      lfc_cutoff  = cfg$de$lfc_cutoff,
      top_n       = cfg$de$top_n_labels
    ) +
      theme(
        legend.position = "none",
        plot.title       = element_text(color = ct_color, face = "bold", size = 23),
        # x-axis: label + text only on bottom row, larger
        axis.title.x = if (is_bottom) element_text(size = 23) else element_blank(),
        axis.text.x  = if (is_bottom) element_text(size = 21) else element_blank(),
        axis.ticks.x = if (is_bottom) element_line()          else element_blank(),
        # y-axis: label + text only on left column, larger
        axis.title.y = if (is_left)   element_text(size = 23) else element_blank(),
        axis.text.y  = if (is_left)   element_text(size = 21) else element_blank(),
        axis.ticks.y = if (is_left)   element_line()          else element_blank()
      )
  })

  # Extract the shared legend from a version that has it
  legend_source <- plot_volcano(
    de_results[[cts_use[1]]]$res,
    title = cts_use[1],
    pval_cutoff = cfg$de$pval_cutoff,
    lfc_cutoff  = cfg$de$lfc_cutoff,
    top_n       = 0   # no labels needed for legend source
  )
  shared_legend <- cowplot::get_legend(legend_source)

  # Assemble grid without legend, then attach shared legend on the right.
  # rel_widths gives the legend ~12% of total width.
  grid_no_legend <- cowplot::plot_grid(
    plotlist = plist,
    ncol     = ncol,
    align    = "hv",
    axis     = "tblr"
  )

  cowplot::plot_grid(
    grid_no_legend,
    shared_legend,
    ncol       = 2,
    rel_widths = c(1, 0.12)
  )
}

plot_top_de_gene_boxplots <- function(res_df, pb_counts, celltype, cfg, top_n = 4) {
  keep <- res_df %>%
    filter(!is.na(pvalue), pvalue < cfg$de$pval_cutoff, abs(log2FoldChange) >= cfg$de$lfc_cutoff) %>%
    slice_head(n = top_n) %>%
    pull(Gene)

  if (length(keep) == 0) return(NULL)

  # compute per-sample library sizes from full pseudobulk matrix
  lib_sizes <- colSums(pb_counts$counts)

  cnt <- as.data.frame(pb_counts$counts[keep, , drop = FALSE]) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Count") %>%
    left_join(
      pb_counts$coldata %>% rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    left_join(
      tibble(Sample = names(lib_sizes), lib_size = lib_sizes),
      by = "Sample"
    ) %>%
    mutate(logCPM = log2((Count / lib_size) * 1e6 + 1))

  pain_col <- cfg$meta$pain_col
  # Ensure Low Pain is always plotted to the left of High Pain
  cnt[[pain_col]] <- factor(cnt[[pain_col]], levels = c("Low Pain", "High Pain"))

  plist <- lapply(keep, function(g) {
    ggplot(filter(cnt, Gene == g), aes(x = .data[[pain_col]], y = logCPM, fill = .data[[pain_col]])) +
      geom_boxplot(outlier.shape = NA, width = 0.55) +
      geom_jitter(width = 0.08, size = 2) +
      scale_fill_manual(values = c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")) +
      labs(title = bquote(italic(.(g)) ~ "(" * .(celltype) * ")"), x = NULL, y = "log2 CPM") +
      theme_classic(base_size = 21) +
      theme(plot.title = element_text(face = "bold", size = 21), legend.position = "none")
  })

  wrap_plots(plist, ncol = min(2, length(plist)))
}

plot_deg_violin_grid <- function(obj, cfg, de_results, n_genes = 16, ncol = 4) {
  # Collect top DEGs across all cell types, pick up to n_genes unique genes
  top_genes <- lapply(names(de_results), function(ct) {
    de_results[[ct]]$res %>%
      filter(!is.na(pvalue), pvalue < 0.05) %>%
      arrange(pvalue) %>%
      slice_head(n = ceiling(n_genes / length(de_results))) %>%
      pull(Gene)
  }) %>% unlist() %>% unique()

  top_genes <- head(top_genes, n_genes)
  top_genes <- intersect(top_genes, rownames(obj[[cfg$rna_assay]]))
  if (length(top_genes) == 0) return(NULL)

  DefaultAssay(obj) <- cfg$rna_assay
  pain_col <- cfg$meta$pain_col

  # Subset to cells with pain group assignment
  cells_use <- rownames(obj@meta.data)[!is.na(obj@meta.data[[pain_col]])]
  sub <- subset(obj, cells = cells_use)

  plist <- lapply(top_genes, function(g) {
    VlnPlot(sub, features = g, group.by = pain_col, pt.size = 0,
            cols = c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")) +
      labs(title = bquote(bold(italic(.(g))))) +
      theme_classic(base_size = 20) +
      theme(
        plot.title = element_text(size = 21),
        legend.position = "none",
        axis.title.x = element_blank()
      )
  })

  wrap_plots(plist, ncol = ncol) +
    plot_annotation(title = "Top DEGs: High/Low Pain (single-cell expression)")
}

# ── Gene selection for DE heatmap: shared high-|LFC| genes ──────────────────
# Filters each cell type's DESeq2 results to pvalue < 0.05, then ranks genes
# by: (1) number of cell types they appear in (most shared first), and within
# ties (2) mean |log2FoldChange| across the cell types they appear in.
# Returns the top `top_n` genes by that combined score.

select_heatmap_genes <- function(res_list, top_n = 40,
                                  pval_cutoff = 0.05, lfc_cutoff = 0.585) {
  # Per-cell-type quota: take top floor(top_n / n_ct) genes per cell type ranked
  # by |LFC|, both pval AND lfc cutoffs enforced. This prevents pan-PBMC shared
  # genes from dominating — each cell type contributes its most distinctive hits.
  n_ct     <- length(res_list)
  quota    <- max(2L, floor(top_n / n_ct))   # min 2 per cell type

  per_ct_genes <- lapply(names(res_list), function(ct) {
    x <- res_list[[ct]]$res
    if (!all(c("Gene", "log2FoldChange", "pvalue") %in% colnames(x))) return(NULL)
    x %>%
      filter(!is.na(pvalue),
             pvalue          <  pval_cutoff,
             abs(log2FoldChange) >= lfc_cutoff) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      slice_head(n = quota) %>%
      pull(Gene)
  })

  # Pool unique genes; if fewer than top_n, backfill with next-best |LFC| genes
  # (still both-cutoff passing) that weren't already selected
  selected <- unique(unlist(per_ct_genes))

  if (length(selected) < top_n) {
    all_sig <- lapply(names(res_list), function(ct) {
      x <- res_list[[ct]]$res
      x %>%
        filter(!is.na(pvalue),
               pvalue              <  pval_cutoff,
               abs(log2FoldChange) >= lfc_cutoff) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        pull(Gene)
    }) %>% unlist() %>% unique()
    backfill <- setdiff(all_sig, selected)
    selected <- c(selected, head(backfill, top_n - length(selected)))
  }

  head(selected, top_n)
}

plot_de_heatmap <- function(res_list, cfg = NULL, top_n = 40,
                            ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC"),
                            pval_cutoff = 0.05, lfc_cutoff = 0.585,
                            genes_use = NULL) {
  ct_order_use <- intersect(ct_order, names(res_list))
  if (length(ct_order_use) == 0) return(NULL)
  
  # ── 1. Gene selection: shared genes appearing in ≥ 2 cell types (UpSet shared) ──
  if (is.null(genes_use) || length(genes_use) == 0) {
    sig_genes_per_ct <- lapply(names(res_list), function(ct) {
      x <- res_list[[ct]]$res
      if (!all(c("Gene", "log2FoldChange", "pvalue") %in% colnames(x))) return(NULL)
      x %>%
        filter(!is.na(pvalue), pvalue < pval_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
        pull(Gene) %>% unique()
    })
    names(sig_genes_per_ct) <- names(res_list)
    sig_genes_per_ct <- Filter(Negate(is.null), sig_genes_per_ct)
    
    all_genes_flat <- unlist(sig_genes_per_ct)
    gene_counts    <- table(all_genes_flat)
    genes_use      <- names(gene_counts[gene_counts >= 2])
    
    if (length(genes_use) == 0) {
      genes_use <- select_heatmap_genes(res_list, top_n = top_n,
                                        pval_cutoff = pval_cutoff, lfc_cutoff = lfc_cutoff)
    }
  }
  if (length(genes_use) == 0) return(NULL)
  
  # ── 2. Gene direction: compute mean LFC across significant cell types ─────────
  gene_mean_lfc <- sapply(genes_use, function(g) {
    lfcs <- sapply(names(res_list), function(ct) {
      x <- res_list[[ct]]$res
      r <- x[!is.na(x$Gene) & x$Gene == g &
               !is.na(x$pvalue) & x$pvalue < pval_cutoff, , drop = FALSE]
      if (nrow(r) == 0) NA_real_ else mean(r$log2FoldChange, na.rm = TRUE)
    })
    mean(lfcs, na.rm = TRUE)
  })
  
  lp_genes <- names(gene_mean_lfc[!is.na(gene_mean_lfc) & gene_mean_lfc < 0])
  hp_genes <- names(gene_mean_lfc[!is.na(gene_mean_lfc) & gene_mean_lfc >= 0])
  lp_genes <- lp_genes[order(-abs(gene_mean_lfc[lp_genes]))]
  hp_genes <- hp_genes[order(-abs(gene_mean_lfc[hp_genes]))]
  gene_col_order <- c(lp_genes, hp_genes)
  
  # ── 3. Build mean log2 CPM: rows = adjacent LP/HP per CT, cols = genes ────────
  row_order_keys <- unlist(lapply(ct_order_use, function(ct) {
    c(paste0(ct, "|||Low Pain"), paste0(ct, "|||High Pain"))
  }))
  
  col_data_list <- list()
  for (ct in ct_order_use) {
    pb      <- res_list[[ct]]$pb
    cnt     <- pb$counts
    coldata <- pb$coldata
    avail   <- intersect(genes_use, rownames(cnt))
    if (length(avail) == 0) next
    lib_sizes <- colSums(cnt)
    log2cpm   <- log2(sweep(cnt[avail, , drop = FALSE], 2, lib_sizes / 1e6, "/") + 1)
    for (grp in c("Low Pain", "High Pain")) {
      samps <- rownames(coldata)[as.character(coldata[["PainGroup"]]) == grp]
      samps <- intersect(samps, colnames(log2cpm))
      if (length(samps) == 0) next
      col_data_list[[paste0(ct, "|||", grp)]] <- rowMeans(log2cpm[, samps, drop = FALSE])
    }
  }
  if (length(col_data_list) == 0) return(NULL)
  
  row_order   <- intersect(row_order_keys, names(col_data_list))
  avail_genes <- intersect(gene_col_order, unique(unlist(lapply(col_data_list, names))))
  if (length(avail_genes) == 0) return(NULL)
  
  expr_mat <- do.call(rbind, lapply(row_order, function(rk) {
    v <- col_data_list[[rk]][avail_genes]; v[is.na(v)] <- 0; v
  }))
  rownames(expr_mat) <- row_order; colnames(expr_mat) <- avail_genes
  
  # ── 4. Per-gene min-max normalisation across conditions ───────────────────────
  expr_norm <- apply(expr_mat, 2, function(col) {
    col_c <- col - mean(col, na.rm = TRUE)
    mx    <- max(abs(col_c), na.rm = TRUE)
    if (mx == 0) return(rep(0, length(col_c)))
    col_c / mx
  })
  rownames(expr_norm) <- row_order; colnames(expr_norm) <- avail_genes
  
  n_rows     <- length(row_order)
  n_genes    <- length(avail_genes)
  n_lp_genes <- length(intersect(lp_genes, avail_genes))
  
  # ── 5. Clean row labels — "CellType\n(Low Pain)" / "(High Pain)" ─────────────
  clean_row_labels <- setNames(
    paste0(gsub("\\|\\|\\|.*$", "", row_order), "\n",
           ifelse(grepl("Low Pain", row_order), "(Low Pain)", "(High Pain)")),
    row_order
  )
  
  # ── 6. Long format — y axis reversed so first CT is at top ───────────────────
  df_long <- as.data.frame(expr_norm) %>%
    rownames_to_column("Row") %>%
    pivot_longer(-Row, names_to = "Gene", values_to = "ExprNorm") %>%
    mutate(
      Row  = factor(Row,  levels = rev(row_order)),
      Gene = factor(Gene, levels = avail_genes)
    )
  
  # ── 7. Main heatmap — no row labels (annotation strip carries identity) ─────
  p_main <- ggplot(df_long, aes(x = Gene, y = Row, fill = ExprNorm)) +
    geom_tile(color = "white", linewidth = 0.2) +
    { if (n_lp_genes > 0 && n_lp_genes < n_genes)
      geom_segment(aes(x = n_lp_genes + 0.5, xend = n_lp_genes + 0.5,
                       y = 0.5, yend = n_rows + 0.5),
                   color = "black", linewidth = 1.5, inherit.aes = FALSE)
    } +
    {
      ct_divs <- seq(2.5, n_rows - 0.5, by = 2)
      if (length(ct_divs) > 0)
        geom_segment(data = data.frame(y = ct_divs),
                     aes(x = 0.5, xend = n_genes + 0.5, y = y, yend = y),
                     color = "grey50", linewidth = 0.5, inherit.aes = FALSE)
    } +
    scale_y_discrete(labels = NULL, expand = expansion(add = 0.5)) +
    scale_x_discrete(position = "bottom") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Min-Max\nNormalized\nExpression") +
    annotate("text", x = max(1, n_lp_genes / 2) + 0.5, y = n_rows + 1.5,
             label = "Low Pain Enriched", fontface = "bold",
             size = 26 / 2.835, hjust = 0.5) +
    { if (n_genes - n_lp_genes > 0)
      annotate("text",
               x = n_lp_genes + (n_genes - n_lp_genes) / 2 + 0.5,
               y = n_rows + 1.5,
               label = "High Pain Enriched", fontface = "bold",
               size = 26 / 2.835, hjust = 0.5)
    } +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 22),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.line.y        = element_blank(),
      plot.margin        = margin(30, 8, 8, 20),
      legend.position    = "none",
      legend.key.height  = unit(5, "cm"),
      legend.key.width   = unit(0.8, "cm"),
      legend.title       = element_text(size = 22),
      legend.text        = element_text(size = 21),
      panel.grid         = element_blank()
    )
  
  # ── 8. Annotation strip drawn inside p_main via annotate() ──────────────────
  # This guarantees pixel-perfect row alignment because the rects share the
  # same y-axis coordinate system as the heatmap tiles. coord_cartesian(clip="off")
  # allows them to render in the left margin area.
  if (!is.null(cfg)) {
    dark_cols  <- cfg$merged_cell_cols
    light_cols <- colorspace::lighten(dark_cols, amount = 0.55)
    ann_colors <- c(
      setNames(light_cols, paste0(names(light_cols), "::Low Pain")),
      setNames(dark_cols,  paste0(names(dark_cols),  "::High Pain"))
    )
    ann_labels <- c(
      setNames(paste0(names(light_cols), "\nLow Pain"),  paste0(names(light_cols), "::Low Pain")),
      setNames(paste0(names(dark_cols),  "\nHigh Pain"), paste0(names(dark_cols),  "::High Pain"))
    )

    # factor(row_order, levels=rev(row_order)): row_order[1] → level n_rows (top, y=n_rows)
    y_positions     <- seq(n_rows, 1, by = -1)
    ann_rect_colors <- sapply(row_order, function(rk) {
      ct <- gsub("\\|\\|\\|.*$", "", rk)
      pg <- gsub("^.*\\|\\|\\|", "", rk)
      if (pg == "Low Pain") unname(light_cols[ct]) else unname(dark_cols[ct])
    })

    # Embed annotation rects in p_main at x = -1 (left of the heatmap data area)
    p_main <- p_main +
      annotate("rect",
               xmin  = -0.7, xmax  = 0.3,
               ymin  = y_positions - 0.45, ymax = y_positions + 0.45,
               fill  = ann_rect_colors, color = "white", linewidth = 0.3) +
      theme(plot.margin = margin(30, 8, 8, 35))   # extra left margin for the strip

    heatmap_body <- p_main

    # Standalone cell-type legend (separate ggplot, no scale conflict)
    ann_legend_df   <- data.frame(fill_key = factor(names(ann_colors), levels = names(ann_colors)),
                                  x = 1, y = 1)
    ann_legend_plot <- ggplot(ann_legend_df, aes(x = x, y = y, fill = fill_key)) +
      geom_tile() +
      scale_fill_manual(values = ann_colors, labels = ann_labels,
                        name = "Cell Type & Pain Group") +
      guides(fill = guide_legend(nrow = 2, ncol = 6, byrow = TRUE)) +
      theme_void() +
      theme(
        legend.position      = "bottom",
        legend.justification = "center",
        legend.text          = element_text(size = 21),
        legend.title         = element_text(size = 22, face = "bold"),
        legend.key.size      = unit(0.65, "cm")
      )
    heatmap_legend <- cowplot::get_legend(ann_legend_plot)

    # Expression gradient legend — thin horizontal bar with −1 / 0 / 1 labels
    expr_legend_plot <- p_main +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-1, 1),
                           breaks = c(-1, 0, 1), labels = c("-1", "0", "1"),
                           name = "Min-Max\nNormalized\nExpression",
                           guide = guide_colorbar(direction      = "horizontal",
                                                  title.position = "top",
                                                  barwidth       = unit(4, "cm"),
                                                  barheight      = unit(0.4, "cm"))) +
      theme(legend.position = "bottom", legend.justification = "center")
    expr_legend <- cowplot::get_legend(expr_legend_plot)

    # Stack: heatmap body on top, both legends side-by-side below
    legends_row <- cowplot::plot_grid(
      expr_legend,
      heatmap_legend,
      nrow = 1, rel_widths = c(1, 3)
    )
    p <- cowplot::plot_grid(heatmap_body, legends_row,
                            ncol = 1, rel_heights = c(1, 0.18))
  } else {
    p <- p_main
  }
  
  p
  
}

plot_gene_overlap_venn <- function(res_list, cfg = NULL, max_sets = 4, pval_cutoff = 0.05, lfc_cutoff = 0.585) {
  suppressPackageStartupMessages(library(ggVennDiagram))

  gene_sets <- lapply(res_list, function(x) {
    x$res %>%
      filter(!is.na(pvalue), pvalue < pval_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
      pull(Gene) %>%
      unique()
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) < 2) return(NULL)

  # ggVennDiagram supports max ~7 sets; keep top N by DEG count
  if (length(gene_sets) > max_sets) {
    top_idx <- order(lengths(gene_sets), decreasing = TRUE)[seq_len(max_sets)]
    gene_sets <- gene_sets[sort(top_idx)]
  }

  # Save original names for color lookup BEFORE shortening
  original_names <- names(gene_sets)

  # Shorten long cell type names for display
  short_names <- gsub("Interferon ", "IFN-", original_names)
  short_names <- gsub(" Monocyte", " Mono", short_names)
  names(gene_sets) <- short_names

  # Look up set edge colors from cell type palette using original names
  set_colors <- rep("black", length(gene_sets))
  if (!is.null(cfg) && !is.null(cfg$cell_cols)) {
    for (i in seq_along(original_names)) {
      if (original_names[i] %in% names(cfg$cell_cols)) {
        set_colors[i] <- cfg$cell_cols[[original_names[i]]]
      }
    }
  }

  p <- ggVennDiagram(gene_sets, label_alpha = 0, set_size = 12,
                     label_size = 10, edge_size = 1.5,
                     set_color = set_colors) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste0("Overlap of DE Genes (top ", length(gene_sets), " cell types)")) +
    theme(
      plot.title = element_text(face = "bold", size = 23),
      plot.margin = margin(10, 30, 10, 30)
    ) + scale_x_continuous(expand = expansion(mult = .2))

  # Attach gene sets as attribute so caller can write CSV
  attr(p, "gene_sets") <- gene_sets
  p
}

run_simple_gsea <- function(res_df, pathways_list = NULL) {
  suppressPackageStartupMessages(library(fgsea))
  if (is.null(pathways_list)) {
    return(NULL)
  }
  stats <- res_df$stat
  names(stats) <- res_df$Gene
  stats <- sort(stats, decreasing = TRUE)
  fgsea(pathways = pathways_list, stats = stats)
}

plot_gsea_bar <- function(gsea_tbl, top_n = 12) {
  if (is.null(gsea_tbl) || nrow(gsea_tbl) == 0) return(NULL)
  df <- gsea_tbl %>%
    arrange(pvalue) %>%
    slice_head(n = top_n) %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES))

  ggplot(df, aes(x = pathway, y = NES)) +
    geom_col() +
    coord_flip() +
    labs(title = "GSEA Functional Analysis", x = NULL, y = "NES") +
    theme_classic(base_size = 21) +
    theme(plot.title = element_text(face = "bold", size = 23))
}

# ── Highlight gene boxplots across cell types split by pain group ────────────

plot_highlight_gene_boxplots <- function(obj, cfg, genes = NULL,
                                          de_results = NULL,
                                          pval_cutoff = NULL) {
  if (is.null(genes)) genes <- cfg$highlight_genes
  genes <- intersect(genes, rownames(obj[[cfg$rna_assay]]))
  if (length(genes) == 0) {
    warning("No highlight genes found in object")
    return(NULL)
  }

  DefaultAssay(obj) <- cfg$rna_assay
  sample_col <- cfg$meta$sample_col
  pain_col   <- cfg$meta$pain_col

  # Use MergedCellType + merged_cell_cols (6 cell types)
  if (!"MergedCellType" %in% colnames(obj@meta.data)) {
    obj <- add_merged_celltype(obj, cfg)
  }
  merged_cols  <- cfg$merged_cell_cols
  ct_order_use <- names(merged_cols)

  cells_use <- rownames(obj@meta.data)[
    !is.na(obj@meta.data[[pain_col]]) &
    !is.na(obj@meta.data[["MergedCellType"]])
  ]
  sub <- subset(obj, cells = cells_use)

  counts <- GetAssayData(sub, assay = cfg$rna_assay, layer = "counts")
  md <- as.data.frame(sub@meta.data[, c("MergedCellType", sample_col, pain_col), drop = FALSE])
  md[] <- lapply(md, as.character)
  colnames(md) <- c("CellType", "Sample", "PainGroup")

  # Pseudobulk per merged cell type
  pb_long <- lapply(ct_order_use, function(ct) {
    idx <- which(md$CellType == ct)
    if (length(idx) < 10) return(NULL)
    sample_ids <- md$Sample[idx]
    ct_counts  <- counts[genes, idx, drop = FALSE]
    pb <- sapply(split(seq_along(sample_ids), sample_ids), function(i) {
      Matrix::rowSums(ct_counts[, i, drop = FALSE])
    })
    if (is.null(dim(pb))) return(NULL)
    lib_sizes      <- colSums(counts[, idx, drop = FALSE])
    lib_per_sample <- sapply(split(seq_along(sample_ids), sample_ids),
                             function(i) sum(lib_sizes[i]))
    lib_per_sample <- lib_per_sample[colnames(pb)]
    cpm     <- sweep(pb, 2, lib_per_sample, "/") * 1e6
    log2cpm <- log2(cpm + 1)
    as.data.frame(as.table(log2cpm)) %>%
      setNames(c("Gene", "Sample", "log2CPM")) %>%
      mutate(CellType = ct)
  }) %>% bind_rows()

  if (nrow(pb_long) == 0) return(NULL)

  sample_pain <- md %>% distinct(Sample, PainGroup)
  pb_long <- left_join(pb_long, sample_pain, by = "Sample")
  pb_long$CellType  <- factor(pb_long$CellType,  levels = ct_order_use)
  pb_long$PainGroup <- factor(pb_long$PainGroup, levels = c("Low Pain", "High Pain"))

  # Color map: all Low Pain (light) first, then all High Pain (dark)
  # This order + guide_legend(nrow=6, ncol=2) gives a 2-col legend:
  # left col = Low Pain, right col = High Pain
  dark_cols  <- merged_cols
  light_cols <- colorspace::lighten(merged_cols, amount = 0.55)
  fill_map <- c(
    setNames(light_cols, paste0(names(light_cols), "::Low Pain")),
    setNames(dark_cols,  paste0(names(dark_cols),  "::High Pain"))
  )
  pb_long$fill_key <- factor(
    paste0(pb_long$CellType, "::", pb_long$PainGroup),
    levels = names(fill_map)
  )

  # Legend labels: just the cell type name (Pain Group shown by shade)
  legend_labels <- c(
    setNames(paste0(names(light_cols), "\nLow Pain"),  paste0(names(light_cols), "::Low Pain")),
    setNames(paste0(names(dark_cols),  "\nHigh Pain"), paste0(names(dark_cols),  "::High Pain"))
  )

  # Significance: look up DESeq2 nominal p-value for each gene × cell type.
  # Falls back to no stars if de_results not provided.
  if (is.null(pval_cutoff)) pval_cutoff <- cfg$de$pval_cutoff

  sig_df <- if (!is.null(de_results)) {
    # Build one row per gene × merged cell type from DESeq2 output
    lapply(names(de_results), function(ct) {
      res <- de_results[[ct]]$res
      if (!all(c("Gene", "pvalue", "log2FoldChange") %in% colnames(res))) return(NULL)
      res %>%
        filter(Gene %in% genes) %>%
        transmute(
          Gene,
          CellType  = ct,
          deseq_p   = pvalue,
          deseq_lfc = log2FoldChange
        )
    }) %>% bind_rows() %>%
      mutate(
        star = if_else(!is.na(deseq_p) & deseq_p < pval_cutoff, "*", ""),
        # blue = gene up in Low Pain (negative LFC); red = up in High Pain
        star_color = case_when(
          is.na(deseq_p) | deseq_p >= pval_cutoff ~ "transparent",
          deseq_lfc <= 0                           ~ "#2166AC",   # Low Pain higher
          TRUE                                     ~ "red4"        # High Pain higher
        )
      )
  } else {
    # No de_results: return empty frame so downstream code works unchanged
    data.frame(Gene = character(0), CellType = character(0),
               star = character(0), star_color = character(0))
  }

  # Build one panel per gene
  plist <- lapply(genes, function(g) {
    df <- pb_long[pb_long$Gene == g, ]
    if (nrow(df) == 0) return(NULL)
    tt <- sig_df[sig_df$Gene == g, ]

    p <- ggplot(df, aes(x = CellType, y = log2CPM, fill = fill_key)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                   outlier.shape = NA) +
      scale_fill_manual(values = fill_map, labels = legend_labels,
                        name = "Cell Type & Pain Group") +
      guides(fill = guide_legend(nrow = 2, ncol = 6, byrow = TRUE)) +
      labs(title = bquote(italic(.(g))), x = NULL, y = "log2 CPM") +
      theme_classic(base_size = 21) +
      theme(
        plot.title   = element_text(face = "bold", size = 23),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
        legend.text  = element_text(size = 19),
        legend.title = element_text(size = 20, face = "bold")
      )

    # Significance stars — one geom_text call per unique color to allow mapping
    if (nrow(tt) > 0) {
      sig_tt <- tt[tt$star != "", ]
      if (nrow(sig_tt) > 0) {
        sig_tt$y <- max(df$log2CPM, na.rm = TRUE) * 1.05
        for (col in unique(sig_tt$star_color)) {
          sub <- sig_tt[sig_tt$star_color == col, ]
          p <- p + geom_text(data = sub, aes(x = CellType, y = y, label = star),
                             inherit.aes = FALSE, size = 9, color = col)
        }
      }
    }
    p
  })
  plist <- Filter(Negate(is.null), plist)
  if (length(plist) == 0) return(NULL)

  # Pagination: 12 per page (4 cols × 3 rows)
  per_page  <- 12
  ncol_page <- 4

  chunks <- split(plist, ceiling(seq_along(plist) / per_page))

  # Build a standalone legend ggplot from fill_map — robust alternative to get_legend
  legend_df <- data.frame(
    fill_key = factor(names(fill_map), levels = names(fill_map)),
    x = 1, y = 1
  )
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = fill_key)) +
    geom_tile() +
    scale_fill_manual(values = fill_map, labels = legend_labels,
                      name = "Cell Type & Pain Group") +
    guides(fill = guide_legend(nrow = 2, ncol = 6, byrow = TRUE)) +
    theme_void() +
    theme(
      legend.position      = "bottom",
      legend.justification = "center",
      legend.text          = element_text(size = 20),
      legend.title         = element_text(size = 21, face = "bold"),
      legend.key.size      = unit(0.65, "cm")
    )
  shared_legend <- cowplot::get_legend(legend_plot)

  pages <- lapply(seq_along(chunks), function(i) {
    chunk      <- chunks[[i]]
    n_in_chunk <- length(chunk)

    # Bottom row: last ncol_page positions
    bottom_start <- n_in_chunk - ((n_in_chunk - 1) %% ncol_page)
    bottom_idx   <- seq(bottom_start, n_in_chunk)
    # Left column: positions 1, 5, 9 …
    left_idx     <- seq(1, n_in_chunk, by = ncol_page)

    chunk <- lapply(seq_along(chunk), function(j) {
      p <- chunk[[j]]
      is_bottom <- j %in% bottom_idx
      is_left   <- j %in% left_idx

      p + theme(
        # x-axis labels only on bottom row
        axis.text.x  = if (is_bottom) element_text(angle = 45, hjust = 1, size = 20)
                       else element_blank(),
        axis.ticks.x = if (is_bottom) element_line() else element_blank(),
        # y-axis label only on left column
        axis.title.y = if (is_left) element_text(size = 21) else element_blank(),
        # no legend on any panel — handled separately below
        legend.position = "none"
      )
    })

    grid <- wrap_plots(chunk, ncol = ncol_page) +
      plot_annotation(title = paste0(
        "Highlight Genes: Merged Cell Types, Pseudobulk log2 CPM (",
        i, "/", length(chunks), ")"
      ))

    # Attach the shared legend below the grid, giving it ~12% of total height
    cowplot::plot_grid(
      grid,
      shared_legend,
      ncol        = 1,
      rel_heights = c(1, 0.15)
    )
  })

  pages
}

# ── DotPlot of highlight genes across clusters ───────────────────────────────

plot_highlight_gene_dotplot <- function(obj, cfg, genes = NULL) {
  if (is.null(genes)) genes <- cfg$highlight_genes
  genes <- intersect(genes, rownames(obj[[cfg$rna_assay]]))
  if (length(genes) == 0) return(NULL)

  DefaultAssay(obj) <- cfg$rna_assay
  ct_col <- cfg$meta$celltype_col
  old_id <- Idents(obj)
  Idents(obj) <- obj@meta.data[[ct_col]]

  p <- DotPlot(obj, features = genes) +
    coord_flip() +
    labs(title = "Highlight Genes Across Cell Types", x = NULL, y = NULL) +
    theme_classic(base_size = 19) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
      axis.text.y = element_text(size = 17)
    )

  Idents(obj) <- old_id
  p
}

# ── Highlight gene boxplots using 6 merged cell types with t-test ────────────

plot_highlight_gene_merged <- function(obj, cfg, genes = NULL, ncol = 3, per_page = 9,
                                       ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  if (is.null(genes)) genes <- cfg$highlight_genes
  genes <- intersect(genes, rownames(obj[[cfg$rna_assay]]))
  if (length(genes) == 0) return(NULL)

  DefaultAssay(obj) <- cfg$rna_assay
  sample_col <- cfg$meta$sample_col
  pain_col   <- cfg$meta$pain_col

  cells_use <- rownames(obj@meta.data)[!is.na(obj@meta.data[[pain_col]]) &
                                         !is.na(obj@meta.data[["MergedCellType"]])]
  sub <- subset(obj, cells = cells_use)
  counts <- GetAssayData(sub, assay = cfg$rna_assay, layer = "counts")
  md <- as.data.frame(sub@meta.data[, c(sample_col, pain_col, "MergedCellType"), drop = FALSE])
  md[] <- lapply(md, as.character)

  all_cts <- unique(md$MergedCellType)
  merged_levels <- intersect(ct_order, all_cts)   # canonical order, present only

  # Build pseudobulk log2 CPM per sample x merged cell type
  pb_long <- lapply(merged_levels, function(ct) {
    idx <- which(md$MergedCellType == ct)
    if (length(idx) < 10) return(NULL)
    sample_ids <- md[[sample_col]][idx]
    ct_counts <- counts[genes, idx, drop = FALSE]
    pb <- sapply(split(seq_along(sample_ids), sample_ids), function(i) {
      Matrix::rowSums(ct_counts[, i, drop = FALSE])
    })
    if (is.null(dim(pb))) return(NULL)
    lib_per_sample <- sapply(split(seq_along(sample_ids), sample_ids), function(i) {
      sum(Matrix::colSums(counts[, idx[i], drop = FALSE]))
    })
    lib_per_sample <- lib_per_sample[colnames(pb)]
    cpm <- sweep(pb, 2, lib_per_sample, "/") * 1e6
    log2cpm <- log2(cpm + 1)
    as.data.frame(as.table(log2cpm)) %>%
      setNames(c("Gene", "Sample", "log2CPM")) %>%
      mutate(CellType = ct)
  }) %>% bind_rows()

  if (nrow(pb_long) == 0) return(NULL)

  sample_pain <- md %>% distinct(.data[[sample_col]], .data[[pain_col]]) %>%
    setNames(c("Sample", "PainGroup"))
  pb_long <- left_join(pb_long, sample_pain, by = "Sample")
  pb_long$CellType  <- factor(pb_long$CellType,  levels = intersect(ct_order, unique(pb_long$CellType)))
  # Low Pain left of High Pain in every dodged boxplot
  pb_long$PainGroup <- factor(pb_long$PainGroup, levels = c("Low Pain", "High Pain"))

  # T-test per gene x cell type
  ttest_df <- pb_long %>%
    group_by(Gene, CellType) %>%
    summarize(
      pval = tryCatch(t.test(log2CPM[PainGroup == "Low Pain"],
                             log2CPM[PainGroup == "High Pain"])$p.value,
                      error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(star = if_else(pval < 0.05, "*", ""))  # single star for any significant difference

  plist <- lapply(genes, function(g) {
    df <- pb_long[pb_long$Gene == g, ]
    if (nrow(df) == 0) return(NULL)
    tt <- ttest_df[ttest_df$Gene == g, ]

    p <- ggplot(df, aes(x = CellType, y = log2CPM, fill = PainGroup)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
      geom_jitter(aes(group = PainGroup),
                  position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
                  size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F"),
                        name = "Pain Group") +
      labs(title = bquote(italic(.(g))), x = NULL, y = "log2 CPM") +
      theme_classic(base_size = 19) +
      theme(
        plot.title = element_text(face = "bold.italic", size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        legend.position = "none"
      )

    # Add significance stars
    if (nrow(tt) > 0) {
      sig_tt <- tt[tt$star != "", ]
      if (nrow(sig_tt) > 0) {
        sig_tt$y <- max(df$log2CPM, na.rm = TRUE) * 1.05
        p <- p + geom_text(data = sig_tt, aes(x = CellType, y = y, label = star),
                           inherit.aes = FALSE, size = 10)
      }
    }
    p
  })
  plist <- Filter(Negate(is.null), plist)
  if (length(plist) == 0) return(NULL)

  chunks <- split(plist, ceiling(seq_along(plist) / per_page))
  lapply(seq_along(chunks), function(i) {
    wrap_plots(chunks[[i]], ncol = ncol) +
      plot_annotation(title = paste0("Highlight Genes: 6 Merged Cell Types (",
                                     i, "/", length(chunks), ")"))
  })
}

# ── UpSet plot of DE gene overlap across cell types ──────────────────────────

build_upset_data <- function(de_results, pval_cutoff = 0.05, lfc_cutoff = 1) {
  gene_sets <- lapply(de_results, function(x) {
    x$res %>%
      filter(!is.na(pvalue), pvalue < pval_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
      pull(Gene) %>%
      unique()
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) < 2) return(NULL)

  all_genes <- unique(unlist(gene_sets))
  mat <- sapply(gene_sets, function(gs) as.integer(all_genes %in% gs))
  rownames(mat) <- all_genes
  list(data = as.data.frame(mat), sets = names(gene_sets), gene_sets = gene_sets)
}

save_upset_plot <- function(upset_data, out_png,
                            width = 12, height = 6, dpi = 300) {
  if (is.null(upset_data)) return(NULL)
  suppressPackageStartupMessages(library(UpSetR))

  png(out_png, width = width, height = height, units = "in",
      res = dpi, type = "cairo")
  print(upset(
    upset_data$data,
    sets = rev(upset_data$sets),
    keep.order = TRUE,
    order.by = "freq",
    text.scale = c(2.0, 1.8, 1.6, 1.4, 2.0, 1.6),
    mb.ratio = c(0.6, 0.4),
    mainbar.y.label = "Shared DE Genes",
    sets.x.label = "DE Genes per Cell Type"
  ))
  invisible(dev.off())

  out_png
}

# ── Bubble plot of DE gene counts per cell type ──────────────────────────────

plot_de_bubble_summary <- function(de_results, pval_cutoff = 0.05, lfc_cutoff = 1,
                                   ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  df <- lapply(names(de_results), function(ct) {
    res <- de_results[[ct]]$res
    up   <- sum(!is.na(res$pvalue) & res$pvalue < pval_cutoff & res$log2FoldChange >= lfc_cutoff, na.rm = TRUE)
    down <- sum(!is.na(res$pvalue) & res$pvalue < pval_cutoff & res$log2FoldChange <= -lfc_cutoff, na.rm = TRUE)
    tibble(CellType = ct,
           Direction = c("Low Pain", "High Pain"),
           Count = c(down, up))
  }) %>% bind_rows()

  if (nrow(df) == 0 || all(df$Count == 0)) return(NULL)

  ct_levels_use <- intersect(ct_order, unique(df$CellType))
  df$CellType  <- factor(df$CellType,  levels = ct_levels_use)
  # Low Pain on bottom row, High Pain on top — mirrors left→right convention on y-axis
  df$Direction <- factor(df$Direction, levels = c("Low Pain", "High Pain"))

  ggplot(df, aes(x = CellType, y = Direction, size = Count, color = Direction)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = Count), size = 3, vjust = -1.2, color = "black") +
    scale_size_area(max_size = 18) +
    scale_color_manual(values = c("High Pain" = "#D95F5F", "Low Pain" = "#6BAED6")) +
    labs(title = "Number of Differential Genes per Cell Type",
         x = NULL, y = NULL) +
    theme_classic(base_size = 21) +
    theme(
      plot.title = element_text(face = "bold", size = 23),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 19),
      legend.position = "right"
    ) +
    guides(size = guide_legend(title = "# DE Genes"))
}

# ── GSEA per cell type using msigdbr ─────────────────────────────────────────

run_gsea_all_celltypes <- function(de_results, species = "Homo sapiens",
                                    category = "H", subcategory = NULL,
                                    min_size = 15, max_size = 500) {
  suppressPackageStartupMessages({
    library(fgsea)
    library(msigdbr)
  })

  pw_df <- msigdbr(species = species, category = category, subcategory = subcategory)
  pathways <- split(pw_df$gene_symbol, pw_df$gs_name)

  gsea_results <- list()
  for (ct in names(de_results)) {
    res <- de_results[[ct]]$res
    if (!"stat" %in% colnames(res)) next
    stats <- res$stat
    names(stats) <- res$Gene
    stats <- stats[!is.na(stats)]
    stats <- sort(stats, decreasing = TRUE)

    gsea_out <- fgsea(pathways = pathways, stats = stats,
                       minSize = min_size, maxSize = max_size)
    gsea_out$CellType <- ct
    gsea_results[[ct]] <- gsea_out
  }

  bind_rows(gsea_results)
}

# ── Leading edge extractor for pathways of interest ─────────────────────────

get_leading_edge <- function(gsea_all, pathway_patterns, pval_col = "pval",
                              sig_cutoff = 0.05) {
  pattern <- paste(pathway_patterns, collapse = "|")

  gsea_all %>%
    filter(grepl(pattern, pathway, ignore.case = TRUE)) %>%
    filter(!is.na(.data[[pval_col]]), .data[[pval_col]] < sig_cutoff) %>%
    select(pathway, CellType, NES, !!sym(pval_col), leadingEdge) %>%
    arrange(pathway, .data[[pval_col]]) %>%
    mutate(
      n_leading = lengths(leadingEdge),
      genes_str = sapply(leadingEdge, paste, collapse = ", ")
    )
}

plot_gsea_heatmap <- function(gsea_all, top_n = 20, pval_col = "pval", sig_cutoff = 0.05,
                              ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  if (is.null(gsea_all) || nrow(gsea_all) == 0) return(NULL)

  # Drop list-columns (e.g. leadingEdge) before any dplyr operations.
  # Coerce to plain data.frame — fgsea returns a data.table where [,logical]
  # row-selects instead of column-selects, returning a logical scalar.
  gsea_all <- as.data.frame(gsea_all)
  gsea_all <- gsea_all[, !vapply(gsea_all, is.list, logical(1)), drop = FALSE]

  # Filter to significant, rank by max |NES|
  sig <- gsea_all %>% filter(.data[[pval_col]] < sig_cutoff)
  if (nrow(sig) == 0) return(NULL)

  top_pw <- sig %>%
    group_by(pathway) %>%
    summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_NES)) %>%
    slice_head(n = top_n) %>%
    pull(pathway)

  ct_levels_use <- intersect(ct_order, unique(gsea_all$CellType))

  df <- gsea_all %>%
    filter(pathway %in% top_pw) %>%
    mutate(
      pathway  = factor(pathway,  levels = rev(top_pw)),
      CellType = factor(CellType, levels = ct_levels_use)
    )

  ggplot(df, aes(x = CellType, y = pathway, fill = NES)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(title = "GSEA: Hallmark Pathways (ranked by |NES|)",
         x = NULL, y = NULL, fill = "NES") +
    theme_minimal(base_size = 19) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 17),
      plot.title = element_text(face = "bold", size = 22)
    )
}

plot_gsea_dotplot <- function(gsea_all, top_n_per_ct = 4, min_per_ct = 2,
                              pval_col = "pval", sig_cutoff = 0.05,
                              ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC"),
                              force_pathways = NULL) {
  if (is.null(gsea_all) || nrow(gsea_all) == 0) return(NULL)

  # Drop list-columns (e.g. leadingEdge) — dplyr count/slice choke on them.
  # Coerce to plain data.frame — fgsea returns a data.table where [,logical]
  # row-selects instead of column-selects, returning a logical scalar.
  gsea_flat <- as.data.frame(gsea_all)
  gsea_flat <- gsea_flat[, !vapply(gsea_flat, is.list, logical(1)), drop = FALSE]

  sig <- gsea_flat %>% filter(.data[[pval_col]] < sig_cutoff)
  if (nrow(sig) == 0) return(NULL)

  ct_levels_use <- intersect(ct_order, unique(gsea_flat$CellType))

  # ── Pathway selection ────────────────────────────────────────────────────────
  # If force_pathways is supplied it is EXCLUSIVE: only those pathways are shown,
  # in exactly the order given. The per-cell-type quota is skipped entirely.
  if (!is.null(force_pathways) && length(force_pathways) > 0) {
    valid_forced   <- intersect(force_pathways, gsea_flat$pathway)
    invalid_forced <- setdiff(force_pathways, gsea_flat$pathway)
    if (length(invalid_forced) > 0) {
      message("  GSEA dotplot: force_pathways not found in results and will be skipped: ",
              paste(invalid_forced, collapse = ", "))
    }
    if (length(valid_forced) == 0) return(NULL)
    # Preserve the user-specified order exactly
    pw_order <- valid_forced
    top_pw   <- valid_forced
  } else {
    # Default: per-cell-type quota logic
    sig_split <- split(sig, sig$CellType)
    top_pw_per_ct <- do.call(rbind, lapply(sig_split, function(x) {
      x <- x[order(-abs(x$NES)), , drop = FALSE]
      head(x, top_n_per_ct)
    }))
    rownames(top_pw_per_ct) <- NULL

    # Log sparse cell types
    ct_tally <- as.data.frame(table(CellType = top_pw_per_ct$CellType),
                               stringsAsFactors = FALSE)
    for (ct in ct_levels_use) {
      n_sig <- ct_tally$Freq[ct_tally$CellType == ct]
      n_sig <- if (length(n_sig) == 0) 0L else n_sig
      if (n_sig < min_per_ct) {
        message("  GSEA dotplot: '", ct, "' has only ", n_sig,
                " significant pathway(s) at p < ", sig_cutoff,
                " (minimum requested: ", min_per_ct, ")")
      }
    }

    top_pw <- unique(top_pw_per_ct$pathway)

    # Order by global max |NES|
    pw_order <- gsea_flat %>%
      filter(pathway %in% top_pw) %>%
      group_by(pathway) %>%
      summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_abs_NES)) %>%
      pull(pathway)

    valid_forced   <- character(0)
    invalid_forced <- character(0)
  }

  # Clean pathway labels: replace underscores with spaces and strip HALLMARK_ prefix
  clean_label <- function(x) gsub("_", " ", gsub("^HALLMARK_", "", x))

  df <- sig %>%
    filter(pathway %in% top_pw) %>%
    mutate(
      neglog10p    = -log10(.data[[pval_col]]),
      pathway_label = clean_label(pathway),
      pathway_label = factor(pathway_label, levels = clean_label(rev(pw_order))),
      CellType      = factor(CellType, levels = ct_levels_use)
    )

  if (!is.null(force_pathways) && length(force_pathways) > 0) {
    subtitle_text <- "Forced pathway list (user-specified order); absent dot = not significant"
  } else {
    subtitle_text <- paste0("Top ", top_n_per_ct,
                            " pathways per cell type by |NES|; absent dot = not significant")
  }

  ggplot(df, aes(x = CellType, y = pathway_label, size = neglog10p, color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    scale_size_continuous(range = c(2, 8)) +
    scale_y_discrete(expand = expansion(add = 0.4)) +
    labs(
      title    = "GSEA: Hallmark Pathways",
      subtitle = subtitle_text,
      x = NULL, y = NULL, size = "-log10(p)", color = "NES"
    ) +
    theme_minimal(base_size = 23) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 23),
      axis.text.y        = element_text(size = 23),
      plot.title         = element_text(face = "bold", size = 23),
      plot.subtitle      = element_text(size = 1, color = "white"),
      legend.text        = element_text(size = 22),
      legend.title       = element_text(size = 23),
      panel.grid.major   = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0, "pt")
    )
}

# ── Balanced GSEA heatmap: per-cell-type quota ───────────────────────────────

plot_gsea_heatmap_balanced <- function(gsea_all, top_n_per_ct = 4,
                                        pval_col = "pval", sig_cutoff = 0.05,
                                        ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  if (is.null(gsea_all) || nrow(gsea_all) == 0) return(NULL)

  # Drop list-columns (e.g. leadingEdge) before any dplyr operations.
  # Coerce to plain data.frame — fgsea returns a data.table where [,logical]
  # row-selects instead of column-selects, returning a logical scalar.
  gsea_all <- as.data.frame(gsea_all)
  gsea_all <- gsea_all[, !vapply(gsea_all, is.list, logical(1)), drop = FALSE]

  sig <- gsea_all %>% filter(.data[[pval_col]] < sig_cutoff)
  if (nrow(sig) == 0) return(NULL)

  # Top N pathways per cell type by |NES|, then pool into a unique set
  top_pw <- sig %>%
    group_by(CellType) %>%
    slice_max(order_by = abs(NES), n = top_n_per_ct, with_ties = FALSE) %>%
    ungroup() %>%
    pull(pathway) %>%
    unique()

  # Order pooled pathways by global max |NES| for a sensible y-axis
  pw_order <- gsea_all %>%
    filter(pathway %in% top_pw) %>%
    group_by(pathway) %>%
    summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_NES)) %>%
    pull(pathway)

  ct_levels_use <- intersect(ct_order, unique(gsea_all$CellType))

  df <- gsea_all %>%
    filter(pathway %in% top_pw) %>%
    mutate(
      # Grey out tiles that are not significant
      NES_plot = ifelse(.data[[pval_col]] < sig_cutoff, NES, NA_real_),
      pathway  = factor(pathway,  levels = rev(pw_order)),
      CellType = factor(CellType, levels = ct_levels_use)
    )

  ggplot(df, aes(x = CellType, y = pathway, fill = NES_plot)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, na.value = "grey92",
                         name = "NES\n(sig. only)") +
    labs(
      title    = "GSEA: Hallmark Pathways (per-cell-type quota)",
      subtitle = paste0("Top ", top_n_per_ct,
                        " pathways per cell type by |NES|; grey = not significant"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 19) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y   = element_text(size = 17),
      plot.title    = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 17, color = "grey40")
    )
}

# ── Merged pseudobulk: collapse to 6 cell types ─────────────────────────────

add_merged_celltype <- function(obj, cfg) {
  ct_col <- cfg$meta$celltype_col
  ct <- as.character(obj@meta.data[[ct_col]])
  merged <- rep(NA_character_, length(ct))
  for (grp_name in names(cfg$merged_celltypes)) {
    merged[ct %in% cfg$merged_celltypes[[grp_name]]] <- grp_name
  }
  obj$MergedCellType <- merged
  obj
}

# ── Pseudobulk PCA ──────────────────────────────────────────────────────────

plot_pseudobulk_pca <- function(obj, cfg, ct_col_use, title_suffix = "",
                                 covariates = NULL, color_map = NULL,
                                 merged_ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  sample_col <- cfg$meta$sample_col
  pain_col   <- cfg$meta$pain_col

  cells_use <- rownames(obj@meta.data)[!is.na(obj@meta.data[[pain_col]]) &
                                         !is.na(obj@meta.data[[ct_col_use]])]
  sub <- subset(obj, cells = cells_use)
  md <- as.data.frame(sub@meta.data[, c(sample_col, pain_col, ct_col_use), drop = FALSE])
  md[] <- lapply(md, as.character)

  counts <- GetAssayData(sub, assay = cfg$rna_assay, layer = cfg$count_layer)

  # Create pseudobulk per sample x cell type
  combos <- unique(paste0(md[[sample_col]], "::", md[[ct_col_use]]))
  pb_list <- list()
  meta_rows <- list()

  for (combo in combos) {
    parts <- strsplit(combo, "::")[[1]]
    samp <- parts[1]; ct <- parts[2]
    idx <- which(md[[sample_col]] == samp & md[[ct_col_use]] == ct)
    if (length(idx) < 10) next
    pb_list[[combo]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    meta_rows[[combo]] <- data.frame(
      SampleCT = combo, Sample = samp, CellType = ct,
      PainGroup = md[[pain_col]][idx[1]],
      stringsAsFactors = FALSE
    )
  }

  if (length(pb_list) < 4) return(NULL)

  pb_mat <- do.call(cbind, pb_list)
  pb_meta <- bind_rows(meta_rows)
  rownames(pb_meta) <- pb_meta$SampleCT

  # Log2 CPM normalize
  lib_sizes <- colSums(pb_mat)
  cpm <- sweep(pb_mat, 2, lib_sizes, "/") * 1e6
  log2cpm <- log2(cpm + 1)

  # Filter low-expression genes
  keep <- rowMeans(log2cpm) > 1
  log2cpm <- log2cpm[keep, ]

  # PCA
  pca <- prcomp(t(log2cpm), scale. = TRUE)
  pca_df <- data.frame(pca$x[, 1:2])
  pca_df$Sample   <- pb_meta$Sample
  pca_df$CellType <- pb_meta$CellType
  pca_df$PainGroup <- pb_meta$PainGroup

  var_exp <- summary(pca)$importance[2, 1:2] * 100

  # Determine color palette
  if (is.null(color_map)) {
    if (ct_col_use == "MergedCellType" && !is.null(cfg$merged_cell_cols)) {
      color_map <- cfg$merged_cell_cols
    } else {
      color_map <- cfg$cell_cols
    }
  }

  # Factor CellType with matching order
  if (ct_col_use == "MergedCellType") {
    ct_order <- intersect(merged_ct_order, unique(pca_df$CellType))
  } else {
    ct_order <- intersect(names(color_map), unique(pca_df$CellType))
  }
  pca_df$CellType <- factor(pca_df$CellType, levels = ct_order)
  color_use <- color_map[ct_order]

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType, shape = PainGroup)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = color_use) +
    scale_shape_manual(values = c("Low Pain" = 16, "High Pain" = 17)) +
    labs(title = paste0("Pseudobulk PCA", title_suffix),
         x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         color = "Cell Type", shape = "Pain Group") +
    theme_classic(base_size = 23) +
    theme(
      plot.title   = element_text(face = "bold", size = 23),
      axis.title   = element_text(size = 23),
      axis.text    = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text  = element_text(size = 21)
    )

  list(plot = p, pca_df = pca_df, log2cpm = log2cpm, pb_meta = pb_meta)
}

# ── T-test cell type proportions ─────────────────────────────────────────────

run_ttest_proportions <- function(obj, cfg) {
  df <- celltype_sample_count_table(obj, cfg)
  ct_levels <- names(cfg$cell_cols)

  results <- lapply(ct_levels, function(ct) {
    sub <- df[df$CellType == ct, ]
    low  <- sub$prop[sub$PainGroup == "Low Pain"]
    high <- sub$prop[sub$PainGroup == "High Pain"]
    if (length(low) < 2 || length(high) < 2) {
      return(tibble(CellType = ct, P.Value = NA_real_,
                    LowMean = mean(low), HighMean = mean(high), method = "t.test"))
    }
    tt <- t.test(low, high)
    tibble(CellType = ct, P.Value = tt$p.value,
           LowMean = mean(low), HighMean = mean(high), method = "t.test")
  }) %>% bind_rows()

  results
}

# ── MASC (Mixed-effects Association of Single Cells) ─────────────────────────

run_masc <- function(obj, cfg) {
  suppressPackageStartupMessages(library(lme4))

  ct_col     <- cfg$meta$celltype_col
  sample_col <- cfg$meta$sample_col
  pain_col   <- cfg$meta$pain_col

  md <- as.data.frame(obj@meta.data[, c(ct_col, sample_col, pain_col), drop = FALSE])
  md[] <- lapply(md, as.character)
  colnames(md) <- c("CellType", "Sample", "PainGroup")
  md <- md[!is.na(md$PainGroup), ]
  md$PainGroup <- factor(md$PainGroup, levels = c("Low Pain", "High Pain"))

  # Sanitize cell type names for R formula compatibility
  ct_original <- md$CellType
  ct_safe <- make.names(md$CellType)
  name_map <- setNames(unique(ct_original), make.names(unique(ct_original)))
  md$CellType_safe <- ct_safe

  cluster <- md$CellType_safe
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, md)

  ct_names <- colnames(designmat)
  results <- list()

  for (i in seq_along(ct_names)) {
    tc <- ct_names[i]
    safe_label <- gsub("cluster", "", tc)
    orig_label <- name_map[[safe_label]]
    message("  MASC: ", orig_label)

    null_fm <- as.formula(paste0("`", tc, "` ~ 1 + (1|Sample)"))
    full_fm <- as.formula(paste0("`", tc, "` ~ PainGroup + (1|Sample)"))
    null_model <- tryCatch(
      glmer(null_fm, data = dataset, family = binomial, nAGQ = 1,
            control = glmerControl(optimizer = "bobyqa")),
      error = function(e) NULL)
    full_model <- tryCatch(
      glmer(full_fm, data = dataset, family = binomial, nAGQ = 1,
            control = glmerControl(optimizer = "bobyqa")),
      error = function(e) NULL)
    if (is.null(null_model) || is.null(full_model)) next
    lrt <- anova(null_model, full_model)
    or <- exp(fixef(full_model)[["PainGroupHigh Pain"]])
    results[[i]] <- tibble(
      CellType = orig_label,
      P.Value = lrt[["Pr(>Chisq)"]][2],
      OR = or,
      method = "MASC"
    )
  }

  bind_rows(results)
}

# ── Unified composition boxplot with legend ──────────────────────────────────

plot_composition_boxplot_with_legend <- function(obj, cfg, prop_res, title = "Cell Composition") {
  df <- celltype_sample_count_table(obj, cfg)
  df$PainGroup <- factor(df$PainGroup, levels = c("Low Pain", "High Pain"))
  df$CellType  <- factor(df$CellType, levels = names(cfg$cell_cols))

  dark_cols  <- cfg$cell_cols
  light_cols <- colorspace::lighten(cfg$cell_cols, amount = 0.55)

  # Build explicit fill for legend
  pain_cols <- c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")

  p <- ggplot(df, aes(x = CellType, y = prop, fill = PainGroup)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
    geom_jitter(aes(group = PainGroup),
                position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
                size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = pain_cols, name = "Pain Group") +
    scale_y_continuous(labels = percent_format()) +
    labs(title = title, x = NULL, y = "% of Total Cells") +
    theme_classic(base_size = 21) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 23),
      legend.position = "right"
    )

  # Add p-values if provided
  if (!is.null(prop_res) && "P.Value" %in% colnames(prop_res)) {
    ann <- prop_res %>%
      mutate(
        CellType = factor(CellType, levels = names(cfg$cell_cols)),
        label = case_when(
          P.Value < 0.001 ~ formatC(P.Value, format = "e", digits = 1),
          P.Value < 0.1   ~ paste0("p=", formatC(P.Value, format = "f", digits = 3)),
          TRUE             ~ paste0("p=", formatC(P.Value, format = "f", digits = 2))
        ),
        star = case_when(P.Value < 0.001 ~ "***", P.Value < 0.01 ~ "**",
                         P.Value < 0.05 ~ "*", P.Value < 0.1 ~ ".", TRUE ~ ""),
        label = paste0(label, star)
      ) %>%
      filter(!is.na(CellType)) %>%
      mutate(y = 0.94 * max(df$prop, na.rm = TRUE))

    if (nrow(ann) > 0) {
      p <- p + geom_text(data = ann, aes(x = CellType, y = y, label = label),
                         inherit.aes = FALSE, size = 3.5, angle = 30, hjust = 0)
    }
  }
  p
}

# ── DE gene count boxplot by pain group ──────────────────────────────────────

plot_de_count_boxplot <- function(de_results, pval_cutoff = 0.05, lfc_cutoff = 0.585,
                                  ct_order = c("CD4+ T", "CD8+ T", "NK", "B", "Monocytes", "DC")) {
  df <- lapply(names(de_results), function(ct) {
    res <- de_results[[ct]]$res
    up   <- sum(!is.na(res$pvalue) & res$pvalue < pval_cutoff & res$log2FoldChange >= lfc_cutoff, na.rm = TRUE)
    down <- sum(!is.na(res$pvalue) & res$pvalue < pval_cutoff & res$log2FoldChange <= -lfc_cutoff, na.rm = TRUE)
    tibble(CellType = ct,
           Direction = c("Low Pain", "High Pain"),
           Count = c(down, up))
  }) %>% bind_rows()

  if (nrow(df) == 0 || all(df$Count == 0)) return(NULL)

  ct_levels_use <- intersect(ct_order, unique(df$CellType))
  df$CellType  <- factor(df$CellType,  levels = ct_levels_use)
  # Low Pain bar dodges left of High Pain
  df$Direction <- factor(df$Direction, levels = c("Low Pain", "High Pain"))

  ggplot(df, aes(x = CellType, y = Count, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.8),
              vjust = -0.3, size = 4.5) +
    scale_fill_manual(values = c("High Pain" = "#D95F5F", "Low Pain" = "#6BAED6"),
                      name = "Pain Group") +
    labs(title = "Number of DE Genes per Cell Type", x = NULL, y = "# DE Genes") +
    theme_classic(base_size = 23) +
    theme(
      plot.title       = element_text(face = "bold", size = 23),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 23),
      axis.text.y      = element_text(size = 22),
      axis.title.y     = element_text(size = 23),
      legend.text      = element_text(size = 22),
      legend.title     = element_text(size = 22),
      legend.position  = "right"
    )
}

# ── Venn with merged cell type groups ────────────────────────────────────────

plot_merged_venn <- function(de_results, cfg, pval_cutoff = 0.05, lfc_cutoff = 0.585,
                             max_sets = 4,
                             force_sets = c("CD8+ T", "NK", "CD4+ T", "Monocytes")) {
  suppressPackageStartupMessages(library(ggVennDiagram))

  # de_results keys are already merged names (CD4+ T, Monocytes, etc.)
  gene_sets <- lapply(de_results, function(x) {
    x$res %>%
      filter(!is.na(pvalue), pvalue < pval_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
      pull(Gene) %>% unique()
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) < 2) return(NULL)

  # If force_sets provided, use those (in order); otherwise fall back to top N by DEG count
  if (!is.null(force_sets) && length(force_sets) > 0) {
    valid_forced <- intersect(force_sets, names(gene_sets))
    if (length(valid_forced) < 2) {
      message("  Merged Venn: fewer than 2 forced cell types have DEGs — falling back to top N")
    } else {
      gene_sets <- gene_sets[valid_forced]
    }
  } else if (length(gene_sets) > max_sets) {
    top_idx <- order(lengths(gene_sets), decreasing = TRUE)[seq_len(max_sets)]
    gene_sets <- gene_sets[sort(top_idx)]
  }

  # Look up edge colors from merged_cell_cols
  set_colors <- rep("black", length(gene_sets))
  if (!is.null(cfg$merged_cell_cols)) {
    for (i in seq_along(gene_sets)) {
      nm <- names(gene_sets)[i]
      if (nm %in% names(cfg$merged_cell_cols)) {
        set_colors[i] <- cfg$merged_cell_cols[[nm]]
      }
    }
  }

  p <- ggVennDiagram(gene_sets, label_alpha = 0, set_size = 10,
                     label_size = 10, edge_size = 3.5,
                     set_color = set_colors) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste0("Overlap of DE Genes (", paste(names(gene_sets), collapse = ", "), ")")) +
    theme(
      plot.title  = element_text(face = "bold", size = 23),
      plot.margin = margin(25, 100, 25, 50),
      legend.text      = element_text(size = 23),
      legend.title     = element_text(size = 23)
    )

  attr(p, "gene_sets") <- gene_sets
  p
}

# Read PNG width/height directly from the IHDR chunk — avoids a magick/png dep.
read_png_dims <- function(f) {
  con <- file(f, "rb")
  on.exit(close(con))
  readBin(con, "raw", n = 16)  # skip 8-byte PNG signature + 8-byte IHDR length/type
  w <- readBin(con, "integer", n = 1, size = 4, endian = "big")
  h <- readBin(con, "integer", n = 1, size = 4, endian = "big")
  c(width = w, height = h)
}

# "00a_qc_violin_by_IMP_sample.png" -> "qc violin by IMP sample"
plots_dir_slide_title <- function(fname) {
  base <- tools::file_path_sans_ext(basename(fname))
  base <- sub("^[0-9]+[a-z]?_", "", base)
  gsub("_", " ", base)
}

# Build a plain-white slide deck where each PNG in `plots_dir` becomes its own
# centered, aspect-preserved slide. No template, no theming. Slide dimensions
# follow the officer default master (10 x 7.5 in for the built-in Office Theme).
build_plots_dir_pptx <- function(plots_dir, out_pptx,
                                 include_title = TRUE) {
  suppressPackageStartupMessages(library(officer))

  if (!dir.exists(plots_dir)) {
    warning("Plots directory not found: ", plots_dir)
    return(NULL)
  }

  png_files <- sort(list.files(plots_dir, pattern = "\\.png$",
                               full.names = TRUE, ignore.case = TRUE))
  if (length(png_files) == 0) {
    warning("No PNG files found in plots directory: ", plots_dir)
    return(NULL)
  }

  doc <- read_pptx()
  sz  <- slide_size(doc)
  slide_width  <- sz$width
  slide_height <- sz$height

  title_h     <- if (include_title) 0.5 else 0
  margin_side <- 0.3
  margin_top  <- if (include_title) 0.15 else 0.2
  margin_bot  <- 0.2

  avail_w <- slide_width  - 2 * margin_side
  avail_h <- slide_height - title_h - margin_top - margin_bot

  for (f in png_files) {
    dims <- tryCatch(read_png_dims(f), error = function(e) c(width = NA, height = NA))
    aspect <- if (!any(is.na(dims)) && dims["height"] > 0) {
      unname(dims["width"] / dims["height"])
    } else {
      avail_w / avail_h
    }

    if (avail_w / avail_h > aspect) {
      h <- avail_h
      w <- aspect * h
    } else {
      w <- avail_w
      h <- w / aspect
    }
    left <- (slide_width - w) / 2
    top  <- title_h + margin_top + (avail_h - h) / 2

    doc <- add_slide(doc, layout = "Blank", master = "Office Theme")

    if (include_title) {
      doc <- ph_with(
        doc,
        value    = plots_dir_slide_title(f),
        location = ph_location(left = margin_side, top = 0.1,
                               width = avail_w, height = title_h)
      )
    }

    doc <- ph_with(
      doc,
      value    = external_img(f, width = w, height = h),
      location = ph_location(left = left, top = top, width = w, height = h)
    )
  }

  print(doc, target = out_pptx)
  message("Wrote slide deck with ", length(png_files), " image(s): ", out_pptx)
  invisible(out_pptx)
}

insert_image_on_slide <- function(doc, slide_index, image_path, box) {
  doc <- on_slide(doc, index = slide_index)
  ph_with(
    x = doc,
    value = external_img(image_path, width = unname(box["width"]), height = unname(box["height"])),
    location = ph_location(
      left = unname(box["left"]),
      top = unname(box["top"]),
      width = unname(box["width"]),
      height = unname(box["height"])
    )
  )
}

build_powerpoint <- function(cfg, fig_paths) {
  if (!file.exists(cfg$template_pptx)) {
    warning("Template PowerPoint not found: ", cfg$template_pptx)
    return(NULL)
  }

  doc <- read_pptx(cfg$template_pptx)

  mapping <- list(
    slide1 = c(fig3A="fig3A_umap_celltype", fig3B="fig3B_adt_featureplot", fig3C="fig3C_canonical_violins",
               fig3D="fig3D_umap_pain", fig3E="fig3E_cell_composition_boxplot"),
    slide2 = c(sfig3A="sfig3A_umap_pool", sfig3B="sfig3B_qc_pool", sfig3C="sfig3C_qc_sample",
               sfig3D="sfig3D_stacked_composition", sfig3E="sfig3E_denovo_marker_heatmap"),
    slide3 = c(fig4A="fig4A_pseudobulk_violins", fig4B="fig4B_gsea", fig4C="fig4C_gene_examples"),
    slide4 = c(sfig4A="sfig4A_de_heatmap", sfig4B="sfig4B_venn")
  )

  for (sl_name in names(mapping)) {
    sl_num <- as.integer(gsub("slide", "", sl_name))
    placements <- cfg$ppt$placements[[sl_name]]
    for (slot in names(mapping[[sl_name]])) {
      key <- mapping[[sl_name]][[slot]]
      img <- fig_paths[[key]]
      if (!is.null(img) && file.exists(img)) {
        doc <- insert_image_on_slide(doc, slide_index = sl_num, image_path = img, box = placements[[slot]])
      }
    }
  }

  out_ppt <- file.path(cfg$output_dir, "IMPACT_auto_filled.pptx")
  print(doc, target = out_ppt)
  out_ppt
}

# ---------------------------- #
# 3) MAIN DRIVER
# ---------------------------- #

run_impact_pipeline <- function(cfg) {
  dir_create2(cfg$output_dir)
  fig_dir <- file.path(cfg$output_dir, "Figures")
  tab_dir <- file.path(cfg$output_dir, "Tables")
  dir_create2(fig_dir)
  dir_create2(tab_dir)

  message2("Reading object: {cfg$input_rds}")
  obj <- readRDS(cfg$input_rds)
  validate_inputs(obj, cfg)
  obj <- standardize_metadata(obj, cfg)

  # Load covariate table (Age, Sex) if provided
  covariates <- NULL
  if (!is.null(cfg$covariate_csv) && file.exists(cfg$covariate_csv)) {
    covariates <- readr::read_csv(cfg$covariate_csv, show_col_types = FALSE)
    covariates <- as.data.frame(covariates)
    message2("Loaded covariates: {nrow(covariates)} samples, columns: {paste(colnames(covariates), collapse=', ')}")
  } else {
    message("No covariate table found — DESeq2 and propeller will run without Age/Sex covariates")
  }

  fig_paths <- list()

  # Figure 3A
  p_fig3A <- plot_umap_celltypes(obj, cfg)
  fig_paths$fig3A_umap_celltype <- save_plot_dual(p_fig3A, "fig3A_umap_celltype", fig_dir, 8, 6)$png

  # Figure 3B — ADT FeaturePlots (9 per page)
  adt_pages <- plot_adt_feature_pages(obj, cfg, per_page = 9, ncol = 3)
  if (!is.null(adt_pages)) {
    for (i in seq_along(adt_pages)) {
      fname <- paste0("fig3B_adt_featureplot_page", i)
      fig_paths[[fname]] <- save_plot_dual(adt_pages[[i]], fname, fig_dir, 12, 12)$png
    }
  }

  # ADT Violin Plots (9 per page)
  adt_vln_pages <- plot_adt_violin_pages(obj, cfg, per_page = 9, ncol = 3)
  if (!is.null(adt_vln_pages)) {
    for (i in seq_along(adt_vln_pages)) {
      fname <- paste0("adt_violin_page", i)
      fig_paths[[fname]] <- save_plot_dual(adt_vln_pages[[i]], fname, fig_dir, 14, 14)$png
    }
  }

  # Figure 3C
  p_fig3C <- plot_canonical_rna_violins(obj, cfg)
  fig_paths$fig3C_canonical_violins <- save_plot_dual(p_fig3C, "fig3C_canonical_violins", fig_dir, 12.5, 12)$png

  # Figure 3D
  p_fig3D <- plot_umap_pain(obj, cfg)
  fig_paths$fig3D_umap_pain <- save_plot_dual(p_fig3D, "fig3D_umap_pain", fig_dir, 8, 6)$png

  # Figure 3E
  prop_res <- run_propeller_limma(obj, cfg, covariates = covariates, transform = "asin")
  readr::write_csv(prop_res, file.path(tab_dir, "celltype_proportions_propeller.csv"))
  p_fig3E <- plot_composition_boxplot_with_legend(obj, cfg, prop_res,
                                                   title = "Cell Composition (propeller/limma)")
  fig_paths$fig3E_cell_composition_boxplot <- save_plot_dual(p_fig3E, "fig3E_cell_composition_boxplot", fig_dir, 16, 6)$png

  # Supplemental Figure 3A
  p_sfig3A <- plot_umap_pool(obj, cfg)
  fig_paths$sfig3A_umap_pool <- save_plot_dual(p_sfig3A, "sfig3A_umap_pool", fig_dir, 8, 6)$png

  # Supplemental Figure 3B
  p_sfig3B <- plot_qc_by_group(obj, cfg$meta$pool_col, "Single-cell QC by Pool")
  if (!is.null(p_sfig3B)) {
    fig_paths$sfig3B_qc_pool <- save_plot_dual(p_sfig3B, "sfig3B_qc_pool", fig_dir, 10, 16)$png
  }

  # Supplemental Figure 3C — QC by sample, colored by pool group
  p_sfig3C <- plot_qc_by_sample_pooled(obj, cfg, "Single-cell QC by Sample")
  if (!is.null(p_sfig3C)) {
    fig_paths$sfig3C_qc_sample <- save_plot_dual(p_sfig3C, "sfig3C_qc_sample", fig_dir, 14, 18)$png
  }

  # Supplemental Figure 3D
  p_sfig3D <- plot_stacked_cell_composition_by_sample(obj, cfg)
  fig_paths$sfig3D_stacked_composition <- save_plot_dual(p_sfig3D, "sfig3D_stacked_composition", fig_dir, 8, 6)$png

  # Supplemental Figure 3E — de novo markers (skipped: requires presto/FindAllMarkers)
  # denovo <- find_denovo_markers(obj, cfg, top_n = 4)
  # readr::write_csv(denovo, file.path(tab_dir, "denovo_markers_top10.csv"))
  # p_sfig3E <- plot_denovo_marker_heatmap(obj, cfg, denovo, top_n = 4)
  # fig_paths$sfig3E_denovo_marker_heatmap <- save_heatmap_direct(p_sfig3E, "sfig3E_denovo_marker_heatmap", fig_dir, 12, 14)$png

  # ══════════════════════════════════════════════════════════════════════════
  # PSEUDOBULK WORKFLOW — Dual: 11 cell types (PCA) + 6 merged (DE/GSEA/Venn)
  # ══════════════════════════════════════════════════════════════════════════

  # Add MergedCellType column
  obj <- add_merged_celltype(obj, cfg)

  # ── PCA: 11 cell types ──────────────────────────────────────────────────
  message2("Generating pseudobulk PCA (all 11 cell types)...")
  pca_11_res <- plot_pseudobulk_pca(obj, cfg, ct_col_use = cfg$meta$celltype_col,
                                     title_suffix = " (11 cell types)")
  if (!is.null(pca_11_res)) {
    fig_paths$pca_11ct <- save_plot_dual(pca_11_res$plot, "pseudobulk_PCA_11celltypes", fig_dir, 10, 8)$png
    readr::write_csv(pca_11_res$pca_df, file.path(tab_dir, "pseudobulk_PCA_11ct_coordinates.csv"))
    readr::write_csv(
      as.data.frame(pca_11_res$log2cpm) %>% rownames_to_column("Gene"),
      file.path(tab_dir, "pseudobulk_log2CPM_11ct.csv")
    )
  }

  # ── PCA: 6 merged cell types ───────────────────────────────────────────
  message2("Generating pseudobulk PCA (6 merged cell types)...")
  pca_6_res <- plot_pseudobulk_pca(obj, cfg, ct_col_use = "MergedCellType",
                                    title_suffix = " (6 merged cell types)")
  if (!is.null(pca_6_res)) {
    fig_paths$pca_6ct <- save_plot_dual(pca_6_res$plot, "pseudobulk_PCA_6celltypes", fig_dir, 10, 8)$png
    readr::write_csv(pca_6_res$pca_df, file.path(tab_dir, "pseudobulk_PCA_6ct_coordinates.csv"))
    readr::write_csv(
      as.data.frame(pca_6_res$log2cpm) %>% rownames_to_column("Gene"),
      file.path(tab_dir, "pseudobulk_log2CPM_6ct.csv")
    )
  }

  # ── Pseudobulk table per 11 cell types ─────────────────────────────────
  message2("Writing pseudobulk tables (11 cell types)...")
  for (ct in names(cfg$cell_cols)) {
    if (!(ct %in% unique(as.character(obj@meta.data[[cfg$meta$celltype_col]])))) next
    tryCatch({
      pb <- pseudo_bulk_counts(obj, cfg, ct, covariates = covariates)
      readr::write_csv(
        as.data.frame(pb$counts) %>% rownames_to_column("Gene"),
        file.path(tab_dir, paste0("pseudobulk_", safe_name(ct), ".csv"))
      )
    }, error = function(e) message("  Skipping pseudobulk table for ", ct, ": ", conditionMessage(e)))
  }

  # ── Cell composition: 3 methods ────────────────────────────────────────

  # Method 1: Propeller/limma (already done above — prop_res)

  # Method 2: T-test
  message2("Running t-test on cell type proportions...")
  ttest_res <- run_ttest_proportions(obj, cfg)
  readr::write_csv(ttest_res, file.path(tab_dir, "celltype_proportions_ttest.csv"))
  p_comp_ttest <- plot_composition_boxplot_with_legend(obj, cfg, ttest_res,
                                                       title = "Cell Composition (t-test)")
  fig_paths$fig3E_composition_ttest <- save_plot_dual(p_comp_ttest, "fig3E_composition_ttest", fig_dir, 14, 8)$png

  # Method 3: MASC
  message2("Running MASC on cell type proportions...")
  masc_res <- tryCatch(run_masc(obj, cfg), error = function(e) {
    message("  !! MASC failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(masc_res)) {
    readr::write_csv(masc_res, file.path(tab_dir, "celltype_proportions_MASC.csv"))
    p_comp_masc <- plot_composition_boxplot_with_legend(obj, cfg, masc_res,
                                                        title = "Cell Composition (MASC)")
    fig_paths$fig3E_composition_masc <- save_plot_dual(p_comp_masc, "fig3E_composition_MASC", fig_dir, 16, 6)$png
  }

  # ── DESeq2 on 6 MERGED cell types ─────────────────────────────────────
  merged_cts <- names(cfg$merged_celltypes)
  de_results <- list()

  # Temporarily swap celltype col to MergedCellType for pseudobulk
  orig_ct_col <- cfg$meta$celltype_col
  cfg$meta$celltype_col <- "MergedCellType"

  for (ct in merged_cts) {
    if (!(ct %in% unique(as.character(obj@meta.data[["MergedCellType"]])))) next
    message2("Running pseudobulk DESeq2 (merged): {ct}")
    one <- tryCatch(
      run_deseq_one_celltype(obj, cfg, ct, covariates = covariates),
      error = function(e) {
        message("  !! SKIPPING ", ct, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(one)) {
      de_results[[ct]] <- one
      readr::write_csv(one$res, file.path(tab_dir, paste0("DESeq2_merged_", safe_name(ct), "_High_vs_LowPain.csv")))
      message2("  Saved DESeq2 CSV for {ct} ({nrow(one$res)} genes)")
      volc <- plot_volcano(one$res, title = ct,
                           pval_cutoff = cfg$de$pval_cutoff,
                           lfc_cutoff = cfg$de$lfc_cutoff,
                           top_n = cfg$de$top_n_labels)
      save_plot_dual(volc, paste0("volcano_merged_", safe_name(ct)), fig_dir, 8, 5)
    }
  }

  # Restore original celltype col
  cfg$meta$celltype_col <- orig_ct_col

  message2("DESeq2 (merged) complete: {length(de_results)}/{length(merged_cts)} cell types succeeded")

  # ── Combined volcano grid (all cell types, shared legend, 2:1 ratio) ────────
  if (length(de_results) > 0) {
    p_volc_grid <- plot_volcano_grid(de_results, cfg)
    if (!is.null(p_volc_grid)) {
      # 2:1 width:height — 18 × 9 inches fits a 3×2 panel neatly
      fig_paths$volcano_grid <- save_plot_dual(
        p_volc_grid, "volcano_merged_grid", fig_dir, width = 21, height = 9
      )$png
      message2("  Saved combined volcano grid")
    }
  }

  # ── Figure 4A pseudobulk boxplots ──────────────────────────────────────
  if (length(de_results) > 0) {
    plist <- lapply(names(de_results)[seq_len(min(4, length(de_results)))], function(ct) {
      plot_top_de_gene_boxplots(de_results[[ct]]$res, de_results[[ct]]$pb, ct, cfg, top_n = 2)
    })
    plist <- plist[!vapply(plist, is.null, logical(1))]
    if (length(plist) > 0) {
      p_fig4A <- wrap_plots(plist, ncol = 1) + plot_annotation(title = "Pseudobulk DESeq2: High/Low Pain")
      fig_paths$fig4A_pseudobulk_violins <- save_plot_dual(p_fig4A, "fig4A_pseudobulk_violins", fig_dir, 9, 12)$png
    }
  }

  # ── Figure 4B — 16-violin grid of top DEGs ─────────────────────────────
  if (length(de_results) > 0) {
    # Use MergedCellType for the violin grid
    cfg_merged <- cfg
    cfg_merged$meta$celltype_col <- "MergedCellType"
    p_fig4B <- plot_deg_violin_grid(obj, cfg_merged, de_results, n_genes = 16, ncol = 4)
    if (!is.null(p_fig4B)) {
      fig_paths$fig4B_deg_violin_grid <- save_plot_dual(p_fig4B, "fig4B_deg_violin_grid", fig_dir, 14, 14)$png
    }
  }

  # ── Figure 4C — GSEA Hallmark (skipped) ────────────────────────────────
  gsea_all <- NULL
  # if (length(de_results) > 0) {
  #   message2("Running GSEA (Hallmark) across {length(de_results)} merged cell types...")
  #   gsea_all <- tryCatch(
  #     run_gsea_all_celltypes(de_results, species = "Homo sapiens", category = "H"),
  #     error = function(e) { message("  !! GSEA failed: ", conditionMessage(e)); NULL }
  #   )
  #   if (!is.null(gsea_all) && nrow(gsea_all) > 0) {
  #     readr::write_csv(
  #       gsea_all %>%
  #         mutate(leading_edge_genes = sapply(leadingEdge, paste, collapse = ", ")) %>%
  #         select(-leadingEdge),
  #       file.path(tab_dir, "GSEA_Hallmark_merged_celltypes.csv")
  #     )
  #     message2("  Saved GSEA CSV ({nrow(gsea_all)} pathway-celltype results)")
  #     p_gsea_heat <- plot_gsea_heatmap_balanced(gsea_all, top_n_per_ct = 4)
  #     if (!is.null(p_gsea_heat)) {
  #       fig_paths$fig4C_gsea_heatmap <- save_plot_dual(p_gsea_heat, "fig4C_gsea_heatmap", fig_dir, 12, 10)$png
  #     }
  #     p_gsea_dot <- plot_gsea_dotplot(gsea_all, top_n_per_ct = 1, min_per_ct = 1,
  #                                     force_pathways = c(
  #                                       "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  #                                       "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  #                                       "HALLMARK_ANGIOGENESIS",
  #                                       "HALLMARK_MYOGENESIS",
  #                                       "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  #                                       "HALLMARK_APOPTOSIS",
  #                                       "HALLMARK_INFLAMMATORY_RESPONSE",
  #                                       "HALLMARK_WNT_BETA_CATENIN_SIGNALING"))
  #     if (!is.null(p_gsea_dot)) {
  #       fig_paths$fig4C_gsea_dotplot <- save_plot_dual(p_gsea_dot, "fig4C_gsea_dotplot", fig_dir, 14, 6)$png
  #     }
  #     le_patterns <- c("TNFA", "TGF_BETA", "ANGIOGENESIS", "WNT")
  #     le_df <- tryCatch(
  #       get_leading_edge(gsea_all, le_patterns, sig_cutoff = 0.05),
  #       error = function(e) { message("  !! Leading edge extraction failed: ", conditionMessage(e)); NULL }
  #     )
  #     if (!is.null(le_df) && nrow(le_df) > 0) {
  #       readr::write_csv(le_df %>% select(-leadingEdge),
  #                        file.path(tab_dir, "GSEA_leading_edge_TNF_TGF_Angio_WNT.csv"))
  #       le_flat <- le_df %>%
  #         group_by(pathway) %>%
  #         summarize(n_sig_celltypes = n(), sig_celltypes = paste(CellType, collapse = ", "),
  #                   all_leading_edge_genes = paste(unique(unlist(leadingEdge)), collapse = ", "),
  #                   .groups = "drop")
  #       readr::write_csv(le_flat,
  #                        file.path(tab_dir, "GSEA_leading_edge_TNF_TGF_Angio_WNT_flat.csv"))
  #     }
  #   }
  # }

  # ── Figure 4D gene examples ────────────────────────────────────────────
  if (length(de_results) > 0) {
    ct_use <- names(de_results)[1]
    p_fig4D <- plot_top_de_gene_boxplots(de_results[[ct_use]]$res, de_results[[ct_use]]$pb, ct_use, cfg, top_n = 4)
    if (!is.null(p_fig4D)) {
      fig_paths$fig4D_gene_examples <- save_plot_dual(p_fig4D, "fig4D_gene_examples", fig_dir, 8, 6)$png
    }
  }

  # ── Figure 5A — highlight gene boxplots (9 per page) ───────────────────
  hl_pages <- plot_highlight_gene_boxplots(obj, cfg,
                                            de_results  = de_results,
                                            pval_cutoff = cfg$de$pval_cutoff)
  if (!is.null(hl_pages)) {
    for (i in seq_along(hl_pages)) {
      fname <- paste0("fig5A_highlight_gene_boxplots_page", i)
      fig_paths[[fname]] <- save_plot_dual(hl_pages[[i]], fname, fig_dir, 24, 14)$png
    }
  }

  # ── Figure 5B — highlight genes DotPlot ────────────────────────────────
  p_hl_dot <- plot_highlight_gene_dotplot(obj, cfg)
  if (!is.null(p_hl_dot)) {
    fig_paths$fig5B_highlight_dotplot <- save_plot_dual(p_hl_dot, "fig5B_highlight_gene_dotplot", fig_dir, 12, 10)$png
  }

  # ── Figure 5B2 — highlight genes, 6 merged cell types + t-test ─────────
  hl_merged_pages <- plot_highlight_gene_merged(obj, cfg)
  if (!is.null(hl_merged_pages)) {
    for (i in seq_along(hl_merged_pages)) {
      fname <- paste0("fig5B2_highlight_merged_page", i)
      fig_paths[[fname]] <- save_plot_dual(hl_merged_pages[[i]], fname, fig_dir, 14, 14)$png
    }
  }

  # ── Figure 5C — DE bubble summary ──────────────────────────────────────
  if (length(de_results) > 0) {
    p_bubble <- plot_de_bubble_summary(de_results,
                                       pval_cutoff = cfg$de$pval_cutoff,
                                       lfc_cutoff = cfg$de$lfc_cutoff)
    if (!is.null(p_bubble)) {
      fig_paths$fig5C_de_bubble <- save_plot_dual(p_bubble, "fig5C_de_bubble_summary", fig_dir, 10, 6)$png
    }
  }

  # ── Figure 5D — DE gene count barplot ──────────────────────────────────
  if (length(de_results) > 0) {
    p_de_bar <- plot_de_count_boxplot(de_results,
                                      pval_cutoff = cfg$de$pval_cutoff,
                                      lfc_cutoff = cfg$de$lfc_cutoff)
    if (!is.null(p_de_bar)) {
      fig_paths$fig5D_de_barplot <- save_plot_dual(p_de_bar, "fig5D_de_count_barplot", fig_dir, 8, 6)$png
    }
  }

  # ── Supplemental Figure 4A — DE heatmap ────────────────────────────────
  p_sfig4A <- plot_de_heatmap(de_results,
                               cfg         = cfg,
                               pval_cutoff = cfg$de$pval_cutoff,
                               lfc_cutoff  = cfg$de$lfc_cutoff)
  if (!is.null(p_sfig4A)) {
    fig_paths$sfig4A_de_heatmap <- save_plot_dual(p_sfig4A, "sfig4A_de_heatmap", fig_dir, 30, 11)$png
  }

  # ── Supplemental Figure 4B — UpSet plot ────────────────────────────────
  if (length(de_results) > 0 && requireNamespace("UpSetR", quietly = TRUE)) {
    upset_data <- build_upset_data(de_results,
                                    pval_cutoff = cfg$de$pval_cutoff,
                                    lfc_cutoff = cfg$de$lfc_cutoff)
    if (!is.null(upset_data)) {
      upset_png <- file.path(fig_dir, "sfig4B_upset_de_overlap.png")
      save_upset_plot(upset_data, out_png = upset_png)
      fig_paths$sfig4B_upset <- upset_png
      upset_df <- lapply(names(upset_data$gene_sets), function(ct) {
        tibble(CellType = ct, Gene = upset_data$gene_sets[[ct]])
      }) %>% bind_rows()
      readr::write_csv(upset_df, file.path(tab_dir, "upset_de_gene_sets.csv"))
    }
  }

  # ── Supplemental Figure 4C — Merged Venn (6 groups) ────────────────────
  if (requireNamespace("ggVennDiagram", quietly = TRUE) && length(de_results) > 0) {
    p_sfig4C <- plot_merged_venn(de_results, cfg,
                                 pval_cutoff = cfg$de$pval_cutoff,
                                 lfc_cutoff = cfg$de$lfc_cutoff)
    if (!is.null(p_sfig4C)) {
      fig_paths$sfig4C_venn <- save_plot_dual(p_sfig4C, "sfig4C_merged_venn", fig_dir, 18, 15)$png
    }
  }

  # Build PPT
  out_ppt <- NULL
  if (isTRUE(cfg$make_pptx)) {
    out_ppt <- build_powerpoint(cfg, fig_paths)
  }

  # Build the plain-white deck from the pre-rendered plots directory
  plots_ppt <- NULL
  if (isTRUE(cfg$make_plots_pptx)) {
    plots_ppt <- build_plots_dir_pptx(cfg$plots_input_dir, cfg$plots_pptx_path)
  }

  list(
    figures = fig_paths,
    de_results = de_results,
    gsea = gsea_all,
    propeller = prop_res,
    ttest = ttest_res,
    masc = masc_res,
    pptx = out_ppt,
    plots_pptx = plots_ppt
  )
}

# ---------------------------- #
# 4) CLI ENTRYPOINT
# ---------------------------- #

if (sys.nframe() == 0) {
  # When sourced interactively, edit cfg above then call:
  #   res <- run_impact_pipeline(cfg)
  #
  # When run from CLI directly:
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    cfg$input_rds <- args[1]
    if (length(args) >= 2) cfg$output_dir <- args[2]
    if (length(args) >= 3) cfg$covariate_csv <- args[3]

    message2("Starting IMPACT pipeline...")
    message2("  Input:      {cfg$input_rds}")
    message2("  Output:     {cfg$output_dir}")
    message2("  Covariates: {cfg$covariate_csv}")

    res <- run_impact_pipeline(cfg)
    message2("Pipeline complete. Results saved to {cfg$output_dir}/")
  } else {
    message2("IMPACT figure pipeline loaded.")
    message2("Usage:")
    message2("  CLI:         Rscript IMPACT_one_shot_figure_pipeline.R <input.rds> [output_dir] [covariate.csv]")
    message2("  Interactive: edit cfg$input_rds, then run res <- run_impact_pipeline(cfg)")
  }
}
