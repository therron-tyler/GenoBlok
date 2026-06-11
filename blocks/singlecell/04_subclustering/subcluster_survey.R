#!/usr/bin/env Rscript
## NPSLE T-cell subcluster survey pipeline
## Merges dims_resolution_scan.R + adt_umap_survey.R
##
## Usage:
##   Rscript subcluster_survey.R \
##     --rds     /path/to/subcluster_final.rds \
##     --outdir  /path/to/output \
##     [--dims     10,15,20]        # max PCs per set; each set = dims 1:N
##     [--res      0.2,0.4,0.6]     # resolutions for per-combo QC pages
##     [--features CD3D,CD4,...]    # RNA canonical markers (comma-separated)
##     [--hto_col  HTO_maxID]       # metadata column for donor split in ADT violins
##
## Output per object: one PDF per dims value → TcellSurvey_<obj>_d<N>.pdf
## PDF page order (per resolution):
##   RNA QC p1 : clustree (res 0.05–1.0, current highlighted) | top-5 DE heatmap
##   RNA QC p2 : seurat UMAP + cells-per-cluster table | canonical violins (paginated)
##   RNA QC p3+: canonical feature plots (paginated, 9 per page)
##   ADT pages : one per T-cell ADT — [cluster UMAP | ADT UMAP | RNA UMAP] / violin HTO×cluster

.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
  library(optparse)
  library(clustree)
  library(RColorBrewer)
  library(gridExtra)
})

# =============================================================================
#  CONSTANTS
# =============================================================================

TCELL_ADT_MARKERS <- c(
  # Core identity / TCR — these pages appear first
  "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8",
  "Hu.CD2",        "Hu.CD5",         "Hu.CD7",
  "Hu.TCR.AB",     "Hu.TCR.Vd2",     "Hu.TCR.Va7.2",
  # Naive / memory
  "Hu.CD45RA",  "Hu.CD45RO",  "Hu.CD62L",
  "Hu.CD27",    "Hu.CD28",    "Hu.CD95",
  "Hu.CD127",   "Hu.CD122",   "HuMs.CD44",
  # Treg
  "Hu.CD25",
  # Effector / cytotoxic
  "Hu.CD57",    "Hu.CD56",    "Hu.CD94",
  "Hu.CX3CR1",  "Hu.KLRG1",   "Hu.CD161",
  # Activation
  "Hu.CD69",    "Hu.CD38-HIT2", "Hu.HLA.DR",
  "Hu.CD26",    "Hu.CD137",     "Hu.CD154",
  # Exhaustion / checkpoint
  "Hu.CD223",   "Hu.CD279",   "Hu.TIGIT",
  "Hu.CD152",   "Hu.CD244",
  # Homing / trafficking
  "Hu.CD183",   "Hu.CD185",   "Hu.CD194",
  "Hu.CD195",   "Hu.CD196",
  # Co-stimulation / misc
  "Hu.CD49d",   "Hu.CD58",    "Hu.CD52"
)

TCELL_RNA_DEFAULT <-
  "CD3D,CD4,CD8A,CCR7,SELL,IL7R,GZMB,GZMK,FOXP3,MKI67,LAG3,PDCD1,ISG15,CX3CR1,NKG7"

# T-cell-relevant subset of ADT→RNA gene map
ADT_RNA_MAP <- c(
  "Hu.CD3-UCHT1"  = "CD3E",   "Hu.CD4-RPA.T4" = "CD4",
  "Hu.CD8"        = "CD8A",   "Hu.CD2"         = "CD2",
  "Hu.CD5"        = "CD5",    "Hu.CD7"         = "CD7",
  "Hu.CD28"       = "CD28",   "Hu.CD27"        = "CD27",
  "Hu.CD25"       = "IL2RA",  "Hu.CD127"       = "IL7R",
  "Hu.CD69"       = "CD69",   "HuMs.CD44"      = "CD44",
  "Hu.CD45RA"     = "PTPRC",  "Hu.CD45RO"      = "PTPRC",
  "Hu.CD62L"      = "SELL",   "Hu.CD161"       = "KLRB1",
  "Hu.CD57"       = "B3GAT1", "Hu.KLRG1"       = "KLRG1",
  "Hu.TCR.AB"     = "TRAC",   "Hu.TCR.Va7.2"   = "TRAV1-2",
  "Hu.TCR.Vd2"    = "TRDV2",  "Hu.CD279"       = "PDCD1",
  "Hu.CD223"      = "LAG3",   "Hu.TIGIT"       = "TIGIT",
  "Hu.CD152"      = "CTLA4",  "Hu.CD244"       = "CD244",
  "Hu.CD137"      = "TNFRSF9","Hu.CD154"       = "CD40LG",
  "Hu.CD183"      = "CXCR3",  "Hu.CD185"       = "CXCR5",
  "Hu.CD194"      = "CCR4",   "Hu.CD195"       = "CCR5",
  "Hu.CD196"      = "CCR6",   "Hu.CD49d"       = "ITGA4",
  "Hu.CD58"       = "CD58",   "Hu.CD52"        = "CD52",
  "Hu.CD56"       = "NCAM1",  "Hu.CD94"        = "KLRD1",
  "Hu.CX3CR1"     = "CX3CR1", "Hu.CD38-HIT2"   = "CD38",
  "Hu.HLA.DR"     = "HLA-DRA","Hu.CD26"        = "DPP4",
  "Hu.CD122"      = "IL2RB",  "Hu.CD95"        = "FAS"
)

CLUSTREE_RES_STEPS <- seq(0.05, 1.0, by = 0.05)  # 20 levels

# =============================================================================
#  CLI
# =============================================================================

option_list <- list(
  make_option("--rds",      type = "character", default = NULL,
              help = "Path to .rds file, named-list .rds, or directory of .rds files"),
  make_option("--outdir",   type = "character", default = NULL,
              help = "Output directory"),
  make_option("--dims",     type = "character", default = "10,15,20",
              help = "Comma-separated max PCs, e.g. '10,15,20' means dims 1:10, 1:15, 1:20 [default: %default]"),
  make_option("--res",      type = "character", default = "0.2,0.4,0.6",
              help = "Comma-separated clustering resolutions for QC pages [default: %default]"),
  make_option("--features", type = "character", default = TCELL_RNA_DEFAULT,
              help = "Comma-separated RNA canonical markers [default: T-cell panel]"),
  make_option("--hto_col",  type = "character", default = "HTO_maxID",
              help = "Metadata column for donor split in ADT violins [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rds) || is.null(opt$outdir))
  stop("--rds and --outdir are required. Run with --help for usage.")

split_csv <- function(x) trimws(strsplit(trimws(x), ",")[[1]])

dims_vec    <- as.integer(split_csv(opt$dims))
resolutions <- as.numeric(split_csv(opt$res))
features    <- split_csv(opt$features)

if (anyNA(dims_vec))    stop("Non-integer value in --dims: ",    opt$dims)
if (anyNA(resolutions)) stop("Non-numeric value in --res: ",     opt$res)

# =============================================================================
#  RDS LOADER
# =============================================================================

load_rds_input <- function(rds_path) {
  if (dir.exists(rds_path)) {
    files <- list.files(rds_path, pattern = "\\.rds$",
                        full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0) stop("No .rds files found in: ", rds_path)
    objs <- lapply(files, readRDS)
    names(objs) <- sub("\\.rds$", "", basename(files), ignore.case = TRUE)
    message("Loaded ", length(objs), " RDS file(s) from directory")
    return(objs)
  }
  if (!file.exists(rds_path)) stop("File not found: ", rds_path)
  obj <- readRDS(rds_path)
  if (is.list(obj) && !inherits(obj, "Seurat")) {
    if (is.null(names(obj)) || any(names(obj) == ""))
      stop("List RDS must be a named list of Seurat objects.")
    message("Loaded named list with ", length(obj), " Seurat object(s)")
    return(obj)
  }
  if (inherits(obj, "Seurat")) {
    nm <- sub("\\.rds$", "", basename(rds_path), ignore.case = TRUE)
    message("Loaded Seurat object (", ncol(obj), " cells) as '", nm, "'")
    return(setNames(list(obj), nm))
  }
  stop("Unsupported RDS content.")
}

# =============================================================================
#  THEME HELPERS
# =============================================================================

umap_base_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      aspect.ratio      = 1,
      axis.line         = element_blank(),
      axis.ticks        = element_blank(),
      axis.text         = element_blank(),
      axis.title        = element_text(size = base_size - 1, color = "grey40"),
      plot.title        = element_text(size = base_size, face = "bold",
                                       margin = margin(b = 4)),
      legend.text       = element_text(size = base_size - 2),
      legend.title      = element_text(size = base_size - 2),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.35, "cm")
    )
}

blank_panel <- function(msg = "No RNA\ncounterpart") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 5, color = "grey55", hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(aspect.ratio    = 1,
          plot.background = element_rect(fill = "grey97", color = "grey85",
                                         linewidth = 0.3))
}

# =============================================================================
#  ADT CLR NORMALIZATION (margin = 2)
# =============================================================================

normalize_adt_clr2 <- function(obj) {
  if (!"ADT" %in% Assays(obj)) {
    warning("No ADT assay — skipping normalization")
    return(obj)
  }
  obj <- NormalizeData(obj, assay = "ADT",
                        normalization.method = "CLR",
                        margin = 2, verbose = FALSE)
  message("  ADT CLR normalization (margin=2) applied")
  obj
}

# =============================================================================
#  CLUSTREE
# =============================================================================

build_clustree_metadata <- function(obj, snn_graph, dims_i) {
  ct_prefix <- paste0("ct_d", dims_i, "_res.")
  message("  Clustree: running FindClusters at ",
          length(CLUSTREE_RES_STEPS), " resolution steps...")
  for (r in CLUSTREE_RES_STEPS) {
    col_name <- paste0(ct_prefix, sprintf("%.2f", r))
    obj <- FindClusters(obj, graph.name  = snn_graph,
                         resolution  = r,
                         cluster.name = col_name,
                         verbose     = FALSE)
  }
  list(obj = obj, prefix = ct_prefix)
}

make_clustree_plot <- function(ct_info, current_res) {
  p <- clustree(ct_info$obj, prefix = ct_info$prefix) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.text   = element_text(size = 8),
      legend.title  = element_text(size = 9)
    )

  # Horizontal reference lines at each resolution level.
  # clustree (via ggraph tree layout) places resolution levels at integer
  # y-positions 0 (lowest res) to n-1 (highest). Adjust yintercept values
  # if your ggraph version uses 1-based or reversed ordering.
  n_lvl    <- length(CLUSTREE_RES_STEPS)
  curr_idx <- which(abs(CLUSTREE_RES_STEPS - current_res) < 1e-6) - 1L  # 0-based

  for (i in 0L:(n_lvl - 1L)) {
    p <- p + geom_hline(
      yintercept = i,
      color      = "grey75",
      linetype   = "dashed",
      linewidth  = 0.25
    )
  }

  p + labs(
    title    = "Clustree  (res 0.05 → 1.0, step 0.05)",
    subtitle = paste0("Current resolution: ", current_res)
  )
}

# =============================================================================
#  CELLS-PER-CLUSTER TABLE
# =============================================================================

cells_table_plot <- function(obj, cluster_col) {
  tbl        <- as.data.frame(table(Cluster = as.character(obj@meta.data[[cluster_col]])))
  tbl$Cluster <- factor(tbl$Cluster,
                         levels = as.character(sort(as.numeric(as.character(tbl$Cluster)))))
  tbl        <- tbl[order(tbl$Cluster), ]
  tbl$Pct    <- paste0(round(100 * tbl$Freq / sum(tbl$Freq), 1), "%")
  names(tbl) <- c("Cluster", "N cells", "% total")

  grob <- tableGrob(
    tbl, rows = NULL,
    theme = ttheme_minimal(
      base_size = 11,
      core    = list(fg_params  = list(fontsize = 11)),
      colhead = list(fg_params  = list(fontsize = 12, fontface = "bold"))
    )
  )
  wrap_elements(grob) +
    plot_annotation(
      title = "Cells per cluster",
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
}

# =============================================================================
#  RNA QC PAGES  — prints directly to the open PDF device
# =============================================================================

print_rna_qc_pages <- function(obj, ct_info, umap_name, cluster_col,
                                features, current_res, dims_i, obj_nm,
                                hto_col = NULL, hto_palette = NULL) {
  lbl      <- paste0(obj_nm, "  |  dims 1:", dims_i, "  |  res ", current_res)
  genes_ok <- intersect(features, rownames(obj[["RNA"]]))
  genes_ms <- setdiff(features,   rownames(obj[["RNA"]]))
  if (length(genes_ms) > 0)
    message("  RNA features missing: ", paste(genes_ms, collapse = ", "))

  DefaultAssay(obj) <- "RNA"
  Idents(obj)       <- cluster_col

  # ── FindAllMarkers for heatmap ─────────────────────────────────────────────
  message("  FindAllMarkers (res ", current_res, ")...")
  markers <- tryCatch(
    FindAllMarkers(obj, assay = "RNA", only.pos = TRUE,
                   min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE),
    error = function(e) { warning("FindAllMarkers failed: ", e$message); NULL }
  )

  # ── RNA QC Page 1: clustree (left) | top-5 heatmap (right) ────────────────
  p_ctree <- make_clustree_plot(ct_info, current_res)

  if (!is.null(markers) && nrow(markers) > 0) {
    top5 <- markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
      ungroup()
    p_heat <- DoHeatmap(obj,
                         features  = unique(top5$gene),
                         group.by  = cluster_col,
                         assay     = "RNA",
                         size      = 3.5,
                         raster    = TRUE) +
      theme(axis.text.y = element_text(size = 7),
            plot.title  = element_text(size = 12, face = "bold")) +
      ggtitle("Top 5 DE markers per cluster  (avg_log2FC)")
  } else {
    p_heat <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "FindAllMarkers\nproduced no results",
               size = 6, color = "grey55") +
      theme_void() +
      ggtitle("Heatmap unavailable")
  }

  page1 <- (p_ctree | p_heat) +
    plot_layout(widths = c(1.2, 1)) +
    plot_annotation(
      title    = paste0(lbl, "  —  RNA QC  [page 1 of 3]"),
      subtitle = "Left: clustree across res 0.05–1.0  |  Right: top 5 markers per cluster",
      theme    = theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size =  9, color = "grey40")
      )
    )
  print(page1)

  # ── RNA QC Page 2: UMAP + cells table + canonical violins ─────────────────
  n_cl   <- length(unique(obj@meta.data[[cluster_col]]))
  p_umap <- DimPlot(obj, reduction = umap_name, group.by = cluster_col,
                     label = TRUE, repel = TRUE, label.size = 4.5, pt.size = 0.5) +
    ggtitle(paste0("Seurat clusters  (res ", current_res,
                   ",  n = ", n_cl, " clusters)")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right",
          plot.title      = element_text(size = 13, face = "bold"),
          aspect.ratio    = 1)

  p_tbl <- cells_table_plot(obj, cluster_col)

  # HTO UMAP — shows donor contribution to each cluster
  if (!is.null(hto_col) && hto_col %in% colnames(obj@meta.data)) {
    p_hto <- DimPlot(obj, reduction = umap_name, group.by = hto_col,
                      label = FALSE, pt.size = 0.5,
                      cols = if (!is.null(hto_palette)) hto_palette else NULL) +
      ggtitle(paste0(hto_col, "  (donor)")) +
      theme_classic(base_size = 12) +
      theme(legend.position = "right",
            plot.title      = element_text(size = 13, face = "bold"),
            aspect.ratio    = 1)
    top_row <- (p_umap | p_hto | p_tbl) + plot_layout(widths = c(1.4, 1.4, 1))
  } else {
    top_row <- (p_umap | p_tbl) + plot_layout(widths = c(1.4, 1))
  }

  if (length(genes_ok) > 0) {
    vln_list <- VlnPlot(obj, features = genes_ok,
                         group.by = cluster_col, pt.size = 0, combine = FALSE)
    vln_list <- lapply(seq_along(vln_list), function(k)
      vln_list[[k]] +
        ggtitle(genes_ok[k]) +
        theme(legend.position = "none",
              plot.title      = element_text(size = 12, face = "bold"),
              axis.title.x    = element_blank(),
              axis.title.y    = element_text(size  =  9),
              axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
              axis.text.y     = element_text(size  =  9)))

    vln_chunks <- split(vln_list, ceiling(seq_along(vln_list) / 9))

    vln_wrap1 <- wrap_plots(vln_chunks[[1]], ncol = 3)
    page2 <- (top_row / vln_wrap1) +
      plot_layout(heights = c(1.3, 1.7)) +
      plot_annotation(
        title    = paste0(lbl, "  —  RNA QC  [page 2 of 3]"),
        subtitle = "Top: cluster UMAP | HTO donor UMAP | cells table  |  Bottom: canonical marker violins by cluster",
        theme    = theme(
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size =  9, color = "grey40")
        )
      )
    print(page2)

    for (k in seq_along(vln_chunks)[-1]) {
      p_cont <- wrap_plots(vln_chunks[[k]], ncol = 3) +
        plot_annotation(
          title = paste0(lbl, "  —  RNA violins (cont. ",
                         k, "/", length(vln_chunks), ")"),
          theme = theme(plot.title = element_text(size = 11, face = "bold"))
        )
      print(p_cont)
    }
  } else {
    page2 <- top_row +
      plot_annotation(title = paste0(lbl, "  —  RNA QC [page 2] — no features found"))
    print(page2)
  }

  # ── RNA QC Page 3+: feature plots ─────────────────────────────────────────
  if (length(genes_ok) > 0) {
    DefaultAssay(obj) <- "RNA"
    fp_list <- FeaturePlot(obj, features = genes_ok, reduction = umap_name,
                            pt.size = 2.0, order = TRUE,
                            combine = FALSE, raster = TRUE)
    fp_list <- lapply(seq_along(fp_list), function(k)
      fp_list[[k]] +
        ggtitle(genes_ok[k]) +
        theme(plot.title = element_text(size = 12, face = "bold")))

    fp_chunks <- split(fp_list, ceiling(seq_along(fp_list) / 9))
    for (k in seq_along(fp_chunks)) {
      pg_sfx <- if (length(fp_chunks) > 1)
        paste0(", pg ", k, "/", length(fp_chunks)) else ""
      p_fp <- wrap_plots(fp_chunks[[k]], ncol = 3) +
        plot_annotation(
          title = paste0(lbl, "  —  RNA Feature Plots  [page 3 of 3", pg_sfx, "]"),
          theme = theme(plot.title = element_text(size = 11, face = "bold"))
        )
      print(p_fp)
    }
  }

  invisible(NULL)
}

# =============================================================================
#  SAMPLE-ORIGIN VIOLIN PAGES  — prints directly to the open PDF device
#  One page per 6 canonical RNA markers; each panel faceted by cluster,
#  x = sample_origin (paired A/B order), fill = burnt-orange (_A) / royal-blue (_B)
# =============================================================================

print_sample_origin_violin_pages <- function(obj, cluster_col, features,
                                              dims_i, current_res, obj_nm,
                                              sample_col      = "sample_origin",
                                              markers_per_page = 6L) {
  if (!sample_col %in% colnames(obj@meta.data)) {
    warning("Column '", sample_col, "' not found — skipping sample-origin violin pages")
    return(invisible(NULL))
  }

  genes_ok <- intersect(features, rownames(obj[["RNA"]]))
  if (length(genes_ok) == 0) {
    message("  No RNA features for sample-origin violin pages — skipping")
    return(invisible(NULL))
  }

  lbl <- paste0(obj_nm, "  |  dims 1:", dims_i, "  |  res ", current_res)

  # Paired A/B x-axis ordering (patients interleaved: NPSLE##_A, NPSLE##_B, ...)
  orig_lvls   <- sort(unique(as.character(obj@meta.data[[sample_col]])))
  a_lvls      <- sort(orig_lvls[grepl("_A$", orig_lvls)])
  sample_ids  <- sub("_A$", "", a_lvls)
  paired_lvls <- as.vector(rbind(paste0(sample_ids, "_A"),
                                  paste0(sample_ids, "_B")))
  paired_lvls <- paired_lvls[paired_lvls %in% orig_lvls]
  unpaired    <- orig_lvls[!orig_lvls %in% paired_lvls]
  all_lvls    <- c(paired_lvls, unpaired)

  fill_colors <- setNames(
    ifelse(grepl("_A$", all_lvls), "#CC5500", "#4169E1"),
    all_lvls
  )

  DefaultAssay(obj) <- "RNA"

  clust_order <- as.character(sort(as.numeric(
    unique(as.character(obj@meta.data[[cluster_col]]))
  )))

  gene_chunks <- split(genes_ok, ceiling(seq_along(genes_ok) / markers_per_page))

  for (pg in seq_along(gene_chunks)) {
    chunk  <- gene_chunks[[pg]]
    pg_sfx <- if (length(gene_chunks) > 1)
      paste0("  (pg ", pg, "/", length(gene_chunks), ")") else ""

    vln_panels <- lapply(chunk, function(gene) {
      vln_df <- FetchData(obj, vars = c(gene, cluster_col, sample_col))
      colnames(vln_df) <- c("expr", "cluster", "sample")
      vln_df$cluster   <- factor(vln_df$cluster, levels = clust_order)
      vln_df$sample    <- factor(vln_df$sample,  levels = all_lvls)

      ggplot(vln_df, aes(x = sample, y = expr, fill = sample)) +
        geom_violin(scale = "width", trim = TRUE,
                    linewidth = 0.12, alpha = 0.85) +
        geom_boxplot(width         = 0.10, fill          = "white",
                     outlier.size  = 0.25, outlier.alpha = 0.30,
                     linewidth     = 0.18) +
        facet_wrap(~ cluster, nrow = 1,
                   labeller = as_labeller(function(x) paste0("C", x))) +
        scale_fill_manual(values = fill_colors, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        labs(title = gene, x = NULL, y = "Log norm.") +
        theme_classic(base_size = 7) +
        theme(
          axis.text.x      = element_text(angle = 55, hjust = 1, size = 5.5),
          axis.text.y      = element_text(size = 6),
          axis.title.y     = element_text(size = 7),
          strip.text       = element_text(size = 7.5, face = "bold"),
          strip.background = element_rect(fill = "grey92", color = "grey70",
                                          linewidth = 0.35),
          panel.spacing    = unit(0.2, "lines"),
          plot.title       = element_text(size = 10, face = "bold",
                                          margin = margin(b = 3))
        )
    })

    p_page <- wrap_plots(vln_panels, ncol = 1) +
      plot_annotation(
        title    = paste0(lbl, "  —  RNA by sample_origin", pg_sfx),
        subtitle = "Canonical markers by cluster  |  x = sample_origin  |  orange = _A   blue = _B",
        theme    = theme(
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size =  9, color = "grey40")
        )
      )
    print(p_page)
  }

  invisible(NULL)
}

# =============================================================================
#  ADT PAGES  — prints directly to the open PDF device
# =============================================================================

print_adt_pages <- function(obj, adt_markers, umap_name, cluster_col,
                             hto_col, hto_palette, dims_i, current_res, obj_nm) {
  adt_present <- intersect(adt_markers, rownames(obj[["ADT"]]))
  adt_missing <- setdiff(adt_markers,   rownames(obj[["ADT"]]))
  if (length(adt_missing) > 0)
    message("  ADT absent from object: ", paste(adt_missing, collapse = ", "))
  message("  Printing ", length(adt_present), " ADT pages (res ", current_res, ")...")

  all_rna   <- rownames(obj[["RNA"]])
  cap_base  <- paste0("dims 1:", dims_i, "  |  res=", current_res,
                       "  |  n_cells=", ncol(obj), "  |  CLR margin=2")

  for (adt in adt_present) {
    rna_gene    <- ADT_RNA_MAP[adt]
    has_rna     <- !is.na(rna_gene) && (rna_gene %in% all_rna)
    adt_display <- sub("^Hu[A-Za-z]*\\.", "", adt)

    # Cluster UMAP (no legend)
    p_cluster <- DimPlot(obj, reduction = umap_name, group.by = cluster_col,
                          label = TRUE, label.size = 3, pt.size = 0.25) +
      ggtitle("Clusters") +
      umap_base_theme() +
      theme(legend.position = "none")

    # ADT UMAP
    DefaultAssay(obj) <- "ADT"
    p_adt <- FeaturePlot(obj, features = adt, reduction = umap_name,
                          pt.size = 0.25, order = TRUE) +
      scale_color_gradientn(
        colors = c("lightgrey", "steelblue2", "navy"),
        limits = c(0, 3), oob = squish,
        name   = "CLR\n(m=2)"
      ) +
      ggtitle(paste0("ADT: ", adt_display)) +
      umap_base_theme() +
      theme(legend.position = "right")

    # RNA counterpart UMAP (or blank)
    if (has_rna) {
      DefaultAssay(obj) <- "RNA"
      p_rna <- FeaturePlot(obj, features = rna_gene, reduction = umap_name,
                            pt.size = 0.25, order = TRUE) +
        scale_color_gradientn(
          colors = c("lightgrey", "tomato2", "firebrick4"),
          name   = "Log\nnorm."
        ) +
        ggtitle(paste0("RNA: ", rna_gene)) +
        umap_base_theme() +
        theme(legend.position = "right")
    } else {
      p_rna <- blank_panel()
    }

    # ADT violin — facet by cluster, x = HTO_maxID, fill = HTO_maxID
    DefaultAssay(obj) <- "ADT"
    vln_df <- FetchData(obj, vars = c(adt, cluster_col, hto_col))
    colnames(vln_df) <- c("expr", "cluster", "donor")
    vln_df$cluster   <- factor(
      vln_df$cluster,
      levels = as.character(sort(as.numeric(unique(as.character(vln_df$cluster)))))
    )
    vln_df$donor <- factor(vln_df$donor)

    # Subset palette to donors present (keeps color consistent across pages)
    donors_here <- levels(vln_df$donor)
    fill_vals   <- hto_palette[donors_here]
    fill_vals[is.na(fill_vals)] <- "grey60"  # fallback for any unexpected HTO

    p_vln <- ggplot(vln_df, aes(x = donor, y = expr, fill = donor)) +
      geom_violin(scale = "width", trim = TRUE,
                  linewidth = 0.15, alpha = 0.85) +
      geom_boxplot(width         = 0.12, fill          = "white",
                   outlier.size  = 0.4,  outlier.alpha = 0.35,
                   linewidth     = 0.25) +
      facet_wrap(~ cluster, nrow = 1,
                 labeller = as_labeller(function(x) paste0("C", x))) +
      scale_fill_manual(values = fill_vals, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
      labs(
        title = paste0(
          adt_display,
          if (has_rna) paste0("  (RNA: ", rna_gene, ")") else "  (no RNA counterpart)",
          "    |    CLR expr by cluster × donor"
        ),
        x = NULL,
        y = "CLR norm. expression (margin=2)"
      ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x      = element_text(angle = 50, hjust = 1, size = 6),
        axis.text.y      = element_text(size = 7),
        axis.title.y     = element_text(size = 8),
        strip.text       = element_text(size = 8.5, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = "grey70",
                                        linewidth = 0.4),
        panel.spacing    = unit(0.25, "lines"),
        plot.title       = element_text(size = 9, face = "bold",
                                        margin = margin(b = 6))
      )

    top_row   <- (p_cluster | p_adt | p_rna) + plot_layout(widths = c(1, 1, 1))
    full_page <- (top_row / p_vln) +
      plot_layout(heights = c(4.2, 6.8)) +
      plot_annotation(
        title   = paste0(obj_nm, "  |  dims 1:", dims_i, "  |  ADT: ", adt_display),
        caption = cap_base,
        theme   = theme(
          plot.title   = element_text(size = 9, color = "grey30", face = "bold"),
          plot.caption = element_text(size = 7, color = "grey55")
        )
      )

    print(full_page)
    DefaultAssay(obj) <- "RNA"
  }
  invisible(NULL)
}

# =============================================================================
#  RESOLUTION GRID PAGE  — one page per dims_i; 2×2 grid at fixed resolutions
# =============================================================================

print_resolution_grid_page <- function(obj_d, snn_name, umap_name, dims_i, obj_nm,
                                        grid_res = c(0.2, 0.4, 0.6, 0.8)) {
  lbl   <- paste0(obj_nm, "  |  dims 1:", dims_i)
  d_tag <- paste0("d", dims_i)

  p_list <- vector("list", length(grid_res))
  for (i in seq_along(grid_res)) {
    r         <- grid_res[i]
    res_tag   <- paste0("r", gsub("\\.", "p", sprintf("%.2f", r)))
    clust_col <- paste0("clust_", d_tag, "_", res_tag)

    if (!clust_col %in% colnames(obj_d@meta.data)) {
      message("    res grid: FindClusters at res ", r, " ...")
      obj_d <- FindClusters(obj_d, graph.name  = snn_name,
                             resolution   = r,
                             cluster.name = clust_col,
                             verbose      = FALSE)
    }

    n_cl <- length(unique(obj_d@meta.data[[clust_col]]))
    p_list[[i]] <- DimPlot(obj_d, reduction = umap_name, group.by = clust_col,
                            label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4) +
      ggtitle(paste0("res ", r, "  (", n_cl, " clusters)")) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none",
            plot.title      = element_text(size = 11, face = "bold"),
            aspect.ratio    = 1)
  }

  p_grid <- wrap_plots(p_list, ncol = 2) +
    plot_annotation(
      title    = paste0(lbl, "  —  Resolution Grid"),
      subtitle = paste0("UMAP at res ", paste(grid_res, collapse = ", ")),
      theme    = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )
  print(p_grid)
  invisible(NULL)
}

# =============================================================================
#  MAIN WORKER
# =============================================================================

run_survey <- function(rds_list, outdir, dims_vec, resolutions, features,
                        hto_col = "HTO_maxID") {
  stopifnot(is.list(rds_list), !is.null(names(rds_list)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  for (obj_nm in names(rds_list)) {
    message("\n════ Object: ", obj_nm, " ════")
    obj <- rds_list[[obj_nm]]

    # Validate PCA
    if (!"pca" %in% names(obj@reductions))
      stop("Object '", obj_nm, "' has no PCA. Run RunPCA() first.")
    max_pcs  <- ncol(obj@reductions$pca)
    dims_use <- dims_vec[dims_vec <= max_pcs]
    if (length(dims_vec) != length(dims_use))
      warning("Dropped dims exceeding available PCs (", max_pcs, "): ",
              paste(dims_vec[dims_vec > max_pcs], collapse = ", "))

    # Validate HTO column
    if (!hto_col %in% colnames(obj@meta.data))
      stop("Column '", hto_col, "' not found in metadata.")

    # ADT CLR margin=2 — done once per object, shared across all dims/res
    obj <- normalize_adt_clr2(obj)

    # Build HTO color palette — Set2 fits 8 donors exactly
    hto_levels  <- sort(unique(as.character(obj@meta.data[[hto_col]])))
    n_hto       <- length(hto_levels)
    pal_nm      <- if (n_hto <= 8) "Set2" else "Set3"
    pal_n       <- min(n_hto, brewer.pal.info[pal_nm, "maxcolors"])
    hto_palette <- setNames(
      colorRampPalette(brewer.pal(max(3L, pal_n), pal_nm))(n_hto),
      hto_levels
    )
    message("  HTO palette: ", pal_nm, " (", n_hto, " donors)")

    obj_outdir <- file.path(outdir, obj_nm)
    dir.create(obj_outdir, showWarnings = FALSE, recursive = TRUE)

    for (dims_i in dims_use) {
      dims_seq  <- seq_len(dims_i)
      d_tag     <- paste0("d", dims_i)
      snn_name  <- paste0("snn_",  d_tag)
      umap_name <- paste0("umap_", d_tag)
      umap_key  <- paste0("umap",  d_tag, "_")

      message("\n  ── dims 1:", dims_i, " ──────────")
      message("  FindNeighbors...")
      obj_d <- FindNeighbors(obj, reduction = "pca", dims = dims_seq,
                              graph.name = snn_name, verbose = FALSE)

      # Clustree metadata (runs FindClusters 20× at 0.05 steps)
      ct_result <- build_clustree_metadata(obj_d, snn_name, dims_i)
      obj_d     <- ct_result$obj  # obj_d now carries ct_d{N}_res.* columns
      ct_info   <- ct_result      # ct_info$obj is a copy: used only for clustree

      message("  RunUMAP (dims 1:", dims_i, ")...")
      obj_d <- RunUMAP(obj_d, reduction = "pca", dims = dims_seq,
                        reduction.name = umap_name, reduction.key = umap_key,
                        seed.use = 42, verbose = FALSE)

      pdf_path <- file.path(obj_outdir,
                             paste0("TcellSurvey_", obj_nm, "_d", dims_i, ".pdf"))
      pdf(pdf_path, width = 17, height = 11)
      message("  PDF open: ", pdf_path)

      for (res_i in resolutions) {
        res_tag   <- paste0("r", gsub("\\.", "p", sprintf("%.2f", res_i)))
        clust_col <- paste0("clust_", d_tag, "_", res_tag)

        message("    res = ", res_i, " ...")
        obj_d <- FindClusters(obj_d, graph.name  = snn_name,
                               resolution  = res_i,
                               cluster.name = clust_col,
                               verbose     = FALSE)
        Idents(obj_d) <- clust_col
        n_cl <- length(unique(obj_d@meta.data[[clust_col]]))
        message("    -> ", n_cl, " clusters")

        # RNA QC pages 1, 2, 3+
        print_rna_qc_pages(
          obj         = obj_d,
          ct_info     = ct_info,
          umap_name   = umap_name,
          cluster_col = clust_col,
          features    = features,
          current_res = res_i,
          dims_i      = dims_i,
          obj_nm      = obj_nm,
          hto_col     = hto_col,
          hto_palette = hto_palette
        )

        # RNA by sample_origin pages — appended after other RNA plots
        print_sample_origin_violin_pages(
          obj         = obj_d,
          cluster_col = clust_col,
          features    = features,
          dims_i      = dims_i,
          current_res = res_i,
          obj_nm      = obj_nm
        )

        # ADT pages — TCELL_ADT_MARKERS order, T-cell markers first
        print_adt_pages(
          obj         = obj_d,
          adt_markers = TCELL_ADT_MARKERS,
          umap_name   = umap_name,
          cluster_col = clust_col,
          hto_col     = hto_col,
          hto_palette = hto_palette,
          dims_i      = dims_i,
          current_res = res_i,
          obj_nm      = obj_nm
        )
      }  # end resolution loop

      # Resolution grid page — once per dims_i, fixed resolutions 0.2/0.4/0.6/0.8
      message("  Resolution grid page (dims 1:", dims_i, ")...")
      print_resolution_grid_page(
        obj_d    = obj_d,
        snn_name = snn_name,
        umap_name = umap_name,
        dims_i   = dims_i,
        obj_nm   = obj_nm
      )

      dev.off()
      message("  PDF saved: ", pdf_path)
    }  # end dims loop
  }  # end object loop

  invisible(NULL)
}

# =============================================================================
#  ENTRY POINT
# =============================================================================

main <- function() {
  cat("Loading RDS:", opt$rds, "\n")
  rds_list <- load_rds_input(opt$rds)

  run_survey(
    rds_list    = rds_list,
    outdir      = opt$outdir,
    dims_vec    = dims_vec,
    resolutions = resolutions,
    features    = features,
    hto_col     = opt$hto_col
  )

  cat("\n====== Done ======\n")
}

main()
