## ============================================================================
##  NPSLE ADT/RNA UMAP Survey Pipeline
##  For dims 1:12 and 1:16 (re-clusters on existing PCA):
##    – Landscape PDF per dim-set, one page per ADT
##    – Page layout: [cluster UMAP | ADT UMAP | RNA UMAP*] / [violin by cluster×sample]
##    * RNA panel shown when a canonical gene counterpart exists; blank otherwise
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(optparse)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
# Usage:
#   Rscript NPSLE_ADT_UMAP_Survey_Pipeline.R \
#     --input_rds  /path/to/object.rds \
#     --out_dir    /path/to/output \
#     --resolution 0.3 \
#     --dims       "1:12,1:16,1:20"
#
# --dims accepts a comma-separated list of colon-range strings, e.g. "1:12,1:16".
# Each range becomes its own dim-set label (dims1_12, dims1_16, …) and PDF.
option_list <- list(
  make_option("--input_rds",   type = "character", default = NULL,
              help = "Path to Seurat RDS file"),
  make_option("--out_dir",     type = "character", default = NULL,
              help = "Output directory"),
  make_option("--resolution",  type = "double",    default = NULL,
              help = "Clustering resolution (default: 0.3)"),
  make_option("--dims",        type = "character", default = NULL,
              help = "Comma-separated dim ranges, e.g. \"1:12,1:16,1:20\" (default: 1:12,1:16)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Parse --dims string into a named list of integer sequences
parse_dims_arg <- function(dims_str) {
  ranges <- trimws(strsplit(dims_str, ",")[[1]])
  out <- lapply(ranges, function(r) {
    parts <- as.integer(strsplit(r, ":")[[1]])
    if (length(parts) != 2 || any(is.na(parts)))
      stop("--dims: cannot parse range '", r, "'. Use colon notation, e.g. 1:12")
    seq(parts[1], parts[2])
  })
  names(out) <- paste0("dims", gsub(":", "_", ranges))
  out
}

# ── 0. Configuration ──────────────────────────────────────────────────────────
cfg <- list(
  input_rds  = file.path(
    "/path/to/data",
    "20260213_GroupComps_NPSLEsamples",
    "NPSLE_CellTypeCog_Hash_dims8_rez3.rds"
  ),
  out_dir    = file.path(
    "/path/to/data",
    "20260213_GroupComps_NPSLEsamples",
    "ADT_UMAP_Survey_Output"
  ),
  resolution = 0.3,
  dims_list  = list(
    dims1_12 = 1:12,
    dims1_16 = 1:16
  )
)

# CLI overrides
if (!is.null(opt$input_rds))  cfg$input_rds  <- opt$input_rds
if (!is.null(opt$out_dir))    cfg$out_dir    <- opt$out_dir
if (!is.null(opt$resolution)) cfg$resolution <- opt$resolution
if (!is.null(opt$dims))       cfg$dims_list  <- parse_dims_arg(opt$dims)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. ADT → canonical RNA gene map ──────────────────────────────────────────
adt_rna_map <- c(
  # T cell lineage / activation
  "Hu.CD3-UCHT1"       = "CD3E",
  "Hu.CD4-RPA.T4"      = "CD4",
  "Hu.CD8"             = "CD8A",
  "Hu.CD2"             = "CD2",
  "Hu.CD5"             = "CD5",
  "Hu.CD7"             = "CD7",
  "Hu.CD28"            = "CD28",
  "Hu.CD27"            = "CD27",
  "Hu.CD25"            = "IL2RA",
  "Hu.CD127"           = "IL7R",
  "Hu.CD69"            = "CD69",
  "Hu.CD44"            = "CD44",     # mapped via HuMs.CD44 below
  "HuMs.CD44"          = "CD44",
  # T cell subsets / memory
  "Hu.CD45RA"          = "PTPRC",
  "Hu.CD45RO"          = "PTPRC",
  "Hu.CD45-HI30"       = "PTPRC",
  "Hu.CD62L"           = "SELL",
  "Hu.CD161"           = "KLRB1",
  "Hu.CD57"            = "B3GAT1",
  "Hu.KLRG1"           = "KLRG1",
  # TCR
  "Hu.TCR.AB"          = "TRAC",
  "Hu.TCR.Va7.2"       = "TRAV1-2",
  "Hu.TCR.Vd2"         = "TRDV2",
  # Checkpoints / co-inhibitory
  "Hu.CD279"           = "PDCD1",
  "Hu.CD274"           = "CD274",
  "Hu.CD223"           = "LAG3",
  "Hu.TIGIT"           = "TIGIT",
  "Hu.CD152"           = "CTLA4",
  "Hu.CD272"           = "BTLA",
  "Hu.CD270"           = "TNFRSF14",
  "Hu.CD85j"           = "LILRB1",
  # Co-stimulatory
  "Hu.CD137"           = "TNFRSF9",
  "Hu.CD134"           = "TNFRSF4",
  "HuMsRt.CD278"       = "ICOS",
  "Hu.CD226-11A8"      = "CD226",
  "Hu.CD244"           = "CD244",
  "Hu.CD155"           = "PVR",
  "Hu.CD112"           = "NECTIN2",
  "Hu.CD352"           = "SLAMF6",
  "Hu.CD319"           = "SLAMF7",
  # NK markers
  "Hu.CD56"            = "NCAM1",
  "Hu.CD314"           = "KLRK1",
  "Hu.CD94"            = "KLRD1",
  "Hu.CD328"           = "SIGLEC7",
  "Hu.CD158"           = "KIR2DL1",
  "Hu.CD158b"          = "KIR2DL2",
  "Hu.CD158e1"         = "KIR3DL1",
  # Monocyte / macrophage
  "Hu.CD14-M5E2"       = "CD14",
  "Hu.CD16"            = "FCGR3A",
  "Hu.CD163"           = "CD163",
  "Hu.CD64"            = "FCGR1A",
  "Hu.CD32"            = "FCGR2A",
  "Hu.CD33"            = "CD33",
  "Hu.CD169"           = "SIGLEC1",
  "Hu.CD204-MSR1"      = "MSR1",
  "Hu.CLEC12A"         = "CLEC12A",
  "Hu.CX3CR1"          = "CX3CR1",
  "Hu.LOX.1"           = "OLR1",
  "Hu.GPR56"           = "ADGRG1",
  # DC
  "Hu.CD11c"           = "ITGAX",
  "Hu.CD1c"            = "CD1C",
  "Hu.CD1d"            = "CD1D",
  "Hu.CD141"           = "THBD",
  "Hu.CD303"           = "CLEC4C",
  "Hu.HLA.DR"          = "HLA-DRA",
  "Hu.HLA.ABC"         = "HLA-A",
  "Hu.HLA.E"           = "HLA-E",
  # B cell
  "Hu.CD19"            = "CD19",
  "Hu.CD20-2H7"        = "MS4A1",
  "Hu.CD22"            = "CD22",
  "Hu.CD21"            = "CR2",
  "Hu.CD24"            = "CD24",
  "Hu.CD38-HIT2"       = "CD38",
  "Hu.CD79b"           = "CD79B",
  "Hu.CD23"            = "FCER2",
  "Hu.CD40"            = "CD40",
  "Hu.CD83"            = "CD83",
  "Hu.CD86"            = "CD86",
  "Hu.CD267"           = "TNFRSF13B",
  "Hu.CD268"           = "TNFRSF13C",
  "Hu.IgM"             = "IGHM",
  "Hu.IgD"             = "IGHD",
  "Hu.Ig.LightChain.k" = "IGKC",
  "Hu.Ig.LightChain.l" = "IGLC2",
  # Chemokine receptors
  "Hu.CD183"           = "CXCR3",
  "Hu.CD185"           = "CXCR5",
  "Hu.CD194"           = "CCR4",
  "Hu.CD195"           = "CCR5",
  "Hu.CD196"           = "CCR6",
  # Integrins / adhesion
  "Hu.CD11a"           = "ITGAL",
  "Hu.CD11b"           = "ITGAM",
  "Hu.CD18"            = "ITGB2",
  "Hu.CD29"            = "ITGB1",
  "Hu.CD49a"           = "ITGA1",
  "Hu.CD49b"           = "ITGA2",
  "Hu.CD49d"           = "ITGA4",
  "HuMs.CD49f"         = "ITGA6",
  "HuMs.integrin.b7"   = "ITGB7",
  "Hu.CD103"           = "ITGAE",
  "Hu.CD54"            = "ICAM1",
  "Hu.CD31"            = "PECAM1",
  "Hu.CD146"           = "MCAM",
  "Hu.CD62P"           = "SELP",
  # Other surface molecules
  "Hu.CD101"           = "CD101",
  "Hu.CD105-43A3"      = "ENG",
  "Hu.CD107a"          = "LAMP1",
  "Hu.CD119"           = "IFNGR1",
  "Hu.CD122"           = "IL2RB",
  "Hu.CD123"           = "IL3RA",
  "Hu.CD124"           = "IL4R",
  "Hu.CD13"            = "ANPEP",
  "Hu.CD154"           = "CD40LG",
  "Hu.CD224"           = "GGT1",
  "Hu.CD26"            = "DPP4",
  "Hu.CD35"            = "CR1",
  "Hu.CD36"            = "CD36",
  "Hu.CD39"            = "ENTPD1",
  "Hu.CD41"            = "ITGA2B",
  "Hu.CD42b"           = "GP1BA",
  "Hu.CD47"            = "CD47",
  "Hu.CD48"            = "CD48",
  "Hu.CD52"            = "CD52",
  "Hu.CD58"            = "CD58",
  "Hu.CD71"            = "TFRC",
  "Hu.CD73"            = "NT5E",
  "Hu.CD81"            = "CD81",
  "Hu.CD82"            = "CD82",
  "Hu.CD88"            = "C5AR1",
  "Hu.CD95"            = "FAS",
  "Hu.CD99"            = "CD99",
  "Hu.FceRIa"          = "FCER1A"
)

# ── 2. Load object ────────────────────────────────────────────────────────────
cat("Loading:", cfg$input_rds, "\n")
npsle        <- readRDS(cfg$input_rds)
all_adts     <- npsle@assays[["ADT"]]@counts@Dimnames[[1]]
all_rna_feat <- rownames(npsle[["RNA"]])

# ── 3. Helper functions ───────────────────────────────────────────────────────
blank_panel <- function(msg = "No RNA\ncounterpart") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 5.5, color = "grey55", hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(aspect.ratio    = 1,
          plot.background = element_rect(fill = "grey97", color = "grey85", linewidth = 0.3))
}

umap_base_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      aspect.ratio  = 1,
      axis.line     = element_blank(),
      axis.ticks    = element_blank(),
      axis.text     = element_blank(),
      axis.title    = element_text(size = base_size - 1, color = "grey40"),
      plot.title    = element_text(size = base_size, face = "bold", margin = margin(b = 4)),
      legend.text   = element_text(size = base_size - 2),
      legend.title  = element_text(size = base_size - 2),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.35, "cm")
    )
}

# ── 4. Main loop over dim sets ────────────────────────────────────────────────
for (dim_label in names(cfg$dims_list)) {
  dims        <- cfg$dims_list[[dim_label]]
  cluster_col <- paste0("clusters_", dim_label)   # e.g. clusters_dims1_12
  umap_name   <- paste0("umap_",     dim_label)   # e.g. umap_dims1_12
  umap_key    <- gsub("[^A-Za-z0-9]", "", dim_label)  # stripped key prefix

  cat("\n========================================\n")
  cat(" Processing:", dim_label, "(dims", min(dims), "–", max(dims), ")\n")
  cat("========================================\n")

  obj <- npsle  # fresh working copy per dim-set

  # Recluster using existing PCA ── do NOT re-run PCA
  cat("  FindNeighbors ... ")
  obj <- FindNeighbors(obj,
                       reduction  = "pca",
                       dims       = dims,
                       graph.name = c("RNA_nn_sv", "RNA_snn_sv"),
                       verbose    = FALSE)
  cat("done\n")

  cat("  FindClusters (res =", cfg$resolution, ") ... ")
  obj <- FindClusters(obj,
                      graph.name   = "RNA_snn_sv",
                      resolution   = cfg$resolution,
                      cluster.name = cluster_col,
                      verbose      = FALSE)
  cat("done\n")

  cat("  RunUMAP ... ")
  obj <- RunUMAP(obj,
                 reduction       = "pca",
                 dims            = dims,
                 reduction.name  = umap_name,
                 reduction.key   = paste0(umap_key, "_"),
                 seed.use        = 42,
                 verbose         = FALSE)
  cat("done\n")

  n_cl <- nlevels(obj[[cluster_col, drop = TRUE]])
  cat("  Clusters detected:", n_cl, "\n")
  cat("  ADTs to process:", length(all_adts), "\n")

  # ── Open PDF ─────────────────────────────────────────────────────────────
  pdf_path <- file.path(cfg$out_dir,
                        paste0("NPSLE_ADT_Survey_", dim_label, ".pdf"))
  pdf(pdf_path, width = 17, height = 11)
  cat("  Writing PDF:", pdf_path, "\n")

  # Page 0: reference cluster UMAP (full legend, no ADT)
  p_ref <- DimPlot(obj,
                   reduction   = umap_name,
                   group.by    = cluster_col,
                   label       = TRUE,
                   label.size  = 4,
                   repel       = TRUE,
                   pt.size     = 0.5) +
    ggtitle(paste0("Seurat Clusters — ", dim_label,
                   "  (res = ", cfg$resolution, ",  n = ", n_cl, ")")) +
    umap_base_theme(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(size = 12, face = "bold"))
  print(p_ref)

  # ── Per-ADT pages ─────────────────────────────────────────────────────────
  for (i in seq_along(all_adts)) {
    adt <- all_adts[i]
    cat(sprintf("  [%d/%d] %s\n", i, length(all_adts), adt))

    rna_gene <- adt_rna_map[adt]
    has_rna  <- !is.na(rna_gene) && (rna_gene %in% all_rna_feat)

    # ── Top-left: cluster UMAP (no legend, small labels) ──────────────────
    p_cluster <- DimPlot(obj,
                         reduction  = umap_name,
                         group.by   = cluster_col,
                         label      = TRUE,
                         label.size = 3,
                         pt.size    = 0.25) +
      ggtitle("Clusters") +
      umap_base_theme() +
      theme(legend.position = "none")

    # ── Top-middle: ADT expression UMAP ───────────────────────────────────
    DefaultAssay(obj) <- "ADT"
    p_adt <- FeaturePlot(obj,
                         features  = adt,
                         reduction = umap_name,
                         pt.size   = 0.25,
                         order     = TRUE) +
      scale_color_gradientn(
        colors = c("lightgrey", "steelblue2", "navy"),
        limits = c(0, 3),
        oob    = squish,
        name   = "CLR\nnorm."
      ) +
      ggtitle(paste0("ADT: ", adt)) +
      umap_base_theme() +
      theme(legend.position = "right")

    # ── Top-right: RNA expression UMAP (or blank placeholder) ─────────────
    if (has_rna) {
      DefaultAssay(obj) <- "RNA"
      p_rna <- FeaturePlot(obj,
                           features  = rna_gene,
                           reduction = umap_name,
                           pt.size   = 0.25,
                           order     = TRUE) +
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

    # ── Bottom: violin per cluster × sample_origin ────────────────────────
    DefaultAssay(obj) <- "ADT"
    vln_df <- FetchData(obj, vars = c(adt, cluster_col, "sample_origin"))
    colnames(vln_df)[1] <- "expr"
    vln_df$cluster      <- factor(vln_df[[cluster_col]])

    # Reorder sample_origin so each patient's A and B sit adjacent:
    # take the A-group order as reference, interleave B next to each A
    orig_lvls  <- levels(vln_df$sample_origin)
    a_lvls     <- orig_lvls[grepl("_A$", orig_lvls)]          # preserves original A ordering
    sample_ids <- sub("_A$", "", a_lvls)                       # e.g. "NPSLE27"
    paired_lvls <- as.vector(rbind(paste0(sample_ids, "_A"),   # A then B for each patient
                                   paste0(sample_ids, "_B")))
    paired_lvls <- paired_lvls[paired_lvls %in% orig_lvls]     # drop any missing batch
    vln_df$sample_origin <- factor(vln_df$sample_origin, levels = paired_lvls)

    # Burnt orange for _A, royal blue for _B
    fill_colors <- setNames(
      ifelse(grepl("_A$", paired_lvls), "#CC5500", "#4169E1"),
      paired_lvls
    )

    p_vln <- ggplot(vln_df,
                    aes(x = sample_origin, y = expr, fill = sample_origin)) +
      geom_violin(scale     = "width",
                  trim      = TRUE,
                  linewidth = 0.15,
                  alpha     = 0.85) +
      geom_boxplot(width         = 0.12,
                   fill          = "white",
                   outlier.size  = 0.4,
                   outlier.alpha = 0.35,
                   linewidth     = 0.25) +
      facet_wrap(~ cluster, nrow = 1,
                 labeller = as_labeller(function(x) paste0("C", x))) +
      scale_fill_manual(values = fill_colors, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
      labs(
        title = paste0(adt,
                       if (has_rna) paste0("  (", rna_gene, ")") else "  (no RNA counterpart)",
                       "    |    expression by cluster × sample"),
        x = NULL,
        y = "CLR norm. expression"
      ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x      = element_text(angle  = 50, hjust = 1, size = 5.5),
        axis.text.y      = element_text(size   = 7),
        axis.title.y     = element_text(size   = 8),
        strip.text       = element_text(size   = 8.5, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = "grey70",
                                        linewidth = 0.4),
        panel.spacing    = unit(0.25, "lines"),
        plot.title       = element_text(size  = 9, face = "bold",
                                        margin = margin(b = 6))
      )

    # ── Assemble page ──────────────────────────────────────────────────────
    top_row   <- (p_cluster | p_adt | p_rna) +
      plot_layout(widths = c(1, 1, 1))

    full_page <- top_row / p_vln +
      plot_layout(heights = c(4.2, 6.8)) +
      plot_annotation(
        title   = paste0(dim_label, "  |  ADT: ", adt),
        caption = paste0("dims ", min(dims), ":", max(dims),
                         "  |  res=", cfg$resolution,
                         "  |  n_cells=", ncol(obj)),
        theme   = theme(
          plot.title   = element_text(size = 8, color = "grey40"),
          plot.caption = element_text(size = 7, color = "grey55")
        )
      )

    print(full_page)
  }

  dev.off()
  cat("  PDF saved:", pdf_path, "\n")
}

cat("\n====== Done ======\n")
