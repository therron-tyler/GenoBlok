#!/usr/bin/env Rscript
# BulkSeq_QC.R — Bulk RNA-seq quality-control visualizations
#
# Usage:
#   Rscript BulkSeq_QC.R <cpm_file> <raw_counts_file> <groupfile> <output_dir>
#
#   cpm_file       : tab-delimited; col1 = Symbol (ENSG IDs, becomes row.names),
#                    col2 = Gene_Symbol, remaining cols = sample CPM values
#   raw_counts_file: same format as cpm_file but integer raw counts.
#                    Pass "none" to skip raw-count panel.
#   groupfile      : CSV with header; col1 = sample_name, remaining cols = metadata
#   output_dir     : output directory (created if absent)
#
# Outputs (all in output_dir):
#   QC_01_filter_grid.png                — 4-panel: CPM dist / sensitivity / legend / raw counts
#   QC_02_hist_all_samples.png           — all-sample overlay, IQR outliers labeled
#   QC_03_hist_by_<col>.png              — per metadata column, lines colored by group
#   QC_04_hist_avg_per_<col1>.png        — per-group averaged lines (first metadata col)
#   QC_05_hist_zoom_outliers_bold.png    — y-clipped to 4000, outlier lines bolded
#   QC_06_correlation_heatmap.png        — square sample×sample Pearson heatmap
#   QC_06_correlation_gene_drivers.csv
#   QC_06_correlation_sample_mean.csv
#   QC_07_PCA_by_<col>.png              — one per metadata column, outliers labeled large
#   QC_07_PCA_coordinates.csv
#   QC_02_hist_sample_bin_data.csv      — per-sample histogram bin counts (log2(CPM+1))

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggrepel)
  library(cowplot)
  library(matrixStats)
  library(scales)
})

set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# 0. ARGUMENTS & QC PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: Rscript BulkSeq_QC.R <cpm_file> <raw_counts_file> <groupfile> <output_dir>",
    "  Pass 'none' as raw_counts_file to skip the raw-count panel.", sep = "\n"
  ))
}
CPM_FILE   <- args[1]
RAW_FILE   <- args[2]
GROUP_FILE <- args[3]
OUT_DIR    <- args[4]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# CPM thresholds shown as vertical dashed lines on all histogram panels
CUTOFF_LIST <- c(5, 10, 20)
# % sample thresholds used in the sensitivity curve (Panel 2 of the 4-panel grid)
PCT_LIST    <- c(10, 20, 25)
# Histogram bin width in log2(CPM+1) space
BIN_WIDTH   <- 0.5
# Max x-axis extent for histograms
MAX_LOG2    <- 16

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
message("=== [1] Loading data ===")

# CPM matrix ─ col1 becomes row.names, col2 = Gene_Symbol, rest = samples
cpm_raw <- read.delim(CPM_FILE, row.names = 1, check.names = FALSE)
gene_sym_map <- data.frame(
  ENSG        = rownames(cpm_raw),
  Gene_Symbol = cpm_raw[["Gene_Symbol"]],
  stringsAsFactors = FALSE
)
cpm <- cpm_raw[, -1, drop = FALSE]   # drop Gene_Symbol column

# Raw counts (optional)
use_raw <- !is.na(RAW_FILE) && RAW_FILE != "none" && file.exists(RAW_FILE)
if (use_raw) {
  raw_raw <- read.delim(RAW_FILE, row.names = 1, check.names = FALSE)
  raw     <- raw_raw[, -1, drop = FALSE]
  message("  Raw counts: ", nrow(raw), " genes x ", ncol(raw), " samples")
} else {
  raw <- NULL
  message("  Raw counts: not provided — panel 4 will be blank")
}

# Group file ─ header row, col1 = sample_name, rest = metadata
gf <- read.csv(GROUP_FILE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(gf)[1] <- "sample_name"
meta_cols        <- setdiff(colnames(gf), "sample_name")

message("  CPM matrix: ", nrow(cpm), " genes x ", ncol(cpm), " samples")
message("  Group file: ", nrow(gf), " samples | metadata: ", paste(meta_cols, collapse = ", "))

# Restrict to samples present in both CPM and group file, in groupfile order
shared_samples <- intersect(gf$sample_name, colnames(cpm))
if (length(shared_samples) == 0)
  stop("No sample names match between CPM file and group file. Check column headers.")

cpm <- cpm[, shared_samples, drop = FALSE]
if (use_raw) {
  raw_shared <- raw[, intersect(shared_samples, colnames(raw)), drop = FALSE]
} else {
  raw_shared <- NULL
}

gf <- gf[match(shared_samples, gf$sample_name), ]
rownames(gf) <- gf$sample_name

message("  Shared samples (CPM + group file): ", length(shared_samples))

# ─────────────────────────────────────────────────────────────────────────────
# 2. AUTO-ASSIGN COLORS AND SHAPES PER METADATA COLUMN
# ─────────────────────────────────────────────────────────────────────────────
message("=== [2] Auto-assigning colors and shapes ===")

META_PALETTES <- list(
  c("#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3","#A65628","#F781BF","#999999"),
  c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
  c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
  c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
)
# Fillable shapes that distinguish well on screen and in print
AUTO_SHAPES <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

col_colors <- list()
col_shapes <- list()
for (i in seq_along(meta_cols)) {
  mc   <- meta_cols[i]
  vals <- sort(unique(as.character(gf[[mc]][!is.na(gf[[mc]])])))
  n    <- length(vals)
  pal  <- META_PALETTES[[(i - 1L) %% length(META_PALETTES) + 1L]]
  col_colors[[mc]] <- setNames(colorRampPalette(pal)(n), vals)
  col_shapes[[mc]] <- setNames(AUTO_SHAPES[((seq_len(n) - 1L) %% length(AUTO_SHAPES)) + 1L], vals)
}

# Per-sample color palette for individual-line histograms
n_samp  <- length(shared_samples)
samp_pal <- setNames(
  colorRampPalette(c("#4575b4","#74add1","#abd9e9","#fdae61","#f46d43","#d73027",
                     "#a50026","#313695","#fee090","#e0f3f8"))(n_samp),
  shared_samples
)

# ─────────────────────────────────────────────────────────────────────────────
# SHARED HELPERS
# ─────────────────────────────────────────────────────────────────────────────

# Build per-sample frequency-polygon data frame from a named list of value vectors.
# group_vec (optional): named character vector mapping sample IDs to a group label.
.make_freq_poly <- function(val_list, breaks, group_vec = NULL) {
  do.call(rbind, lapply(names(val_list), function(sid) {
    vals <- val_list[[sid]]
    vals <- vals[vals > 0 & vals <= max(breaks)]
    h    <- hist(vals, breaks = breaks, plot = FALSE, right = TRUE)
    df   <- data.frame(Sample = sid, bin_mid = h$mids, count = h$counts,
                       stringsAsFactors = FALSE)
    if (!is.null(group_vec)) df$Group <- group_vec[sid]
    df
  }))
}

# Build per-group averaged frequency-polygon data frame.
.make_freq_poly_avg <- function(val_list, group_vec, breaks) {
  grps <- unique(group_vec)
  do.call(rbind, lapply(grps, function(g) {
    sids <- names(group_vec)[group_vec == g]
    sids <- sids[sids %in% names(val_list)]
    if (length(sids) == 0) return(NULL)
    # Average log2(CPM+1) per gene across samples in this group
    mat  <- do.call(cbind, val_list[sids])
    avg  <- if (ncol(mat) > 1) rowMeans(mat) else mat[, 1]
    avg  <- avg[avg > 0 & avg <= max(breaks)]
    h    <- hist(avg, breaks = breaks, plot = FALSE, right = TRUE)
    data.frame(Group = g, bin_mid = h$mids, count = h$counts, stringsAsFactors = FALSE)
  }))
}

# Shared x-axis for log2(CPM+1) histograms
breaks_h    <- seq(0, MAX_LOG2, by = BIN_WIDTH)
.log2_breaks <- seq(BIN_WIDTH / 2, MAX_LOG2 - BIN_WIDTH / 2, by = BIN_WIDTH)
.log2_labels <- paste0("(", seq(0, MAX_LOG2 - BIN_WIDTH, by = BIN_WIDTH),
                        ",", seq(BIN_WIDTH, MAX_LOG2, by = BIN_WIDTH), "]")

# Colour palette for CUTOFF_LIST vlines (shared across all histogram panels)
cutoff_colors <- setNames(
  colorRampPalette(c("#d7191c","#fdae61","#1a9641","#2c7bb6","#7B2D8B"))(length(CUTOFF_LIST)),
  as.character(CUTOFF_LIST)
)
# Reusable geom layers for CUTOFF vlines
cpm_vlines <- lapply(CUTOFF_LIST, function(t) {
  geom_vline(xintercept = log2(t + 1),
             color = cutoff_colors[as.character(t)],
             linetype = "dashed", linewidth = 0.85)
})

# Build CPM log2(CPM+1) list for all shared samples (used in every histogram section)
cpm_mat  <- as.matrix(cpm)
cpm_list <- setNames(lapply(shared_samples, function(s) log2(cpm_mat[, s] + 1)),
                     shared_samples)

# IQR outlier detection on per-sample median log2(CPM+1) — shared across sections
med_all <- sapply(shared_samples, function(s) median(cpm_list[[s]]))
iqr_all <- IQR(med_all)
q1_all  <- quantile(med_all, 0.25)
q3_all  <- quantile(med_all, 0.75)
iqr_outliers <- shared_samples[
  med_all < q1_all - 1.5 * iqr_all | med_all > q3_all + 1.5 * iqr_all
]
if (length(iqr_outliers) > 0)
  message("  IQR outliers: ", paste(iqr_outliers, collapse = ", "))

# ─────────────────────────────────────────────────────────────────────────────
# 3b. 4-PANEL QC HISTOGRAM GRID
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3b] 4-panel QC histogram grid ===")

poly_all <- .make_freq_poly(cpm_list, breaks_h)

# ── Panel 1: CPM distribution — all samples, per-sample palette ──────────────
cutoff_annot_p1 <- lapply(seq_along(CUTOFF_LIST), function(i) {
  t      <- CUTOFF_LIST[i]
  min_s  <- ceiling(0.20 * n_samp)
  n_pass <- sum(rowSums(cpm_mat > t) >= min_s)
  annotate("text",
           x = log2(t + 1) + 0.10, y = Inf,
           vjust = 1.8 + (i - 1) * 3.8,
           label = paste0("t=", t, " CPM: ", n_pass, " genes (≥20%)"),
           color = cutoff_colors[[as.character(t)]], size = 2.4, hjust = 0)
})

p_p1 <- ggplot(poly_all, aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.55, linewidth = 0.40) +
  cpm_vlines + cutoff_annot_p1 +
  scale_color_manual(values = samp_pal, guide = "none") +
  scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                     guide = guide_axis(angle = 90)) +
  labs(title    = paste0("CPM distribution — all ", n_samp, " samples"),
       subtitle = "log2(CPM+1) | bin 0.5 | zeros excluded",
       x = "log2(CPM+1) bin", y = "Gene count") +
  theme_cowplot(12)

# ── Panel 2: Sensitivity curve ────────────────────────────────────────────────
cpm_sweep <- seq(0, 30, by = 0.25)
sens_df <- do.call(rbind, lapply(PCT_LIST, function(pct) {
  min_s   <- ceiling(pct / 100 * n_samp)
  n_genes <- sapply(cpm_sweep, function(t) sum(rowSums(cpm_mat > t) >= min_s))
  data.frame(log2_cpm  = log2(cpm_sweep + 1),
             n_genes   = n_genes,
             pct_label = paste0(pct, "% samples"),
             stringsAsFactors = FALSE)
}))
pct_line_colors <- setNames(
  colorRampPalette(c("#d7191c","#1a9641","#2c7bb6"))(length(PCT_LIST)),
  paste0(PCT_LIST, "% samples")
)
combo_pts <- do.call(rbind, lapply(PCT_LIST, function(pct) {
  min_s <- ceiling(pct / 100 * n_samp)
  do.call(rbind, lapply(CUTOFF_LIST, function(cut) {
    data.frame(log2_cpm  = log2(cut + 1),
               n_genes   = sum(rowSums(cpm_mat > cut) >= min_s),
               pct_label = paste0(pct, "% samples"),
               pt_label  = as.character(sum(rowSums(cpm_mat > cut) >= min_s)),
               stringsAsFactors = FALSE)
  }))
}))

p_p2 <- ggplot(sens_df, aes(x = log2_cpm, y = n_genes,
                              color = pct_label, group = pct_label)) +
  geom_line(linewidth = 1) +
  geom_vline(data = data.frame(xint = log2(CUTOFF_LIST + 1)),
             aes(xintercept = xint), color = "grey65", linetype = "dashed",
             linewidth = 0.7, inherit.aes = FALSE) +
  annotate("text", x = log2(CUTOFF_LIST + 1) + 0.06, y = Inf,
           vjust = 1.6, hjust = 0, label = paste0("CPM=", CUTOFF_LIST),
           color = "grey40", size = 2.4) +
  geom_point(data = combo_pts, aes(x = log2_cpm, y = n_genes, color = pct_label),
             size = 3.5, shape = 16, show.legend = FALSE) +
  geom_label_repel(data = combo_pts,
                   aes(x = log2_cpm, y = n_genes, label = pt_label, color = pct_label),
                   size = 2.8, show.legend = FALSE, box.padding = 0.35,
                   max.overlaps = 30, segment.size = 0.3, fill = "white",
                   alpha = 0.9, label.size = 0.2) +
  scale_color_manual(values = pct_line_colors, name = "% samples\nrequired") +
  scale_x_continuous(breaks = log2(c(0, 1, 2, 5, 10, 15, 20, 25, 30) + 1),
                     labels = paste0("log2(", c(0, 1, 2, 5, 10, 15, 20, 25, 30), "+1)"),
                     guide = guide_axis(angle = 45)) +
  labs(title    = "Genes retained vs. CPM threshold",
       subtitle = paste0(length(PCT_LIST), " PCT lines x ", length(CUTOFF_LIST),
                         " CUTOFF | ", n_samp, " samples"),
       x = "log2(CPM+1)", y = "Genes retained") +
  theme_cowplot(12)

# ── Panel 3: Sample color legend ─────────────────────────────────────────────
n_per_col  <- 30L
legend_df  <- data.frame(
  Sample  = names(samp_pal),
  col_num = ceiling(seq_len(n_samp) / n_per_col),
  row_num = rev(((seq_len(n_samp) - 1L) %% n_per_col) + 1L),
  stringsAsFactors = FALSE
)
n_leg_cols <- max(legend_df$col_num)
p_p3 <- ggplot(legend_df, aes(x = col_num, y = row_num, color = Sample)) +
  geom_point(size = 2.0, shape = 15) +
  geom_text(aes(x = col_num + 0.12, label = Sample),
            hjust = 0, size = 1.8, color = "grey15") +
  scale_color_manual(values = samp_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
  theme_void() +
  labs(title = paste0("Sample color legend (n=", n_samp, ")"))

# ── Panel 4: Raw count distribution ──────────────────────────────────────────
if (!is.null(raw_shared) && ncol(raw_shared) > 0) {
  raw_ids   <- colnames(raw_shared)
  raw_mat   <- as.matrix(raw_shared)
  raw_list  <- setNames(lapply(raw_ids, function(s) log2(raw_mat[, s] + 1)), raw_ids)
  poly_raw  <- .make_freq_poly(raw_list, breaks_h)

  raw_annot <- lapply(seq_along(CUTOFF_LIST), function(i) {
    t      <- CUTOFF_LIST[i]
    n_pass <- sum(rowSums(raw_mat > t) >= ceiling(0.20 * length(raw_ids)))
    annotate("text", x = log2(t + 1) + 0.10, y = Inf,
             vjust = 1.8 + (i - 1) * 3.8,
             label = paste0("t=", t, " counts: ", n_pass, " genes (>=20%)"),
             color = cutoff_colors[[as.character(t)]], size = 2.4, hjust = 0)
  })

  p_p4 <- ggplot(poly_raw, aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
    geom_line(alpha = 0.45, linewidth = 0.40) +
    lapply(CUTOFF_LIST, function(t)
      geom_vline(xintercept = log2(t + 1),
                 color = cutoff_colors[as.character(t)],
                 linetype = "dashed", linewidth = 0.85)) +
    raw_annot +
    scale_color_manual(values = samp_pal[raw_ids], guide = "none") +
    scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                       guide = guide_axis(angle = 90)) +
    labs(title    = paste0("Raw count distribution — ", length(raw_ids), " samples"),
         subtitle = "log2(count+1) | bin 0.5 | zeros excluded",
         x = "log2(count+1) bin", y = "Gene count") +
    theme_cowplot(12)
} else {
  p_p4 <- ggplot() +
    labs(title = "Raw counts not provided") + theme_void()
}

grid_4 <- plot_grid(p_p1, p_p2, p_p3, p_p4, ncol = 2, nrow = 2)
ggsave(file.path(OUT_DIR, "QC_01_filter_grid.png"),
       grid_4, width = 16, height = 10, dpi = 150, bg = "white")
message("  Saved: QC_01_filter_grid.png")

# ─────────────────────────────────────────────────────────────────────────────
# 3c. STANDALONE ALL-SAMPLES HISTOGRAM — IQR OUTLIERS LABELED
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3c] Standalone all-samples histogram ===")

# Label points at mode bin of each outlier's polygon
.label_pts <- function(poly_df, label_ids) {
  if (length(label_ids) == 0) return(NULL)
  do.call(rbind, lapply(label_ids, function(s) {
    df_s <- poly_df[poly_df$Sample == s, ]
    if (nrow(df_s) == 0) return(NULL)
    row         <- df_s[which.max(df_s$count), , drop = FALSE]
    row$display <- s
    row
  }))
}

label_pts_3c <- .label_pts(poly_all, iqr_outliers)

p_3c <- ggplot(poly_all, aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.55, linewidth = 0.40) +
  cpm_vlines +
  {if (!is.null(label_pts_3c) && nrow(label_pts_3c) > 0)
    geom_label_repel(data = label_pts_3c,
                     aes(x = bin_mid, y = count, label = display),
                     inherit.aes = FALSE,
                     size = 3.2, box.padding = 0.45, max.overlaps = 60,
                     segment.size = 0.30, segment.color = "grey50",
                     fill = "white", alpha = 0.92, label.size = 0.25,
                     fontface = "bold")
  else NULL} +
  scale_color_manual(values = samp_pal, guide = "none") +
  scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                     guide = guide_axis(angle = 90)) +
  labs(title    = paste0("CPM distribution — all ", n_samp, " samples"),
       subtitle = paste0("log2(CPM+1) | bin 0.5 | zeros excluded | ",
                         length(iqr_outliers), " IQR outliers labeled"),
       x = "log2(CPM+1) bin", y = "Gene count") +
  theme_cowplot(12)

hist_3c <- plot_grid(p_3c, p_p3, ncol = 2, rel_widths = c(3, 0.55 * n_leg_cols))
ggsave(file.path(OUT_DIR, "QC_02_hist_all_samples.png"),
       hist_3c,
       width  = 16 + 3 * n_leg_cols,
       height = 10, dpi = 150, bg = "white")
message("  Saved: QC_02_hist_all_samples.png")

# Export per-sample histogram bin counts for follow-up analysis
hist_export <- data.frame(
  Sample  = poly_all$Sample,
  bin_lo  = poly_all$bin_mid - BIN_WIDTH / 2,
  bin_mid = poly_all$bin_mid,
  bin_hi  = poly_all$bin_mid + BIN_WIDTH / 2,
  count   = poly_all$count,
  stringsAsFactors = FALSE
)
write.csv(hist_export,
          file.path(OUT_DIR, "QC_02_hist_sample_bin_data.csv"),
          row.names = FALSE)
message("  Saved: QC_02_hist_sample_bin_data.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 3c-VARIANT. GROUPED HISTOGRAMS BY EACH METADATA COLUMN
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3c-variant] Grouped histograms by metadata column ===")

for (mc in meta_cols) {
  message("  Grouping by: ", mc)

  grp_vec_mc  <- setNames(as.character(gf[shared_samples, mc]), shared_samples)
  poly_mc     <- .make_freq_poly(cpm_list, breaks_h, grp_vec_mc)
  grp_pal_mc  <- col_colors[[mc]]
  label_mc    <- .label_pts(poly_mc, iqr_outliers)

  # Group color legend
  grp_vals_mc <- names(grp_pal_mc)
  leg_mc_df   <- data.frame(Group = grp_vals_mc,
                              x = 1L,
                              y = rev(seq_along(grp_vals_mc)),
                              stringsAsFactors = FALSE)
  p_leg_mc <- ggplot(leg_mc_df, aes(x = x, y = y, color = Group)) +
    geom_point(size = 5, shape = 15) +
    geom_text(aes(x = x + 0.08, label = Group), hjust = 0, size = 3.5, color = "grey15") +
    scale_color_manual(values = grp_pal_mc, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +
    scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
    theme_void() +
    labs(title = mc)

  p_mc <- ggplot(poly_mc, aes(x = bin_mid, y = count, group = Sample, color = Group)) +
    geom_line(alpha = 0.65, linewidth = 0.50) +
    cpm_vlines +
    {if (!is.null(label_mc) && nrow(label_mc) > 0)
      geom_label_repel(data = label_mc,
                       aes(x = bin_mid, y = count, label = display),
                       inherit.aes = FALSE,
                       size = 3.2, box.padding = 0.45, max.overlaps = 60,
                       segment.size = 0.30, segment.color = "grey50",
                       fill = "white", alpha = 0.92, label.size = 0.25,
                       fontface = "bold")
    else NULL} +
    scale_color_manual(values = grp_pal_mc, name = mc) +
    scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                       guide = guide_axis(angle = 90)) +
    labs(title    = paste0("CPM distribution grouped by: ", mc),
         subtitle = paste0(n_samp, " individual sample lines | log2(CPM+1) | ",
                           "IQR outliers labeled"),
         x = "log2(CPM+1) bin", y = "Gene count") +
    theme_cowplot(12) +
    theme(legend.position = "none")

  safe_mc    <- gsub("[^A-Za-z0-9_]", "_", mc)
  hist_mc    <- plot_grid(p_mc, p_leg_mc, ncol = 2, rel_widths = c(4, 1))
  ggsave(file.path(OUT_DIR, paste0("QC_03_hist_by_", safe_mc, ".png")),
         hist_mc, width = 14, height = 7, dpi = 150, bg = "white")
  message("    Saved: QC_03_hist_by_", safe_mc, ".png")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3d. PER-GROUP AVERAGED HISTOGRAM (first metadata column as "subject" grouping)
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3d] Per-group averaged histogram — ", meta_cols[1], " ===")

mc1       <- meta_cols[1]
grp_vec_1 <- setNames(as.character(gf[shared_samples, mc1]), shared_samples)
poly_avg  <- .make_freq_poly_avg(cpm_list, grp_vec_1, breaks_h)
grp_pal_1 <- col_colors[[mc1]]

# Outlier detection on per-group mean log2(CPM+1)
avg_vals  <- setNames(
  sapply(names(grp_pal_1), function(g) {
    sids <- names(grp_vec_1)[grp_vec_1 == g]
    sids <- sids[sids %in% names(cpm_list)]
    if (length(sids) == 0) return(NA_real_)
    mat <- do.call(cbind, cpm_list[sids])
    mean(if (ncol(mat) > 1) rowMeans(mat) else mat[, 1])
  }),
  names(grp_pal_1)
)
avg_vals    <- avg_vals[!is.na(avg_vals)]
iqr_avg     <- IQR(avg_vals)
q1_avg      <- quantile(avg_vals, 0.25)
q3_avg      <- quantile(avg_vals, 0.75)
avg_outliers <- names(avg_vals)[avg_vals < q1_avg - 1.5 * iqr_avg |
                                 avg_vals > q3_avg + 1.5 * iqr_avg]

label_avg <- if (length(avg_outliers) > 0) {
  do.call(rbind, lapply(avg_outliers, function(g) {
    df_g <- poly_avg[poly_avg$Group == g, ]
    if (nrow(df_g) == 0) return(NULL)
    row         <- df_g[which.max(df_g$count), , drop = FALSE]
    row$display <- g
    row
  }))
} else { NULL }

# Group color legend for mc1
grp_vals_1 <- names(grp_pal_1)
leg_1_df   <- data.frame(Group = grp_vals_1,
                          x = 1L, y = rev(seq_along(grp_vals_1)),
                          stringsAsFactors = FALSE)
p_leg_avg  <- ggplot(leg_1_df, aes(x = x, y = y, color = Group)) +
  geom_point(size = 5, shape = 15) +
  geom_text(aes(x = x + 0.08, label = Group), hjust = 0, size = 3.5, color = "grey15") +
  scale_color_manual(values = grp_pal_1, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
  theme_void() +
  labs(title = mc1)

p_avg <- ggplot(poly_avg, aes(x = bin_mid, y = count, group = Group, color = Group)) +
  geom_line(alpha = 0.80, linewidth = 0.80) +
  cpm_vlines +
  {if (!is.null(label_avg) && nrow(label_avg) > 0)
    geom_label_repel(data = label_avg,
                     aes(x = bin_mid, y = count, label = display),
                     inherit.aes = FALSE,
                     size = 3.5, box.padding = 0.45, max.overlaps = 50,
                     segment.size = 0.30, fill = "white",
                     alpha = 0.92, label.size = 0.25, fontface = "bold")
  else NULL} +
  scale_color_manual(values = grp_pal_1, name = mc1) +
  scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                     guide = guide_axis(angle = 90)) +
  labs(title    = paste0("CPM distribution — averaged per ", mc1),
       subtitle = paste0("One line per ", mc1, " value (mean log2(CPM+1) across members) | ",
                         "IQR outlier groups labeled"),
       x = "log2(CPM+1) bin", y = "Gene count") +
  theme_cowplot(12) +
  theme(legend.position = "none")

safe_mc1  <- gsub("[^A-Za-z0-9_]", "_", mc1)
hist_avg  <- plot_grid(p_avg, p_leg_avg, ncol = 2, rel_widths = c(4, 1))
ggsave(file.path(OUT_DIR, paste0("QC_04_hist_avg_per_", safe_mc1, ".png")),
       hist_avg, width = 14, height = 7, dpi = 150, bg = "white")
message("  Saved: QC_04_hist_avg_per_", safe_mc1, ".png")

# ─────────────────────────────────────────────────────────────────────────────
# 3e. ZOOMED OVERLAY HISTOGRAM — y <= 4000, IQR OUTLIERS BOLDED
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3e] Zoomed histogram, outliers bolded ===")

poly_bg   <- poly_all[!poly_all$Sample %in% iqr_outliers, ]
poly_bold <- poly_all[ poly_all$Sample %in% iqr_outliers, ]

p_3e <- ggplot() +
  geom_line(data = poly_bg,
            aes(x = bin_mid, y = count, group = Sample, color = Sample),
            alpha = 0.50, linewidth = 0.40) +
  cpm_vlines +
  {if (nrow(poly_bold) > 0)
    geom_line(data = poly_bold,
              aes(x = bin_mid, y = count, group = Sample, color = Sample),
              alpha = 1.00, linewidth = 1.40)
  else NULL} +
  {if (!is.null(label_pts_3c) && nrow(label_pts_3c) > 0)
    geom_label_repel(data = label_pts_3c,
                     aes(x = bin_mid, y = count, label = display),
                     inherit.aes = FALSE,
                     size = 3.5, box.padding = 0.45, max.overlaps = 60,
                     segment.size = 0.30, segment.color = "grey50",
                     fill = "white", alpha = 0.92, label.size = 0.25,
                     fontface = "bold")
  else NULL} +
  coord_cartesian(ylim = c(0, 4000)) +
  scale_color_manual(values = samp_pal, guide = "none") +
  scale_x_continuous(breaks = .log2_breaks, labels = .log2_labels,
                     guide = guide_axis(angle = 90)) +
  labs(title    = paste0("CPM distribution — zoomed y<=4000 | ",
                         length(iqr_outliers), " IQR outliers bolded"),
       subtitle = paste0("All ", n_samp, " samples | log2(CPM+1) | bin 0.5 | zeros excluded"),
       x = "log2(CPM+1) bin", y = "Gene count") +
  theme_cowplot(12)

hist_3e <- plot_grid(p_3e, p_p3, ncol = 2, rel_widths = c(3, 0.55 * n_leg_cols))
ggsave(file.path(OUT_DIR, "QC_05_hist_zoom_outliers_bold.png"),
       hist_3e,
       width  = 16 + 3 * n_leg_cols,
       height = 10, dpi = 150, bg = "white")
message("  Saved: QC_05_hist_zoom_outliers_bold.png")

# ─────────────────────────────────────────────────────────────────────────────
# 6B. SAMPLE-SAMPLE PEARSON CORRELATION HEATMAP
# ─────────────────────────────────────────────────────────────────────────────
message("=== [6B] Sample-sample Pearson correlation heatmap ===")

log2_cpm_mat <- log2(cpm_mat + 1)
cor_m        <- cor(log2_cpm_mat, method = "pearson")
hc_m         <- hclust(as.dist(1 - cor_m), method = "complete")
cor_colors_h <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Annotation data frame — all metadata columns, values as character
ann_df <- gf[shared_samples, meta_cols, drop = FALSE]
for (mc in meta_cols) ann_df[[mc]] <- as.character(ann_df[[mc]])
rownames(ann_df) <- shared_samples

# Auto annotation colors: one named vector per metadata column
ann_colors_list <- lapply(meta_cols, function(mc) {
  vals <- sort(unique(ann_df[[mc]][!is.na(ann_df[[mc]])]))
  col_colors[[mc]][vals]
})
names(ann_colors_list) <- meta_cols

# Square PNG: scale pixel size with sample count
px_per_samp <- max(10L, min(40L, as.integer(3600 / n_samp)))
img_dim     <- max(2000L, as.integer(n_samp * px_per_samp + 800L))
fsize_col   <- max(14, min(36, as.integer(36 * 40 / n_samp)))

# Build annotation legend panel separately (placed at the top of the figure).
# Skip columns with > 25 unique values (e.g. sort_order) — too many levels to key.
leg_cols <- meta_cols[sapply(meta_cols, function(mc) length(unique(ann_df[[mc]])) <= 25)]

if (length(leg_cols) > 0) {
  ann_leg_plots <- lapply(leg_cols, function(mc) {
    vals <- sort(unique(ann_df[[mc]]))
    df   <- data.frame(x = 0, y = -seq_along(vals),
                       lab = vals, col = ann_colors_list[[mc]][vals],
                       stringsAsFactors = FALSE)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(aes(color = I(col)), shape = 15, size = 3) +
      geom_text(aes(x = 0.15, label = lab), hjust = 0, size = 2.4, color = "grey15") +
      scale_x_continuous(limits = c(-0.3, 3.5), expand = expansion(0)) +
      labs(title = mc) +
      theme_void() +
      theme(plot.title  = element_text(size = 9, face = "bold", hjust = 0),
            plot.margin = margin(2, 8, 2, 4))
  })
  p_ann_legend <- plot_grid(plotlist = ann_leg_plots, nrow = 1)
  max_n_lev    <- max(sapply(leg_cols, function(mc) length(unique(ann_df[[mc]]))))
  ann_prop     <- min(0.20, max(0.06, max_n_lev * 0.008))
} else {
  p_ann_legend <- NULL
  ann_prop     <- 0
}

# Draw heatmap without annotation legend (color-scale bar is kept via legend = TRUE)
ph_obj <- pheatmap(cor_m,
         cluster_rows      = hc_m,
         cluster_cols      = hc_m,
         annotation_col    = ann_df,
         annotation_row    = ann_df,
         annotation_colors = ann_colors_list,
         color             = cor_colors_h,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         fontsize_col      = fsize_col,
         fontsize_row      = fsize_col,
         fontsize          = max(10, min(14, as.integer(14 * 40 / n_samp))),
         main              = paste0("Sample-sample Pearson | log2(CPM+1) | ",
                                    nrow(cpm_mat), " genes | ", n_samp, " samples"),
         border_color      = NA,
         annotation_legend = FALSE,
         silent            = TRUE)

# Compose: annotation legend strip at top, heatmap filling the rest
p_heat   <- ggdraw() + draw_grob(ph_obj$gtable)
combined <- if (!is.null(p_ann_legend)) {
  plot_grid(p_ann_legend, p_heat, ncol = 1, rel_heights = c(ann_prop, 1 - ann_prop))
} else {
  p_heat
}

png(file.path(OUT_DIR, "QC_06_correlation_heatmap.png"),
    width = img_dim, height = img_dim, res = 150, type = "cairo")
print(combined)
dev.off()
message("  Saved: QC_06_correlation_heatmap.png")

# Per-gene variance CSV (top 500 drivers)
gene_var  <- sort(rowVars(log2_cpm_mat), decreasing = TRUE)
top_500   <- names(gene_var)[seq_len(min(500L, length(gene_var)))]
driver_df <- data.frame(
  ENSG             = top_500,
  Gene_Symbol      = gene_sym_map$Gene_Symbol[match(top_500, gene_sym_map$ENSG)],
  log2cpm_variance = gene_var[top_500],
  stringsAsFactors = FALSE
)
write.csv(driver_df, file.path(OUT_DIR, "QC_06_correlation_gene_drivers.csv"),
          row.names = FALSE)
message("  Saved: QC_06_correlation_gene_drivers.csv")

# Per-sample mean Pearson correlation
cor_nd        <- cor_m
diag(cor_nd)  <- NA
samp_mean_cor <- sort(rowMeans(cor_nd, na.rm = TRUE))
write.csv(data.frame(Sample       = names(samp_mean_cor),
                     mean_pearson = round(samp_mean_cor, 4),
                     stringsAsFactors = FALSE),
          file.path(OUT_DIR, "QC_06_correlation_sample_mean.csv"), row.names = FALSE)
message("  Saved: QC_06_correlation_sample_mean.csv")

# ─────────────────────────────────────────────────────────────────────────────
# PCA — ONE PNG PER METADATA COLUMN, OUTLIERS LABELED LARGE
# ─────────────────────────────────────────────────────────────────────────────
message("=== [PCA] ===")

# Remove zero-variance genes before prcomp
nz_genes     <- rowVars(log2_cpm_mat) > 0
pca_mat      <- t(log2_cpm_mat[nz_genes, , drop = FALSE])

pca_res  <- prcomp(pca_mat, scale. = TRUE, center = TRUE)
pct_var  <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_base <- data.frame(
  Sample = rownames(pca_res$x),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  stringsAsFactors = FALSE
)

# Outlier detection: 1.5×IQR on PC1 and PC2 independently
.iqr_out <- function(x) {
  q  <- quantile(x, c(0.25, 0.75))
  iq <- IQR(x)
  x < q[1] - 1.5 * iq | x > q[2] + 1.5 * iq
}
pca_out_idx <- .iqr_out(pca_base$PC1) | .iqr_out(pca_base$PC2)
pca_out     <- pca_base$Sample[pca_out_idx]
message("  PCA outliers: ", length(pca_out),
        if (length(pca_out) > 0) paste0(" — ", paste(pca_out, collapse = ", ")) else "")

pca_label_df <- pca_base[pca_out_idx, , drop = FALSE]

# Save PCA coordinates with full metadata
write.csv(cbind(pca_base, gf[pca_base$Sample, meta_cols, drop = FALSE]),
          file.path(OUT_DIR, "QC_07_PCA_coordinates.csv"), row.names = FALSE)
message("  Saved: QC_07_PCA_coordinates.csv")

# One PNG per metadata column
for (mc in meta_cols) {
  grp_pca   <- as.character(gf[pca_base$Sample, mc])
  pca_df_mc <- cbind(pca_base, Group = grp_pca, stringsAsFactors = FALSE)

  grp_pal_pca   <- col_colors[[mc]]
  grp_shape_pca <- col_shapes[[mc]]

  # Ensure all levels present in pca_df_mc are in the palette/shape vectors
  present_grps <- unique(grp_pca[!is.na(grp_pca)])
  missing_fill  <- setdiff(present_grps, names(grp_pal_pca))
  if (length(missing_fill) > 0) {
    extra <- setNames(rep("#888888", length(missing_fill)), missing_fill)
    grp_pal_pca <- c(grp_pal_pca, extra)
  }
  missing_shape <- setdiff(present_grps, names(grp_shape_pca))
  if (length(missing_shape) > 0) {
    extra_s <- setNames(rep(21L, length(missing_shape)), missing_shape)
    grp_shape_pca <- c(grp_shape_pca, extra_s)
  }

  p_pca <- ggplot(pca_df_mc,
                  aes(x = PC1, y = PC2, fill = Group, shape = Group)) +
    geom_point(size = 5, color = "grey25", stroke = 0.4) +
    {if (nrow(pca_label_df) > 0)
      geom_label_repel(data = cbind(pca_label_df,
                                    Group = as.character(gf[pca_label_df$Sample, mc])),
                       aes(x = PC1, y = PC2, label = Sample),
                       inherit.aes = FALSE,
                       size = 5.5, fontface = "bold",
                       box.padding   = 0.55,
                       point.padding = 0.40,
                       max.overlaps  = 80,
                       segment.size  = 0.45,
                       segment.color = "grey40",
                       fill          = "white",
                       alpha         = 0.92,
                       label.size    = 0.30)
    else NULL} +
    scale_fill_manual(values  = grp_pal_pca,   name = mc) +
    scale_shape_manual(values = grp_shape_pca, name = mc) +
    labs(title    = paste0("PCA — colored by: ", mc),
         subtitle = paste0(n_samp, " samples | ", sum(nz_genes),
                           " genes | outliers labeled"),
         x = paste0("PC1: ", round(pct_var[1] * 100, 1), "% variance"),
         y = paste0("PC2: ", round(pct_var[2] * 100, 1), "% variance")) +
    theme_linedraw(base_size = 14) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title     = element_text(size = 7),
          legend.text      = element_text(size = 6),
          legend.key.size  = unit(0.35, "cm"),
          legend.position  = "right",
          axis.text        = element_text(size = 12),
          axis.title       = element_text(size = 16)) +
    guides(fill  = guide_legend(ncol = 2, override.aes = list(size = 2.5)),
           shape = guide_legend(ncol = 2, override.aes = list(size = 2.5)))

  safe_mc_pca <- gsub("[^A-Za-z0-9_]", "_", mc)
  ggsave(file.path(OUT_DIR, paste0("QC_07_PCA_by_", safe_mc_pca, ".png")),
         p_pca, width = 10, height = 7, dpi = 150, bg = "white")
  message("  Saved: QC_07_PCA_by_", safe_mc_pca, ".png")
}

message("=== Done. All outputs in: ", OUT_DIR, " ===")
