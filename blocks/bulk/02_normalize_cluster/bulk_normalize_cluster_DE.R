###############################################################################
# PROMOTER Bulk RNA-seq Analysis Pipeline
# CD11b-sorted myeloid cells | Continuous pain (QST-PPT) | DESeq2 + GSEA
#
# Inputs:
#   - PROMOTOR_BulkSeq_20260216_final_count.csv   raw integer counts
#   - CPM_PROMOTOR_BulkSeq_20260216_final_count.txt  CPM (Symbol + Gene_Symbol first)
#   - example_pain_metadata.csv    study_id / qst_ppt_tr_avg_v1 / etc.
#   - groupfile_by_Day/RunOrder/Visit/Sample.csv  (no header, 4 cols)
#
# NOTE: PROMO24 V1 has 3 technical replicates (S3, S21, S93) — kept separate.
#
# Pipeline sections:
#   1.  Load data
#   2.  Build sample metadata (colData)
#   3.  Filter lowly expressed genes (CPM > 2)
#   4.  VST normalisation
#   5.  Variance-based variable gene selection
#   6.  Hierarchical clustering of samples (Pearson | Ward.D2)
#   7.  K-means clustering of variable genes + heatmaps (k = 3,4,5,6)
#   8.  DESeq2 helper function
#   9.  DESeq2 ~ qst_ppt_tr_avg_v1 (V1 samples, continuous PPT)
#   10. Example gene visualisation (top DE genes, box+jitter)
#   11. GSEA (Hallmark gene sets, fgsea, dotplot)
#   12. Summary
###############################################################################
.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggrepel)
  library(cowplot)
  library(fgsea)
  library(msigdbr)
  library(matrixStats)
  library(scales)
  library(cluster)
  library(ggVennDiagram)
  library(ggalluvial)
})
# clusterProfiler loaded separately — ggtree/enrichplot dependency can be tricky.
# Script continues even if unavailable; run_ora_cp() returns NULL gracefully.
HAS_CP <- requireNamespace("clusterProfiler", quietly = TRUE)
if (HAS_CP) suppressPackageStartupMessages(library(clusterProfiler)) else
  message("  NOTE: clusterProfiler not found — CP ORA comparison will be skipped")

set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# 0. PATHS  ── edit these
# ─────────────────────────────────────────────────────────────────────────────
RAW_COUNTS <- "PROMOTOR_BulkSeq_20260216_final_count.txt"
CPM_FILE   <- "CPM_PROMOTOR_BulkSeq_20260216_final_count.txt"
PAIN_META  <- "example_pain_metadata.csv"
GF_DAY     <- "groupfile_by_Day.csv"
GF_RUN     <- "groupfile_by_RunOrder.csv"
GF_VISIT   <- "groupfile_by_Visit.csv"
GF_SAMPLE  <- "groupfile_by_Sample.csv"
OUT_DIR    <- paste0("PROMOTER_pipeline_output_", format(Sys.Date(), "%Y%m%d"))
dir.create(OUT_DIR, showWarnings = FALSE)

DIR_QC      <- file.path(OUT_DIR, "QC")
DIR_HVG     <- file.path(OUT_DIR, "HVG")
DIR_OUTLIER <- file.path(OUT_DIR, "Outlier_QC")
DIR_CLUST   <- file.path(OUT_DIR, "Gene_Clustering")
DIR_ORA     <- file.path(OUT_DIR, "ORA")
DIR_DESEQ2  <- file.path(OUT_DIR, "DESeq2")
DIR_GSEA    <- file.path(OUT_DIR, "GSEA")
for (d in c(DIR_QC, DIR_HVG, DIR_OUTLIER, DIR_CLUST, DIR_ORA, DIR_DESEQ2, DIR_GSEA))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
CPM_THRESH      <- 5         # CPM threshold for low-expression filter
DESEQ_MIN_COUNT <- CPM_THRESH
                             # DESeq2 raw-count floor — kept in sync with CPM_THRESH
MIN_PCT_SAMPLES <- 0.20    # min fraction of V1 PROMO samples a gene must be expressed in
                             # used by both the CPM filter (Section 3) and the
                             # DESeq2 internal filter inside run_deseq2()
CUTOFF_LIST     <- c(5, 10, 20)
                             # candidate CPM thresholds shown in the QC histogram (Section 3b)
PCT_LIST        <- c(10, 20, 25)
                             # % sample thresholds shown as lines in Panel 2 sensitivity curve
N_VAR_GENES    <- 2000     # top variable genes to use for clustering/heatmap
K_VALUES       <- 6          # k-means k values to test (set to 6 only — speeds up the loop)
N_EXAMPLE_GENES <- 12      # top DE genes to plot individually
GSEA_TOP_N     <- 30       # pathways shown in GSEA dotplot
LFC_THRESH     <- 0.585    # post-hoc |LFC| filter (not passed to results())
ALPHA          <- 0.05
# Subject IDs used in outlier diagnostics (6C, 6D, 6E) and HVG influence analysis.
# These samples remain in all analyses; this list flags them for visualization only.
HVG_EXCLUDE_SUBJECTS <- c("PROMO02", "PROMO04")

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
message("=== [1] Loading data ===")

read_gf <- function(path) {
  read.csv(path, header = FALSE,
           col.names = c("SampleID","Group","Color","Num"),
           stringsAsFactors = FALSE)
}
gf_day   <- read_gf(GF_DAY)
gf_run   <- read_gf(GF_RUN)
gf_visit <- read_gf(GF_VISIT)
gf_samp  <- read_gf(GF_SAMPLE)

raw <- read.delim(RAW_COUNTS, row.names = 1, check.names = FALSE)
message("  Raw count matrix: ", nrow(raw), " genes x ", ncol(raw), " samples")

cpm_raw <- read.delim(CPM_FILE, row.names = 1, check.names = FALSE)
gene_sym_map <- data.frame(
  ENSG        = rownames(cpm_raw),
  Gene_Symbol = cpm_raw[["Gene_Symbol"]],
  stringsAsFactors = FALSE
)
cpm <- cpm_raw[ , -1, drop = FALSE]

# Pain metadata — strip underscore from study_id (PROMO_13 → PROMO13)
pain_raw <- read.csv(PAIN_META, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
pain_raw$SubjectID <- gsub("_", "", pain_raw$study_id)
message("  Pain metadata: ", nrow(pain_raw), " subjects")
message("  PPT range: ", round(min(pain_raw$qst_ppt_tr_avg_v1, na.rm=TRUE), 2),
        " – ", round(max(pain_raw$qst_ppt_tr_avg_v1, na.rm=TRUE), 2))

shared_cols <- intersect(colnames(raw), colnames(cpm))
raw <- raw[ , shared_cols]
cpm <- cpm[ , shared_cols]
message("  Shared samples (raw + CPM): ", length(shared_cols))

# Drop excluded subjects from all matrices before any analysis
excl_pattern <- paste(HVG_EXCLUDE_SUBJECTS, collapse = "|")
shared_cols  <- shared_cols[!grepl(excl_pattern, shared_cols)]
raw          <- raw[, shared_cols, drop = FALSE]
cpm          <- cpm[, shared_cols, drop = FALSE]
message("  Excluded: ", paste(HVG_EXCLUDE_SUBJECTS, collapse = ", "),
        " — ", length(shared_cols), " samples remaining")

# ─────────────────────────────────────────────────────────────────────────────
# 2. BUILD SAMPLE METADATA
# ─────────────────────────────────────────────────────────────────────────────
message("=== [2] Building colData ===")

dedup <- function(df) df[!duplicated(df$SampleID), ]

meta <- data.frame(SampleID = shared_cols, stringsAsFactors = FALSE) %>%
  left_join(dedup(gf_day)  [, c("SampleID","Group")] %>% rename(Day      = Group), by = "SampleID") %>%
  left_join(dedup(gf_run)  [, c("SampleID","Group")] %>% rename(RunOrder = Group), by = "SampleID") %>%
  left_join(dedup(gf_visit)[, c("SampleID","Group","Color")] %>%
              rename(Visit = Group, VisitColor = Color),             by = "SampleID") %>%
  left_join(dedup(gf_samp) [, c("SampleID","Group","Color")] %>%
              rename(SampleGroup = Group, SampleColor = Color),      by = "SampleID")

meta$SubjectID <- str_extract(meta$SampleID, "PROMO\\d+|PerlmanNorm\\d+")
meta$IsPerlman <- grepl("PerlmanNorm", meta$SampleID)
meta$Visit     <- factor(meta$Visit, levels = c("Visit1","Visit2","Visit3"))
meta$DayNum    <- as.integer(str_extract(meta$Day, "\\d+"))

meta <- meta %>%
  left_join(
    pain_raw %>% select(SubjectID, qst_ppt_tr_avg_v1, qst_ppt_tr_avg_v2,
                        promis_global07_v1, promis_global07_v2, age, sex_at_birth.factor) %>%
      rename(Sex = sex_at_birth.factor, Age = age, PROMIS = promis_global07_v1, PROMIS_v2 = promis_global07_v2),
    by = "SubjectID"
  )

# Numerical PROMIS change from Visit 1 to Visit 2 (raw, unscaled), as V2 - V1.
# Negative = PROMIS dropped (less pain / treatment response); positive = worsened.
# NA for subjects missing either v1 or v2 — used to subset to v1+v2 samples in 11.6.
meta$promis_diff <- meta$PROMIS_v2 - meta$PROMIS

# Center and scale PPT to improve GLM convergence (as DESeq2 recommends)
# Original values are retained in qst_ppt_raw for plotting
meta$qst_ppt_raw    <- meta$qst_ppt_tr_avg_v1
meta$qst_ppt_tr_avg_v1 <- as.numeric(scale(meta$qst_ppt_tr_avg_v1))
meta$qst_ppt_tr_avg_v2 <- as.numeric(scale(meta$qst_ppt_tr_avg_v2))

rownames(meta) <- meta$SampleID

message("  colData: ", nrow(meta), " samples | ",
        sum(!meta$IsPerlman), " PROMO + ", sum(meta$IsPerlman), " Perlman")
write.csv(meta, file.path(DIR_QC, "sample_metadata.csv"), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 3. FILTER LOWLY EXPRESSED GENES
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3] Filtering low-expression genes ===")

# Gene expression filter: a gene must have CPM > CPM_THRESH
# in >= MIN_PCT_SAMPLES of V1 PROMO samples.  Using the V1 PROMO group as the
# reference population keeps the filter anchored to the primary analysis group
# and avoids dilution from the larger cross-visit / Perlman sample set.
n_v1_promo  <- sum(meta$Visit == "Visit1" & !meta$IsPerlman, na.rm = TRUE)
min_samples <- ceiling(MIN_PCT_SAMPLES * n_v1_promo)

v1_promo_ids_filt <- rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman &
                                     rownames(meta) %in% colnames(cpm)]
cpm_v1_filt       <- as.matrix(cpm[, v1_promo_ids_filt])
keep              <- rowSums(cpm_v1_filt > CPM_THRESH) >= min_samples
raw_filt          <- raw[keep, ]
cpm_filt          <- cpm[keep, ]

message("  CPM>", CPM_THRESH, " in >= ", min_samples,
        " V1 PROMO samples (", round(MIN_PCT_SAMPLES*100), "% of ", n_v1_promo, ")")
message("  Genes: ", nrow(raw), " → ", nrow(raw_filt), " retained")

# Optional: filter artifact/uninformative gene classes before DE and clustering.
# Toggle each flag TRUE/FALSE in the PARAMETERS section (top of script) or here.
# Filtered from both raw_filt and cpm_filt so all downstream steps are consistent.
#
# FILTER_RIBOSOMAL : RPL* / RPS* ribosomal protein genes
#   These are extremely highly expressed housekeeping genes. They dominate DE
#   results in bulk RNA-seq not because of biology but because of their sheer
#   abundance, and are rarely interpretable in an immunology context.
#
# FILTER_RP_LNCRNA : RP*-\d+ style names (e.g. RP11-206L10, RP4-758J18)
#   These are BAC clone-based identifiers for poorly annotated genomic loci,
#   mostly lncRNAs and pseudogenes. They are NOT ribosomal proteins despite
#   the RP prefix. Filter if you are not studying lncRNA biology.
#
# FILTER_MITO     : MT-* mitochondrial genes
# FILTER_PSEUDO   : genes ending in -PS or containing "pseudogene" annotation
#   (requires BioMart lookup — left FALSE by default)

FILTER_RIBOSOMAL <- FALSE   # RPL*, RPS* — kept for unbiased HVG selection
FILTER_RP_LNCRNA <- FALSE    # RP\d+-* BAC clone lncRNAs
FILTER_MITO      <- FALSE   # MT-* — kept for unbiased HVG selection

exclude_genes <- c()

if (FILTER_RIBOSOMAL) {
  ribo <- gene_sym_map$ENSG[grepl("^RPL|^RPS", gene_sym_map$Gene_Symbol)]
  exclude_genes <- union(exclude_genes, ribo)
  message("  Ribosomal protein genes flagged for removal (RPL*/RPS*): ", length(ribo))
}
if (FILTER_RP_LNCRNA) {
  rp_lnc <- gene_sym_map$ENSG[grepl("^RP\\d+-", gene_sym_map$Gene_Symbol)]
  exclude_genes <- union(exclude_genes, rp_lnc)
  message("  BAC clone lncRNA genes flagged for removal (RP##-*): ", length(rp_lnc))
}
if (FILTER_MITO) {
  mito <- gene_sym_map$ENSG[grepl("^MT-", gene_sym_map$Gene_Symbol)]
  exclude_genes <- union(exclude_genes, mito)
  message("  Mitochondrial genes flagged for removal (MT-*): ", length(mito))
}

if (length(exclude_genes) > 0) {
  n_before <- nrow(raw_filt)
  keep_genes <- !rownames(raw_filt) %in% exclude_genes
  raw_filt   <- raw_filt[keep_genes, ]
  cpm_filt   <- cpm_filt[keep_genes, ]
  message("  After artifact gene removal: ", n_before, " → ", nrow(raw_filt), " genes")
}

# ─────────────────────────────────────────────────────────────────────────────────
# 3b. GENE FILTER VISUALIZATION
# Per-sample frequency polygons: each sample drawn as its own coloured line.
# y-axis = linear gene count; x-axis = raw CPM values, bin width = 0.5.
# Zeros excluded via open left bracket (x > 0 only).
# CUTOFF_LIST vertical dashed lines shown on both panels.
# ─────────────────────────────────────────────────────────────────────────────────
message("=== [3b] Gene filter visualization ===")

# Palette for CUTOFF_LIST vertical lines
cutoff_colors <- setNames(
  colorRampPalette(c("#d7191c", "#fdae61", "#1a9641", "#2c7bb6", "#7B2D8B"))(length(CUTOFF_LIST)),
  as.character(CUTOFF_LIST)
)

# Per-cutoff annotation: genes passing + top-5 gene symbols
mean_cpm_v1 <- rowMeans(cpm_v1_filt)
cutoff_info <- lapply(CUTOFF_LIST, function(t) {
  n_pass <- sum(rowSums(cpm_v1_filt > t) >= min_samples)
  top_g  <- names(sort(mean_cpm_v1[mean_cpm_v1 > t], decreasing = TRUE))[1:5]
  syms   <- gene_sym_map$Gene_Symbol[match(top_g, gene_sym_map$ENSG)]
  list(threshold = t, n_pass = n_pass, top_symbols = syms[!is.na(syms)])
})
names(cutoff_info) <- as.character(CUTOFF_LIST)

# Helper: build per-sample frequency polygon data frame from a named list of value vectors
.make_freq_poly <- function(val_list, breaks, visit_vec, cumulative = FALSE) {
  do.call(rbind, lapply(names(val_list), function(sid) {
    vals <- val_list[[sid]]
    vals <- vals[vals > 0 & vals <= max(breaks)]  # exclude zeros and values beyond plot range
    h    <- hist(vals, breaks = breaks, plot = FALSE, right = TRUE)
    counts <- if (cumulative) rev(cumsum(rev(h$counts))) else h$counts
    data.frame(Sample = sid, Visit = visit_vec[sid],
               bin_mid = h$mids, count = counts,
               stringsAsFactors = FALSE)
  }))
}

# ---- Panel 1: CPM expression distribution (V1 PROMO only) ------------------
breaks_expr <- seq(0, 16, by = 0.5)
cpm_list_v1 <- setNames(
  lapply(v1_promo_ids_filt, function(sid) log2(cpm[, sid] + 1)),
  v1_promo_ids_filt
)
visit_v1_map  <- setNames(rep("Visit1", length(v1_promo_ids_filt)), v1_promo_ids_filt)
poly_expr_df  <- .make_freq_poly(cpm_list_v1, breaks_expr, visit_v1_map)

# Unique-per-sample palette (legend shown in panel 3 of this grid)
n_v1_samp   <- length(v1_promo_ids_filt)
samp_pal_v1 <- setNames(
  colorRampPalette(c("#4575b4","#74add1","#abd9e9","#fdae61","#f46d43","#d73027",
                     "#a50026","#313695","#fee090","#e0f3f8"))(n_v1_samp),
  v1_promo_ids_filt
)

# Bin interval labels for x-axis: (lower, upper] for each 0.5-wide bin in log2(CPM+1) space
.log2_bin_breaks <- seq(0.25, 15.75, by = 0.5)
.log2_bin_labels <- paste0("(", seq(0, 15.5, by = 0.5), ",", seq(0.5, 16, by = 0.5), "]")

# Cutoff vlines at log2(CPM+1) positions
cpm_vline_layers <- lapply(CUTOFF_LIST, function(t) {
  geom_vline(xintercept = log2(t + 1),
             color      = cutoff_colors[as.character(t)],
             linetype   = "dashed", linewidth = 0.85)
})
cpm_annot_layers <- mapply(function(ci, idx) {
  annotate("text",
           x = log2(ci$threshold + 1) + 0.1, y = Inf,
           vjust = 1.8 + (idx - 1) * 3.8,
           label = paste0("t=", ci$threshold, " CPM: ", ci$n_pass, " genes"),
           color = cutoff_colors[[as.character(ci$threshold)]],
           size = 2.5, hjust = 0)
}, cutoff_info, seq_along(cutoff_info), SIMPLIFY = FALSE)

p_expr_hist <- ggplot(poly_expr_df,
                      aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.60, linewidth = 0.45) +
  cpm_vline_layers +
  cpm_annot_layers +
  scale_color_manual(values = samp_pal_v1, guide = "none") +
  scale_x_continuous(
    breaks = .log2_bin_breaks,
    labels = .log2_bin_labels,
    guide  = guide_axis(angle = 90)
  ) +
  labs(
    title    = paste0("Gene expression distribution — V1 PROMO samples (zeros excluded)"),
    subtitle = paste0(n_v1_samp, " samples | log2(CPM+1) per gene | bin width: 0.5"),
    x = "log2(CPM+1) bin", y = "Gene count"
  ) +
  theme_cowplot(12)

# ---- Panel 2: sensitivity curve — genes retained vs. log2(CPM+1) threshold ---
# 3 continuous lines (one per PCT_LIST value) sweeping over CPM thresholds.
# 9 labeled points: each of the 3 CUTOFF_LIST values marked on each PCT line.
pct_line_colors <- setNames(
  colorRampPalette(c("#d7191c", "#1a9641", "#2c7bb6"))(length(PCT_LIST)),
  paste0(PCT_LIST, "% samples")
)

cpm_sweep <- seq(0, 30, by = 0.25)

sens2_df <- do.call(rbind, lapply(PCT_LIST, function(pct) {
  min_s   <- ceiling(pct / 100 * n_v1_promo)
  n_genes <- sapply(cpm_sweep, function(t) sum(rowSums(cpm_v1_filt > t) >= min_s))
  data.frame(log2_cpm  = log2(cpm_sweep + 1),
             n_genes   = n_genes,
             pct_label = paste0(pct, "% samples"),
             stringsAsFactors = FALSE)
}))

combo_pts <- do.call(rbind, lapply(PCT_LIST, function(pct) {
  min_s <- ceiling(pct / 100 * n_v1_promo)
  do.call(rbind, lapply(CUTOFF_LIST, function(cut) {
    n_g <- sum(rowSums(cpm_v1_filt > cut) >= min_s)
    data.frame(log2_cpm  = log2(cut + 1),
               n_genes   = n_g,
               pct_label = paste0(pct, "% samples"),
               pt_label  = as.character(n_g),
               stringsAsFactors = FALSE)
  }))
}))

cutoff_vline_df <- data.frame(
  xint  = log2(CUTOFF_LIST + 1),
  label = paste0("CPM=", CUTOFF_LIST)
)

p_sensitivity <- ggplot(sens2_df,
                        aes(x = log2_cpm, y = n_genes,
                            color = pct_label, group = pct_label)) +
  geom_line(linewidth = 1) +
  geom_vline(data = cutoff_vline_df,
             aes(xintercept = xint),
             color = "grey65", linetype = "dashed", linewidth = 0.7,
             inherit.aes = FALSE) +
  annotate("text",
           x     = log2(CUTOFF_LIST + 1) + 0.06,
           y     = Inf, vjust = 1.6, hjust = 0,
           label = paste0("CPM=", CUTOFF_LIST),
           color = "grey40", size = 2.5) +
  geom_point(data = combo_pts,
             aes(x = log2_cpm, y = n_genes, color = pct_label),
             size = 3.5, shape = 16, show.legend = FALSE) +
  geom_label_repel(data = combo_pts,
                   aes(x = log2_cpm, y = n_genes,
                       label = pt_label, color = pct_label),
                   size = 2.8, show.legend = FALSE,
                   box.padding = 0.35, max.overlaps = 30,
                   segment.size = 0.3, segment.color = "grey55",
                   fill = "white", alpha = 0.9, label.size = 0.2) +
  scale_color_manual(values = pct_line_colors, name = "% samples\nrequired") +
  scale_x_continuous(
    breaks = log2(c(0, 1, 2, 5, 10, 15, 20, 25, 30) + 1),
    labels = paste0("log2(", c(0, 1, 2, 5, 10, 15, 20, 25, 30), "+1)"),
    guide  = guide_axis(angle = 45)
  ) +
  labs(title    = "Genes retained vs. CPM threshold — 9 (CUTOFF × PCT) combos",
       subtitle = paste0(length(PCT_LIST), " PCT lines × ", length(CUTOFF_LIST),
                         " CUTOFF points | ", n_v1_promo, " V1 PROMO samples"),
       x = "log2(CPM + 1)", y = "Genes retained") +
  theme_cowplot(12)

# ---- Panel 3: per-sample color legend (V1 PROMO) ----------------------------
n_per_col_leg <- 20L
legend_df <- data.frame(
  Sample  = names(samp_pal_v1),
  col_num = ceiling(seq_len(n_v1_samp) / n_per_col_leg),
  row_num = rev(((seq_len(n_v1_samp) - 1L) %% n_per_col_leg) + 1L),
  stringsAsFactors = FALSE
)
p_sample_legend <- ggplot(legend_df, aes(x = col_num, y = row_num, color = Sample)) +
  geom_point(size = 3.5, shape = 15) +
  geom_text(aes(x = col_num + 0.12, label = Sample),
            hjust = 0, size = 3.2, color = "grey15") +
  scale_color_manual(values = samp_pal_v1, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
  theme_void() +
  labs(title = paste0("Sample color legend — V1 PROMO (n = ", n_v1_samp, ")"))

# ---- Panel 4: RAW COUNT distribution — V1 PROMO samples + all CUTOFF_LIST lines ---
# Uses unnormalized integer counts (not CPM) — relevant for DESeq2 which requires raw counts.
# Thresholds (CUTOFF_LIST) applied as raw-count cutoffs at the same numeric values as CPM.
raw_v1_mat_p4 <- as.matrix(raw[, v1_promo_ids_filt])
raw_list_p4 <- setNames(
  lapply(v1_promo_ids_filt, function(sid) log2(raw[, sid] + 1)),
  v1_promo_ids_filt
)
breaks_raw  <- seq(0, 16, by = 0.5)
poly_raw_df <- .make_freq_poly(raw_list_p4, breaks_raw, visit_v1_map)

n_p4_samp   <- length(v1_promo_ids_filt)
samp_pal_p4 <- samp_pal_v1

raw_vline_layers <- lapply(CUTOFF_LIST, function(t) {
  geom_vline(xintercept = log2(t + 1),
             color      = cutoff_colors[as.character(t)],
             linetype   = "dashed", linewidth = 0.85)
})
raw_annot_layers <- lapply(seq_along(CUTOFF_LIST), function(i) {
  t      <- CUTOFF_LIST[i]
  n_pass <- sum(rowSums(raw_v1_mat_p4 > t) >= min_samples)
  annotate("text",
           x = log2(t + 1) + 0.1, y = Inf,
           vjust = 1.8 + (i - 1) * 3.8,
           label = paste0("t=", t, " counts: ", n_pass, " genes"),
           color = cutoff_colors[[as.character(t)]],
           size = 2.5, hjust = 0)
})

p_raw_hist <- ggplot(poly_raw_df,
                     aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.45, linewidth = 0.40) +
  raw_vline_layers +
  raw_annot_layers +
  scale_color_manual(values = samp_pal_p4, guide = "none") +
  scale_x_continuous(
    breaks = .log2_bin_breaks,
    labels = .log2_bin_labels,
    guide  = guide_axis(angle = 90)
  ) +
  labs(
    title    = paste0("Raw count distribution — V1 PROMO samples (zeros excluded)"),
    subtitle = paste0(n_p4_samp, " samples | log2(raw count+1) per gene | bin width: 0.5"),
    x = "log2(raw count+1) bin", y = "Gene count"
  ) +
  theme_cowplot(12)

filter_grid <- plot_grid(p_expr_hist,     p_sensitivity,
                          p_sample_legend, p_raw_hist,
                          ncol = 2, nrow = 2)
ggsave(file.path(DIR_QC, "00_filter_visualization.png"),
       filter_grid, width = 16, height = 10, dpi = 150, bg = "white")
message("  Saved: 00_filter_visualization.png")

# Per-cutoff summary table
cutoff_summary <- do.call(rbind, lapply(cutoff_info, function(ci) {
  data.frame(
    threshold      = ci$threshold,
    n_genes_pass   = ci$n_pass,
    pct_genes_pass = round(100 * ci$n_pass / nrow(cpm), 1),
    top5_genes     = paste(ci$top_symbols, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))
write.csv(cutoff_summary, file.path(DIR_QC, "00_filter_cutoff_summary.csv"),
          row.names = FALSE)
message("  Saved: 00_filter_cutoff_summary.csv")

# ─────────────────────────────────────────────────────────────────────────────────
# 3c. STANDALONE HISTOGRAM: ALL SAMPLES (freq polygon, Panel 1 method, all visits)
# Per-sample colorRampPalette (identical to Panel 1). Zeros excluded, bin = 0.5.
# Outliers detected via 1.5×IQR on per-sample median log2(CPM+1); HVG_EXCLUDE_SUBJECTS
# and PROMO26 V1/V2 always labeled. Color key embedded as Panel-3-style swatch grid.
# ─────────────────────────────────────────────────────────────────────────────────
message("=== [3c] Standalone histogram — all samples ===")

all_samp_ids_3c <- shared_cols
n_all_3c        <- length(all_samp_ids_3c)

# Visit annotation for display
visit_all_map_3c <- setNames(
  ifelse(is.na(meta[all_samp_ids_3c, "Visit"]), "Unknown",
         as.character(meta[all_samp_ids_3c, "Visit"])),
  all_samp_ids_3c
)

# Per-sample colorRampPalette — exact Panel 1 palette, expanded to all samples
samp_pal_3c <- setNames(
  colorRampPalette(c("#4575b4","#74add1","#abd9e9","#fdae61","#f46d43","#d73027",
                     "#a50026","#313695","#fee090","#e0f3f8"))(n_all_3c),
  all_samp_ids_3c
)

# Frequency polygon data — same .make_freq_poly helper and breaks as Panel 1
cpm_list_3c <- setNames(
  lapply(all_samp_ids_3c, function(sid) log2(cpm[, sid] + 1)),
  all_samp_ids_3c
)
poly_3c_df <- .make_freq_poly(cpm_list_3c, seq(0, 16, by = 0.5), visit_all_map_3c)

# Outlier detection: flag samples whose median log2(CPM+1) lies outside 1.5×IQR fence
med_3c   <- sapply(all_samp_ids_3c, function(sid) median(cpm_list_3c[[sid]]))
q1_3c    <- quantile(med_3c, 0.25)
q3_3c    <- quantile(med_3c, 0.75)
iqr_3c   <- IQR(med_3c)
auto_out_3c <- all_samp_ids_3c[med_3c < q1_3c - 1.5*iqr_3c | med_3c > q3_3c + 1.5*iqr_3c]

# PROMO26 V1 and V2 — always labeled
promo26_3c <- rownames(meta)[
  !is.na(meta$SubjectID) & meta$SubjectID == "PROMO26" &
  !is.na(meta$Visit)     & meta$Visit %in% c("Visit1", "Visit2") &
  rownames(meta) %in% all_samp_ids_3c
]

# HVG_EXCLUDE_SUBJECTS — always labeled
hvg_excl_3c <- rownames(meta)[
  !is.na(meta$SubjectID) & meta$SubjectID %in% HVG_EXCLUDE_SUBJECTS &
  rownames(meta) %in% all_samp_ids_3c
]

label_samps_3c <- unique(c(auto_out_3c, promo26_3c, hvg_excl_3c))

# Label placed at the mode (highest-count bin) of each flagged sample's polygon
label_pts_3c <- do.call(rbind, lapply(label_samps_3c, function(sid) {
  df_s        <- poly_3c_df[poly_3c_df$Sample == sid, ]
  row         <- df_s[which.max(df_s$count), , drop = FALSE]
  row$display <- sid
  row
}))

# CUTOFF vlines — same cutoff_colors palette as Section 3b
cpm_vline_3c <- lapply(CUTOFF_LIST, function(t) {
  geom_vline(xintercept = log2(t + 1),
             color      = cutoff_colors[as.character(t)],
             linetype   = "dashed", linewidth = 0.85)
})

# Main frequency polygon panel
p_hist_3c <- ggplot(poly_3c_df,
                    aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.55, linewidth = 0.40) +
  cpm_vline_3c +
  geom_label_repel(
    data        = label_pts_3c,
    aes(x = bin_mid, y = count, label = display),
    inherit.aes = FALSE,
    size        = 2.5, box.padding = 0.40, max.overlaps = 50,
    segment.size = 0.30, segment.color = "grey50",
    fill = "white", alpha = 0.90, label.size = 0.2
  ) +
  scale_color_manual(values = samp_pal_3c, guide = "none") +
  scale_x_continuous(
    breaks = .log2_bin_breaks,
    labels = .log2_bin_labels,
    guide  = guide_axis(angle = 90)
  ) +
  labs(
    title    = paste0("Gene expression distribution — all ", n_all_3c,
                      " samples (zeros excluded)"),
    subtitle = paste0("log2(CPM+1) per gene | bin width: 0.5 | ",
                      "outliers & PROMO26 V1/V2 labeled"),
    x = "log2(CPM+1) bin", y = "Gene count"
  ) +
  theme_cowplot(12)

# Panel-3-style swatch color key (exact Panel 3 code, adapted for n_all_3c samples)
n_per_col_3c <- 20L
legend_3c_df <- data.frame(
  Sample  = names(samp_pal_3c),
  col_num = ceiling(seq_len(n_all_3c) / n_per_col_3c),
  row_num = rev(((seq_len(n_all_3c) - 1L) %% n_per_col_3c) + 1L),
  stringsAsFactors = FALSE
)
p_legend_3c <- ggplot(legend_3c_df, aes(x = col_num, y = row_num, color = Sample)) +
  geom_point(size = 3.5, shape = 15) +
  geom_text(aes(x = col_num + 0.12, label = Sample),
            hjust = 0, size = 2.8, color = "grey15") +
  scale_color_manual(values = samp_pal_3c, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
  theme_void() +
  labs(title = paste0("Sample color key (n = ", n_all_3c, ")"))

n_key_cols_3c   <- max(legend_3c_df$col_num)
hist_3c_combined <- plot_grid(p_hist_3c, p_legend_3c,
                               ncol = 2,
                               rel_widths = c(3, 0.55 * n_key_cols_3c))

ggsave(file.path(DIR_QC, "00c_expr_hist_all_samples.png"),
       hist_3c_combined,
       width  = 16 + 3 * n_key_cols_3c,
       height = 10, dpi = 150, bg = "white")
message("  Saved: 00c_expr_hist_all_samples.png")

# ─────────────────────────────────────────────────────────────────────────────────
# 3d. STANDALONE HISTOGRAM: V1/V2-AVERAGED PER SUBJECT
# Per-gene rowMeans of log2(CPM+1) across Visit1 + Visit2 for each subject.
# Subjects with only one visit use that visit directly.
# Note: PROMO24 V1 has 3 technical replicates — all averaged into one subject line.
# One line per subject. PROMO07 always labeled; outliers auto-labeled (1.5×IQR).
# ─────────────────────────────────────────────────────────────────────────────────
message("=== [3d] V1/V2-averaged histogram — one line per subject ===")

# All unique subject IDs present in shared_cols
subj_ids_3d <- unique(meta$SubjectID[!is.na(meta$SubjectID) &
                                      rownames(meta) %in% shared_cols])

# Per-subject: average log2(CPM+1) across V1 + V2 (or any available visit if neither exists)
avg_cpm_list_3d <- lapply(subj_ids_3d, function(subj) {
  v12_samps <- rownames(meta)[
    !is.na(meta$SubjectID) & meta$SubjectID == subj &
    !is.na(meta$Visit)     & meta$Visit %in% c("Visit1", "Visit2") &
    rownames(meta) %in% shared_cols
  ]
  if (length(v12_samps) == 0) {
    v12_samps <- rownames(meta)[
      !is.na(meta$SubjectID) & meta$SubjectID == subj &
      rownames(meta) %in% shared_cols
    ]
  }
  if (length(v12_samps) == 1)
    log2(cpm[, v12_samps, drop = TRUE] + 1)
  else
    rowMeans(log2(cpm[, v12_samps, drop = FALSE] + 1))
})
names(avg_cpm_list_3d) <- subj_ids_3d

n_subj_3d <- length(subj_ids_3d)

# Per-subject colorRampPalette (Panel 1 method)
samp_pal_3d <- setNames(
  colorRampPalette(c("#4575b4","#74add1","#abd9e9","#fdae61","#f46d43","#d73027",
                     "#a50026","#313695","#fee090","#e0f3f8"))(n_subj_3d),
  subj_ids_3d
)

# Frequency polygon data (same helper + breaks as Panel 1)
subj_visit_map_3d <- setNames(rep("Averaged", n_subj_3d), subj_ids_3d)
poly_3d_df <- .make_freq_poly(avg_cpm_list_3d, seq(0, 16, by = 0.5), subj_visit_map_3d)

# Outlier detection: 1.5×IQR fence on per-subject median log2(CPM+1)
med_3d      <- sapply(avg_cpm_list_3d, median)
q1_3d       <- quantile(med_3d, 0.25)
q3_3d       <- quantile(med_3d, 0.75)
iqr_3d      <- IQR(med_3d)
auto_out_3d <- subj_ids_3d[med_3d < q1_3d - 1.5*iqr_3d | med_3d > q3_3d + 1.5*iqr_3d]

# PROMO07 always labeled
label_subjs_3d <- unique(c(auto_out_3d, "PROMO07"))
label_subjs_3d <- label_subjs_3d[label_subjs_3d %in% subj_ids_3d]

# Label placed at mode bin of each flagged subject's polygon
label_pts_3d <- do.call(rbind, lapply(label_subjs_3d, function(subj) {
  df_s        <- poly_3d_df[poly_3d_df$Sample == subj, ]
  row         <- df_s[which.max(df_s$count), , drop = FALSE]
  row$display <- subj
  row
}))

# CUTOFF vlines
cpm_vline_3d <- lapply(CUTOFF_LIST, function(t) {
  geom_vline(xintercept = log2(t + 1),
             color      = cutoff_colors[as.character(t)],
             linetype   = "dashed", linewidth = 0.85)
})

# Main frequency polygon panel
p_hist_3d <- ggplot(poly_3d_df,
                    aes(x = bin_mid, y = count, group = Sample, color = Sample)) +
  geom_line(alpha = 0.60, linewidth = 0.45) +
  cpm_vline_3d +
  geom_label_repel(
    data        = label_pts_3d,
    aes(x = bin_mid, y = count, label = display),
    inherit.aes = FALSE,
    size        = 2.5, box.padding = 0.40, max.overlaps = 50,
    segment.size = 0.30, segment.color = "grey50",
    fill = "white", alpha = 0.90, label.size = 0.2
  ) +
  scale_color_manual(values = samp_pal_3d, guide = "none") +
  scale_x_continuous(
    breaks = .log2_bin_breaks,
    labels = .log2_bin_labels,
    guide  = guide_axis(angle = 90)
  ) +
  labs(
    title    = paste0("Gene expression distribution — V1/V2-averaged per subject (",
                      n_subj_3d, " subjects, zeros excluded)"),
    subtitle = paste0("Per-gene mean log2(CPM+1) across Visit1 & Visit2 | ",
                      "bin width: 0.5 | outliers & PROMO07 labeled"),
    x = "log2(CPM+1) bin", y = "Gene count"
  ) +
  theme_cowplot(12)

# Panel-3-style swatch color key (exact Panel 3 code, adapted for subjects)
n_per_col_3d <- 20L
legend_3d_df <- data.frame(
  Sample  = names(samp_pal_3d),
  col_num = ceiling(seq_len(n_subj_3d) / n_per_col_3d),
  row_num = rev(((seq_len(n_subj_3d) - 1L) %% n_per_col_3d) + 1L),
  stringsAsFactors = FALSE
)
p_legend_3d <- ggplot(legend_3d_df, aes(x = col_num, y = row_num, color = Sample)) +
  geom_point(size = 3.5, shape = 15) +
  geom_text(aes(x = col_num + 0.12, label = Sample),
            hjust = 0, size = 2.8, color = "grey15") +
  scale_color_manual(values = samp_pal_3d, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
  theme_void() +
  labs(title = paste0("Subject color key (n = ", n_subj_3d, ")"))

n_key_cols_3d    <- max(legend_3d_df$col_num)
hist_3d_combined <- plot_grid(p_hist_3d, p_legend_3d,
                               ncol = 2,
                               rel_widths = c(3, 0.55 * n_key_cols_3d))

ggsave(file.path(DIR_QC, "00d_expr_hist_v1v2_avg_per_subject.png"),
       hist_3d_combined,
       width  = 16 + 3 * n_key_cols_3d,
       height = 10, dpi = 150, bg = "white")
message("  Saved: 00d_expr_hist_v1v2_avg_per_subject.png")

# ─────────────────────────────────────────────────────────────────────────────────
# 3e. ZOOMED OVERLAY HISTOGRAM: ALL SAMPLES, y = 0–4000, PROMO07 BOLDED
# Same overlaid frequency polygon as Section 3c. coord_cartesian clips at 4000
# without dropping data. PROMO07 samples drawn on top with a heavier line.
# ─────────────────────────────────────────────────────────────────────────────────
message("=== [3e] Zoomed overlay histogram — y 0-4000, PROMO07 bolded ===")

promo07_samp_ids <- rownames(meta)[
  !is.na(meta$SubjectID) & meta$SubjectID == "PROMO07" &
  rownames(meta) %in% all_samp_ids_3c
]

poly_3e_bg   <- poly_3c_df[!poly_3c_df$Sample %in% promo07_samp_ids, ]
poly_3e_bold <- poly_3c_df[ poly_3c_df$Sample %in% promo07_samp_ids, ]

p_hist_3e <- ggplot() +
  geom_line(data = poly_3e_bg,
            aes(x = bin_mid, y = count, group = Sample, color = Sample),
            alpha = 0.55, linewidth = 0.40) +
  cpm_vline_3c +
  geom_line(data = poly_3e_bold,
            aes(x = bin_mid, y = count, group = Sample, color = Sample),
            alpha = 1.00, linewidth = 1.30) +
  geom_label_repel(
    data        = label_pts_3c,
    aes(x = bin_mid, y = count, label = display),
    inherit.aes = FALSE,
    size        = 2.5, box.padding = 0.40, max.overlaps = 50,
    segment.size = 0.30, segment.color = "grey50",
    fill = "white", alpha = 0.90, label.size = 0.2
  ) +
  coord_cartesian(ylim = c(0, 4000)) +
  scale_color_manual(values = samp_pal_3c, guide = "none") +
  scale_x_continuous(
    breaks = .log2_bin_breaks,
    labels = .log2_bin_labels,
    guide  = guide_axis(angle = 90)
  ) +
  labs(
    title    = paste0("Gene expression distribution — all ", n_all_3c,
                      " samples, zoomed y ≤ 4000 (zeros excluded)"),
    subtitle = paste0("log2(CPM+1) per gene | bin width: 0.5 | ",
                      "PROMO07 bolded | outliers & PROMO26 V1/V2 labeled"),
    x = "log2(CPM+1) bin", y = "Gene count"
  ) +
  theme_cowplot(12)

hist_3e_combined <- plot_grid(p_hist_3e, p_legend_3c,
                               ncol = 2,
                               rel_widths = c(3, 0.55 * n_key_cols_3c))

ggsave(file.path(DIR_QC, "00e_expr_hist_zoom_promo07bold.png"),
       hist_3e_combined,
       width  = 16 + 3 * n_key_cols_3c,
       height = 10, dpi = 150, bg = "white")
message("  Saved: 00e_expr_hist_zoom_promo07bold.png")


# ─────────────────────────────────────────────────────────────────────────────
# 4. VST NORMALISATION (all samples, blind — for clustering only)
# ─────────────────────────────────────────────────────────────────────────────
message("=== [4] VST normalisation for clustering ===")

dds_full <- DESeqDataSetFromMatrix(
  countData = round(raw_filt),
  colData   = meta,
  design    = ~ 1
)
dds_full <- estimateSizeFactors(dds_full)
vst_full <- vst(dds_full, blind = TRUE)
vst_mat  <- assay(vst_full)


# ─────────────────────────────────────────────────────────────────────────────
# 5. VARIABLE GENE SELECTION — TWO PARALLEL APPROACHES
#
# Both approaches select PROMO-only samples and return the top N_VAR_GENES.
# Outputs (mean-variance plots, spot-checks, PCA, heatmaps, k-means) are
# generated for EACH approach with a suffix in the filename:
#   _VST   : variance computed on DESeq2 VST-normalised values
#   _rawCPM: variance computed on raw CPM values (before log transformation)
#
# VST inflates variance for lowly expressed genes because the dispersion model
# assigns relatively larger stabilised values to genes with few counts.
# rawCPM variance is the most direct measure: no transformation applied before
# computing spread. High-expression genes will tend to dominate since absolute
# variance scales with mean in count data — this is visible in the mean-variance
# plot and is precisely what VST is designed to correct.
# The two gene lists will partially overlap but will differ at the margins.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [5] Variable gene selection — VST and rawCPM approaches ===")

promo_cols_vst <- rownames(meta)[!meta$IsPerlman & rownames(meta) %in% colnames(vst_mat)]
promo_cols_cpm <- rownames(meta)[!meta$IsPerlman & rownames(meta) %in% colnames(cpm_filt)]

# ── VST approach ──────────────────────────────────────────────────────────────
gene_vars_vst        <- rowVars(vst_mat[, promo_cols_vst])
names(gene_vars_vst) <- rownames(vst_mat)
top_var_vst          <- names(sort(gene_vars_vst, decreasing = TRUE))[1:N_VAR_GENES]

# ── CPM coefficient-of-variation approach — CV = sd/mean per gene, all PROMO ──
# CV normalises variability by mean expression, so high-CV genes are variable
# relative to their own mean rather than in absolute terms.
cpm_promo_mat    <- as.matrix(cpm_filt[, promo_cols_cpm])
gene_vars_rowvar <- setNames(
  apply(cpm_promo_mat, 1, function(x) {
    m <- mean(x, na.rm = TRUE)
    if (m == 0) return(0)
    sd(x, na.rm = TRUE) / m
  }),
  rownames(cpm_filt)
)
top_var_rowvar <- names(sort(gene_vars_rowvar, decreasing = TRUE))[1:N_VAR_GENES]

# ── min-max LFC approach — log2(CPM+1) range across V1 PROMO samples ─────────
#
# For each gene: if expressed (log2(CPM+1) > 0) in >= MIN_PCT_SAMPLES of V1 PROMO
# samples, compute LFC = max(log2(CPM+1)) - min(log2(CPM+1)) across those samples.
# Genes with LFC >= 1 are selected as HVGs (1 log2 unit = 2-fold range).
# Uses ALL V1 PROMO samples (including PROMO02/04).
v1_promo_lfc     <- rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman &
                                    rownames(meta) %in% colnames(cpm_filt)]
log2cpm_v1_lfc   <- log2(as.matrix(cpm_filt[, v1_promo_lfc]) + 1)
v1_min_n_lfc     <- ceiling(MIN_PCT_SAMPLES * length(v1_promo_lfc))
expressed_v1_lfc <- rowSums(log2cpm_v1_lfc > 0) >= v1_min_n_lfc

lfc_vals        <- apply(log2cpm_v1_lfc[expressed_v1_lfc, , drop = FALSE], 1,
                          function(x) max(x) - min(x))
gene_vars_lfc            <- setNames(rep(0, nrow(cpm_filt)), rownames(cpm_filt))
gene_vars_lfc[names(lfc_vals)] <- lfc_vals

top_var_lfc_all <- names(gene_vars_lfc)[gene_vars_lfc >= 1]
top_var_lfc <- top_var_lfc_all
message("  minMaxLFC HVG: V1 PROMO samples used: ", length(v1_promo_lfc))
message("  minMaxLFC HVG: genes expressed in >= ", v1_min_n_lfc, " V1 samples: ",
        sum(expressed_v1_lfc))
message("  minMaxLFC HVG: genes with LFC >= 1: ", sum(gene_vars_lfc >= 1),
        " | selected for analysis: ", length(top_var_lfc))

# ── P10-P90 LFC approach — robust percentile range across V1 PROMO samples ──
#
# Identical expressed-gene filter as min-max LFC. LFC = P90 - P10 of
# log2(CPM+1) across V1 PROMO samples. Reduces influence of single outlier
# samples relative to max-min. Same LFC >= 1 threshold.
lfc_vals_pct <- apply(log2cpm_v1_lfc[expressed_v1_lfc, , drop = FALSE], 1,
                      function(x) quantile(x, 0.9, names = FALSE) -
                                  quantile(x, 0.1, names = FALSE))
gene_vars_lfc_pct <- setNames(rep(0, nrow(cpm_filt)), rownames(cpm_filt))
gene_vars_lfc_pct[names(lfc_vals_pct)] <- lfc_vals_pct

top_var_lfc_pct_all <- names(gene_vars_lfc_pct)[gene_vars_lfc_pct >= 1]
top_var_lfc_pct <- top_var_lfc_pct_all
message("  P10-P90 LFC HVG: genes with LFC >= 1: ", sum(gene_vars_lfc_pct >= 1),
        " | selected: ", length(top_var_lfc_pct))

# ── P5-P95 LFC approach — broader percentile range across V1 PROMO samples ──
#
# Wider window than P10-P90: captures signal from samples near the distribution
# tails while still trimming the single most extreme values.
# By construction: P10-P90 LFC <= P5-P95 LFC <= min-max LFC for any gene,
# so P10-P90 set ⊆ P5-P95 set ⊆ min-max set.
lfc_vals_p5p95 <- apply(log2cpm_v1_lfc[expressed_v1_lfc, , drop = FALSE], 1,
                        function(x) quantile(x, 0.95, names = FALSE) -
                                    quantile(x, 0.05, names = FALSE))
gene_vars_lfc_p5p95 <- setNames(rep(0, nrow(cpm_filt)), rownames(cpm_filt))
gene_vars_lfc_p5p95[names(lfc_vals_p5p95)] <- lfc_vals_p5p95

top_var_lfc_p5p95_all <- names(gene_vars_lfc_p5p95)[gene_vars_lfc_p5p95 >= 1]
top_var_lfc_p5p95     <- top_var_lfc_p5p95_all
message("  P5-P95 LFC HVG: genes with LFC >= 1: ", sum(gene_vars_lfc_p5p95 >= 1),
        " | selected: ", length(top_var_lfc_p5p95))

# Mean CPM for display (shared by all approaches)
gene_means_cpm <- rowMeans(cpm_filt[, promo_cols_cpm])

message("  VST approach — min variance in top-", N_VAR_GENES, ": ",
        round(min(gene_vars_vst[top_var_vst]), 4))
message("  CPM CV — min CV in top-", N_VAR_GENES, ": ",
        round(min(gene_vars_rowvar[top_var_rowvar]), 4))
message("  LFC HVG count (LFC>=1): ", sum(gene_vars_lfc >= 1),
        " | selected: ", length(top_var_lfc))
message("  P10-P90 LFC HVG count (LFC>=1): ", sum(gene_vars_lfc_pct >= 1),
        " | selected: ", length(top_var_lfc_pct))
message("  P5-P95 LFC HVG count (LFC>=1): ", sum(gene_vars_lfc_p5p95 >= 1),
        " | selected: ", length(top_var_lfc_p5p95))
message("  Genes in common VST ∩ rowVar:        ", length(intersect(top_var_vst, top_var_rowvar)))
message("  Genes in common VST ∩ LFC:           ", length(intersect(top_var_vst, top_var_lfc)))
message("  Genes in common VST ∩ LFC_PCT:       ", length(intersect(top_var_vst, top_var_lfc_pct)))
message("  Genes in common VST ∩ LFC_P5P95:     ", length(intersect(top_var_vst, top_var_lfc_p5p95)))
message("  Genes in common rowVar ∩ LFC:        ", length(intersect(top_var_rowvar, top_var_lfc)))
message("  Genes in common rowVar ∩ LFC_PCT:    ", length(intersect(top_var_rowvar, top_var_lfc_pct)))
message("  Genes in common rowVar ∩ LFC_P5P95:  ", length(intersect(top_var_rowvar, top_var_lfc_p5p95)))
message("  Genes in common LFC ∩ LFC_PCT:       ", length(intersect(top_var_lfc, top_var_lfc_pct)))
message("  Genes in common LFC ∩ LFC_P5P95:     ", length(intersect(top_var_lfc, top_var_lfc_p5p95)))
message("  Genes in common LFC_PCT ∩ LFC_P5P95: ", length(intersect(top_var_lfc_pct, top_var_lfc_p5p95)))

# Save combined variance table
var_result <- data.frame(
  ENSG            = names(gene_vars_vst),
  Gene_Symbol     = gene_sym_map$Gene_Symbol[match(names(gene_vars_vst), gene_sym_map$ENSG)],
  variance_vst    = gene_vars_vst,
  cv_cpm          = gene_vars_rowvar[names(gene_vars_vst)],
  lfc_minmax      = gene_vars_lfc[names(gene_vars_vst)],
  lfc_p5p95       = gene_vars_lfc_p5p95[names(gene_vars_vst)],
  lfc_p10p90      = gene_vars_lfc_pct[names(gene_vars_vst)],
  mean_cpm        = gene_means_cpm[names(gene_vars_vst)],
  top_vst         = names(gene_vars_vst) %in% top_var_vst,
  top_rowvar      = names(gene_vars_vst) %in% top_var_rowvar,
  top_lfc         = names(gene_vars_vst) %in% top_var_lfc,
  top_lfc_p5p95   = names(gene_vars_vst) %in% top_var_lfc_p5p95,
  top_lfc_pct     = names(gene_vars_vst) %in% top_var_lfc_pct
) %>% arrange(desc(variance_vst))
write.csv(var_result, file.path(DIR_HVG, "variance_variable_genes_all.csv"), row.names = FALSE)

# Export individual HVG lists for each method
write.csv(
  data.frame(ENSG = top_var_vst,
             Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_var_vst, gene_sym_map$ENSG)],
             variance_vst = gene_vars_vst[top_var_vst],
             stringsAsFactors = FALSE) %>% arrange(desc(variance_vst)),
  file.path(DIR_HVG, "HVG_list_VST.csv"), row.names = FALSE
)
write.csv(
  data.frame(ENSG = top_var_rowvar,
             Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_var_rowvar, gene_sym_map$ENSG)],
             cv_cpm = gene_vars_rowvar[top_var_rowvar],
             stringsAsFactors = FALSE) %>% arrange(desc(cv_cpm)),
  file.path(DIR_HVG, "HVG_list_CPM_CV.csv"), row.names = FALSE
)
write.csv(
  data.frame(ENSG = top_var_lfc,
             Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_var_lfc, gene_sym_map$ENSG)],
             lfc_minmax = gene_vars_lfc[top_var_lfc],
             stringsAsFactors = FALSE) %>% arrange(desc(lfc_minmax)),
  file.path(DIR_HVG, "HVG_list_minMaxLFC.csv"), row.names = FALSE
)
write.csv(
  data.frame(ENSG = top_var_lfc_pct,
             Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_var_lfc_pct, gene_sym_map$ENSG)],
             lfc_p10p90 = gene_vars_lfc_pct[top_var_lfc_pct],
             stringsAsFactors = FALSE) %>% arrange(desc(lfc_p10p90)),
  file.path(DIR_HVG, "HVG_list_P10P90LFC.csv"), row.names = FALSE
)
write.csv(
  data.frame(ENSG = top_var_lfc_p5p95,
             Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_var_lfc_p5p95, gene_sym_map$ENSG)],
             lfc_p5p95 = gene_vars_lfc_p5p95[top_var_lfc_p5p95],
             stringsAsFactors = FALSE) %>% arrange(desc(lfc_p5p95)),
  file.path(DIR_HVG, "HVG_list_P5P95LFC.csv"), row.names = FALSE
)
message("  Saved: variance_variable_genes_all.csv, HVG_list_VST.csv, ",
        "HVG_list_CPM_CV.csv, HVG_list_minMaxLFC.csv, HVG_list_P10P90LFC.csv, HVG_list_P5P95LFC.csv")

# ── Shared infrastructure: short label helper + sample annotation ─────────────
# Helper: shorten sample IDs for axis labels
#   PROMO13_V1_CD11b_S1  →  PROMO13_V1  |  PerlmanNorm102_...  →  Perlman102
# make.unique() handles PROMO24's 3 V1 tech reps: V1, V1.1, V1.2
short_label <- function(x) {
  labels <- ifelse(
    grepl("PerlmanNorm", x),
    paste0("Perlman", str_extract(x, "(?<=PerlmanNorm)\\d+")),
    str_extract(x, "PROMO\\d+_V\\d+")
  )
  make.unique(labels, sep = ".")
}

# min-max normalisation helper (CPM-based, shared by both approaches)
min_max_norm <- function(x, min_val = -1, max_val = 1) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  min_val + (x - rng[1]) / (rng[2] - rng[1]) * (max_val - min_val)
}

# Build sample annotation data frame (all metadata columns; gene-set-independent)
samp_ids_all <- colnames(vst_mat)
ann_col <- data.frame(
  Visit      = as.character(meta[samp_ids_all, "Visit"]),
  Day        = meta[samp_ids_all, "Day"],
  DayNum     = meta[samp_ids_all, "DayNum"],
  RunOrder   = meta[samp_ids_all, "RunOrder"],
  IsControl  = ifelse(meta[samp_ids_all, "IsPerlman"], "Perlman", "PROMO"),
  PPT_scaled = meta[samp_ids_all, "qst_ppt_tr_avg_v1"],
  PROMIS     = meta[samp_ids_all, "PROMIS"],
  Age        = meta[samp_ids_all, "Age"],
  Sex        = as.character(meta[samp_ids_all, "Sex"]),
  row.names  = samp_ids_all,
  stringsAsFactors = FALSE
)
ann_col$Sex[ann_col$Sex %in% c("NA", "")] <- NA

ann_colors <- list(
  Visit      = c(Visit1 = "#d84b4b", Visit2 = "#4bd84b", Visit3 = "#4b4bd8"),
  IsControl  = c(Perlman = "#E69F00", PROMO = "#AAAAAA"),
  PPT_scaled = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  PROMIS     = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  DayNum     = colorRampPalette(c("#ffffcc","#800026"))(100),
  Age        = colorRampPalette(c("#fff7ec","#7f2704"))(100),
  Sex        = c(Female = "#e78ac3", Male = "#66c2a5")
)
ann_colors_promo           <- ann_colors
ann_colors_promo$IsControl <- NULL

# ─────────────────────────────────────────────────────────────────────────────
# 5c. MEAN-VARIANCE COMPARISON GRID — VST | rowVar CPM | MAD²
#
# Three panels side-by-side to compare HVG selection methods:
#   Left  — VST variance: no y-axis log (shows raw stabilised variance)
#   Mid   — rowVar on raw CPM (all PROMO incl. 02/04): log2(var+1) on y
#   Right — MAD² on raw CPM (excl. 02/04): log2(var+1) on y
# All panels: x = log2(mean CPM) across PROMO samples.
# Red dots = top-N_VAR_GENES selected by each method.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [5c] Mean-variance comparison grid ===")

make_mv_panel <- function(gene_means, gene_vars, top_genes, x_label, y_label,
                           y_transform = identity, title_str, label_n = 20) {
  df <- data.frame(
    ENSG     = names(gene_vars),
    mean_cpm = gene_means[names(gene_vars)],
    variance = gene_vars,
    is_top   = names(gene_vars) %in% top_genes
  ) %>% filter(!is.na(variance), !is.na(mean_cpm), mean_cpm > 0, variance >= 0)
  df$y_val <- y_transform(df$variance)
  df <- df %>% filter(is.finite(y_val))

  top_labels <- df %>%
    filter(is_top) %>%
    arrange(desc(y_val)) %>%
    slice_head(n = label_n) %>%
    mutate(Gene_Symbol = gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)])

  ggplot(df, aes(x = log2(mean_cpm + 0.01), y = y_val,
                 color = is_top, alpha = is_top)) +
    geom_point(size = 1.8, stroke = 0) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(x = log2(mean_cpm + 0.01), y = y_val, label = Gene_Symbol),
      inherit.aes = FALSE,
      color = "#8B0000", size = 2.2, max.overlaps = 25,
      box.padding = 0.3, min.segment.length = 0
    ) +
    scale_color_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#d84b4b"),
                       labels = c("FALSE" = "Other", "TRUE" = "Top HVGs")) +
    scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 0.75), guide = "none") +
    labs(title = title_str,
         x = "log2(mean CPM)", y = y_label, color = NULL) +
    theme_cowplot(10) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          plot.title  = element_text(size = 10, face = "bold"))
}

p_mv_vst <- make_mv_panel(
  gene_means  = gene_means_cpm,
  gene_vars   = gene_vars_vst,
  top_genes   = top_var_vst,
  y_transform = identity,
  y_label     = "Variance (VST)",
  title_str   = "VST variance"
)

p_mv_rowvar <- make_mv_panel(
  gene_means  = gene_means_cpm,
  gene_vars   = gene_vars_rowvar,
  top_genes   = top_var_rowvar,
  y_transform = identity,
  y_label     = "CV (sd/mean CPM)",
  title_str   = "CPM Coefficient of Variation (all PROMO)"
)

p_mv_lfc <- make_mv_panel(
  gene_means  = gene_means_cpm,
  gene_vars   = gene_vars_lfc,
  top_genes   = top_var_lfc,
  y_transform = identity,
  y_label     = "min-max LFC (log2 range)",
  title_str   = "min-max LFC (V1 PROMO, LFC≥1)"
)

p_mv_lfc_pct <- make_mv_panel(
  gene_means  = gene_means_cpm,
  gene_vars   = gene_vars_lfc_pct,
  top_genes   = top_var_lfc_pct,
  y_transform = identity,
  y_label     = "P10-P90 LFC (log2 range)",
  title_str   = "P10-P90 LFC (V1 PROMO, LFC≥1)"
)

p_mv_lfc_p5p95 <- make_mv_panel(
  gene_means  = gene_means_cpm,
  gene_vars   = gene_vars_lfc_p5p95,
  top_genes   = top_var_lfc_p5p95,
  y_transform = identity,
  y_label     = "P5-P95 LFC (log2 range)",
  title_str   = "P5-P95 LFC (V1 PROMO, LFC≥1)"
)

mv_compare_grid <- plot_grid(p_mv_vst, p_mv_rowvar, p_mv_lfc, p_mv_lfc_p5p95, p_mv_lfc_pct,
                              ncol = 5, nrow = 1)
ggsave(file.path(DIR_HVG, "00_mean_variance_comparison_grid.png"),
       mv_compare_grid, width = 30, height = 7, dpi = 150, bg = "white")
message("  Saved: 00_mean_variance_comparison_grid.png")

# Focused comparator: all three LFC methods side by side
lfc_comparator_grid <- plot_grid(p_mv_lfc, p_mv_lfc_p5p95, p_mv_lfc_pct, ncol = 3, nrow = 1)
ggsave(file.path(DIR_HVG, "00_mean_variance_LFC_comparator.png"),
       lfc_comparator_grid, width = 21, height = 7, dpi = 150, bg = "white")
message("  Saved: 00_mean_variance_LFC_comparator.png")
message("  Overlap VST    \u2229 rowVar:      ", length(intersect(top_var_vst, top_var_rowvar)))
message("  Overlap VST    \u2229 LFC:         ", length(intersect(top_var_vst, top_var_lfc)))
message("  Overlap VST    \u2229 LFC_PCT:     ", length(intersect(top_var_vst, top_var_lfc_pct)))
message("  Overlap VST    \u2229 LFC_P5P95:   ", length(intersect(top_var_vst, top_var_lfc_p5p95)))
message("  Overlap rowVar \u2229 LFC:         ", length(intersect(top_var_rowvar, top_var_lfc)))
message("  Overlap LFC    \u2229 LFC_PCT:     ", length(intersect(top_var_lfc, top_var_lfc_pct)))
message("  Overlap LFC    \u2229 LFC_P5P95:   ", length(intersect(top_var_lfc, top_var_lfc_p5p95)))
message("  Overlap LFC_PCT \u2229 LFC_P5P95:  ", length(intersect(top_var_lfc_pct, top_var_lfc_p5p95)))

# ─────────────────────────────────────────────────────────────────────────────
# 5d. VENN DIAGRAM — min-max LFC vs P10-P90 LFC HVGs
#
# Compares the two LFC-based HVG sets and spot-checks representative genes from
# each region of the Venn diagram using CPM vs PPT scatter plots.
# ─────────────────────────────────────────────────────────────────────────────
# [5d] Venn diagram deferred — generated after section 6F so both LFC gene lists
# use the PROMO02/04-excluded cohort, matching the k-means gene universe.

# Spot-check scatter plots: top N genes by mean CPM from each Venn region
.venn_spotcheck <- function(gene_ids, n_show = 9) {
  v1_ids_sc <- rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman &
                               !is.na(meta$qst_ppt_raw) &
                               rownames(meta) %in% colnames(cpm_filt)]
  valid <- gene_ids[gene_ids %in% rownames(cpm_filt)]
  if (length(valid) == 0) return(NULL)
  mean_cpm_set <- rowMeans(cpm_filt[valid, v1_ids_sc, drop = FALSE])
  top_ids <- names(sort(mean_cpm_set, decreasing = TRUE))[seq_len(min(n_show, length(mean_cpm_set)))]

  plot_list <- lapply(top_ids, function(g) {
    sym     <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    if (length(sym) == 0 || is.na(sym)) sym <- g
    lfc_mm  <- round(gene_vars_lfc[g],     2)
    lfc_pct <- round(gene_vars_lfc_pct[g], 2)
    df <- data.frame(
      PPT = meta[v1_ids_sc, "qst_ppt_raw"],
      CPM = as.numeric(cpm_filt[g, v1_ids_sc]),
      stringsAsFactors = FALSE
    )
    ggplot(df, aes(PPT, CPM)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      labs(title    = sym,
           subtitle = paste0("mmLFC=", lfc_mm, " | p10p90=", lfc_pct),
           x = "PPT", y = "CPM") +
      theme_cowplot(9) +
      theme(plot.title    = element_text(face = "bold", size = 9),
            plot.subtitle = element_text(size = 7))
  })
  plot_grid(plotlist = plot_list, ncol = 3)
}

# Spotcheck calls deferred to after section 6F (same reason as Venn above).


# ─────────────────────────────────────────────────────────────────────────────
# 6B. SAMPLE–SAMPLE CORRELATION HEATMAP (gene-set-independent, runs once)
#
# Input: ALL unfiltered genes, log2(CPM+1), PROMO samples only.
# This is a global QC view of transcriptional similarity between samples
# and does not depend on which 2000 genes were selected.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [6B] Sample-sample Pearson correlation — pre- and post-filter ===")

promo_ids_6b   <- colnames(cpm)[colnames(cpm) %in% rownames(meta) &
                                 !meta[colnames(cpm), "IsPerlman"]]
promo_ids_6b   <- promo_ids_6b[!is.na(promo_ids_6b)]
promo_short_6b <- short_label(promo_ids_6b)

# Shared annotation — all metadata cols except IsControl (all are PROMO here)
ann_promo <- ann_col[promo_ids_6b, , drop = FALSE]
rownames(ann_promo) <- promo_short_6b
ann_promo$IsControl <- NULL

# Helper: build and save one correlation heatmap + gene-contribution CSV
save_cor_heatmap <- function(expr_mat, sample_short, ann_df, title_str, filename,
                              save_contrib = TRUE, out_dir = OUT_DIR) {
  cor_m   <- cor(expr_mat, method = "pearson")
  hc_m    <- hclust(as.dist(1 - cor_m), method = "complete")
  colnames(cor_m) <- sample_short
  rownames(cor_m) <- sample_short
  hc_m$labels     <- sample_short
  cor_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  png(file.path(out_dir, filename), width = 6900, height = 6700, res = 300, type = "cairo")
  pheatmap(cor_m,
           cluster_rows = hc_m, cluster_cols = hc_m,
           annotation_col = ann_df,
           annotation_colors = ann_colors_promo,
           color = cor_colors,
           show_rownames = FALSE, show_colnames = TRUE,
           fontsize_col = 12, fontsize_row = 22, fontsize = 15,
           main = title_str,
           legend = TRUE, annotation_legend = TRUE, border_color = NA)
  dev.off()
  message("  Saved: ", filename)

  if (save_contrib) {
    # Per-gene CPM variance across samples — identifies genes driving
    # sample-to-sample correlation differences. CPM is already normalised so
    # no additional scaling is applied.
    gene_cor_var <- rowVars(expr_mat)
    names(gene_cor_var) <- rownames(expr_mat)
    gene_cor_var <- sort(gene_cor_var, decreasing = TRUE)

    top_driver_genes <- head(names(gene_cor_var), 500)
    driver_df <- data.frame(
      ENSG        = top_driver_genes,
      Gene_Symbol = gene_sym_map$Gene_Symbol[match(top_driver_genes, gene_sym_map$ENSG)],
      cpm_variance = gene_cor_var[top_driver_genes],
      stringsAsFactors = FALSE
    )
    csv_name <- sub("\\.png$", "_gene_drivers.csv", filename)
    write.csv(driver_df, file.path(out_dir, csv_name), row.names = FALSE)
    message("  Saved: ", csv_name)

    # Per-sample mean absolute deviation from median correlation — flags outlier samples
    cor_no_diag       <- cor_m
    diag(cor_no_diag) <- NA
    samp_mean_cor     <- rowMeans(cor_no_diag, na.rm = TRUE)
    samp_cor_df <- data.frame(
      Sample      = names(samp_mean_cor),
      mean_pearson = round(samp_mean_cor, 4),
      stringsAsFactors = FALSE
    ) %>% arrange(mean_pearson)
    csv_samp <- sub("\\.png$", "_sample_mean_cor.csv", filename)
    write.csv(samp_cor_df, file.path(out_dir, csv_samp), row.names = FALSE)
    message("  Saved: ", csv_samp)

    invisible(list(cor_m = cor_m, gene_drivers = driver_df,
                   sample_mean_cor = samp_cor_df))
  }
}

# ── Pre-filter: all genes, raw CPM ───────────────────────────────────────────
save_cor_heatmap(
  expr_mat     = as.matrix(cpm[, promo_ids_6b]),
  sample_short = promo_short_6b,
  ann_df       = ann_promo,
  title_str    = paste0("PROMO samples | Pearson correlation | ",
                        nrow(cpm), " genes (unfiltered) raw CPM | Complete HC"),
  filename     = "01b_correlation_PROMO_prefilter.png",
  out_dir      = DIR_QC
)

# ── Post-filter: filtered genes, raw CPM ─────────────────────────────────────
promo_ids_filt   <- intersect(promo_ids_6b, colnames(cpm_filt))
promo_short_filt <- short_label(promo_ids_filt)

cor_heatmap_result <- save_cor_heatmap(
  expr_mat     = as.matrix(cpm_filt[, promo_ids_filt]),
  sample_short = promo_short_filt,
  ann_df       = ann_promo[promo_short_filt, , drop = FALSE],
  title_str    = paste0("PROMO samples | Pearson correlation | ",
                        nrow(cpm_filt), " genes (filtered) raw CPM | Complete HC"),
  filename     = "01b_correlation_PROMO_postfilter.png",
  save_contrib = TRUE,
  out_dir      = DIR_QC
)

# ── Gene-driver dotplot: top 40 genes by z-score variance ────────────────────
if (!is.null(cor_heatmap_result) && is.list(cor_heatmap_result)) {
  driver_top <- cor_heatmap_result$gene_drivers %>%
    slice_head(n = 40) %>%
    arrange(cpm_variance) %>%
    mutate(Gene_Symbol = factor(Gene_Symbol, levels = Gene_Symbol))

  p_driver <- ggplot(driver_top, aes(x = cpm_variance, y = Gene_Symbol, fill = cpm_variance)) +
    geom_col() +
    scale_fill_gradient(low = "#AAAAAA", high = "#B2182B") +
    labs(title    = "Top genes driving sample-sample correlation variance",
         subtitle = paste0("Raw CPM variance across ", length(promo_ids_filt),
                           " PROMO samples | higher = gene varies more across samples"),
         x = "CPM variance", y = NULL, fill = "CPM var") +
    theme_cowplot(11) +
    theme(axis.text.y = element_text(size = 8), legend.position = "none")

  ggsave(file.path(DIR_QC, "01b_correlation_gene_drivers.png"),
         p_driver, width = 8, height = 9, dpi = 150, bg = "white")
  message("  Saved: 01b_correlation_gene_drivers.png")

  # Sample-level mean correlation bar chart (shows outlier samples)
  samp_cor <- cor_heatmap_result$sample_mean_cor %>%
    arrange(mean_pearson) %>%
    mutate(Sample = factor(Sample, levels = Sample))

  p_samp_cor <- ggplot(samp_cor, aes(x = mean_pearson, y = Sample,
                                      fill = mean_pearson)) +
    geom_col() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = median(samp_cor$mean_pearson)) +
    geom_vline(xintercept = median(samp_cor$mean_pearson),
               linetype = "dashed", color = "grey40") +
    labs(title    = "Mean Pearson correlation per sample vs. all other PROMO samples",
         subtitle = "Low values = transcriptionally distinct outliers",
         x = "Mean Pearson r", y = NULL) +
    theme_cowplot(9) +
    theme(axis.text.y = element_text(size = 6), legend.position = "none")

  ggsave(file.path(DIR_QC, "01b_correlation_sample_mean_cor.png"),
         p_samp_cor, width = 8, height = max(4, length(promo_ids_filt) * 0.14),
         dpi = 150, bg = "white")
  message("  Saved: 01b_correlation_sample_mean_cor.png")
}

# ─────────────────────────────────────────────────────────────────────────────
# 6C-ii. PATHWAY SETS + HELPER FUNCTIONS (used by 6E, section 9, section 11.5)
# ─────────────────────────────────────────────────────────────────────────────
message("=== [6C-ii] Loading pathway sets ===")

# Load pathway sets here so they are available to downstream sections too
hallmark_sets   <- msigdbr(species = "Homo sapiens", category = "H")
pathways_h      <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

reactome_sets   <- msigdbr(species = "Homo sapiens", category = "C2",
                            subcategory = "CP:REACTOME")
pathways_reactome <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)

gobp_sets       <- msigdbr(species = "Homo sapiens", category = "C5",
                            subcategory = "GO:BP")
pathways_gobp   <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)

message("  Pathways loaded — Hallmark: ", length(pathways_h),
        " | REACTOME: ", length(pathways_reactome),
        " | GO:BP: ", length(pathways_gobp))

save_gsea_tbl <- function(gsea_res, filename, out_dir = OUT_DIR, leading_edge = FALSE) {
  df <- as.data.frame(gsea_res)
  if (leading_edge && "leadingEdge" %in% names(gsea_res)) {
    df$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ", ")
  } else {
    df <- df[, !vapply(df, is.list, logical(1)), drop = FALSE]
  }
  write.csv(df %>% arrange(pval), file.path(out_dir, filename), row.names = FALSE)
  message("  Saved: ", file.path(basename(out_dir), filename))
  invisible(df)
}
# Generic dotplot helper (handles any pathway prefix; used here and in section 11)
plot_gsea_dotplot_generic <- function(gsea_tbl, top_n = 25, pval_col = "pval",
                                       sig_cutoff = 0.05, title = "GSEA",
                                       strip_prefix = NULL) {
  if (is.null(gsea_tbl) || nrow(gsea_tbl) == 0) return(NULL)
  df <- as.data.frame(gsea_tbl)
  df <- df[, !vapply(df, is.list, logical(1)), drop = FALSE]
  sig <- df %>% filter(.data[[pval_col]] < sig_cutoff)
  if (nrow(sig) == 0) {
    message("    No significant pathways at p<", sig_cutoff,
            " — showing top ", top_n, " by pval")
    sig <- df %>% arrange(.data[[pval_col]]) %>% slice_head(n = top_n)
  }
  plot_df <- sig %>%
    mutate(abs_NES = abs(NES)) %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = top_n) %>%
    mutate(
      neglog10p     = -log10(.data[[pval_col]]),
      pathway_label = if (!is.null(strip_prefix))
                        gsub(strip_prefix, "", pathway, perl = TRUE)
                      else pathway,
      pathway_label = gsub("_", " ", pathway_label),
      pathway_label = factor(pathway_label, levels = pathway_label[order(NES)])
    )
  ggplot(plot_df, aes(x = NES, y = pathway_label, size = neglog10p, color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    scale_size_continuous(range = c(3, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title    = title,
         subtitle = "NES: positive = enriched in outliers | size = -log10(p)",
         x = "NES", y = NULL, size = "-log10(p)", color = "NES") +
    theme_cowplot(12) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.y   = element_text(size = 9),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3))
}

n_load <- 25
make_loading_bar <- function(df, pc_col, abs_col, pc_label) {
  top_df <- df %>% arrange(desc(.data[[abs_col]])) %>% slice_head(n = n_load) %>%
    arrange(.data[[pc_col]]) %>%
    mutate(Gene_Symbol = factor(Gene_Symbol, levels = Gene_Symbol))
  ggplot(top_df, aes(.data[[pc_col]], Gene_Symbol, fill = .data[[pc_col]] > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                      guide = "none") +
    labs(title = paste0("Top ", n_load, " ", pc_label, " loadings"),
         x = paste(pc_label, "loading"), y = NULL) +
    theme_cowplot(9) +
    theme(axis.text.y = element_text(size = 7))
}

# Per-sample gene contribution bar: shows which genes drive an outlier's PC position.
# contrib_vec: named numeric vector (names = ENSG); contribution = centered_expr * loading.
make_contrib_bar <- function(contrib_vec, pc_label, sample_id, pc_score, n_top = 25) {
  n_show  <- min(n_top, length(contrib_vec))
  top_idx <- order(abs(contrib_vec), decreasing = TRUE)[seq_len(n_show)]
  top_c   <- contrib_vec[top_idx]
  syms    <- gene_sym_map$Gene_Symbol[match(names(top_c), gene_sym_map$ENSG)]
  syms    <- ifelse(is.na(syms), names(top_c), syms)
  df <- data.frame(
    Gene_Symbol  = factor(syms, levels = syms[order(top_c)]),
    Contribution = as.numeric(top_c),
    stringsAsFactors = FALSE
  )
  ggplot(df, aes(x = Contribution, y = Gene_Symbol, fill = Contribution > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
    labs(title    = paste0(sample_id, " — ", pc_label, " gene contributions"),
         subtitle = paste0(pc_label, " score = ", round(pc_score, 2),
                           " | top ", n_show, " genes by |contribution|"),
         x = paste0(pc_label, " contribution (centred × loading)"), y = NULL) +
    theme_cowplot(9) +
    theme(axis.text.y = element_text(size = 7))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6E. PROMO36 / PROMO39 — CORRELATION SCATTER, OUTLIER GSEA, PCA DRIVERS
#
# PROMO36_V1, PROMO39_V1, and PROMO39_V2 stand out in the sample–sample
# Pearson correlation matrix (Section 6B).  Reference mean is restricted to
# Visit1 and Visit2 PROMO samples only (no Visit3, no Perlman).
# Analyses mirror Sections 6C, 6C-ii, 6D, and 6D-ii.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [6E] PROMO36/PROMO39 — correlation scatter, GSEA, PCA drivers ===")

OUTLIER_SHORT_6E <- c("PROMO36_V1", "PROMO39_V1", "PROMO39_V2")
promo_all_6e     <- promo_ids_6b[
  !is.na(meta[promo_ids_6b, "Visit"]) &
  meta[promo_ids_6b, "Visit"] %in% c("Visit1", "Visit2")
]
outlier_ids_6e   <- promo_all_6e[
  str_extract(promo_all_6e, "PROMO\\d+_V\\d+") %in% OUTLIER_SHORT_6E
]
rest_ids_6e      <- setdiff(promo_all_6e, outlier_ids_6e)

if (length(outlier_ids_6e) == 0) {
  message("  WARNING: none of ", paste(OUTLIER_SHORT_6E, collapse = ", "),
          " found in promo_ids_6b — skipping Section 6E")
} else {
  message("  Outlier samples (n = ", length(outlier_ids_6e), "): ",
          paste(str_extract(outlier_ids_6e, "PROMO\\d+_V\\d+"), collapse = ", "))
  message("  Reference samples (n = ", length(rest_ids_6e), ")")

  log2cpm_all_6e <- log2(as.matrix(cpm_filt[, promo_all_6e]) + 1)
  mean_rest_6e   <- rowMeans(log2cpm_all_6e[, rest_ids_6e, drop = FALSE])

  # ── 6E. Scatter: each outlier vs. mean of all other PROMO samples ───────────
  scatter_list_6e <- lapply(outlier_ids_6e, function(oid) {
    short_id <- str_extract(oid, "PROMO\\d+_V\\d+")
    out_expr <- log2cpm_all_6e[, oid]
    scat_df  <- data.frame(
      gene      = names(out_expr),
      outlier   = as.numeric(out_expr),
      mean_rest = mean_rest_6e[names(out_expr)],
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(outlier), !is.na(mean_rest)) %>%
      mutate(
        Gene_Symbol = gene_sym_map$Gene_Symbol[match(gene, gene_sym_map$ENSG)],
        residual    = outlier - mean_rest
      )

    # Write scatter data for this outlier sample
    write.csv(scat_df,
              file.path(DIR_OUTLIER, paste0("scatter_data_", short_id, "_vs_rest.csv")),
              row.names = FALSE)
    message("  Saved: scatter_data_", short_id, "_vs_rest.csv")

    r_val  <- round(cor(scat_df$outlier, scat_df$mean_rest, method = "pearson"), 3)
    lab_df <- scat_df %>%
      filter(outlier > 11) %>%
      mutate(dir = ifelse(residual >= 0, "up", "dn"))

    ggplot(scat_df, aes(mean_rest, outlier)) +
      geom_point(alpha = 0.25, size = 0.5, color = "#AAAAAA") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "#555555", linewidth = 0.7) +
      geom_text_repel(data = lab_df,
                      aes(label = Gene_Symbol, color = dir),
                      size = 2.3, max.overlaps = Inf, box.padding = 0.3,
                      show.legend = FALSE) +
      scale_color_manual(values = c(up = "#B2182B", dn = "#2166AC")) +
      labs(title    = paste(short_id, "vs. mean of other PROMO samples"),
           subtitle = paste0("r = ", r_val,
                             " | log2(CPM+1) | identity line dashed\n",
                             "Labeled: outlier log2(CPM+1) > 11 | ",
                             "Red = higher in ", short_id,
                             " | Blue = lower in ", short_id),
           x = "Mean log2(CPM+1) — other PROMO",
           y = paste0("log2(CPM+1) — ", short_id)) +
      theme_cowplot(11)
  })

  scatter_combined_6e <- plot_grid(plotlist = scatter_list_6e,
                                    ncol = length(scatter_list_6e))
  ggsave(file.path(DIR_OUTLIER, "00h_outlier_scatter_PROMO36_PROMO39.png"),
         scatter_combined_6e,
         width = 7 * max(length(scatter_list_6e), 1), height = 7,
         dpi = 150, bg = "white")
  message("  Saved: 00h_outlier_scatter_PROMO36_PROMO39.png")

  # ── 6E-ii. GSEA from pooled residuals ───────────────────────────────────────
  message("=== [6E-ii] Outlier GSEA — PROMO36/PROMO39 vs. rest ===")

  resid_list_6e <- lapply(outlier_ids_6e, function(oid) {
    resid <- as.numeric(log2cpm_all_6e[, oid]) - mean_rest_6e
    setNames(resid, rownames(log2cpm_all_6e))
  })
  pooled_resid_6e <- Reduce("+", resid_list_6e) / length(resid_list_6e)

  resid_sym_6e <- setNames(
    pooled_resid_6e,
    gene_sym_map$Gene_Symbol[match(names(pooled_resid_6e), gene_sym_map$ENSG)]
  )
  resid_sym_6e <- resid_sym_6e[!is.na(names(resid_sym_6e)) & !duplicated(names(resid_sym_6e))]
  resid_sym_6e <- sort(resid_sym_6e, decreasing = TRUE)
  message("  Ranked gene list: ", length(resid_sym_6e), " genes")

  gsea_6e_h  <- fgsea(pathways_h,        stats = resid_sym_6e, minSize = 15, maxSize = 500)
  gsea_6e_r  <- fgsea(pathways_reactome, stats = resid_sym_6e, minSize = 15, maxSize = 500)
  gsea_6e_bp <- fgsea(pathways_gobp,     stats = resid_sym_6e, minSize = 15, maxSize = 500)

  save_gsea_tbl(gsea_6e_h,  "GSEA_outlier_PROMO36_PROMO39_Hallmark.csv", out_dir = DIR_OUTLIER)
  save_gsea_tbl(gsea_6e_r,  "GSEA_outlier_PROMO36_PROMO39_Reactome.csv", out_dir = DIR_OUTLIER)
  save_gsea_tbl(gsea_6e_bp, "GSEA_outlier_PROMO36_PROMO39_GOBP.csv",     out_dir = DIR_OUTLIER)

  dp_6e_h  <- plot_gsea_dotplot_generic(gsea_6e_h,  top_n = 25,
    title = "Outlier GSEA (PROMO36/39 vs. PROMO rest) — Hallmark",
    strip_prefix = "^HALLMARK_")
  dp_6e_r  <- plot_gsea_dotplot_generic(gsea_6e_r,  top_n = 30,
    title = "Outlier GSEA (PROMO36/39 vs. PROMO rest) — REACTOME",
    strip_prefix = "^REACTOME_")
  dp_6e_bp <- plot_gsea_dotplot_generic(gsea_6e_bp, top_n = 30,
    title = "Outlier GSEA (PROMO36/39 vs. PROMO rest) — GO:BP",
    strip_prefix = "^GOBP_")

  for (dp_pair in list(
    list(dp = dp_6e_h,  fn = "00h2_GSEA_PROMO36_PROMO39_Hallmark_dotplot.png",  w = 10, h = 7),
    list(dp = dp_6e_r,  fn = "00h2_GSEA_PROMO36_PROMO39_Reactome_dotplot.png",  w = 14, h = 12),
    list(dp = dp_6e_bp, fn = "00h2_GSEA_PROMO36_PROMO39_GOBP_dotplot.png",      w = 14, h = 12)
  )) {
    if (!is.null(dp_pair$dp)) {
      ggsave(file.path(DIR_OUTLIER, dp_pair$fn), dp_pair$dp,
             width = dp_pair$w, height = dp_pair$h, dpi = 150, bg = "white")
      message("  Saved: ", dp_pair$fn)
    }
  }

  # Per-gene residual table: each outlier's log2(CPM+1) minus reference mean
  resid_df_6e <- data.frame(
    ENSG        = rownames(log2cpm_all_6e),
    Gene_Symbol = gene_sym_map$Gene_Symbol[match(rownames(log2cpm_all_6e),
                                                   gene_sym_map$ENSG)],
    mean_rest_log2cpm = round(mean_rest_6e, 4),
    stringsAsFactors = FALSE
  )
  for (oid in outlier_ids_6e) {
    short_id <- str_extract(oid, "PROMO\\d+_V\\d+")
    resid_col <- as.numeric(log2cpm_all_6e[, oid]) - mean_rest_6e
    resid_df_6e[[paste0("resid_", short_id)]] <- round(resid_col, 4)
  }
  resid_df_6e$pooled_mean_resid <- round(pooled_resid_6e, 4)
  resid_df_6e <- resid_df_6e %>% arrange(desc(abs(pooled_mean_resid)))
  write.csv(resid_df_6e,
            file.path(DIR_OUTLIER, "outlier_PROMO36_PROMO39_gene_residuals.csv"),
            row.names = FALSE)
  message("  Saved: outlier_PROMO36_PROMO39_gene_residuals.csv")

  # ── 6E-iii. PCA — all PROMO samples, all visits, all filtered genes ──────────
  message("=== [6E-iii] PCA — all PROMO, all visits — PROMO36/39 drivers ===")

  pca_6e     <- prcomp(t(log2cpm_all_6e), center = TRUE, scale. = FALSE)
  pct_var_6e <- round(100 * pca_6e$sdev^2 / sum(pca_6e$sdev^2), 1)

  is_outlier_6e <- str_extract(rownames(pca_6e$x), "PROMO\\d+_V\\d+") %in% OUTLIER_SHORT_6E
  pca_6e_df <- data.frame(
    PC1        = pca_6e$x[, 1],
    PC2        = pca_6e$x[, 2],
    ShortID    = short_label(rownames(pca_6e$x)),
    Visit      = as.character(meta[rownames(pca_6e$x), "Visit"]),
    Sex        = as.character(meta[rownames(pca_6e$x), "Sex"]),
    is_outlier = is_outlier_6e,
    stringsAsFactors = FALSE
  )
  pca_6e_df$Visit <- ifelse(is.na(pca_6e_df$Visit), "Unknown", pca_6e_df$Visit)
  pca_6e_df$Sex   <- ifelse(is.na(pca_6e_df$Sex) | pca_6e_df$Sex == "NA",
                              "Unknown", pca_6e_df$Sex)
  pca_6e_df$label <- ifelse(pca_6e_df$is_outlier, pca_6e_df$ShortID, NA)

  p_pca_6e <- ggplot(pca_6e_df, aes(PC1, PC2, color = is_outlier, shape = Visit)) +
    geom_point(aes(size = is_outlier), alpha = 0.85) +
    geom_text_repel(aes(label = label), na.rm = TRUE, size = 3,
                    box.padding = 0.4, max.overlaps = 20, color = "#8B0000") +
    scale_color_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#d84b4b"),
                       labels = c("FALSE" = "Other",
                                  "TRUE"  = paste(OUTLIER_SHORT_6E, collapse = "/"))) +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    scale_shape_manual(values = c(Visit1 = 16, Visit2 = 15, Visit3 = 17, Unknown = 8)) +
    labs(title    = "PCA — all filtered genes | all PROMO samples",
         subtitle = paste0(nrow(log2cpm_all_6e), " genes | log2(CPM+1) centered | ",
                           "all visits"),
         x = paste0("PC1 (", pct_var_6e[1], "%)"),
         y = paste0("PC2 (", pct_var_6e[2], "%)"),
         color = NULL, shape = "Visit") +
    theme_cowplot(12)

  load_df_6e <- data.frame(
    ENSG        = rownames(pca_6e$rotation),
    Gene_Symbol = gene_sym_map$Gene_Symbol[match(rownames(pca_6e$rotation),
                                                   gene_sym_map$ENSG)],
    PC1 = pca_6e$rotation[, 1],
    PC2 = pca_6e$rotation[, 2]
  ) %>%
    mutate(PC1_abs = abs(PC1), PC2_abs = abs(PC2))

  p_pc1_bar_6e <- make_loading_bar(load_df_6e, "PC1", "PC1_abs", "PC1")
  p_pc2_bar_6e <- make_loading_bar(load_df_6e, "PC2", "PC2_abs", "PC2")

  pca_6e_grid <- plot_grid(p_pca_6e, p_pc1_bar_6e, p_pc2_bar_6e,
                             ncol = 3, rel_widths = c(1.5, 1, 1))
  ggsave(file.path(DIR_OUTLIER, "00h3_PCA_allgenes_PROMO36_PROMO39_drivers.png"),
         pca_6e_grid, width = 19, height = 7, dpi = 150, bg = "white")
  message("  Saved: 00h3_PCA_allgenes_PROMO36_PROMO39_drivers.png")

  write.csv(load_df_6e %>% arrange(desc(PC1_abs)),
            file.path(DIR_OUTLIER, "PCA_PROMO36_PROMO39_PC1_PC2_loadings.csv"),
            row.names = FALSE)
  message("  Saved: PCA_PROMO36_PROMO39_PC1_PC2_loadings.csv")

  # Per-gene PC score contributions for each outlier sample
  centered_mat_6e <- sweep(t(log2cpm_all_6e), 2, pca_6e$center, "-")

  contrib_plot_list_6e <- lapply(outlier_ids_6e, function(oid) {
    short_id    <- str_extract(oid, "PROMO\\d+_V\\d+")
    contrib_pc1 <- setNames(centered_mat_6e[oid, ] * pca_6e$rotation[, 1],
                             colnames(centered_mat_6e))
    contrib_pc2 <- setNames(centered_mat_6e[oid, ] * pca_6e$rotation[, 2],
                             colnames(centered_mat_6e))
    list(
      pc1 = make_contrib_bar(contrib_pc1, "PC1", short_id, pca_6e$x[oid, 1]),
      pc2 = make_contrib_bar(contrib_pc2, "PC2", short_id, pca_6e$x[oid, 2])
    )
  })

  contrib_panels_6e <- unlist(lapply(contrib_plot_list_6e, function(x) list(x$pc1, x$pc2)),
                               recursive = FALSE)
  contrib_grid_6e <- plot_grid(plotlist = contrib_panels_6e,
                                ncol = 2, nrow = length(outlier_ids_6e))
  ggsave(file.path(DIR_OUTLIER, "00h4_PCA_PROMO36_PROMO39_gene_contributions.png"),
         contrib_grid_6e,
         width  = 12,
         height = 5.5 * length(outlier_ids_6e),
         dpi = 150, bg = "white")
  message("  Saved: 00h4_PCA_PROMO36_PROMO39_gene_contributions.png")

  contrib_rows_6e <- do.call(rbind, lapply(outlier_ids_6e, function(oid) {
    short_id <- str_extract(oid, "PROMO\\d+_V\\d+")
    c1 <- setNames(centered_mat_6e[oid, ] * pca_6e$rotation[, 1], colnames(centered_mat_6e))
    c2 <- setNames(centered_mat_6e[oid, ] * pca_6e$rotation[, 2], colnames(centered_mat_6e))
    data.frame(
      Sample      = short_id,
      ENSG        = names(c1),
      Gene_Symbol = gene_sym_map$Gene_Symbol[match(names(c1), gene_sym_map$ENSG)],
      contrib_PC1 = as.numeric(c1),
      contrib_PC2 = as.numeric(c2),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(contrib_rows_6e %>% arrange(Sample, desc(abs(contrib_PC1))),
            file.path(DIR_OUTLIER, "PCA_PROMO36_PROMO39_gene_contributions.csv"),
            row.names = FALSE)
  message("  Saved: PCA_PROMO36_PROMO39_gene_contributions.csv")

  # ── 6E-iv. Heatmap of top 100 driver genes ───────────────────────────────────
  message("  Generating heatmap of top driver genes")

  top_driver_genes_6e <- head(names(sort(abs(pooled_resid_6e), decreasing = TRUE)), 100)
  valid_rows_6e       <- top_driver_genes_6e[top_driver_genes_6e %in% rownames(cpm_filt)]

  if (length(valid_rows_6e) > 0) {
    # Min-max normalize [-1, 1] per row on raw CPM values (no log2)
    cpm_raw_6e  <- as.matrix(cpm_filt[valid_rows_6e, promo_all_6e])
    hm_mat_6e   <- t(apply(cpm_raw_6e, 1, min_max_norm))
    row_syms_6e <- gene_sym_map$Gene_Symbol[match(rownames(hm_mat_6e), gene_sym_map$ENSG)]
    rownames(hm_mat_6e) <- ifelse(is.na(row_syms_6e), rownames(hm_mat_6e), row_syms_6e)

    is_out_col_6e <- str_extract(colnames(hm_mat_6e), "PROMO\\d+_V\\d+") %in% OUTLIER_SHORT_6E

    # Place PROMO36/39 columns on the right (next to row labels); cluster non-outliers only
    non_out_cols <- colnames(hm_mat_6e)[!is_out_col_6e]
    out_cols     <- colnames(hm_mat_6e)[ is_out_col_6e]
    if (length(non_out_cols) > 1) {
      d_non_out   <- as.dist(1 - cor(hm_mat_6e[, non_out_cols, drop = FALSE],
                                     use = "pairwise.complete.obs"))
      hc_non_out  <- hclust(d_non_out, method = "complete")
      non_out_cols <- non_out_cols[hc_non_out$order]
    }
    hm_mat_6e_ord <- hm_mat_6e[, c(non_out_cols, out_cols), drop = FALSE]

    ann_hm_6e <- data.frame(
      Outlier = ifelse(colnames(hm_mat_6e_ord) %in% out_cols, "PROMO36/39", "Other"),
      Visit   = as.character(meta[colnames(hm_mat_6e_ord), "Visit"]),
      row.names = colnames(hm_mat_6e_ord)
    )
    ann_hm_colors_6e <- list(
      Outlier = c("PROMO36/39" = "#d84b4b", "Other" = "#AAAAAA"),
      Visit   = c(Visit1 = "#d84b4b", Visit2 = "#4bd84b", Visit3 = "#4b4bd8")
    )

    png(file.path(DIR_OUTLIER, "00h5_heatmap_top_driver_genes_PROMO36_PROMO39.png"),
        width = 4800, height = 5000, res = 200, type = "cairo")
    pheatmap(hm_mat_6e_ord,
             annotation_col    = ann_hm_6e,
             annotation_colors = ann_hm_colors_6e,
             scale             = "none",
             cluster_cols      = FALSE,
             clustering_distance_rows = "correlation",
             show_colnames = FALSE,
             show_rownames = TRUE,
             fontsize_row  = 22,
             border_color  = NA,
             main = paste0("Top 100 driver genes — PROMO36/39 outliers\n",
                           "raw CPM min-max [-1,1] per row | V1+V2 PROMO samples"))
    dev.off()
    message("  Saved: 00h5_heatmap_top_driver_genes_PROMO36_PROMO39.png")
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# 5d. VENN DIAGRAM — min-max LFC vs P5-P95 LFC vs P10-P90 LFC HVGs (3-way)
#
# By construction these sets are nested: P10-P90 ⊆ P5-P95 ⊆ min-max.
# The Venn cleanly separates:
#   all-three   = bulk-distribution variable (all three methods agree)
#   mm+P5P95    = P5-P95 captures but P10-P90 doesn't (moderate-outlier-driven)
#   mm-only     = single-outlier-driven (max/min pulled by one extreme sample)
# The two smaller sets (P10-P90-only, P5-P95-only) are always empty.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [5d] Venn diagram — 3-way LFC HVG comparison (clean cohort) ===")

# All-pairwise and triple overlaps
overlap_all3      <- Reduce(intersect, list(top_var_lfc, top_var_lfc_p5p95, top_var_lfc_pct))
overlap_mm_p5p95  <- setdiff(intersect(top_var_lfc, top_var_lfc_p5p95), top_var_lfc_pct)
lfc_only_hvg      <- setdiff(top_var_lfc, union(top_var_lfc_p5p95, top_var_lfc_pct))
overlap_hvg       <- intersect(top_var_lfc, top_var_lfc_pct)  # kept for downstream spotcheck

message("  All-3 overlap (min-max ∩ P5-P95 ∩ P10-P90): ", length(overlap_all3))
message("  min-max ∩ P5-P95 only (not P10-P90):        ", length(overlap_mm_p5p95))
message("  min-max only (not P5-P95 or P10-P90):       ", length(lfc_only_hvg),
        " — outlier-driven")
message("  (P5-P95-only and P10-P90-only are always 0 by nested set property)")

venn_list_clean <- list(
  "min-max LFC" = top_var_lfc,
  "P5-P95 LFC"  = top_var_lfc_p5p95,
  "P10-P90 LFC" = top_var_lfc_pct
)
p_venn <- ggVennDiagram(venn_list_clean, label_alpha = 0) +
  scale_fill_gradient(low = "#EFF3FF", high = "#2166AC") +
  scale_color_manual(values = c("#d84b4b", "#E08214", "#2166AC")) +
  labs(title    = "HVG overlap — min-max LFC vs P5-P95 LFC vs P10-P90 LFC",
       subtitle = paste0("min-max: n=", length(top_var_lfc),
                         " | P5-P95: n=", length(top_var_lfc_p5p95),
                         " | P10-P90: n=", length(top_var_lfc_pct),
                         " | PROMO02/04 excluded at load")) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"))
ggsave(file.path(DIR_HVG, "00_venn_LFC_3way.png"),
       p_venn, width = 8, height = 7, dpi = 150, bg = "white")
message("  Saved: 00_venn_LFC_3way.png")

venn_groups <- list(
  list(ids = overlap_all3,     fn = "00_spotcheck_venn_all3.png",       desc = "all-3 overlap"),
  list(ids = overlap_mm_p5p95, fn = "00_spotcheck_venn_mm_p5p95.png",   desc = "min-max+P5-P95 only"),
  list(ids = lfc_only_hvg,     fn = "00_spotcheck_venn_lfc_only.png",   desc = "min-max only")
)
for (vg in venn_groups) {
  g <- .venn_spotcheck(vg$ids)
  if (!is.null(g)) {
    ggsave(file.path(DIR_HVG, vg$fn), g,
           width = 12, height = 12, dpi = 150, bg = "white")
    message("  Saved: ", vg$fn, " (", vg$desc, ")")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 5e. IFNAR SCATTER — V1 PROMO raw CPM vs PPT (if IFNAR1/IFNAR2 in any HVG list)
#
# Scans every HVG set (VST, CPM CV, min-max LFC, P5-P95 LFC, P10-P90 LFC).
# For any gene whose symbol matches ^IFNAR (i.e. IFNAR1, IFNAR2) and that
# appears in at least one HVG list, emit a raw CPM vs PPT scatter for V1 PROMO
# samples with a linear fit + 95% CI. Subtitle lists which HVG sets selected
# the gene, so the connection to the upstream HVG selection is explicit.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [5e] IFNAR scatter (PPT vs raw CPM) — V1 PROMO ===")

ifnar_hvg_sets <- list(
  VST       = top_var_vst,
  CPM_CV    = top_var_rowvar,
  minMaxLFC = top_var_lfc,
  P5P95LFC  = top_var_lfc_p5p95,
  P10P90LFC = top_var_lfc_pct
)

ifnar_ensg_all <- gene_sym_map$ENSG[grepl("^IFNAR", gene_sym_map$Gene_Symbol)]
ifnar_in_any   <- character(0)
ifnar_membership <- list()  # ENSG -> char vector of HVG set names that include it
for (eid in ifnar_ensg_all) {
  hits <- names(ifnar_hvg_sets)[vapply(ifnar_hvg_sets,
                                        function(s) eid %in% s, logical(1))]
  if (length(hits) > 0) {
    ifnar_in_any <- c(ifnar_in_any, eid)
    ifnar_membership[[eid]] <- hits
  }
}

if (length(ifnar_in_any) == 0) {
  message("  No IFNAR* gene found in any HVG list — skipping scatter")
} else {
  v1_ids_ifnar <- rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman &
                                  !is.na(meta$qst_ppt_raw) &
                                  rownames(meta) %in% colnames(cpm_filt)]
  message("  IFNAR genes in HVG lists (n = ", length(ifnar_in_any), "): ",
          paste(gene_sym_map$Gene_Symbol[match(ifnar_in_any, gene_sym_map$ENSG)],
                collapse = ", "))

  for (eid in ifnar_in_any) {
    sym  <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == eid][1]
    sets <- ifnar_membership[[eid]]
    df_ifnar <- data.frame(
      PPT     = meta[v1_ids_ifnar, "qst_ppt_raw"],
      CPM     = as.numeric(cpm_filt[eid, v1_ids_ifnar]),
      ShortID = short_label(v1_ids_ifnar),
      stringsAsFactors = FALSE
    )
    fit  <- lm(CPM ~ PPT, data = df_ifnar)
    r_val <- round(cor(df_ifnar$PPT, df_ifnar$CPM,
                       use = "pairwise.complete.obs"), 3)
    p_val <- signif(summary(fit)$coefficients[2, 4], 3)

    p_ifnar <- ggplot(df_ifnar, aes(PPT, CPM)) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      geom_point(color = "#2c7bb6", alpha = 0.85, size = 2.2) +
      geom_text_repel(aes(label = ShortID), size = 2.4,
                      box.padding = 0.3, max.overlaps = 25,
                      color = "grey25") +
      labs(
        title    = paste0(sym, " — raw CPM vs PPT (V1 PROMO)"),
        subtitle = paste0("HVG sets: ", paste(sets, collapse = ", "),
                          " | n = ", nrow(df_ifnar),
                          " | Pearson r = ", r_val,
                          " | slope p = ", p_val),
        x = "PPT (raw qst_ppt_tr_avg_v1)", y = "Raw CPM"
      ) +
      theme_cowplot(12)

    fn_ifnar <- paste0("00_IFNAR_scatter_", sym, "_V1_PROMO.png")
    ggsave(file.path(DIR_HVG, fn_ifnar),
           p_ifnar, width = 7, height = 5.5, dpi = 150, bg = "white")
    message("  Saved: ", fn_ifnar)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 5b / 6A / 7.  GENE-SET LOOP — runs once for VST, once for LFC
#
# For each approach:
#   (i)   Mean-variance plot with top-20 labels
#   (ii)  Spot-check: top-16 variable genes as CPM scatter vs PPT
#   (iii) PCA of PROMO samples + per-sample variance contribution (SSD)
#   (iv)  Gene × sample heatmap (CPM min-max, HC column order)
#   (v)   K-means heatmaps (V1 PROMO only, k = 4:10)
#
# Display is always CPM min-max [-1,1].
# Clustering distances use the approach-specific scaled expression matrix:
#   VST approach    → z-scored VST values
#   rawCPM approach → z-scored log2(CPM+1) for clustering distance
#                     (raw CPM too skewed for Pearson correlation)
# ─────────────────────────────────────────────────────────────────────────────
gene_sets <- list(
  VST = list(
    top_genes  = top_var_vst,
    gene_vars  = gene_vars_vst,
    label      = "VST",
    var_ylab   = "Variance (VST)",
    clust_expr = vst_mat[, promo_cols_vst]
  ),
  LFC = list(
    top_genes  = top_var_lfc,
    gene_vars  = gene_vars_lfc,
    label      = "LFC",
    var_ylab   = "min-max LFC (log2 range, V1 PROMO)",
    # log2(CPM+1) for clustering distance — LFC method is V1-based but
    # we cluster across all PROMO for the full heatmap view
    clust_expr = log2(cpm_filt[, promo_cols_cpm] + 1)
  ),
  LFC_P5P95 = list(
    top_genes  = top_var_lfc_p5p95,
    gene_vars  = gene_vars_lfc_p5p95,
    label      = "LFC_P5P95",
    var_ylab   = "P5-P95 LFC (log2 range, V1 PROMO)",
    clust_expr = log2(cpm_filt[, promo_cols_cpm] + 1)
  ),
  LFC_PCT = list(
    top_genes  = top_var_lfc_pct,
    gene_vars  = gene_vars_lfc_pct,
    label      = "LFC_PCT",
    var_ylab   = "P10-P90 LFC (log2 range, V1 PROMO)",
    clust_expr = log2(cpm_filt[, promo_cols_cpm] + 1)
  )
)

v1_promo_ids <- rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman &
                                !is.na(meta$qst_ppt_raw)]

sil_summary_list <- list()  # accumulates per-k silhouette stats across all gene sets
# Accumulates k-means cluster assignments across all HVG methods and clustering modes.
# Used by section 7e (Jaccard z-vs-raw comparison) after the gene_sets loop.
# Key format: "<HVG method>__<file_tag>__k<k>"; value: named int vector (ENSG -> cluster).
kmeans_clusters_all <- list()

for (gs_name in names(gene_sets)) {
  gs            <- gene_sets[[gs_name]]
  top_var_genes <- gs$top_genes
  gene_vars_gs  <- gs$gene_vars
  suf           <- gs$label   # file suffix: "VST" or "logCPM"
  message(sprintf("\n--- Gene set: %s ---", suf))

  # ── CPM min-max normalised matrix for display ──────────────────────────────
  cpm_var  <- cpm_filt[top_var_genes, ]
  norm_mat <- t(apply(cpm_var, 1, min_max_norm))
  colnames(norm_mat) <- colnames(cpm_var)

  # ── Scaled expression for clustering distance (approach-specific) ──────────
  expr_clust <- gs$clust_expr[top_var_genes[top_var_genes %in% rownames(gs$clust_expr)], ]
  scaled_clust <- t(scale(t(expr_clust)))
  scaled_clust <- scaled_clust[complete.cases(scaled_clust), ]

  # ── (i) Mean-variance plot ─────────────────────────────────────────────────
  mv_df <- data.frame(
    ENSG     = names(gene_means_cpm),
    mean_cpm = gene_means_cpm,
    variance = gene_vars_gs[names(gene_means_cpm)],
    is_top   = names(gene_means_cpm) %in% top_var_genes,
    row.names = names(gene_means_cpm)
  ) %>% filter(!is.na(variance), mean_cpm > 0, variance > 0) %>%
    mutate(log2_var = log2(variance))

  # Top-20 by variance within HVG set — used for labeling in the loop pass.
  # A second, post-GSEA pass (see section 11b) replaces these labels with
  # leading-edge genes for the rawCPM approach.
  mv_top20 <- mv_df %>%
    filter(is_top) %>% arrange(desc(variance)) %>% slice_head(n = 20) %>%
    mutate(Gene_Symbol = gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)])

  mv_plot <- ggplot(mv_df, aes(x = log2(mean_cpm + 0.01), y = log2_var,
                                color = is_top, alpha = is_top)) +
    geom_point(size = 0.5) +
    geom_text_repel(data = mv_top20,
                    aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
                    color = "#8B0000", size = 2.5, max.overlaps = 20,
                    box.padding = 0.3, inherit.aes = FALSE) +
    scale_color_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#d84b4b"),
                       labels = c("FALSE" = "Other", "TRUE" = if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT"))
                         paste0("HVGs (LFC≥1, n=", length(top_var_genes), ")")
                         else paste0("Top ", N_VAR_GENES))) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
    labs(title = paste("Mean-variance relationship —", suf, "approach"),
         subtitle = if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT"))
           paste0("Red = all HVGs (LFC≥1, n=", length(top_var_genes),
                  ") | filter: CPM>", CPM_THRESH, " | top-20 labelled")
         else
           paste0("Red = top-", N_VAR_GENES,
                  " | filter: CPM>", CPM_THRESH, " | top-20 labelled"),
         x = "log2(mean CPM)", y = paste0("log2(", gs$var_ylab, ")"), color = NULL) +
    theme_cowplot(12) + theme(legend.position = "top")

  ggsave(file.path(DIR_CLUST, paste0("00_mean_variance_", suf, ".png")),
         mv_plot, width = 7, height = 5, dpi = 150, bg = "white")
  message("  Saved: 00_mean_variance_", suf, ".png")

  # Persist mv_df for post-GSEA enhanced labeling
  if (suf == "LFC")       mv_df_lfc       <- mv_df
  if (suf == "LFC_P5P95") mv_df_lfc_p5p95 <- mv_df
  if (suf == "LFC_PCT")   mv_df_lfc_pct   <- mv_df

  # ── (ii) Spot-check: top 16 genes, CPM vs PPT ─────────────────────────────
  spot_genes <- head(top_var_genes, 16)
  spot_list  <- lapply(spot_genes, function(g) {
    sym <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    vr  <- round(gene_vars_gs[g], 3)
    df  <- data.frame(
      PPT     = meta[v1_promo_ids, "qst_ppt_raw"],
      CPM     = as.numeric(cpm_filt[g, v1_promo_ids]),
      ShortID = short_label(v1_promo_ids),
      stringsAsFactors = FALSE
    )
    # Label top-3 and bottom-3 by CPM; PROMO24 replicates always labeled in purple
    # (sanity check: technical replicates should cluster tightly regardless of gene)
    n_label <- 3
    df_label <- rbind(
      df %>% arrange(desc(CPM)) %>% slice_head(n = n_label) %>% mutate(end = "high"),
      df %>% arrange(CPM)       %>% slice_head(n = n_label) %>% mutate(end = "low")
    )
    df_promo24 <- df %>% filter(grepl("PROMO24", ShortID)) %>% mutate(end = "promo24")
    df_label   <- bind_rows(df_label, df_promo24) %>%
      distinct(ShortID, .keep_all = TRUE)

    ggplot(df, aes(PPT, CPM)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_point(data = df %>% filter(grepl("PROMO24", ShortID)),
                 color = "#7B2D8B", size = 2.8, alpha = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      geom_text_repel(data = df_label,
                      aes(label = ShortID,
                          color = end),
                      size = 2, max.overlaps = 20,
                      box.padding = 0.25, show.legend = FALSE) +
      scale_color_manual(values = c(high = "#B2182B", low = "#2166AC",
                                    promo24 = "#7B2D8B")) +
      labs(title = sym, subtitle = paste0("var=", vr),
           x = "PPT", y = "CPM") +
      theme_cowplot(9) +
      theme(plot.title    = element_text(face = "bold", size = 9),
            plot.subtitle = element_text(size = 7))
  })
  png(file.path(DIR_CLUST, paste0("00b_spotcheck_", suf, ".png")),
      width = 4*400, height = 4*380, res = 120, type = "cairo")
  print(plot_grid(plotlist = spot_list, ncol = 4))
  dev.off()
  message("  Saved: 00b_spotcheck_", suf, ".png")

  # ── (iii) PCA + per-sample SSD ────────────────────────────────────────────
  promo_ids_pca <- intersect(promo_cols_vst, colnames(scaled_clust))
  sc_promo      <- scaled_clust[, promo_ids_pca]

  pca_res <- prcomp(t(sc_promo), center = FALSE, scale. = FALSE)
  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  pca_df <- data.frame(
    PC1     = pca_res$x[,1], PC2 = pca_res$x[,2],
    PPT     = meta[rownames(pca_res$x), "qst_ppt_raw"],
    Visit   = as.character(meta[rownames(pca_res$x), "Visit"]),
    ShortID = short_label(rownames(pca_res$x)),
    stringsAsFactors = FALSE
  )
  pca_df$dist_centre <- sqrt(scale(pca_df$PC1)^2 + scale(pca_df$PC2)^2)
  pca_df$is_outlier  <- pca_df$dist_centre > 3

  pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = PPT, label = ShortID)) +
    geom_point(aes(shape = Visit), size = 3, alpha = 0.85) +
    geom_text_repel(data = subset(pca_df, is_outlier),
                    size = 3, max.overlaps = 20, box.padding = 0.4,
                    color = "black") +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = median(pca_df$PPT, na.rm = TRUE),
                          na.value = "grey70") +
    labs(title = paste("PCA of PROMO samples —", suf, "gene selection"),
         subtitle = if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT"))
           paste0(length(top_var_genes), " HVGs (LFC≥1) | black labels = >3 SD from centroid")
         else
           paste0("Top-", N_VAR_GENES, " variable genes | black labels = >3 SD from centroid"),
         x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)"),
         color = "PPT", shape = "Visit") +
    theme_cowplot(12)
  ggsave(file.path(DIR_CLUST, paste0("00c_PCA_", suf, ".png")),
         pca_plot, width = 8, height = 6, dpi = 150, bg = "white")
  message("  Saved: 00c_PCA_", suf, ".png")

  # Per-sample SSD
  gene_means_sc  <- rowMeans(sc_promo)
  ssd_per_sample <- colSums((sc_promo - gene_means_sc)^2)
  ssd_df <- data.frame(
    ShortID = short_label(names(ssd_per_sample)),
    SSD     = ssd_per_sample,
    PPT     = meta[names(ssd_per_sample), "qst_ppt_raw"],
    stringsAsFactors = FALSE
  ) %>% arrange(desc(SSD))
  ssd_thresh     <- mean(ssd_df$SSD) + 2 * sd(ssd_df$SSD)
  ssd_df$is_high <- ssd_df$SSD > ssd_thresh
  ssd_df$ShortID <- factor(ssd_df$ShortID, levels = ssd_df$ShortID)

  ssd_plot <- ggplot(ssd_df, aes(ShortID, SSD, fill = is_high)) +
    geom_col() +
    geom_hline(yintercept = ssd_thresh, linetype = "dashed", color = "navy") +
    scale_fill_manual(values = c("FALSE"="#AAAAAA","TRUE"="#d84b4b"),
                      labels = c("FALSE"="Normal","TRUE"=">mean+2SD"), name = NULL) +
    labs(title = paste("Per-sample variance contribution —", suf),
         subtitle = if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT"))
           paste0("SSD from gene-wise mean | ", length(top_var_genes), " HVGs (LFC≥1) | dashed = mean+2SD")
         else
           paste0("SSD from gene-wise mean | top-", N_VAR_GENES, " genes | dashed = mean+2SD"),
         x = NULL, y = "Sum of squared deviations") +
    theme_cowplot(11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "top")
  ggsave(file.path(DIR_CLUST, paste0("00d_sample_variance_contribution_", suf, ".png")),
         ssd_plot, width = 14, height = 5, dpi = 150, bg = "white")
  message("  Saved: 00d_sample_variance_contribution_", suf, ".png")
  message("  High-leverage samples: ",
          paste(as.character(ssd_df$ShortID[ssd_df$is_high]), collapse = ", "))

  # ── (iv) 6A: Gene × sample heatmap ────────────────────────────────────────
  samp_ids  <- colnames(scaled_clust)
  samp_short <- short_label(samp_ids)
  ann_col_disp <- ann_col[samp_ids, , drop = FALSE]
  rownames(ann_col_disp) <- samp_short

  # Column order: correlation HC on scaled expression
  cor_samp  <- cor(scaled_clust, method = "pearson")
  hc_samp   <- hclust(as.dist(1 - cor_samp), method = "complete")
  hc_samp$labels <- samp_short
  hc_col_order   <- hc_samp$order

  norm_display   <- norm_mat[, samp_ids[hc_col_order]]
  colnames(norm_display) <- samp_short[hc_col_order]
  ann_col_hc     <- ann_col_disp[samp_short[hc_col_order], , drop = FALSE]

  # Row clustering on scaled expression (correlation-based)
  sc_var    <- scaled_clust[rownames(norm_display)[rownames(norm_display) %in% rownames(scaled_clust)], ]
  cor_genes <- cor(t(sc_var), method = "pearson")
  hc_genes  <- hclust(as.dist(1 - cor_genes), method = "complete")

  png(file.path(DIR_CLUST, paste0("01a_heatmap_genes_x_samples_", suf, ".png")),
      width = 2200, height = 1800, res = 150, type = "cairo")
  pheatmap(norm_display,
           cluster_rows = hc_genes, cluster_cols = FALSE,
           annotation_col = ann_col_hc,
           annotation_colors = ann_colors,
           color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           show_rownames = FALSE, show_colnames = TRUE,
           fontsize_col = 8,
           main = if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT")) {
             sprintf("%d HVGs (%s) | CPM min-max | HC cols", length(top_var_genes), suf)
           } else {
             sprintf("Top-%d variable genes (%s) | CPM min-max | HC cols", N_VAR_GENES, suf)
           },
           fontsize = 8, legend = TRUE, annotation_legend = TRUE, border_color = NA)
  dev.off()
  message("  Saved: 01a_heatmap_genes_x_samples_", suf, ".png")

  # ── (v) 7: K-means on V1 PROMO samples ────────────────────────────────────
  v1_ids    <- intersect(rownames(meta)[meta$Visit == "Visit1" & !meta$IsPerlman],
                         colnames(norm_mat))
  # Caveat A fix: subset to V1 PROMO samples FIRST, then min-max each gene over
  # ONLY those samples. The matrix k-means clusters must be scaled to span
  # [-1,1] across the exact samples being clustered — not inherited from a
  # min-max computed over the full (all-visit + Perlman) sample set.
  cpm_v1_raw <- as.matrix(cpm_filt[top_var_genes, v1_ids])
  norm_v1    <- t(apply(cpm_v1_raw, 1, min_max_norm))
  colnames(norm_v1) <- v1_ids
  v1_ids_km  <- v1_ids      # PROMO02/04 already dropped in section 6F
  norm_v1_km <- norm_v1

  # log2(CPM+1) min-max counterpart — same V1-subset-first scaling, but on the
  # log-transformed values. Used only to RE-DISPLAY the raw-CPM k-means result
  # (heatmap 02g, gene PCA 02h) in log space, so each cluster can be traced
  # before vs after the log2 transform. Clustering itself stays on raw CPM.
  log2_v1_raw   <- log2(as.matrix(cpm_filt[top_var_genes, v1_ids]) + 1)
  norm_v1_log2  <- t(apply(log2_v1_raw, 1, min_max_norm))
  colnames(norm_v1_log2) <- v1_ids

  # Column dendrogram uses all remaining V1 PROMO samples
  cor_v1    <- cor(norm_v1, method = "pearson")
  hc_v1     <- hclust(as.dist(1 - cor_v1), method = "complete")
  hc_v1$labels <- short_label(v1_ids)

  # Annotation for V1 PROMO columns
  ann_v1 <- ann_col[v1_ids, , drop = FALSE]
  rownames(ann_v1) <- short_label(v1_ids)
  ann_v1$IsControl <- NULL  # all are PROMO in V1 k-means

  # Pre-compute min-max CPM for ALL expressed genes in the k-means sample set
  # (used for per-cluster GSEA ranking)
  cpm_v1_km_all <- as.matrix(cpm_filt[, v1_ids_km])
  all_genes_norm_km <- t(apply(cpm_v1_km_all, 1, min_max_norm))

  # ── Z-score matrix for the parallel clustering pass ──────────────────────
  # Raw CPM on V1 PROMO samples (cpm_v1_raw, built above) then gene-wise
  # center + scale: each row has mean = 0, sd = 1 across the V1 PROMO cohort.
  # No log2 transform — scale() is applied directly to raw CPM.
  # Genes with sd = 0 (constant expression) produce 0/0 = NaN — drop them and
  # warn, so the matrix passed to kmeans() is always finite.
  z_v1_full <- t(scale(t(cpm_v1_raw), center = TRUE, scale = TRUE))
  z_finite  <- is.finite(rowSums(z_v1_full))
  if (any(!z_finite))
    message("  Z-score: dropped ", sum(!z_finite),
            " gene(s) with sd = 0 across V1 PROMO (", suf, ")")
  z_v1 <- z_v1_full[z_finite, , drop = FALSE]
  colnames(z_v1) <- v1_ids

  # Heatmap display transform for the z-score pass:
  #   raw CPM  →  gene-wise z-score  →  per-gene min-max [-1, 1].
  # The min-max rescale keeps the divergent colormap centered consistently
  # with the raw-CPM-min-max pass while preserving the z-score's gene-wise
  # standardization (so a peak gene has the same colormap saturation
  # regardless of the gene's absolute CPM magnitude).
  z_minmax_v1 <- t(apply(z_v1, 1, min_max_norm))

  # ── Clustering methods (parallel passes) ──────────────────────────────────
  # Each entry produces its own full set of artifacts (heatmap, PCAs,
  # silhouette, violin, ORA, leading-edge boxplots). file_tag is appended to
  # every output filename so the two passes never collide.
  #
  # Heatmap display (one per method, written as 02_heatmap_kmeans_*):
  #   minmax_raw  → raw CPM min-max [-1, 1]
  #   zscore_raw  → z-score → min-max [-1, 1]
  # The PCA-alt panels (02h) still use the log2(CPM+1) min-max transform as
  # a separate diagnostic; only the k-means HEATMAPS were switched off log2.
  norm_v1_log2_z <- norm_v1_log2[rownames(z_v1), , drop = FALSE]
  cluster_methods <- list(
    list(
      file_tag     = "",
      method_label = "raw CPM min-max",
      cluster_mat  = norm_v1_km,
      sil_dist_mat = norm_v1_km,
      pca_main     = list(mat = norm_v1_km,    label = "row min-max"),
      pca_alt      = list(mat = norm_v1_log2,  label = "log2 min-max"),
      displays = list(
        list(suf2 = "",  values = norm_v1,
             breaks = seq(-1, 1, length.out = 101),
             transform_label = "raw CPM min-max [-1,1]")
      )
    ),
    list(
      file_tag     = "_ZSCORE",
      method_label = "z-score raw CPM",
      cluster_mat  = z_v1,
      sil_dist_mat = z_v1,
      pca_main     = list(mat = z_v1,           label = "z-score"),
      pca_alt      = list(mat = norm_v1_log2_z, label = "log2 min-max"),
      displays = list(
        list(suf2 = "",  values = z_minmax_v1,
             breaks = seq(-1, 1, length.out = 101),
             transform_label = "z-score → min-max [-1,1]")
      )
    )
  )

  kmeans_results <- list()
  for (k in K_VALUES) {
    for (cm in cluster_methods) {
      method_tag <- cm$file_tag
      message("  k=", k, " (", suf, method_tag, " — ", cm$method_label, ")")
      set.seed(40)
      km <- kmeans(cm$cluster_mat, centers = k, nstart = 25, iter.max = 100)
      kmeans_results[[paste0("k", k, method_tag)]] <- km
      kmeans_clusters_all[[paste0(suf, "__", method_tag, "__k", k)]] <- km$cluster

      row_order   <- order(km$cluster)
      cluster_pal <- setNames(
        colorRampPalette(brewer.pal(max(k, 8), "Set2"))(k),
        paste0("C", 1:k)
      )

      # Each method now writes a single 02_heatmap_kmeans_*.png on the
      # transform declared in its displays list.
      for (disp_idx in seq_along(cm$displays)) {
        disp        <- cm$displays[[disp_idx]]
        mat_display <- disp$values[row_order, , drop = FALSE]
        colnames(mat_display) <- short_label(v1_ids)
        row_ann <- data.frame(
          GeneCluster = factor(paste0("C", km$cluster[row_order]),
                               levels = paste0("C", 1:k)),
          row.names   = rownames(mat_display)
        )
        hm_title <- if (suf %in% c("LFC", "LFC_P5P95", "LFC_PCT")) {
          sprintf("%d HVGs (%s%s) | k=%d | %s | V1 PROMO",
                  length(top_var_genes), suf, method_tag, k, disp$transform_label)
        } else {
          sprintf("Top-%d genes (%s%s) | k=%d | %s | V1 PROMO",
                  N_VAR_GENES, suf, method_tag, k, disp$transform_label)
        }
        hm_fn <- sprintf("02_heatmap_kmeans_k%d_%s%s.png", k, suf, method_tag)
        png(file.path(DIR_CLUST, hm_fn),
            width = 2400, height = 2400, res = 150, type = "cairo")
        pheatmap(mat_display,
                 cluster_rows = FALSE, cluster_cols = hc_v1,
                 annotation_col = ann_v1,
                 annotation_row = row_ann,
                 annotation_colors = c(ann_colors_promo, list(GeneCluster = cluster_pal)),
                 color  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
                 breaks = disp$breaks,
                 show_rownames = FALSE, show_colnames = TRUE,
                 fontsize_col = 8,
                 main = hm_title,
                 fontsize = 8, border_color = NA, annotation_legend = TRUE)
        dev.off()
        message("  Saved: ", hm_fn)
      }

      # ── Gene-space PCA on the clustering matrix (02f) ──────────────────────
      # Each point is ONE GENE projected by PCA of the exact matrix k-means
      # clustered (cm$pca_main$mat). Smeared colors over one cloud = no real
      # cluster structure.
      pca_genes  <- prcomp(cm$pca_main$mat, center = TRUE, scale. = FALSE)
      pca_g_pct  <- round(100 * pca_genes$sdev^2 / sum(pca_genes$sdev^2), 1)
      pca_g_df   <- data.frame(
        PC1     = pca_genes$x[, 1],
        PC2     = pca_genes$x[, 2],
        Cluster = factor(paste0("C", km$cluster[rownames(pca_genes$x)]),
                         levels = paste0("C", 1:k)),
        stringsAsFactors = FALSE
      )
      p_pca_genes <- ggplot(pca_g_df, aes(PC1, PC2, color = Cluster)) +
        geom_point(size = 0.8, alpha = 0.6) +
        scale_color_manual(values = cluster_pal, name = "k-means\ncluster") +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        labs(
          title    = sprintf("Gene PCA colored by k-means cluster — k=%d | %s%s | %s",
                             k, suf, method_tag, cm$method_label),
          subtitle = paste0(nrow(pca_g_df), " genes | PCA of ", cm$pca_main$label,
                            " matrix used for clustering\n",
                            "Smeared colors over one cloud = no real cluster structure (poor fit)"),
          x = paste0("PC1 (", pca_g_pct[1], "%)"),
          y = paste0("PC2 (", pca_g_pct[2], "%)")
        ) +
        theme_cowplot(12)
      ggsave(file.path(DIR_CLUST,
                       sprintf("02f_gene_PCA_kmeans_k%d_%s%s.png", k, suf, method_tag)),
             p_pca_genes, width = 8, height = 6, dpi = 150, bg = "white")
      message("  Saved: 02f_gene_PCA_kmeans_k", k, "_", suf, method_tag, ".png")

      # ── Gene PCA on the alternate display matrix (02h) ─────────────────────
      # Same clusters colored on the alt-transform PCA. Structure that changes
      # between 02f and 02h is transform-driven, not biological.
      pca_genes_log2 <- prcomp(cm$pca_alt$mat, center = TRUE, scale. = FALSE)
      pca_gl_pct     <- round(100 * pca_genes_log2$sdev^2 / sum(pca_genes_log2$sdev^2), 1)
      pca_gl_df      <- data.frame(
        PC1     = pca_genes_log2$x[, 1],
        PC2     = pca_genes_log2$x[, 2],
        Cluster = factor(paste0("C", km$cluster[rownames(pca_genes_log2$x)]),
                         levels = paste0("C", 1:k)),
        stringsAsFactors = FALSE
      )
      p_pca_genes_log2 <- ggplot(pca_gl_df, aes(PC1, PC2, color = Cluster)) +
        geom_point(size = 0.8, alpha = 0.6) +
        scale_color_manual(values = cluster_pal, name = "k-means\ncluster") +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        labs(
          title    = sprintf("Gene PCA (%s) colored by k-means cluster — k=%d | %s%s | %s",
                             cm$pca_alt$label, k, suf, method_tag, cm$method_label),
          subtitle = paste0(nrow(pca_gl_df), " genes | PCA of ", cm$pca_alt$label,
                            " | clusters from ", cm$method_label, "\n",
                            "Compare to 02f: structure that changes here is transform-driven"),
          x = paste0("PC1 (", pca_gl_pct[1], "%)"),
          y = paste0("PC2 (", pca_gl_pct[2], "%)")
        ) +
        theme_cowplot(12)
      ggsave(file.path(DIR_CLUST,
                       sprintf("02h_gene_PCA_log2_kmeans_k%d_%s%s.png", k, suf, method_tag)),
             p_pca_genes_log2, width = 8, height = 6, dpi = 150, bg = "white")
      message("  Saved: 02h_gene_PCA_log2_kmeans_k", k, "_", suf, method_tag, ".png")

      # ── Silhouette plot (requires k >= 2; undefined for k=1) ──────────────
      if (k >= 2) {
        dist_genes_km <- dist(cm$sil_dist_mat)
        sil_obj <- silhouette(km$cluster, dist_genes_km)
        sil_df  <- as.data.frame(sil_obj[, 1:3])
        colnames(sil_df) <- c("Cluster", "Neighbor", "SilWidth")
        sil_df$Gene    <- rownames(cm$sil_dist_mat)
        sil_df$Cluster <- factor(paste0("C", sil_df$Cluster), levels = paste0("C", 1:k))
        sil_df <- sil_df %>%
          group_by(Cluster) %>%
          arrange(desc(SilWidth), .by_group = TRUE) %>%
          mutate(gene_idx = row_number()) %>%
          ungroup()

        avg_sil <- sil_df %>%
          group_by(Cluster) %>%
          summarise(avg_sil = round(mean(SilWidth), 3), .groups = "drop")
        avg_sil_label <- paste(paste0(avg_sil$Cluster, "=", avg_sil$avg_sil), collapse = " | ")
        sil_summary_list[[paste0(suf, method_tag, "_k", k)]] <- data.frame(
          gs = paste0(suf, method_tag), k = k,
          Cluster  = avg_sil$Cluster,
          mean_sil = avg_sil$avg_sil,
          stringsAsFactors = FALSE
        )

        p_sil <- ggplot(sil_df, aes(x = gene_idx, y = SilWidth, fill = Cluster)) +
          geom_col(width = 1) +
          facet_wrap(~Cluster, scales = "free_x", nrow = 1) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
          scale_fill_manual(values = cluster_pal, guide = "none") +
          labs(
            title    = sprintf("Silhouette plot — k=%d | %s%s | %s",
                               k, suf, method_tag, cm$method_label),
            subtitle = paste0("Mean silhouette by cluster: ", avg_sil_label,
                              " (PROMO02/04 excluded at load from cohort)"),
            x = "Gene (sorted by silhouette width within cluster)",
            y = "Silhouette width"
          ) +
          theme_cowplot(11) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                strip.text = element_text(face = "bold"))
        ggsave(file.path(DIR_CLUST,
                         sprintf("02b_silhouette_k%d_%s%s.png", k, suf, method_tag)),
               p_sil,
               width  = max(8, k * 2.5), height = 4,
               dpi = 150, bg = "white")
        message("  Saved: 02b_silhouette_k", k, "_", suf, method_tag, ".png")
      } else {
        message("  k=1: silhouette undefined — skipping")
      }

      # ── Per-cluster absolute expression violin (diagnostic) ────────────────
      # Shows log2(CPM+1) per cluster to reveal whether a cluster is genuinely
      # lowly expressed or only appears equal after scaling. Uses
      # names(km$cluster) directly so the z-score pass (which may have dropped
      # zero-sd genes) shows exactly the genes that were clustered.
      clust_gene_ids   <- intersect(names(km$cluster), rownames(cpm_filt))
      clust_log2cpm_v1 <- log2(as.matrix(cpm_filt[clust_gene_ids, v1_ids_km]) + 1)
      clust_expr_long  <- do.call(rbind, lapply(seq_len(k), function(ci) {
        gene_ids <- names(km$cluster)[km$cluster == ci]
        valid    <- gene_ids[gene_ids %in% rownames(clust_log2cpm_v1)]
        if (length(valid) == 0) return(NULL)
        data.frame(
          Cluster = paste0("C", ci),
          log2cpm = as.vector(clust_log2cpm_v1[valid, ]),
          stringsAsFactors = FALSE
        )
      }))
      clust_expr_long$Cluster <- factor(clust_expr_long$Cluster,
                                        levels = paste0("C", seq_len(k)))
      p_clust_expr <- ggplot(clust_expr_long, aes(x = Cluster, y = log2cpm, fill = Cluster)) +
        geom_violin(alpha = 0.70, trim = TRUE) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.4, outlier.alpha = 0.4) +
        scale_fill_manual(values = cluster_pal, guide = "none") +
        labs(
          title    = sprintf("Per-cluster log2(CPM+1) — k=%d | %s%s | %s",
                             k, suf, method_tag, cm$method_label),
          subtitle = paste0(
            "Absolute expression before scaling\n",
            "Clusters with low log2(CPM+1) appear identical on the heatmap despite lower expression"
          ),
          x = "Cluster", y = "log2(CPM+1) across V1 PROMO samples"
        ) +
        theme_cowplot(12)
      ggsave(file.path(DIR_CLUST,
                       sprintf("02e_cluster_expression_violin_k%d_%s%s.png", k, suf, method_tag)),
             p_clust_expr,
             width = max(6, k * 1.4), height = 5,
             dpi = 150, bg = "white")
      message("  Saved: 02e_cluster_expression_violin_k", k, "_", suf, method_tag, ".png")

      # ── Per-cluster ORA — clusterProfiler only; facet-grid output ──────────
      # Foreground: gene symbols in each cluster (for the current method).
      # Background/universe: all expressed gene symbols in the cleaned V1 cohort.
      # Output: individual CSVs per cluster + one combined facet-grid PNG per k.
      # k=1 is skipped (all HVGs in one cluster — no contrast possible).
      if (k >= 2) {
        message("  Running per-cluster ORA — clusterProfiler (k=", k, ", ",
                suf, method_tag, ")")

        ora_bg_syms <- unique(na.omit(
          gene_sym_map$Gene_Symbol[gene_sym_map$ENSG %in% rownames(all_genes_norm_km)]
        ))

        run_ora_cp <- function(fg_syms, bg_syms, pathway_list, min_gs = 5, max_gs = 500) {
          if (!HAS_CP) return(NULL)
          t2g <- do.call(rbind, lapply(names(pathway_list), function(pw) {
            data.frame(term = pw, gene = pathway_list[[pw]], stringsAsFactors = FALSE)
          }))
          res <- tryCatch(
            clusterProfiler::enricher(
              gene          = fg_syms,
              universe      = bg_syms,
              TERM2GENE     = t2g,
              minGSSize     = min_gs,
              maxGSSize     = max_gs,
              pAdjustMethod = "BH",
              pvalueCutoff  = 1,
              qvalueCutoff  = 1
            ),
            error = function(e) NULL
          )
          if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
          df <- as.data.frame(res)
          data.frame(
            Description = df$ID,
            Count       = df$Count,
            GeneRatio   = sapply(df$GeneRatio, function(x) eval(parse(text = x))),
            BgRatio     = sapply(df$BgRatio,   function(x) eval(parse(text = x))),
            pvalue      = df$pvalue,
            p.adjust    = df$p.adjust,
            geneID      = df$geneID,
            stringsAsFactors = FALSE
          )
        }

        # Single-cluster panel builder (for facet grid assembly)
        make_ora_panel <- function(df, cluster_label, top_n = 15, strip_prefix = NULL) {
          if (is.null(df) || nrow(df) == 0) return(NULL)
          plot_df <- df %>%
            arrange(p.adjust) %>%
            slice_head(n = top_n) %>%
            mutate(
              neglog10padj = pmin(-log10(p.adjust), 10),
              pw_label     = if (!is.null(strip_prefix))
                               gsub(strip_prefix, "", Description, perl = TRUE)
                             else Description,
              pw_label     = gsub("_", " ", pw_label),
              pw_label     = factor(pw_label, levels = rev(pw_label))
            )
          ggplot(plot_df, aes(x = GeneRatio, y = pw_label,
                               size = Count, color = neglog10padj)) +
            geom_point() +
            scale_color_gradient(low = "#AAAAAA", high = "#B2182B",
                                 limits = c(0, 10), oob = scales::squish,
                                 name = "-log10(p.adj)") +
            scale_size_continuous(range = c(2, 7), name = "Genes") +
            labs(title = cluster_label, x = "Gene Ratio", y = NULL) +
            theme_cowplot(9) +
            theme(
              plot.title  = element_text(face = "bold", size = 10),
              axis.text.y = element_text(size = 7),
              panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
              legend.position = "bottom",
              legend.text  = element_text(size = 7),
              legend.title = element_text(size = 8)
            )
        }

        cp_h_list  <- list()
        cp_bp_list <- list()

        for (clust_idx in seq_len(k)) {
          fg_ensg <- names(km$cluster)[km$cluster == clust_idx]
          fg_syms <- unique(na.omit(
            gene_sym_map$Gene_Symbol[match(fg_ensg, gene_sym_map$ENSG)]
          ))
          clust_label <- sprintf("k%d_C%d_%s%s", k, clust_idx, suf, method_tag)

          cp_h  <- run_ora_cp(fg_syms, ora_bg_syms, pathways_h,    min_gs = 5, max_gs = 500)
          cp_bp <- run_ora_cp(fg_syms, ora_bg_syms, pathways_gobp, min_gs = 5, max_gs = 500)

          cp_h_list[[paste0("C", clust_idx)]]  <- cp_h
          cp_bp_list[[paste0("C", clust_idx)]] <- cp_bp

          if (!is.null(cp_h)  && nrow(cp_h)  > 0)
            write.csv(cp_h,  file.path(DIR_ORA,
              sprintf("ORA_cluster_%s_Hallmark_CP.csv", clust_label)), row.names = FALSE)
          if (!is.null(cp_bp) && nrow(cp_bp) > 0)
            write.csv(cp_bp, file.path(DIR_ORA,
              sprintf("ORA_cluster_%s_GOBP_CP.csv", clust_label)), row.names = FALSE)
        }

        # Assemble all cluster panels into one facet-grid figure per pathway set
        for (ora_spec in list(
          list(cp_list = cp_h_list,  pw_type = "Hallmark", strip = "^HALLMARK_", tag = "Hallmark"),
          list(cp_list = cp_bp_list, pw_type = "GO:BP",    strip = "^GOBP_",     tag = "GOBP")
        )) {
          panels <- lapply(names(ora_spec$cp_list), function(cn) {
            make_ora_panel(ora_spec$cp_list[[cn]], cluster_label = cn,
                            top_n = 15, strip_prefix = ora_spec$strip)
          })
          panels <- panels[!sapply(panels, is.null)]
          if (length(panels) > 0) {
            n_cols_fg <- min(3, length(panels))
            n_rows_fg <- ceiling(length(panels) / n_cols_fg)
            fn_fg     <- sprintf("02c_ORA_CP_facet_k%d_%s%s_%s.png",
                                 k, suf, method_tag, ora_spec$tag)
            ggsave(file.path(DIR_ORA, fn_fg),
                   plot_grid(plotlist = panels, ncol = n_cols_fg),
                   width  = 10 * n_cols_fg,
                   height =  7 * n_rows_fg,
                   dpi = 150, bg = "white")
            message("  Saved: ", fn_fg)
          }
        }
        message("  Per-cluster ORA complete for k=", k, " (", suf, method_tag, ")")

        # ── Leading-edge boxplots for selected Hallmark pathways ──────────────
        # For each target pathway, emit one boxplot per cluster where the pathway
        # is significant (p.adjust < ALPHA). enricher above uses pvalueCutoff = 1,
        # so cp_h_list contains every tested pathway — the ALPHA filter here is
        # what actually defines "significant cluster" for this pathway.
        hallmark_box_targets <- c(
          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        )
        for (target_pw in hallmark_box_targets) {
          hit_clusters <- character(0)
          for (cn in names(cp_h_list)) {
            df_c <- cp_h_list[[cn]]
            if (is.null(df_c) || nrow(df_c) == 0) next
            row_match <- df_c[df_c$Description == target_pw, , drop = FALSE]
            if (nrow(row_match) == 0) next
            if (!is.na(row_match$p.adjust[1]) && row_match$p.adjust[1] < ALPHA)
              hit_clusters <- c(hit_clusters, cn)
          }
          if (length(hit_clusters) == 0) {
            message("  ", target_pw, " not significant in any cluster for ",
                    suf, method_tag, " — skipping boxplot")
            next
          }
          for (cn in hit_clusters) {
            df_c       <- cp_h_list[[cn]]
            row_match  <- df_c[df_c$Description == target_pw, , drop = FALSE]
            padj_cn    <- row_match$p.adjust[1]
            genes_cn   <- strsplit(row_match$geneID[1], "/", fixed = TRUE)[[1]]
            # Filter leading-edge symbols to genes that are HVGs for the CURRENT
            # method (top_var_genes is the in-scope HVG set for this suf — see
            # the gene_sets loop at line ~1930). The cluster foreground was built
            # from top_var_genes upstream, so this should normally drop zero
            # genes; kept as a defensive check so a boxplot can never display a
            # gene that wasn't in the same HVG set as the clustered heatmap.
            ensg_all <- gene_sym_map$ENSG[match(genes_cn, gene_sym_map$Gene_Symbol)]
            keep_hvg <- !is.na(ensg_all) &
                        ensg_all %in% top_var_genes &
                        ensg_all %in% rownames(cpm_filt)
            n_drop   <- sum(!keep_hvg)
            if (n_drop > 0)
              message("  ", target_pw, " in ", cn, " (", suf, method_tag,
                      "): dropped ", n_drop,
                      " leading-edge gene(s) not in HVG set for this method")
            ensg_all <- ensg_all[keep_hvg]
            sym_all  <- genes_cn[keep_hvg]
            # Take first 6 from the HVG-filtered list (order preserved from geneID)
            take_ix  <- seq_len(min(6, length(ensg_all)))
            ensg_6   <- ensg_all[take_ix]
            sym_6    <- sym_all [take_ix]
            if (length(ensg_6) == 0) {
              message("  No HVG leading-edge genes for ", target_pw, " in ", cn,
                      " (", suf, method_tag, ") — skipping boxplot")
              next
            }
            cpm_box <- as.matrix(cpm_filt[ensg_6, v1_ids_km, drop = FALSE])
            rownames(cpm_box) <- sym_6
            box_df <- data.frame(
              Gene     = rep(rownames(cpm_box), times = ncol(cpm_box)),
              Sample   = rep(short_label(colnames(cpm_box)), each = nrow(cpm_box)),
              log2cpm  = as.vector(log2(cpm_box + 1)),
              stringsAsFactors = FALSE
            )
            box_df$Gene <- factor(box_df$Gene, levels = rownames(cpm_box))

            # Per-gene 1.5*IQR outlier flag on log2(CPM+1)
            box_df <- box_df %>%
              group_by(Gene) %>%
              mutate(
                .q1  = quantile(log2cpm, 0.25, na.rm = TRUE),
                .q3  = quantile(log2cpm, 0.75, na.rm = TRUE),
                .iqr = .q3 - .q1,
                is_outlier = log2cpm < (.q1 - 1.5 * .iqr) |
                             log2cpm > (.q3 + 1.5 * .iqr)
              ) %>%
              ungroup() %>%
              select(-.q1, -.q3, -.iqr)

            # Per-gene 5% / 95% percentile markers on log2(CPM+1)
            pctl_df <- box_df %>%
              group_by(Gene) %>%
              summarise(
                p05 = quantile(log2cpm, 0.05, na.rm = TRUE),
                p95 = quantile(log2cpm, 0.95, na.rm = TRUE),
                .groups = "drop"
              ) %>%
              mutate(x_num = as.numeric(Gene))

            pw_label <- gsub("_", " ", gsub("^HALLMARK_", "", target_pw))
            p_box_le <- ggplot(box_df, aes(x = Gene, y = log2cpm)) +
              geom_boxplot(fill = "#88A8C9", alpha = 0.7, outlier.shape = NA) +
              geom_segment(
                data = pctl_df,
                aes(x = x_num - 0.4, xend = x_num + 0.4, y = p05, yend = p05),
                color = "#7B2D8B", linetype = "dashed", linewidth = 0.5,
                inherit.aes = FALSE
              ) +
              geom_segment(
                data = pctl_df,
                aes(x = x_num - 0.4, xend = x_num + 0.4, y = p95, yend = p95),
                color = "#7B2D8B", linetype = "dashed", linewidth = 0.5,
                inherit.aes = FALSE
              ) +
              geom_jitter(data = box_df %>% filter(!is_outlier),
                          width = 0.18, size = 0.9, alpha = 0.55, color = "grey25") +
              geom_point(data = box_df %>% filter(is_outlier),
                         size = 1.6, color = "#B2182B") +
              ggrepel::geom_text_repel(
                data = box_df %>% filter(is_outlier),
                aes(label = Sample),
                size = 2.6, color = "grey15",
                box.padding = 0.35, point.padding = 0.2,
                max.overlaps = 50,
                segment.size = 0.25, segment.color = "grey50",
                min.segment.length = 0
              ) +
              labs(
                title    = paste0(pw_label, " — leading-edge genes (first 6)"),
                subtitle = paste0(
                  "Cluster ", cn, " (k=", k, ", ", suf, method_tag,
                  " — ", cm$method_label, ")",
                  " | p.adj = ", signif(padj_cn, 3),
                  " | log2(CPM+1) across ", length(v1_ids_km), " V1 PROMO samples",
                  " | dashed = per-gene 5%/95% pctile | red = 1.5×IQR outliers (heatmap ID)"
                ),
                x = NULL, y = "log2(CPM + 1)"
              ) +
              theme_cowplot(11) +
              theme(axis.text.x = element_text(angle = 30, hjust = 1))

            fn_box <- sprintf("02c_leadingedge_boxplot_%s_%s_k%d_%s%s.png",
                              target_pw, cn, k, suf, method_tag)
            ggsave(file.path(DIR_ORA, fn_box), p_box_le,
                   width = 8, height = 5, dpi = 150, bg = "white")
            message("  Saved: ", fn_box)
          }
        }
      }
    }  # end cluster_methods loop
  }    # end k loop

  # Save cluster membership for largest k — one CSV per clustering method
  for (cm in cluster_methods) {
    km_final <- kmeans_results[[paste0("k", max(K_VALUES), cm$file_tag)]]
    if (is.null(km_final)) next
    cluster_df <- data.frame(
      ENSG        = names(km_final$cluster),
      Gene_Symbol = gene_sym_map$Gene_Symbol[match(names(km_final$cluster), gene_sym_map$ENSG)],
      Cluster     = km_final$cluster,
      variance    = gene_vars_gs[names(km_final$cluster)],
      stringsAsFactors = FALSE
    ) %>% arrange(Cluster, desc(variance))
    write.csv(cluster_df,
              file.path(DIR_CLUST, sprintf("kmeans_k%d_gene_clusters_%s%s.csv",
                                          max(K_VALUES), suf, cm$file_tag)),
              row.names = FALSE)
  }

  # ── Z-score-per-patient CSVs + cluster-mean line plots + correlation matrix
  #
  # For each HVG set we have two parallel k-means passes (raw CPM min-max and
  # z-score). The shared gene universe is names(km_z$cluster) — z_v1 drops
  # zero-sd genes that the raw pass keeps. We export the z-score matrix once
  # per cluster mode (rows sorted by that mode's cluster), then derive:
  #   - per-cluster mean z-score per sample → faceted line plot in heatmap col order
  #   - per-cluster mean z-score vector per patient → correlated against pain
  #     metadata (qst_ppt_v1, promis_v1, global07_ge4_v1, age, sex_male)
  km_raw_g <- kmeans_results[[paste0("k", max(K_VALUES))]]
  km_z_g   <- kmeans_results[[paste0("k", max(K_VALUES), "_ZSCORE")]]
  if (!is.null(km_raw_g) && !is.null(km_z_g)) {
    shared_genes_kg <- intersect(names(km_raw_g$cluster), names(km_z_g$cluster))
    if (length(shared_genes_kg) > 0) {
      samp_order_kg     <- v1_ids[hc_v1$order]
      samp_short_ord_kg <- short_label(samp_order_kg)
      z_export_mat_kg   <- z_v1[shared_genes_kg, samp_order_kg, drop = FALSE]
      colnames(z_export_mat_kg) <- samp_short_ord_kg

      raw_clusters_kg <- km_raw_g$cluster[shared_genes_kg]
      z_clusters_kg   <- km_z_g$cluster[shared_genes_kg]
      k_used_kg <- max(K_VALUES)

      # Build per-patient pain metadata table once (V1 numeric cols + sex_male).
      # Same SubjectID across PROMO24 tech reps gets identical values — that's
      # expected and only mildly weights PROMO24 in the correlation.
      pain_lookup_kg <- pain_raw
      rownames(pain_lookup_kg) <- pain_lookup_kg$SubjectID
      subj_per_samp_kg <- meta[samp_order_kg, "SubjectID"]
      sex_bin_kg <- ifelse(
        is.na(pain_lookup_kg[subj_per_samp_kg, "sex_at_birth.factor"]), NA_integer_,
        as.integer(pain_lookup_kg[subj_per_samp_kg, "sex_at_birth.factor"] == "Male")
      )
      # Column name → display label. Numeric V1 metrics from pain metadata plus
      # sex_male (binary from sex_at_birth.factor) and sort_batch (integer code).
      # read.csv mangles "Sort batch" → "Sort.batch"; reference the mangled name.
      meta_cols_kg <- c(
        qst_ppt_v1         = "qst_ppt_tr_avg_v1",
        promis_v1          = "promis_global07_v1",
        global07_ge4       = "global07_ge4_v1",
        age                = "age",
        fsq_thinking_v1    = "fsq_thinking_v1",
        fsq_v1             = "fsq_v1",
        wpi_v1             = "wpi_v1",
        sss_v1             = "sss_v1",
        fatigue_tscore_v1  = "promis_cat_v10_fatigue_4_tscore_v1",
        paininterf_tscore_v1 = "promis_cat_v11_pain_interference_4_tscore_v1",
        ptsjc28_v1         = "ptsjc28_sum.correct_v1",
        pttjc28_v1         = "pttjc28_sum.correct_v1",
        CDAI_v1            = "CDAI_v1",
        sort_batch         = "Sort.batch"
      )
      meta_cols_present_kg <- meta_cols_kg[meta_cols_kg %in% colnames(pain_lookup_kg)]
      missing_meta_kg      <- setdiff(meta_cols_kg, colnames(pain_lookup_kg))
      if (length(missing_meta_kg) > 0)
        message("  Correlation matrix: metadata columns not found in PAIN_META, skipped: ",
                paste(missing_meta_kg, collapse = ", "))

      meta_mat_kg <- as.data.frame(
        lapply(meta_cols_present_kg,
               function(cn) as.numeric(pain_lookup_kg[subj_per_samp_kg, cn])),
        stringsAsFactors = FALSE
      )
      names(meta_mat_kg) <- names(meta_cols_present_kg)
      meta_mat_kg$sex_male <- sex_bin_kg
      rownames(meta_mat_kg) <- samp_order_kg

      cluster_pal_kg <- setNames(
        colorRampPalette(brewer.pal(max(k_used_kg, 8), "Set2"))(k_used_kg),
        paste0("C", seq_len(k_used_kg))
      )

      for (mode_tag in c("", "_ZSCORE")) {
        clusters_use_kg <- if (mode_tag == "") raw_clusters_kg else z_clusters_kg
        mode_label_kg   <- if (mode_tag == "") "raw CPM k-means" else "z-score k-means"

        # ── CSV: gene rows × patient columns, both cluster labels included ───
        df_zcsv <- data.frame(
          ENSG           = shared_genes_kg,
          Gene_Symbol    = gene_sym_map$Gene_Symbol[match(shared_genes_kg, gene_sym_map$ENSG)],
          Cluster_raw    = raw_clusters_kg,
          Cluster_zscore = z_clusters_kg,
          stringsAsFactors = FALSE
        )
        df_zcsv <- cbind(df_zcsv, as.data.frame(z_export_mat_kg))
        df_zcsv <- df_zcsv[order(clusters_use_kg), , drop = FALSE]
        fn_zcsv <- sprintf("kmeans_k%d_zscore_per_patient_%s%s.csv",
                           k_used_kg, suf, mode_tag)
        write.csv(df_zcsv, file.path(DIR_CLUST, fn_zcsv), row.names = FALSE)
        message("  Saved: ", fn_zcsv)

        # ── Per-cluster mean z-score per sample (matrix used by line + corr)
        cluster_means_mat_kg <- do.call(rbind, lapply(seq_len(k_used_kg), function(ci) {
          g_ids <- shared_genes_kg[clusters_use_kg == ci]
          if (length(g_ids) == 0) return(rep(NA_real_, length(samp_order_kg)))
          colMeans(z_v1[g_ids, samp_order_kg, drop = FALSE])
        }))
        rownames(cluster_means_mat_kg) <- paste0("C", seq_len(k_used_kg))
        colnames(cluster_means_mat_kg) <- samp_order_kg

        # ── Faceted line plot: one panel per cluster, sample order = heatmap col order
        line_df_kg <- do.call(rbind, lapply(seq_len(k_used_kg), function(ci) {
          vals <- cluster_means_mat_kg[ci, ]
          if (all(is.na(vals))) return(NULL)
          data.frame(
            Cluster    = paste0("C", ci),
            SampleID   = samp_order_kg,
            SampleLab  = samp_short_ord_kg,
            MeanZScore = as.numeric(vals),
            stringsAsFactors = FALSE
          )
        }))
        if (!is.null(line_df_kg) && nrow(line_df_kg) > 0) {
          line_df_kg$Cluster   <- factor(line_df_kg$Cluster,
                                          levels = paste0("C", seq_len(k_used_kg)))
          line_df_kg$SampleLab <- factor(line_df_kg$SampleLab, levels = samp_short_ord_kg)

          p_line_kg <- ggplot(line_df_kg, aes(x = SampleLab, y = MeanZScore,
                                              group = Cluster, color = Cluster)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
                       linewidth = 0.4) +
            geom_line(linewidth = 0.8) +
            geom_point(size = 1.8) +
            facet_wrap(~ Cluster, ncol = 2, scales = "free_y") +
            scale_color_manual(values = cluster_pal_kg, guide = "none") +
            labs(
              title    = sprintf("Per-cluster mean z-score per sample — %s | k=%d | %s%s",
                                  mode_label_kg, k_used_kg, suf, mode_tag),
              subtitle = "Sample order = k-means heatmap column order | y = mean z-score across cluster genes",
              x = NULL, y = "Mean z-score (across cluster genes)"
            ) +
            theme_cowplot(11) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                              size = 7),
                  strip.text  = element_text(face = "bold"))
          fn_line_kg <- sprintf("02i_cluster_mean_zscore_per_sample_k%d_%s%s.png",
                                k_used_kg, suf, mode_tag)
          ggsave(file.path(DIR_CLUST, fn_line_kg), p_line_kg,
                 width  = max(14, length(samp_order_kg) * 0.30),
                 height = 2 + 1.8 * ceiling(k_used_kg / 2),
                 dpi = 150, bg = "white")
          message("  Saved: ", fn_line_kg)

          # Save the underlying matrix
          line_csv_df <- data.frame(
            Cluster   = line_df_kg$Cluster,
            SampleID  = line_df_kg$SampleID,
            SampleLab = as.character(line_df_kg$SampleLab),
            MeanZScore = line_df_kg$MeanZScore,
            stringsAsFactors = FALSE
          )
          write.csv(line_csv_df,
                    file.path(DIR_CLUST,
                              sprintf("kmeans_k%d_cluster_mean_zscore_per_sample_%s%s.csv",
                                      k_used_kg, suf, mode_tag)),
                    row.names = FALSE)
        }

        # ── Correlation matrix: per-cluster mean-z vector vs metadata columns
        cor_mat_kg <- matrix(NA_real_,
                              nrow = k_used_kg, ncol = ncol(meta_mat_kg),
                              dimnames = list(rownames(cluster_means_mat_kg),
                                              colnames(meta_mat_kg)))
        p_mat_kg <- cor_mat_kg
        for (ci in seq_len(k_used_kg)) {
          for (mc in seq_len(ncol(meta_mat_kg))) {
            x_v <- cluster_means_mat_kg[ci, ]
            y_v <- meta_mat_kg[, mc]
            ok  <- !is.na(x_v) & !is.na(y_v)
            if (sum(ok) >= 3 && sd(x_v[ok]) > 0 && sd(y_v[ok]) > 0) {
              ct <- suppressWarnings(cor.test(x_v[ok], y_v[ok], method = "pearson"))
              cor_mat_kg[ci, mc] <- as.numeric(ct$estimate)
              p_mat_kg  [ci, mc] <- ct$p.value
            }
          }
        }

        cor_long_kg <- data.frame(
          Cluster  = rownames(cor_mat_kg)[row(cor_mat_kg)],
          Metadata = colnames(cor_mat_kg)[col(cor_mat_kg)],
          r        = as.vector(cor_mat_kg),
          p        = as.vector(p_mat_kg),
          stringsAsFactors = FALSE
        )
        cor_long_kg$lab <- ifelse(
          is.na(cor_long_kg$r), "",
          sprintf("%.2f%s", cor_long_kg$r,
                  ifelse(!is.na(cor_long_kg$p) & cor_long_kg$p < 0.05, "*", ""))
        )
        cor_long_kg$Cluster  <- factor(cor_long_kg$Cluster,
                                        levels = rownames(cor_mat_kg))
        cor_long_kg$Metadata <- factor(cor_long_kg$Metadata,
                                        levels = colnames(cor_mat_kg))

        p_corr_kg <- ggplot(cor_long_kg, aes(x = Metadata, y = Cluster, fill = r)) +
          geom_tile(color = "white") +
          geom_text(aes(label = lab), size = 3) +
          scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, limits = c(-1, 1),
                                na.value = "grey85") +
          labs(
            title    = sprintf("Per-cluster mean z-score vs pain metadata — %s | k=%d | %s%s",
                                mode_label_kg, k_used_kg, suf, mode_tag),
            subtitle = paste0("Pearson r per (cluster, metadata) | * = p<0.05 | ",
                              "n = ", length(samp_order_kg), " V1 PROMO samples"),
            x = NULL, y = NULL, fill = "r"
          ) +
          theme_cowplot(11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        fn_corr_kg <- sprintf("02j_cluster_vs_metadata_corr_k%d_%s%s.png",
                              k_used_kg, suf, mode_tag)
        ggsave(file.path(DIR_CLUST, fn_corr_kg), p_corr_kg,
               width  = max(8, 1.5 + 0.55 * ncol(meta_mat_kg)),
               height = 1.8 + 0.6 * k_used_kg,
               dpi = 150, bg = "white")
        message("  Saved: ", fn_corr_kg)

        write.csv(cor_long_kg[, c("Cluster", "Metadata", "r", "p")],
                  file.path(DIR_CLUST,
                            sprintf("kmeans_k%d_cluster_vs_metadata_corr_%s%s.csv",
                                    k_used_kg, suf, mode_tag)),
                  row.names = FALSE)
      }
    }
  }

}   # end gene_sets loop

# ─────────────────────────────────────────────────────────────────────────────
# 7b. GLOBAL SILHOUETTE SUMMARY — mean silhouette ± SD across clusters per k
#
# One figure per gene-set approach: x = k value, y = mean silhouette averaged
# across all clusters (error bars = SD across clusters).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [7b] Global silhouette summary ===")

if (length(sil_summary_list) > 0) {
  sil_summary_df <- do.call(rbind, sil_summary_list)
  rownames(sil_summary_df) <- NULL

  sil_global <- sil_summary_df %>%
    group_by(gs, k) %>%
    summarise(
      mean_sil_global = mean(mean_sil),
      sd_sil_global   = sd(mean_sil),
      se_sil_global   = sd(mean_sil) / sqrt(n()),
      n_clusters      = n(),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean_sil_global - sd_sil_global,
      ymax = mean_sil_global + sd_sil_global
    )

  for (gs_name_sil in unique(sil_global$gs)) {
    df_sil <- sil_global %>% filter(gs == gs_name_sil)
    p_sil_global <- ggplot(df_sil, aes(x = k, y = mean_sil_global)) +
      geom_line(color = "#2166AC", linewidth = 1) +
      geom_point(size = 3, color = "#2166AC") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax),
                    width = 0.25, color = "#2166AC", linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
      scale_x_continuous(breaks = df_sil$k) +
      labs(
        title    = paste0("Global silhouette — ", gs_name_sil, " gene selection"),
        subtitle = "Mean silhouette per k \u00b1 SD across clusters | higher = better separation",
        x = "Number of clusters (k)",
        y = "Mean silhouette width"
      ) +
      theme_cowplot(12) +
      theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
    ggsave(file.path(DIR_CLUST, paste0("02d_silhouette_summary_", gs_name_sil, ".png")),
           p_sil_global, width = 7, height = 5, dpi = 150, bg = "white")
    message("  Saved: 02d_silhouette_summary_", gs_name_sil, ".png")
  }

  write.csv(sil_summary_df,
            file.path(DIR_CLUST, "silhouette_summary_all_k_all_gs.csv"),
            row.names = FALSE)
  message("  Saved: silhouette_summary_all_k_all_gs.csv")
} else {
  message("  sil_summary_list is empty — skipping global silhouette summary")
}

# ─────────────────────────────────────────────────────────────────────────────
# 7e. JACCARD COMPARISON — z-score vs raw-CPM k-means cluster overlap
#
# For each HVG method (VST / LFC / LFC_P5P95 / LFC_PCT), the gene_sets loop
# ran k-means twice on the same HVG set: once on raw-CPM min-max values
# (file_tag "") and once on z-scored raw CPM ("_ZSCORE"). This block quantifies
# how stable cluster membership is between the two transforms.
#
# Pairing: greedy maximum-Jaccard — for each raw cluster, take its best
# (highest Jaccard) z-score cluster, drop that z-cluster from the pool, and
# move to the next raw cluster. Ties broken by lower z-cluster ID.
#
# Outputs:
#   kmeans_jaccard_z_vs_raw.csv     — per-(HVG, paired-cluster) overlap stats
#   02g_kmeans_jaccard_boxplot.png  — overlap vs different counts (x = HVG)
# ─────────────────────────────────────────────────────────────────────────────
message("=== [7e] Jaccard cluster overlap — z-score vs raw CPM ===")

K_JACCARD <- max(K_VALUES)

# Greedy max-Jaccard pairing. raw_vec / z_vec: named int vectors (ENSG -> cluster).
greedy_match_clusters <- function(raw_vec, z_vec) {
  raw_ids <- sort(unique(raw_vec))
  z_ids   <- sort(unique(z_vec))
  raw_sets <- lapply(setNames(raw_ids, raw_ids),
                     function(i) names(raw_vec)[raw_vec == i])
  z_sets   <- lapply(setNames(z_ids,   z_ids),
                     function(j) names(z_vec)[z_vec   == j])

  available_z <- as.character(z_ids)
  rows <- list()
  for (i in as.character(raw_ids)) {
    if (length(available_z) == 0) {
      rows[[i]] <- data.frame(
        raw_cluster = as.integer(i), matched_z_cluster = NA_integer_,
        raw_n = length(raw_sets[[i]]), z_n = NA_integer_,
        overlap = 0L, union_n = NA_integer_,
        jaccard = NA_real_, diff_n = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }
    jaccs <- vapply(available_z, function(j) {
      o <- length(intersect(raw_sets[[i]], z_sets[[j]]))
      u <- length(union    (raw_sets[[i]], z_sets[[j]]))
      if (u == 0) 0 else o / u
    }, numeric(1))
    best_j <- available_z[which.max(jaccs)]
    o      <- length(intersect(raw_sets[[i]], z_sets[[best_j]]))
    u      <- length(union    (raw_sets[[i]], z_sets[[best_j]]))
    rows[[i]] <- data.frame(
      raw_cluster       = as.integer(i),
      matched_z_cluster = as.integer(best_j),
      raw_n             = length(raw_sets[[i]]),
      z_n               = length(z_sets[[best_j]]),
      overlap           = o,
      union_n           = u,
      jaccard           = if (u == 0) NA_real_ else o / u,
      diff_n            = u - o,
      stringsAsFactors  = FALSE
    )
    available_z <- setdiff(available_z, best_j)
  }
  do.call(rbind, rows)
}

hvg_methods_7e <- c("VST", "LFC", "LFC_P5P95", "LFC_PCT")
jacc_rows <- list()
for (hvg in hvg_methods_7e) {
  # method_tag "" gives "<hvg>____k<k>"; "_ZSCORE" gives "<hvg>___ZSCORE__k<k>"
  raw_key <- paste0(hvg, "____k", K_JACCARD)
  z_key   <- paste0(hvg, "___ZSCORE__k", K_JACCARD)
  if (!raw_key %in% names(kmeans_clusters_all) ||
      !z_key   %in% names(kmeans_clusters_all)) {
    message("  Skipping ", hvg, " — missing raw or z-score k=", K_JACCARD,
            " cluster assignment")
    next
  }
  pair_df            <- greedy_match_clusters(kmeans_clusters_all[[raw_key]],
                                              kmeans_clusters_all[[z_key]])
  pair_df$hvg_method <- hvg
  jacc_rows[[hvg]]   <- pair_df
}

if (length(jacc_rows) == 0) {
  message("  No HVG method had both raw and z-score cluster results — skipping")
} else {
  jacc_all <- do.call(rbind, jacc_rows)
  jacc_all <- jacc_all[, c("hvg_method", "raw_cluster", "matched_z_cluster",
                            "raw_n", "z_n", "overlap", "union_n",
                            "jaccard", "diff_n")]
  write.csv(jacc_all,
            file.path(DIR_CLUST, "kmeans_jaccard_z_vs_raw.csv"),
            row.names = FALSE)
  message("  Saved: kmeans_jaccard_z_vs_raw.csv")

  box_long <- rbind(
    data.frame(hvg_method = jacc_all$hvg_method, category = "Overlap",
               count = jacc_all$overlap, stringsAsFactors = FALSE),
    data.frame(hvg_method = jacc_all$hvg_method, category = "Different",
               count = jacc_all$diff_n, stringsAsFactors = FALSE)
  )
  present_methods     <- intersect(hvg_methods_7e, unique(box_long$hvg_method))
  box_long$hvg_method <- factor(box_long$hvg_method, levels = present_methods)
  box_long$category   <- factor(box_long$category,   levels = c("Overlap", "Different"))

  p_jacc <- ggplot(box_long, aes(x = hvg_method, y = count, fill = category)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
                 alpha = 0.75, outlier.shape = NA) +
    geom_jitter(aes(color = category),
                position = position_jitterdodge(jitter.width = 0.15,
                                                 dodge.width = 0.8),
                size = 1.6, alpha = 0.9) +
    scale_fill_manual(values  = c("Overlap"   = "#1a9641",
                                  "Different" = "#d7191c")) +
    scale_color_manual(values = c("Overlap"   = "#0f6627",
                                  "Different" = "#8B0000"), guide = "none") +
    labs(
      title    = paste0("k-means cluster stability: z-score vs raw CPM (k=",
                        K_JACCARD, ")"),
      subtitle = paste0("Greedy max-Jaccard pairing | each point = one paired ",
                        "cluster | Overlap = |raw ∩ z|, ",
                        "Different = |raw ∪ z| - |raw ∩ z|"),
      x = "HVG method", y = "Gene count per paired cluster", fill = NULL
    ) +
    theme_cowplot(12) +
    theme(legend.position = "top")
  ggsave(file.path(DIR_CLUST, "02g_kmeans_jaccard_boxplot.png"),
         p_jacc, width = 9, height = 6, dpi = 150, bg = "white")
  message("  Saved: 02g_kmeans_jaccard_boxplot.png")

  jacc_summary <- aggregate(jacc_all$jaccard,
                            by = list(hvg_method = jacc_all$hvg_method),
                            FUN = function(x) round(mean(x, na.rm = TRUE), 3))
  for (rr in seq_len(nrow(jacc_summary))) {
    message("  Mean Jaccard (", jacc_summary$hvg_method[rr], "): ",
            jacc_summary$x[rr])
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 7f. SANKEY — raw CPM k-means clusters → z-score k-means clusters (per HVG)
#
# One alluvial per HVG method (VST, LFC, LFC_P5P95, LFC_PCT). Genes flow from
# their raw-CPM cluster (left axis) to their z-score cluster (right axis), so
# reassignment between the two heatmap modes is visible at a glance.
# Reads from kmeans_clusters_all built during the gene_sets loop.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [7f] Sankey: raw CPM cluster -> z-score cluster (per HVG method) ===")

K_SANKEY <- max(K_VALUES)
for (hvg in hvg_methods_7e) {
  raw_key_sk <- paste0(hvg, "____k", K_SANKEY)
  z_key_sk   <- paste0(hvg, "___ZSCORE__k", K_SANKEY)
  if (!raw_key_sk %in% names(kmeans_clusters_all) ||
      !z_key_sk   %in% names(kmeans_clusters_all)) {
    message("  Skipping Sankey for ", hvg, " - missing raw or z-score assignments")
    next
  }
  raw_vec_sk <- kmeans_clusters_all[[raw_key_sk]]
  z_vec_sk   <- kmeans_clusters_all[[z_key_sk]]
  shared_sk  <- intersect(names(raw_vec_sk), names(z_vec_sk))
  if (length(shared_sk) == 0) {
    message("  Skipping Sankey for ", hvg, " - no shared genes"); next
  }
  cluster_levels_sk <- paste0("C", seq_len(K_SANKEY))
  sank_df <- data.frame(
    ENSG  = shared_sk,
    raw_c = factor(paste0("C", raw_vec_sk[shared_sk]), levels = cluster_levels_sk),
    z_c   = factor(paste0("C", z_vec_sk[shared_sk]),   levels = cluster_levels_sk),
    stringsAsFactors = FALSE
  )

  sank_counts <- sank_df %>%
    dplyr::group_by(raw_c, z_c) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  flow_pal_sk <- setNames(
    colorRampPalette(brewer.pal(max(K_SANKEY, 8), "Set2"))(K_SANKEY),
    cluster_levels_sk
  )

  p_sank <- ggplot(sank_counts,
                   aes(axis1 = raw_c, axis2 = z_c, y = n)) +
    geom_alluvium(aes(fill = raw_c), width = 1/12,
                  colour = "grey30", alpha = 0.85) +
    geom_stratum(width = 1/12, fill = "white", colour = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("raw CPM k-means", "z-score k-means"),
                     expand = c(0.05, 0.05)) +
    scale_fill_manual(values = flow_pal_sk, name = "Raw\ncluster") +
    labs(
      title    = sprintf("Gene cluster reassignment - %s (k=%d, n=%d shared genes)",
                          hvg, K_SANKEY, length(shared_sk)),
      subtitle = "Raw CPM k-means clusters (left) -> z-score min-max k-means clusters (right)",
      x = NULL, y = "Number of genes"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x     = element_text(face = "bold", size = 11),
          legend.position = "right")

  fn_sank <- sprintf("02k_sankey_raw_vs_zscore_k%d_%s.png", K_SANKEY, hvg)
  ggsave(file.path(DIR_CLUST, fn_sank), p_sank,
         width = 10, height = 7, dpi = 150, bg = "white")
  write.csv(sank_counts,
            file.path(DIR_CLUST,
                      sprintf("02k_sankey_counts_k%d_%s.csv", K_SANKEY, hvg)),
            row.names = FALSE)
  message("  Saved: ", fn_sank)

  # Per-gene reassignment table (every shared gene, both cluster labels)
  reassign_df <- data.frame(
    ENSG        = shared_sk,
    Gene_Symbol = gene_sym_map$Gene_Symbol[match(shared_sk, gene_sym_map$ENSG)],
    Cluster_raw    = as.integer(raw_vec_sk[shared_sk]),
    Cluster_zscore = as.integer(z_vec_sk[shared_sk]),
    Reassigned     = as.integer(raw_vec_sk[shared_sk]) != as.integer(z_vec_sk[shared_sk]),
    stringsAsFactors = FALSE
  ) %>% arrange(Cluster_raw, Cluster_zscore)
  write.csv(reassign_df,
            file.path(DIR_CLUST,
                      sprintf("02k_sankey_gene_reassignment_k%d_%s.csv",
                              K_SANKEY, hvg)),
            row.names = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7c. LFC OUTLIER-FRACTION DIAGNOSTIC
#
# Explains why min-max LFC silhouette tends to be higher than P10-P90 despite
# P10-P90 being a subset:
#
# min-max-only genes are outlier-driven — a single extreme sample determines
# max or min. In k-means these genes create artificial cohesion: if all outlier
# samples for a set of genes happen to share the same cluster, those genes will
# be geometrically tight (high silhouette) regardless of biological coherence.
#
# Visualizations:
#   Panel A — scatter of P10-P90 LFC vs min-max LFC, colored by HVG membership.
#             Points above the identity line = min-max inflated by outlier(s).
#   Panel B — histogram of outlier fraction = (min-max − P10-P90) / min-max
#             for min-max-only genes. Values > 0.5 = >50% of the range driven
#             by a single extreme sample.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [7c] LFC outlier-fraction diagnostic ===")

lfc_diag_df <- data.frame(
  ENSG       = rownames(cpm_filt),
  lfc_minmax = gene_vars_lfc[rownames(cpm_filt)],
  lfc_p5p95  = gene_vars_lfc_p5p95[rownames(cpm_filt)],
  lfc_p10p90 = gene_vars_lfc_pct[rownames(cpm_filt)],
  stringsAsFactors = FALSE
) %>%
  filter(lfc_minmax >= 1) %>%
  mutate(
    category = case_when(
      ENSG %in% top_var_lfc & ENSG %in% top_var_lfc_pct  ~ "min-max & P10-P90",
      ENSG %in% top_var_lfc & ENSG %in% top_var_lfc_p5p95 ~ "min-max & P5-P95 (not P10-P90)",
      ENSG %in% top_var_lfc ~ "min-max only",
      TRUE ~ "other"
    ),
    outlier_excess = lfc_minmax - lfc_p10p90,
    outlier_frac   = ifelse(lfc_minmax > 0, outlier_excess / lfc_minmax, 0)
  ) %>%
  filter(category != "other")

cat_pal_diag <- c(
  "min-max & P10-P90"              = "#2166AC",
  "min-max & P5-P95 (not P10-P90)" = "#E08214",
  "min-max only"                   = "#d84b4b"
)

p_lfc_scatter <- ggplot(lfc_diag_df,
                         aes(x = lfc_p10p90, y = lfc_minmax, color = category)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = cat_pal_diag, name = NULL) +
  labs(
    title    = "HVG decomposition: min-max LFC vs P10-P90 LFC",
    subtitle = paste0(
      "Identity line = P10-P90 equals min-max (no outlier inflation)\n",
      "Points above line: min-max range inflated by 1-2 extreme samples\n",
      "Higher min-max silhouette explanation: outlier-driven genes create\n",
      "artificial k-means cohesion when extreme samples co-cluster"
    ),
    x = "P10-P90 LFC (robust range)",
    y = "min-max LFC (full range)"
  ) +
  theme_cowplot(11) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 8, color = "grey30", lineheight = 1.3))

p_outlier_hist <- ggplot(
  lfc_diag_df %>% filter(category == "min-max only"),
  aes(x = outlier_frac)
) +
  geom_histogram(binwidth = 0.05, fill = "#d84b4b", alpha = 0.75, color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.52, y = Inf, vjust = 1.5, hjust = 0,
           label = ">0.5: outlier-driven", color = "grey30", size = 3) +
  labs(
    title    = "Outlier fraction — min-max-only HVGs",
    subtitle = paste0(
      "Fraction of LFC unexplained by P10-P90\n",
      ">0.5: majority of range driven by a single extreme sample"
    ),
    x = "(min-max LFC − P10-P90 LFC) / min-max LFC",
    y = "Gene count"
  ) +
  theme_cowplot(11)

lfc_diag_grid <- plot_grid(p_lfc_scatter, p_outlier_hist, ncol = 2)
ggsave(file.path(DIR_HVG, "00_lfc_outlier_fraction_diagnostic.png"),
       lfc_diag_grid, width = 14, height = 6, dpi = 150, bg = "white")
message("  Saved: 00_lfc_outlier_fraction_diagnostic.png")


# ─────────────────────────────────────────────────────────────────────────────
# 8. DESEQ2 HELPER FUNCTION
#
# IMPORTANT — lfcThreshold is NOT passed to results()/lfcShrink().
# Passing it there triggers a TREAT-style test (H0: |LFC| <= threshold)
# which is highly conservative and gives padj~1 at clinical sample sizes.
# The LFC threshold is applied post-hoc in the Significant column only.
# ─────────────────────────────────────────────────────────────────────────────
run_deseq2 <- function(counts, col_data, design_formula, contrast_or_name,
                        label, alpha = ALPHA,
                        results_type = c("contrast", "name", "apeglm"),
                        out_dir = OUT_DIR) {
  results_type <- match.arg(results_type)

  dds  <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData   = col_data,
    design    = design_formula
  )
  # Per-sample raw-count filter: gene must have count >= DESEQ_MIN_COUNT in
  # at least MIN_PCT_SAMPLES of the samples in this DE comparison.
  # This mirrors the CPM filter (Section 3) but operates on raw counts —
  # which is what DESeq2 actually models — and uses the same percentage
  # threshold so the two filters stay in sync.
  deseq_min_n <- ceiling(MIN_PCT_SAMPLES * ncol(dds))
  keep <- rowSums(counts(dds) >= DESEQ_MIN_COUNT) >= deseq_min_n
  message("  DESeq2 filter: count >= ", DESEQ_MIN_COUNT,
          " in >= ", deseq_min_n, " samples (",
          round(MIN_PCT_SAMPLES*100), "% of ", ncol(dds),
          ") | ", nrow(dds), " \u2192 ", sum(keep), " genes")
  dds  <- dds[keep, ]
  dds  <- DESeq(dds, quiet = TRUE)

  res <- switch(results_type,
    contrast = results(dds, contrast        = contrast_or_name,
                       alpha = alpha,
                       independentFiltering = TRUE),
    name     = results(dds, name            = contrast_or_name,
                       alpha = alpha),
    apeglm   = lfcShrink(dds, coef         = contrast_or_name,
                         type = "apeglm")
  )

  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSG") %>%
    left_join(gene_sym_map, by = "ENSG") %>%
    mutate(
      Gene        = Gene_Symbol,
      # For continuous variables LFC is a gradient (not a group ratio),
      # so significance is defined by pvalue only — no LFC threshold.
      # Lowered nominal p-value cutoff to 0.05.
      Significant = !is.na(pvalue) & pvalue < 0.05
    ) %>%
    arrange(pvalue)

  out_file <- file.path(out_dir, paste0("DESeq2_", label, ".csv"))
  write.csv(res_df, out_file, row.names = FALSE)
  message("  Saved: ", basename(out_file),
          " | Sig genes (pvalue<0.05): ",
          sum(res_df$Significant, na.rm=TRUE))

  # Volcano — x = LFC (gradient direction), y = -log10(pvalue)
  plot_df    <- res_df %>% filter(!is.na(pvalue))
  top_labels <- plot_df %>% filter(Significant) %>% slice_min(pvalue, n = 20)
  if (nrow(top_labels) == 0) top_labels <- plot_df %>% slice_min(pvalue, n = 10)

  vp <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue), color = Significant)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("TRUE" = "#d84b4b", "FALSE" = "#888888")) +
    geom_text_repel(data = top_labels, aes(label = Gene_Symbol),
                    size = 2.5, max.overlaps = 20, box.padding = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "navy") +
    labs(title = label,
         x = "log2 Fold Change per SD of Predictor",
         y = "-log10(pvalue)") +
    theme_cowplot(12) + theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0("Volcano_", label, ".png")),
         vp, width = 7, height = 6, dpi = 150, bg = "white")

  list(dds = dds, res = res, res_df = res_df)
}

# ─────────────────────────────────────────────────────────────────────────────
# 9. DESEQ2 — BASELINE V1 ~ qst_ppt_tr_avg_v1 (CONTINUOUS PPT)
#
# Uses V1 PROMO samples only. PPT is centered+scaled (see section 2).
#
# HOW TO INTERPRET:
#   The log2FoldChange column is the estimated change in log2 expression
#   per ONE STANDARD DEVIATION increase in PPT (because we scaled the variable).
#   Higher PPT = less pain sensitive.
#
#   Positive LFC: gene is higher in subjects with high PPT (less pain sensitive)
#   Negative LFC: gene is higher in subjects with low PPT (more pain sensitive)
#
#   This is fundamentally a regression coefficient, NOT a group comparison.
#   The Wald test p-value tests H0: coefficient = 0 (no linear relationship).
#   padj is BH-corrected across all tested genes.
#
#   Because PPT is continuous there is no single "fold change" between two
#   groups — LFC here means: for every 1 SD (~1.3 units) increase in PPT,
#   expression changes by 2^LFC fold. Use it directionally, not as a ratio.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [9] DESeq2: Baseline V1 ~ qst_ppt_tr_avg_v1 (continuous PPT) ===")

meta_v1 <- meta %>%
  filter(Visit == "Visit1", !IsPerlman, !is.na(qst_ppt_tr_avg_v1)) %>%
  mutate(Visit = droplevels(Visit))

# Pass raw (unfiltered) counts — DESeq2's internal count filter (DESEQ_MIN_COUNT
# in >= MIN_PCT_SAMPLES of samples) defines its own gene universe independently
# of CPM_THRESH. Changing CPM_THRESH should not affect DE results.
counts_v1 <- raw[ , rownames(meta_v1)]
message("  V1 samples with PPT data: ", nrow(meta_v1))
message("  PPT (scaled) mean: ", round(mean(meta_v1$qst_ppt_tr_avg_v1), 3),
        " | SD: ", round(sd(meta_v1$qst_ppt_tr_avg_v1), 3))

de_ppt <- run_deseq2(
  counts           = counts_v1,
  col_data         = meta_v1,
  design_formula   = ~ qst_ppt_tr_avg_v1,
  contrast_or_name = "qst_ppt_tr_avg_v1",
  label            = "V1_continuous_PPT",
  results_type     = "name",
  out_dir          = DIR_DESEQ2
)

# ─────────────────────────────────────────────────────────────────────────────
# 10. EXAMPLE GENE VISUALISATION
#
# For top DE genes: scatter plot of VST expression vs raw PPT value.
# Shows the continuous relationship directly.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [10] Example gene visualisation ===")

top_genes <- de_ppt$res_df %>%
  filter(!is.na(pvalue)) %>%
  filter(Significant) %>%
  filter(ENSG %in% rownames(cpm_filt)) %>%
  slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
  pull(ENSG)

if (length(top_genes) > 0) {
  # CPM expression for V1 samples present in cpm_filt (excludes HVG-excluded subjects)
  cpm_v1 <- cpm_filt[, intersect(rownames(meta_v1), colnames(cpm_filt))]

  plot_list <- lapply(top_genes, function(g) {
    sym <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    lfc <- round(de_ppt$res_df$log2FoldChange[de_ppt$res_df$ENSG == g], 3)
    pnom  <- signif(de_ppt$res_df$pvalue[de_ppt$res_df$ENSG == g], 3)
    df  <- data.frame(
      PPT  = meta_v1[colnames(cpm_v1), "qst_ppt_raw"],
      expr = as.numeric(cpm_v1[g, ]),
      stringsAsFactors = FALSE
    )
    ggplot(df, aes(PPT, expr)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      labs(title = sym,
           subtitle = paste0("LFC/SD=", lfc, "  p=", pnom),
           x = "PPT (raw)", y = "CPM") +
      theme_cowplot(10) +
      theme(plot.title    = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 7))
  })

  n_cols <- min(4, length(plot_list))
  n_rows <- ceiling(length(plot_list) / n_cols)

  png(file.path(DIR_DESEQ2, "DEseq_example_genes_V1.png"),
      width = n_cols * 400, height = n_rows * 380, res = 120, type = "cairo")
  print(plot_grid(plotlist = plot_list, ncol = n_cols))
  dev.off()
  message("  Saved: DEseq_example_genes_V1.png (",
          length(plot_list), " genes)")
} else {
  message("  No significant genes at pvalue<0.05 — showing top genes by pvalue.")
  message("  Plotting top genes by nominal pvalue as fallback")

  # Fallback: top genes by unadjusted p-value, restricted to cpm_filt universe
  top_genes_nominal <- de_ppt$res_df %>%
    filter(!is.na(pvalue)) %>%
    filter(ENSG %in% rownames(cpm_filt)) %>%
    slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
    pull(ENSG)

  cpm_v1 <- cpm_filt[, intersect(rownames(meta_v1), colnames(cpm_filt))]

  plot_list <- lapply(top_genes_nominal, function(g) {
    sym <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    lfc <- round(de_ppt$res_df$log2FoldChange[de_ppt$res_df$ENSG == g], 3)
    pnom <- signif(de_ppt$res_df$pvalue[de_ppt$res_df$ENSG == g], 3)
    df  <- data.frame(
      PPT  = meta_v1[colnames(cpm_v1), "qst_ppt_raw"],
      expr = as.numeric(cpm_v1[g, ]),
      stringsAsFactors = FALSE
    )
    ggplot(df, aes(PPT, expr)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      labs(title = sym,
           subtitle = paste0("LFC/SD=", lfc, "  pnom=", pnom, " (nominal)"),
           x = "PPT (raw)", y = "CPM") +
      theme_cowplot(10) +
      theme(plot.title    = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 7, color = "#888888"))
  })

  n_cols <- min(4, length(plot_list))
  n_rows <- ceiling(length(plot_list) / n_cols)

  png(file.path(DIR_DESEQ2, "DEseq_example_genes_V1_nominal.png"),
      width = n_cols * 400, height = n_rows * 380, res = 120, type = "cairo")
  print(plot_grid(plotlist = plot_list, ncol = n_cols))
  dev.off()
  message("  Saved: DEseq_example_genes_V1_nominal.png (top by pvalue, unadjusted)")
}

# ─────────────────────────────────────────────────────────────────────────────
# 11. GSEA — Hallmark gene sets, ranked by DESeq2 Wald stat
#
# Ranking on the Wald statistic (stat column) is preferred over LFC alone
# because it incorporates both effect size and precision.
# Positive stat → positively associated with PPT (less pain sensitive).
# Negative stat → negatively associated with PPT (more pain sensitive).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [11] GSEA (Hallmark gene sets) ===")

run_simple_gsea <- function(res_df, pathways_list = NULL, min_size = 15, max_size = 500) {
  if (is.null(pathways_list)) return(NULL)
  stats        <- res_df$stat
  names(stats) <- res_df$Gene_Symbol
  stats        <- stats[!is.na(stats) & !is.na(names(stats))]
  stats        <- stats[!duplicated(names(stats))]
  stats        <- sort(stats, decreasing = TRUE)
  fgsea(pathways = pathways_list, stats = stats,
        minSize = min_size, maxSize = max_size)
}

plot_gsea_dotplot_bulk <- function(gsea_tbl, top_n = GSEA_TOP_N,
                                   pval_col = "pval", sig_cutoff = 0.05,
                                   title = "GSEA: Hallmark Pathways") {
  if (is.null(gsea_tbl) || nrow(gsea_tbl) == 0) return(NULL)

  # fgsea returns a data.table — coerce and drop list columns
  df <- as.data.frame(gsea_tbl)
  df <- df[ , !vapply(df, is.list, logical(1)), drop = FALSE]

  # Select top_n by |NES| among significant
  sig <- df %>% filter(.data[[pval_col]] < sig_cutoff)
  if (nrow(sig) == 0) {
    message("  No significant pathways at p < ", sig_cutoff,
            " — showing top ", top_n, " by pval regardless")
    sig <- df %>% arrange(.data[[pval_col]]) %>% slice_head(n = top_n)
  }

  plot_df <- sig %>%
    mutate(abs_NES = abs(NES)) %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = top_n) %>%
    mutate(
      neglog10p     = -log10(.data[[pval_col]]),
      pathway_label = gsub("_", " ", gsub("^HALLMARK_", "", pathway)),
      pathway_label = factor(pathway_label,
                             levels = pathway_label[order(NES)])
    )

  ggplot(plot_df, aes(x = NES, y = pathway_label,
                      size = neglog10p, color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0) +
    scale_size_continuous(range = c(3, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title    = title,
         subtitle = paste0("Color/position = NES | size = -log10(p) | ",
                           "positive NES = enriched with high PPT (less pain sensitive)"),
         x = "Normalized Enrichment Score (NES)",
         y = NULL,
         size  = "-log10(p)",
         color = "NES") +
    theme_cowplot(14) +
    theme(
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      axis.text.y      = element_text(size = 12),
      legend.text      = element_text(size = 11),
      legend.title     = element_text(size = 12),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
}

# Pathway sets loaded in section 6C-ii; reuse them here
message("  Running V1 GSEA on Hallmark, REACTOME, and GO:BP pathway sets")

gsea_v1_h  <- run_simple_gsea(de_ppt$res_df, pathways_list = pathways_h)
gsea_v1_r  <- run_simple_gsea(de_ppt$res_df, pathways_list = pathways_reactome)
gsea_v1_bp <- run_simple_gsea(de_ppt$res_df, pathways_list = pathways_gobp)

# Save tables
save_gsea_tbl(gsea_v1_h,  "GSEA_V1_Hallmark_PPT.csv", out_dir = DIR_GSEA)
save_gsea_tbl(gsea_v1_r,  "GSEA_V1_Reactome_PPT.csv", out_dir = DIR_GSEA)
save_gsea_tbl(gsea_v1_bp, "GSEA_V1_GOBP_PPT.csv",     out_dir = DIR_GSEA)

# Back-compat alias so downstream sections that reference gsea_res still work
gsea_res <- gsea_v1_h

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {

  # ── Summary table (no list columns) — Hallmark only ─────────────────────────
  gsea_df  <- as.data.frame(gsea_res)
  gsea_out <- gsea_df[ , !vapply(gsea_df, is.list, logical(1)), drop = FALSE] %>%
    arrange(pval)
  write.csv(gsea_out, file.path(DIR_GSEA, "GSEA_Hallmark_PPT.csv"), row.names = FALSE)
  message("  GSEA results saved: GSEA_Hallmark_PPT.csv")

  n_sig <- sum(gsea_out$pval < 0.05, na.rm = TRUE)
  message("  Significant pathways (p<0.05): ", n_sig)

  # ── Leading edge CSV — one row per pathway, genes as comma-separated string ─
  sig_gsea <- gsea_df[!is.na(gsea_df$pval) & gsea_df$pval < 0.05, ]

  if (nrow(sig_gsea) > 0) {
    leading_edge_df <- data.frame(
      pathway        = sig_gsea$pathway,
      pathway_label  = gsub("_", " ", gsub("^HALLMARK_", "", sig_gsea$pathway)),
      NES            = sig_gsea$NES,
      pval           = sig_gsea$pval,
      padj           = sig_gsea$padj,
      n_leading_edge = lengths(sig_gsea$leadingEdge),
      genes_str      = sapply(sig_gsea$leadingEdge, paste, collapse = ", "),
      stringsAsFactors = FALSE
    ) %>% arrange(pval)

    write.csv(leading_edge_df,
              file.path(DIR_GSEA, "GSEA_Hallmark_PPT_leading_edge.csv"),
              row.names = FALSE)
    message("  Leading edge genes saved: GSEA_Hallmark_PPT_leading_edge.csv")
    message("  Significant pathways with leading edge: ", nrow(leading_edge_df))
  } else {
    message("  No significant pathways — leading edge CSV not written")
  }

  # ── Section 11b: Four-color annotated mean-variance plot ────────────────────
  #
  # Overlays four biologically-annotated gene classes on the HVG scatter:
  #
  #   RED    — all top-N_VAR_GENES HVGs (MAD² ranked)
  #   BLUE   — HVG + DESeq2 significant + GSEA leading-edge + negative LFC
  #              (top 20 by MAD² variance; negative LFC = upregulated in low PPT
  #              / high pain-sensitivity patients)
  #   GREEN  — HVG + DESeq2 significant + NOT in any GSEA leading-edge (top 5)
  #   PURPLE — HVG + GSEA leading-edge + NOT DESeq2 significant (top 5)
  #
  # Categories are hierarchical: blue > green > purple > red.
  # Only colored (blue/green/purple) genes are labelled with gene symbols.
  # ─────────────────────────────────────────────────────────────────────────────
  message("=== [11b] Building 4-color annotated mean-variance plot ===")

  if (exists("mv_df_lfc")) {

    # ── Resolve leading-edge ENSG set ─────────────────────────────────────────
    le_ensg_set <- if (nrow(sig_gsea) > 0) {
      unique(gene_sym_map$ENSG[
        gene_sym_map$Gene_Symbol %in% unlist(sig_gsea$leadingEdge)])
    } else {
      character(0)
    }
    message("  Leading-edge unique ENSGs: ", length(le_ensg_set))

    # ── Join DESeq2 annotations (nominal p<0.05 for visualisation breadth) ──────
    # The strict Significant flag (pvalue<0.001) is retained in the DE CSV.
    # At pvalue<0.001 only 6 genes overlap DE+LE; switching to p<0.05 gives
    # 114 overlap genes (22 neg-LFC / 92 pos-LFC), making the categories visible.
    de_ann <- de_ppt$res_df %>%
      select(ENSG, log2FoldChange, pvalue, stat, Significant)

    mv_ann <- mv_df_lfc %>%
      left_join(de_ann, by = "ENSG") %>%
      mutate(
        Gene_Symbol = gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)],
        in_le    = ENSG %in% le_ensg_set,
        de_nom   = !is.na(pvalue) & pvalue < 0.05,
        lfc_neg  = !is.na(log2FoldChange) & log2FoldChange < 0
      )

    message("  HVGs with nominal DE (p<0.05):       ",
            sum(mv_ann$is_top & mv_ann$de_nom, na.rm = TRUE))
    message("  HVGs in GSEA leading-edge:           ",
            sum(mv_ann$is_top & mv_ann$in_le,  na.rm = TRUE))
    message("  HVGs in both, neg LFC (blue pool):   ",
            sum(mv_ann$is_top & mv_ann$de_nom & mv_ann$in_le & mv_ann$lfc_neg, na.rm = TRUE))

    # ── Assign categories (hierarchical: blue > green > purple > red) ─────────
    blue_ensg <- mv_ann %>%
      filter(is_top, de_nom, in_le, lfc_neg) %>%
      arrange(desc(variance)) %>% slice_head(n = 20) %>% pull(ENSG)

    green_ensg <- mv_ann %>%
      filter(is_top, de_nom, !in_le, !ENSG %in% blue_ensg) %>%
      arrange(desc(variance)) %>% slice_head(n = 5) %>% pull(ENSG)

    purple_ensg <- mv_ann %>%
      filter(is_top, in_le, !de_nom, !ENSG %in% blue_ensg) %>%
      arrange(desc(variance)) %>% slice_head(n = 5) %>% pull(ENSG)

    message("  Blue  (DE\u2229GSEA, low PPT): ", length(blue_ensg))
    message("  Green (DE only):         ", length(green_ensg))
    message("  Purple (GSEA only):      ", length(purple_ensg))

    # Factor levels control legend order and drawing order (bg → fg)
    cat_levels <- c("Other", "HVG",
                    "GSEA only — top 5",
                    "DE only — top 5",
                    "DE \u2229 GSEA, low PPT — top 20")

    mv_ann <- mv_ann %>%
      mutate(gene_cat = factor(case_when(
        ENSG %in% blue_ensg   ~ "DE \u2229 GSEA, low PPT — top 20",
        ENSG %in% green_ensg  ~ "DE only — top 5",
        ENSG %in% purple_ensg ~ "GSEA only — top 5",
        is_top                ~ "HVG",
        TRUE                  ~ "Other"
      ), levels = cat_levels)) %>%
      # Draw background genes first so colored dots sit on top
      arrange(gene_cat)

    cat_colors <- c(
      "Other"                           = "#BBBBBB",
      "HVG"                             = "#d84b4b",
      "GSEA only — top 5"               = "#7B2D8B",
      "DE only — top 5"                 = "#1a9641",
      "DE \u2229 GSEA, low PPT — top 20" = "#2166AC"
    )
    cat_alpha  <- c(
      "Other"                           = 0.25,
      "HVG"                             = 0.65,
      "GSEA only — top 5"               = 1.0,
      "DE only — top 5"                 = 1.0,
      "DE \u2229 GSEA, low PPT — top 20" = 1.0
    )
    cat_size   <- c(
      "Other"                           = 0.9,
      "HVG"                             = 1.4,
      "GSEA only — top 5"               = 2.5,
      "DE only — top 5"                 = 2.5,
      "DE \u2229 GSEA, low PPT — top 20" = 2.5
    )

    # Genes to label (blue, green, purple only)
    labeled_df <- mv_ann %>%
      filter(!gene_cat %in% c("Other", "HVG")) %>%
      select(ENSG, mean_cpm, log2_var, gene_cat, Gene_Symbol)

    # Per-category repel colors match point colors
    repel_colors <- c(
      "GSEA only — top 5"               = "#7B2D8B",
      "DE only — top 5"                 = "#1a9641",
      "DE \u2229 GSEA, low PPT — top 20" = "#2166AC"
    )

    mv_4color <- ggplot(mv_ann,
                        aes(x = log2(mean_cpm + 0.01), y = log2_var,
                            color = gene_cat, alpha = gene_cat,
                            size  = gene_cat)) +
      geom_point(stroke = 0) +
      # Separate repel call per color category so label color matches point color
      geom_text_repel(
        data  = filter(labeled_df, gene_cat == "DE \u2229 GSEA, low PPT \u2014 top 20"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors[["DE \u2229 GSEA, low PPT \u2014 top 20"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      geom_text_repel(
        data  = filter(labeled_df, gene_cat == "DE only \u2014 top 5"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors[["DE only \u2014 top 5"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      geom_text_repel(
        data  = filter(labeled_df, gene_cat == "GSEA only \u2014 top 5"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors[["GSEA only \u2014 top 5"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      scale_color_manual(values = cat_colors, name = NULL) +
      scale_alpha_manual(values = cat_alpha,  guide = "none") +
      scale_size_manual( values = cat_size,   guide = "none") +
      guides(color = guide_legend(
        override.aes = list(size = 3, alpha = 1),
        nrow = 3
      )) +
      labs(
        title    = "Mean-variance relationship — min-max LFC",
        subtitle = paste0(
          "HVGs (min-max LFC>=1) | CPM>", CPM_THRESH,
          "\nBlue: DE\u2229GSEA low-PPT (nom p<0.05)  |  Green: DE only  |  Purple: GSEA only"
        ),
        x = "log2(mean CPM)",
        y = "log2(min-max LFC)"
      ) +
      theme_cowplot(12) +
      theme(
        legend.position  = "top",
        legend.text      = element_text(size = 9),
        plot.subtitle    = element_text(size = 9, color = "grey30"),
        plot.margin      = margin(t = 8, r = 15, b = 5, l = 5)
      )

    ggsave(file.path(DIR_GSEA, "00_mean_variance_rawCPM_leadingedge.png"),
           mv_4color, width = 8, height = 6, dpi = 150, bg = "white")
    message("  Saved: 00_mean_variance_rawCPM_leadingedge.png")

  } else {
    message("  mv_df_lfc not found — skipping 4-color mean-variance plot")
  }

  # ── Section 11c: Bridge plot — DESeq2 × GSEA for min-max LFC HVGs ─────────────
  #
  # Scatter of log2(LFC) vs DESeq2 Wald stat for all HVGs (LFC>=1).
  # Bridges three analysis layers in one view:
  #   X-axis : transcriptional variability (min-max LFC rank)
  #   Y-axis : PPT association magnitude + direction (Wald stat)
  #   Color  : GSEA leading-edge membership and DE direction
  #
  # Genes top-right: high variability + positive PPT association (less pain).
  # Genes top-left : high variability + negative PPT association (more pain).
  # Dashed lines at Wald ±2 mark a rough "noteworthy DE" zone.
  # ─────────────────────────────────────────────────────────────────────────────
  message("=== [11c] Bridge plot: DESeq2 \u00d7 GSEA \u00d7 min-max LFC HVGs ===")

  if (exists("mv_df_lfc") && !is.null(de_ppt$res_df)) {

    de_full <- de_ppt$res_df %>%
      select(ENSG, Gene_Symbol, log2FoldChange, stat, pvalue, Significant)

    le_ensg_bridge <- if (exists("le_ensg_set")) le_ensg_set else character(0)

    bridge_df <- mv_df_lfc %>%
      filter(is_top) %>%
      left_join(de_full, by = "ENSG") %>%
      mutate(
        Gene_Symbol = coalesce(
          Gene_Symbol,
          gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)]
        ),
        in_le   = ENSG %in% le_ensg_bridge,
        de_nom  = !is.na(pvalue) & pvalue < 0.05,
        lfc_neg = !is.na(log2FoldChange) & log2FoldChange < 0,
        bridge_cat = factor(case_when(
          in_le &  de_nom &  lfc_neg ~ "LE + DE, low PPT",
          in_le &  de_nom & !lfc_neg ~ "LE + DE, high PPT",
          in_le & !de_nom            ~ "LE only",
         !in_le &  de_nom            ~ "DE only (nom p<0.05)",
          TRUE                       ~ "HVG only"
        ), levels = c("HVG only", "DE only (nom p<0.05)",
                      "LE only", "LE + DE, high PPT", "LE + DE, low PPT"))
      ) %>%
      arrange(bridge_cat)

    bridge_colors <- c(
      "HVG only"             = "#CCCCCC",
      "DE only (nom p<0.05)" = "#1a9641",
      "LE only"              = "#9C6EBA",
      "LE + DE, high PPT"    = "#E08214",
      "LE + DE, low PPT"     = "#2166AC"
    )
    bridge_size <- c(
      "HVG only"             = 0.8,
      "DE only (nom p<0.05)" = 1.6,
      "LE only"              = 1.6,
      "LE + DE, high PPT"    = 2.2,
      "LE + DE, low PPT"     = 2.2
    )
    bridge_alpha <- c(
      "HVG only"             = 0.3,
      "DE only (nom p<0.05)" = 0.8,
      "LE only"              = 0.8,
      "LE + DE, high PPT"    = 1.0,
      "LE + DE, low PPT"     = 1.0
    )

    # Label top genes per category ranked by |stat| × log2(MAD²)
    label_bridge <- bridge_df %>%
      filter(bridge_cat != "HVG only", !is.na(stat)) %>%
      mutate(score = abs(stat) * log2_var) %>%
      group_by(bridge_cat) %>%
      slice_max(score, n = 6) %>%
      ungroup()

    bridge_plot <- ggplot(bridge_df,
                          aes(x = log2_var, y = stat,
                              color = bridge_cat,
                              alpha = bridge_cat,
                              size  = bridge_cat)) +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                 color = "grey60", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "solid",
                 color = "grey40", linewidth = 0.3) +
      geom_point(stroke = 0) +
      geom_text_repel(
        data        = label_bridge,
        aes(label   = Gene_Symbol),
        size        = 2.4, max.overlaps = 25,
        box.padding = 0.35, seed = 42,
        show.legend = FALSE
      ) +
      scale_color_manual(values = bridge_colors, name = NULL) +
      scale_alpha_manual(values = bridge_alpha,  guide = "none") +
      scale_size_manual( values = bridge_size,   guide = "none") +
      guides(color = guide_legend(
        override.aes = list(size = 3, alpha = 1), nrow = 3
      )) +
      labs(
        title    = "DESeq2 \u00d7 GSEA bridge \u2014 min-max LFC HVGs",
        subtitle = paste0(
          "X = variability (log2 min-max LFC) | Y = DESeq2 Wald stat per SD of PPT\n",
          "+stat = high PPT (less pain sensitive) | \u2212stat = low PPT (more pain sensitive)\n",
          "Dashed = Wald \u00b12 | LE = Hallmark leading-edge genes | nom p<0.05"
        ),
        x = "log2(min-max LFC)",
        y = "DESeq2 Wald statistic"
      ) +
      theme_cowplot(12) +
      theme(
        legend.position = "top",
        legend.text     = element_text(size = 9),
        plot.subtitle   = element_text(size = 8, color = "grey30", lineheight = 1.3),
        plot.margin     = margin(t = 8, r = 15, b = 5, l = 5)
      )

    ggsave(file.path(DIR_GSEA, "00e_bridge_DESeq2_GSEA_HVG.png"),
           bridge_plot, width = 8, height = 7, dpi = 150, bg = "white")
    message("  Saved: 00e_bridge_DESeq2_GSEA_HVG.png")

  } else {
    message("  Skipping bridge plot \u2014 mv_df_lfc or de_ppt not available")
  }

  # \u2500\u2500 Section 11d: 4-color annotated mean-variance plot \u2014 P10-P90 LFC HVGs \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  message("=== [11d] Building 4-color annotated mean-variance plot (P10-P90 LFC) ===")

  if (exists("mv_df_lfc_pct")) {

    le_ensg_set_d <- if (nrow(sig_gsea) > 0) {
      unique(gene_sym_map$ENSG[
        gene_sym_map$Gene_Symbol %in% unlist(sig_gsea$leadingEdge)])
    } else {
      character(0)
    }

    de_ann_d <- de_ppt$res_df %>%
      select(ENSG, log2FoldChange, pvalue, stat, Significant)

    mv_ann_d <- mv_df_lfc_pct %>%
      left_join(de_ann_d, by = "ENSG") %>%
      mutate(
        Gene_Symbol = gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)],
        in_le   = ENSG %in% le_ensg_set_d,
        de_nom  = !is.na(pvalue) & pvalue < 0.05,
        lfc_neg = !is.na(log2FoldChange) & log2FoldChange < 0
      )

    blue_ensg_d <- mv_ann_d %>%
      filter(is_top, de_nom, in_le, lfc_neg) %>%
      arrange(desc(variance)) %>% slice_head(n = 20) %>% pull(ENSG)
    green_ensg_d <- mv_ann_d %>%
      filter(is_top, de_nom, !in_le, !ENSG %in% blue_ensg_d) %>%
      arrange(desc(variance)) %>% slice_head(n = 5) %>% pull(ENSG)
    purple_ensg_d <- mv_ann_d %>%
      filter(is_top, in_le, !de_nom, !ENSG %in% blue_ensg_d) %>%
      arrange(desc(variance)) %>% slice_head(n = 5) %>% pull(ENSG)

    cat_levels_d <- c("Other", "HVG",
                      "GSEA only \u2014 top 5",
                      "DE only \u2014 top 5",
                      "DE \u2229 GSEA, low PPT \u2014 top 20")

    mv_ann_d <- mv_ann_d %>%
      mutate(gene_cat = factor(case_when(
        ENSG %in% blue_ensg_d   ~ "DE \u2229 GSEA, low PPT \u2014 top 20",
        ENSG %in% green_ensg_d  ~ "DE only \u2014 top 5",
        ENSG %in% purple_ensg_d ~ "GSEA only \u2014 top 5",
        is_top                  ~ "HVG",
        TRUE                    ~ "Other"
      ), levels = cat_levels_d)) %>%
      arrange(gene_cat)

    cat_colors_d <- c(
      "Other"                           = "#BBBBBB",
      "HVG"                             = "#d84b4b",
      "GSEA only \u2014 top 5"          = "#7B2D8B",
      "DE only \u2014 top 5"            = "#1a9641",
      "DE \u2229 GSEA, low PPT \u2014 top 20" = "#2166AC"
    )
    cat_alpha_d  <- c(
      "Other"                           = 0.25,
      "HVG"                             = 0.65,
      "GSEA only \u2014 top 5"          = 1.0,
      "DE only \u2014 top 5"            = 1.0,
      "DE \u2229 GSEA, low PPT \u2014 top 20" = 1.0
    )
    cat_size_d   <- c(
      "Other"                           = 0.9,
      "HVG"                             = 1.4,
      "GSEA only \u2014 top 5"          = 2.5,
      "DE only \u2014 top 5"            = 2.5,
      "DE \u2229 GSEA, low PPT \u2014 top 20" = 2.5
    )

    labeled_df_d <- mv_ann_d %>%
      filter(!gene_cat %in% c("Other", "HVG")) %>%
      select(ENSG, mean_cpm, log2_var, gene_cat, Gene_Symbol)

    repel_colors_d <- c(
      "GSEA only \u2014 top 5"          = "#7B2D8B",
      "DE only \u2014 top 5"            = "#1a9641",
      "DE \u2229 GSEA, low PPT \u2014 top 20" = "#2166AC"
    )

    mv_4color_d <- ggplot(mv_ann_d,
                          aes(x = log2(mean_cpm + 0.01), y = log2_var,
                              color = gene_cat, alpha = gene_cat,
                              size  = gene_cat)) +
      geom_point(stroke = 0) +
      geom_text_repel(
        data  = filter(labeled_df_d, gene_cat == "DE \u2229 GSEA, low PPT \u2014 top 20"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors_d[["DE \u2229 GSEA, low PPT \u2014 top 20"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      geom_text_repel(
        data  = filter(labeled_df_d, gene_cat == "DE only \u2014 top 5"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors_d[["DE only \u2014 top 5"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      geom_text_repel(
        data  = filter(labeled_df_d, gene_cat == "GSEA only \u2014 top 5"),
        aes(x = log2(mean_cpm + 0.01), y = log2_var, label = Gene_Symbol),
        inherit.aes = FALSE,
        color = repel_colors_d[["GSEA only \u2014 top 5"]],
        size  = 2.5, max.overlaps = 30, box.padding = 0.35, seed = 42
      ) +
      scale_color_manual(values = cat_colors_d, name = NULL) +
      scale_alpha_manual(values = cat_alpha_d,  guide = "none") +
      scale_size_manual( values = cat_size_d,   guide = "none") +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 3)) +
      labs(
        title    = "Mean-variance relationship \u2014 P10-P90 LFC",
        subtitle = paste0(
          "HVGs (P10-P90 LFC>=1) | CPM>", CPM_THRESH,
          "\nBlue: DE\u2229GSEA low-PPT (nom p<0.05)  |  Green: DE only  |  Purple: GSEA only"
        ),
        x = "log2(mean CPM)",
        y = "log2(P10-P90 LFC)"
      ) +
      theme_cowplot(12) +
      theme(
        legend.position  = "top",
        legend.text      = element_text(size = 9),
        plot.subtitle    = element_text(size = 9, color = "grey30"),
        plot.margin      = margin(t = 8, r = 15, b = 5, l = 5)
      )

    ggsave(file.path(DIR_GSEA, "00_mean_variance_P10P90LFC_leadingedge.png"),
           mv_4color_d, width = 8, height = 6, dpi = 150, bg = "white")
    message("  Saved: 00_mean_variance_P10P90LFC_leadingedge.png")

  } else {
    message("  mv_df_lfc_pct not found \u2014 skipping P10-P90 4-color mean-variance plot")
  }

  # \u2500\u2500 Section 11e: Bridge plot \u2014 DESeq2 \u00d7 GSEA for P10-P90 LFC HVGs \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  message("=== [11e] Bridge plot: DESeq2 \u00d7 GSEA \u00d7 P10-P90 LFC HVGs ===")

  if (exists("mv_df_lfc_pct") && !is.null(de_ppt$res_df)) {

    de_full_e <- de_ppt$res_df %>%
      select(ENSG, Gene_Symbol, log2FoldChange, stat, pvalue, Significant)

    le_ensg_e <- if (exists("le_ensg_set")) le_ensg_set else character(0)

    bridge_df_e <- mv_df_lfc_pct %>%
      filter(is_top) %>%
      left_join(de_full_e, by = "ENSG") %>%
      mutate(
        Gene_Symbol = coalesce(
          Gene_Symbol,
          gene_sym_map$Gene_Symbol[match(ENSG, gene_sym_map$ENSG)]
        ),
        in_le   = ENSG %in% le_ensg_e,
        de_nom  = !is.na(pvalue) & pvalue < 0.05,
        lfc_neg = !is.na(log2FoldChange) & log2FoldChange < 0,
        bridge_cat = factor(case_when(
          in_le &  de_nom &  lfc_neg ~ "LE + DE, low PPT",
          in_le &  de_nom & !lfc_neg ~ "LE + DE, high PPT",
          in_le & !de_nom            ~ "LE only",
         !in_le &  de_nom            ~ "DE only (nom p<0.05)",
          TRUE                       ~ "HVG only"
        ), levels = c("HVG only", "DE only (nom p<0.05)",
                      "LE only", "LE + DE, high PPT", "LE + DE, low PPT"))
      ) %>%
      arrange(bridge_cat)

    bridge_colors_e <- c(
      "HVG only"             = "#CCCCCC",
      "DE only (nom p<0.05)" = "#1a9641",
      "LE only"              = "#9C6EBA",
      "LE + DE, high PPT"    = "#E08214",
      "LE + DE, low PPT"     = "#2166AC"
    )
    bridge_size_e <- c(
      "HVG only"             = 0.8,
      "DE only (nom p<0.05)" = 1.6,
      "LE only"              = 1.6,
      "LE + DE, high PPT"    = 2.2,
      "LE + DE, low PPT"     = 2.2
    )
    bridge_alpha_e <- c(
      "HVG only"             = 0.3,
      "DE only (nom p<0.05)" = 0.8,
      "LE only"              = 0.8,
      "LE + DE, high PPT"    = 1.0,
      "LE + DE, low PPT"     = 1.0
    )

    label_bridge_e <- bridge_df_e %>%
      filter(bridge_cat != "HVG only", !is.na(stat)) %>%
      mutate(score = abs(stat) * log2_var) %>%
      group_by(bridge_cat) %>%
      slice_max(score, n = 6) %>%
      ungroup()

    bridge_plot_e <- ggplot(bridge_df_e,
                            aes(x = log2_var, y = stat,
                                color = bridge_cat,
                                alpha = bridge_cat,
                                size  = bridge_cat)) +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                 color = "grey60", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "solid",
                 color = "grey40", linewidth = 0.3) +
      geom_point(stroke = 0) +
      geom_text_repel(
        data        = label_bridge_e,
        aes(label   = Gene_Symbol),
        size        = 2.4, max.overlaps = 25,
        box.padding = 0.35, seed = 42,
        show.legend = FALSE
      ) +
      scale_color_manual(values = bridge_colors_e, name = NULL) +
      scale_alpha_manual(values = bridge_alpha_e,  guide = "none") +
      scale_size_manual( values = bridge_size_e,   guide = "none") +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 3)) +
      labs(
        title    = "DESeq2 \u00d7 GSEA bridge \u2014 P10-P90 LFC HVGs",
        subtitle = paste0(
          "X = variability (log2 P10-P90 LFC) | Y = DESeq2 Wald stat per SD of PPT\n",
          "+stat = high PPT (less pain sensitive) | \u2212stat = low PPT (more pain sensitive)\n",
          "Dashed = Wald \u00b12 | LE = Hallmark leading-edge genes | nom p<0.05"
        ),
        x = "log2(P10-P90 LFC)",
        y = "DESeq2 Wald statistic"
      ) +
      theme_cowplot(12) +
      theme(
        legend.position = "top",
        legend.text     = element_text(size = 9),
        plot.subtitle   = element_text(size = 8, color = "grey30", lineheight = 1.3),
        plot.margin     = margin(t = 8, r = 15, b = 5, l = 5)
      )

    ggsave(file.path(DIR_GSEA, "00e_bridge_DESeq2_GSEA_P10P90LFC.png"),
           bridge_plot_e, width = 8, height = 7, dpi = 150, bg = "white")
    message("  Saved: 00e_bridge_DESeq2_GSEA_P10P90LFC.png")

  } else {
    message("  Skipping P10-P90 bridge plot \u2014 mv_df_lfc_pct or de_ppt not available")
  }

  # Dotplots — one per pathway collection
  dp_v1_h  <- plot_gsea_dotplot_generic(gsea_v1_h,  top_n = GSEA_TOP_N,
    title = "V1 GSEA — Hallmark | V1 ~ PPT (continuous)",
    strip_prefix = "^HALLMARK_")
  dp_v1_r  <- plot_gsea_dotplot_generic(gsea_v1_r,  top_n = GSEA_TOP_N,
    title = "V1 GSEA — REACTOME | V1 ~ PPT (continuous)",
    strip_prefix = "^REACTOME_")
  dp_v1_bp <- plot_gsea_dotplot_generic(gsea_v1_bp, top_n = GSEA_TOP_N,
    title = "V1 GSEA — GO:BP | V1 ~ PPT (continuous)",
    strip_prefix = "^GOBP_")

  for (dp_pair in list(
    list(dp = dp_v1_h,  fn = "04_GSEA_V1_Hallmark_PPT_dotplot.png",  w = 10, h = 8),
    list(dp = dp_v1_r,  fn = "04_GSEA_V1_Reactome_PPT_dotplot.png",  w = 14, h = 12),
    list(dp = dp_v1_bp, fn = "04_GSEA_V1_GOBP_PPT_dotplot.png",      w = 14, h = 12)
  )) {
    if (!is.null(dp_pair$dp)) {
      ggsave(file.path(DIR_GSEA, dp_pair$fn), dp_pair$dp,
             width = dp_pair$w, height = dp_pair$h, dpi = 150, bg = "white")
      message("  Saved: ", dp_pair$fn)
    }
  }
} else {
  message("  GSEA returned no results — check stat column in DE output")
}

# ─────────────────────────────────────────────────────────────────────────────
# 11.5  DESEQ2 + GSEA — V2 SAMPLES ~ qst_ppt_tr_avg_v2
#
# Uses Visit-2-specific PPT values (qst_ppt_tr_avg_v2) as the predictor,
# matching V2 expression to concurrent pain phenotype at that timepoint.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [11.5] DESeq2: V2 ~ qst_ppt_tr_avg_v2 (continuous PPT) ===")

meta_v2 <- meta %>%
  filter(Visit == "Visit2", !IsPerlman, !is.na(qst_ppt_tr_avg_v2)) %>%
  mutate(Visit = droplevels(Visit))

message("  V2 samples with V2 PPT data: ", nrow(meta_v2))

if (nrow(meta_v2) >= 5) {
  counts_v2 <- raw[, rownames(meta_v2)]

  de_ppt_v2 <- run_deseq2(
    counts           = counts_v2,
    col_data         = meta_v2,
    design_formula   = ~ qst_ppt_tr_avg_v2,
    contrast_or_name = "qst_ppt_tr_avg_v2",
    label            = "V2_continuous_PPT",
    results_type     = "name",
    out_dir          = DIR_DESEQ2
  )

  # ── V2 example gene scatter plots ──
  message("=== [11.5a] V2 example gene visualisation ===")
  cpm_v2 <- cpm_filt[, intersect(rownames(meta_v2), colnames(cpm_filt))]

  top_genes_v2 <- de_ppt_v2$res_df %>%
    filter(!is.na(pvalue)) %>%
    filter(Significant) %>%
    filter(ENSG %in% rownames(cpm_filt)) %>%
    slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
    pull(ENSG)

  if (length(top_genes_v2) == 0) {
    top_genes_v2 <- de_ppt_v2$res_df %>%
      filter(!is.na(pvalue)) %>%
      filter(ENSG %in% rownames(cpm_filt)) %>%
      slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
      pull(ENSG)
    v2_suffix <- "_nominal"
    message("  No significant V2 genes — using top by nominal p-value")
  } else {
    v2_suffix <- ""
  }

  plot_list_v2 <- lapply(top_genes_v2, function(g) {
    sym  <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    lfc  <- round(de_ppt_v2$res_df$log2FoldChange[de_ppt_v2$res_df$ENSG == g], 3)
    pnom <- signif(de_ppt_v2$res_df$pvalue[de_ppt_v2$res_df$ENSG == g], 3)
    df   <- data.frame(
      PPT  = meta_v2[colnames(cpm_v2), "qst_ppt_raw"],
      expr = as.numeric(cpm_v2[g, ]),
      stringsAsFactors = FALSE
    )
    ggplot(df, aes(PPT, expr)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      labs(title    = sym,
           subtitle = paste0("LFC/SD=", lfc, "  p=", pnom),
           x = "PPT (raw)", y = "CPM") +
      theme_cowplot(10) +
      theme(plot.title    = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 7))
  })

  if (length(plot_list_v2) > 0) {
    n_cols_v2 <- min(4, length(plot_list_v2))
    n_rows_v2 <- ceiling(length(plot_list_v2) / n_cols_v2)
    out_fn_v2 <- paste0("DEseq_example_genes_V2", v2_suffix, ".png")
    png(file.path(DIR_DESEQ2, out_fn_v2),
        width = n_cols_v2 * 400, height = n_rows_v2 * 380, res = 120, type = "cairo")
    print(plot_grid(plotlist = plot_list_v2, ncol = n_cols_v2))
    dev.off()
    message("  Saved: ", out_fn_v2, " (", length(plot_list_v2), " genes)")
  }

  message("=== [11.5b] GSEA on V2 DESeq2 results ===")
  gsea_v2_h  <- run_simple_gsea(de_ppt_v2$res_df, pathways_list = pathways_h)
  gsea_v2_r  <- run_simple_gsea(de_ppt_v2$res_df, pathways_list = pathways_reactome)
  gsea_v2_bp <- run_simple_gsea(de_ppt_v2$res_df, pathways_list = pathways_gobp)

  save_gsea_tbl(gsea_v2_h,  "GSEA_V2_Hallmark_PPT.csv", out_dir = DIR_GSEA)
  save_gsea_tbl(gsea_v2_r,  "GSEA_V2_Reactome_PPT.csv", out_dir = DIR_GSEA)
  save_gsea_tbl(gsea_v2_bp, "GSEA_V2_GOBP_PPT.csv",     out_dir = DIR_GSEA)

  for (dp_pair in list(
    list(gs = gsea_v2_h,  fn = "04_GSEA_V2_Hallmark_PPT_dotplot.png",  w = 10, h = 8,
         title = "V2 GSEA — Hallmark | V2 ~ V2 PPT (continuous)", pfx = "^HALLMARK_"),
    list(gs = gsea_v2_r,  fn = "04_GSEA_V2_Reactome_PPT_dotplot.png",  w = 14, h = 12,
         title = "V2 GSEA — REACTOME | V2 ~ V2 PPT (continuous)", pfx = "^REACTOME_"),
    list(gs = gsea_v2_bp, fn = "04_GSEA_V2_GOBP_PPT_dotplot.png",      w = 14, h = 12,
         title = "V2 GSEA — GO:BP | V2 ~ V2 PPT (continuous)",   pfx = "^GOBP_")
  )) {
    dp <- plot_gsea_dotplot_generic(dp_pair$gs, top_n = GSEA_TOP_N,
                                     title = dp_pair$title,
                                     strip_prefix = dp_pair$pfx)
    if (!is.null(dp)) {
      ggsave(file.path(DIR_GSEA, dp_pair$fn), dp,
             width = dp_pair$w, height = dp_pair$h, dpi = 150, bg = "white")
      message("  Saved: ", dp_pair$fn)
    }
  }
} else {
  message("  Skipping V2 DESeq2/GSEA — fewer than 5 V2 samples with PPT data")
  de_ppt_v2 <- NULL
}

# ─────────────────────────────────────────────────────────────────────────────
# 11.6  DESEQ2 + GSEA — V1 SAMPLES ~ PROMIS DIFFERENCE (V2 - V1)
#
# Evaluates baseline V1 expression against the change in PROMIS score
# from Visit 1 to Visit 2. Restricted to V1 PROMO samples whose subject has
# both a v1 and a v2 PROMIS value (promis_diff is non-NA).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [11.6] DESeq2: V1 Baseline ~ promis_diff (V2 - V1) ===")

meta_promis_diff <- meta %>%
  filter(Visit == "Visit1", !IsPerlman, !is.na(promis_diff)) %>%
  mutate(Visit = droplevels(Visit))

message("  V1 samples with both V1 and V2 PROMIS data: ", nrow(meta_promis_diff))

if (nrow(meta_promis_diff) >= 5) {
  counts_promis <- raw[, rownames(meta_promis_diff)]

  # Scale the difference so LFC represents change per 1 Standard Deviation
  meta_promis_diff$promis_diff_scaled <- as.numeric(scale(meta_promis_diff$promis_diff))

  de_promis <- run_deseq2(
    counts           = counts_promis,
    col_data         = meta_promis_diff,
    design_formula   = ~ promis_diff_scaled,
    contrast_or_name = "promis_diff_scaled",
    label            = "V1_continuous_PROMIS_diff",
    results_type     = "name",
    out_dir          = DIR_DESEQ2
  )

  # ── PROMIS Diff example gene scatter plots ──
  message("=== [11.6a] PROMIS Diff example gene visualisation ===")
  cpm_promis <- cpm_filt[, intersect(rownames(meta_promis_diff), colnames(cpm_filt))]

  top_genes_promis <- de_promis$res_df %>%
    filter(!is.na(pvalue)) %>%
    filter(Significant) %>%
    filter(ENSG %in% rownames(cpm_filt)) %>%
    slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
    pull(ENSG)

  if (length(top_genes_promis) == 0) {
    top_genes_promis <- de_promis$res_df %>%
      filter(!is.na(pvalue)) %>%
      filter(ENSG %in% rownames(cpm_filt)) %>%
      slice_min(pvalue, n = N_EXAMPLE_GENES) %>%
      pull(ENSG)
    promis_suffix <- "_nominal"
    message("  No significant PROMIS diff genes — using top by nominal p-value")
  } else {
    promis_suffix <- ""
  }

  plot_list_promis <- lapply(top_genes_promis, function(g) {
    sym  <- gene_sym_map$Gene_Symbol[gene_sym_map$ENSG == g]
    lfc  <- round(de_promis$res_df$log2FoldChange[de_promis$res_df$ENSG == g], 3)
    pnom <- signif(de_promis$res_df$pvalue[de_promis$res_df$ENSG == g], 3)
    df   <- data.frame(
      PROMIS_diff = meta_promis_diff[colnames(cpm_promis), "promis_diff"],
      expr        = as.numeric(cpm_promis[g, ]),
      stringsAsFactors = FALSE
    )
    ggplot(df, aes(PROMIS_diff, expr)) +
      geom_point(color = "#2c7bb6", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#d7191c",
                  fill = "#d7191c", alpha = 0.15) +
      labs(title    = sym,
           subtitle = paste0("LFC/SD=", lfc, "  p=", pnom),
           x = "PROMIS (V2 - V1)", y = "CPM (V1 Baseline)") +
      theme_cowplot(10) +
      theme(plot.title    = element_text(face = "bold", size = 10),
            plot.subtitle = element_text(size = 7))
  })

  if (length(plot_list_promis) > 0) {
    n_cols_p <- min(4, length(plot_list_promis))
    n_rows_p <- ceiling(length(plot_list_promis) / n_cols_p)
    out_fn_p <- paste0("DEseq_example_genes_PROMIS_diff", promis_suffix, ".png")
    png(file.path(DIR_DESEQ2, out_fn_p),
        width = n_cols_p * 400, height = n_rows_p * 380, res = 120, type = "cairo")
    print(plot_grid(plotlist = plot_list_promis, ncol = n_cols_p))
    dev.off()
    message("  Saved: ", out_fn_p, " (", length(plot_list_promis), " genes)")
  }

  # ── GSEA for PROMIS Diff ──
  message("=== [11.6b] GSEA on PROMIS diff DESeq2 results ===")
  gsea_promis_h  <- run_simple_gsea(de_promis$res_df, pathways_list = pathways_h)
  save_gsea_tbl(gsea_promis_h, "GSEA_PROMIS_diff_Hallmark.csv", out_dir = DIR_GSEA)

  dp_promis <- plot_gsea_dotplot_generic(gsea_promis_h, top_n = GSEA_TOP_N,
                                   title = "V1 GSEA — Hallmark | V1 ~ PROMIS Diff",
                                   strip_prefix = "^HALLMARK_")
  if (!is.null(dp_promis)) {
    ggsave(file.path(DIR_GSEA, "04_GSEA_PROMIS_diff_Hallmark_dotplot.png"), dp_promis,
           width = 10, height = 8, dpi = 150, bg = "white")
    message("  Saved: 04_GSEA_PROMIS_diff_Hallmark_dotplot.png")
  }
} else {
  message("  Skipping PROMIS Diff DESeq2/GSEA — fewer than 5 V1 samples with V1 & V2 PROMIS data")
  de_promis <- NULL
}

# ─────────────────────────────────────────────────────────────────────────────
# 12. SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
message("=== [12] Summary ===")

summary_lines <- c(
  "PROMOTER Bulk RNA-seq Pipeline Summary",
  paste("Run date:", Sys.time()),
  "",
  "── Input ──────────────────────────────────────────────────",
  paste("  Total samples in matrix:              ", ncol(raw)),
  paste("  Genes raw:                            ", nrow(raw)),
  paste("  Genes after CPM>", CPM_THRESH,
        " in >=", round(MIN_PCT_SAMPLES*100), "% samples:", nrow(raw_filt)),
  paste("  Pain metadata subjects:               ", nrow(pain_raw)),
  paste("  V1 samples with PPT data:             ", nrow(meta_v1)),
  "",
  "── Variable genes ─────────────────────────────────────────",
  "  5 HVG methods compared:",
  "    1. VST variance (rowVars on VST values, all PROMO)",
  "    2. CPM rowVar (row variance on raw CPM, all PROMO)",
  "    3. min-max LFC (log2 CPM+1 range, V1 PROMO, LFC>=1)",
  "    4. P5-P95 LFC (5th-95th pct LFC, V1 PROMO, LFC>=1)",
  "    5. P10-P90 LFC (10th-90th pct LFC, V1 PROMO, LFC>=1)",
  paste("  Top variable genes (VST/rowVar cap):  ", N_VAR_GENES),
  paste("  min-max LFC HVG count (LFC>=1):       ", sum(gene_vars_lfc >= 1)),
  paste("  P5-P95 LFC HVG count (LFC>=1):        ", sum(gene_vars_lfc_p5p95 >= 1)),
  paste("  P10-P90 LFC HVG count (LFC>=1):       ", sum(gene_vars_lfc_pct >= 1)),
  paste("  K-means k values:                     ", paste(K_VALUES, collapse=",")),
  "",
  "── DESeq2 ─────────────────────────────────────────────────",
  "  V1 model: ~ qst_ppt_tr_avg_v1 (centered + scaled)",
  "  V2 model: ~ qst_ppt_tr_avg_v2 (centered + scaled)",
  paste0("  Internal filter: count >= ", DESEQ_MIN_COUNT, " in >= ",
         round(MIN_PCT_SAMPLES*100), "% of samples"),
  "  LFC interpretation: change per 1 SD increase in PPT",
  "  +LFC = higher expression with higher PPT (less pain sensitive)",
  paste("  Significant genes (pvalue<0.001):",
        sum(de_ppt$res_df$Significant, na.rm=TRUE)),
  "",
  "── Output files ───────────────────────────────────────────",
  "  sample_metadata.csv",
  "  variance_variable_genes_all.csv    (all 5 HVG methods)",
  "  HVG_list_VST.csv / HVG_list_CPM_CV.csv / HVG_list_minMaxLFC.csv / HVG_list_P5P95LFC.csv / HVG_list_P10P90LFC.csv",
  paste0("  kmeans_k", max(K_VALUES), "_gene_clusters.csv"),
  "  00_filter_visualization.png          (2x2 grid: incl. raw count panel)",
  "  00_mean_variance_comparison_grid.png (1x5: VST, CPM var, min-max LFC, P5-P95 LFC, P10-P90 LFC)",
  "  00_mean_variance_LFC_comparator.png  (3-panel: min-max LFC vs P5-P95 LFC vs P10-P90 LFC)",
  "  00_venn_LFC_3way.png                 (3-way Venn: min-max / P5-P95 / P10-P90)",
  "  00_lfc_outlier_fraction_diagnostic.png (explains min-max vs P10-P90 silhouette paradox)",
  "  01a_heatmap_genes_x_samples.png",
  "  01b_correlation_PROMO_prefilter.png  (all genes, raw CPM Pearson relative scale)",
  "  01b_correlation_PROMO_postfilter.png (filtered genes, raw CPM Pearson relative scale)",
  "  02_heatmap_kmeans_k{k}.png",
  "  02b_silhouette_k{k}_{method}.png     (k>=2 only)",
  "  02c_ORA_CP_facet_k{k}_{method}_{Hallmark|GOBP}.png  (CP ORA facet grid per k)",
  "  02e_cluster_expression_violin_k{k}.png  (absolute log2CPM per cluster — heatmap diagnostic)",
  "  03_example_genes_PPT_scatter.png",
  "  DESeq2_V1_continuous_PPT.csv + Volcano_V1_continuous_PPT.png",
  "  GSEA_Hallmark_PPT.csv",
  "  GSEA_Hallmark_PPT_leading_edge.csv (one row per sig pathway, genes as comma-sep string)",
  "  04_GSEA_Hallmark_PPT_dotplot.png",
  "  00_mean_variance_LFC_leadingedge.png     (min-max LFC: leading-edge HVGs labelled)",
  "  00_mean_variance_P10P90LFC_leadingedge.png (P10-P90 LFC: leading-edge HVGs labelled)",
  "  00e_bridge_DESeq2_GSEA_HVG.png        (min-max LFC bridge plot)",
  "  00e_bridge_DESeq2_GSEA_P10P90LFC.png  (P10-P90 LFC bridge plot)"
)
writeLines(summary_lines, file.path(OUT_DIR, "pipeline_summary.txt"))
cat(paste(summary_lines, collapse="\n"), "\n")

message("\n=== PIPELINE COMPLETE — outputs in: ", OUT_DIR, " ===")
