###############################################################################
# PROMOTER — V1 → paired-V2 GENE-CLUSTER ROBUSTNESS
# Spin-off of bulk_normalize_cluster_DE.R
#
# Question: is the V1 gene clustering (k=6, LFC_P5P95 z-score pass) reproducible
# at Visit 2? Restricted to subjects with a PAIRED V1 + V2 RNA-seq sample, so
# both visits use the identical subject set (apples-to-apples).
#
# Two robustness methods (both produced):
#   (A) PROJECT  — keep the V1 cluster centroids fixed; assign each gene to its
#                  nearest V1 centroid using its V2 z-score profile. Cluster
#                  identities stay aligned automatically (shared centroids).
#   (B) RECLUSTER— run a fresh k=6 k-means on the V2 z-score matrix, then align
#                  the V2 labels to V1 by maximum gene overlap (Hungarian).
#                  Adjusted Rand Index (label-invariant) is the headline metric.
#
# Sankey (ggalluvial) per method: V1 cluster -> V2 cluster, flows = genes,
# each cluster its own color.
#
# Emulated parent plots (V1 baseline + V2): silhouette (02b), gene-PCA scatter
# (02f), per-cluster ORA facet grid (02c), and the Z-score-per-patient CSV +
# cluster-mean line plot (02i) + per-cluster vs metadata correlation matrix.
#
# Inputs (same working dir as parent, on the cluster):
#   - kmeans_k6_gene_clusters_LFC_P5P95_ZSCORE.csv   V1 gene->cluster map
#   - CPM_PROMOTOR_BulkSeq_20260216_final_count.txt  CPM (Gene_Symbol col first)
#   - example_pain_metadata.csv
###############################################################################
.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(RColorBrewer)
  library(cowplot)
  library(pheatmap)      # clustered heatmaps (parent 02_heatmap_kmeans style)
  library(cluster)       # silhouette()
  library(ggalluvial)    # Sankey / alluvial
  library(ggrepel)
  library(scales)
  library(msigdbr)       # ORA pathway sets
})
HAS_CP   <- requireNamespace("clusterProfiler", quietly = TRUE)
if (HAS_CP) suppressPackageStartupMessages(library(clusterProfiler)) else
  message("  NOTE: clusterProfiler not found — per-cluster ORA will be skipped")
HAS_CLUE <- requireNamespace("clue", quietly = TRUE)  # optimal label alignment

set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# 0. PATHS  ── edit these
# ─────────────────────────────────────────────────────────────────────────────
CLUSTER_CSV <- "kmeans_k6_gene_clusters_LFC_P5P95_ZSCORE.csv"
CPM_FILE    <- "CPM_PROMOTOR_BulkSeq_20260216_final_count.txt"
PAIN_META   <- "example_pain_metadata.csv"

SUF      <- "LFC_P5P95"   # label tag carried through filenames (matches parent)
KM_SEED  <- 40            # parent uses set.seed(40) right before kmeans()

OUT_DIR  <- paste0("gene_cluster_robustness_", format(Sys.Date(), "%Y%m%d"))
DIR_SANKEY    <- file.path(OUT_DIR, "Sankey")
DIR_SIL       <- file.path(OUT_DIR, "Silhouette")
DIR_PCA       <- file.path(OUT_DIR, "GenePCA")
DIR_HEATMAP   <- file.path(OUT_DIR, "Heatmaps")
DIR_CLUSTMETA <- file.path(OUT_DIR, "ClusterMeta")
DIR_ORA       <- file.path(OUT_DIR, "ORA")
for (d in c(OUT_DIR, DIR_SANKEY, DIR_SIL, DIR_PCA, DIR_HEATMAP, DIR_CLUSTMETA, DIR_ORA))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# helpers (copied from parent so the spin-off stays self-contained)
# ─────────────────────────────────────────────────────────────────────────────
min_max_norm <- function(x, min_val = -1, max_val = 1) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  min_val + (x - rng[1]) / (rng[2] - rng[1]) * (max_val - min_val)
}

# Adjusted Rand Index — label-invariant agreement between two clusterings.
adj_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) sum(choose(x, 2))
  idx  <- comb2(as.vector(tab))
  e    <- comb2(rowSums(tab)) * comb2(colSums(tab)) / choose(sum(tab), 2)
  maxi <- (comb2(rowSums(tab)) + comb2(colSums(tab))) / 2
  if ((maxi - e) == 0) return(NA_real_)
  (idx - e) / (maxi - e)
}

# Align `new` integer labels (1..K) to `ref` labels by maximising gene overlap.
# Returns a length-K vector: aligned_label_for[new_cluster_i].
align_labels <- function(ref, new, K) {
  m <- table(factor(new, levels = seq_len(K)), factor(ref, levels = seq_len(K)))
  m <- matrix(as.numeric(m), nrow = K, ncol = K)
  if (HAS_CLUE) {
    as.integer(clue::solve_LSAP(m, maximum = TRUE))   # m[i, mapping[i]] maximised
  } else {                                            # greedy fallback
    mapping <- integer(K); avail <- rep(TRUE, K); rows <- rep(TRUE, K)
    for (step in seq_len(K)) {
      mm <- m; mm[!rows, ] <- -1; mm[, !avail] <- -1
      ix <- which(mm == max(mm), arr.ind = TRUE)[1, ]
      mapping[ix[1]] <- ix[2]; rows[ix[1]] <- FALSE; avail[ix[2]] <- FALSE
    }
    mapping
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
message("=== [1] Load data ===")

clust_v1 <- read.csv(CLUSTER_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("ENSG", "Gene_Symbol", "Cluster") %in% colnames(clust_v1)))
K        <- max(clust_v1$Cluster)
message("  V1 clustering: ", nrow(clust_v1), " genes | k = ", K)

cpm_raw <- read.delim(CPM_FILE, row.names = 1, check.names = FALSE)
gene_sym_map <- data.frame(
  ENSG        = rownames(cpm_raw),
  Gene_Symbol = cpm_raw[["Gene_Symbol"]],
  stringsAsFactors = FALSE
)
cpm <- cpm_raw[, -1, drop = FALSE]   # drop Gene_Symbol column → numeric CPM matrix
cpm <- as.matrix(cpm)
message("  CPM matrix: ", nrow(cpm), " genes x ", ncol(cpm), " samples")

pain_raw <- read.csv(PAIN_META, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
pain_raw$SubjectID <- gsub("_", "", pain_raw$study_id)   # PROMO_13 → PROMO13
rownames(pain_raw) <- pain_raw$SubjectID

# ─────────────────────────────────────────────────────────────────────────────
# 2. IDENTIFY PAIRED V1 + V2 RNA-seq SAMPLES (from CPM column names)
#    Sample IDs look like PROMO37_V2_CD11b_S12 → SubjectID + Visit.
#    Pairing for CLUSTERING is defined by having a V2 EXPRESSION column, not by
#    metadata v2 values (which are clinical follow-up and may not coincide).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [2] Identify paired V1+V2 expression samples ===")

samp_info <- data.frame(
  col     = colnames(cpm),
  Subject = str_extract(colnames(cpm), "PROMO\\d+"),
  Visit   = str_match(colnames(cpm), "_V(\\d)_")[, 2],
  stringsAsFactors = FALSE
)
samp_info <- samp_info[!is.na(samp_info$Subject) & samp_info$Visit %in% c("1", "2"), ]

subj_visits   <- split(samp_info$Visit, samp_info$Subject)
paired_subj   <- sort(names(subj_visits)[sapply(subj_visits,
                       function(v) all(c("1", "2") %in% v))])
message("  Paired subjects (V1 + V2 expression sample): ", length(paired_subj))
if (length(paired_subj) < 5)
  stop("Fewer than 5 paired subjects — cannot cluster V2 robustly.")

# Cross-check against metadata v2 fields (informational only)
meta_has_v2 <- pain_raw$SubjectID[
  !is.na(pain_raw$promis_global07_v2) | !is.na(pain_raw$qst_ppt_tr_avg_v2)]
expr_only <- setdiff(paired_subj, meta_has_v2)
meta_only <- setdiff(intersect(meta_has_v2, names(subj_visits)), paired_subj)
if (length(expr_only) > 0)
  message("  NOTE: V2 expression but no v2 metadata value: ",
          paste(expr_only, collapse = ", "))
if (length(meta_only) > 0)
  message("  NOTE: v2 metadata value but no paired V2 expression sample (excluded): ",
          paste(meta_only, collapse = ", "))

# ─────────────────────────────────────────────────────────────────────────────
# 3. BUILD PER-SUBJECT CPM (average technical replicates within Subject×Visit)
#    so each subject contributes exactly one V1 and one V2 column, in the same
#    subject order across both visits (e.g. PROMO24 V1 has 3 tech reps).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [3] Build paired per-subject CPM matrices ===")

genes_use <- intersect(clust_v1$ENSG, rownames(cpm))
message("  Clustering genes present in CPM matrix: ", length(genes_use),
        " / ", nrow(clust_v1))

subj_visit_mat <- function(visit) {
  mats <- sapply(paired_subj, function(s) {
    cols <- samp_info$col[samp_info$Subject == s & samp_info$Visit == visit]
    rowMeans(cpm[genes_use, cols, drop = FALSE])
  })
  colnames(mats) <- paired_subj
  mats
}
cpm_v1 <- subj_visit_mat("1")
cpm_v2 <- subj_visit_mat("2")

# ─────────────────────────────────────────────────────────────────────────────
# 4. GENE-WISE Z-SCORE PER VISIT (across the paired subjects), drop zero-sd genes
#    Mirrors the parent's z-score clustering pass: raw CPM → gene-wise center+
#    scale. Genes constant within either visit produce NaN — drop so both
#    matrices share an identical, finite gene universe.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [4] Z-score (gene-wise, per visit) ===")

z_v1_full <- t(scale(t(cpm_v1), center = TRUE, scale = TRUE))
z_v2_full <- t(scale(t(cpm_v2), center = TRUE, scale = TRUE))
genes_z   <- genes_use[is.finite(rowSums(z_v1_full)) & is.finite(rowSums(z_v2_full))]
dropped   <- setdiff(genes_use, genes_z)
if (length(dropped) > 0)
  message("  Dropped ", length(dropped), " gene(s) with sd = 0 at V1 or V2 across paired subjects")

z_v1 <- z_v1_full[genes_z, , drop = FALSE]
z_v2 <- z_v2_full[genes_z, , drop = FALSE]

v1_lab <- setNames(clust_v1$Cluster[match(genes_z, clust_v1$ENSG)], genes_z)
message("  Genes clustered in both visits: ", length(genes_z))

cluster_pal <- setNames(
  colorRampPalette(brewer.pal(max(K, 8), "Set2"))(K),
  paste0("C", seq_len(K))
)

# ─────────────────────────────────────────────────────────────────────────────
# 5. METHOD A — PROJECT genes onto FIXED V1 centroids using V2 profiles
#    V1 centroid_c = mean V1 z-profile (across subjects) of cluster c's genes.
#    Each gene's V2 z-profile (same subject order) is assigned to the nearest
#    V1 centroid. Tests whether a gene's across-subject pattern is preserved.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [5] Method A — project to fixed V1 centroids ===")

centroids_v1 <- t(sapply(seq_len(K), function(c) colMeans(z_v1[names(v1_lab)[v1_lab == c], , drop = FALSE])))
rownames(centroids_v1) <- paste0("C", seq_len(K))

v2_proj_lab <- setNames(
  apply(z_v2, 1, function(g) which.min(colSums((t(centroids_v1) - g)^2))),
  genes_z
)
proj_retention <- mean(v2_proj_lab == v1_lab)
message(sprintf("  Overall projection retention (gene keeps its V1 cluster at V2): %.1f%%",
                100 * proj_retention))

# ─────────────────────────────────────────────────────────────────────────────
# 6. METHOD B — RE-CLUSTER V2 then align labels to V1 (max overlap)
# ─────────────────────────────────────────────────────────────────────────────
message("=== [6] Method B — re-cluster V2 + align labels ===")

set.seed(KM_SEED)
km_v2       <- kmeans(z_v2, centers = K, nstart = 25, iter.max = 100)
v2_raw_lab  <- setNames(km_v2$cluster, genes_z)
map_v2      <- align_labels(ref = v1_lab, new = v2_raw_lab, K = K)
v2_recl_lab <- setNames(map_v2[v2_raw_lab], genes_z)

ari      <- adj_rand_index(v1_lab, v2_raw_lab)   # label-invariant
recl_agreement <- mean(v2_recl_lab == v1_lab)
message(sprintf("  Adjusted Rand Index (V1 vs V2 re-cluster): %.3f", ari))
message(sprintf("  Post-alignment label agreement: %.1f%%", 100 * recl_agreement))

# ─────────────────────────────────────────────────────────────────────────────
# 7. SANKEY DIAGRAMS (ggalluvial) — one per method
# ─────────────────────────────────────────────────────────────────────────────
message("=== [7] Sankey diagrams ===")

make_sankey <- function(v1, v2, method_label, fn, headline) {
  df <- as.data.frame(table(
    V1 = factor(paste0("C", v1), levels = paste0("C", seq_len(K))),
    V2 = factor(paste0("C", v2), levels = paste0("C", seq_len(K)))
  ))
  df <- df[df$Freq > 0, ]
  p <- ggplot(df, aes(axis1 = V1, axis2 = V2, y = Freq)) +
    geom_alluvium(aes(fill = V1), alpha = 0.75, width = 1/8, knot.pos = 0.3) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.2) +
    scale_x_discrete(limits = c("Visit 1", "Visit 2"), expand = c(.07, .07)) +
    scale_fill_manual(values = cluster_pal, name = "Cluster", drop = FALSE) +
    labs(
      title    = sprintf("Gene cluster mapping V1 → V2 (%s) | %s", method_label, SUF),
      subtitle = paste0(headline, " | ", length(genes_z), " genes | ",
                        length(paired_subj), " paired subjects | flows colored by V1 cluster"),
      x = NULL, y = "Genes"
    ) +
    theme_cowplot(12) + theme(legend.position = "right")
  ggsave(file.path(DIR_SANKEY, fn), p, width = 9, height = 7, dpi = 150, bg = "white")
  message("  Saved: ", fn)
}

make_sankey(v1_lab, v2_proj_lab, "project to fixed V1 centroids",
            "sankey_A_project_to_V1_centroids.png",
            sprintf("Retention %.1f%%", 100 * proj_retention))
make_sankey(v1_lab, v2_recl_lab, "re-cluster V2 + align",
            "sankey_B_recluster_aligned.png",
            sprintf("ARI %.3f", ari))

# ─────────────────────────────────────────────────────────────────────────────
# 8. SILHOUETTE PLOTS  (parent 02b style) — V1 baseline + V2 re-cluster
# ─────────────────────────────────────────────────────────────────────────────
message("=== [8] Silhouette plots ===")

silhouette_plot <- function(z_mat, labels, visit_tag, mode_label) {
  if (K < 2) { message("  k<2 — silhouette undefined"); return(invisible()) }
  sil_obj <- silhouette(as.integer(labels), dist(z_mat))
  sil_df  <- as.data.frame(sil_obj[, 1:3])
  colnames(sil_df) <- c("Cluster", "Neighbor", "SilWidth")
  sil_df$Cluster <- factor(paste0("C", sil_df$Cluster), levels = paste0("C", seq_len(K)))
  sil_df <- sil_df %>% group_by(Cluster) %>%
    arrange(desc(SilWidth), .by_group = TRUE) %>%
    mutate(gene_idx = row_number()) %>% ungroup()
  avg <- sil_df %>% group_by(Cluster) %>%
    summarise(a = round(mean(SilWidth), 3), .groups = "drop")
  lab <- paste(paste0(avg$Cluster, "=", avg$a), collapse = " | ")

  p <- ggplot(sil_df, aes(gene_idx, SilWidth, fill = Cluster)) +
    geom_col(width = 1) +
    facet_wrap(~Cluster, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = cluster_pal, guide = "none") +
    labs(title    = sprintf("Silhouette — %s | k=%d | %s", visit_tag, K, SUF),
         subtitle = paste0(mode_label, " | mean silhouette by cluster: ", lab),
         x = "Gene (sorted by silhouette width within cluster)", y = "Silhouette width") +
    theme_cowplot(11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(face = "bold"))
  fn <- sprintf("02b_silhouette_%s_%s.png", visit_tag, SUF)
  ggsave(file.path(DIR_SIL, fn), p, width = max(8, K * 2.5), height = 4, dpi = 150, bg = "white")
  message("  Saved: ", fn)
}
silhouette_plot(z_v1, v1_lab,      "V1_paired",       "V1 labels (from input CSV), distance on V1 z-score")
silhouette_plot(z_v2, v2_recl_lab, "V2_recluster",    "V2 re-cluster (aligned), distance on V2 z-score")

# ─────────────────────────────────────────────────────────────────────────────
# 9. GENE-PCA SCATTER  (parent 02f style) — each point is one gene
# ─────────────────────────────────────────────────────────────────────────────
message("=== [9] Gene-PCA scatter ===")

gene_pca_plot <- function(z_mat, labels, visit_tag, mode_label) {
  pca <- prcomp(z_mat, center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    Cluster = factor(paste0("C", labels[rownames(pca$x)]),
                                     levels = paste0("C", seq_len(K))))
  p <- ggplot(df, aes(PC1, PC2, color = Cluster)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = cluster_pal, name = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title    = sprintf("Gene PCA colored by cluster — %s | k=%d | %s", visit_tag, K, SUF),
         subtitle = paste0(nrow(df), " genes | ", mode_label,
                           "\nSmeared colors over one cloud = weak cluster structure"),
         x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_cowplot(12)
  fn <- sprintf("02f_gene_PCA_%s_%s.png", visit_tag, SUF)
  ggsave(file.path(DIR_PCA, fn), p, width = 8, height = 6, dpi = 150, bg = "white")
  message("  Saved: ", fn)
}
gene_pca_plot(z_v1, v1_lab,      "V1_paired",    "PCA of V1 z-score matrix | V1 cluster labels")
gene_pca_plot(z_v2, v2_recl_lab, "V2_recluster", "PCA of V2 z-score matrix | V2 re-cluster (aligned)")

# ─────────────────────────────────────────────────────────────────────────────
# 9B. CLUSTERED HEATMAPS  (parent 02_heatmap_kmeans style)
#   Display transform mirrors the parent's z-score pass exactly:
#     gene-wise z-score  →  per-gene min-max [-1, 1]  →  RdBu, breaks -1..1.
#   Rows NOT clustered — ordered by cluster (order(labels)); columns clustered by
#   correlation distance (1 - Pearson, complete linkage), matching the parent.
#   annotation_row = GeneCluster (shared cluster_pal); annotation_col = per-visit
#   subject metadata (V1 heatmap → v1 fields, V2 heatmaps → v2 fields).
# ─────────────────────────────────────────────────────────────────────────────
message("=== [9B] Clustered heatmaps (parent style) ===")

# Per-visit subject-annotation column specs (display name -> metadata column).
# Names kept identical across visits so one ann_colors list serves both.
hm_ann_spec_v1 <- c(PPT = "qst_ppt_tr_avg_v1", PROMIS = "promis_global07_v1",
                    Global07_ge4 = "global07_ge4_v1", Age = "age")
hm_ann_spec_v2 <- c(PPT = "qst_ppt_tr_avg_v2", PROMIS = "promis_global07_v2",
                    Global07_ge4 = "global07_ge4_V2", Age = "age")

# Parent color scheme: continuous diverging for PPT/PROMIS, sequential for Age,
# discrete for Sex and the global07 ≥4 binary flag.
ann_colors_sub <- list(
  PPT          = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  PROMIS       = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  Age          = colorRampPalette(c("#fff7ec", "#7f2704"))(100),
  Sex          = c(Female = "#e78ac3", Male = "#66c2a5"),
  Global07_ge4 = c(`0` = "#bdbdbd", `1` = "#de2d26")
)

build_ann_col <- function(spec, subj_order) {
  df <- data.frame(row.names = subj_order)
  for (nm in names(spec)) {
    col <- spec[[nm]]
    if (!col %in% colnames(pain_raw)) {
      message("  Heatmap annotation col not found, skipped: ", col); next
    }
    if (nm == "Global07_ge4") {
      v <- suppressWarnings(as.numeric(pain_raw[subj_order, col]))
      df[[nm]] <- factor(ifelse(is.na(v), NA_character_, as.character(as.integer(v))),
                         levels = c("0", "1"))
    } else {
      df[[nm]] <- as.numeric(pain_raw[subj_order, col])
    }
  }
  sx <- pain_raw[subj_order, "sex_at_birth.factor"]
  df$Sex <- factor(ifelse(sx %in% c("Male", "Female"), as.character(sx), NA_character_),
                   levels = c("Female", "Male"))
  df
}

make_cluster_heatmap <- function(z_mat, labels, ann_col, visit_tag, mode_label) {
  # z-score (already gene-wise) → per-gene min-max [-1, 1] for display
  disp        <- t(apply(z_mat, 1, min_max_norm))
  row_order   <- order(labels[rownames(disp)])
  mat_display <- disp[row_order, , drop = FALSE]
  row_ann <- data.frame(
    GeneCluster = factor(paste0("C", labels[rownames(mat_display)]),
                         levels = paste0("C", seq_len(K))),
    row.names   = rownames(mat_display)
  )
  # Column dendrogram: correlation distance on the z-score matrix (parent uses
  # 1 - Pearson, complete linkage). Subjects already carry short IDs (PROMO##).
  hc_cols <- hclust(as.dist(1 - cor(z_mat, method = "pearson")), method = "complete")

  hm_title <- sprintf("%d genes | %s | k=%d | %s | %d paired subj | %s",
                      nrow(mat_display), visit_tag, K, mode_label,
                      length(paired_subj), SUF)
  fn <- sprintf("02_heatmap_kmeans_%s_%s.png", visit_tag, SUF)
  png(file.path(DIR_HEATMAP, fn), width = 2400, height = 2400, res = 150, type = "cairo")
  pheatmap(mat_display,
           cluster_rows      = FALSE, cluster_cols = hc_cols,
           annotation_col    = ann_col,
           annotation_row    = row_ann,
           annotation_colors = c(ann_colors_sub, list(GeneCluster = cluster_pal)),
           color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks            = seq(-1, 1, length.out = 101),
           show_rownames     = FALSE, show_colnames = TRUE, fontsize_col = 9,
           main              = hm_title,
           fontsize = 8, border_color = NA, annotation_legend = TRUE)
  dev.off()
  message("  Saved: ", fn)
}

ann_v1_hm <- build_ann_col(hm_ann_spec_v1, paired_subj)
ann_v2_hm <- build_ann_col(hm_ann_spec_v2, paired_subj)
make_cluster_heatmap(z_v1, v1_lab,      ann_v1_hm, "V1_paired",     "V1 labels (input CSV)")
make_cluster_heatmap(z_v2, v2_recl_lab, ann_v2_hm, "V2_recluster",  "V2 re-cluster (aligned)")
make_cluster_heatmap(z_v2, v2_proj_lab, ann_v2_hm, "V2_project",    "V2 project to fixed V1 centroids")

# ─────────────────────────────────────────────────────────────────────────────
# 10. Z-SCORE-PER-PATIENT CSV + CLUSTER-MEAN LINE PLOT + CORRELATION MATRIX
#     (parent 02i + correlation block). V1 uses v1 metadata; V2 uses v2 metadata.
# ─────────────────────────────────────────────────────────────────────────────
message("=== [10] Cluster-mean per-patient + metadata correlation ===")

# Metadata column specs per visit (display name = metadata column)
meta_spec_v1 <- c(qst_ppt_v1 = "qst_ppt_tr_avg_v1", promis_v1 = "promis_global07_v1",
                  global07_ge4_v1 = "global07_ge4_v1", age = "age")
meta_spec_v2 <- c(qst_ppt_v2 = "qst_ppt_tr_avg_v2", promis_v2 = "promis_global07_v2",
                  global07_ge4_V2 = "global07_ge4_V2", age = "age")

build_meta_mat <- function(spec, subj_order) {
  present <- spec[spec %in% colnames(pain_raw)]
  miss    <- setdiff(spec, colnames(pain_raw))
  if (length(miss) > 0)
    message("  Metadata cols not found, skipped: ", paste(miss, collapse = ", "))
  mm <- as.data.frame(lapply(present, function(cn) as.numeric(pain_raw[subj_order, cn])))
  names(mm) <- names(present)
  sx <- pain_raw[subj_order, "sex_at_birth.factor"]
  mm$sex_male <- ifelse(is.na(sx), NA_integer_, as.integer(sx == "Male"))
  rownames(mm) <- subj_order
  mm
}

cluster_meta_block <- function(z_mat, labels, meta_mat, visit_tag, mode_label) {
  subj_order <- colnames(z_mat)
  # per-cluster mean z-score per subject
  cmat <- do.call(rbind, lapply(seq_len(K), function(ci) {
    g <- names(labels)[labels == ci]
    if (length(g) == 0) return(rep(NA_real_, length(subj_order)))
    colMeans(z_mat[g, subj_order, drop = FALSE])
  }))
  rownames(cmat) <- paste0("C", seq_len(K)); colnames(cmat) <- subj_order

  # gene × subject z-score CSV with cluster labels
  zcsv <- data.frame(
    ENSG = rownames(z_mat),
    Gene_Symbol = gene_sym_map$Gene_Symbol[match(rownames(z_mat), gene_sym_map$ENSG)],
    Cluster = labels[rownames(z_mat)], stringsAsFactors = FALSE
  )
  zcsv <- cbind(zcsv, as.data.frame(z_mat))
  zcsv <- zcsv[order(zcsv$Cluster), ]
  write.csv(zcsv, file.path(DIR_CLUSTMETA,
            sprintf("zscore_per_patient_%s_%s.csv", visit_tag, SUF)), row.names = FALSE)

  # faceted cluster-mean line plot (subject order = hierarchical-cluster order)
  hc_order <- subj_order[hclust(dist(t(z_mat)), method = "complete")$order]
  line_df  <- do.call(rbind, lapply(seq_len(K), function(ci) {
    data.frame(Cluster = paste0("C", ci), Subject = subj_order,
               MeanZScore = as.numeric(cmat[ci, ]), stringsAsFactors = FALSE)
  }))
  line_df$Cluster <- factor(line_df$Cluster, levels = paste0("C", seq_len(K)))
  line_df$Subject <- factor(line_df$Subject, levels = hc_order)
  p_line <- ggplot(line_df, aes(Subject, MeanZScore, group = Cluster, color = Cluster)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    facet_wrap(~Cluster, ncol = 2, scales = "free_y") +
    scale_color_manual(values = cluster_pal, guide = "none") +
    labs(title    = sprintf("Per-cluster mean z-score per subject — %s | k=%d | %s", visit_tag, K, SUF),
         subtitle = paste0(mode_label, " | y = mean z-score across cluster genes"),
         x = NULL, y = "Mean z-score") +
    theme_cowplot(11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(DIR_CLUSTMETA, sprintf("02i_cluster_mean_zscore_per_subject_%s_%s.png", visit_tag, SUF)),
         p_line, width = max(12, length(subj_order) * 0.30),
         height = 2 + 1.8 * ceiling(K / 2), dpi = 150, bg = "white")

  # per-cluster mean-z vs metadata correlation tile
  cor_mat <- matrix(NA_real_, K, ncol(meta_mat),
                    dimnames = list(rownames(cmat), colnames(meta_mat)))
  p_mat <- cor_mat
  for (ci in seq_len(K)) for (mc in seq_len(ncol(meta_mat))) {
    x <- cmat[ci, ]; y <- meta_mat[subj_order, mc]; ok <- !is.na(x) & !is.na(y)
    if (sum(ok) >= 3 && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok]))
      cor_mat[ci, mc] <- ct$estimate; p_mat[ci, mc] <- ct$p.value
    }
  }
  cl <- data.frame(Cluster = rownames(cor_mat)[row(cor_mat)],
                   Metadata = colnames(cor_mat)[col(cor_mat)],
                   r = as.vector(cor_mat), p = as.vector(p_mat), stringsAsFactors = FALSE)
  cl$lab <- ifelse(is.na(cl$r), "", sprintf("%.2f%s", cl$r,
                   ifelse(!is.na(cl$p) & cl$p < 0.05, "*", "")))
  cl$Cluster  <- factor(cl$Cluster,  levels = rownames(cor_mat))
  cl$Metadata <- factor(cl$Metadata, levels = colnames(cor_mat))
  p_corr <- ggplot(cl, aes(Metadata, Cluster, fill = r)) +
    geom_tile(color = "white") + geom_text(aes(label = lab), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey85") +
    labs(title    = sprintf("Per-cluster mean z-score vs metadata — %s | k=%d | %s", visit_tag, K, SUF),
         subtitle = paste0("Pearson r | * = p<0.05 | n = ", length(subj_order), " paired subjects"),
         x = NULL, y = NULL) +
    theme_cowplot(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(DIR_CLUSTMETA, sprintf("cluster_metadata_correlation_%s_%s.png", visit_tag, SUF)),
         p_corr, width = 8, height = 5, dpi = 150, bg = "white")
  write.csv(cl, file.path(DIR_CLUSTMETA,
            sprintf("cluster_metadata_correlation_%s_%s.csv", visit_tag, SUF)), row.names = FALSE)
  message("  ", visit_tag, ": z-score CSV + line plot + correlation written")
}

cluster_meta_block(z_v1, v1_lab,      build_meta_mat(meta_spec_v1, paired_subj),
                   "V1_paired",    "V1 labels (input CSV)")
cluster_meta_block(z_v2, v2_recl_lab, build_meta_mat(meta_spec_v2, paired_subj),
                   "V2_recluster", "V2 re-cluster (aligned)")

# ─────────────────────────────────────────────────────────────────────────────
# 11. PER-CLUSTER ORA  (parent 02c style; clusterProfiler) — V1 + V2
# ─────────────────────────────────────────────────────────────────────────────
message("=== [11] Per-cluster ORA ===")

if (HAS_CP) {
  hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
  pathways_h    <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
  gobp_sets     <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  pathways_gobp <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)

  ora_bg_syms <- unique(na.omit(
    gene_sym_map$Gene_Symbol[gene_sym_map$ENSG %in% rownames(cpm)]))

  run_ora_cp <- function(fg_syms, bg_syms, pathway_list, min_gs = 5, max_gs = 500) {
    t2g <- do.call(rbind, lapply(names(pathway_list), function(pw)
      data.frame(term = pw, gene = pathway_list[[pw]], stringsAsFactors = FALSE)))
    res <- tryCatch(clusterProfiler::enricher(
      gene = fg_syms, universe = bg_syms, TERM2GENE = t2g,
      minGSSize = min_gs, maxGSSize = max_gs, pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1), error = function(e) NULL)
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    df <- as.data.frame(res)
    data.frame(Description = df$ID, Count = df$Count,
               GeneRatio = sapply(df$GeneRatio, function(x) eval(parse(text = x))),
               pvalue = df$pvalue, p.adjust = df$p.adjust, stringsAsFactors = FALSE)
  }
  make_ora_panel <- function(df, cluster_label, top_n = 15, strip_prefix = NULL) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    pdf <- df %>% arrange(p.adjust) %>% slice_head(n = top_n) %>%
      mutate(neglog10padj = pmin(-log10(p.adjust), 10),
             pw_label = if (!is.null(strip_prefix)) gsub(strip_prefix, "", Description, perl = TRUE) else Description,
             pw_label = gsub("_", " ", pw_label),
             pw_label = factor(pw_label, levels = rev(pw_label)))
    ggplot(pdf, aes(GeneRatio, pw_label, size = Count, color = neglog10padj)) +
      geom_point() +
      scale_color_gradient(low = "#AAAAAA", high = "#B2182B", limits = c(0, 10),
                           oob = scales::squish, name = "-log10(p.adj)") +
      scale_size_continuous(range = c(2, 7), name = "Genes") +
      labs(title = cluster_label, x = "Gene Ratio", y = NULL) +
      theme_cowplot(9) +
      theme(plot.title = element_text(face = "bold", size = 10),
            axis.text.y = element_text(size = 7), legend.position = "bottom")
  }

  ora_for_visit <- function(labels, visit_tag) {
    for (spec in list(list(pw = pathways_h,    strip = "^HALLMARK_", tag = "Hallmark"),
                      list(pw = pathways_gobp, strip = "^GOBP_",     tag = "GOBP"))) {
      panels <- lapply(seq_len(K), function(ci) {
        fg <- names(labels)[labels == ci]
        fg_syms <- unique(na.omit(gene_sym_map$Gene_Symbol[match(fg, gene_sym_map$ENSG)]))
        res <- run_ora_cp(fg_syms, ora_bg_syms, spec$pw)
        if (!is.null(res) && nrow(res) > 0)
          write.csv(res, file.path(DIR_ORA,
            sprintf("ORA_%s_C%d_%s_%s.csv", visit_tag, ci, SUF, spec$tag)), row.names = FALSE)
        make_ora_panel(res, paste0("C", ci), strip_prefix = spec$strip)
      })
      panels <- panels[!sapply(panels, is.null)]
      if (length(panels) > 0) {
        nc <- min(3, length(panels))
        fn <- sprintf("02c_ORA_facet_%s_%s_%s.png", visit_tag, SUF, spec$tag)
        ggsave(file.path(DIR_ORA, fn), plot_grid(plotlist = panels, ncol = nc),
               width = 10 * nc, height = 7 * ceiling(length(panels) / nc),
               dpi = 150, bg = "white")
        message("  Saved: ", fn)
      }
    }
  }
  ora_for_visit(v1_lab,      "V1_paired")
  ora_for_visit(v2_recl_lab, "V2_recluster")
} else {
  message("  clusterProfiler unavailable — ORA skipped")
}

# ─────────────────────────────────────────────────────────────────────────────
# 12. ROBUSTNESS SUMMARY — per-gene assignments + contingency tables + metrics
# ─────────────────────────────────────────────────────────────────────────────
message("=== [12] Robustness summary ===")

assign_df <- data.frame(
  ENSG               = genes_z,
  Gene_Symbol        = gene_sym_map$Gene_Symbol[match(genes_z, gene_sym_map$ENSG)],
  V1_cluster         = v1_lab[genes_z],
  V2_project         = v2_proj_lab[genes_z],
  V2_recluster_aligned = v2_recl_lab[genes_z],
  stays_project      = v1_lab[genes_z] == v2_proj_lab[genes_z],
  stays_recluster    = v1_lab[genes_z] == v2_recl_lab[genes_z],
  stringsAsFactors = FALSE
)
write.csv(assign_df, file.path(OUT_DIR, sprintf("gene_cluster_assignments_V1V2_%s.csv", SUF)),
          row.names = FALSE)

write.csv(as.data.frame.matrix(table(V1 = v1_lab, V2_project = v2_proj_lab)),
          file.path(OUT_DIR, "contingency_A_project.csv"))
write.csv(as.data.frame.matrix(table(V1 = v1_lab, V2_recluster = v2_recl_lab)),
          file.path(OUT_DIR, "contingency_B_recluster.csv"))

per_clust <- data.frame(
  Cluster = paste0("C", seq_len(K)),
  n_genes_V1 = as.integer(table(factor(v1_lab, levels = seq_len(K)))),
  retention_project = sapply(seq_len(K), function(c)
    mean(v2_proj_lab[v1_lab == c] == c)),
  retention_recluster = sapply(seq_len(K), function(c)
    mean(v2_recl_lab[v1_lab == c] == c))
)
write.csv(per_clust, file.path(OUT_DIR, sprintf("per_cluster_retention_%s.csv", SUF)), row.names = FALSE)

summary_lines <- c(
  "PROMOTER V1→V2 Cluster Robustness Summary",
  paste("Run date:", Sys.time()),
  paste("Paired subjects:", length(paired_subj)),
  paste("Genes clustered in both visits:", length(genes_z)),
  paste("k:", K),
  "",
  sprintf("Method A (project to fixed V1 centroids) — overall retention: %.1f%%", 100 * proj_retention),
  sprintf("Method B (re-cluster V2 + align)         — Adjusted Rand Index: %.3f", ari),
  sprintf("Method B post-alignment label agreement: %.1f%%", 100 * recl_agreement),
  "",
  "Per-cluster retention:",
  paste(capture.output(print(per_clust, row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path(OUT_DIR, "robustness_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
message("=== DONE — outputs in ", OUT_DIR, " ===")
