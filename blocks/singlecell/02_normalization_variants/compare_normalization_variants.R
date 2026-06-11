#!/usr/bin/env Rscript
# =============================================================================
# compare_normalization_variants.R
#
# Author: Tyler Therron
# Purpose:
#   Mirror singlecell_figure_engine.R but run the IMPACT cohort
#   through four normalization / integration variants on the same dims & res:
#     - merged   : NormalizeData + FindVariableFeatures + ScaleData + PCA
#     - sctv2    : SCTransform v2 + PCA
#     - rpca     : per-pool SCTransform v2 + Seurat RPCA integration
#     - harmony  : SCTransform v2 + Harmony on hash_run
#
#   IMP12 is excluded up-front by default. CellType-dependent figures from the
#   v5-2 script are skipped because each variant produces its own clusters.
#   UMAP / QC / ADT / clustree / de novo cluster marker figures are produced per
#   variant per (dim_spec x resolution) combo, then a cross-variant comparison
#   panel at each combo.
#
# pca_dims now accepts a list of dim specifications separated by ';'. Within a
# spec, comma-separated ranges/integers are unioned. Bare integer N is shorthand
# for 1:N (backwards-compatible). Examples:
#   --pca_dims 12                       -> [ 1:12 ]
#   --pca_dims 1:9,11:12                -> [ {1..9, 11, 12} ]    (skip PC10)
#   --pca_dims "1:9,11:12;1:10;1:12"    -> three specs
#
# Variants can be dropped two ways (they compose):
#   --variant merged
#   --variant all --exclude_variants sctv2,rpca,harmony
#
# Each (dim_spec x resolution) combo also gets a PC-overlay facet UMAP: for
# each PC used in that spec, the cell's PC (or Harmony) score is rendered as
# a feature plot on the UMAP, paged into groups of four.
#
# Example:
#   Rscript compare_normalization_variants.R \
#     --input_rds /path/to/WHL1-5_raw_merged.rds \
#     --out_dir   /path/to/output \
#     --exclude_variants sctv2,rpca,harmony \
#     --pca_dims "1:9,11:12;1:10;1:11;1:12" \
#     --final_res 0.3 \
#     --resolutions 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0
# =============================================================================

.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(sctransform)
  library(glmGamPoi)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(RColorBrewer)
  library(colorspace)
  library(scales)
  library(clustree)
  library(ggrepel)
  library(Matrix)
  library(optparse)
  library(glue)
})

# =============================================================================
## CLI
# =============================================================================
option_list <- list(
  make_option("--input_rds",        type = "character", default = NULL,
              help = "Raw pre-norm merged WHL1-5 RDS [REQUIRED]"),
  make_option("--out_dir",          type = "character", default = NULL,
              help = "Output base directory [REQUIRED]"),
  make_option("--variant",          type = "character", default = "all",
              help = "Comma-list of variants or 'all'. Options: merged, sctv2, rpca, harmony [default: all]"),
  make_option("--exclude_variants", type = "character", default = "",
              help = "Comma-list of variants to drop after --variant resolution. Composes with --variant. [default: none]"),
  make_option("--exclude_samples",  type = "character", default = "IMP12",
              help = "Comma-list of HTO_maxID values to drop before any processing [default: IMP12]"),
  make_option("--sample_col",       type = "character", default = "HTO_maxID",
              help = "Metadata column for patient IMP IDs [default: HTO_maxID]"),
  make_option("--pool_col",         type = "character", default = "hash_run",
              help = "Metadata column for WHL pool / batch [default: hash_run]"),
  make_option("--n_hvg",            type = "integer",   default = 3000L,
              help = "Variable features for SCT-based variants [default: 3000]"),
  make_option("--n_hvg_logn",       type = "integer",   default = 2000L,
              help = "Variable features for the merged log-norm variant [default: 2000]"),
  make_option("--n_pcs",            type = "integer",   default = 30L,
              help = "Number of PCs to compute [default: 30]. Auto-bumped if a dim spec asks for more."),
  make_option("--pca_dims",         type = "character", default = "20",
              help = paste0(
                "Dim specs separated by ';'. Within a spec, comma-list of ",
                "ranges/integers (e.g. '1:9,11:12'). Bare int N means 1:N. ",
                "Example: '1:9,11:12;1:10;1:12'. [default: 20]"
              )),
  make_option("--final_res",        type = "double",    default = 0.5,
              help = "Clustering resolution to mark as 'final' in clustree/comparison [default: 0.5]"),
  make_option("--resolutions",      type = "character",
              default = "0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0",
              help = "Comma-list of resolutions for the clustering sweep / clustree"),
  make_option("--rpca_k_anchor",    type = "integer",   default = 5L,
              help = "k.anchor for FindIntegrationAnchors (RPCA) [default: 5]"),
  make_option("--skip_stages",      type = "character", default = "",
              help = paste0(
                "Comma-list of figure stages to skip. Useful for exploratory ",
                "dim-spec sweeps where you don't want to wait for slow stages. ",
                "Options: elbow, clustree, umaps, qc, adt, composition, ",
                "markers, pc_overlay, compare. [default: none]"
              )),
  make_option("--seed",             type = "integer",   default = 42L)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds) || is.null(opt$out_dir))
  stop("--input_rds and --out_dir are required.")

ALL_VARIANTS <- c("merged", "sctv2", "rpca", "harmony")
requested_variants <- if (grepl("^all$", trimws(opt$variant), ignore.case = TRUE))
  ALL_VARIANTS else trimws(strsplit(opt$variant, ",")[[1]])
excluded_variants <- if (nzchar(trimws(opt$exclude_variants)))
  trimws(strsplit(opt$exclude_variants, ",")[[1]]) else character(0)
unknown_req  <- setdiff(requested_variants, ALL_VARIANTS)
unknown_excl <- setdiff(excluded_variants, ALL_VARIANTS)
if (length(unknown_req) > 0)
  stop("Unknown variant(s) in --variant: ", paste(unknown_req, collapse = ", "),
       ". Allowed: ", paste(ALL_VARIANTS, collapse = ", "))
if (length(unknown_excl) > 0)
  stop("Unknown variant(s) in --exclude_variants: ",
       paste(unknown_excl, collapse = ", "),
       ". Allowed: ", paste(ALL_VARIANTS, collapse = ", "))
active_variants <- setdiff(requested_variants, excluded_variants)
if (length(active_variants) == 0)
  stop("No active variants after --variant / --exclude_variants resolution.")

ALL_STAGES <- c("elbow", "clustree", "umaps", "qc", "adt",
                "composition", "markers", "pc_overlay", "compare")
skip_stages <- if (nzchar(trimws(opt$skip_stages)))
  trimws(strsplit(opt$skip_stages, ",")[[1]]) else character(0)
unknown_stages <- setdiff(skip_stages, ALL_STAGES)
if (length(unknown_stages) > 0)
  stop("Unknown --skip_stages value(s): ",
       paste(unknown_stages, collapse = ", "),
       ". Allowed: ", paste(ALL_STAGES, collapse = ", "))

# -----------------------------------------------------------------------------
# pca_dims parsing
# -----------------------------------------------------------------------------
# parse_dim_spec("12")          -> 1:12        (bare-integer shorthand)
# parse_dim_spec("1:9,11:12")   -> c(1:9, 11:12)
# parse_dim_spec("5,10")        -> c(5, 10)    (literal, NOT shorthand)
parse_dim_spec <- function(spec) {
  spec <- trimws(spec)
  if (!nzchar(spec)) stop("Empty pca_dims spec.")
  if (!grepl("[:,]", spec)) {
    n <- suppressWarnings(as.integer(spec))
    if (is.na(n) || n < 1)
      stop("Invalid --pca_dims bare integer: '", spec, "' (must be >= 1).")
    return(seq_len(n))
  }
  pieces <- trimws(strsplit(spec, ",")[[1]])
  pieces <- pieces[nzchar(pieces)]
  dims <- unlist(lapply(pieces, function(p) {
    if (grepl(":", p)) {
      rng <- suppressWarnings(as.integer(trimws(strsplit(p, ":")[[1]])))
      if (length(rng) != 2 || any(is.na(rng)) || rng[1] < 1 || rng[2] < rng[1])
        stop("Invalid --pca_dims range: '", p, "'")
      return(seq.int(rng[1], rng[2]))
    }
    v <- suppressWarnings(as.integer(p))
    if (is.na(v) || v < 1)
      stop("Invalid --pca_dims value: '", p, "'")
    v
  }))
  sort(unique(as.integer(dims)))
}

parse_dim_specs <- function(s) {
  s <- trimws(s)
  if (!nzchar(s)) stop("--pca_dims is empty.")
  specs <- trimws(strsplit(s, ";")[[1]])
  specs <- specs[nzchar(specs)]
  if (length(specs) == 0) stop("--pca_dims parsed to zero specs.")
  out <- lapply(specs, parse_dim_spec)
  # Drop duplicates (same dim set after sort) while preserving first-seen order
  keys <- vapply(out, function(d) paste(d, collapse = ","), character(1))
  out[!duplicated(keys)]
}

# format_dim_tag(c(1:9, 11:12)) -> "1-9_11-12"
# format_dim_tag(1:10)          -> "1-10"
# format_dim_tag(c(1,3,5))      -> "1_3_5"
format_dim_tag <- function(dims) {
  dims <- sort(unique(as.integer(dims)))
  if (length(dims) == 0) return("")
  if (length(dims) == 1) return(as.character(dims))
  d <- diff(dims)
  starts <- c(1L, which(d != 1L) + 1L)
  ends   <- c(starts[-1] - 1L, length(dims))
  parts  <- mapply(function(s, e) {
    if (s == e) as.character(dims[s])
    else paste0(dims[s], "-", dims[e])
  }, starts, ends, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  paste(parts, collapse = "_")
}

pca_dims_specs <- parse_dim_specs(opt$pca_dims)
max_dim_req    <- max(unlist(pca_dims_specs))
if (max_dim_req > opt$n_pcs) {
  message("  NOTE: max(pca_dims_specs)=", max_dim_req,
          " > --n_pcs=", opt$n_pcs, " -- bumping n_pcs to ", max_dim_req)
  opt$n_pcs <- max_dim_req
}

# =============================================================================
## CONFIG
# =============================================================================
cfg <- list(
  input_rds       = opt$input_rds,
  out_dir         = opt$out_dir,
  exclude_samples = trimws(strsplit(opt$exclude_samples, ",")[[1]]),
  sample_col      = opt$sample_col,
  pool_col        = opt$pool_col,

  n_hvg           = opt$n_hvg,
  n_hvg_logn      = opt$n_hvg_logn,
  n_pcs           = opt$n_pcs,
  pca_dims_specs  = pca_dims_specs,
  final_res       = opt$final_res,
  resolutions     = as.numeric(trimws(strsplit(opt$resolutions, ",")[[1]])),
  rpca_k_anchor   = opt$rpca_k_anchor,
  skip_stages     = skip_stages,

  regress_vars     = c("percent.mt", "S.Score", "G2M.Score"),
  exclude_patterns = c("^MT-", "^RPS", "^RPL",
                       "^IGHV", "^IGLV", "^IGKV",
                       "^TRAV", "^TRBV", "^TRDV", "^TRGV"),

  rna_assay    = "RNA",
  adt_assay    = "ADT",
  count_layer  = "counts",

  low_pain  = c("IMP38","IMP35","IMP29","IMP11","IMP13","IMP21","IMP24",
                "IMP37","IMP42","IMP30","IMP34"),
  high_pain = c("IMP16","IMP18","IMP19","IMP33","IMP36","IMP39","IMP46",
                "IMP25","IMP23","IMP15"),

  pool_groups = list(
    "WHL1"   = c("IMP38","IMP35","IMP29","IMP25","IMP23","IMP15"),
    "WHL2/3" = c("IMP11","IMP19","IMP21","IMP30","IMP33","IMP37","IMP39"),
    "WHL4/5" = c("IMP13","IMP16","IMP18","IMP24","IMP34","IMP36","IMP42","IMP46")
  ),
  pool_cols = c("WHL1" = "#E41A1C", "WHL2/3" = "#377EB8", "WHL4/5" = "#4DAF4A"),

  adt_markers = c(
    "Hu.CD3-UCHT1"  = "CD3",
    "Hu.CD4-RPA.T4" = "CD4",
    "Hu.CD8"        = "CD8",
    "Hu.CD335"      = "NKp46",
    "Hu.CD20-2H7"   = "CD20",
    "Hu.CD14-M5E2"  = "CD14",
    "Hu.CD16"       = "CD16",
    "Hu.CD1c"       = "CD1c",
    "Hu.CD303"      = "CD303"
  ),

  umap_w = 10,
  umap_h = 8,

  # Point size for the PC-overlay facet plot (user requested 'large').
  pc_overlay_pt_size = 1.8
)

set.seed(opt$seed)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

message("== IMPACT Normalization Variant Pipeline ==")
message("  input_rds        : ", cfg$input_rds)
message("  out_dir          : ", cfg$out_dir)
message("  variants         : ", paste(active_variants, collapse = ", "),
        if (length(excluded_variants))
          paste0("   (excluded: ", paste(excluded_variants, collapse = ", "), ")")
        else "")
message("  excluded samples : ", paste(cfg$exclude_samples, collapse = ", "))
message("  pca_dims specs   : ",
        paste(vapply(cfg$pca_dims_specs, format_dim_tag, character(1)),
              collapse = " | "))
message("  n_pcs (compute)  : ", cfg$n_pcs)
message("  final_res        : ", cfg$final_res)
message("  resolutions      : ", paste(cfg$resolutions, collapse = ", "))
message("  n_hvg (SCT)      : ", cfg$n_hvg, "  n_hvg (log-norm): ", cfg$n_hvg_logn)
if (length(cfg$skip_stages) > 0)
  message("  skip_stages      : ", paste(cfg$skip_stages, collapse = ", "))

# =============================================================================
## HELPERS
# =============================================================================

save_png <- function(p, name, dir, width, height, dpi = 300) {
  path <- file.path(dir, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, units = "in", dpi = dpi,
         bg = "white", limitsize = FALSE)
  message("  Saved: ", basename(path))
  invisible(path)
}

# DoHeatmap + ggsave frequently produces a blank PNG; render via a real device.
save_png_device <- function(p, name, dir, width, height, dpi = 300) {
  path <- file.path(dir, paste0(name, ".png"))
  png(path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  print(p)
  dev.off()
  message("  Saved: ", basename(path))
  invisible(path)
}

should_skip <- function(stage, cfg) isTRUE(stage %in% cfg$skip_stages)

# Pulls the per-combo metadata stamped onto obj before each figure stage.
umap_subtitle <- function(obj, show_res = TRUE) {
  cm <- attr(obj, "variant_meta")
  if (is.null(cm)) return("")
  if (show_res)
    paste0("res = ", cm$res, "  |  dims = ", cm$dim_tag)
  else
    paste0("dims = ", cm$dim_tag)
}

add_clustree_res_marker <- function(p, resolutions, final_res) {
  res_sorted <- sort(unique(resolutions))
  idx <- match(final_res, res_sorted)
  if (is.na(idx)) return(p)
  p +
    geom_hline(yintercept = idx - 1, linetype = "dashed",
               colour = "red", linewidth = 0.9) +
    labs(caption = paste0("red dashed line = final resolution (", final_res, ")"))
}

filter_hvg <- function(genes, patterns) {
  keep <- rep(TRUE, length(genes))
  for (pat in patterns) keep <- keep & !grepl(pat, genes)
  genes[keep]
}

safe_counts <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
}

ensure_mito_ribo <- function(obj) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt <- grep("^MT-", rownames(obj), value = TRUE)
    obj[["percent.mt"]] <- if (length(mt) > 0)
      PercentageFeatureSet(obj, pattern = "^MT-") else 0
  }
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    rb <- grep("^RP[SL][[:digit:]]|^RPSA", rownames(obj), value = TRUE)
    obj[["percent.ribo"]] <- if (length(rb) > 0)
      PercentageFeatureSet(obj, features = rb) else 0
  }
  obj
}

make_sample_pal <- function(samples) {
  samples <- sort(unique(samples))
  base_cols <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"),
                 brewer.pal(8, "Dark2"))
  setNames(colorRampPalette(base_cols)(length(samples)), samples)
}

make_pool_pal <- function(pools, cfg) {
  pools <- sort(unique(as.character(pools)))
  hit   <- pools %in% names(cfg$pool_cols)
  pal   <- character(length(pools))
  pal[hit]  <- cfg$pool_cols[pools[hit]]
  if (any(!hit)) {
    pal[!hit] <- colorRampPalette(brewer.pal(8, "Dark2"))(sum(!hit))
  }
  setNames(pal, pools)
}

assign_pool_group <- function(sample_ids, cfg) {
  pool <- rep(NA_character_, length(sample_ids))
  for (grp in names(cfg$pool_groups))
    pool[sample_ids %in% cfg$pool_groups[[grp]]] <- grp
  pool
}

assign_pain_group <- function(obj, cfg) {
  ids <- as.character(obj@meta.data[[cfg$sample_col]])
  pg <- dplyr::case_when(
    ids %in% cfg$low_pain  ~ "Low Pain",
    ids %in% cfg$high_pain ~ "High Pain",
    TRUE ~ NA_character_
  )
  obj$PainGroup <- factor(pg, levels = c("Low Pain", "High Pain"))
  message(glue("  PainGroup assigned: ",
               "{sum(pg == 'Low Pain',  na.rm = TRUE)} Low, ",
               "{sum(pg == 'High Pain', na.rm = TRUE)} High, ",
               "{sum(is.na(pg))} NA"))
  obj
}

# =============================================================================
## LOAD + PREP (run once, shared across variants)
# =============================================================================
load_and_prep <- function(cfg) {
  message("\n[load] Reading: ", cfg$input_rds)
  obj <- readRDS(cfg$input_rds)
  DefaultAssay(obj) <- cfg$rna_assay
  message("  Cells: ", ncol(obj), "  Genes: ", nrow(obj))

  if (length(cfg$exclude_samples) > 0) {
    keep_cells <- rownames(obj@meta.data)[
      !as.character(obj@meta.data[[cfg$sample_col]]) %in% cfg$exclude_samples
    ]
    n_before <- ncol(obj)
    obj <- subset(obj, cells = keep_cells)
    message(glue("  Excluded {paste(cfg$exclude_samples, collapse=',')}: ",
                 "{n_before} -> {ncol(obj)} cells"))
  }

  tryCatch(obj <- JoinLayers(obj, assay = cfg$rna_assay),
           error = function(e) NULL)

  obj <- ensure_mito_ribo(obj)
  obj <- assign_pain_group(obj, cfg)

  obj$PoolGroup <- assign_pool_group(
    as.character(obj@meta.data[[cfg$sample_col]]), cfg
  )

  message("  Cell cycle scoring (NormalizeData -> CellCycleScoring)...")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(
    obj,
    s.features   = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident    = FALSE,
    verbose      = FALSE
  )
  message("  Phase: ",
          paste(names(table(obj$Phase)), table(obj$Phase),
                sep = "=", collapse = " | "))

  print(table(obj@meta.data[[cfg$sample_col]], useNA = "ifany"))
  obj
}

# =============================================================================
## PCA + (per-dim-spec) UMAP/CLUSTER HELPERS
# =============================================================================

run_pca_phase <- function(obj, cfg, features, plot_dir, label,
                          reduction_pca = "pca") {
  obj <- RunPCA(obj, features = features, npcs = cfg$n_pcs,
                reduction.name = reduction_pca, verbose = FALSE)
  if (!should_skip("elbow", cfg)) {
    p_elbow <- ElbowPlot(obj, ndims = cfg$n_pcs, reduction = reduction_pca) +
      ggtitle(paste0("PCA Elbow -- ", label)) +
      theme_cowplot(13)
    save_png(p_elbow, "01_elbow_pca", plot_dir, 8, 5)
  }
  obj
}

# Runs UMAP + neighbors + the full resolution sweep for ONE dim_spec, writing
# a clustree to dim_plot_dir. Adds reductions and metadata columns onto obj
# under names tagged with the dim_spec so multiple specs coexist.
run_umap_cluster_dimspec <- function(obj, cfg, dim_vec, tag,
                                     reduction_for_umap, graph_name_base,
                                     dim_plot_dir, label) {
  # Seurat's `[[<-` sanitizes graph AND reduction names through make.names()
  # (warning + silent rename), so a tag like "1-9_11-12" gets stored as
  # "1.9_11.12" and the original lookup in FindClusters / DimPlot fails.
  # Pre-sanitize here so the names we pass in match the names Seurat actually
  # stores. The filesystem tag is kept as-is for readability.
  graph_name     <- make.names(paste0(graph_name_base, "_dims_", tag))
  snn_name       <- make.names(paste0(graph_name, "_snn"))
  nn_name        <- make.names(paste0(graph_name, "_nn"))
  graph_prefix   <- paste0(snn_name, "_res.")
  reduction_umap <- make.names(paste0("umap_", tag))

  # Seurat's internal name sanitization for reduction slots is not guaranteed
  # to match make.names() across versions (a mismatch here is what produced the
  # 'umap_<tag>' not-found error downstream). Capture the reductions before the
  # call and read back whatever Seurat actually stored, so downstream DimPlot /
  # FeaturePlot lookups always reference the real slot name.
  reducs_before  <- Reductions(obj)
  obj <- RunUMAP(obj, reduction = reduction_for_umap, dims = dim_vec,
                 reduction.name = reduction_umap, verbose = FALSE)
  new_reducs <- setdiff(Reductions(obj), reducs_before)
  if (!reduction_umap %in% Reductions(obj)) {
    if (length(new_reducs) >= 1) {
      message("  NOTE: UMAP reduction stored as '", new_reducs[length(new_reducs)],
              "' (requested '", reduction_umap, "'); using stored name.")
      reduction_umap <- new_reducs[length(new_reducs)]
    } else {
      stop("RunUMAP did not add a reduction for dims tag '", tag, "'.")
    }
  }
  obj <- FindNeighbors(obj, reduction = reduction_for_umap, dims = dim_vec,
                       graph.name = c(nn_name, snn_name), verbose = FALSE)
  for (res in cfg$resolutions) {
    obj <- FindClusters(obj, graph.name = snn_name, resolution = res,
                        verbose = FALSE)
    message("    [dims=", tag, "  res=", res, "] -> ",
            length(unique(obj@meta.data[[paste0(graph_prefix, res)]])),
            " clusters")
  }

  if (!should_skip("clustree", cfg)) {
    ct_p <- tryCatch({
      p <- clustree(obj@meta.data, prefix = graph_prefix) +
        ggtitle(paste0("Clustree -- ", label, "  |  dims=", tag)) +
        theme(legend.position = "right")
      add_clustree_res_marker(p, cfg$resolutions, cfg$final_res)
    }, error = function(e) {
      message("  Clustree failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(ct_p))
      save_png(ct_p, "02_clustree", dim_plot_dir, 12, 14)
  }

  list(
    obj            = obj,
    dim_tag        = tag,
    dim_vec        = dim_vec,
    graph_prefix   = graph_prefix,
    reduction_umap = reduction_umap,
    graph_name     = graph_name
  )
}

# =============================================================================
## VARIANT PREP (shared once per variant; the dim-spec loop runs after)
##
## Each prep returns:
##   $obj          fully prepped Seurat object (HVG/SCT/integration/PCA all done)
##   $variant_meta list with:
##     variant            : variant id
##     reduction_for_umap : reduction passed to RunUMAP/FindNeighbors per spec
##                          ("pca" for merged/sctv2/rpca, "harmony" for harmony)
##     reduction_pca      : reduction used as the "PC" source for overlay plots
##     graph_name_base    : prefix for graph/cluster column names
##     marker_assay/layer : assay+layer for FindAllMarkers and DoHeatmap
##     label              : human-readable label
## =============================================================================

process_merged_prep <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[merged] log-normalize -> ScaleData -> PCA")
  DefaultAssay(obj) <- cfg$rna_assay
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = cfg$n_hvg_logn, verbose = FALSE)
  cleaned <- filter_hvg(VariableFeatures(obj), cfg$exclude_patterns)
  VariableFeatures(obj) <- cleaned
  message("  HVGs: ", cfg$n_hvg_logn, " -> ", length(cleaned),
          " after stripping MT/ribo/TCR/Ig")
  write.csv(data.frame(gene = cleaned),
            file.path(data_dir, "merged_hvg_list.csv"), row.names = FALSE)
  obj <- ScaleData(obj, features = cleaned,
                   vars.to.regress = cfg$regress_vars, verbose = FALSE)
  obj <- run_pca_phase(obj, cfg, features = cleaned,
                       plot_dir = plot_dir, label = "Merged (log-norm)",
                       reduction_pca = "pca")
  list(obj = obj,
       variant_meta = list(
         variant            = "merged",
         reduction_for_umap = "pca",
         reduction_pca      = "pca",
         graph_name_base    = "RNA",
         marker_assay       = cfg$rna_assay,
         marker_layer       = "data",
         label              = "Merged (log-norm)"
       ))
}

run_sct <- function(obj, cfg) {
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
          " after stripping MT/ribo/TCR/Ig")
  list(obj = obj, hvg = cleaned_hvg)
}

process_sctv2_prep <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[sctv2] SCTransform v2 -> PCA")
  DefaultAssay(obj) <- cfg$rna_assay
  sct <- run_sct(obj, cfg)
  obj <- sct$obj
  write.csv(data.frame(gene = sct$hvg),
            file.path(data_dir, "sctv2_hvg_list.csv"), row.names = FALSE)
  obj <- run_pca_phase(obj, cfg, features = sct$hvg,
                       plot_dir = plot_dir, label = "SCTransform v2",
                       reduction_pca = "pca")
  list(obj = obj,
       variant_meta = list(
         variant            = "sctv2",
         reduction_for_umap = "pca",
         reduction_pca      = "pca",
         graph_name_base    = "SCT",
         marker_assay       = "SCT",
         marker_layer       = "data",
         label              = "SCTransform v2"
       ))
}

# Integration runs once with the full PC range (1:n_pcs); the dim-spec loop
# selects subsets downstream. Re-running integration per spec would multiply
# the most expensive step by length(specs); the user-selected dim_spec only
# changes UMAP/Neighbors/Clusters anyway.
process_rpca_prep <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[rpca] SCTransform v2 per ", cfg$pool_col,
          " -> RPCA integration (dims=1:", cfg$n_pcs, ") -> PCA")
  DefaultAssay(obj) <- cfg$rna_assay

  obj_list <- SplitObject(obj, split.by = cfg$pool_col)
  message("  Pools: ",
          paste(names(obj_list), vapply(obj_list, ncol, integer(1)),
                sep = "=", collapse = " | "))

  obj_list <- lapply(seq_along(obj_list), function(i) {
    nm <- names(obj_list)[i]
    message("  [pool=", nm, "] SCTransform v2 (cells=", ncol(obj_list[[i]]), ")")
    x <- obj_list[[i]]
    DefaultAssay(x) <- cfg$rna_assay
    SCTransform(
      x,
      vst.flavor          = "v2",
      method              = "glmGamPoi",
      vars.to.regress     = cfg$regress_vars,
      variable.features.n = cfg$n_hvg,
      verbose             = FALSE
    )
  }) |> setNames(names(obj_list))

  features <- SelectIntegrationFeatures(
    obj_list, nfeatures = cfg$n_hvg, verbose = FALSE
  )
  features <- filter_hvg(features, cfg$exclude_patterns)
  message("  Integration features: ", length(features),
          " after stripping MT/ribo/TCR/Ig")
  write.csv(data.frame(gene = features),
            file.path(data_dir, "rpca_integration_features.csv"),
            row.names = FALSE)

  obj_list <- PrepSCTIntegration(
    obj_list, anchor.features = features, verbose = FALSE
  )
  obj_list <- lapply(obj_list, function(x) {
    RunPCA(x, features = features, npcs = cfg$n_pcs, verbose = FALSE)
  })

  min_cells <- min(vapply(obj_list, ncol, integer(1)))
  k_filter  <- min(200, max(5, min_cells - 1))
  full_dims <- seq_len(cfg$n_pcs)

  anchors <- FindIntegrationAnchors(
    object.list          = obj_list,
    normalization.method = "SCT",
    anchor.features      = features,
    reduction            = "rpca",
    dims                 = full_dims,
    k.anchor             = cfg$rpca_k_anchor,
    k.filter             = k_filter,
    verbose              = TRUE
  )

  obj_int <- IntegrateData(
    anchorset            = anchors,
    normalization.method = "SCT",
    dims                 = full_dims,
    verbose              = TRUE
  )

  meta_keep <- setdiff(colnames(obj@meta.data), colnames(obj_int@meta.data))
  if (length(meta_keep) > 0) {
    common_cells <- intersect(colnames(obj_int), colnames(obj))
    obj_int <- AddMetaData(
      obj_int, metadata = obj@meta.data[common_cells, meta_keep, drop = FALSE]
    )
  }

  if (cfg$adt_assay %in% names(obj@assays) &&
      !cfg$adt_assay %in% names(obj_int@assays)) {
    message("  Re-attaching ", cfg$adt_assay, " assay onto integrated object")
    common_cells <- intersect(colnames(obj_int), colnames(obj))
    adt_sub <- tryCatch(
      subset(obj[[cfg$adt_assay]], cells = common_cells),
      error = function(e) {
        adt_counts <- safe_counts(obj, cfg$adt_assay)
        CreateAssayObject(counts = adt_counts[, common_cells, drop = FALSE])
      }
    )
    obj_int[[cfg$adt_assay]] <- adt_sub
  }

  DefaultAssay(obj_int) <- "integrated"
  obj_int <- run_pca_phase(obj_int, cfg, features = NULL,
                           plot_dir = plot_dir, label = "SCTv2 + RPCA",
                           reduction_pca = "pca")

  list(obj = obj_int,
       variant_meta = list(
         variant            = "rpca",
         reduction_for_umap = "pca",
         reduction_pca      = "pca",
         graph_name_base    = "integrated",
         marker_assay       = "SCT",
         marker_layer       = "data",
         label              = "SCTv2 + RPCA"
       ))
}

# Harmony correction is run once with dims.use = 1:n_pcs. Each dim_spec then
# picks a subset of corrected components for UMAP/clustering.
process_harmony_prep <- function(obj, cfg, plot_dir, data_dir) {
  message("\n[harmony] SCTransform v2 -> PCA -> RunHarmony(",
          cfg$pool_col, "), dims.use=1:", cfg$n_pcs)
  DefaultAssay(obj) <- cfg$rna_assay
  sct <- run_sct(obj, cfg)
  obj <- sct$obj
  write.csv(data.frame(gene = sct$hvg),
            file.path(data_dir, "harmony_hvg_list.csv"), row.names = FALSE)
  obj <- run_pca_phase(obj, cfg, features = sct$hvg,
                       plot_dir = plot_dir, label = "SCTv2 + Harmony",
                       reduction_pca = "pca")
  obj <- RunHarmony(
    obj,
    group.by.vars  = cfg$pool_col,
    reduction      = "pca",
    dims.use       = seq_len(cfg$n_pcs),
    reduction.save = "harmony",
    verbose        = TRUE
  )
  list(obj = obj,
       variant_meta = list(
         variant            = "harmony",
         reduction_for_umap = "harmony",
         # Overlay Harmony-corrected components (these drive clustering for
         # the harmony variant); naming will say "Harmony N".
         reduction_pca      = "harmony",
         graph_name_base    = "SCT_harmony",
         marker_assay       = "SCT",
         marker_layer       = "data",
         label              = "SCTv2 + Harmony"
       ))
}

# =============================================================================
## FIGURE STAGES (per-combo; combo metadata stamped onto obj as attribute)
## =============================================================================

stage_umaps <- function(obj, cfg, plot_dir, label) {
  if (should_skip("umaps", cfg)) {
    message("  [skip] umaps"); return(invisible(NULL))
  }
  cm <- attr(obj, "variant_meta")
  reduc <- cm$reduction_umap

  samples  <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  pools    <- sort(unique(as.character(obj@meta.data[[cfg$pool_col]])))
  samp_pal <- make_sample_pal(samples)
  pool_pal <- make_pool_pal(pools, cfg)
  pain_pal <- c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")

  p_clust <- DimPlot(obj, reduction = reduc, label = TRUE,
                     repel = TRUE, label.size = 5, raster = FALSE) +
    NoLegend() +
    labs(title = paste0(label, " -- clusters"),
         subtitle = umap_subtitle(obj)) +
    theme_cowplot(13)
  save_png(p_clust, "03a_umap_clusters", plot_dir, cfg$umap_w, cfg$umap_h)

  p_pool <- DimPlot(obj, reduction = reduc, group.by = cfg$pool_col,
                    cols = pool_pal, pt.size = 0.3, raster = FALSE) +
    labs(title = paste0(label, " -- by ", cfg$pool_col),
         subtitle = umap_subtitle(obj, show_res = FALSE)) +
    theme_cowplot(13)
  save_png(p_pool, "03b_umap_by_pool", plot_dir, cfg$umap_w + 2, cfg$umap_h)

  p_imp <- DimPlot(obj, reduction = reduc, group.by = cfg$sample_col,
                   cols = samp_pal, pt.size = 0.3, raster = FALSE) +
    labs(title = paste0(label, " -- by IMP sample"),
         subtitle = umap_subtitle(obj, show_res = FALSE)) +
    theme_cowplot(13)
  save_png(p_imp, "03c_umap_by_IMP_sample", plot_dir,
           cfg$umap_w + 4, cfg$umap_h)

  if ("PainGroup" %in% colnames(obj@meta.data) &&
      any(!is.na(obj@meta.data$PainGroup))) {
    p_pain <- DimPlot(obj, reduction = reduc, group.by = "PainGroup",
                      cols = pain_pal, pt.size = 0.3, raster = FALSE,
                      na.value = "grey90") +
      labs(title = paste0(label, " -- by Pain Group"),
           subtitle = umap_subtitle(obj, show_res = FALSE)) +
      theme_cowplot(13)
    save_png(p_pain, "03d_umap_by_pain_group", plot_dir,
             cfg$umap_w + 2, cfg$umap_h)
  }

  invisible(NULL)
}

stage_qc <- function(obj, cfg, plot_dir, label) {
  if (should_skip("qc", cfg)) {
    message("  [skip] qc"); return(invisible(NULL))
  }
  qc_feats <- intersect(
    c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_Ribo"),
    colnames(obj@meta.data)
  )
  if (length(qc_feats) == 0) {
    message("  No QC features found -- skipping QC stage.")
    return(invisible(NULL))
  }
  n_feat <- length(qc_feats)

  plist <- VlnPlot(obj, features = qc_feats, group.by = "seurat_clusters",
                   pt.size = 0, combine = FALSE)
  for (i in seq_along(plist))
    plist[[i]] <- plist[[i]] +
      ggtitle(qc_feats[i]) +
      theme(legend.position = "none",
            plot.title   = element_text(face = "bold", size = 17),
            axis.title.x = element_blank(),
            axis.text.x  = if (i < n_feat) element_blank()
                           else element_text(angle = 45, hjust = 1, size = 10),
            axis.ticks.x = if (i < n_feat) element_blank() else element_line())
  p <- wrap_plots(plist, ncol = 1) +
    plot_annotation(title = paste0("QC by cluster -- ", label))
  save_png(p, "04a_qc_violin_by_cluster", plot_dir,
           width  = max(12, length(levels(obj$seurat_clusters)) * 0.45),
           height = 3.2 * n_feat)

  samples <- sort(unique(as.character(obj@meta.data[[cfg$sample_col]])))
  pools_per_sample <- obj@meta.data |>
    as.data.frame() |>
    dplyr::transmute(Sample = as.character(.data[[cfg$sample_col]]),
                     Pool   = as.character(.data[[cfg$pool_col]])) |>
    dplyr::count(Sample, Pool) |>
    dplyr::group_by(Sample) |>
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  pool_pal <- make_pool_pal(pools_per_sample$Pool, cfg)
  sample_cols <- setNames(
    ifelse(pools_per_sample$Pool %in% names(pool_pal),
           pool_pal[pools_per_sample$Pool], "grey60"),
    pools_per_sample$Sample
  )
  sample_order <- pools_per_sample$Sample[
    order(pools_per_sample$Pool, pools_per_sample$Sample)
  ]
  obj@meta.data[[cfg$sample_col]] <- factor(
    as.character(obj@meta.data[[cfg$sample_col]]), levels = sample_order
  )
  plist2 <- VlnPlot(obj, features = qc_feats, group.by = cfg$sample_col,
                    pt.size = 0, combine = FALSE, cols = sample_cols)
  for (i in seq_along(plist2))
    plist2[[i]] <- plist2[[i]] +
      ggtitle(qc_feats[i]) +
      theme(legend.position = "none",
            plot.title   = element_text(face = "bold", size = 17),
            axis.title.x = element_blank(),
            axis.text.x  = if (i < n_feat) element_blank()
                           else element_text(angle = 45, hjust = 1, size = 10),
            axis.ticks.x = if (i < n_feat) element_blank() else element_line())
  p2 <- wrap_plots(plist2, ncol = 1) +
    plot_annotation(title = paste0("QC by sample (colored by pool) -- ", label))
  save_png(p2, "04b_qc_violin_by_IMP_sample", plot_dir,
           width = max(14, length(samples) * 0.5), height = 3.2 * n_feat)

  invisible(NULL)
}

stage_adt_featureplots <- function(obj, cfg, plot_dir, label,
                                   per_page = 9, ncol = 3) {
  if (should_skip("adt", cfg)) {
    message("  [skip] adt"); return(invisible(NULL))
  }
  if (!(cfg$adt_assay %in% names(obj@assays))) {
    message("  ", cfg$adt_assay, " assay not present -- skipping ADT FeaturePlots.")
    return(invisible(NULL))
  }
  cm <- attr(obj, "variant_meta")
  reduc <- cm$reduction_umap
  prev_assay <- DefaultAssay(obj)
  DefaultAssay(obj) <- cfg$adt_assay

  feats <- names(cfg$adt_markers)
  feats <- feats[feats %in% rownames(obj[[cfg$adt_assay]])]
  if (length(feats) == 0) {
    message("  No requested ADT markers found -- skipping.")
    DefaultAssay(obj) <- prev_assay
    return(invisible(NULL))
  }

  chunks <- split(feats, ceiling(seq_along(feats) / per_page))
  adt_limits <- c(0, 3)

  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    plist <- lapply(seq_along(chunk), function(j) {
      f       <- chunk[[j]]
      display <- cfg$adt_markers[[f]]
      FeaturePlot(obj, features = f, reduction = reduc, raster = FALSE) +
        scale_color_gradientn(
          colours = c("lightgrey", "#4B0082"),
          limits  = adt_limits, oob = scales::squish,
          name    = "Surface\nMarker\nIntensity"
        ) +
        ggtitle(display) +
        theme_cowplot(11) +
        theme(plot.title = element_text(face = "bold", size = 15),
              legend.position = "right")
    })
    panel <- wrap_plots(plist, ncol = ncol) +
      plot_annotation(title = paste0("ADT FeaturePlots (", i, "/",
                                     length(chunks), ") -- ", label),
                      subtitle = umap_subtitle(obj, show_res = FALSE))
    save_png(panel, paste0("05_adt_featureplots_p", i), plot_dir,
             width = 4.5 * ncol, height = 4.5 * ceiling(length(chunk) / ncol))
  }

  DefaultAssay(obj) <- prev_assay
  invisible(NULL)
}

stage_cluster_composition <- function(obj, cfg, plot_dir, data_dir, label) {
  if (should_skip("composition", cfg)) {
    message("  [skip] composition"); return(invisible(NULL))
  }
  clust_vec <- as.character(obj$seurat_clusters)
  samp_vec  <- as.character(obj@meta.data[[cfg$sample_col]])
  keep      <- !is.na(clust_vec) & !is.na(samp_vec)
  tbl       <- as.data.frame(
    table(seurat_clusters = clust_vec[keep],
          sample_id       = samp_vec[keep], useNA = "no"),
    stringsAsFactors = FALSE
  )
  colnames(tbl) <- c("seurat_clusters", cfg$sample_col, "n")
  tbl$n <- as.integer(tbl$n)

  prop_df <- tbl |>
    dplyr::group_by(.data[[cfg$sample_col]]) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()
  write.csv(prop_df, file.path(data_dir, "cluster_composition_by_sample.csv"),
            row.names = FALSE)

  clust_lvls <- levels(factor(suppressWarnings(
    as.integer(as.character(prop_df$seurat_clusters)))))
  if (anyNA(clust_lvls) || length(clust_lvls) == 0)
    clust_lvls <- sort(unique(as.character(prop_df$seurat_clusters)))
  prop_df$seurat_clusters <- factor(prop_df$seurat_clusters, levels = clust_lvls)
  clust_pal <- setNames(
    colorRampPalette(c(brewer.pal(8, "Set2"),
                       brewer.pal(8, "Set1")))(length(clust_lvls)),
    clust_lvls
  )

  samp_pool <- obj@meta.data |>
    as.data.frame() |>
    dplyr::transmute(Sample = as.character(.data[[cfg$sample_col]]),
                     Pool   = as.character(.data[[cfg$pool_col]])) |>
    dplyr::count(Sample, Pool) |>
    dplyr::group_by(Sample) |>
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  samp_order <- samp_pool$Sample[order(samp_pool$Pool, samp_pool$Sample)]
  prop_df[[cfg$sample_col]] <- factor(
    prop_df[[cfg$sample_col]], levels = samp_order
  )

  p <- ggplot(prop_df,
              aes(x = .data[[cfg$sample_col]], y = prop,
                  fill = seurat_clusters)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = clust_pal, name = "Cluster") +
    labs(title = paste0("Cluster Composition by IMP Sample -- ", label),
         subtitle = umap_subtitle(obj),
         x = NULL, y = "Proportion") +
    theme_cowplot(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  save_png(p, "06_cluster_composition_by_IMP_sample", plot_dir,
           width = max(14, length(samp_order) * 0.55), height = 7)
  invisible(NULL)
}

stage_markers <- function(obj, cfg, plot_dir, data_dir, label) {
  if (should_skip("markers", cfg)) {
    message("  [skip] markers"); return(invisible(NULL))
  }
  cm <- attr(obj, "variant_meta")
  assay_use <- cm$marker_assay
  if (!(assay_use %in% names(obj@assays))) {
    message("  Assay '", assay_use, "' missing -- skipping markers.")
    return(invisible(NULL))
  }
  DefaultAssay(obj) <- assay_use

  if (assay_use == "SCT")
    obj <- tryCatch(PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE),
                    error = function(e) {
                      message("  PrepSCTFindMarkers failed: ",
                              conditionMessage(e)); obj
                    })

  Idents(obj) <- obj$seurat_clusters

  markers <- FindAllMarkers(
    obj, assay = assay_use, only.pos = TRUE,
    min.pct = 0.05, logfc.threshold = 0.2,
    test.use = "wilcox", verbose = FALSE
  )
  markers_sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, ]
  write.csv(markers,     file.path(data_dir, "cluster_markers_all.csv"),
            row.names = FALSE)
  write.csv(markers_sig, file.path(data_dir, "cluster_markers_sig.csv"),
            row.names = FALSE)
  message("  Significant cluster markers: ", nrow(markers_sig))

  if (nrow(markers_sig) == 0) {
    message("  No significant markers -- skipping heatmap.")
    return(invisible(markers))
  }

  top5 <- markers_sig |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = 5, with_ties = FALSE) |>
    dplyr::pull(gene) |>
    unique()
  cl_levels <- levels(obj$seurat_clusters)
  Idents(obj) <- factor(as.character(obj$seurat_clusters), levels = cl_levels)

  obj <- tryCatch(
    ScaleData(obj, features = top5, assay = assay_use, verbose = FALSE),
    error = function(e) { message("  ScaleData failed: ",
                                  conditionMessage(e)); obj }
  )
  n_cl <- length(cl_levels)

  p_heat <- DoHeatmap(obj, features = top5, assay = assay_use,
                      slot = "scale.data", size = 5, angle = 45) +
    scale_fill_gradientn(colors = c("#2166AC", "white", "#B2182B"),
                         name = "Scaled\nexpression") +
    theme(axis.text.y = element_text(size = 13, face = "italic")) +
    ggtitle(paste0("De novo top-5 markers per cluster -- ", label))
  save_png_device(
    p_heat, "07_denovo_top5_marker_heatmap", plot_dir,
    width  = max(12, n_cl * 0.9),
    height = min(40, max(10, length(top5) * 0.34))
  )

  invisible(markers_sig)
}

# NEW: facet-grid UMAP where each panel is a FeaturePlot of one PC's (or
# Harmony component's) cell-score, paged into groups of `per_page`. Points are
# enlarged via cfg$pc_overlay_pt_size.
stage_pc_overlay_facet <- function(obj, cfg, plot_dir, label,
                                   per_page = 4, ncol = 2) {
  if (should_skip("pc_overlay", cfg)) {
    message("  [skip] pc_overlay"); return(invisible(NULL))
  }
  cm <- attr(obj, "variant_meta")
  reduction_pca  <- cm$reduction_pca
  reduction_umap <- cm$reduction_umap
  dim_vec        <- cm$dim_vec

  if (!(reduction_pca %in% names(obj@reductions))) {
    message("  Reduction '", reduction_pca,
            "' missing -- skipping PC overlay.")
    return(invisible(NULL))
  }
  emb <- Embeddings(obj, reduction = reduction_pca)
  if (max(dim_vec) > ncol(emb)) {
    message("  WARN: dim_vec includes component ", max(dim_vec),
            " but '", reduction_pca, "' has only ", ncol(emb),
            " -- clipping.")
    dim_vec <- dim_vec[dim_vec <= ncol(emb)]
  }
  if (length(dim_vec) == 0) return(invisible(NULL))

  prefix <- if (identical(reduction_pca, "harmony")) "Harmony" else "PC"
  meta_cols <- paste0("__", prefix, dim_vec, "_score")
  for (i in seq_along(dim_vec)) {
    obj@meta.data[[meta_cols[i]]] <- emb[, dim_vec[i]]
  }

  chunks <- split(seq_along(dim_vec),
                  ceiling(seq_along(dim_vec) / per_page))

  for (k in seq_along(chunks)) {
    idxs <- chunks[[k]]
    plist <- lapply(idxs, function(j) {
      f         <- meta_cols[j]
      pc_label  <- paste0(prefix, " ", dim_vec[j])
      vals      <- obj@meta.data[[f]]
      vlim      <- max(abs(range(vals, na.rm = TRUE)))
      FeaturePlot(obj, features = f, reduction = reduction_umap,
                  raster = FALSE, pt.size = cfg$pc_overlay_pt_size,
                  order = TRUE) +
        scale_color_gradient2(
          low      = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
          midpoint = 0,
          limits   = c(-vlim, vlim), oob = scales::squish,
          name     = pc_label
        ) +
        ggtitle(pc_label) +
        theme_cowplot(11) +
        theme(plot.title = element_text(face = "bold", size = 15),
              legend.position = "right")
    })
    panel <- wrap_plots(plist, ncol = ncol) +
      plot_annotation(
        title    = paste0("PC contribution overlay (", k, "/",
                          length(chunks), ") -- ", label),
        subtitle = umap_subtitle(obj)
      )
    save_png(panel, paste0("08_pc_overlay_p", k), plot_dir,
             width  = 6   * ncol,
             height = 5.5 * ceiling(length(idxs) / ncol))
  }

  invisible(NULL)
}

make_variant_figures <- function(obj, cfg, plot_dir, data_dir, label) {
  message("\n[figures] ", label)
  stage_umaps(obj, cfg, plot_dir, label)
  stage_qc(obj, cfg, plot_dir, label)
  stage_adt_featureplots(obj, cfg, plot_dir, label)
  stage_cluster_composition(obj, cfg, plot_dir, data_dir, label)
  stage_pc_overlay_facet(obj, cfg, plot_dir, label)
  stage_markers(obj, cfg, plot_dir, data_dir, label)
}

# =============================================================================
## CROSS-VARIANT COMPARISON (per combo)
# =============================================================================
make_variant_comparison <- function(results, cfg, comp_dir, dim_tag, res_val) {
  if (should_skip("compare", cfg)) {
    message("  [skip] compare"); return(invisible(NULL))
  }
  message("\n[compare] dims=", dim_tag, "  res=", res_val,
          "  variants=", paste(names(results), collapse = ","))
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  variant_labels <- c(
    merged  = "Merged (log-norm)",
    sctv2   = "SCTv2",
    rpca    = "SCTv2 + RPCA",
    harmony = "SCTv2 + Harmony"
  )

  samples <- sort(unique(unlist(lapply(results, function(o)
    as.character(o@meta.data[[cfg$sample_col]])))))
  pools   <- sort(unique(unlist(lapply(results, function(o)
    as.character(o@meta.data[[cfg$pool_col]])))))
  samp_pal <- make_sample_pal(samples)
  pool_pal <- make_pool_pal(pools, cfg)
  pain_pal <- c("Low Pain" = "#6BAED6", "High Pain" = "#D95F5F")

  subtitle_str <- paste0("dims = ", dim_tag, "  |  res = ", res_val)

  make_panel <- function(group_col, palette, na_value = "grey90",
                         title, file_tag) {
    plist <- lapply(names(results), function(v) {
      o    <- results[[v]]
      cm   <- attr(o, "variant_meta")
      DimPlot(o, reduction = cm$reduction_umap, group.by = group_col,
              cols = palette, pt.size = 0.3, raster = FALSE,
              na.value = na_value) +
        labs(title = variant_labels[v]) +
        theme_cowplot(11) +
        theme(plot.title = element_text(face = "bold", size = 14))
    })
    p <- wrap_plots(plist, ncol = 2) +
      plot_layout(guides = "collect") +
      plot_annotation(title = title, subtitle = subtitle_str)
    save_png(p, file_tag, comp_dir,
             width = cfg$umap_w * 2 + 2, height = cfg$umap_h * 2)
  }

  make_panel(cfg$pool_col, pool_pal, "grey90",
             paste0("UMAP comparison -- by ", cfg$pool_col),
             "00_compare_umap_by_pool")
  make_panel(cfg$sample_col, samp_pal, "grey90",
             "UMAP comparison -- by IMP sample",
             "01_compare_umap_by_IMP_sample")
  if (all(vapply(results, function(o)
        "PainGroup" %in% colnames(o@meta.data), logical(1)))) {
    make_panel("PainGroup", pain_pal, "grey90",
               "UMAP comparison -- by Pain Group",
               "02_compare_umap_by_pain_group")
  }

  summary_df <- tibble::tibble(
    variant = names(results),
    label   = variant_labels[names(results)],
    dim_tag = dim_tag,
    res     = res_val,
    n_cells = vapply(results, ncol, integer(1)),
    n_clusters = vapply(results, function(o)
      length(levels(droplevels(o$seurat_clusters))), integer(1))
  )
  write.csv(summary_df, file.path(comp_dir, "variant_summary.csv"),
            row.names = FALSE)
  message("  Variant summary written.")
  invisible(NULL)
}

# =============================================================================
## MAIN
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Loading + prep")
message(paste(rep("=", 60), collapse = ""))
obj_base <- load_and_prep(cfg)

variant_results <- list()  # named by variant; each: list(obj, dim_results, variant_meta)

for (v in active_variants) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("VARIANT: ", v)
  message(paste(rep("=", 60), collapse = ""))

  v_root <- file.path(cfg$out_dir, paste0("variant_", v))
  v_plot <- file.path(v_root, "plots")   # variant-shared (elbow)
  v_data <- file.path(v_root, "data")    # variant-shared (HVGs, etc.)
  v_ckpt <- file.path(v_root, paste0("ckpt_", v, ".rds"))
  for (d in c(v_root, v_plot, v_data))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  prep <- switch(
    v,
    merged  = process_merged_prep(obj_base,  cfg, v_plot, v_data),
    sctv2   = process_sctv2_prep(obj_base,   cfg, v_plot, v_data),
    rpca    = process_rpca_prep(obj_base,    cfg, v_plot, v_data),
    harmony = process_harmony_prep(obj_base, cfg, v_plot, v_data)
  )
  obj_v <- prep$obj
  vm    <- prep$variant_meta

  dim_results <- list()
  for (spec in cfg$pca_dims_specs) {
    tag    <- format_dim_tag(spec)
    d_root <- file.path(v_root, paste0("dims_", tag))
    d_plot <- file.path(d_root, "plots")
    dir.create(d_plot, recursive = TRUE, showWarnings = FALSE)

    message("\n  [", v, "  dims=", tag, "] UMAP + neighbors + cluster sweep")
    dr <- run_umap_cluster_dimspec(
      obj_v, cfg, spec, tag,
      reduction_for_umap = vm$reduction_for_umap,
      graph_name_base    = vm$graph_name_base,
      dim_plot_dir       = d_plot,
      label              = vm$label
    )
    obj_v <- dr$obj
    dr$obj <- NULL  # don't duplicate the whole object inside dim_results
    dim_results[[tag]] <- dr
  }

  saveRDS(obj_v, v_ckpt)
  message("  Checkpoint: ", basename(v_ckpt))

  # Per-(dim_spec x resolution) figures
  for (tag in names(dim_results)) {
    dr <- dim_results[[tag]]
    for (r in cfg$resolutions) {
      combo_root <- file.path(v_root, paste0("dims_", tag),
                              paste0("res_", format(r, nsmall = 1)))
      combo_plot <- file.path(combo_root, "plots")
      combo_data <- file.path(combo_root, "data")
      for (d in c(combo_plot, combo_data))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)

      cluster_col <- paste0(dr$graph_prefix, r)
      if (!cluster_col %in% colnames(obj_v@meta.data)) {
        message("  [skip] Missing cluster column: ", cluster_col)
        next
      }
      obj_combo <- obj_v
      obj_combo$seurat_clusters <- factor(obj_combo@meta.data[[cluster_col]])
      Idents(obj_combo) <- obj_combo$seurat_clusters

      attr(obj_combo, "variant_meta") <- list(
        variant        = vm$variant,
        dim_tag        = tag,
        dim_vec        = dr$dim_vec,
        graph_prefix   = dr$graph_prefix,
        reduction_umap = dr$reduction_umap,
        reduction_pca  = vm$reduction_pca,
        marker_assay   = vm$marker_assay,
        marker_layer   = vm$marker_layer,
        res            = r,
        label          = vm$label
      )

      message("\n  [", v, "  dims=", tag, "  res=", r,
              "] figures (", length(levels(obj_combo$seurat_clusters)),
              " clusters)")
      make_variant_figures(
        obj_combo, cfg, combo_plot, combo_data,
        paste0(vm$label, "  |  dims=", tag, "  |  res=", r)
      )
    }
  }

  variant_results[[v]] <- list(
    obj          = obj_v,
    dim_results  = dim_results,
    variant_meta = vm
  )
}

# Cross-variant comparison per combo
if (length(variant_results) > 1) {
  ref_tags <- names(variant_results[[1]]$dim_results)
  for (tag in ref_tags) {
    for (r in cfg$resolutions) {
      comp_root <- file.path(
        cfg$out_dir, "variant_comparison",
        paste0("dims_", tag), paste0("res_", format(r, nsmall = 1))
      )
      objs_for_combo <- list()
      for (v in names(variant_results)) {
        vr <- variant_results[[v]]
        if (!(tag %in% names(vr$dim_results))) next
        dr <- vr$dim_results[[tag]]
        cluster_col <- paste0(dr$graph_prefix, r)
        if (!cluster_col %in% colnames(vr$obj@meta.data)) next
        o <- vr$obj
        o$seurat_clusters <- factor(o@meta.data[[cluster_col]])
        Idents(o) <- o$seurat_clusters
        attr(o, "variant_meta") <- list(
          variant        = v,
          dim_tag        = tag,
          dim_vec        = dr$dim_vec,
          graph_prefix   = dr$graph_prefix,
          reduction_umap = dr$reduction_umap,
          reduction_pca  = vr$variant_meta$reduction_pca,
          marker_assay   = vr$variant_meta$marker_assay,
          marker_layer   = vr$variant_meta$marker_layer,
          res            = r,
          label          = vr$variant_meta$label
        )
        objs_for_combo[[v]] <- o
      }
      if (length(objs_for_combo) > 1) {
        tryCatch(
          make_variant_comparison(objs_for_combo, cfg, comp_root, tag, r),
          error = function(e)
            message("  [compare dims=", tag, " res=", r, "] failed: ",
                    conditionMessage(e))
        )
      }
    }
  }
}

message("\n", paste(rep("=", 60), collapse = ""))
message("Pipeline complete. Outputs under: ", cfg$out_dir)
for (v in names(variant_results)) {
  message("  - variant_", v, "/")
  for (tag in names(variant_results[[v]]$dim_results))
    message("      dims_", tag, "/  (", length(cfg$resolutions), " resolutions)")
}
if (length(variant_results) > 1)
  message("  - variant_comparison/dims_<tag>/res_<r>/")
message(paste(rep("=", 60), collapse = ""))
