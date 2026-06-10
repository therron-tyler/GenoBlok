#!/usr/bin/env Rscript

.libPaths(c("/path/to/your/R_library", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
#  T CELL ADT PANEL
#  Curated from the NPSLE CITEseq panel. Grouped by functional category.
#  B cell, myeloid, platelet, and isotype controls excluded.
# =============================================================================
TCELL_ADT_MARKERS <- c(
  # Core identity / TCR
  "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8",
  "Hu.CD2",       "Hu.CD5",         "Hu.CD7",
  "Hu.TCR.AB",    "Hu.TCR.Vd2",     "Hu.TCR.Va7.2",
  # Naive / Memory
  "Hu.CD45RA",  "Hu.CD45RO",  "Hu.CD62L",
  "Hu.CD27",    "Hu.CD28",    "Hu.CD95",
  "Hu.CD127",   "Hu.CD122",   "HuMs.CD44",
  # Treg
  "Hu.CD25",
  # Effector / cytotoxic
  "Hu.CD57",    "Hu.CD56",    "Hu.CD94",
  "Hu.CX3CR1",  "Hu.KLRG1",  "Hu.CD161",
  # Activation
  "Hu.CD69",   "Hu.CD38-HIT2", "Hu.HLA.DR",
  "Hu.CD26",   "Hu.CD137",     "Hu.CD154",
  # Exhaustion / checkpoint
  "Hu.CD223",  "Hu.CD279",  "Hu.TIGIT",
  "Hu.CD152",  "Hu.CD244",
  # Homing / trafficking
  "Hu.CD183",  "Hu.CD185",  "Hu.CD194",
  "Hu.CD195",  "Hu.CD196",
  # Co-stimulation / misc
  "Hu.CD49d",  "Hu.CD58",   "Hu.CD52"
)

# =============================================================================
#  CLI PARSER
# =============================================================================
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  get_flag <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx == length(args)) stop("Missing value after ", flag)
    args[[idx + 1]]
  }

  rds_path   <- get_flag("--rds")
  outdir     <- get_flag("--outdir")
  dims_s     <- get_flag("--dims",     "10,20,30")
  res_s      <- get_flag("--res",      "0.2,0.5,0.8")
  features_s <- get_flag("--features",
    "CD3D,CD4,CD8A,CCR7,SELL,IL7R,GZMB,GZMK,FOXP3,MKI67,LAG3,PDCD1,ISG15,CX3CR1,NKG7")

  if (is.null(rds_path) || is.null(outdir)) {
    stop(
      "Usage: Rscript SingleCell_RDS_list_DimRezScanner_UMAP.R \\\n",
      "  --rds <file.rds | rds_dir/> \\\n",
      "  --outdir <dir> \\\n",
      "  [--dims 7,8,11,12] \\\n",
      "  [--res 0.2,0.4,0.6] \\\n",
      "  [--features CD3D,CD4,...]"
    )
  }

  split_clean <- function(x) {
    x <- trimws(x)
    if (nchar(x) == 0) character(0) else trimws(strsplit(x, ",")[[1]])
  }

  dims_vec    <- as.integer(split_clean(dims_s))
  resolutions <- as.numeric(split_clean(res_s))
  features    <- split_clean(features_s)

  if (anyNA(dims_vec))    stop("Non-integer dims in --dims: ",    dims_s)
  if (anyNA(resolutions)) stop("Non-numeric resolutions in --res: ", res_s)

  list(rds_path = rds_path, outdir = outdir,
       dims_vec = dims_vec, resolutions = resolutions, features = features)
}

# =============================================================================
#  LOADER  —  single .rds Seurat, named list .rds, or directory of .rds files
# =============================================================================
load_rds_input <- function(rds_path) {
  if (dir.exists(rds_path)) {
    files <- list.files(rds_path, pattern = "\\.rds$",
                        full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0) stop("No .rds files in: ", rds_path)
    objs <- lapply(files, readRDS)
    names(objs) <- sub("\\.rds$", "", basename(files), ignore.case = TRUE)
    message("Loaded ", length(objs), " RDS file(s) from directory")
    return(objs)
  }
  if (!file.exists(rds_path)) stop("Path not found: ", rds_path)
  obj <- readRDS(rds_path)
  if (is.list(obj) && !inherits(obj, "Seurat")) {
    if (is.null(names(obj)) || any(names(obj) == ""))
      stop("RDS contains unnamed list. Provide a named list or single Seurat object.")
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
#  HELPER: paginate any list of ggplot panels → one JPEG per page
#  Returns vector of output file paths.
# =============================================================================
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
                paste0(title, "  (page ", pg, "/", length(chunks), ")")
              else title
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

# =============================================================================
#  MAIN WORKER
# =============================================================================
run_umap_grid <- function(
    rds_list,
    outdir,
    dims_vec        = c(10, 20, 30),
    resolutions     = c(0.2, 0.5, 0.8),
    features        = c("CD3D","CD4","CD8A","CCR7","SELL","IL7R",
                        "GZMB","GZMK","FOXP3","MKI67","LAG3","PDCD1",
                        "ISG15","CX3CR1","NKG7"),
    adt_markers     = TCELL_ADT_MARKERS,
    group_by        = NULL,
    label           = TRUE,
    repel           = TRUE,
    dpi             = 300,
    # DimPlot grid layout
    max_per_page    = 9,       # <=9 keeps panels square per page
    panel_size      = 6,       # each DimPlot: panel_size x panel_size inches
    # RNA FeaturePlot
    feature_ncol    = 3,
    feature_pt_size = 3,
    feature_panel_w = 5.5,
    feature_panel_h = 5.0,
    # ADT FeaturePlot
    adt_ncol        = 4,
    adt_pt_size     = 3,
    adt_panel_w     = 5.0,
    adt_panel_h     = 4.5,
    # VlnPlot by cluster
    vln_ncol        = 3,
    vln_panel_w     = 5.0,
    vln_panel_h     = 3.5,
    drop_idents     = c("Negative", "Doublet")
) {
  stopifnot(is.list(rds_list), !is.null(names(rds_list)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  grid_df <- expand.grid(dims = dims_vec, resolution = resolutions,
                         stringsAsFactors = FALSE)

  results <- lapply(names(rds_list), function(nm) {
    message("\n-- Processing: ", nm, " --")
    obj <- rds_list[[nm]]
    DefaultAssay(obj) <- "RNA"

    if (!"pca" %in% names(obj@reductions))
      stop("Object '", nm, "' has no PCA reduction. Run RunPCA() first.")

    max_pcs  <- ncol(obj@reductions$pca)
    bad_dims <- dims_vec[dims_vec > max_pcs]
    if (length(bad_dims) > 0) {
      warning("dims (", paste(bad_dims, collapse=","), ") exceed available PCs (",
              max_pcs, "). Clamping.")
      dims_vec <<- dims_vec[dims_vec <= max_pcs]
    }

    for (bad in drop_idents) {
      if (bad %in% levels(Idents(obj))) {
        obj <- subset(obj, idents = bad, invert = TRUE)
        message("  Dropped ident: ", bad)
      }
    }

    # Check ADT assay once per object
    has_adt <- "ADT" %in% names(obj@assays)
    if (has_adt) {
      adt_present <- intersect(adt_markers, rownames(obj@assays[["ADT"]]))
      adt_missing <- setdiff(adt_markers,   rownames(obj@assays[["ADT"]]))
      message("  ADT assay: ", length(adt_present), " markers present, ",
              length(adt_missing), " missing")
    } else {
      adt_present <- character(0)
      message("  No ADT assay — ADT plots will be skipped")
    }

    ds_dir <- file.path(outdir, nm)
    dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)

    jobs <- lapply(seq_len(nrow(grid_df)), function(i) {
      dims_i   <- grid_df$dims[i]
      res_i    <- grid_df$resolution[i]
      dims_use <- seq_len(dims_i)

      res_tag   <- gsub("\\.", "p", as.character(res_i))
      tag       <- paste0("d", dims_i, "_r", res_tag)
      umap_name <- paste0("umap_", tag)
      snn_name  <- paste0("snn_",  tag)

      message("  dims=1:", dims_i, "  res=", res_i)

      obj2 <- obj
      obj2 <- FindNeighbors(obj2, reduction = "pca", dims = dims_use,
                             graph.name = snn_name, verbose = FALSE)
      obj2 <- FindClusters(obj2, graph.name = snn_name,
                            resolution = res_i, verbose = FALSE)
      obj2 <- RunUMAP(obj2, reduction = "pca", dims = dims_use,
                      reduction.name = umap_name, verbose = FALSE)

      n_clust   <- length(unique(Idents(obj2)))
      combo_lbl <- paste0("dims 1:", dims_i, "  |  res ", res_i,
                          "  |  ", n_clust, " clusters")
      message("    -> ", n_clust, " clusters")

      # ── 1. DimPlot ────────────────────────────────────────────────────────
      p_dim <- DimPlot(obj2, reduction = umap_name, group.by = group_by,
                       label = label, repel = repel) +
        ggtitle(combo_lbl) +
        theme(plot.title = element_text(size = 11))
      ggsave(
        file.path(ds_dir, paste0("UMAP_dims", dims_i, "_res", res_i, ".jpeg")),
        p_dim, width = panel_size, height = panel_size, dpi = dpi
      )

      # ── 2. RNA FeaturePlot ────────────────────────────────────────────────
      genes_present <- intersect(features, rownames(obj2))
      genes_missing <- setdiff(features,   rownames(obj2))
      rna_fp_files  <- NA_character_

      if (length(genes_present) > 0) {
        DefaultAssay(obj2) <- "RNA"
        fp_list <- FeaturePlot(obj2, features = genes_present,
                               reduction = umap_name, pt.size = feature_pt_size,
                               order = TRUE, combine = FALSE, raster = TRUE)
        fp_list <- lapply(seq_along(fp_list), function(k)
          fp_list[[k]] +
            ggtitle(genes_present[[k]]) +
            theme(plot.title = element_text(size = 11, face = "bold")))

        rna_fp_files <- save_panel_pages(
          fp_list,
          base_file = file.path(ds_dir,
                        paste0("FeaturePlots_RNA_dims", dims_i, "_res", res_i, ".jpeg")),
          ncol      = feature_ncol,
          panel_w   = feature_panel_w,
          panel_h   = feature_panel_h,
          dpi       = dpi,
          title     = paste0("RNA Feature Plots  |  ", combo_lbl),
          subtitle  = if (length(genes_missing) > 0)
                        paste0("Missing: ", paste(genes_missing, collapse = ", "))
                      else NULL
        )
      } else {
        warning("No RNA features found for '", nm, "' (", tag, "). Skipping.")
      }

      # ── 3. ADT FeaturePlot ────────────────────────────────────────────────
      adt_fp_files <- NA_character_

      if (has_adt && length(adt_present) > 0) {
        DefaultAssay(obj2) <- "ADT"
        adt_fp_list <- FeaturePlot(obj2, features = adt_present,
                                   reduction = umap_name, pt.size = adt_pt_size,
                                   order = TRUE, combine = FALSE, raster = TRUE)
        # Strip "Hu." / "HuMs." / "HuMsRt." prefixes for cleaner panel titles
        adt_display <- sub("^Hu[A-Za-z]*\\.", "", adt_present)
        adt_fp_list <- lapply(seq_along(adt_fp_list), function(k)
          adt_fp_list[[k]] +
            ggtitle(adt_display[[k]]) +
            theme(plot.title = element_text(size = 10, face = "bold")))

        adt_fp_files <- save_panel_pages(
          adt_fp_list,
          base_file = file.path(ds_dir,
                        paste0("FeaturePlots_ADT_dims", dims_i, "_res", res_i, ".jpeg")),
          ncol      = adt_ncol,
          panel_w   = adt_panel_w,
          panel_h   = adt_panel_h,
          dpi       = dpi,
          title     = paste0("ADT Feature Plots  |  ", combo_lbl)
        )
        DefaultAssay(obj2) <- "RNA"
      }

      # ── 4. VlnPlot by cluster (RNA markers) ──────────────────────────────
      # One panel per gene, violins split by seurat_clusters.
      # Lets you correlate cluster identity with marker expression directly.
      vln_files <- NA_character_

      if (length(genes_present) > 0) {
        DefaultAssay(obj2) <- "RNA"
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

        vln_files <- save_panel_pages(
          vln_list,
          base_file = file.path(ds_dir,
                        paste0("VlnPlot_clusters_dims", dims_i, "_res", res_i, ".jpeg")),
          ncol      = vln_ncol,
          panel_w   = vln_panel_w,
          panel_h   = vln_panel_h,
          dpi       = dpi,
          title     = paste0("RNA Violin by Cluster  |  ", combo_lbl)
        )
      }

      list(
        dims = dims_i, resolution = res_i, tag = tag,
        plot         = p_dim,
        rna_fp_files = rna_fp_files,
        adt_fp_files = adt_fp_files,
        vln_files    = vln_files
      )
    })  # end jobs lapply

    # ── Paginated DimPlot grid (rows=dims, cols=res) ───────────────────────
    plot_lookup <- setNames(lapply(jobs, `[[`, "plot"),
                            vapply(jobs, `[[`, character(1), "tag"))
    ordered <- list()
    for (d in dims_vec)
      for (r in resolutions) {
        key <- paste0("d", d, "_r", gsub("\\.", "p", as.character(r)))
        if (!is.null(plot_lookup[[key]])) ordered[[key]] <- plot_lookup[[key]]
      }

    page_ncol  <- min(3L, length(resolutions))
    chunks     <- split(ordered, ceiling(seq_along(ordered) / max_per_page))
    grid_files <- character(length(chunks))

    for (pg in seq_along(chunks)) {
      ch       <- chunks[[pg]]
      n_cols   <- min(page_ncol, length(ch))
      n_rows   <- ceiling(length(ch) / n_cols)
      pg_label <- if (length(chunks) > 1)
        paste0("  (page ", pg, "/", length(chunks), ")") else ""
      pg_plot  <- wrap_plots(ch, ncol = n_cols) +
        plot_annotation(
          title    = paste0(nm, "  |  UMAP parameter grid", pg_label),
          subtitle = paste0("Rows = dims (", paste(dims_vec, collapse = ", "),
                            ")  |  Cols = res (", paste(resolutions, collapse = ", "), ")")
        )
      suffix    <- if (length(chunks) > 1) paste0("_page", pg) else ""
      grid_file <- file.path(ds_dir, paste0("UMAP_grid", suffix, ".jpeg"))
      ggsave(grid_file, pg_plot,
             width  = n_cols * panel_size,
             height = n_rows * panel_size,
             dpi    = dpi)
      grid_files[pg] <- grid_file
      message("  Grid page ", pg, "/", length(chunks), " saved: ", grid_file)
    }

    list(jobs = jobs, grid_files = grid_files)
  })  # end object lapply

  names(results) <- names(rds_list)
  results
}

# =============================================================================
#  RUN
# =============================================================================
main <- function() {
  a        <- parse_args()
  rds_list <- load_rds_input(a$rds_path)

  umap_out <- run_umap_grid(
    rds_list    = rds_list,
    outdir      = a$outdir,
    dims_vec    = a$dims_vec,
    resolutions = a$resolutions,
    features    = a$features
  )

  index_file <- file.path(a$outdir, "umap_grid_index.txt")
  con <- file(index_file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(umap_out))
    for (gf in umap_out[[nm]]$grid_files)
      writeLines(paste0(nm, "\t", gf), con)

  message("\nDone. Index: ", index_file)
  invisible(umap_out)
}

main()
