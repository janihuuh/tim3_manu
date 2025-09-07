# sc_utils.R
# Lightweight, path-free utilities for Seurat/SingleR/scds and enrichment helpers
# Dependencies: Seurat, SingleCellExperiment, SingleR, scds, dplyr, tibble, purrr,
#               data.table, ggplot2, uwot, stringr
# Suggests: stats

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SingleR)
  library(scds)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(uwot)
  library(stringr)
  library(data.table)
})

# ----------------------------- Small helpers ---------------------------------

extract_pair1 <- function(x) sub("\\_.*", "", x)
extract_pair2 <- function(x) sub(".*\\_", "", x)
extract_name  <- function(x) sub("\\_.*", "", x)
extract_file  <- function(x) sub(".*\\/", "", x)

extract_cluster_number <- function(strs) {
  # from "0 T cells" -> "0"
  vapply(strs, function(s) strsplit(s, " +")[[1]][1], character(1))
}

extract_coarse_phenotype <- function(strs) {
  # from "0 T cells" -> "T"
  vapply(strs, function(s) strsplit(s, " +")[[1]][2], character(1))
}

ordered_timepoints <- function(x, levels_vec) {
  factor(vapply(x, function(s) strsplit(s, "_")[[1]][2], character(1)),
         levels = levels_vec)
}

facets_nice <- ggplot2::theme(
  strip.background = element_rect(fill = "grey96"),
  strip.text = element_text(colour = "black")
)

# ----------------------------- QC --------------------------------------------

get_qc <- function(seurat_object,
                   min_mito = 0, max_mito = 15,
                   min_ribo = 5, max_ribo = 50,
                   min_features = 300, max_features = 5e3,
                   min_counts = 1e3, max_counts = 3e4) {
  md <- seurat_object@meta.data %>% as.data.frame()
  md$barcode <- colnames(seurat_object)

  mito_out   <- md %>% filter(percent.mt   > max_mito | percent.mt   < min_mito) %>% pull(barcode)
  ribo_out   <- md %>% filter(percent.ribo > max_ribo | percent.ribo < min_ribo) %>% pull(barcode)
  feat_out   <- md %>% filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode)
  umi_out    <- md %>% filter(nCount_RNA   > max_counts | nCount_RNA   < min_counts)   %>% pull(barcode)

  outliers   <- c(mito_out, ribo_out, feat_out, umi_out)
  keep       <- setdiff(colnames(seurat_object), outliers)
  subset(seurat_object, cells = keep)
}

# ----------------------------- Preprocess ------------------------------------

preprocess_seurat <- function(orig_object, cells_to_use = NULL,
                              hvg_n = 2000, clip_max = 10,
                              remove_genes = character(),
                              high_var_sd_thresh = 10,
                              npcs = 50,
                              umap_lr = 1) {

  object <- if (!is.null(cells_to_use)) subset(orig_object, cells = cells_to_use) else orig_object

  # normalize + HVGs
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = hvg_n, clip.max = clip_max)

  # prune overly high variance and unwanted genes
  hvg <- VariableFeatures(object)
  too_hvg <- HVFInfo(object) %>%
    rownames_to_column("gene") %>%
    filter(variance.standardized > high_var_sd_thresh) %>%
    pull("gene")

  hvg <- setdiff(hvg, union(too_hvg, remove_genes))
  VariableFeatures(object) <- hvg

  object <- ScaleData(object, features = hvg)
  object <- RunPCA(object, features = hvg, npcs = npcs)

  nPCs <- sum(object[["pca"]]@stdev > 2)
  message(sprintf("Selected PCs with stdev > 2: %d", nPCs))

  # UMAP
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = umap_lr)
  object
}

# ----------------------------- Doublets --------------------------------------

get_doublets <- function(seurat_object) {
  sce <- SingleCellExperiment(list(counts = seurat_object@assays$RNA@counts),
                              colData = seurat_object@meta.data)
  sce <- cxds(sce)
  sce <- bcds(sce)
  sce <- cxds_bcds_hybrid(sce)

  md <- colData(sce)
  seurat_object$cxds_doublet_score <- md$cxds_score
  seurat_object$bcds_doublet_score <- md$bcds_score
  seurat_object$hybrid_doublet_score <- md$hybrid_score

  # simple [0,1] scaling (avoid division by zero)
  scale01 <- function(v) if (diff(range(v)) == 0) rep(0, length(v)) else (v - min(v)) / diff(range(v))
  seurat_object$cxds_doublet_score_norm <- scale01(md$cxds_score)
  seurat_object$bcds_doublet_score_norm <- scale01(md$bcds_score)

  seurat_object
}

# ----------------------------- Clustering ------------------------------------

find_neighbors_all_dims <- function(seurat_object, reduction = "pca") {
  dims <- seq_len(ncol(seurat_object[[reduction]]@cell.embeddings))
  FindNeighbors(seurat_object, reduction = reduction, dims = dims)
}

get_clustering <- function(seurat_object,
                           reduction = "pca",
                           resolutions = c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3),
                           verbose = FALSE) {
  seurat_object <- find_neighbors_all_dims(seurat_object, reduction = reduction)
  FindClusters(seurat_object, resolution = resolutions, verbose = verbose)
}

plot_clustering_overview <- function(seurat_object,
                                     col_prefix = "SNN_res.") {
  cols <- grep(sprintf("^%s", col_prefix), colnames(seurat_object@meta.data), value = TRUE)
  # order by numeric resolution extracted from column names
  res_vals <- as.numeric(sub(paste0("^", col_prefix), "", cols))
  cols <- cols[order(res_vals)]

  n_clusters <- vapply(cols, function(cn) length(levels(seurat_object@meta.data[[cn]])), numeric(1))
  data.frame(resolution = sort(res_vals), nClusters = n_clusters) |>
    ggplot(aes(resolution, nClusters)) +
    geom_point(shape = 21) + theme_bw()
}

# ----------------------------- Latent slots ----------------------------------

put_latents_seurat <- function(seurat_object, latent, run_umap = TRUE, umap_n_neighbors = 15, umap_min_dist = 0.1) {
  stopifnot(ncol(seurat_object) == nrow(latent))
  rownames(latent) <- colnames(seurat_object)

  latent_dr <- CreateDimReducObject(key = "latent_", embeddings = as.matrix(latent))
  seurat_object[["latent"]] <- latent_dr

  if (run_umap) {
    latent_umap <- uwot::umap(latent, n_neighbors = umap_n_neighbors, min_dist = umap_min_dist) %>%
      as.data.frame() %>% `colnames<-`(c("UMAP1", "UMAP2"))
    rownames(latent_umap) <- colnames(seurat_object)
    seurat_object[["latent_umap"]] <- CreateDimReducObject(key = "latent_umap_", embeddings = as.matrix(latent_umap))
  }

  seurat_object
}

get_latent_umap <- function(seurat_object, n_neighbors = 15, min_dist = 0.1) {
  emb <- seurat_object[["latent"]]@cell.embeddings
  umap_df <- uwot::umap(emb, n_neighbors = n_neighbors, min_dist = min_dist)
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  seurat_object[["latent_umap"]] <- CreateDimReducObject(key = "latent_umap_", embeddings = as.matrix(umap_df))
  seurat_object
}

get_latent_clustering <- function(seurat_object,
                                  resolutions = c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3),
                                  verbose = FALSE) {
  get_clustering(seurat_object, reduction = "latent", resolutions = resolutions, verbose = verbose)
}

# ----------------------------- Seurat repair ---------------------------------

fix_seurat <- function(seurat_object) {
  md         <- seurat_object@meta.data
  counts     <- seurat_object@assays$RNA@counts
  scaled     <- seurat_object@assays$RNA@scale.data
  latent     <- seurat_object[["latent"]]
  latent_um  <- seurat_object[["latent_umap"]]

  rownames(md) <- md$barcode %||% colnames(seurat_object) # fallback if barcode absent
  old_idents <- Idents(seurat_object)

  new_obj <- CreateSeuratObject(counts = counts)
  new_obj@meta.data <- md
  new_obj@assays$RNA@counts     <- counts
  new_obj@assays$RNA@scale.data <- scaled

  if (!is.null(latent))    new_obj[["latent"]]      <- latent
  if (!is.null(latent_um)) new_obj[["latent_umap"]] <- latent_um
  Idents(new_obj) <- old_idents
  new_obj
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------- SCVI I/O --------------------------------------

write_scvi_inputs <- function(seurat_object,
                              out_dir,
                              split_by = "orig.ident",
                              exclude_genes = character()) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  Idents(seurat_object) <- seurat_object[[split_by]][,1]

  splits <- levels(Idents(seurat_object))
  invisible(lapply(splits, function(idv) {
    so <- subset(seurat_object, idents = idv)
    m  <- as.data.frame(so@assays$RNA@counts)
    if (length(exclude_genes)) {
      m <- m[!rownames(m) %in% exclude_genes, , drop = FALSE]
    }
    data.table::fwrite(m, file.path(out_dir, paste0(idv, ".csv")), sep = ",", quote = FALSE)
  }))
}

# ----------------------------- SingleR wrappers ------------------------------

run_singler <- function(seurat_object,
                        cluster = NULL,
                        method = c("cell", "cluster"),
                        sample_n = NULL,
                        refs = c("hpca", "blueprint")) {

  method <- match.arg(method)
  use_hpca <- "hpca" %in% refs
  use_blue <- "blueprint" %in% refs

  if (!is.null(sample_n)) {
    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(seq_len(ncol(seurat_object)), sample_n)])
  }

  sce <- as.SingleCellExperiment(seurat_object)
  hpca.se   <- if (use_hpca) SingleR::HumanPrimaryCellAtlasData() else NULL
  blueprint <- if (use_blue) SingleR::BlueprintEncodeData() else NULL

  if (method == "cell") {
    res <- list()
    if (use_hpca) res$hpca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.fine)
    if (use_blue) res$blue <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)

    if (is.null(sample_n)) {
      if (use_hpca) seurat_object$singler_hpca_pred <- res$hpca$first.labels
      if (use_blue) seurat_object$singler_blueprint_pred <- res$blue$first.labels
      return(seurat_object)
    } else {
      out <- tibble(
        barcode = colnames(seurat_object),
        cluster = seurat_object$seurat_clusters %||% Idents(seurat_object) %>% as.character()
      )
      if (use_hpca) out$singler_hpca_pred <- res$hpca$labels
      if (use_blue) out$singler_blueprint_pred <- res$blue$labels
      return(out)
    }
  }

  if (method == "cluster") {
    if (is.null(cluster)) cluster <- Idents(seurat_object) %>% as.character()
    res <- list()
    if (use_hpca) res$hpca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,
                                               labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    if (use_blue) res$blue <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1,
                                               labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    out <- tibble(cluster = rownames(res[[1]]))
    if (use_hpca) out$singler_hpca_pred <- res$hpca$labels
    if (use_blue) out$singler_blueprint_pred <- res$blue$labels
    return(out)
  }
}

# ----------------------------- Enrichment ------------------------------------

# expected term object format:
# a data frame with columns: gene (character), set (term name)
# universe_df and genes_df should contain at least: gene (character)
# genes_df also needs: cluster, direction

get_hypergeometric <- function(genes_df, universe_df, term_df) {
  # for each set, compute overlap and phyper p-value + ratios
  term_list <- split(term_df$gene, term_df$set)
  gene_univ <- intersect(universe_df$gene, unique(term_df$gene))
  gene_hits <- intersect(genes_df$gene, gene_univ)

  tibble::tibble(
    ID = names(term_list),
    Description = names(term_list),
    set_size = vapply(term_list, length, integer(1)),
    overlap = vapply(term_list, function(gs) length(intersect(gs, gene_hits)), integer(1))
  ) %>%
    mutate(
      universe_size = length(gene_univ),
      query_size = length(gene_hits),
      pvalue = phyper(q = overlap - 1, m = set_size, n = universe_size - set_size, k = query_size, lower.tail = FALSE),
      GeneRatio = paste0(overlap, "/", query_size),
      BgRatio   = paste0(set_size, "/", universe_size)
    ) %>%
    arrange(pvalue)
}

get_enrichment_term <- function(genes_df, universe_df, term_df, term_name,
                                directions = c("up", "down"),
                                clusters = NULL) {
  if (is.null(clusters)) clusters <- unique(genes_df$cluster)
  res <- list(); i <- 1
  for (cl in clusters) {
    for (dirn in directions) {
      degs <- genes_df %>% filter(.data$cluster == cl, .data$direction == dirn)
      if (nrow(degs) == 0) next
      tmp <- get_hypergeometric(degs, universe_df, term_df)
      if (nrow(tmp) == 0) next
      res[[i]] <- tmp %>% mutate(cluster = cl, direction = dirn, term = term_name)
      i <- i + 1
    }
  }
  if (!length(res)) return(tibble())
  bind_rows(res)
}

get_enrichment_full <- function(genes_df,
                                universe_df,
                                hallmark = NULL,
                                tf = NULL,
                                reactome = NULL,
                                kegg = NULL,
                                directions = c("up", "down"),
                                clusters = NULL) {
  pieces <- list()
  if (!is.null(hallmark)) pieces$hallmark <- get_enrichment_term(genes_df, universe_df, hallmark, "hallmark", directions, clusters)
  if (!is.null(tf))       pieces$tf       <- get_enrichment_term(genes_df, universe_df, tf,       "tf",       directions, clusters)
  if (!is.null(reactome)) pieces$reactome <- get_enrichment_term(genes_df, universe_df, reactome, "reactome", directions, clusters)
  if (!is.null(kegg))     pieces$kegg     <- get_enrichment_term(genes_df, universe_df, kegg,     "kegg",     directions, clusters)

  df <- bind_rows(pieces)
  if (!nrow(df)) return(df)

  # compute odds ratio from ratios robustly
  split_ratio <- function(x) {
    as.numeric(strsplit(x, "/")[[1]])
  }
  or <- pmap_dbl(df[,c("GeneRatio","BgRatio")], function(GeneRatio, BgRatio) {
    a_b <- split_ratio(GeneRatio); c_d <- split_ratio(BgRatio)
    a <- a_b[1]; b <- a_b[2]
    c <- c_d[1]; d <- c_d[2]
    # guard against zeros
    if (any(c(b, d) == 0)) return(NA_real_)
    (a / b) / (c / d)
  })
  df$or <- or
  df
}
