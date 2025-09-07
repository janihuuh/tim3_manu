# ------------------------------------------------------------------------------
# TIM-3 manuscript: end-to-end scRNA-seq preprocessing + labeling + scVI I/O
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

# If you placed the helpers into R/sc_utils.R:
# source("R/sc_utils.R")

# --------------------------- Project settings ---------------------------------
project        <- "tim3"
data_dir       <- "data/scRNAseq"
results_dir    <- file.path("results", project)
qc_dir_before  <- file.path(results_dir, "qc", "before_1")
qc_dir_after   <- file.path(results_dir, "qc", "after_1")
scvi_in_dir    <- file.path("results", "scvi", "input_files")       # python inputs
scvi_out_dir   <- file.path("results", "scvi", "results")           # python outputs
scvi_latent_fp <- file.path(scvi_out_dir, paste0(project, "_latent.csv"))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir_before, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir_after,  recursive = TRUE, showWarnings = FALSE)
dir.create(scvi_in_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(scvi_out_dir,  recursive = TRUE, showWarnings = FALSE)

# Helper to name samples from folder path (avoid custom getSeuratName)
get_seurat_name <- function(p) {
  # For 10x runs, parent dir of the matrix folder is typically the sample ID
  bn <- basename(p)
  if (bn %in% c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")) {
    return(basename(dirname(p)))
  }
  bn
}

# Discover 10x folders robustly
all_dirs <- list.dirs(data_dir, recursive = TRUE, full.names = TRUE)
tenx_dirs <- Filter(function(d) {
  any(file.exists(file.path(d, "matrix.mtx"))) ||
    any(file.exists(file.path(d, "matrix.mtx.gz"))) ||
    any(file.exists(file.path(d, "barcodes.tsv.gz")))
}, all_dirs)

stopifnot(length(tenx_dirs) > 0)

# ------------------------------ 0) Create Seurat ------------------------------
scrnaseq_objects <- lapply(tenx_dirs, function(x) {
  message(get_seurat_name(x))
  m <- Read10X(data.dir = x)
  CreateSeuratObject(counts = m, project = get_seurat_name(x), min.cells = 3, min.features = 200)
})

tim3_seurat <- scrnaseq_objects[[1]]  # start with one object; merge later if desired

# ------------------------------ 1) Basic QC -----------------------------------
tim3_seurat <- PercentageFeatureSet(tim3_seurat, pattern = "^MT-",  col.name = "percent.mt")
tim3_seurat <- PercentageFeatureSet(tim3_seurat, pattern = "^RP",   col.name = "percent.ribo")

# Optional: cell cycle gene percentage if you have 'cycle.genes' defined
if (exists("cycle.genes")) {
  tim3_seurat <- PercentageFeatureSet(tim3_seurat, features = cycle.genes, col.name = "percent.cycle")
}

tim3_seurat@meta.data$barcode <- colnames(tim3_seurat)

# Optional: your own QC plotting util, only if available
if (exists("plotQC")) {
  tim3_seurat %>% plotQC(folder = qc_dir_before)
}

# Use the cleaned utility (from sc_utils.R); fallback: keep all cells if not available
if (exists("get_qc")) {
  tim3_seurat <- get_qc(tim3_seurat)
}

# Doublets (scds via helper)
if (exists("get_doublets")) {
  tim3_seurat <- get_doublets(tim3_seurat)
  # Use your preferred threshold; 1.8 was your prior heuristic
  tim3_seurat <- subset(tim3_seurat, subset = hybrid_doublet_score < 1.8)
}

# ------------------ 2) SingleR predictions & rare-type pruning ----------------
if (exists("run_singler")) {
  tim3_seurat <- run_singler(tim3_seurat, method = "cell")
  md <- tim3_seurat@meta.data

  if ("singler_hpca_pred" %in% names(md)) {
    keep_hpca <- md %>% group_by(singler_hpca_pred) %>% summarise(n = n(), .groups = "drop") %>% filter(n >= 10) %>% pull(singler_hpca_pred)
    tim3_seurat$singler_hpca_pred <- ifelse(tim3_seurat$singler_hpca_pred %in% keep_hpca, tim3_seurat$singler_hpca_pred, "rare")
  }
  if ("singler_blueprint_pred" %in% names(md)) {
    keep_blue <- md %>% group_by(singler_blueprint_pred) %>% summarise(n = n(), .groups = "drop") %>% filter(n >= 10) %>% pull(singler_blueprint_pred)
    tim3_seurat$singler_blueprint_pred <- ifelse(tim3_seurat$singler_blueprint_pred %in% keep_blue, tim3_seurat$singler_blueprint_pred, "rare")
  }
}

# ------------------------------ 3) Seurat prep --------------------------------
# Build gene exclusion lists (TCR/BCR & unwanted)
clonality_genes <- grep("^(TRAV|TRBV|TRGV|TRDV|IGLV|IGLC|IGLL|IGKV|IGHV|IGKC|IGH|IGK)", rownames(tim3_seurat), value = TRUE)
unwanted_genes  <- grep("^(LINC|AC|AL|MT-|RP)", rownames(tim3_seurat), value = TRUE)
remove_genes    <- unique(c(clonality_genes, unwanted_genes))

if (exists("preprocess_seurat")) {
  tim3_seurat <- preprocess_seurat(tim3_seurat, cells_to_use = colnames(tim3_seurat), remove_genes = remove_genes)
} else {
  # minimal fallback
  tim3_seurat <- NormalizeData(tim3_seurat)
  tim3_seurat <- FindVariableFeatures(tim3_seurat)
  tim3_seurat <- ScaleData(tim3_seurat)
  tim3_seurat <- RunPCA(tim3_seurat, npcs = 50)
  tim3_seurat <- RunUMAP(tim3_seurat, dims = 1:20, learning.rate = 1)
}

if (exists("get_clustering")) {
  tim3_seurat <- get_clustering(tim3_seurat, reduction = "pca")
} else {
  tim3_seurat <- FindNeighbors(tim3_seurat, dims = 1:20)
  tim3_seurat <- FindClusters(tim3_seurat, resolution = c(seq(0.1,1,0.1), seq(1.2,2,0.2), 2.5, 3))
}

# ------------------------------ 4) scVI I/O -----------------------------------
# Make safe file names for orig.ident to avoid slashes
tim3_seurat$orig.ident <- gsub("\\/", "_", tim3_seurat$orig.ident)

if (exists("write_scvi_inputs")) {
  write_scvi_inputs(tim3_seurat, out_dir = scvi_in_dir, split_by = "orig.ident", exclude_genes = clonality_genes)
} else {
  # Basic fallback writer by orig.ident
  Idents(tim3_seurat) <- tim3_seurat$orig.ident
  for (idv in levels(Idents(tim3_seurat))) {
    so <- subset(tim3_seurat, idents = idv)
    m  <- as.data.frame(so@assays$RNA@counts)
    m  <- m[!rownames(m) %in% clonality_genes, , drop = FALSE]
    fwrite(m, file.path(scvi_in_dir, paste0(idv, ".csv")), sep = ",", quote = FALSE)
  }
}

# --- Run your Python scVI pipeline separately; then read latent back in ----
if (file.exists(scvi_latent_fp)) {
  latents <- data.table::fread(scvi_latent_fp) %>% as.matrix()
  if (exists("put_latents_seurat")) {
    tim3_seurat <- put_latents_seurat(tim3_seurat, latent = latents)
  }
  if (exists("get_latent_clustering")) {
    tim3_seurat <- get_latent_clustering(tim3_seurat)
  } else {
    # quick latent-based clustering if helpers not loaded
    if (!is.null(tim3_seurat[["latent"]])) {
      dims_lat <- seq_len(ncol(tim3_seurat[["latent"]]@cell.embeddings))
      tim3_seurat <- FindNeighbors(tim3_seurat, reduction = "latent", dims = dims_lat)
      tim3_seurat <- FindClusters(tim3_seurat, resolution = c(seq(0.1,1,0.1), seq(1.2,2,0.2), 2.5, 3))
    }
  }
  if (exists("fix_seurat")) {
    tim3_seurat <- fix_seurat(tim3_seurat)
  }
}

# --------------------------- 5) Decide on clustering --------------------------
if (exists("plot_clustering_overview")) {
  p <- plot_clustering_overview(tim3_seurat)
  ggsave(file.path(qc_dir_after, "scatter_clustering.png"), p, width = 5, height = 4)
}

# Choose a working resolution (edit as needed)
if ("RNA_snn_res.0.7" %in% colnames(tim3_seurat@meta.data)) {
  Idents(tim3_seurat) <- tim3_seurat$RNA_snn_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
} else if ("SNN_res.0.7" %in% colnames(tim3_seurat@meta.data)) {
  Idents(tim3_seurat) <- tim3_seurat$SNN_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
}

tim3_seurat$cluster <- Idents(tim3_seurat)

# DE
all_markers <- FindAllMarkers(tim3_seurat, test.use = "t")
fwrite(all_markers, file.path(results_dir, "de.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Save object
saveRDS(tim3_seurat, file.path(results_dir, "tim3_seurat.rds"))