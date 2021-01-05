

# Autoencoder -------------------------------------------------------------

activation_fns <- c("relu", "sigmoid", "softmax", "softplus", "softsign", "tanh", "selu", "elu", "exponential")



# Clustering --------------------------------------------------------------

hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

# De analysis -------------------------------------------------------------

de_methods <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")

de_df_columns <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene")

# Dimensional reduction ---------------------------------------------------

dim_red_methods <- c("pca", "umap", "tsne")


# Gene sets ---------------------------------------------------------------

gene_set_methods <- c("mean", "gsva", "ssgsea", "zscore", "plage")



# Hotspot analysis --------------------------------------------------------

hotspot_list_slots <-
  c("df", "max_hotspots", "mtr_name", "n_start", "sample", "smooth_span", "suggestion",
    "threshold_icld", "threshold_stw", "threshold_stpv", "tot_wss")


# Seurat analysis ---------------------------------------------------------

seurat_process_fns <- c("SCTransform","NormalizeData", "FindVariableFeatures", "ScaleData",
                        "RunPCA", "FindNeighbors", "FindClusters", "RunTSNE", "RunUMAP" )

seurat_methods <- c("spatial", "single_cell")

seurat_coords_from_opts <- c("umap", "tsne")


# Trajectory analysis -----------------------------------------------------

empty_ctdf <- data.frame(barcodes = character(0),
                         sample = character(0),
                         x = numeric(0),
                         y = numeric(0),
                         projection_length = numeric(0),
                         trajectory_part = character(0),
                         stringsAsFactors = FALSE)

empty_segment_df <- data.frame(x = numeric(0),
                               y = numeric(0),
                               xend = numeric(0),
                               yend = numeric(0),
                               part = character(0),
                               stringsAsFactors = FALSE)

trajectory_df_colnames <- c("trajectory_part", "trajectory_order", "trajectory_part_order")



# Version  ----------------------------------------------------------------

current_spata_version <- list(major = 1, minor = 1, patch = 0)
