

# Autoencoder -------------------------------------------------------------

activation_fns <- c("relu", "sigmoid", "softmax", "softplus", "softsign", "tanh", "selu", "elu", "exponential")

# Version  ----------------------------------------------------------------

current_spata_version <- list(major = 1, minor = 1, patch = 0)


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


# Dimensional reduction ---------------------------------------------------

dim_red_methods <- c("pca", "umap", "tsne")


# Gene sets ---------------------------------------------------------------

gene_set_methods <- c("mean", "gsva", "ssgsea", "zscore", "plage")


# De analysis -------------------------------------------------------------

de_methods <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")

de_df_columns <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene")



