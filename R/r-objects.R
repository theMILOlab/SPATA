
# Autoencoder -------------------------------------------------------------

activation_fns <- c("relu", "sigmoid", "softmax", "softplus", "softsign", "tanh", "selu", "elu", "exponential")


# Classes -----------------------------------------------------------------

numeric_classes <- c("numeric", "integer", "double")


# Clustering --------------------------------------------------------------

hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

# De analysis -------------------------------------------------------------

de_methods <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")

de_df_columns <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene")


# Data structure ----------------------------------------------------------

spata_df_vars <- c("barcodes", "sample")

coords_df_vars <- c(spata_df_vars, "x", "y")


# Dimensional reduction ---------------------------------------------------

dim_red_methods <- c("pca", "umap", "tsne")

# Gene sets ---------------------------------------------------------------

gene_set_methods <- c("mean", "gsva", "ssgsea", "zscore", "plage")


# Hotspot analysis --------------------------------------------------------

center_coords <- c("center_x", "center_y")

hotspot_list_slots <-
  c("df", "mtr_name", "sample", "smooth_span", "suggestion",
    "threshold_qntl", "threshold_stpv", "threshold_stw")

hotspot_plot_types <- c("encircle", "expression", "density_2d", "density_filled")

pr_list_slots <-
  list("hotspot" = hotspot_list_slots)

pr_methods <- c("hotspot")



# Information -------------------------------------------------------------

directory_options <- c("cell_data_set", "seurat_object", "spata_object")

default_colors <- viridis::viridis(15)

default_instructions_object <-
  methods::new(Class = "default_instructions",
               average_genes = FALSE,
               clrp = "milo",
               clrsp = "inferno",
               colors = default_colors,
               complete = TRUE,
               display_image = FALSE,
               display_labels = TRUE,
               display_points = FALSE,
               display_residuals = TRUE,
               display_title = TRUE,
               display_trajectory_parts = TRUE,
               max_adj_pval = 0.5,
               method_de = "wilcox",
               method_dr = "umap",
               method_gs = "mean",
               method_hclust = "complete",
               method_ovl = "classic",
               method_padj = "fdr",
               normalize = TRUE,
               n_highest_lfc = 50,
               n_lowest_pval = 50,
               n_pcs = 10,
               position = "fill",
               pt_alpha = 0.9,
               pt_clr = "lightgrey",
               pt_clrp = "milo",
               pt_clrsp = "inferno",
               pt_size = 2,
               show_rownames = FALSE,
               show_colnames = FALSE,
               smooth = TRUE,
               smooth_span = 0.02,
               uniform_genes = "discard",
               verbose = TRUE)


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
