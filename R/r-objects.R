

seurat_process_fns <- c("SCTransform","NormalizeData", "FindVariableFeatures", "ScaleData",
                        "RunPCA", "FindNeighbors", "FindClusters", "RunTSNE", "RunUMAP" )


seurat_methods <- c("spatial", "single_cell")

seurat_coords_from_opts <- c("umap", "tsne")
