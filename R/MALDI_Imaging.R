#' @title Initiate a spata-object from MALDI Imaging
#'
#' @description Creates, saves and returns an object of class spata
#' from scratch. Several samples can be stored in one object, though we recommand to stick
#' to one. (See details for more.)
#'
#' @param coordinates Data.frame containing "unique identifier" meaning a unique description of
#' each pixel/spot as well as  x and y coordinates of the position of the spot
#' @param intensity_matrix Matrix of intensities with rownames for Proteins/Metabolites/Lipids, and "unique identifier"
#' as colnames
#' @param gene_set_path Character value (or NULL). Specifies the path to a
#' .RDS-file containing a data.frame that is to be used as input for slot @@used_genesets.
#'
#'  Must have the character-variables
#'
#'    \itemize{
#'     \item{\emph{'ont'}: The gene set names.}
#'     \item{\emph{'gene'}: The belonging gene names.}
#'     }
#'
#' If set to NULL the default gene-set data.frame will used. Run \code{?gsdf} to get more information.
#'
#' @param sample_names Character vector. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'
#' @param file_name Character value. The name-suffix for the file name under which the
#' spata-object is stored. Is prefixed with \emph{'spata-obj-'} and suffixed with \emph{'.RDS'}.
#'
#' @inherit verbose params
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the spata-object is processed. For all of these arguments apply the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete spata-object. Use \code{validateSpataObject()} after
#' initiating it in order to see which slots are valid and which are not.
#'
#' Handling more than one sample:
#'
#' Several samples can be stored in one object. If so, the count-matrices will be combined to one matrix which is given to the seurat-object that is temporarily
#' initiated in order to perform the pre processing steps. Sample related unambiguity with respect to the barcode's belonging is maintained
#' by suffixing the barcode-sequences with the respective sample name specified in \code{sample_names}. The meta.data data.frame of the
#' seurat-object is joined with a variable called \emph{sample} denoting the sample-belonging of every barcode which can be used as input
#' for pre processing functions.
#'
#' For now we recommend to stick to one sample per spata-object.
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export
#'
#'
#'
#'



#Input files
#coordinates<- file with "unique identifier" and x and y coordinates
#Intensity matrix <- Intensity matrix of metabolites as rownames and "unique identifier" as colnames
#If requested an Image (png)


initiateSpataObject_MALDI<-function(coordinates,
                                    intensity_matrix,
                                    geneSets=NULL,
                                    sample_names=NULL,
                                    file_name=NULL,
                                    PCA_Comp=30,
                                    NN=50,
                                    verbose=T){

  if(base::is.null(sample_names)){sample_names="sample_1"}
  if(base::is.null(file_name)){file_name="MALDI_SPATAobj.RDS"}
  if(!methods::is(intensity_matrix, "matrix")) stop("intensity_matrix: Is not from class 'matrix' ")
  if(base::is.null(geneSets)){geneSets=SPATA::gsdf}
  string_contain=length(stringr::str_detect(rownames(intensity_matrix), "_") %in% TRUE) <= 1

  if(string_contain) {
    message(" rownames of intensity matrix contain a '_'. This is not allowed and will be changed into '-' ")
    rownames(intensity_matrix) <- str_replace(rownames(intensity_matrix), "_", "-")
  }

  #check inputs
  if (base::length(base::unique(coordinates$barcode %in% colnames(intensity_matrix)))!=1) stop("barcodes from coordinates are not similar to colnames of the intensity_matrix")

  #Add sample name to coordinates
  coordinates <- coordinates %>% dplyr::mutate(sample=sample_names) %>% dplyr::select("barcodes","sample","x","y")

  #Change Coordinates for image in SPATA

  if (verbose==T) { message("---------- Create SPATA object from MALDI  ------------- ") }

  #Crease empty SPATA object
  obj=methods::new(Class="spata")

  obj@coordinates <- coordinates
  obj@data@norm_exp <- intensity_matrix
  obj@data@counts=as(intensity_matrix,"dgCMatrix")
  obj@samples <- unique(obj@coordinates$sample)
  obj@fdata <- obj@coordinates %>% dplyr::select(barcodes, sample) %>% dplyr::mutate(segment = "")
  obj@used_genesets=geneSets
  obj@trajectories=list(sample=list())
  names(obj@trajectories)=obj@samples

  ## Run Analysis ##


  if (verbose==T) { message("---------- Run PCA Analysis 1/4 ------------- ")}
  pca=base::as.data.frame(pcaMethods::pca(t(intensity_matrix),method="svd", nPcs=PCA_Comp)@scores)
  if (verbose==T) { message("---------- Select Eigenvalues Analysis 2/4 ------------- ")}
  if (verbose==T) { message("---------- Run SNN-Cluster Analysis 3/4 ------------- ")}

  nearest <- RANN::nn2(pca, k = NN, treetype = "bd", searchtype = "priority")
  nearest$nn.idx <- nearest$nn.idx
  nearest$nn.dists <- nearest$nn.dists
  edges = reshape::melt(t(nearest$nn.idx[, 1:NN]))
  base::colnames(edges) = c("B", "A", "C")
  edges = edges[, c("A", "B", "C")]
  edges$B = edges$C
  edges$C = 1
  edges = base::unique(base::transform(edges, A = pmin(A, B), B = pmax(A, B)))
  names(edges) <- c("V1", "V2", "weight")
  edges$V1 <- base::rownames(pca)[edges$V1]
  edges$V2 <- base::rownames(pca)[edges$V2]
  g <- igraph::graph.data.frame(edges, directed = F)
  graph.out = igraph::cluster_louvain(g)
  clust.assign = base::factor(graph.out$membership, levels = base::sort(base::unique(graph.out$membership)))
  Cluster_out=base::data.frame(ID=base::rownames(pca), Cluster=clust.assign); base::rownames(Cluster_out)=Cluster_out$ID

  obj@fdata$cluster=Cluster_out[obj@fdata$barcodes, "Cluster"]

  if (verbose==T) { message("---------- Run Dimensional Reduction 4/4 ------------- ")}

  ts=Rtsne::Rtsne(pca, perplexity=30)
  obj@dim_red@TSNE <- base::data.frame(barcodes=base::rownames(pca), sample=sample_names, tsne1=ts$Y[,1], tsne2=ts$Y[,2])
  #plot(obj@dim_red@TSNE$tsne1, obj@dim_red@TSNE$tsne2)
  umap=umap::umap(pca)
  obj@dim_red@UMAP <- base::data.frame(barcodes=base::rownames(pca), sample=sample_names, umap1=umap$layout[,1], umap2=umap$layout[,2])

  if (verbose==T) { message("---------- Analysis Done... export data  ------------- ")}

  base::saveRDS(obj, file_name)

  return(obj)

}


