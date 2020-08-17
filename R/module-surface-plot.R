#' @title UI of surface plot module
#'
#' @param id The namespace id.
#'

moduleSurfacePlotUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::column(width = 12,
                  shiny::wellPanel(
                    shiny::fluidRow(
                      shiny::column(width = 4,
                                    shiny::uiOutput(ns("sample_opts")),
                                    shiny::uiOutput(ns("aes_clr_opts")),
                                    shiny::uiOutput(ns("aes_clr_opts_detailed")),
                                    shiny::conditionalPanel(
                                      condition = "input.aes_clr_opts == 'gene_set'", ns = ns,
                                      shiny::selectInput(ns("method_gs"),
                                                         label = "Gene set method",
                                                         choices = c("Mean" = "mean",
                                                                     "Gene Set Variation Analysis" = "gsva",
                                                                     "Gene Set Enrichment Analysis" = "ssgsea",
                                                                     "Z-Score" = "zscore",
                                                                     "Plage" = "plage" ))
                                    ),
                                    shiny::selectInput(ns("clrsp"),
                                                       label = "Color spectrum",
                                                       choices = c("Magma" = "magma",
                                                                   "Inferno" = "inferno",
                                                                   "Cividis" = "cividis",
                                                                   "Plasma" = "plasma",
                                                                   "Viridis" = "viridis")),
                                    shiny::HTML("<br>")
                      ),
                      shiny::column(width = 8,
                                    shiny::plotOutput(ns("surface_plot"), dblclick = ns("surface_plot_dblclick"))
                      )

                    ),
                    shiny::HTML("<br>"),
                    shiny::fluidRow(
                      shiny::column(width = 4, align = "center",
                                    shiny::actionButton(ns("update_plot"), label = "Plot & Update"),
                                    shinyWidgets::dropdownButton(
                                      shiny::sliderInput(ns("pt_size"), label = "Size of points", min = 1, max = 5, step = 0.01, value = 2.75),
                                      shiny::sliderInput(ns("pt_alpha"), label = "Transparency of points", min = 0.01, max = 0.99, step = 0.01, value = 0.1),
                                      shinyWidgets::materialSwitch(ns("perform_smoothing"), value = FALSE, label = "Smooth values:", status = "success"),
                                      shiny::conditionalPanel(condition = "input.perform_smoothing == 1", ns = ns,
                                                              shiny::sliderInput(ns("span_smoothing"),
                                                                                 label = "Smoothing Degree",
                                                                                 min = 0.01, max = 0.1, step = 0.001, value = 0.015)
                                      ),
                                      shinyWidgets::materialSwitch(ns("perform_normalization"), value = TRUE, label = "Normalize values:", status = "success"),
                                      label = NULL,
                                      up = T,
                                      circle = F,
                                      icon = shiny::icon("gear"),
                                      inline = T
                                    )
                      ),
                      shiny::column(width = 8, align = "center",
                                    shinyWidgets::checkboxGroupButtons(
                                      inputId = ns("display_add_ons"),
                                      label = NULL,
                                      choices = c("Legend" = "legend", "Image" = "image", "Title" = "title", "Coordinates" = "coords", "Segmentation" = "segmentation"),
                                      direction = "horizontal", justified = F, individual = F)
                      )
                    )
                  )
    )
  )




}


#' @title Server of surface plot module
#'
#' @param id  The namespace id.
#' @param object A valid spata-object.
#' @param final_plot The final plot that is to be displayed. (See details.).
#' @param reactive_object A valid (reactive) spata-object.
#'
#' @return A reactive list with several slots:
#'  \enumerate{
#'   \item $assembled_plot() The surface plot as a ggplot-object.
#'   \item $dblclick() A list containing information regarding the double clicked position in the plot.
#'   \item $current_setting() A list with information about the settings of \code{assembled_plot} (e.g. sample, color_to, smooth, smoothing_span ...)}
#'
#' @details The argument \code{final_plot} takes a ggplot object as input which is going to be displayed as the final plot. This allows to
#' adjust the output of \code{$assembled_plot()} outside of the module. If no further adjustment is needed determine \code{final_plot} as:
#' \code{shiny::reactive(*module_return_variable*()$assembled_plot())}
#'

moduleSurfacePlotServer <- function(id, object, final_plot, reactive_object){

  shiny::moduleServer(
    id = id,
    module = function(input,
                      output,
                      session){

      # Reactive values -----------------------------------------------------------

      return_plot <- shiny::reactiveVal(list())

      current <- shiny::reactiveValues(

        sample = samples(object)[1],
        color_code = "gene_set",
        gene_set = character(),
        method_gs = character(),
        genes = character(),
        feature = character(),
        clrsp = character(),
        smooth = logical(),
        span = numeric()

      )

      reset_select_gene_sets <- shiny::reactiveVal(value = 0)
      reset_select_genes <- shiny::reactiveVal(value = 0)

      all_features <- getFeatureNames(object) %>% base::unname()
      all_gene_sets <- getGeneSets(object = object)
      all_genes <- getGenes(object = object)

      # Render UIs and Outputs --------------------------------------------------

      # select outputs
      output$sample_opts <- shiny::renderUI({

        ns <- session$ns

        shiny::selectInput(ns("sample_opts"),
                           label = "Choose sample:",
                           choices = samples(object),
                           selected = samples(object)[1])

      })

      output$aes_clr_opts <- shiny::renderUI({

        ns <- session$ns

        shiny::selectInput(inputId = ns("aes_clr_opts"),
                           label = "Choose colour code:",
                           choices = c("Gene set" = "gene_set", "Genes" = "genes", "Feature" = "feature"),
                           selected = "features")

      })

      select_gene_sets <- shiny::eventReactive(reset_select_gene_sets(),{

        ns <- session$ns

        shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                  label = "Choose gene set:",
                                  choices = all_gene_sets,
                                  selected = all_gene_sets[1],
                                  options = list(`live-search` = TRUE),
                                  multiple = F)

      })
      select_genes <- shiny::eventReactive(reset_select_genes(),{

        ns <- session$ns

        shiny::tagList(
          shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                    label = "Choose gene(s):",
                                    choices = all_genes,
                                    #selected = all_genes()[1],
                                    options = shinyWidgets::pickerOptions(
                                      liveSearch = TRUE,
                                      actionsBox = TRUE),
                                    multiple = TRUE),
          shiny::checkboxInput(ns("reset_select_genes"),
                               label = "Automatic reset",
                               value = TRUE))

      })
      select_features <- shiny::reactive({

        ns <- session$ns

        shiny::selectInput(inputId = ns("aes_clr_opts_detailed"),
                           label = "Choose feature:",
                           choices = all_features,
                           #selected = all_features()[1],
                           multiple = F
        )

      })

      output$aes_clr_opts_detailed <- shiny::renderUI({

        if(input$aes_clr_opts == "gene_set"){

          return(select_gene_sets())

        } else if(input$aes_clr_opts == "genes"){

          return(select_genes())

        } else if(input$aes_clr_opts == "feature"){

          return(select_features())

        }

      })


      # plot output
      output$surface_plot <- shiny::renderPlot({

        shiny::req(final_plot())

        final_plot()

      })


      # Plot add-ons ------------------------------------------------------------

      #----- Image add-on -----#

      image_add_on <- shiny::reactive({

        ## set up background
        if("image" %in% input$display_add_ons){

          ## extract image info
          img_info <-
            image(object, of_sample = current$sample) %>%
            grDevices::as.raster() %>%
            magick::image_read() %>%
            magick::image_info()

          st_image <-
            grDevices::as.raster(image(object, of_sample = current$sample)) %>%
            magick::image_read() %>% magick::image_flip()

          image_add_on <-
            ggplot2::annotation_raster(raster = st_image,
                                       xmin = 0, ymin = 0,
                                       xmax = img_info$width,
                                       ymax = img_info$height)


        } else {

          image_add_on <- NULL

        }


      })


      #----- Geom point add-on -----#

      # sample coordinates
      sample_coords <- shiny::reactive({

        sample_coords <-
          coordsSpatial(object = object, of_sample = current$sample)

        print("sample_coords")
        return(sample_coords)

      })

      # rna_assay
      rna_assay <- shiny::reactive({

        rna_assay <-
          exprMtr(object = object, of_sample = current$sample)

        print("update rna_assay")
        return(rna_assay)

      })

      # gene_vls
      gene_vls <- shiny::reactive({

        genes <- current$genes

        # compute mean if neccessary
        if(length(genes) > 1){
          rna_assay <- base::colMeans(rna_assay()[genes,])
        } else {
          rna_assay <- rna_assay()[genes,]
        }


        # convert to data frame
        gene_vls <-
          hlpr_normalize_vctr(rna_assay) %>%
          as.data.frame() %>%
          magrittr::set_colnames(value = "expr_score") %>%
          tibble::rownames_to_column(var = "barcodes")

        print("update gene vls")
        return(gene_vls)

      })

      # geneset_vls
      geneset_vls <- shiny::reactive({

        shiny::req(current$gene_set)

        gene_set_df <- object@used_genesets

        genes <-
          gene_set_df %>%
          dplyr::filter(ont == current$gene_set) %>%
          dplyr::filter(gene %in% base::rownames(rna_assay())) %>%
          dplyr::pull(gene)


        if(current$method_gs == "mean"){

          geneset_vls <-
            base::colMeans(rna_assay()[genes, ]) %>%
            hlpr_normalize_vctr() %>%
            base::as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")

        } else if(current$method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {

          shiny::showNotification(
            ui = stringr::str_c("Calculating gene set score according to method: '", current$method_gs, "'. This might take a few moments.", sep = ""),
            type = "message")

          geneset_vls <-
            GSVA::gsva(expr = rna_assay()[genes,], gset.idx.list = gene_set_df, mx.diff = 1, parallel.sz = 2, method = current$method_gs, verbose = F) %>%
            t() %>%
            hlpr_normalize_vctr() %>%
            as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")

        }

        print("update geneset_vls")
        return(geneset_vls)


      })

      # fdata
      fdata <- shiny::reactive({

        fdata <-
          featureData(object = object, of_sample = current$sample)[, c("barcodes", current$feature)]

        return(fdata)

      })

      # joined data.frame
      joined_df <- shiny::reactive({

        if(current$color_code == "genes"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = gene_vls(), by = "barcodes")

        } else if(current$color_code == "gene_set"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = geneset_vls(), by = "barcodes")

        } else if(current$color_code == "feature"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = fdata(), by = c("barcodes"))

        }

        print("update joined df")
        return(joined_df)

      })

      # variable
      variable <- shiny::reactive({

        if(current$color_code %in% c("genes", "gene_set")){

          variable <- "expr_score"

        } else if(current$color_code == "feature") {

          variable <- current$feature

        }

        return(variable)

      })

      # smoothed_df
      smoothed_df <- shiny::reactive({

        if(base::isTRUE(current$smooth)){

          smoothed_df <-
            hlpr_smooth_shiny(coords_df = joined_df(),
                              variable = variable(),
                              smooth_span = current$span)

          return(smoothed_df)

        } else {

          return(joined_df())

        }

      })

      #normalized_df
      normalized_df <- shiny::reactive({

        if(base::isTRUE(current$normalize) & current$color_code != "feature"){
          print("normalizing")
          normalized_df <-
            purrr::imap_dfr(.x = smoothed_df(),
                            .f = hlpr_normalize_imap,
                            aspect = "Gene",
                            verbose = TRUE,
                            subset = variable()
            )

          return(normalized_df)


        } else {

          return(smoothed_df())

        }


      })

      # geom_point_add_on
      geom_point_add_on <- shiny::reactive({

        color <- dplyr::pull(.data = normalized_df(), variable())

        add_on <-
          list(
            ggplot2::geom_point(data = normalized_df(),
                                mapping = ggplot2::aes(x = x, y = y, color = color),
                                size = input$pt_size,
                                alpha = (1-input$pt_alpha))
          )

        return(add_on)

      })


      #----- Scale color add-on -----#
      scale_color_add_on <- shiny::reactive({

        if(current$color_code %in% c("genes", "gene_set")){

          add_on <-
            ggplot2::scale_color_viridis_c(option = current$clrsp)

        } else {

          feature <- featureData(object, of_sample = current$sample)[, current$feature]

          if(is.numeric(feature)){

            add_on <-
              ggplot2::scale_color_viridis_c(option = current$clrsp)

          } else {

            add_on <-
              ggplot2::scale_color_manual(values = c("#C4432A", "#2C6CA3", "#478C3D", "#F7E70A", "#FFA500", "#56D9ED", "#C934BD",
                                                     "#3A389C", "#64DB74", "#C9B972", "#4F1211", "#CD4F39", "#00868B", "#8B7355",
                                                     "#CAFF70", "#525252","#FFD700", "#1C86EE", "#EEAEEE", "#8B2252"))

          }

        }

        print("update colour add on")
        return(add_on)

      })


      #----- Theme add-ons -----#
      coords_add_on <- shiny::reactive({

        if("coords" %in% input$display_add_ons){

          add_on <-
            list(ggplot2::theme_bw(),
                 ggplot2::theme(
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()
                 ))

        } else {

          add_on <-
            list(ggplot2::theme_void())

        }

        return(add_on)

      })

      legend_add_on <- shiny::reactive({

        if("legend" %in% input$display_add_ons){

          if(current$color_code %in% c("gene_set", "genes")){

            legend_title = "Expr.\nscore"

          } else {

            legend_title = current$feature

          }

          add_on <-
            list(ggplot2::labs(color = legend_title))

        } else {

          add_on <-
            list(ggplot2::theme(legend.position = "none"))

        }

        print("update legend add on ")
        return(add_on)


      })

      title_add_on <- shiny::reactive({

        if("title" %in% input$display_add_ons){

          if(current$color_code == "genes"){

            genes <- current$genes

            if(length(genes) > 5){

              genes <- c(genes[1:5], stringr::str_c("... +", (length(genes)-5), sep = " "))

            }

            genes_string <- stringr::str_c(genes, collapse = ", ")

            plot_title <- stringr::str_c("Genes:", genes_string, sep = " ")

          } else if(current$color_code == "gene_set"){

            gene_set <- current$gene_set

            gene_set_string <- stringr::str_c(gene_set, " (", current$method_gs, ")", sep = "")

            plot_title <- stringr::str_c("Gene set:", gene_set_string, sep = " ")

          } else {

            plot_title <- stringr::str_c("Feature:", current$feature, sep = " ")

          }

          add_on <- ggplot2::labs(title = plot_title)

        } else {

          add_on <- NULL

        }

        return(add_on)


      })

      segmentation_add_on <- reactive({

        if("segmentation" %in% input$display_add_ons){

          if(nrow(segmentation_df()) == 0){

            showNotification(ui = stringr::str_c("Sample", current$sample, "has not been segmented so far.", sep = " "))
            return(list())

          } else {

            segm_data <-
              coordinates(object = object, of_sample = current$sample) %>%
              dplyr::left_join(y = segmentation_df(), by = "barcodes")

            segm_layer <-
              list(
                ggplot2::geom_point(data = segm_data,
                                    mapping = ggplot2::aes(x = x, y = y, fill = segment),
                                    size = shiny::isolate(input$pt_size),
                                    shape = 21,
                                    color = ggplot2::alpha("white", 0.01)
                ),
                ggsci::scale_fill_igv()
              )

          }

        } else {

          return(list())

        }

      })

      segmentation_df <- reactive({

        segm_df <-
          featureData(object = reactive_object()) %>%
          dplyr::filter(segment != "") %>%
          dplyr::select(barcodes, segment)

        return(segm_df)

      })

      print(3)


      # Assembled plot ----------------------------------------------------------

      assembled_plot <- shiny::reactive({

        shiny::req(input$update_plot)

        ggplot2::ggplot() +
          image_add_on() +
          geom_point_add_on() +
          scale_color_add_on() +
          segmentation_add_on() +
          coords_add_on() +
          legend_add_on() +
          title_add_on() +
          ggplot2::coord_equal()

      })


      # Observe events ----------------------------------------------------------

      # update plot by updating reactive values
      oe <- shiny::observeEvent(input$update_plot, {

        current$sample = input$sample_opts
        current$color_code = input$aes_clr_opts

        if(current$color_code == "genes"){

          current$genes = input$aes_clr_opts_detailed

        } else if(current$color_code == "gene_set"){

          current$gene_set = input$aes_clr_opts_detailed
          current$method_gs = input$method_gs

        } else if(current$color_code == "feature"){

          current$feature = input$aes_clr_opts_detailed

        }

        current$size = input$pt_size
        current$clrsp = input$clrsp
        current$smooth = input$perform_smoothing
        current$span = input$span_smoothing
        current$normalize = input$perform_normalization

        if(base::isTRUE(input$reset_select_genes) &&
           current$color_code == "genes"){
          reset_select_genes((reset_select_genes() + 1))
        }

      })


      # Return values -----------------------------------------------------------

      return_list <- shiny::reactive({

        list(
          assembled_plot = shiny::reactive({assembled_plot()}),
          dblclick = shiny::reactive({input$surface_plot_dblclick}),
          current_setting = shiny::reactive({current}),
          smoothed_df = shiny::reactive({smoothed_df()}),
          variable = shiny::reactive({variable()})
        )

      })


      base::return(return_list)

    })


}
