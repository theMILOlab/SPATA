##----- server part of spatial segmentation mini app

shiny_server_spatial_trajectories <- function(input, output, session){

# Reactive values ---------------------------------------------------------

    spata_obj <- shiny::reactiveVal(value = object)
    plot_obj <- shiny::reactiveVal(value = object)

    vertices_df <-
      shiny::reactiveVal(value = data.frame(x = numeric(0),
                                            y = numeric(0)))

    segment_trajectory_df <- shiny::reactiveVal(
      value = data.frame(x = numeric(0),
                         y = numeric(0),
                         xend = numeric(0),
                         yend = numeric(0),
                         part = character(0))
    )

    compiled_trajectory_df <- shiny::reactiveVal(
      value = data.frame(barcodes = character(0),
                         sample = character(0),
                         x = numeric(0),
                         y = numeric(0),
                         projection_length = numeric(0),
                         trajectory_part = character(0),
                         stringsAsFactors = F)
    )

    current <- reactiveValues(

      sample = samples(object)[1],
      color_code = character(),
      gene_set = character(),
      method_gs = character(),
      genes = character(),
      feature = character(),
      clrsp = character(),
      smooth = logical(),
      span = numeric()

    )


# Reactive expressions ----------------------------------------------------

    ### reactive plot hierarchy

    ## 1. add ons to assemble the plot

    # image add on  ----------

    image_add_on <- shiny::reactive({

      ## set up background
      if("image" %in% input$display_add_ons){

        ## extract image info
        image <-
          image(object = plot_obj(), of_sample = current$sample) %>%
          grDevices::as.raster() %>%
          magick::image_read()

        image_info <- magick::image_info(image = image)

        image_flipped <- magick::image_flip(image = image)

        image_add_on <-
          ggplot2::annotation_raster(raster = image_flipped,
                                     xmin = 0, ymin = 0,
                                     xmax = image_info$width,
                                     ymax = image_info$height)


      } else {

        image_add_on <- NULL

      }


    })

    # geom point add on  ----------

    # sample coordinates
    sample_coords <- shiny::reactive({

      sample_coords <-
        coordinates(object = spata_obj(), of_sample = current$sample)

      print("sample_coords")
      return(sample_coords)

    })

    # rna_assay
    rna_assay <- shiny::reactive({

      rna_assay <-
        exprMtr(object = spata_obj(), of_sample = current$sample)

      print("update rna_assay")
      return(rna_assay)

    })

    # gene_vls
    gene_vls <- shiny::reactive({

      genes <- current$genes

      ## compute mean if neccessary
      if(base::length(genes) > 1){
        rna_assay <- base::colMeans(rna_assay()[genes,])
      } else {
        rna_assay <- rna_assay()[genes,]
      }


      ## convert to data frame
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

      gene_set_df <- spata_obj()@used_genesets

      genes <-
        gene_set_df %>%
        dplyr::filter(ont == current$gene_set) %>%
        dplyr::filter(gene %in% rownames(rna_assay())) %>%
        dplyr::pull(gene)


      if(current$method_gs == "mean"){

        geneset_vls <-
          base::colMeans(rna_assay()[genes, ]) %>%
          hlpr_normalize_vctr() %>%
          as.data.frame() %>%
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
        featureData(object = spata_obj(), of_sample = current$sample)[, c("barcodes", current$feature)]

      print("update featuer data")
      print(head(fdata))
      return(fdata)

    })

    # join coordinates with specific data.frame
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


    # geom_point_add_on
    geom_point_add_on <- shiny::reactive({

      color <- dplyr::pull(.data = smoothed_df(), variable())

      add_on <-
        list(
          ggplot2::geom_point(data = smoothed_df(),
                              mapping = ggplot2::aes(x = x, y = y, color = color),
                              size = input$pt_size,
                              alpha = (1-input$pt_alpha))
        )

      return(add_on)


    })


    # scale_colour_add_on  ----------

    scale_color_add_on <- shiny::reactive({

      if(current$color_code %in% c("genes", "gene_set")){

        add_on <-
          ggplot2::scale_color_viridis_c(option = current$clrsp)

      } else {

        feature <- featureData(spata_obj(), of_sample = current$sample)[, current$feature]

        if(base::is.numeric(feature)){

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


    # theme_add_ons  ----------

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


    ## 2. additional layer
    trajectory_segment_layer <- reactive({

      new_layer <- list()

      ## update geom_point layer
      if(nrow(vertices_df()) >= 1){

        new_layer[[1]] <-
          ggplot2::geom_point(data = vertices_df(),
                              mapping = ggplot2::aes(x = x, y = y),
                              size = 3.5, color = "black")

      }

      ## update geom_segment layer
      if(base::nrow(segment_trajectory_df()) >= 1){

        new_layer[[2]] <-
          ggplot2::geom_segment(data = segment_trajectory_df(),
                                mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                size = 1.25, color = "black",
                                arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
          )

      }

      return(new_layer)

    })

    trajectory_color_layer <- reactive({

      if(!base::nrow(compiled_trajectory_df()) == 0){

        joined_traj_df <-
          dplyr::left_join(x = compiled_trajectory_df(),
                           y = dplyr::select(smoothed_df(), -x, -y),
                           by = "barcodes")

        color <- dplyr::pull(.data = joined_traj_df, variable())

        add_on_layer <-
          list(
            ggplot2::geom_point(data = smoothed_df(),
                                mapping = ggplot2::aes(x = x, y = y),
                                color = "lightgrey", size = input$pt_size),
            ggplot2::geom_point(data = joined_traj_df,
                                mapping = ggplot2::aes(x = x, y = y, color = color),
                                size = input$pt_size)

          )

      } else {

        add_on_layer <- list()

      }

      return(add_on_layer)

    })


    ## 3 assembled plot
    assembled_plot <- shiny::reactive({

      shiny::req(input$update_plot)

        ggplot2::ggplot() +
          image_add_on() +
          geom_point_add_on() +
          scale_color_add_on() +
          coords_add_on() +
          legend_add_on() +
          title_add_on() +
          ggplot2::coord_equal()+
            trajectory_color_layer() +
            trajectory_segment_layer()

    })




# Render UIs and Outputs --------------------------------------------------

    ## select outputs
    output$sample_opts <- shiny::renderUI({
      #print("##### samples opts #####")

      shiny::selectInput("sample_opts",
                         label = "Choose sample:",
                         choices = samples(plot_obj()),
                         selected = samples(plot_obj())[1])

    })

    output$aes_clr_opts <- shiny::renderUI({

      shiny::selectInput(inputId = "aes_clr_opts",
                         label = "Choose colour code:",
                         choices = c("Gene set" = "gene_set", "Genes" = "genes", "Feature" = "feature"),
                         selected = "features")

    })

    output$aes_clr_opts_detailed <- shiny::renderUI({

      print("##### color opts detailed #####")

      if(input$aes_clr_opts == "gene_set"){

        gene_sets <- getGeneSets(object = plot_obj())

        shinyWidgets::pickerInput(inputId = "aes_clr_opts_detailed",
                                  label = "Choose gene set:",
                                  choices = gene_sets,
                                  selected = gene_sets[1],
                                  options = list(`live-search` = TRUE),
                                  multiple = F
        )

      } else if(input$aes_clr_opts == "genes"){


        genes <- getGenes(object = plot_obj(), in_sample = current$sample)

        shinyWidgets::pickerInput(inputId = "aes_clr_opts_detailed",
                                  label = "Choose gene(s):",
                                  choices = genes,
                                  selected = genes[1],
                                  options = list(`live-search` = TRUE),
                                  multiple = T
        )

      } else if(input$aes_clr_opts == "feature"){

        features <- featureNames(object = plot_obj()) %>% base::unname()

        shiny::selectInput(inputId = "aes_clr_opts_detailed",
                           label = "Choose feature:",
                           choices = features,
                           selected = features[1],
                           multiple = F
        )

      }


    })


    ## plot output
    output$surface_plot <- shiny::renderPlot({

      shiny::req(assembled_plot())

      assembled_plot()

      })




# Observe events and reactive events --------------------------------------

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

      current$clrsp = input$clrsp
      current$smooth = input$perform_smoothing
      current$span = input$span_smoothing

    })

    # 1. add trajectory vertice consecutively
    oe <- shiny::observeEvent(input$traj_vertices, {

      # 1. prolong and update data.frame
      vrtcs_list <- input$traj_vertices
      new_df <- dplyr::add_row(.data = vertices_df(),
                               x = vrtcs_list$x,
                               y = vrtcs_list$y)

      vertices_df(new_df)

      # 2. update trajectory df
      n_vrt <- nrow(vertices_df())

      if(n_vrt >= 2){

        stdf <-
          segment_trajectory_df() %>%
          dplyr::add_row(
            x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
            y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
            xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
            yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
            part = stringr::str_c("part", n_vrt-1 , sep = "_")
          )

        print(stdf)
        segment_trajectory_df(stats::na.omit(stdf))

      } else {

        segment_trajectory_df(data.frame(
                                x = numeric(0),
                                y = numeric(0),
                                xend = numeric(0),
                                yend = numeric(0),
                                part = character(0),
                                stringsAsFactors = F))

      }

    })

    # 3.1
    oe <- shiny::observeEvent(input$highlight_trajectory, {

      checkpoint(evaluate = nrow(segment_trajectory_df()) >= 1, case_false = "insufficient_n_vertices2", e_n)

        print(head(segment_trajectory_df()))

        compiled_trajectory_df <-
          hlpr_compile_trajectory(segment_trajectory_df = segment_trajectory_df(),
                                  trajectory_width = input$trajectory_width,
                                  object = spata_obj(),
                                  sample = current$sample)

        print(head(compiled_trajectory_df))

        compiled_trajectory_df(compiled_trajectory_df)

    })


    # 3.2 reset current vertices
    oe <- shiny::observeEvent(input$reset_trajectory, {

      vertices_df(data.frame(x = numeric(0),
                             y = numeric(0)))

      segment_trajectory_df(data.frame(x = numeric(0),
                               y = numeric(0),
                               xend = numeric(0),
                               yend = numeric(0),
                               part = character(0),
                               stringsAsFactors = F))

      compiled_trajectory_df(data.frame(barcodes = character(0),
                                        sample = character(0),
                                        x = numeric(0),
                                        y = numeric(0),
                                        projection_length = numeric(0),
                                        trajectory_part = character(0),
                                        stringsAsFactors = F))

    })

    ##--- 4. save the highlighted trajectory
    oe <- shiny::observeEvent(input$save_trajectory, {

      ## control
      checkpoint(evaluate = base::nrow(compiled_trajectory_df()) > 0, case_false = "insufficient_n_vertices2")
      checkpoint(evaluate = shiny::isTruthy(x = input$name_trajectory), case_false = "invalid_trajectory_name")
      checkpoint(evaluate = !input$name_trajectory %in% getTrajectoryNames(spata_obj(), current$sample),
                    case_false = "occupied_trajectory_name")

        ## save trajectory
        spata_obj <- spata_obj()

        spata_obj@trajectories[[current$sample]][[input$name_trajectory]] <-
          methods::new("spatialTrajectory",
                       compiled_trajectory_df = compiled_trajectory_df(),
                       segment_trajectory_df = segment_trajectory_df(),
                       comment = input$comment_trajectory,
                       name = input$name_trajectory,
                       sample = current$sample)

        spata_obj(spata_obj)


        ## control
        check <- trajectory(spata_obj(), trajectory_name = input$name_trajectory, of_sample = current$sample)

        ## feedback and reset
        if(base::identical(check@compiled_trajectory_df, compiled_trajectory_df())){

          shiny::showNotification(ui = "Trajectory has been stored.", type = "message", duration = 7)


          vertices_df(data.frame(x = numeric(0),
                                 y = numeric(0)))

          segment_trajectory_df(data.frame(x = numeric(0),
                                           y = numeric(0),
                                           xend = numeric(0),
                                           yend = numeric(0),
                                           part = character(0),
                                           stringsAsFactors = F))

          compiled_trajectory_df(data.frame(barcodes = character(0),
                                            sample = character(0),
                                            x = numeric(0),
                                            y = numeric(0),
                                            projection_length = numeric(0),
                                            trajectory_part = character(0),
                                            stringsAsFactors = F))

        } else {

          shiny::showNotification(ui = "Could not save trajectory.")

        }

    })

    ##--- 5. close application and return spata object
    oe <- shiny::observeEvent(input$close_app, {

      stopApp(returnValue = spata_obj())

    })

  }
