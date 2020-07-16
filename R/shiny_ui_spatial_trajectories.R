##----- ui part of spatial segmentation mini app

shiny_ui_spatial_trajectories <- shiny::fluidPage(

    ##----- title
    shiny::titlePanel(title = "Spatial Trajectories"),

    shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0, 10), color = "red"),

    ##----- main panel
    shiny::mainPanel(
      shiny::fluidRow(
        shiny::column(width = 4,
                      shiny::wellPanel(
                        shiny::tags$h3(shiny::strong("Instructions")),
                        shiny::HTML("<br>"),
                        shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                        shiny::HTML("<br>"),
                        shiny::helpText("2. Determine the vertices of the trajectory by 'double - clicking' the position on the plot."),
                        shiny::HTML("<br>"),
                        shiny::helpText("3. Highlight or reset the trajectory by clicking the respective button below."),
                        shiny::sliderInput("trajectory_width", label = "Determine width of trajectory", value = 20, min = 5, max = 100, step = 1),
                        shiny::HTML("<br>"),
                        shiny::splitLayout(
                          shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                          shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                          cellWidths = c("50%", "50%")
                        ),
                        shiny::HTML("<br>"),
                        shiny::helpText("4. Enter the name you want to give the trajectory as well as a 'guiding comment' and click the 'Save'-button."),
                        shiny::splitLayout(
                          shiny::actionButton("save_trajectory", "Save Segment", width = "100%"),
                          shiny::textInput("name_trajectory", label = NULL, placeholder = "Name trajectory", value = ""),
                          cellWidths = c("50%", "50%")
                        ),
                        shiny::textInput("comment_trajectory", label = NULL, placeholder = "A guiding comment.", value = ""),
                        shiny::HTML("<br>"),
                        shiny::helpText("5. If you are done click on 'Close application'."),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(width = 12, align = "center",
                                        shiny::actionButton("close_app", label = "Close application", width = "50%")
                          )
                        )
                      )),
        shiny::column(width = 8,
                      shiny::wellPanel(
                        shiny::tags$h3(shiny::strong("Surface Plot")),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(width = 4,
                                        shiny::uiOutput("sample_opts"),
                                        shiny::uiOutput("aes_clr_opts"),
                                        shiny::uiOutput("aes_clr_opts_detailed"),
                                        shiny::conditionalPanel(
                                          condition = "input.aes_clr_opts == 'gene_set'",
                                          shiny::selectInput("method_gs",
                                                             label = "'Gene set method'",
                                                             choices = c("Mean" = "mean",
                                                                         "Gene Set Variation Analysis" = "gsva",
                                                                         "Gene Set Enrichment Analysis" = "ssgsea",
                                                                         "Z-Score" = "zscore",
                                                                         "Plage" = "plage" ))
                                          ),
                                        shiny::selectInput("clrsp",
                                                           label = "Color spectrum",
                                                           choices = c("Magma" = "magma",
                                                                       "Inferno" = "inferno",
                                                                       "Cividis" = "cividis",
                                                                       "Plasma" = "plasma",
                                                                       "Viridis" = "viridis")),
                                        shiny::HTML("<br>")
                                        ),
                          shiny::column(width = 8,
                                        shiny::plotOutput("surface_plot", dblclick = "traj_vertices")
                                        )

                       ),
                       shiny::HTML("<br>"),
                       shiny::fluidRow(
                         shiny::column(width = 4, align = "center",
                         shiny::actionButton("update_plot", label = "Plot & Update", width = "50%"),
                         shinyWidgets::dropdownButton(
                           shiny::sliderInput("pt_size", label = "Size of points", min = 1, max = 5, step = 0.01, value = 2.75),
                           shiny::sliderInput("pt_alpha", label = "Transparency of points", min = 0.01, max = 0.99, step = 0.01, value = 0.1),
                           shinyWidgets::materialSwitch("perform_smoothing", value = F, label = "Smooth values:", status = "success"),
                           shiny::conditionalPanel(condition = "input.perform_smoothing == 1",
                               shiny::sliderInput("span_smoothing",
                                                  label = "Smoothing Degree",
                                                  min = 0.01, max = 0.1, step = 0.001, value = 0.015)
                           ),
                           label = NULL,
                           up = T,
                           circle = F,
                           icon = shiny::icon("gear"),
                           inline = T
                         )
                       ),
                       shiny::column(width = 8, align = "center",
                                     shinyWidgets::checkboxGroupButtons(
                                       inputId = "display_add_ons",
                                       label = NULL,
                                       choices = c("Legend" = "legend", "Image" = "image", "Title" = "title", "Coordinates" = "coords"),
                                       direction = "horizontal", justified = F, individual = F)
                                     )
                      )
                    )
                )

      )
    )

  )
