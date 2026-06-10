# Tyler Therron, MS
# 20250213


#  ======================================. Libraries ==================
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(readr)
library(shiny)
library(shiny.fluent)
library(htmlwidgets) 
library(webshot) 
library(future)
library(DT)

#  ======================================. Libraries ==================
options(shiny.maxRequestSize = 1000 * 1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  title = "Myeloid Reporter Levels and Gene Expression in Mice",
  tags$head(
    tags$style(HTML("
      body {
        font-size: 16px; /* Lock the font size */
      }
      * {
        font-size: inherit !important;
      }
      .png-image {
        width: clamp(250px, 30vw, 540px) !important;
        height: auto !important;
        margin-left: clamp(60px, 8vw, 160px) !important;
        margin-bottom: clamp(-280px, -14vw, -100px) !important;
      }
      .jpg-image {
        width: clamp(60px, 7vw, 130px) !important;
        height: auto !important;
      }
      .colorkey-image {
        width: clamp(35px, 4vw, 75px) !important;
        height: auto !important;
      }
      .color-key-row {
        margin-top: 0; /* Adjust this value for more/less spacing */
        display: flex;
        margin-bottom: 0; /* Adjust this value for more/less spacing */
        height: 100%;
        width: 100%;
        margin-left: 20px;
      }
      .color-key-row-left-margin {
        margin-left: 1px; /* Adjust this value for more left margin space */
      }
      .flow-data-row {
        margin-top: 0px;
        display: flex;
        justify-content: left;
        align-items: left;
        text-align: left;
        height: 15px;
      }
      .flow-data-row-left-margin {
        margin-left: -5%; /* Adjust this value for more left margin space */
      }
      .flow-key-row-left-margin {
        margin-left: -1%;
      }
      .arrow-row-left-margin {
        margin-left: 150px; /* Adjust this value for more left margin space */
      }
      .arrow-row {
        display: flex;
        justify-content: left;
        align-items: left;
        text-align: left;
        margin-top: -50px; /* Adjust this value for more/less spacing */
        height: 50px;
      }
      .arrow-row-left-margin2 {
        margin-left: 150px; /* Adjust this value for more left margin space */
      }
      .arrow-row2 {
        display: flex;
        justify-content: left;
        align-items: left;
        text-align: left;
        margin-top: -110px; /* Adjust this value for more/less spacing */
        height: 50px;
      }
      .arrow-column {
        display: flex;
        justify-content: center;
        align-items: center;
        text-align: center;
      }
      .arrow-image {
        width: clamp(180px, 46vw, 700px);
        height: auto;
      }
      .arrow-image2 {
        width: clamp(180px, 46vw, 700px);
        height: auto;
      }
      .vert-arrow-image {
        width: auto;
        height: min(15vh, 200px);
        min-height: 80px;
        margin-right: -10px;
        margin-left: 30px;
      }
      .logo-img {
        max-width: 120px; /* Adjust as needed */
        float: right;
        margin-top: 10px;
        margin-bottom: -10%;
        border-radius: 20px;
      }
      @media only screen and (max-width: 1000px) {
        .col-md-3 {
          flex: 0 0 100%;
          max-width: 100%;
        }
      }
    ")),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no")  # Disable user scaling
  ),
  fluidRow(
    column(6,
      titlePanel(div(p("Myeloid Reporter Levels and Gene Expression in Mice", style = "font-size: 24px !important; margin-top: 2.5%")))
    ),
    column(6, 
           # Logo image on the right
           tags$div(
             tags$img(src = "3.png", class = "logo-img")
           )
          )
  ),
     sidebarLayout(
        sidebarPanel(
          div(p("For each myeloid reporter mouse model, 
                visualize the fluorescence intensity measured by flow 
                cytometry and the matching expression of the reporter gene by 
                single-cell RNA-seq in the joint synovium.", style = "font-size: 13px !important; font-weight: bold;")),
          br(),
          selectInput("gene_choices", div(p("Select Mouse Model to Visualize:", style = "font-size: 12px !important")), choices = NULL),
          uiOutput('select_pri'), br(),
          div(p("Please cite:", style = "font-size: 10px !important; 
                font-weight: bold !important; margin-bottom: -5%")),
          tags$div(uiOutput("citation1"), 
                   style = "font-size: 9px !important; text-align: left; padding-top: 20px;"), width = 3
        ),
        mainPanel(
          br(), br(),
          fluidRow(uiOutput("combined_images")),
          width = 9
        )
      )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # Track current tab so it persists when mouse model changes
  current_tab <- reactiveVal("Monocytes/Macrophages")

  observeEvent(input$main_tabs, {
    current_tab(input$main_tabs)
  }, ignoreInit = TRUE)

  observe({
    genes_path <- "./genes"
    genes_choice_list <- list.files(c(path = genes_path,
                                      full.names = FALSE,
                                      recursive = FALSE))
    updateSelectInput(session, "gene_choices",
                      choices = genes_choice_list)
  })
  
  #VI
  output$combined_images <- renderUI({
    req(input$gene_choices)
    
    # Directories and file processing (same as before)
    selected_gene_dir <- file.path("./genes", input$gene_choices)
    color_key_path <- file.path("./Color_Keys")
    horz_arrow_image_path <- file.path("./horz_arrow")
    vert_arrow_image_path <- file.path("./vert_arrow")
    
    # Process PNG images (unchanged)
    all_png_files <- list.files(path = selected_gene_dir, pattern = "\\.png$", full.names = TRUE)
    png_6cells <- all_png_files[grepl("_6cells\\.png$", all_png_files)]
    png_other  <- all_png_files[!grepl("_6cells\\.png$", all_png_files)]
    
    # Process JPEG images
    jpg_files <- list.files(path = selected_gene_dir, pattern = "\\.jpg$", full.names = TRUE)
    jpg_files <- jpg_files[ !grepl("(Eosinophil\\.jpg$|DC\\.jpg$|Neutrophil\\.jpg$|CD11b-\\.jpg$)", jpg_files) ]
    
    # Define desired order and reorder jpg_files using endsWith()
    order_vec <- c("Macrophage.jpg", "CM.jpg", "NCM.jpg", "TRMC.jpg",
                   "CD11b-.jpg", "CX3CR1+ MHCII-.jpg", "CX3CR1+ MHCII+.jpg",
                   "CX3CR1- MHCII-.jpg", "CX3CR1- MHCII+.jpg")
    ordered_jpg_files <- unlist(lapply(order_vec, function(pat) {
      matching <- jpg_files[sapply(basename(jpg_files), function(x) endsWith(x, pat))]
      return(matching)
    }))
    ordered_jpg_files <- unique(ordered_jpg_files)
    jpg_files <- ordered_jpg_files
    
    # Split JPEG images into top (first 5) and bottom (the rest)
    n <- length(jpg_files)
    top_jpg <- jpg_files[1:min(4, n)]
    bottom_jpg <- if(n > 4) jpg_files[5:n] else character(0)
    
    # Process color key 
    order_vec_key <- c("Macrophage.png", "CM.png", "NCM.png", "TRMC.png",
                       "CD11b-.png", "CX3CR1+ MHCII-.png", "CX3CR1+ MHCII+.png",
                       "CX3CR1- MHCII-.png", "CX3CR1- MHCII+.png")
    
    color_key_PNGs <- list.files(path = color_key_path, pattern = "\\.png$", full.names = TRUE)
    
    ordered_key_files <- unlist(lapply(order_vec_key, function(pat) {
      # Use basename to ignore directory parts and match the ending exactly
      matching <- color_key_PNGs[sapply(basename(color_key_PNGs), function(x) endsWith(x, pat))]
      return(matching)
    }))
    
    ordered_key_files <- unique(ordered_key_files)  # In case of duplicates
    
    n_keys <- length(ordered_key_files)
    top_keys <- ordered_key_files[1:min(4, n_keys)]
    bottom_keys <- if(n_keys > 4) ordered_key_files[5:n_keys] else character(0)
    
    # Aroows
    horz_arrow_images <- list.files(path = horz_arrow_image_path, pattern = "\\.png$", full.names = TRUE)
    vert_arrow_images <- list.files(path = vert_arrow_image_path, pattern = "\\.png$", full.names = TRUE)
    
    
    # Build the layout:
    # Outer fluidRow: Left column for vertical arrow; right column for JPEG images
    combined_images <- tagList(
      tabsetPanel(
        id = "main_tabs",
        selected = isolate(current_tab()),
        tabPanel("Monocytes/Macrophages",
      fluidRow(
        # Left column for vertical arrow (this image will span the height of the JPEG block)
        style = "display: flex; align-items: center;", 
        column(width = 2,
               div(
                 imageOutput("vert_arrow_image_1"),
                 style = "width:100%; height:500px; object-fit: contain; margin-left: 0; margin-bottom: 30%"
               )
        ),
        # Right column with nested fluidRows for JPEG images
        column(width = 10,
               fluidRow(
                 class = "flow-data-row flow-data-row-left-margin",
                 style = "margin-top: 10%; margin-bottom: 10%; padding: 0",
                 lapply(seq_along(top_jpg), function(i) {
                   div(
                     imageOutput(paste0("gene_image_jpg_", i)),
                     style = "display:inline-block; width:18%; margin: .01%; vertical-align: top; margin-left: 4%;"
                   )
                 })
               ),
               fluidRow(
                 div(
                   class = "flow-key-row-left-margin",
                   style = "margin-top: 1%; margin-bottom: 10%; padding: 0",
                   lapply(seq_along(top_keys), function(i) {
                     div(
                       imageOutput(paste0("color_key_png_", i, "_tab1")),
                       style = "display:inline-block; width:18%; margin: .01%; vertical-align: top; margin-left: 4%"
                     )
                   })
                 )
               )
        )
      ),
      # Horizontal arrow row (optional)
      fluidRow(
        column(12, div(imageOutput("horz_arrow_image_1"), class = "arrow-image")),
        style = "margin-top: clamp(-470px, -35vw, -180px); margin-left: 9%; margin-bottom: 0px; padding: 0px;",
        class = "arrow-row arrow-row-left-margin"
      ),
      fluidRow(
        style = "display: flex; flex-wrap: nowrap; margin-top: 5vw;",
        column(8,
               lapply(seq_along(png_other), function(i) {
                 div(
                   imageOutput(paste0("gene_image_png_other_", i)),
                   style = "text-align:left; margin-left: 30%;"
                 )
               })
        ),
        column(4,
               class = "color-key-row color-key-row-left-margin",
               div(
                 imageOutput("legend_vlnplot_tab1"),
                 style = "margin-left: 3%; width:1000px; height:1000px"
               )
        )
      )
     ),
     tabPanel("Macrophage Subsets",
              fluidRow(
                # Left column for vertical arrow (this image will span the height of the JPEG block)
                style = "display: flex; align-items: center;", 
                column(width = 2,
                       div(
                         imageOutput("vert_arrow_image_2"),
                         style = "width:100%; height:500px; object-fit: contain; 
                         margin-left: 0; margin-bottom: 30%"
                       )
                ),
                # Right column with nested fluidRows for JPEG images
                column(width = 10,
                         fluidRow(
                           class = "flow-data-row flow-data-row-left-margin",
                           style = "margin-top: 10%; margin-bottom: 10%; padding: 0;",
                           lapply(seq_along(bottom_jpg), function(i) {
                             div(
                               imageOutput(paste0("gene_image_jpg_", i + length(top_jpg))),
                               style = "display:inline-block; width:18%; margin: .01%; vertical-align: top; margin-left: 4%;"
                             )
                           })
                         ),
                       fluidRow(
                         div(
                           class = "flow-key-row-left-margin",
                           style = "margin-top: 1%; margin-bottom: 10%; padding: 0;",
                           lapply(seq_along(bottom_keys), function(i) {
                             div(
                               imageOutput(paste0("color_key_png_", i, "_tab2")),
                               style = "display:inline-block; width:18%; margin: .01%; vertical-align: top; margin-left: 4%"
                             )
                           })
                         )
                       )
                )
              ),
              # Horizontal arrow row (optional)
              fluidRow(
                column(12, div(imageOutput("horz_arrow_image_2"), class = "arrow-image2")),
                style = "margin-top: clamp(-470px, -35vw, -180px); margin-left: 9%; margin-bottom: 0px; padding: 0px;",
                class = "arrow-row2 arrow-row-left-margin2"
              ),
           fluidRow(
             style = "display: flex; flex-wrap: nowrap; margin-top: 5vw;",
             column(8,
                    lapply(seq_along(png_6cells), function(i) {
                      div(
                        imageOutput(paste0("gene_image_png_6cells_", i)),
                        style = "text-align:left; margin-left: 22%;"
                      )
                    })
             ),
             column(4,
                    class = "color-key-row color-key-row-left-margin",
                    div(
                      imageOutput("legend_vlnplot_tab2"),
                      style = "margin-left: 14%;"
                    )
             )
           )
        )
    )
  )
    
    # Render the JPEG images
    lapply(seq_along(jpg_files), function(i) {
      output[[paste0("gene_image_jpg_", i)]] <- renderImage({
        list(src = jpg_files[i], alt = basename(jpg_files[i]), class = "jpg-image")
      }, deleteFile = FALSE)
    })
    
    # Render the PNG images (unchanged)
    lapply(seq_along(png_other), function(i) {
      output[[paste0("gene_image_png_other_", i)]] <- renderImage({
        list(src = png_other[i], 
             alt = input$gene_choices,
             class = "png-image"
             )
      }, deleteFile = FALSE)
    })
    
    lapply(seq_along(png_6cells), function(i) {
      output[[paste0("gene_image_png_6cells_", i)]] <- renderImage({
        list(src = png_6cells[i], 
             alt = input$gene_choices, 
             class = "png-image")
      }, deleteFile = FALSE)
    })
    
    lapply(seq_along(color_key_PNGs), function(i) {
      output[[paste0("color_key_png_", i)]] <- renderImage({
        list(src = color_key_PNGs[i], alt = basename(color_key_PNGs[i]), class = "colorkey-image")
      }, deleteFile = FALSE)
    })
    
    lapply(seq_along(top_keys), function(i) {
      output[[paste0("color_key_png_", i, "_tab1")]] <- renderImage({
        list(src = top_keys[i],
             alt = basename(top_keys[i]),
             class = "colorkey-image")
      }, deleteFile = FALSE)
    })
    
    lapply(seq_along(bottom_keys), function(i) {
      output[[paste0("color_key_png_", i, "_tab2")]] <- renderImage({
        list(src = bottom_keys[i],
             alt = basename(bottom_keys[i]),
             class = "colorkey-image")
      }, deleteFile = FALSE)
    })
    
    output$horz_arrow_image_1 <- renderImage({
      list(src = horz_arrow_images,
           alt = basename(horz_arrow_images),
           style = "width: clamp(180px, 46vw, 700px); height: clamp(30px, 3.5vw, 60px); object-fit: fill;"
           )
    }, deleteFile = FALSE)

    output$vert_arrow_image_1 <- renderImage({
      list(src = vert_arrow_images,
           alt = basename(vert_arrow_images),
           style = "width: clamp(80px, 10vw, 160px); height: clamp(120px, 14vh, 220px); object-fit: contain;"
           )
    }, deleteFile = FALSE)

    output$horz_arrow_image_2 <- renderImage({
      list(src = horz_arrow_images,
           alt = basename(horz_arrow_images),
           style = "width: clamp(180px, 46vw, 700px); height: clamp(30px, 3.5vw, 60px); object-fit: fill;"
           )
    }, deleteFile = FALSE)

    output$vert_arrow_image_2 <- renderImage({
      list(src = vert_arrow_images,
           alt = basename(vert_arrow_images),
           style = "width: clamp(80px, 10vw, 160px); height: clamp(120px, 14vh, 220px); object-fit: contain;")
    }, deleteFile = FALSE)

    output$legend_vlnplot_tab1 <- renderImage({
      list(
        src = file.path(color_key_path, "Legend_VlnPlot.png"),
        alt = "Legend VlnPlot",
        style = "width: clamp(120px, 14vw, 240px); height: auto; object-fit: contain;"
      )
    }, deleteFile = FALSE)

    output$legend_vlnplot_tab2 <- renderImage({
      req(input$gene_choices != "Ms4a3-cre-Ai3-YFP")
      list(
        src = file.path(color_key_path, "Legend_VlnPlot_6cells.png"),
        alt = "Legend VlnPlot 6cells",
        style = "width: clamp(120px, 14vw, 240px); height: auto; object-fit: contain;"
      )
    }, deleteFile = FALSE)
    
    output$citation1 <- renderUI({
      HTML("Wang Y, Dowling S, Rodriguez V, Maciuch J, 
      Mayer M, Therron T, Gurra G, Shah C, Cuda CM, Makinde HKM, 
      Ginhoux F, Vohringer D, Harrington C, Lawrence T, 
      Grainger JR, Winter DR, Perlman H. <b>Comprehensive analysis of myeloid reporter mice.</b>")
    })
    
    output$citation2 <- renderUI({
      HTML("Wang Y, Dowling S, Rodriguez V, Maciuch J, 
      Mayer M, Therron T, Gurra G, Shah C, Cuda CM, Makinde HKM, 
      Ginhoux F, Vohringer D, Harrington C, Lawrence T, 
      Grainger JR, Winter DR, Perlman H. <b>Comprehensive analysis of myeloid reporter mice.</b>")
    })
    
    combined_images
  })
  
  

}

# Run the application
shinyApp(ui = ui, server = server)
