# app.R
library(shiny)
library(dplyr)
library(DT)
library(NGLVieweR)



df <- read.delim(
  "../data/ted_clean/regions_with_parch_fnote.tsv",
  sep             = "\t",
  header          = TRUE,
  stringsAsFactors = FALSE
) %>% filter(tipo == "link" | tipo == "internal_link" ) %>%
    mutate(across(where(is.character), ~ gsub("[<>]", "x", .)))%>%
    mutate(arch_length = nchar(lite_cat_20)) %>%
  arrange(
    arch_length,
    desc(afd_dc),
    desc(lite_dc)
  )





# 2) columnas de interés y filtro inicial
cols <- c(
  'acc','feature','tipo','start','end','length','lite_dc','afd_dc','lite_cat_20', 'lite_cat_10'
  #,
  #'lite_cat_10','afd_cat_20','fcr','ncpr','isop','fn','fp','fe','kappa',
  #'omega','delta'
)

df <- df %>%
  filter(lite_dc >= 0.7) %>%
  select(all_of(cols))


ui <- fluidPage(
  titlePanel("Disordered linkers"),
  
  # filtros arriba
  fluidRow(
    column(2,
      sliderInput("lite_dc", "Disorder window (MobiDB Lite)",
                  min = 0.7, max = 1, value = c(0.7,1), step = 0.01)
    ),
    column(2,
    div(style = "display: flex; align-items: center;",
      # le damos al textInput un ancho fijo menor
      div(style = "flex: 1; margin-right: 10px;",
        textInput("lite_cat_20", "Protein Architecture", value = "A1A"),
        checkboxInput("exact_arch", "Exact match", value = FALSE)
      ),
      # el checkbox al lado
      div(style = "white-space: nowrap;",
        #checkboxInput("exact_arch", "Exact match", value = FALSE)
      )
    )
    ),
    column(3,
      sliderInput("length", "Min Length",
                  min = floor(min(df$length, na.rm=TRUE)),
                  max = 300,
                  value = c(floor(min(df$length, na.rm=TRUE)), max(df$length, na.rm=TRUE)),
                  step = 1)
    ),
    column(2, downloadButton("download_csv", "Download"))
  ),
  # checkboxes inline para columnas visibles
  fluidRow(
    column(12,
      checkboxGroupInput(
        "cols_visible", "Columnas visibles",
        choices  = cols,
        selected = c('acc','feature','length','lite_dc','afd_dc','lite_cat_20', 'tipo',
                     'lite_cat_10'
                     #,'fcr','ncpr','fe','kappa','omega'
                     ),
        inline   = TRUE
      )
    )
  ),
  tags$hr(),
  
  # layout: tabla a la izquierda, visor NGL a la derecha
  fluidRow( style='margin-left:5px; width:95%',
    column(8,
    
      DTOutput("tabla"),
     
      
    ),
    column(4,
      NGLVieweROutput("nglviewer", width = "95%", height = "600px")
    )
  )
)

server <- function(input, output, session) {
  # reactive con el df filtrado (sin formatear)
#   base_data <- reactive({
#     df %>%
#       filter(
#         lite_dc >= input$lite_dc[1],
#         lite_dc <= input$lite_dc[2],
#         # manejo del filtro de arquitectura con exact/contains
#         ( input$lite_cat_20 == "" ) |
#         ( input$exact_arch  & lite_cat_20 == input$lite_cat_20 ) |
#         ( !input$exact_arch & grepl(input$lite_cat_20, lite_cat_20, fixed = TRUE) ),
#         length  >= input$length
#       )
#   })

base_data <- reactive({
  # 1) Partir el texto por comas y limpiar espacios
  pats <- if (input$lite_cat_20 == "") {
    character(0)
  } else {
    trimws(unlist(strsplit(input$lite_cat_20, ",")))
  }

  df %>%
    filter(
      # tus otros filtros
      lite_dc >= input$lite_dc[1],
      lite_dc <= input$lite_dc[2],
      length  >= input$length[1],
      length  <= input$length[2],
      # 2) nuevo filtro de arquitectura múltiple
      (
        # sin patrón ingresado → todo pasa
        length(pats) == 0
        # exact match → pertenece a alguno de los patrones
        | ( input$exact_arch & lite_cat_20 %in% pats )
        # contains match → contiene alguno de los patrones
        | ( !input$exact_arch & 
            rowSums(sapply(pats, function(p) grepl(p, lite_cat_20, fixed = TRUE))) > 0
          )
      )
    )
})

  
  # reactive con el df formateado para la tabla
  table_data <- reactive({
    base_data() %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      mutate(
        feature = sprintf(
          '<a href="https://ted.cathdb.info/uniprot/%s" target="_blank">%s</a>',
          acc, feature
        ),
        acc = sprintf(
          '<a href="https://mobidb.org/%s" target="_blank">%s</a>',
          acc, acc
        )
      ) %>%
      select(all_of(input$cols_visible))
  })
  
  # renderizamos la tabla con DT
  output$tabla <- renderDT({
    datatable(
      table_data(),
      escape    = FALSE,
      filter    = "top",
      selection = "single",
      options   = list(pageLength = 10, scrollX = TRUE)
    )
  }, server = TRUE)
  
  # descarga el PDB de AlphaFold al seleccionar una fila
  pdb_file <- eventReactive(input$tabla_rows_selected, {
    sel <- input$tabla_rows_selected
    
    req(sel)
    row <- base_data()[sel, ]
    print(row)
    tmp <- file.path(tempdir(), "latest_af.pdb")
    url <- sprintf(
      "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v4.pdb",
      toupper(row$acc)
    )
    download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
    list(path = tmp, start = row$start, end = row$end)
  })

    output$download_csv <- downloadHandler(
    filename = function() {
      paste0("linkers_filtrados_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Tomamos el df filtrado, seleccionamos columnas visibles y redondeamos
      df_out <- base_data() %>%
        select(all_of(input$cols_visible)) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3)))
      write.csv(df_out, file, row.names = FALSE)
    }
  )

  
  # renderizamos el visor NGLVieweR
  output$nglviewer <- renderNGLVieweR({
    info <- pdb_file()
    req(info)
    NGLVieweR(info$path) %>%
      stageParameters(backgroundColor = "white") %>%
      addRepresentation("cartoon", param = list(color = "grey")) %>%
      addRepresentation("cartoon",
        param = list(sele = paste0(info$start, "-", info$end),
                     color = "red")
      )
  })
}

shinyApp(ui, server)
