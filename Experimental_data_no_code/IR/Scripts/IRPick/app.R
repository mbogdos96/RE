library(shiny)
library(plotly)

ui <- fluidPage(
  titlePanel("IR Plot App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose an IR data file (CSV format):"),
      actionButton("clear", "Clear Selection")
    ),
    mainPanel(
      plotlyOutput("ir_plot"),
      verbatimTextOutput("wavenumbers")
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, header = FALSE)
    colnames(df) <- c("Wavenumber", "Transmission")
    df
  })
  
  selected_peaks <- reactiveVal(NULL)
  
  observeEvent(event_data("plotly_click"), {
    x <- event_data("plotly_click")$x
    x <- round(x)
    
    current_peaks <- selected_peaks()
    
    if (is.null(current_peaks)) {
      current_peaks <- x
    } else {
      current_peaks <- unique(c(current_peaks, x))
    }
    
    selected_peaks(current_peaks)
  })
  
  observeEvent(input$clear, {
    selected_peaks(NULL)
  })
  
  observeEvent(input$file, {
    selected_peaks(NULL)
  })
  
  output$ir_plot <- renderPlotly({
    x <- data()$Wavenumber
    y <- data()$Transmission
    
    plot_data <- data.frame(x = x, y = y)
    
    # Highlight selected peaks
    plot_data$peak <- ifelse(x %in% selected_peaks(), y, NA)
    
    p <- plot_ly(data = plot_data, x = ~x, y = ~y, type = 'scatter', mode = 'lines', name = 'IR Spectrum') %>%
      add_trace(y = ~peak, name = 'Selected Peaks', marker = list(color = 'red'))
    
    p %>%
      layout(
        xaxis = list(title = 'Wavenumber [cm^-1]', autorange = 'reversed'), # Reverse the x-axis
        yaxis = list(title = 'Transmission'),
        showlegend = TRUE
      )
  })
  
  output$wavenumbers <- renderText({
    if (!is.null(selected_peaks())) {
      selected_wavenumbers <- round(selected_peaks())
      paste("Selected Wavenumbers: ", paste(selected_wavenumbers, collapse = ", "))
    }
  })
}

shinyApp(ui, server)
