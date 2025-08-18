library(shiny)
library(ggplot2)
library(DT)
library(shinythemes)
library(scales)
library(reactable)
library(yarrr)
library(ggpirate)
library(shinycssloaders)
library(shinyjs)
library(shinyThings)

# Define UI for application
ui <- navbarPage(
  title = div(
    style = "font-size: 20px; font-weight: bold;",  # Adjust font size and weight
    "SUPEROMICS: Supermere, Extracellular Vesicle, and Exomere Omics Explorer"
  ),  # Navbar title
  theme = shinytheme("cerulean"),         
  
  # Home Tab 
  tabPanel(
  theme = shinytheme("cerulean"),         
    "Home",
    sidebarLayout(
      fluid = TRUE,
      sidebarPanel(
        width = 3,  
        fluid = TRUE,
        shinyThings::radioSwitchButtons("dataset", "Dataset:",
                    choices = c("Proteomics", "RNA-seq", "Lipidomics"),
                    #inline = TRUE,
                    selected = "Proteomics"), 
        uiOutput("dynamicUI"),
        shinyThings::radioSwitchButtons("grouping", "Sample Grouping:", 
                    choices = c("Individual", "Grouped"),
                    selected = "Grouped"),
        div(
          style = "margin-top: 8px; font-weight: bold; font-size: 14px;",
          "Plot Faceting:"
        ),
        checkboxInput("facet_isolation", "Isolation Method", value = FALSE),
        checkboxInput("facet_growth", "Growth Conditions", value = FALSE),
        downloadButton("downloadPlot", "Download Plot"),
        style = "padding: 20px;"  # Add padding for better layout
      ),
      mainPanel(
        fluidRow(
          column(12, withSpinner(plotOutput("plot", height = "400px")))  
        ),
        fluidRow(
          column(12, withSpinner(reactableOutput("datatable")))
        ),
        style = "padding: 20px;"  
      )
    )
  ),
  
  # About Tab
  tabPanel(
    "About",
    fluidPage(
      h3("About This App"),
      div(
        p("This app was developed to accompany the publication ", 
          tags$em("'A comprehensive analysis of supermere, exomere, and extracellular vesicle isolation and cargo in colorectal cancer'"), 
          " [DOI:XXXXXXXX]."),
        p("It provides an interactive platform to explore EV and NVEP Omics datasets generated from DiFi cells across various data types, including:"),
        tags$ul(
          tags$li(tags$b("Proteomics:"), " Supports searching for individual proteins using Gene Name nomenclature. Data is reported as median normalized counts."),
          tags$li(tags$b("RNA-seq:"), " Provides an overview of all host genome small RNA types (sRNA type) and browsing expression for each individual sRNA by type (e.g., miRNA, lncRNA, snRNA). Data is reported as reads per million total reads."),
          tags$li(tags$b("Lipidomics:"), " Enables an overview of all lipids categorized by Category, Class, Subclass, or Species. Data is reported as log2-normalized counts.")
        )
      ),
      div(
        h4("Features:"),
        tags$ul(
          tags$li(tags$b("Sample Grouping:"), " Choose 'Individual' for sample-level dot plots or 'Grouped' for pirate plots aggregated by particle type."),
          tags$li(tags$b("Faceting Options:"), " Customize plots by Isolation Method (UC vs FPLC) and Growth Conditions (2D vs 3D).")
        )
      ),
      div(
        h4("Resources:"),
        tags$ul(
          tags$li("Refer to the publication for detailed explanations on methods used for data acquisition and preprocessing."),
          tags$li("Raw data available at publicly accessible databases: proteomics at ", tags$a(href = "https://www.peptideatlas.org/PASS/PXD066872", "PeptideAtlas"), 
                  ", RNA-seq at ", tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234567", "GEO"), 
                  ", and lipidomics at ", tags$a(href = "https://www.lipidmaps.org/data/lipidomics/", "LIPID MAPS"),
          tags$li("Access the source code and datasets on ", tags$a(href = "https://github.com/rskks/GOFiltering/tree/main/browse", "GitHub"), ".")
        )
      )
    )
  ),
  
  # Contact Tab
  tabPanel(
    "Contact",
    fluidPage(
      h3("Contact Us"),
      p("Please refer to the original ", tags$a(href = "https://github.com", "publication"), " for the corresponding author's contact information."),
      p("If you have questions or suggestions about the website, you can email ", tags$a(href = "mailto:ostutanov@gmail.com", "Oleg Tutanov"), ".")
    )
  )
)
  

# Define server logic
server <- function(input, output, session) {
  plotOutput("plot") %>% withSpinner(color = "#007FFF")
  shinyThings::updateRadioSwitchButtons(session = session, 
                                        inputId = "dataset",
                                        selected = "Proteomics")
  
  # Load the data
  protein_data <- reactive({
    file_path <- "protein_data.csv"
     cat("Loading data from:", file_path, "\n")  # Debugging output
    if (file.exists(file_path)) {
      cat("File exists. Reading the data...\n")
      data <- read.csv(file_path)
      cat("Data loaded. Converting to data frame...\n")
      data <- as.data.frame(data, check.names = FALSE)
      rownames(data) <- data$protein
      data <- data[,-1]  # Remove the $protein column
      
      # Convert columns to numeric and log non-numeric values
      for (col in names(data)) {
        suppressWarnings({
          numeric_values <- as.numeric(data[[col]])
          non_numeric_indices <- which(is.na(numeric_values))  # Identify problematic rows
          
          if (length(non_numeric_indices) > 0) {
            cat("Non-numeric values found in column '", col, "' at rows: ", 
                paste(non_numeric_indices, collapse = ", "), "\n", sep = "")
          }
          
          # Replace non-numeric values with 0
          numeric_values[is.na(numeric_values)] <- 0
          data[[col]] <- numeric_values
        })
      }
      
      data <- round(data, 2)  # Round data to 2 decimal places
      cat("Data processing complete. Returning data...\n")
      return(data)
    } else {
      stop("File not found: ", file_path)
    }
  })
  protein_meta <- reactive({
    file_path <- "protein_meta2_fixed.csv"
    if (file.exists(file_path)) {
      data <- read.csv(file_path, row.names = 1)
      #data <- t(data)
      data <- as.data.frame(data, check.names = FALSE)
      return(data)
    } else {
      stop("File not found: ", file_path)
    }
  })
  
  # RNA-seq data (updated for consistency)
  rna_data <- reactive({
    req(input$rnaseq_type) # Ensure input exists
    
    # Map input to file paths
    file_path <- switch(input$rnaseq_type,
                        "sRNA type" = "rna_data_cat.csv",
                        "miRNA" = "rna_data_mirna.csv",
                        "lncRNA" = "rna_data_lncrna.csv",
                        "snRNA" = "rna_data_snrna.csv",
                        "tRNA" = "rna_data_trna.csv",
                        "snoRNA" = "rna_data_snorna.csv",
                        "rRNA" = "rna_data_rrna.csv",
                        "yRNA" = "rna_data_yrna.csv"
                        )
    
    if (file.exists(file_path)) {
      data <- as.data.frame(read.csv(file_path, stringsAsFactors = FALSE))
      
      # Debugging output
      cat("Loaded data from:", file_path, "\n")
      #cat("Column names:\n", colnames(data), "\n")
      #cat("First few rows:\n")
      #print(head(data))
      
      # Conditional logic for row name processing
      if ("rna" %in% colnames(data)) {
        if (input$rnaseq_type != "lncRNA") {
          # Process row names (skip for lncRNA)
          data$rna <- gsub("[|:;].*", "", data$rna)
          rownames(data) <- make.unique(data$rna)  # Ensure unique row names
          data <- data[, -1]  # Remove the 'rna' column
          cat("Row names set successfully.\n")
        } else {
          cat("Skipping row name processing for lncRNA.\n")
          rownames(data) <- make.unique(data$rna)  # Ensure unique row names
          data <- data[, -1]  # Remove the 'rna' column
        }
      } else {
        stop("Column 'rna' not found in the data.")
      }
      
      return(data)
    } else {
      showNotification("Selected RNA-seq file not found.", type = "error")
      # Return an empty data frame to avoid issues downstream
      return(data.frame())
    }
  })
  
  # Lipid data (updated for consistency)
  lipid_data <- reactive({
    req(input$lipidomics_type) # Ensure input exists
    
    file_path <- switch(input$lipidomics_type,
                        "Category" = "lipid_data_cat.csv",
                        "Class" = "lipid_data_cl.csv",
                        "Subclass" = "lipid_data_subcl.csv",
                        "Species" = "lipid_data_spe.csv")
    
    if (file.exists(file_path)) {
      data <- as.data.frame(read.csv(file_path, stringsAsFactors = FALSE))
      
      # Debugging output
      cat("Loaded data from:", file_path, "\n")
      #cat("Column names:\n", colnames(data), "\n")
      #cat("First few rows:\n")
      #print(head(data))
      
      if ("lipid" %in% colnames(data)) {
        rownames(data) <- data$lipid
        data <- data[,-1]  # Remove the 'lipid' column
        cat("Row names set successfully.\n")
      } else {
        stop("Column 'lipid' not found in the data.")
      }
      
      return(data)
    } else {
      showNotification("Selected lipidomics file not found.", type = "error")
      # Return an empty data frame to avoid issues downstream
      return(data.frame())
    }
  })
  

  # Dynamic UI rendering for Proteomics, RNA-seq, or Lipidomics (with sub-types for Lipidomics)
  output$dynamicUI <- renderUI({
    if (input$dataset == "Lipidomics") {
      # Lipidomics-specific UI
      tagList(
        selectInput("lipidomics_type", "Lipidomics Type:",
                    choices = c("Category", "Class", "Subclass", "Species"),
                    selected = "Category"),
        selectizeInput("lipid", "Lipid:", choices = NULL) # Initially empty
      )
    } else if (input$dataset == "Proteomics") {
      selectizeInput("protein", "Protein:", choices = NULL) # Placeholder
    } else if (input$dataset == "RNA-seq") {
      # RNAseq-specific UI
      tagList(
        selectInput("rnaseq_type", "RNAseq Dataset:",
                    choices = c("sRNA type", "miRNA", "lncRNA", "snRNA", "tRNA", "snoRNA", "rRNA", "yRNA"),
                    selected = "sRNA type"),
        selectizeInput("rna", "RNA:", choices = NULL) # Placeholder
      )
    }
  })
  
  # Update lipidomics dropdown when lipidomics_type changes
  observeEvent(input$lipidomics_type, {
    req(input$lipidomics_type) # Ensure lipidomics type is selected
    data <- lipid_data()
    req(data) # Ensure lipid_data is loaded
    
    updateSelectizeInput(session, "lipid", choices = rownames(data), server = TRUE)
  })
  
  # Update proteomics dropdown when protein_data changes
  observe({
    data <- protein_data()
    req(data) # Ensure protein_data is loaded
    
    updateSelectizeInput(session, "protein", choices = rownames(data), server = TRUE)
  })
  
  # Update RNA-seq dropdown when rna_data changes
  observe({
    data <- rna_data()
    req(data) # Ensure rna_data is loaded
    
    updateSelectizeInput(session, "rna", choices = rownames(data), server = TRUE)
  })
  
  # Observe event to update inputs on dataset switching
  observeEvent(input$dataset, {
    # Reset the dropdown menu for proteins
    if (input$dataset == "Proteomics") {
      updateSelectInput(session, "protein", 
                        choices = rownames(protein_data()),
                        selected = NULL)
    }
    
    # Reset the dropdown menu for RNA
    if (input$dataset == "RNA-seq") {
      updateSelectInput(session, "rna", 
                        choices = rownames(rna_data()),
                        selected = NULL)
    }
    
    # Reset the dropdown menu for lipids
    if (input$dataset == "Lipidomics") {
      lipid_choices <- rownames(lipid_data())
      updateSelectInput(session, "lipid", 
                        choices = lipid_choices,
                        selected = NULL)
    }
  })
  
  custom_theme <- theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),  # X-axis tick labels
    axis.text.y = element_text(size = 13),  # Y-axis tick labels
    axis.title.x = element_text(size = 17),  # X-axis label font size
    axis.title.y = element_text(size = 13),  # Y-axis label font size
    plot.title = element_text(size = 19, face = "bold"),  # Title font size
    strip.text = element_text(size = 15, face = "bold"),  # Facet labels font size
    legend.text = element_text(size = 14),  # Legend text size
    legend.title = element_text(size = 15, face = "bold")  # Legend title font size
  )
  
  
  
  apply_faceting <- function(plot, input) {
    if (input$facet_isolation & input$facet_growth) {
      plot + facet_grid(rows = vars(Isolation), cols = vars(Growth))
    } else if (input$facet_isolation) {
      plot + facet_grid(cols = vars(Isolation))
    } else if (input$facet_growth) {
      plot + facet_grid(cols = vars(Growth))
    } else {
      plot
    }
  }
  
  prepare_total_data <- function(data, meta_data) {
    data.frame(
      Condition = colnames(data),
      Value = as.numeric(data),
      Individual = factor(unlist(meta_data["Individual", colnames(data)])),
      Particle = factor(unlist(meta_data["Particle", colnames(data)])),
      Isolation = factor(unlist(meta_data["Isolation", colnames(data)])),
      Growth = factor(unlist(meta_data["Growth", colnames(data)]))
    )
  }
  
  custom_labels <- function(total_data) {
    function(x) {
      sapply(x, function(i) {
        paste(
          unique(total_data$Particle[total_data$Individual == i]), 
          unique(total_data$Isolation[total_data$Individual == i]), 
          unique(total_data$Growth[total_data$Individual == i]), 
          i
        )
      })
    }
  }
  
  # Reactive plot generation
  plot_reactive <- reactive({
    data <- switch(input$dataset,
                   "Proteomics" = protein_data(),
                   "RNA-seq" = rna_data(),
                   "Lipidomics" = lipid_data())
    
    target_input <- switch(input$dataset,
                           "Proteomics" = input$protein,
                           "RNA-seq" = input$rna,
                           "Lipidomics" = input$lipid)
    
    y_label <- switch(input$dataset,
                      "Proteomics" = "normalized counts",
                      "RNA-seq" = "reads per million total reads",
                      "Lipidomics" = "log-2 normalized counts")
    
    if (!is.null(target_input) && target_input %in% rownames(data)) {
      selected_data <- data[target_input, , drop = FALSE]
      meta_data <- protein_meta()
      total_data <- prepare_total_data(selected_data, meta_data)
      
      if (input$grouping == "Individual") {
        p <- ggplot(total_data, aes(x = Individual, y = Value, color = Particle)) +
          geom_point(position = position_jitter(width = 0.1), size = 3) +
          labs(
            y = y_label,
            x = NULL,
            title = paste("Expression of", target_input)
          ) +
          scale_x_discrete(labels = custom_labels(total_data)) +
          custom_theme
      } else {
        p <- ggplot(total_data, aes(x = Particle, y = Value)) +
          geom_pirate(aes(colour = Particle), bars = FALSE,
                      points_params = list(shape = 19, alpha = 0.2),
                      lines_params = list(size = 0.8)) +
          labs(
            y = y_label,
            title = paste("Grouped Expression for", target_input)) +
          custom_theme
      }
      
      p <- apply_faceting(p, input)
      return(p)
    } else {
      ggplot() +
        labs(title = "No data available") +
        custom_theme
    }
  })
  
  # Render plot in the UI
  output$plot <- renderPlot({
    print(plot_reactive())
  })
  
  # Download handler for the plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot_", input$dataset, "_", input$protein, ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_reactive(), device = "png", width = 16, height = 9)
    }
  )
  
  
  data <- reactive({
    switch(input$dataset,
           "Proteomics" = protein_data(),
           "RNA-seq" = rna_data(),
           "Lipidomics" = lipid_data())
  })
  
  # Render the reactable table
  output$datatable <- renderReactable({
    df <- data()
    
    if (nrow(df) == 0 || ncol(df) == 0) {
      cat("Data is empty. Returning a placeholder table.\n")
      return(reactable(data.frame(Message = "No data available."), 
                       pagination = FALSE, searchable = FALSE))
    }
    
    cat("Rendering table with valid data.\n")
    reactable(
      df,
      defaultColDef = colDef(format = colFormat(digits = 2)),
      searchable = TRUE,
      highlight = TRUE,
      pagination = TRUE,
      resizable = TRUE,
      wrap = FALSE,
      defaultPageSize = 10,
      pageSizeOptions = c(5, 10, 20),
      theme = reactable::reactableTheme(
        headerStyle = list(backgroundColor = '#f5f5f5', color = '#333'),
        cellStyle = list(backgroundColor = '#fff'),
        highlightColor = '#007ba7'
      )
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
