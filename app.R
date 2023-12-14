# Loading necessary libraries for UI components
library(shiny)
library(ggplot2)
library(shinydashboard)
library(tidyverse)
library(readxl) # For readxl
library(ggpubr)
library(cowplot)
library(DT) # For enabling scrolling in DataTable
library(plotly) # For interactive plotting
library(shinyWidgets)
library(shinyBS)
library(RSQLite)
library(shinycssloaders)

# Define UI layout for Shiny app
ui <- dashboardPage(
    dashboardHeader(title = "Fowler Lab"), # Dashboard header with title
    dashboardSidebar(
        # Sidebar with navigation items
        sidebarMenu(
            menuItem("Home", tabName = "home", icon = icon("home")),
            menuItem("About", tabName = "about", icon = icon("info-circle"))
        )
    ),
    dashboardBody(
        # Main body of the dashboard
        tags$style(".small-box.bg-green { background-color: #808000 !important; color: #FFFFFF !important; }"),
        tags$style(".small-box.bg-maroon { background-color: #E8A11C !important; color: #FFFFFF !important; }"),
        tags$style(".small-box.bg-blue { background-color: #695B73 !important; color: #FFFFFF !important; }"),
        tabItems(
            # Home tab
            tabItem("home",
                fluidRow(
                    box(
                        title = "Welcome to the Gene-wise Prediction Explorer!",
                        HTML("This interactive tool allows you to explore gene-wise prediction using a dataset of your choice. Select the dataset, predictor, and the gene you want to investigate. The app will provide you with visualizations of predictor score distributions, helping you gain a deeper understanding of predictor performance at the gene level.<i> The density plot is estimated with scores from all SNVs.</i>"),
                        width = 12,
                        collapsible = TRUE
                    )
                ),
                fluidRow(
                    column(width = 2, id="datasetColumn", selectInput("dataset", "Select a Dataset:",
                        choices = c("One Star Missense", "ClinVar 2023 Dataset without training variants", "ClinVar 2023 Dataset", "ClinGen SVI Calibration Dataset")),
                        bsTooltip('datasetColumn', 'Choose a dataset. Dataset descriptions can be found on the about page.', placement = "top", trigger = "hover")),
                    column(width = 2, id="predictorColumn", selectInput("predictor", "Select a Predictor:",
                        choices = c('REVEL', 'BayesDel')), bsTooltip('predictorColumn', 'Choose between exploring the REVEL or BayesDel predictors.', placement = "top", trigger = "hover")),
                    column(width = 2, id="geneColumn", selectizeInput("selectedGene", "Select a Gene:", 
                        choices = NULL,
                        options = list(
                            placeholder = "Search for a gene", # Placeholder text
                            maxOptions = 10 # Maximum number of displayed options
                        )),
                        bsTooltip('geneColumn', 'Select a gene from the drop-down list, or type to search for a gene.', placement = "top", trigger = "hover")),
                        tags$head(tags$style(HTML(".small-box {height: 100px}"))),
                        valueBoxOutput("concordant", width = 2),
                        valueBoxOutput("discordant", width = 2),
                        valueBoxOutput("indeterminate", width = 2)
                ),
                fluidRow(
                    column(width = 6,
                        box(
                            title = "Aggregate Genome Plot",
                            plotlyOutput("render_wg_plot") %>% withSpinner(), # Plot output area
                            width = NULL,
                            collapsible = TRUE,
                            downloadButton('downloadWGSummaryTable', "Download Summary Table")
                        )
                    ),
                    column(width = 6,
                        box(
                            title = "Gene-wise Plot",
                            plotlyOutput("render_gene_wise_plot") %>% withSpinner(), # Plot output area
                            width = NULL,
                            collapsible = TRUE,
                            downloadButton('downloadSummaryTable', "Download Summary Table"),
                            id = "geneWisePlot"
                        ),
                        bsTooltip('geneWisePlot', 'Hover over colored bars to view summary statistics. Click to zoom in.', placement = "top", trigger = "hover")
                    )
                ),
                fluidRow(
                    box(
                        width = 12,
                        title = "Data Table",
                        collapsible = TRUE,
                        fluidRow(
                            column(6,
                                pickerInput(
                                    inputId = "categories",
                                    label = "Select Categories",
                                    choices = c("Concordant", "Discordant", "Indeterminate"),
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)
                                )
                            ),
                            column(6,
                                pickerInput(
                                    inputId = "intervals",
                                    label = "Select Interval(s)",
                                    choices = NULL,
                                    multiple = TRUE,
                                    options = list(`actions-box` = TRUE)
                                )
                            )
                        ),
                        dataTableOutput("data_table"),
                        downloadButton('downloadTable', "Download Table")
                    )
                )
            ),
            # About tab
            tabItem("about",
                fluidRow(
                    box(
                        title = "Welcome to the Gene-wise Prediction Explorer!",
                        HTML("This interactive tool allows you to explore gene-wise prediction using a dataset of your choice. Select the dataset, predictor, and the gene you want to investigate. The app will provide you with visualizations of predictor score distributions, helping you gain a deeper understanding of predictor performance at the gene level.<i> The density plot is estimated with scores from all SNVs.</i>"),
                        width = 12,
                        collapsible = TRUE
                    )
                ),
                fluidRow(
                    box(
                        title = "About the Data",
                        HTML("<p><strong>ClinGen SVI Calibration Dataset:</strong><br>
                             The training (ClinVar 2019) and test (ClinVar 2020) data sets used by the ClinGen SVI for predictor calibrations were downloaded from their supplement. The datasets were combined to create a comprehensive dataset of 20,948 variants across 2,711 genes.</p>
                             <p><strong>ClinVar 2023 Dataset:</strong><br>
                             A VCF file containing all variants present in ClinVar was downloaded from the ClinVar database in August 2023. This VCF file was annotated with REVEL predictor scores, BayesDel predictor scores, and gnomAD allele frequencies using Annovar. Variants meeting the following criteria were retained: 1+ stars, non-VUS, and missense. Genes without any pathogenic variants and variants with allele frequencies exceeding 0.01 were excluded. Exome allele frequency was used unless only genome allele frequency was available. The final filtered dataset consisted of 57,628 variants across 2,971 genes.</p>
                             <p><strong>ClinVar 2023 Dataset without training variants:</strong><br>
                             REVEL and BayesDel training variants were inferred by identifying variants that were missing from ClinGen SVI calibration dataset when cross-referenced with the December 2020 ClinVar download. Variants present in the ClinVar 2020 download, but not the ClinGen SVI calibration dataset were flagged as training variants and excluded from the ClinVar 2023 Dataset. This dataset consisted of 39,517 variants across 2,820 genes.</p>"),
                        width = 12,
                        collapsible = TRUE
                    )
                ),
                fluidRow(
                    box(
                        title = "Fowler Lab",
                        HTML("The <a href='https://fowlerlab.gs.washington.edu/'>Fowler Lab</a> is in the <a href='https://www.gs.washington.edu/'>Department of Genome Sciences</a> at the University of Washington, located in Seattle, WA."),
                        width = 12,
                        collapsible = TRUE
                    )
                )
            )
        )
    )
)

# Define server function
server <- function(input, output, session) {
  # Load pre-filtered variant data
  # agg_bd_data_clinvar23 <- read_csv("agg_bd_data.csv")
  # agg_revel_data_clinvar23 <- read_csv("new_agg_revel_data.csv")

  # Create a dictionary-like structure using named lists
  fileDictionary <- list(
    "ClinVar 2023 Dataset" = list(var1 = "Value1", var2 = "Value2"), # Placeholder variables until cleaned data is generated
    "One Star Missense" = list(var1 = "revel_one_star_missense.csv", var2 = "bd_one_star_missense.csv"),
    "ClinVar 2023 Dataset without training variants" = list(var1 = "revel_no_training_var.csv", var2 = "bd_no_training_var.csv"),
    "ClinGen SVI Calibration Dataset" = list(var1 = "Value1", var2 = "Value2") # Placeholder variables until cleaned data is generated
  )

  # Load and process one_star_missense data
  # Previously used clinvar_2023_annovar_37_parsed_one_star_missesne_no_vus.csv
  one_star_missense <- read_csv('clinvar_2023_annoavar_37_parsed_one_star_missense_no_vus_AF_filtered_w_geneid_cleanedup_nogenedup.csv') %>%
    mutate(
      REVEL = as.numeric(as.character(REVEL)),
      BayesDel_noAF_score = as.numeric(as.character(BayesDel_noAF_score))
    )

  # Load and process no_training_variants data
  # Previously used clinvar23_no_training_var_AF_fixed.csv
  no_training_variants <- read_csv("clinvar23_no_training_var_AF_fixed_w_geneid_cleaned_up.csv") %>%
    mutate(REVEL = as.numeric(as.character(REVEL)))

  # Load and combine ClinGen data
  clinvar19 <- read_excel('Clinvar_2019_dataset_080823.xlsx')
  clinvar20 <- read_excel('Clinvar_2020_dataset_080823.xlsx')
  clingen <- rbind(clinvar19, clinvar20) %>%
    mutate(REVEL_score = as.numeric(as.character(REVEL_score)))

  # Load and process clinvar23 data
  clinvar23 <- read_csv('clinvar_2023_annoavar_37_parsed_one_star_missense_no_vus_AF_filtered.csv') %>%
    mutate(
      REVEL = as.numeric(as.character(REVEL)),
      BayesDel_noAF_score = as.numeric(as.character(BayesDel_noAF_score))
    )

  datasetDictionary <- list(
    "ClinVar 2023 Dataset" = clinvar23,
    "One Star Missense" = one_star_missense,
    "ClinVar 2023 Dataset without training variants" = no_training_variants,
    "ClinGen SVI Calibration Dataset" = clingen
  )

  # Create lists of gene names for each dataframe
  gene_ntr <- unique(no_training_variants$Gene)
  gene_osm <- unique(one_star_missense$Gene)
  gene_clingen <- unique(clingen$genename)
  gene_clinvar23 <- unique(clinvar23$Gene)

  # Special list of ACMG genes
  ACMG_genes <- c('ACTA2', 'CSF1R', 'ACTC1', 'ACVRL1', 'APC', 'APOB', 'ATP7B', 'BAG3', 'BMPR1A', 'BRCA1', 'BRCA2', 'BTD', 'CACNA1S', 'CASQ2', 'COL3A1', 'DES', 'DSG2', 'DSC2', 'DSP', 'ENG', 'FBN1', 'FLNC', 'GAA', 'GLA', 'HFE', 'HNF1A', 'KCNH2', 'KCNQ1', 'LDLR', 'LMNA', 'MAX', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'MYBPC3', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NF2', 'OTC', 'PALB2', 'PCSK9', 'PKP2', 'PMS2', 'PRKAG2', 'PTEN', 'RB1', 'RBM20', 'RET', 'RPE65', 'RYR1', 'RYR2', 'SCN5A', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD3', 'SMAD4', 'STK11', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMEM43', 'TNNC1', 'TNNI3', 'TNNT2', 'TP53', 'TPM1', 'TRDN', 'TSC1', 'TSC2', 'TTN', 'TTR', 'VHL', 'WT1')

  # Retrieve filtered REVEL data based on the selected gene
  filtered_revel_data <- reactive({
    con <- dbConnect(RSQLite::SQLite(), dbname = "data")
    query <- glue::glue("
      SELECT REVEL
      FROM gene_wise_revel_scores
      WHERE Gene = '{input$selectedGene}'
    ")
    data <- dbGetQuery(con, query)
    return(data)
  })

  # Retrieve filtered BayesDel data based on the selected gene
  filtered_bd_data <- reactive({
    con <- dbConnect(RSQLite::SQLite(), dbname = "data")
    query <- glue::glue("
      SELECT BayesDel_nsfp33a_noAF
      FROM gene_wise_bd_scores
      WHERE Gene = '{input$selectedGene}'
    ")
    data <- dbGetQuery(con, query)
    return(data)
  })

  # Retrieve all gene names from genomic_coordinates
  genes <- reactive({
    con <- dbConnect(RSQLite::SQLite(), dbname = "data")
    query <- 'SELECT DISTINCT Gene FROM genomic_coordinates;'
    data <- dbGetQuery(con, query)
    return(data)
  })

  # Populate the gene selection dropdown at startup
  observe({
    query_results <- genes()$Gene
    query_results <- query_results[2:length(query_results)] # Excludes column header
    updateSelectizeInput(inputId = "selectedGene", choices = query_results, selected = "BRCA1")
  })

  # Update interval choices based on the selected predictor
  observe({
    if (input$predictor == 'REVEL') {
      interval_choices <- c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic')
    } else if (input$predictor == 'BayesDel') {
      interval_choices <- c("Moderate Benign", "Supporting Benign", "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic")
    } else {
      interval_choices <- character(0)
    }
    # Update the choices for the "intervals" input
    updatePickerInput(session, "intervals", choices = interval_choices)
  })

  # Reactively create a filtered data table based on user selections
  # Filtering logic based on selection of gene, selection of interval, and selection of whether gene is concordant, discordant, or indeterminate
  filteredTable <- reactive({
    if (input$predictor == 'REVEL') {
      my_dict <- list(
        'Strong Pathogenic' = 'pr_variant_pass_strong_path',
        'Moderate Pathogenic' = 'pr_variant_pass_moderate_path',
        'Supporting Pathogenic' = 'pr_variant_pass_supporting_path',
        'Very Strong Benign' = 'pr_variant_pass_very_strong_ben',
        'Strong Benign' = 'pr_variant_pass_strong_ben',
        'Moderate Benign' = 'pr_variant_pass_moderate_ben',
        'Supporting Benign' = 'pr_variant_pass_supporting_ben')
      data <- filter(datasetDictionary[[input$dataset]], Gene == input$selectedGene) # Filter for gene
      # Join DataFrames on Gene, Chr, and GeneID
      df <- read_csv(fileDictionary[[input$dataset]]$var1)
      data <- inner_join(data, df, by = c("Gene", "Chr", "GeneID"))
      # Define the REVEL intervals
      intervals <- list(
        'Very Strong Benign' = c(0, 0.003),
        'Strong Benign' = c(0.003, 0.016),
        'Moderate Benign' = c(0.016, 0.183),
        'Supporting Benign' = c(0.183, 0.290),
        "No Data" = c(0.290, 0.644),
        'Supporting Pathogenic' = c(0.644, 0.773),
        'Moderate Pathogenic' = c(0.773, 0.932),
        'Strong Pathogenic' = c(0.932, 1)
        )
        # Create a logical vector for each interval and combine them
        filter_conditions <- lapply(input$intervals, function(interval) {
        bounds <- intervals[[interval]]
        data[['REVEL']] >= bounds[1] & data[['REVEL']] < bounds[2]
        })
        combined_condition <- Reduce("|", filter_conditions)
        # Apply the combined filter to the dataframe
        data <- data[combined_condition, ]
        # Now proceed with the filtering as before
        column_filter_conditions <- lapply(my_dict[input$intervals], function(column) {
        sapply(input$categories, function(category) {
            grepl(category, data[[column]])
        }) %>% apply(1, any) # Check if any category is present in each row
        })
        combined_column_condition <- Reduce("|", column_filter_conditions)
        # Apply the combined filter to the dataframe
        data <- data[combined_column_condition, ]
    } else if (input$predictor == 'BayesDel') {
       my_dict <- list(
         'Strong Pathogenic' = 'pr_variant_pass_strong_path',
         'Moderate Pathogenic' = 'pr_variant_pass_moderate_path',
         'Supporting Pathogenic' = 'pr_variant_pass_supporting_path',
         'Moderate Benign' = 'pr_variant_pass_moderate_ben',
         'Supporting Benign' = 'pr_variant_pass_supporting_ben')
       data <- filter(datasetDictionary[[input$dataset]], Gene == input$selectedGene) # Filter for gene
       # Define the BayesDel intervals
       intervals <- list(
        'Moderate Benign' = c(-1.30, -0.36),
        'Supporting Benign' = c(-0.36, -0.18),
        "No Data" = c(-0.18, 0.13),
        'Supporting Pathogenic' = c(0.13, 0.27),
        'Moderate Pathogenic' = c(0.27, 0.5),
        'Strong Pathogenic' = c(0.5, 1)
       )
       # Join DataFrames on Gene, Chr, and GeneID
       df <- read_csv(fileDictionary[[input$dataset]]$var2)
       data <- inner_join(data, df, by = c("Gene", "Chr", "GeneID"))
       # Create a logical vector for each interval and combine them
       filter_conditions <- lapply(input$intervals, function(interval) {
       bounds <- intervals[[interval]]
       data[['BayesDel_noAF_score']] > bounds[1] & data[['BayesDel_noAF_score']] <= bounds[2]
       })
       combined_condition <- Reduce("|", filter_conditions)
       # Apply the combined filter to the dataframe
       data <- data[combined_condition, ]
       # Now proceed with the filtering as before
       column_filter_conditions <- lapply(my_dict[input$intervals], function(column) {
       sapply(input$categories, function(category) {
         grepl(category, data[[column]])
       }) %>% apply(1, any) # Check if any category is present in each row
       })
       combined_column_condition <- Reduce("|", column_filter_conditions)
       # Apply the combined filter to the dataframe
       data <- data[combined_column_condition, ]
    }
    return(data)
  })

  # Render the data table with options for scrolling
  output$data_table <- renderDataTable({
    datatable(filteredTable(), options = list(scrollX = TRUE))
  })

  # Render ValueBox for concordant genes
  output$concordant <- renderValueBox({
    rect_data <- interval_categories()
    data <- rect_data$data
    if (input$predictor == 'REVEL') {
        conditions = c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic')
    } else if (input$predictor == 'BayesDel') {
        conditions = c("Moderate Benign", "Supporting Benign", '', "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic")
    }
    filter <- conditions[data == "Concordant"]
    concordant <- as.character(filter[complete.cases(filter)])
    if (length(concordant) > 0) {
        result <- paste(concordant, collapse = ", ")
    } else {
        result <- "None"
    }
    valueBox(
        value = tags$p("Concordant", style = "font-size: 64%;"),
        subtitle =tags$p(result, style = "font-size: 82%;"),
        color = "green",
        icon=icon('check')
    )
  })

  # Render ValueBox for discordant genes
  output$discordant <- renderValueBox({
    rect_data <- interval_categories()
    data <- rect_data$data
    if (input$predictor == 'REVEL') {
        conditions = c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic')
    } else if (input$predictor == 'BayesDel') {
        conditions = c("Moderate Benign", "Supporting Benign", '', "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic")
    }
    filter <- conditions[data == "Discordant"]
    discordant <- as.character(filter[complete.cases(filter)])
    if (length(discordant) > 0) {
        result <- paste(discordant, collapse = ", ")
    } else {
        result <- "None"
    }
    valueBox(
        value = tags$p("Discordant", style = "font-size: 64%;"),
        subtitle =tags$p(result, style = "font-size: 82%;"),
        color = "maroon",
        icon=icon('minus')
    )
  })

  # Render ValueBox for indeterminate genes
  output$indeterminate <- renderValueBox({
    rect_data <- interval_categories()
    data <- rect_data$data
    if (input$predictor == 'REVEL') {
        conditions = c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic')
    } else if (input$predictor == 'BayesDel') {
        conditions = c("Moderate Benign", "Supporting Benign", '', "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic")
    }
    filter <- conditions[data == "Indeterminate"]
    indeterminate <- as.character(filter[complete.cases(filter)])
    if (length(indeterminate) > 0) {
        result <- paste(indeterminate, collapse = ", ")
    } else {
        result <- "None"
    }
    valueBox(
        value = tags$p("Indeterminate", style = "font-size: 64%;"),
        subtitle =tags$p(result, style = "font-size: 82%;"),
        color = "blue"
    )
  })
  
  # Reactively identify interval categories
  interval_categories <- reactive({
    if (input$predictor == 'REVEL') {
        data_subset <- read_csv(fileDictionary[[input$dataset]]$var1)
        data_subset <- filter(data_subset, Gene == input$selectedGene)
        
        v_strong_benign <- data_subset$pr_variant_pass_very_strong_ben
        strong_benign <- data_subset$pr_variant_pass_strong_ben
        moderate_benign <- data_subset$pr_variant_pass_moderate_ben
        supporting_benign <- data_subset$pr_variant_pass_supporting_ben
        strong_path <- data_subset$pr_variant_pass_strong_path
        moderate_path <- data_subset$pr_variant_pass_moderate_path
        supporting_path <- data_subset$pr_variant_pass_supporting_path

        v_strong_benign_controls <- data_subset$control_variants_very_strong_benign
        strong_benign_controls <- data_subset$control_variants_strong_benign
        moderate_benign_controls <- data_subset$control_variants_moderate_benign
        supporting_benign_controls <- data_subset$control_variants_supporting_benign
        strong_path_controls <- data_subset$control_variants_strong_pathogenic
        moderate_path_controls <- data_subset$control_variants_moderate_pathogenic
        supporting_path_controls <- data_subset$control_variants_supporting_pathogenic

        rect_data <- data.frame(
            data = c(v_strong_benign, strong_benign, moderate_benign, supporting_benign, "No Data", supporting_path, moderate_path, strong_path),
            categories = c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic'),
            x_left = c(0, 0.003, 0.016, 0.183, 0.290, 0.644, 0.773, 0.932), # Left side of the rectangles
            x_right = c(0.003, 0.016, 0.183, 0.290, 0.644, 0.773, 0.932, 1), # Right side of the rectangles
            interval = c("", "", "", "", "", "", "", ""), # Labels for intervals
            pass_fail = factor(c(v_strong_benign, strong_benign, moderate_benign, supporting_benign, "No Data", supporting_path, moderate_path, strong_path), levels = c("Concordant", "Discordant", "Indeterminate", "No Data")),
            number_control_variants = c(v_strong_benign_controls, strong_benign_controls, moderate_benign_controls, supporting_benign_controls, "No Data", supporting_path_controls, moderate_path_controls, strong_path_controls)
        )
    } else {
        data_subset <- read_csv(fileDictionary[[input$dataset]]$var2)
        data_subset <- filter(data_subset, Gene == input$selectedGene)
        
        moderate_benign <- data_subset$pr_variant_pass_moderate_ben
        supporting_benign <- data_subset$pr_variant_pass_supporting_ben
        strong_path <- data_subset$pr_variant_pass_strong_path
        moderate_path <- data_subset$pr_variant_pass_moderate_path
        supporting_path <- data_subset$pr_variant_pass_supporting_path

        moderate_benign_controls <- data_subset$control_variants_moderate_benign
        supporting_benign_controls <- data_subset$control_variants_supporting_benign
        strong_path_controls <- data_subset$control_variants_strong_pathogenic
        moderate_path_controls <- data_subset$control_variants_moderate_pathogenic
        supporting_path_controls <- data_subset$control_variants_supporting_pathogenic

        rect_data <- data.frame(
            data = c(moderate_benign, supporting_benign, "No Data", supporting_path, moderate_path, strong_path),
            categories = c('Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic'),
            x_left = c(-1.30, -0.36, -0.18, 0.13, 0.27, 0.5), # Left side of the rectangles
            x_right = c(-0.36, -0.18, 0.13, 0.27, 0.5, 0.76), # Right side of the rectangles
            interval = c("", "", "", "", "", ""), # Labels for intervals
            pass_fail = factor(c(moderate_benign, supporting_benign, "No Data", supporting_path, moderate_path, strong_path), levels = c("Concordant", "Discordant", "Indeterminate", "No Data")),
            number_control_variants = c(moderate_benign_controls, supporting_benign_controls, "No Data", supporting_path_controls, moderate_path_controls, strong_path_controls)
      )
    }
    return(rect_data)
  })

# Reactively identify interval categories
#   interval_categories <- reactive({
#     if (input$predictor == 'REVEL') {
#       data_subset <- read_csv(fileDictionary[[input$dataset]]$var1)
#       data_subset <- filter(data_subset, Gene == input$selectedGene)
#       v_strong_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_very_strong_benign']), 'variant_pass'][1, 1]
#       strong_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_strong_benign']), 'variant_pass'][1, 1]
#       moderate_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_moderate_benign']), 'variant_pass'][1, 1]
#       supporting_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_supporting_benign']), 'variant_pass'][1, 1]
#       strong_path <- data_subset[complete.cases(data_subset[, 'pass_fail_strong_pathogenic']), 'variant_pass'][1, 1]
#       moderate_path <- data_subset[complete.cases(data_subset[, 'pass_fail_moderate_pathogenic']), 'variant_pass'][1, 1]
#       supporting_path <- data_subset[complete.cases(data_subset[, 'pass_fail_supporting_pathogenic']), 'variant_pass'][1, 1]

#       v_strong_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_very_strong_benign']), 'control_variants_very_strong_benign'][1, 1]
#       strong_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_strong_benign']), 'control_variants_strong_benign'][1, 1]
#       moderate_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_benign']), 'control_variants_moderate_benign'][1, 1]
#       supporting_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_benign']), 'control_variants_supporting_benign'][1, 1]
#       strong_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_strong_pathogenic']), 'control_variants_strong_pathogenic'][1, 1]
#       moderate_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_pathogenic']), 'control_variants_moderate_pathogenic'][1, 1]
#       supporting_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_pathogenic']), 'control_variants_supporting_pathogenic'][1, 1]

#       v_strong_benign_percent <- data_subset[complete.cases(data_subset[, 'control_variants_very_strong_benign']), 'percent'][1, 1]
#       strong_benign_percent <- data_subset[complete.cases(data_subset[, 'control_variants_strong_benign']), 'percent'][1, 1]
#       moderate_benign_percent <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_benign']), 'percent'][1, 1]
#       supporting_benign_percent <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_benign']), 'percent'][1, 1]
#       strong_path_percent <- data_subset[complete.cases(data_subset[, 'control_variants_strong_pathogenic']), 'percent'][1, 1]
#       moderate_path_percent <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_pathogenic']), 'percent'][1, 1]
#       supporting_path_percent <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_pathogenic']), 'percent'][1, 1]

#       rect_data <- data.frame(
#         data = c(v_strong_benign[[1, 1]], strong_benign[[1, 1]], moderate_benign[[1, 1]], supporting_benign[[1, 1]], "No Data", supporting_path[[1, 1]], moderate_path[[1, 1]], strong_path[[1, 1]]),
#         categories = c('Very Strong Benign', 'Strong Benign', 'Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic'),
#         x_left = c(0, 0.003, 0.016, 0.183, 0.290, 0.644, 0.773, 0.932), # Left side of the rectangles
#         x_right = c(0.003, 0.016, 0.183, 0.290, 0.644, 0.773, 0.932, 1), # Right side of the rectangles
#         interval = c("", "", "", "", "", "", "", ""), # Labels for intervals
#         pass_fail = factor(c(v_strong_benign[[1, 1]], strong_benign[[1, 1]], moderate_benign[[1, 1]], supporting_benign[[1, 1]], "No Data", supporting_path[[1, 1]], moderate_path[[1, 1]], strong_path[[1, 1]]), levels = c("Concordant", "Discordant", "Indeterminate", "No Data")),
#         number_control_variants = c(v_strong_benign_controls[[1, 1]], strong_benign_controls[[1, 1]], moderate_benign_controls[[1, 1]], supporting_benign_controls[[1, 1]], "No Data", supporting_path_controls[[1, 1]], moderate_path_controls[[1, 1]], strong_path_controls[[1, 1]]),
#         percent_misclassification = c(v_strong_benign_percent[[1, 1]], strong_benign_percent[[1, 1]], moderate_benign_percent[[1, 1]], supporting_benign_percent[[1, 1]], "No Data", supporting_path_percent[[1, 1]], moderate_path_percent[[1, 1]], strong_path_percent[[1, 1]])
#       )

#     } else if (input$predictor == 'BayesDel') {
#       data_subset <- read_csv(fileDictionary[[input$dataset]]$var2)
#       data_subset <- filter(data_subset, Gene == input$selectedGene)
#       moderate_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_moderate_benign']), 'variant_pass'][1, 1]
#       supporting_benign <- data_subset[complete.cases(data_subset[, 'pass_fail_supporting_benign']), 'variant_pass'][1, 1]
#       strong_path <- data_subset[complete.cases(data_subset[, 'pass_fail_strong_pathogenic']), 'variant_pass'][1, 1]
#       moderate_path <- data_subset[complete.cases(data_subset[, 'pass_fail_moderate_pathogenic']), 'variant_pass'][1, 1]
#       supporting_path <- data_subset[complete.cases(data_subset[, 'pass_fail_supporting_pathogenic']), 'variant_pass'][1, 1]

#       moderate_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_benign']), 'control_variants_moderate_benign'][1, 1]
#       supporting_benign_controls <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_benign']), 'control_variants_supporting_benign'][1, 1]
#       strong_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_strong_pathogenic']), 'control_variants_strong_pathogenic'][1, 1]
#       moderate_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_moderate_pathogenic']), 'control_variants_moderate_pathogenic'][1, 1]
#       supporting_path_controls <- data_subset[complete.cases(data_subset[, 'control_variants_supporting_pathogenic']), 'control_variants_supporting_pathogenic'][1, 1]

#       rect_data <- data.frame(
#         data = c(moderate_benign[[1, 1]], supporting_benign[[1, 1]], "No Data", supporting_path[[1, 1]], moderate_path[[1, 1]], strong_path[[1, 1]]),
#         categories = c('Moderate Benign', 'Supporting Benign', '', 'Supporting Pathogenic', 'Moderate Pathogenic', 'Strong Pathogenic'),
#         x_left = c(-1.30, -0.36, -0.18, 0.13, 0.27, 0.5), # Left side of the rectangles
#         x_right = c(-0.36, -0.18, 0.13, 0.27, 0.5, 0.76), # Right side of the rectangles
#         interval = c("", "", "", "", "", ""), # Labels for intervals
#         pass_fail = factor(c(moderate_benign[[1, 1]], supporting_benign[[1, 1]], "No Data", supporting_path[[1, 1]], moderate_path[[1, 1]], strong_path[[1, 1]]), levels = c("Concordant", "Discordant", "Indeterminate", "No Data")),
#         number_control_variants = c(moderate_benign_controls[[1, 1]], supporting_benign_controls[[1, 1]], "No Data", supporting_path_controls[[1, 1]], moderate_path_controls[[1, 1]], strong_path_controls[[1, 1]])
#       )
#     } 
#     return(rect_data)
#   })

  # Whole Genome plot based on selected dataset
  wg_plot <- reactive({
    if (input$dataset == 'ClinVar 2023 Dataset') {
      data_subset_wg <- clinvar23
    } else if (input$dataset == 'One Star Missense') {
        data_subset_wg <- one_star_missense
    } else if (input$dataset == 'ClinVar 2023 Dataset without training variants') {
       data_subset_wg <- no_training_variants
    } else if (input$dataset == 'ClinGen SVI Calibration Dataset') {
        data_subset_wg <- clingen
    }
    if (input$predictor == 'REVEL') {
      phist = gghistogram(data_subset_wg, x = 'REVEL', fill = 'Significance', rug = FALSE,position = 'stack', palette = c('Benign' = 'royalblue3', 'Likely_benign' = 'steelblue1','Likely_pathogenic' = 'coral', 'Pathogenic' = 'firebrick3','Pathogenic/Likely_pathogenic' = 'coral', 'Benign/Likely_benign' = 'steelblue1'), alpha = 0.9, color = NA) + theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=14))+ geom_vline(xintercept = 0.932, linetype = 'twodash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.773, linetype = 'longdash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.644, linetype = 'dashed',color = '#990000',size = 0.5)+ geom_vline(xintercept = 0.290, linetype = 'dashed', color = '#0033CC',size = 0.5)+ geom_vline(xintercept = 0.183, linetype = 'longdash',color = '#0033CC',size = 0.5)+ geom_vline(xintercept = 0.016, linetype = 'twodash',color = '#0033CC',size = 0.5) + geom_vline(xintercept = 0.003, linetype = 'solid',color = '#0033CC',size = 0.5)+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(0,1))+ ggtitle('REVEL scores') +
      theme(
        axis.line = element_line(size = 0.2),  # Adjust the axis line thickness
        axis.text = element_text(size = 9)  # Adjust the axis text size
      ) + xlab("REVEL Scores") + ylab("Variant Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0)))+scale_x_continuous(expand = c(0,0), limits = c(0,1))
      phist_plotly <- ggplotly(phist, tooltip = FALSE)
    } else {
      phist = gghistogram(data_subset_wg, x = 'BayesDel_noAF_score', fill = 'Significance', rug = FALSE,position = 'stack', palette = c('Benign' = 'royalblue3', 'Likely_benign' = 'steelblue1','Likely_pathogenic' = 'coral', 'Pathogenic' = 'firebrick3','Pathogenic/Likely_pathogenic' = 'coral', 'Benign/Likely_benign' = 'steelblue1'), alpha = 0.9, color = NA, bins = 30) + theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=14))+ geom_vline(xintercept = 0.5, linetype = 'twodash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.27, linetype = 'longdash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.13, linetype = 'dashed',color = '#990000',size = 0.5)+ geom_vline(xintercept = -0.36, linetype = 'dashed', color = '#0033CC',size = 0.5)+ geom_vline(xintercept = -0.18, linetype = 'longdash',color = '#0033CC',size = 0.5)+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(-1.30,0.76)) +
      theme(
        axis.line = element_line(size = 0.2),  # Adjust the axis line thickness
        axis.text = element_text(size = 9)  # Adjust the axis text size
      ) + ggtitle('BayesDel scores') + xlab("BayesDel Scores") + ylab("Variant Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0)))+scale_x_continuous(expand = c(0,0))
      phist_plotly <- ggplotly(phist, tooltip = FALSE)
    }
    phist_plotly
  })

  # Render the whole genome distribution plot
  output$render_wg_plot <- renderPlotly({
    wg_plot()
  })

  # Gene-wise distribution plot based on selected dataset
  gene_wise_plot <- reactive({
    if (input$dataset == 'One Star Missense') {
        data_subset <- one_star_missense
    } else if (input$dataset == 'ClinVar 2023 Dataset without training variants') {
       data_subset <- no_training_variants
    } else if (input$dataset == 'ClinGen SVI Calibration Dataset') {
        data_subset <- clingen
    } else if (input$dataset == 'ClinVar 2023 Dataset') {
        data_subset <- clinvar23
    }
    data_subset <- filter(data_subset,data_subset$Gene == input$selectedGene)
    rect_data <- interval_categories()
    data <- rect_data$data
    number_control_variants <- rect_data$number_control_variants
    # percent_misclassification <- rect_data$percent_misclassification # This can be uncommented once percent_misclassification is added to analysis dataframe
    if (input$predictor == 'REVEL') {
      phist <- gghistogram(data_subset, x = 'REVEL', fill = 'Significance', rug = FALSE,position = 'stack', palette = c('Benign' = 'royalblue3', 'Likely_benign' = 'steelblue1','Likely_pathogenic' = 'coral', 'Pathogenic' = 'firebrick3','Pathogenic/Likely_pathogenic' = 'coral', 'Benign/Likely_benign' = 'steelblue1'), alpha = 0.9, color = NA) + theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=14))+ geom_vline(xintercept = 0.932, linetype = 'twodash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.773, linetype = 'longdash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.644, linetype = 'dashed',color = '#990000',size = 0.5)+ geom_vline(xintercept = 0.290, linetype = 'dashed', color = '#0033CC',size = 0.5)+ geom_vline(xintercept = 0.183, linetype = 'longdash',color = '#0033CC',size = 0.5)+ geom_vline(xintercept = 0.016, linetype = 'twodash',color = '#0033CC',size = 0.5) + geom_vline(xintercept = 0.003, linetype = 'solid',color = '#0033CC',size = 0.5)+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(0,1)) + ggtitle('REVEL scores', input$selectedGene) +
      theme(
        axis.line = element_line(size = 0.2),  # Adjust the axis line thickness
        axis.text = element_text(size = 9)  # Adjust the axis text size
      ) + xlab("REVEL Scores") + ylab("Variant Count") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      # Retrieve the data for the density plot
      density_data <- filtered_revel_data()  # Replace with the correct call to retrieve REVEL scores
      # Add density plot to the existing plot
      pdensity <- ggdensity(density_data, x = 'REVEL', fill = '#999999', color = NA) +
      scale_y_continuous(expand = expansion(mult = c(0, 0)))+scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
      theme(
        axis.line = element_line(size = 0.2),  # Adjust the axis line thickness
        axis.text = element_text(size = 9)  # Adjust the axis text size
      ) + xlab("REVEL Scores")
      phist <- ggplotly(phist, tooltip = NULL)
      pdensity <- ggplotly(pdensity, tooltip = NULL)
      # Create a combined interactive plot with dual y-axes
      phist <- subplot(phist, pdensity, shareX = TRUE, nrows=2, heights=c(0.01, 0.99))
      combined_label <- paste0('Pass/Fail: ', rect_data$pass_fail, '<br>', 'Category: ', rect_data$categories, '<br>', 'Number Control Variants: ', number_control_variants)

      # Create the rectangles and text annotations
      ann <- ggplot(rect_data, aes(xmin = x_left, xmax = x_right)) +
        geom_rect(aes(ymin = 1.2, ymax = 1.4, fill = pass_fail, text = combined_label), alpha = 0.5) +
        geom_text(aes(x = (x_left + x_right) / 2, y = 1.3, label = interval), vjust = 0, fontface = "bold") +
        scale_fill_manual(values = c("Concordant" = "#808000", "Discordant" = "#E8A11C", "Indeterminate" = "#695B73", "No Data" = "white"), na.value = "#695B73", guide=FALSE) +
        theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=9), axis.title=element_blank())+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(0,1))+
        theme(
          axis.line = element_blank(),  # Remove axis lines
          axis.ticks = element_blank(), # Remove axis ticks
          axis.text = element_blank(),
        )
      # Convert 'ann' to an interactive plot with hovers
      interactive_ann <- ggplotly(ann, tooltip = c("text")) %>%
        layout(
          hoverdistance = 2 # Increase the hoverdistance for less sensitivity
        )

      # Create a plotly subplot with 'interactive_ann' and 'phist'
      subplot <- subplot(
        interactive_ann,
        phist,
        shareX = TRUE, shareY = TRUE,
        heights = c(0.1, 0.9),
        nrows = 2
      )

      subplot %>%
        layout(
            yaxis2 = list(side = "right", title = "All SNV Density"),
            yaxis3 = list(overlaying = "y2", side = "left", title = "Variant Count")
        )

    } else {
        phist = gghistogram(data_subset, x = 'BayesDel_noAF_score', fill = 'Significance', rug = FALSE,position = 'stack', palette = c('Benign' = 'royalblue3', 'Likely_benign' = 'steelblue1','Likely_pathogenic' = 'coral', 'Pathogenic' = 'firebrick3','Pathogenic/Likely_pathogenic' = 'coral', 'Benign/Likely_benign' = 'steelblue1'), alpha = 0.9, color = NA, bins = 30) + theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=14))+ geom_vline(xintercept = 0.5, linetype = 'twodash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.27, linetype = 'longdash', color = '#990000',size = 0.5) + geom_vline(xintercept = 0.13, linetype = 'dashed',color = '#990000',size = 0.5)+ geom_vline(xintercept = -0.36, linetype = 'dashed', color = '#0033CC',size = 0.5)+ geom_vline(xintercept = -0.18, linetype = 'longdash',color = '#0033CC',size = 0.5)+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(-1.30,0.76)) + ggtitle('BayesDel scores', input$selectedGene) +
        theme(
          axis.line = element_line(size = 0.2),  # Adjust the axis line thickness
          axis.text = element_text(size = 9)  # Adjust the axis text size
        ) + xlab("BayesDel Scores") + ylab("Variant Count") + scale_y_continuous(expand = expansion(mult = c(0, 0)))
        phist <- ggplotly(phist, tooltip = NULL)

        # Retrieve the data for the density plot
        density_data <- filtered_bd_data()
        # Add density plot to the existing plot
        pdensity <- ggdensity(density_data, x = 'BayesDel_nsfp33a_noAF', fill = '#999999', color = NA) + xlab('BayesDel Scores') +
        scale_y_continuous(expand = expansion(mult = c(0, 0)))+scale_x_continuous(expand = c(0,0), limits = c(-1.30,0.76))
        pdensity <- ggplotly(pdensity, tooltip = NULL)
        # Create a combined interactive plot with dual y-axes
        phist <- subplot(phist, pdensity, shareX = TRUE, nrows=2, heights=c(0.01, 0.99))
        combined_label <- paste0('Pass/Fail: ', rect_data$pass_fail, '<br>', 'Category: ', rect_data$categories, '<br>', 'Number Control Variants: ', number_control_variants)

        # Create the rectangles and text annotations
        ann <- ggplot(rect_data, aes(xmin = x_left, xmax = x_right)) +
          geom_rect(aes(ymin = 1.2, ymax = 1.4, fill = pass_fail, text = combined_label), alpha = 0.5) +
          geom_text(aes(x = (x_left + x_right) / 2, y = 1.3, label = interval), vjust = 0, fontface = "bold") +
          scale_fill_manual(values = c("Concordant" = "#808000", "Discordant" = "#E8A11C", "Indeterminate" = "#695B73", "No Data" = "white"), na.value = "#695B73", guide=FALSE) +
          theme_half_open(11, rel_small = 1) + theme(legend.position = "none", axis.text=element_text(size=14), axis.title=element_blank())+xlim(c(0,1))+scale_x_continuous(expand = c(0,0), limits = c(0,1))+
          theme(
            axis.line = element_blank(),  # Remove axis lines
            axis.ticks = element_blank(), # Remove axis ticks
            axis.text = element_blank(),
          ) + scale_x_continuous(expand = c(0,0), limits = c(-1.30,0.76))

        # Convert 'ann' to an interactive plot with hovers
        interactive_ann <- ggplotly(ann, tooltip = c("text")) %>%
        layout(
          hoverdistance = 2 # Increase the hoverdistance for less sensitivity
        )

        # Create a plotly subplot with 'interactive_ann' and 'phist'
        subplot <- subplot(
            interactive_ann,
            phist,
            shareX = TRUE, shareY = TRUE,
            heights = c(0.1, 0.9),
            nrows = 2
        )
        subplot %>%
            layout(
                yaxis2 = list(side = "right", title = "All SNV Density"),
                yaxis3 = list(overlaying = "y2", side = "left", title = "Variant Count")
            )
    }
  })

  # Render the gene-wise distribution plot
  output$render_gene_wise_plot <- renderPlotly({
    gene_wise_plot() %>% layout(dragmode = "zoom", dragmode = FALSE)
  })

  # Reactively create a summary table for filters
  filter_summary_table <- reactive({
    filter_interval_categories <- interval_categories()
    if (input$predictor == 'REVEL') {
      filter_interval_categories <- filter_interval_categories[, c("categories", "pass_fail", "number_control_variants")]
    } else {
      filter_interval_categories <- filter_interval_categories[, c("categories", "pass_fail", "number_control_variants")]
    }
    return(filter_interval_categories)
  })

  # Download handlers for various tables
  output$downloadTable <- downloadHandler(
    filename = function() { "filtered_table.csv" },
    content = function(fname) {
      write.csv(filteredTable(), fname)
    }
  )

  output$downloadSummaryTable <- downloadHandler(
    filename = function() { "gene_wise_summary_table.csv" },
    content = function(fname) {
      write.csv(filter_summary_table(), fname)
    }
  )

  output$downloadWGSummaryTable <- downloadHandler(
    filename = function() { "WG_summary_table.csv" },
    content = function(fname) {
      write.csv(filter_summary_table(), fname)
    }
  )
}

# Run the app
shinyApp(ui, server)
