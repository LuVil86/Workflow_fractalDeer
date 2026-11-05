library(shiny)
library(bslib)
library(ggplot2)
library(gridExtra)
library(fs)
library(shinyFiles) # Package requis pour l'exploration de dossiers locaux

# Source des fonctions existantes
source("./getNSDValues.R")
source("./refineClusteringByNSD.R")

# Interface Utilisateur (UI)
ui <- page_sidebar(
  title = "Analyse et Raffinement NSD Multi-fichiers",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  sidebar = sidebar(
    width=400,
    title = "Configuration",
    
    # 1. Bouton d'exploration de dossier natif
    shinyDirButton("folder_btn", " Sélectionner un dossier", 
                   title = "Choisir le dossier contenant vos fichiers .csv", 
                   class = "btn-secondary w-100 mb-2"),
    
    # Text output discret pour afficher le dossier actuellement sélectionné
    verbatimTextOutput("current_folder", placeholder = TRUE),
    
    # 2. Choix dynamique du fichier (mis à jour selon le dossier sélectionné)
    selectInput("selected_stepSize", "Choisir la stepSize à analyser :", 
                choices = NULL, width = "100%"),
    
    # 3. Sélection du paramètre de raffinement
    selectInput("parameter", "Paramètre de raffinement :",
                choices = c("Ratio Mean NSD" = "ratio_meanNSD",
                            "Ratio End NSD" = "ratio_endNSD",
                            "Ratio Cumulative NSD" = "ratio_cumulativeNSD")),
    
    hr(),
    actionButton("run_analysis", "Lancer l'analyse", class = "btn-primary w-100"),
    
    # Section de téléchargement conditionnelle
    uiOutput("download_ui")
  ),
  
  # Zone principale avec les trois onglets
  navset_card_pill(
    nav_panel(
      title = "1. NSD Initial (getNSDValues)",
      card(
        card_header("Graphique du comportement initial"),
        plotOutput("plot_initial", height = "600px")
      )
    ),
    nav_panel(
      title = "2. Raffinement (refineClusteringByNSD)",
      card(
        card_header("Graphique de comparaison après raffinement"),
        plotOutput("plot_refined", height = "700px")
      )
    ),
    nav_panel(
      title = "3. NSD Final (getNSDValues)",
      card(
        card_header("Graphique du comportement final (après raffinement)"),
        plotOutput("plot_final", height = "600px")
      )
    )
  )
)

# Logique Serveur (Server)
server <- function(input, output, session) {
  
  # Stockage du chemin du fichier raffiné
  refined_file_path <- reactiveVal(NULL)
  
  # Initialisation de la racine du système de fichiers (Donne accès à tout le système local)
  volumes <- c(Home = fs::path_home(), Root = "/", `Dossier Actuel` = getwd())
  shinyDirChoose(input, "folder_btn", roots = volumes, session = session)
  
  # Conversion réactive du choix de l'explorateur en chemin absolu textuel
  selected_dir_path <- reactive({
    req(input$folder_btn)
    # Permet de reconstruire le chemin complet depuis la racine choisie
    path <- parseDirPath(volumes, input$folder_btn)
    if (length(path) == 0) return(NULL)
    return(path)
  })
  
  # Afficher le chemin du dossier sélectionné sous le bouton
  output$current_folder <- renderText({
    path <- selected_dir_path()
    if (is.null(path)) "Aucun dossier sélectionné" else path
  })
  
  # Surveillance du dossier sélectionné pour mettre à jour la liste des fichiers CSV
  observe({
    dir_path <- selected_dir_path()
    req(dir_path)
    
    if (dir_exists(dir_path)) {
      # Lister uniquement les fichiers .csv valides (exclut les "REFINED")
      all_files <- read_csv(list.files(dir_path, pattern = "results_Chisq_meanNSD", full.names = TRUE))
      
      stepSizes <- sort(all_files[,"stepSize"]$stepSize)
      print(stepSizes)
      # Mettre à jour le sélecteur de fichiers
      updateSelectInput(session, "selected_stepSize", choices = stepSizes)
    } else {
      updateSelectInput(session, "selected_stepSize", choices = character(0))
    }
  })

  # Déterminer le chemin complet du fichier sélectionné
  target_file <- reactive({
    dir_path <- selected_dir_path()
    all_files <- read_csv(list.files(dir_path, pattern = "results_Chisq_meanNSD", full.names = TRUE))
    fileFullName <- all_files[which(all_files$stepSize==input$selected_stepSize), "fullPath"]
    fileInputName <-str_split_i(fileFullName, pattern = "/", i = -1)
    req(dir_path, fileInputName)
    if (fileInputName == "") return(NULL)
    return(file.path(dir_path, fileInputName))
  })
  
  # Déclencher les calculs au clic sur le bouton
  data_analysis <- eventReactive(input$run_analysis, {
    file_path <- target_file()
    dir_path <- selected_dir_path()
    if (is.null(file_path) || !file.exists(file_path)) {
      showNotification("Fichier introuvable ou invalide.", type = "error")
      return(NULL)
    }
    
    showModal(modalDialog("Calculs et génération des graphiques en cours...", footer = NULL))
    on.exit(removeModal())
    
    # --- ÉTAPE 1 : Premier graphique ---
    png(tf1 <- tempfile(fileext = ".png"))
    dev.control("enable")
    getNSDValues(file_path, display = "behaviour", save.plot = FALSE, mutate.df = TRUE, show.breakpoints = TRUE)
    p_init <- recordPlot()
    dev.off()
    unlink(tf1)
    
    # --- ÉTAPE 2 : Exécution du raffinement ---
    png(tf2 <- tempfile(fileext = ".png"))
    dev.control("enable")
    refineClusteringByNSD(file_path,
                          parameter = input$parameter,
                          display.plot = TRUE,
                          save.plot = FALSE,
                          save.data.frame = TRUE, 
                          show.breakpoints = TRUE,
                          save.BIC = FALSE)
    p_ref <- recordPlot()
    dev.off()
    unlink(tf2)
    
    # --- ÉTAPE 3 : Reconstruction du chemin du fichier raffiné ---
    plotTitle<-strsplit(file_path, split=".csv")[[1]][1]
    plotTitle<-gsub("_"," ",tail(strsplit(plotTitle, split="/")[[1]],n=1))
    
    generated_csv <- file.path(dir_path, paste0(plotTitle, " REFINED ", toupper(input$parameter), ".csv"))
    refined_file_path(generated_csv)
    
    # --- ÉTAPE 4 : Troisième graphique ---
    p_final <- NULL
    if (file.exists(generated_csv)) {
      png(tf3 <- tempfile(fileext = ".png"))
      dev.control("enable")
      getNSDValues(generated_csv, display = "behaviour", save.plot = FALSE, mutate.df = FALSE, show.breakpoints = TRUE)
      p_final <- recordPlot()
      dev.off()
      unlink(tf3)
    } else {
      showNotification("Le fichier raffiné n'a pas pu être localisé.", type = "warning")
    }
    
    list(plot_initial = p_init, plot_refined = p_ref, plot_final = p_final)
  })
  
  # Génération dynamique des boutons de téléchargement
  output$download_ui <- renderUI({
    res <- data_analysis()
    if (is.null(res)) return(NULL)
    
    tagList(
      hr(),
      h5("Téléchargements"),
      downloadButton("download_csv", "Télécharger le CSV final", class = "btn-outline-success w-100 mb-2"),
      downloadButton("download_plot", "Télécharger le Plot Final (PNG)", class = "btn-outline-info w-100")
    )
  })
  
  # Téléchargement du fichier CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      outStepSize <- input$selected_stepSize
      
      paste0("REFINED_", toupper(input$parameter), "_",outStepSize )
    },
    content = function(file) {
      req(selected_dir_path())
      if (file.exists(selected_dir_path())) {
        file.copy(selected_dir_path(), file)
      } else {
        showNotification("Erreur : Fichier introuvable.", type = "error")
      }
    }
  )
  
  # Téléchargement du Graphique Final
  output$download_plot <- downloadHandler(
    filename = function() {
      plot_name <- strsplit(input$selected_file, split = "\\.")[]
      paste0("PLOT_FINAL_", toupper(input$parameter), "_", plot_name, ".png")
    },
    content = function(file) {
      res <- data_analysis()
      req(res$plot_final)
      
      png(file, width = 800, height = 600, res = 100)
      replayPlot(res$plot_final)
      dev.off()
    }
  )
  
  # Rendu graphique des trois onglets
  output$plot_initial <- renderPlot({
    res <- data_analysis()
    if (is.null(res)) return(NULL)
    replayPlot(res$plot_initial)
  })
  
  output$plot_refined <- renderPlot({
    res <- data_analysis()
    if (is.null(res)) return(NULL)
    replayPlot(res$plot_refined)
  })
  
  output$plot_final <- renderPlot({
    res <- data_analysis()
    if (is.null(res) || is.null(res$plot_final)) return(NULL)
    replayPlot(res$plot_final)
  })
}

# Lancement de l'application
shinyApp(ui = ui, server = server)