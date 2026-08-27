
generate_intervals_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h2("Generate Species Intervals From Trees"),
    
    sidebarLayout(
      sidebarPanel(
        fileInput(
          ns("tree_file"),
          label = span(
            "Upload timetree",
            tooltip(
              bsicons::bs_icon(
                "info-circle",
                title = "About tree format"
              ),
              paste(
                "Upload a rooted, time-scaled Newick tree.",
                "Tip labels must contain binomial species names (eg. 
                Homo_sapiens).",
                "If there are duplicates, add a unique identifier (eg. 
                Homo_sapiens_1)",
                "Unclassied species should be Genus_sp._unique_identifier",
                "Accepted formats: .nwk, .tree, .tre, and .treefile."
              ),
              placement = "right"
            )
          ),
          accept = c(".nwk", ".tree", ".tre", ".treefile")
        ),
        
        numericInput(
          ns("monophyly_threshold"),
          "Minimum monophyly (%)",
          value = 95,
          min = 0,
          max = 100,
          step = 1
        ),
        
        radioButtons(
          ns("unclassified_handling"),
          "Unclassified tips",
          choices = c(
            "Ignore" = "ignore",
            "Count as non-monophyly" = "include"
          ),
          selected = "ignore"
        ),
        
        selectInput(
          ns("merge_species_1"),
          "Keep species name",
          choices = NULL
        ),
        
        selectInput(
          ns("merge_species_2"),
          "Merge this species into it",
          choices = NULL
        ),
        
        actionButton(
          ns("merge_species"),
          "Merge species"
        ),
        
        div(
          class = "d-grid gap-2 mt-3",
          
          actionButton(
            ns("add_tree"),
            "Add tree intervals"
          ),
          
          downloadButton(
            ns("download_intervals"),
            "Download generated intervals"
          )
        )
        ),
      
      mainPanel(
        accordion(
          
          accordion_panel(
            "Tree figure",
            plotOutput(
              ns("tree_plot"),
              height = "700px"
            )
          ),
          
          accordion_panel(
            "Tree summary",
            tableOutput(ns("tree_summary"))
          ),
          
          accordion_panel(
            "Generated intervals",
            tableOutput(ns("interval_preview"))
          ),
          
          accordion_panel(
            "Recognized tip names",
            tableOutput(ns("tip_species_preview"))
          ),
          
          accordion_panel(
            "Species groups",
            tableOutput(ns("species_groups_preview"))
          ),
          
          open = "Tree summary",
          multiple = TRUE
        )
      )
    )
  )
}


generate_intervals_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    saved_intervals <- reactiveVal(tibble())
    
    merge_pair <- reactiveVal(NULL)
    
    tree <- reactive({
      req(input$tree_file)
      
      ape::read.tree(
        input$tree_file$datapath
      )
    })
    
    get_species_mrca <- function(current_tree, species_tips) {
      
      tip_numbers <- which(
        current_tree$tip.label %in% species_tips
      )
      
      if (length(tip_numbers) != length(species_tips)) {
        stop("Not all species tips were matched in the tree.")
      }
      
      ape::getMRCA(current_tree, tip_numbers)
    }
    
    extract_species_name <- function(tip_label) {
      
      cleaned_label <- tip_label |>
        str_replace_all(" ", "_")
      
      match <- str_match(
        cleaned_label,
        "^([^_]+)_(sp\\.?|[^_]+)(?:_|$)"
      )
      
      species <- if_else(
        is.na(match[, 1]),
        NA_character_,
        paste0(
          match[, 2],
          "_",
          if_else(match[, 3] == "sp", "sp.", match[, 3])
        )
      )
      
      species
    }
    
    
    tip_species <- reactive({
      current_tree <- tree()
      
      tibble(
        tip_index = seq_along(current_tree$tip.label),
        tip_label = current_tree$tip.label,
        species = extract_species_name(current_tree$tip.label)
      ) |>
        mutate(
          status = case_when(
            is.na(species) ~ "Unrecognized",
            str_detect(species, "_sp\\.$") ~ "Unclassified",
            TRUE ~ "Recognized"
          )
        )
    })
    
    observe({
      species_choices <- tip_species() |>
        filter(!is.na(species)) |>
        distinct(species) |>
        arrange(species) |>
        pull(species)
      
      req(length(species_choices) > 0)
      
      updateSelectInput(
        session,
        "merge_species_1",
        choices = species_choices
      )
      
      updateSelectInput(
        session,
        "merge_species_2",
        choices = species_choices
      )
    })
    
    
    observeEvent(input$merge_species, {
      req(input$merge_species_1, input$merge_species_2)
      
      if (input$merge_species_1 == input$merge_species_2) {
        showNotification(
          "Please select two different species.",
          type = "error"
        )
        return()
      }
      
      merge_pair(
        c(
          input$merge_species_1,
          input$merge_species_2
        )
      )
      
      showNotification(
        paste(
          input$merge_species_2,
          "merged into",
          input$merge_species_1
        ),
        type = "message"
      )
    })
    
    
    merged_tip_species <- reactive({
      tip_data <- tip_species()
      selected_pair <- merge_pair()
      
      if (is.null(selected_pair)) {
        return(tip_data)
      }
      
      tip_data |>
        mutate(
          species = if_else(
            species %in% selected_pair,
            selected_pair[[1]],
            species
          )
        )
    })
    
    species_groups <- reactive({
      current_tree <- tree()
      tip_data <- merged_tip_species()
      
      groups <- tip_data |>
        filter(!is.na(species)) |>
        group_by(species) |>
        summarise(
          n_tips = n(),
          tips = list(tip_label),
          tip_indices = list(tip_index),
          all_unclassified = all(status == "Unclassified"),
          .groups = "drop"
        ) |>
        mutate(
          crown_node = map_int(
            tip_indices,
            ~ if (length(.x) > 1) {
              ape::getMRCA(current_tree, .x)
            } else {
              NA_integer_
            }
          ),
          
          monophyly_percent = pmap_dbl(
            list(species, tips, crown_node, all_unclassified),
            function(species_name, species_tips,
                     crown_node, unclassified) {
              
              if (unclassified) return(NA_real_)
              if (length(species_tips) == 1) return(100)
              
              descendant_tips <- ape::extract.clade(
                current_tree,
                crown_node
              )$tip.label
              
              descendant_data <- tip_data |>
                filter(tip_label %in% descendant_tips)
              
              if (input$unclassified_handling == "ignore") {
                descendant_data <- descendant_data |>
                  filter(status != "Unclassified")
              }
              
              100 * sum(
                descendant_data$species == species_name,
                na.rm = TRUE
              ) / nrow(descendant_data)
            }
          ),
          
          status = case_when(
            all_unclassified ~ "Unclassified",
            n_tips == 1 ~ "Singleton",
            monophyly_percent >= input$monophyly_threshold ~
              "Meets threshold",
            TRUE ~ "Below threshold"
          )
        )
      groups
    })
    
    
    output$tree_summary <- renderTable({
      current_tree <- tree()
      
      tibble(
        `Tree property` = c(
          "File",
          "Number of tips",
          "Number of internal nodes",
          "Rooted",
          "Has branch lengths"
        ),
        `Value` = c(
          input$tree_file$name,
          ape::Ntip(current_tree),
          ape::Nnode(current_tree),
          ape::is.rooted(current_tree),
          !is.null(current_tree$edge.length)
        )
      )
    },
    rownames = FALSE)
    
    
    intervals <- reactive({
      current_tree <- tree()
      
      validate(
        need(
          ape::is.ultrametric(current_tree, tol = 1e-4),
          "The uploaded tree must be ultrametric."
        )
      )
      
      depths <- ape::node.depth.edgelength(current_tree)
      
      tree_height <- max(
        depths[seq_len(ape::Ntip(current_tree))]
      )
      
      node_age <- tree_height - depths
      
      get_parent <- function(node) {
        parent <- current_tree$edge[
          current_tree$edge[, 2] == node,
          1
        ]
        
        if (length(parent) == 0) {
          return(NA_integer_)
        }
        
        parent[1]
      }
      
      species_groups() |>
        filter(status == "Meets threshold") |>
        mutate(
          stem_node = map_int(crown_node, get_parent),
          crown_age = node_age[crown_node],
          stem_age = node_age[stem_node],
          crown_age = if_else(
            abs(crown_age) < 1e-10,
            0,
            crown_age
          ),
          stem_age = if_else(
            abs(stem_age) < 1e-10,
            0,
            stem_age
          ),
          speciation_interval = stem_age - crown_age
        ) |>
        filter(!is.na(stem_age)) |>
        select(
          species,
          crown_age,
          stem_age,
          speciation_interval,
          n_tips
        )
    })
    
    observeEvent(input$add_tree, {
      
      message("STEP 1: Add tree button clicked")
      
      showNotification(
        "STEP 1: Button clicked",
        type = "message",
        duration = 5
      )
      
      req(input$tree_file)
      
      current_tree <- tree()
      
      ultrametric_result <- ape::is.ultrametric(
        current_tree,
        tol = 1e-4
      )
      
      message(
        "STEP 2: Ultrametric = ",
        ultrametric_result
      )
      
      showNotification(
        paste("STEP 2: Ultrametric =", ultrametric_result),
        type = if (ultrametric_result) "message" else "error",
        duration = 5
      )
      
      if (!ultrametric_result) {
        return()
      }
      
      new_intervals <- tryCatch(
        {
          intervals()
        },
        error = function(e) {
          
          message(
            "STEP 3 ERROR: ",
            conditionMessage(e)
          )
          
          showNotification(
            paste(
              "STEP 3: intervals() failed:",
              conditionMessage(e)
            ),
            type = "error",
            duration = NULL
          )
          
          NULL
        }
      )
      
      if (is.null(new_intervals)) {
        return()
      }
      
      message(
        "STEP 3: Calculated interval rows = ",
        nrow(new_intervals)
      )
      
      showNotification(
        paste(
          "STEP 3: Calculated",
          nrow(new_intervals),
          "intervals"
        ),
        type = if (nrow(new_intervals) > 0) "message" else "warning",
        duration = 5
      )
      
      if (nrow(new_intervals) == 0) {
        return()
      }
      
      new_intervals <- new_intervals |>
        mutate(
          tree_file = input$tree_file$name,
          .before = 1
        )
      
      saved_intervals(
        bind_rows(
          saved_intervals(),
          new_intervals
        )
      )
      
      message(
        "STEP 4: Total saved rows = ",
        nrow(saved_intervals())
      )
      
      showNotification(
        paste(
          "STEP 4:",
          nrow(new_intervals),
          "intervals added successfully"
        ),
        type = "message",
        duration = 8
      )
    })
    
    
    output$interval_preview <- renderTable({
      validate(
        need(
          nrow(saved_intervals()) > 0,
          "Upload a tree and click 'Add tree intervals'."
        )
      )
      
      saved_intervals()
    })
    
    output$tip_species_preview <- renderTable({
      tip_species()
    },
    rownames = FALSE)
    
    output$species_groups_preview <- renderTable({
      species_groups() |>
        mutate(
          monophyly_percent = round(monophyly_percent, 1)
        ) |>
        select(
          species,
          n_tips,
          monophyly_percent,
          status
        )
    },
    rownames = FALSE)
    
    output$tree_plot <- renderPlot({
      current_tree <- tree()
      tip_data <- merged_tip_species()
      
      plot_tree <- current_tree
      
      plot_tree$tip.label <- if_else(
        is.na(tip_data$species),
        tip_data$tip_label,
        tip_data$species
      )
      
      get_parent <- function(node) {
        parent <- current_tree$edge[
          current_tree$edge[, 2] == node,
          1
        ]
        
        if (length(parent) == 0) {
          return(NA_integer_)
        }
        
        parent[1]
      }
      
      eligible_groups <- species_groups() |>
        filter(status == "Meets threshold") |>
        mutate(
          stem_node = map_int(crown_node, get_parent)
        )
      
      crown_nodes <- eligible_groups$crown_node |>
        unique() |>
        na.omit() |>
        as.integer()
      
      stem_nodes <- eligible_groups$stem_node |>
        unique() |>
        na.omit() |>
        as.integer()
      
      par(mar = c(1, 1, 1, 1))
      
      plot(
        plot_tree,
        use.edge.length = FALSE,
        direction = "rightwards",
        show.tip.label = TRUE,
        cex = 0.6,
        no.margin = TRUE
      )
      
      depths <- ape::node.depth.edgelength(current_tree)
      
      tree_height <- max(
        depths[seq_len(ape::Ntip(current_tree))]
      )
      
      node_age <- tree_height - depths
      
      crown_labels <- paste0(
        "C\n",
        round(node_age[crown_nodes], 2),
        " MY"
      )
      
      stem_labels <- paste0(
        "S\n",
        round(node_age[stem_nodes], 2),
        " MY"
      )
      
      if (length(crown_nodes) > 0) {
        ape::nodelabels(
          text = crown_labels,
          node = crown_nodes,
          frame = "circle",
          bg = "white",
          cex = 0.6
        )
      }
      
      if (length(stem_nodes) > 0) {
        ape::nodelabels(
          text = stem_labels,
          node = stem_nodes,
          frame = "rect",
          bg = "white",
          cex = 0.6
        )
      }
      
      legend(
        "topleft",
        legend = c("C = crown node", "S = stem node"),
        bty = "n",
        cex = 0.8
      )
    })
    
    
    output$download_intervals <- downloadHandler(
      filename = function() {
        "combined_tree_intervals.csv"
      },
      content = function(file) {
        req(nrow(saved_intervals()) > 0)
        
        write_csv(
          saved_intervals(),
          file
        )
      }
    )
    
    return(saved_intervals)
  })
}

