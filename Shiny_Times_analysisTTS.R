library(shiny)
library(tidyverse)
library(bslib)
library(multimode)
library(genefilter)
library(shinycssloaders)

#---Helper HRM Python files-BootmodeV3------------------------------------------

run_hrm_python <- function(
    values,
    outer_iterations,
    inner_iterations
) {
  
  # Temporary data file: one number per row
  data_file <- tempfile(fileext = ".txt")
  
  writeLines(
    as.character(na.omit(values)),
    data_file
  )
  
  data_file <- normalizePath(
    data_file,
    winslash = "/",
    mustWork = TRUE
  )
  
  # Temporary CSV where BOOTMODE saves its results
  output_file <- tempfile(fileext = ".csv")
  
  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  # Answers to BOOTMODE's prompts
  response_file <- tempfile(fileext = ".txt")
  
  writeLines(
    c(
      paste0('"', data_file, '"'),
      paste0('"', output_file, '"'),
      "",
      "y",
      as.character(outer_iterations),
      as.character(inner_iterations)
    ),
    response_file
  )
  
  response_file <- normalizePath(
    response_file,
    winslash = "/",
    mustWork = TRUE
  )
  
  python_script <- normalizePath(
    "bootmode_v3.py",
    winslash = "/",
    mustWork = TRUE
  )
  
  python_commands <- if (.Platform$OS.type == "windows") {
    c("python", "python3")
  } else {
    c("python3", "python")
  }

  python_locations <- unname(Sys.which(python_commands))
  python_exe <- python_locations[python_locations != ""][1]
  
  if (length(python_exe) == 0 || is.na(python_exe)) {
    stop("Python 3 could not be found on the server.")
  }
  
  python_exe <- normalizePath(
    python_exe,
    winslash = "/",
    mustWork = TRUE
  )
  
  on.exit(
    unlink(c(
      data_file,
      response_file,
      output_file
    )),
    add = TRUE
  )
  
  result <- system2(
    command = python_exe,
    args = shQuote(python_script),
    stdin = response_file,
    stdout = TRUE,
    stderr = TRUE
  )
  
  status <- attr(result, "status")
  
  if (!is.null(status) && status != 0) {
    stop(paste(result, collapse = "\n"))
  }
  
  if (!file.exists(output_file)) {
    stop(
      paste(
        "BOOTMODE did not create its results CSV.",
        paste(result, collapse = "\n")
      )
    )
  }
  
  readr::read_csv(
    output_file,
    show_col_types = FALSE
  )
}
#---UI--------------------------------------------------------------------------
tts_analysis_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h2("Time to Speciation Interval Analysis"),
    
    sidebarLayout(
      sidebarPanel(
        fileInput(
          ns("file"),
          label = span(
            "Upload CSV",
            tooltip(
              bsicons::bs_icon(
                "info-circle",
                title = "About file format"
              ),
              "Requires a .csv file with headers: species, crown_age, and stem_age.",
              placement = "right"
            )
          ),
          accept = ".csv"
        ),
        
        numericInput(
          ns("binwidth"),
          label = span(
            "Time step (Myr)",
            tooltip(
              bsicons::bs_icon(
                "info-circle",
                title = "About binwidth"
              ),
              "Determines the width of the time bins.",
              "Each species is counted in every bin overlapped by its crown-to-stem interval.",
              "For example, 5–8 Myr with 1 Myr bins contributes to 5–6, 6–7, and 7–8.",
              "An interval ending after 8 Myr, such as 5.1–8.1, also contributes to 8–9.",
              placement = "right"
            )
          ),
          value = 1,
          min = 0.01,
          step = 0.05
        ),
        
        numericInput(
          ns("bootstrap"),
          label = span(
            "Bootstrap replicates",
            tooltip(
              bsicons::bs_icon(
                "info-circle",
                title = "About bootstrap"
              ),
              "Number of bootstrap replicates used for the half-range mode (HRM).",
              "HRM resamples species intervals and provides a 95% confidence interval.",
              placement = "right"
            )
          ),
          value = 500,
          min = 10,
          step = 1
        ),
        
        numericInput(
          ns("bootmode_outer"),
          "BOOTMODE outer bootstrap iterations",
          value = 100,
          min = 100,
          step = 100
        ),
        
        numericInput(
          ns("bootmode_inner"),
          "BOOTMODE inner bootstrap iterations",
          value = 500,
          min = 10,
          step = 10
        ),
        
        radioButtons(
          ns("y_type"),
          "Y-axis",
          choices = c("Density", "Histogram"),
          selected = "Density"
        ),
        
        sliderInput(
          ns("x_limits"),
          label = span(
            "X-axis limits",
            tooltip(
              bsicons::bs_icon(
                "info-circle",
                title = "About slider bar"
              ),
              "Changes the visible range of the plots.",
              "Does not change the mode calculations.",
              placement = "right"
            )
          ),
          min = 0,
          max = 500,
          value = c(0, 1),
          step = 1
        ),
        
        div(
          style = "display: flex; flex-direction: column; gap: 10px;",
          
          downloadButton(
            ns("download_data"),
            "Download interval CSV"
          ),
          
          downloadButton(
            ns("download_plot"),
            "Download plot"
          ),
          
          downloadButton(
            ns("download_bootmode"),
            "Download BOOTMODE stats"
          )
        )
      ),
      
      mainPanel(
        accordion(
          accordion_panel(
            "TTS interval plot",
            plotOutput(
              ns("tts_plot"),
              height = "500px"
            )
          ),
          
          accordion_panel(
            "Data preview",
            div(
              style = "overflow-x: auto; width: 100%;",
              tableOutput(ns("preview"))
            )
          ),
          
          accordion_panel(
            "TTS estimate",
            div(
              style = "overflow-x: auto; width: 100%;",
              
              p(
                "Calculating bootstrap estimates may take a few moments."
              ),
              
              shinycssloaders::withSpinner(
                tableOutput(ns("tts_mode_table")),
                type = 4
              )
            )
          ),
          
          accordion_panel(
            "Methods",
            tableOutput(ns("methods_table"))
          ),
          
          open = "TTS interval plot",
          multiple = TRUE
        )
      )
    )
  )
}


#-------------------------------------------------------------------------------


tts_analysis_server <- function(id, generated_intervals = NULL) {
  moduleServer(id, function(input, output, session) {
  
  tts_data <- reactive({
    #requires data input from user
    data <- if (!is.null(input$file)) {
      
      read_csv(
        input$file$datapath,
        show_col_types = FALSE
      )
      
    } else {
      
      req(!is.null(generated_intervals))
      
      generated_data <- generated_intervals()
      
      req(
        !is.null(generated_data),
        nrow(generated_data) > 0
      )
      
      generated_data
    }
    
    lower_names <- tolower(names(data))
    
    species_options <- c("species", "sp.", "sp")
    crown_options <- c(
      "crown_age",
      "crown_time",
      "crown age (mya)",
      "crown age",
      "crown_age_mya",
      "crown"
    )
    stem_options <- c(
      "stem_age",
      "stem_time",
      "stem age (mya)",
      "stem age",
      "stem_age_mya",
      "stem"
    )
    
    species_match <- match(species_options, lower_names, nomatch = 0)
    crown_match <- match(crown_options, lower_names, nomatch = 0)
    stem_match <- match(stem_options, lower_names, nomatch = 0)
    
    species_col <- names(data)[species_match[species_match > 0][1]]
    crown_col <- names(data)[crown_match[crown_match > 0][1]]
    stem_col <- names(data)[stem_match[stem_match > 0][1]]
    
    
    if (any(is.na(c(species_col, crown_col, stem_col)))) {
      stop(
        paste(
          "CSV must contain:",
          "species or Species;",
          "crown, crown_time, crown_age, or Crown Age (mya);",
          "stem, stem_time, stem_age, or Stem Age (mya)."
        )
      )
    }
    
    data <- data |>
      rename(
        species = all_of(species_col),
        crown_age = all_of(crown_col),
        stem_age = all_of(stem_col)
      ) |>
      mutate(
        crown_age = as.numeric(crown_age),
        stem_age = as.numeric(stem_age)
      )
    
    validate(
      need(
        all(is.finite(data$crown_age)) &&
          all(is.finite(data$stem_age)),
        "Crown and stem ages must be numeric."
      ),
      need(
        all(data$crown_age >= 0),
        "Crown ages cannot be negative."
      ),
      need(
        all(data$stem_age > data$crown_age),
        "Each stem age must be older than its crown age."
      )
    )
    
    data |>
      mutate(row_id = row_number())
  
  })
  
#---Objects that calculate species bins-----------------------------------------
  

  # Reproduce the Excel time-point frequency calculation
  tts_point_frequencies <- reactive({
    data <- tts_data()
    time_step <- input$binwidth
    
    req(
      nrow(data) > 0,
      is.finite(time_step),
      time_step > 0
    )
    
    max_time <- ceiling(
      max(data$stem_age, na.rm = TRUE) / time_step
    ) * time_step
    
    time_points <- seq(
      from = time_step,
      to = max_time,
      by = time_step
    ) |>
      round(digits = 12)
    
    frequencies <- vapply(
      time_points,
      function(time) {
        sum(
          time > data$crown_age &
            time < data$stem_age,
          na.rm = TRUE
        )
      },
      integer(1)
    )
    
    tibble(
      time_point = time_points,
      frequency = frequencies,
      cumulative_frequency = cumsum(frequencies)
    )
  })
  
  
  # Reproduce the Excel column named Final
  bootmode_final <- reactive({
    point_data <- tts_point_frequencies()
    
    final_values <- rep(
      point_data$time_point,
      times = point_data$frequency
    )
    
    validate(
      need(
        length(final_values) > 0,
        "No time points occur inside the supplied crown–stem intervals."
      )
    )
    
    tibble(Final = final_values)
  })
  
  observeEvent(
    list(tts_data(), input$binwidth),
    {
      data <- tts_data()
      
      x_max <- ceiling(
        max(data$stem_age, na.rm = TRUE)
      )
      
      updateSliderInput(
        session,
        "x_limits",
        min = 0,
        max = x_max,
        value = c(0, x_max),
        step = 1
      )
    }
  )

#-------------------------------------------------------------------------------
  
  density_peak <- reactive({
    density_values <- bootmode_final() |>
      pull(Final)
    
    validate(
      need(
        length(unique(density_values)) > 1,
        "At least two distinct occupied time points are required."
      )
    )
    
    d <- density(
      density_values,
      na.rm = TRUE
    )
    
    peak_time <- d$x[which.max(d$y)]
    
    message(
      sprintf(
        "Preliminary smooth-density TTS peak: %.2f Mya",
        peak_time
      )
    )
    
    tibble(
      peak_time = peak_time,
      max_density = max(d$y)
    )
  })
  
  tts_mode <- reactive({
    tts_point_frequencies() |>
      filter(
        frequency > 0,
        frequency == max(frequency, na.rm = TRUE)
      )
  })
  
  hrm_result <- reactive({
    bootmode_values <- bootmode_final() |>
      pull(Final)
    
    run_hrm_python(
      values = bootmode_values,
      outer_iterations = as.integer(input$bootmode_outer),
      inner_iterations = as.integer(input$bootmode_inner)
    )
  })
  
  bootstrap_hrm <- reactive({
    data <- tts_data()
    time_step <- input$binwidth
    n_boot <- as.integer(input$bootstrap)
    
    req(
      nrow(data) > 0,
      is.finite(time_step),
      time_step > 0,
      is.finite(n_boot),
      n_boot > 0
    )
    
    set.seed(7)
    
    bootstrap_modes <- replicate(n_boot, {
      
      # Resample whole species intervals
      boot_data <- data |>
        slice_sample(
          n = nrow(data),
          replace = TRUE
        )
      
      max_time <- ceiling(
        max(boot_data$stem_age, na.rm = TRUE) / time_step
      ) * time_step
      
      time_points <- seq(
        from = time_step,
        to = max_time,
        by = time_step
      ) |>
        round(digits = 12)
      
      frequencies <- vapply(
        time_points,
        function(time) {
          sum(
            time > boot_data$crown_age &
              time < boot_data$stem_age,
            na.rm = TRUE
          )
        },
        integer(1)
      )
      
      final_values <- rep(
        time_points,
        times = frequencies
      )
      
      if (length(final_values) == 0) {
        return(NA_real_)
      }
      
      genefilter::half.range.mode(final_values)
    })
    
    bootstrap_modes <- bootstrap_modes[
      is.finite(bootstrap_modes)
    ]
    
    validate(
      need(
        length(bootstrap_modes) > 0,
        "No valid R bootstrap estimates were produced."
      )
    )
    
    tibble(
      hrm_mode = mean(bootstrap_modes),
      ci_lower = unname(
        quantile(bootstrap_modes, 0.025)
      ),
      ci_upper = unname(
        quantile(bootstrap_modes, 0.975)
      ),
      standard_error = sd(bootstrap_modes)
    )
  })
  
  locmode_result <- reactive({
    x <- bootmode_final() |>
      pull(Final)
    
    fit <- multimode::locmodes(
      data = x,
      mod0 = 1,
      display = FALSE
    )
    
    tibble(
      peak_time = fit$locations[1],
      max_density = fit$fvalue[1]
    )
  })
  
  output$preview <- renderTable({
    head(tts_data(), 10)
  })
  
  tts_plot_object <- reactive({
    final_data <- bootmode_final()
    point_data <- tts_point_frequencies() |>
      filter(frequency > 0)
    
    mode_points <- tts_mode()
    peak <- density_peak()
    
    if (input$y_type == "Density") {
      
      ggplot(final_data, aes(x = Final)) +
        geom_density(linewidth = 1, fill = "turquoise") +
        geom_vline(
          xintercept = peak$peak_time,
          linetype = "dashed",
          linewidth = 1
        ) +
        coord_cartesian(xlim = input$x_limits) +
        theme_classic() +
        labs(
          x = "Interval times (my)",
          y = "Density",
          title = "TTS interval density through time"
          )+
        theme_classic(base_size = 12) +
        theme(
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(
            size = 12,
            hjust = 0.5
          )
        )
      
    } else {
      
      ggplot(
        point_data,
        aes(x = time_point, y = frequency)
      ) +
        geom_col(
          width = input$binwidth,
          fill = "turquoise",
          color = NA
        ) +
        geom_vline(
          data = mode_points,
          aes(xintercept = time_point),
          linetype = "dashed",
          linewidth = 1
        ) +
        coord_cartesian(xlim = input$x_limits) +
        theme_classic() +
        labs(
          x = "Time before present (mya)",
          y = "Number of species intervals",
          title = "TTS interval occupancy through time"
          )+
        theme_classic(base_size = 12) +
        theme(
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(
            size = 12,
            hjust = 0.5
          )
        )
    }
  })
  
  output$tts_plot <- renderPlot({
    tts_plot_object()
  })
  
  output$tts_mode_table <- renderTable({
    mode_bin <- tts_mode()
    peak <- density_peak()
    hrm <- bootstrap_hrm()
    
    bootmode <- hrm_result() |>
      slice(1)
    
    n_modes <- nrow(mode_bin)
    
    results <- mode_bin |>
      mutate(
        `Binned mode` = if (n_modes == 1) {
          "Unique mode"
        } else {
          paste0(
            "Tied mode ",
            seq_len(n_modes),
            " of ",
            n_modes
          )
        }
      ) |>
      transmute(
        `Binned mode` = `Binned mode`,
        `Modal time point (Mya)` = time_point,
        `Species interval frequency` = frequency,
        
        `Smooth density peak (Mya)` =
          round(peak$peak_time, 2),
        
        `R bootstrap HRM estimate (Mya)` =
          round(hrm$hrm_mode, 2),
        
        `R bootstrap HRM 95% CI lower` =
          round(hrm$ci_lower, 2),
        
        `R bootstrap HRM 95% CI upper` =
          round(hrm$ci_upper, 2),
        
        `R bootstrap HRM SE` =
          round(hrm$standard_error, 2),
        
        `BOOTMODE occupancy HRM estimate (Mya)` =
          round(bootmode$mode, 2),
        
        `BOOTMODE 95% CI lower` =
          round(bootmode$ci_low_2_5_percent, 2),
        
        `BOOTMODE 95% CI upper` =
          round(bootmode$ci_high_97_5_percent, 2),
        
        `BOOTMODE bootstrap SE` =
          round(bootmode$bootstrap_std_error, 2)
      )
    
    if (n_modes > 1) {
      results[2:n_modes, 5:ncol(results)] <- NA
    }
    
    results
  },
  rownames = FALSE,
  na = "")
  
  
  output$methods_table <- renderTable({
    tibble(
      `Method component` = c(
        "Input species intervals",
        "Time step",
        "Occupancy mode",
        "Smooth density",
        "BOOTMODE input",
        "BOOTMODE outer iterations",
        "BOOTMODE inner iterations",
        "R bootstrap unit",
        "R bootstrap replicates",
        "Confidence intervals"
      ),
      `Summary` = c(
        nrow(tts_data()),
        paste0(input$binwidth, " Myr"),
        paste(
          "Time point or tied time points with the greatest",
          "crown–stem interval occupancy"
        ),
        paste(
          "Kernel density calculated from time points repeated",
          "according to crown–stem interval occupancy"
        ),
        paste(
          "Time points repeated according to",
          "crown–stem interval occupancy"
        ),
        input$bootmode_outer,
        input$bootmode_inner,
        paste(
          "Species intervals resampled, followed by reconstruction",
          "of the time-point occupancy distribution"
        ),
        input$bootstrap,
        "95% percentile intervals"
      )
    )
  },
  rownames = FALSE)
  
  
  output$download_data <- downloadHandler(
    filename = function() {
      "tts_intervals.csv"
    },
    content = function(file) {
      write_csv(tts_data(), file)
    }
  )
  
  output$download_bootmode <- downloadHandler(
    filename = function() {
      paste0(
        "bootmode_statistics_",
        Sys.Date(),
        ".csv"
      )
    },
    
    content = function(file) {
      bootmode_stats <- hrm_result() |>
        select(-source_file)
      
      write_csv(
        bootmode_stats,
        file
      )
    }
  )
  
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(
        "tts_",
        tolower(input$y_type),
        "_plot.png"
      )
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = tts_plot_object(),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  )
  
  })
}
