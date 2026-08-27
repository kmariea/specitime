library(shiny)
library(tidyverse)
library(bslib)
library(ape)
library(multimode)
library(genefilter)
library(shinycssloaders)

source("Shiny_Times_analysisTTS.R")
source("generate_intervals.R")
source("home_page.R")

ui <- page_navbar(
  id = "main_nav",
  title = "Time to Speciation App",
  
  theme = bs_theme(
    version = 5,
    preset = "minty"
  ),
  
  nav_panel(
    "Home",
    value = "home",
    home_ui("home")
  ),
  
  nav_panel(
    "Generate intervals",
    value = "generate",
    generate_intervals_ui("trees")
  ),
  
  nav_panel(
    "Analyze intervals",
    value = "analysis",
    tts_analysis_ui("analysis")
  )
)

server <- function(input, output, session) {
  
  home_server("home")
  
  generated_intervals <- generate_intervals_server("trees")
  
  tts_analysis_server(
    "analysis",
    generated_intervals = generated_intervals
  )
}

shinyApp(ui, server)

