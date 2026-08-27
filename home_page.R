home_ui <- function(id) {
  ns <- NS(id)
  
  page_fillable(
    div(
      class = "text-center py-5",
      
      h1("Time to Speciation"),
      
      p(
        class = "lead",
        "Generate crown–stem intervals from a time-calibrated phylogenetic tree
        and estimate the time to speciation.", br(),
        "This work is based on Tree of Life Reveals Clock-Like Speciation and 
        Diversification (Hedges et al., 2015).")
    ),
    
    layout_columns(
      card(
        card_header("1. Generate intervals"),
        
        p(
          "Upload a time-calibrated phylogenetic tree and automatically generate
          the crown–stem intervals needed for the analysis.",
          br(),
          "A crown node contains all of the population-level data for a species
          while a stem node is the most recent common ancestor (MRCA) to a species",
          br(),
          "Automatically generating the intervals requires specific formating of
          the tip labels (genus_species_uniqueID), which should be contained in
          a time-calibrated phylogentic tree of newick or .tre format.",
          br(),
          "The output is a .csv file with crown and stem intervals of each species
          with a crown node. Multiple trees can be analysed and added to a single
          .csv file."
        ),
        
        card_footer(
          actionLink(
            ns("go_to_generate"),
            "Generate intervals →"
          )
        )
      ),
      
      card(
        card_header("2. Analyze intervals"),
        
        p(
          "Estimate the time to speciation from generated
          or previously saved crown–stem intervals."
        ),
        
        card_footer(
          actionLink(
            ns("go_to_analysis"),
            "Analyze intervals →"
          )
        )
      )
    )
  )
}

home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$go_to_generate, {
      nav_select(
        "main_nav",
        selected = "generate",
        session = session$rootScope()
      )
    })
    
    observeEvent(input$go_to_analysis, {
      nav_select(
        "main_nav",
        selected = "analysis",
        session = session$rootScope()
      )
    })
  })
}