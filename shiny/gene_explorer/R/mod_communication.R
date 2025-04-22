# Module COMMUNICATION
#
# Facilitates communication between modules.

require(shiny)

mod_communication_server <- function(id, state, data) {
  
  server <- function(input, output, session) {
    
    shiny::observe({
      gid <- NULL
      if(!is.null(state$gene_input))
        gid <- state$gene_input
      if(!is.null(state$gene_cor))
        gid <- state$gene_cor
      state$gene_plot <- gid
    })

  }
  
  shiny::moduleServer(id, server)
}