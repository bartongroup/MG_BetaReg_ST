# GENE EXPLORER
#
# Main server
#

server <- function(input, output, session) {
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # Initialise app state, reactive object for communication between modules
  # Note: an undefined reactive value defaults to NULL, so we can start with empty list
  app_state <- new_app_state()
  
  # server logic: modules
  mod_global_input_server("global_input", DATA, app_state)
  mod_dot_plot_server("dot_plot", DATA, app_state)
  mod_spatial_plot_server("spatial_plot", DATA, app_state)
}
