# Module COMMUNICATION
#
# Facilitates communication between modules.


mod_communication_server <- function(id, data_set, state) {
  
  server <- function(input, output, session) {
    
  
    moduleServer(id, server)
  }
}