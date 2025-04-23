# GENE SET ENRICHMENT EXPLORER
#
# Global settings and variables
#


library(shiny)
#library(shinydashboard)
library(jsonlite)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

options(dplyr.summarise.inform = FALSE)
select <- dplyr::select
rename <- dplyr::rename

# Source all R files

files_R <- list.files(c("../R", "modules"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

# Read config json file

config <- jsonlite::read_json("../shiny_config.json")

# Read all data, make it available to all modules

shd <- sh_read_data(c("dexset", "features", "gse"), config, with_progress = FALSE)
 
DATA <- list(
  dex = shd$dexset,
  names = names(shd$dexset),
  features = shd$features,
  gse = shd$gse
)

DATA$id2name <- set_names(DATA$features$name, DATA$features$id)
DATA$name2id <- set_names(DATA$features$id, DATA$features$name)

