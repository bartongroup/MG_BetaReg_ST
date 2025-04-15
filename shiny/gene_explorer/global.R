# GENE EXPLORER
#
# Global settings and variables
#


library(shiny)
library(shinydashboard)
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

GENES <- jsonlite::read_json("ge_genes.json", simplifyVector = TRUE)

# Read all data, make it available to all modules

DATA <- list()
DATA$`1` <- qs2::qs_read("../data/seu_1.qs")
DATA$`2` <- qs2::qs_read("../data/seu_2.qs")
