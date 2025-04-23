library(dexdash)
library(readr)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
options(dplyr.summarise.inform = FALSE)

# Source all R files

files_R <- list.files(c("../R"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)


# Read config json file

config <- jsonlite::read_json("../shiny_config.json")

# Read all data, make it available to all modules

shd <- sh_read_data(c("dexset", "features", "fterms"), config, with_progress = FALSE)
dexdash::run_app(shd$dexset, shd$features, shd$fterms, title = "Differential expression vs MPC",
                 x_variable = "cell_type", colour_variable = "group")
