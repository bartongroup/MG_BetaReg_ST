# GENE EXPLORER
#
# Global settings and variables
#

require(dplyr)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

options(dplyr.summarise.inform = FALSE)
select <- dplyr::select
rename <- dplyr::rename

# Source all R files

files_R <- list.files("R", pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

