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

# Read Seurat data

DATA <- list()
DATA$`1` <- qs2::qs_read("../data/seu_1.qs")
DATA$`2` <- qs2::qs_read("../data/seu_2.qs")

# Open correlation files

COR <- list()
COR$`1` <- arrow::open_dataset("../data/cor_1.parquet")
COR$`2` <- arrow::open_dataset("../data/cor_2.parquet")

# Read config json file

all_genes <- SeuratObject::Features(DATA$`1`)
GENES <- read_genes_for_pathways("genes_of_interest.xlsx", all_genes)

# Read gene descriptions

GENE_INFO <- qs2::qs_read("../data/gene_info.qs")
