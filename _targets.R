suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(qs2)
})

options(tidyverse.quiet = TRUE, dplyr.summarise.inform = FALSE)

# attach R packages
required_packages <- read.delim("R/packages.txt", header = FALSE, col.names = "name")$name
tar_option_set(packages = required_packages, format = "qs")

# Create dirs if necessary
for (d in c("data", "tab", "fig", "doc/images")) if (!dir.exists(d)) dir.create(d)

# for interactive session only
if (interactive()) sapply(required_packages, library, character.only = TRUE)

# load all functions from .R files
files_R <- list.files(c("R", "targets"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

# prevent other packages from stealing from tidyverse
filter <- dplyr::filter
rename <- dplyr::rename
select <- dplyr::select
slice <- dplyr::slice
count <- dplyr::count

# Add session info
sessinfo <- list(
  tar_force(session_info, sessioninfo::session_info(), force = TRUE)
)

# Targets
c(
  sessinfo,
  targets_main()
)

