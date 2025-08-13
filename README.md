# BetaReg_ST

This repository contains software for processing and analysing Visium HD data from 8 wpc pancreatic samples.  
It is written in R and uses the [`targets`](https://docs.ropensci.org/targets/) package to build a fully reproducible analysis pipeline.  

The pipeline definitions can be found in the following files:

- `R/tar_main.R`  
- `R/tar_8um.R`  
- `R/tar_2um.R`  

## Running the `targets` Pipeline

Some stages of the pipeline are memory-intensive and require long runtimes. These steps should be executed on a high-memory cluster node with minor adjustments to the workflow.

### Set up the R environment

Create and restore the analysis environment using [`renv`](https://rstudio.github.io/renv/):

```
install.packages("renv")
renv::restore()
```

### Prepare input data for the external script

Run the following command to generate the input files:

```
targets::tar_make(tidyselect::starts_with("sav_external_cor"))
```

This will produce files named `ext/seu_?.rds`, where `?` is the sample number (e.g., `1` and `2` in our case).

### Run the correlation script on the cluster

Copy the `.rds` files and the script `scripts/cscore_correlation_all.R` to the cluster.  
Run the script for each file:

```
cscore_correlation_all.R seu_1.rds cor_1.parquet
cscore_correlation_all.R seu_2.rds cor_2.parquet
```

**Notes:**

- It is recommended to submit these jobs via `qsub`.
- Request **at least 100 GB of RAM**.
- On our cluster, these jobs took **12–16 hours** to complete.

### Retrieve and link the results

When the jobs finish, copy the generated `.parquet` files back to your laptop into the `shiny/data` directory (for use in the Shiny apps). Also, create symbolic links (rather than duplicating these large files) in the `ext` directory.

### Complete the pipeline

Resume the pipeline:

```
targets::tar_make()
```


### Generate the report and run Shiny apps

Render the Quarto report:

```
quarto::quarto_render("doc/report.qmd")
```

Launch the Shiny applications:


```
shiny::runApp("shiny/de")
shiny::runApp("shiny/gsea")
shiny::runApp("shiny/gene_explorer")
```

