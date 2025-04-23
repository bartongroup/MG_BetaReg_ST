# BetaReg_ST



## Running targets pipeline

Part of the pipeline requires large memory and a lot of time, so it needs to be executed on the cluster node. This needs some tweaks. First, create an R environment using *renv*:

```
install.packages("renv")
renv::restore()
```

Next, run a part of the pipeline to create input data for the external script:

```
targets::tar_make(tidyselect::starts_with("sav_external_cor"))
```

This will create files `ext/seu_?.rds`, which `?` corresponds to the sample number (`1` and `2` in our case). Copy them to the cluster together with the script `scripts/cscore_correlation_all.R`. Then run the script for each of the files:

```
cscore_correlation_all.R seu_1.rds cor_1.parquet
cscore_correlation_all.R seu_2.rds cor_2.parquet
```

It's best to use `qsub` to submit these commands to the cluster. Make sure you request at least 100 GB of RAM. When this is ready (took 12-16 hours on our cluster), copy the parquet files back to the laptop, to the `shiny/data` directory. They will be used by the Shiny app. Also, create a link (do not make a second copy of these massive files) of each parquet file in the `ext` directory. Now you can continue running the pipeline:

```
targets::tar_make()
```

Once this is done, you can render the Quarto report with

```
quarto::quarto_render("doc/report.qmd")
```

and run the Shiny apps

```
shiny::runApp("shiny/de")
shiny::runApp("shiny/gsea")
shiny::runApp("shiny/gene_explorer")
```