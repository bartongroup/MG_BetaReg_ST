okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#666666")

th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )


thm <- th +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

# Minimal dependencies: rlang (installed with ggplot2)
# Written by ChatGPT 5
# Warning: this does not work
gg_isolate_env <- function(p) {
  stopifnot(inherits(p, "ggplot"))

  # Recursively set environments of any quosures to emptyenv()
  strip_quosure_envs <- function(x) {
    if (rlang::is_quosure(x)) {
      return(rlang::quo_set_env(x, rlang::empty_env()))
    }
    if (is.list(x)) {
      return(lapply(x, strip_quosure_envs))
    }
    x
  }

  # 1) Kill the plot-level environment
  p$plot_env <- rlang::empty_env()

  # 2) Strip envs from top-level aes()
  if (length(p$mapping)) {
    p$mapping <- strip_quosure_envs(p$mapping)
  }

  # 3) Strip envs from each layer's aes()
  if (length(p$layers)) {
    p$layers <- lapply(p$layers, function(layer) {
      if (length(layer$mapping)) {
        layer$mapping <- strip_quosure_envs(layer$mapping)
      }
      layer
    })
  }

  # 4) Strip envs from facet quosures (wrap, grid, etc.)
  if (!is.null(p$facet) && !is.null(p$facet$params)) {
    params <- p$facet$params
    for (nm in c("facets", "rows", "cols")) {
      if (!is.null(params[[nm]])) {
        params[[nm]] <- strip_quosure_envs(params[[nm]])
      }
    }
    p$facet$params <- params
  }

  # 5) (Optional) If you store custom data as quosures in labels/scales, clean them too
  if (length(p$labels)) p$labels <- strip_quosure_envs(p$labels)
  if (!is.null(p$scales) && length(p$scales$scales)) {
    p$scales$scales <- lapply(p$scales$scales, function(s) {
      # some scales may carry expressions/quosures in labels/trans breaks
      s[] <- strip_quosure_envs(unclass(s))
      class(s) <- attr(s, "class")
      s
    })
  }

  p
}

#' Save plot to PNG in doc/images
#'
#' Wrapper around ggsave() that writes a ggplot to doc/images/<name>[_<sample>].png.
#'
#' @param gg A ggplot object.
#' @param name Base file name (without extension).
#' @param sample Optional sample name to append to the file.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution in dots per inch.
#'
#' @return Invisibly, the file path written.
gs <- function(gg, name, sample = NULL, width, height, dpi = 240, path = "doc/images") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if(is.null(sample)) {
    fname <- file.path(path, paste0(name, ".png"))
  } else {
    fname <- file.path(path, paste0(name, "_", sample, ".png"))
  }

  ggsave(filename = fname, plot = gg, device = "png", limitsize = FALSE,
         width = width, height = height, dpi = dpi)
}

gp <- function(gg, name, sample = NULL, width, height, path = "doc/images") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if(is.null(sample)) {
    fname <- file.path(path, paste0(name, ".pdf"))
  } else {
    fname <- file.path(path, paste0(name, "_", sample, ".pdf"))
  }

  ggsave(filename = fname, plot = gg, device = "pdf", limitsize = FALSE,
         width = width, height = height)
}


#' Save R object to RDS in doc/images
#'
#' Writes a serialised object to doc/images/<name>[_<sample>].rds.
#'
#' @param gg An R object to serialise (often a ggplot).
#' @param name Base file name (without extension).
#' @param sample Optional sample name to append to the file.
#'
#' @return Invisibly, the file path written.
rs <- function(gg, name, sample = NULL) {
  if(is.null(sample)) {
    fname <- file.path("doc", "images", paste0(name, "rds"))
  } else {
    fname <- file.path("doc", "images", paste0(name, "_", sample, ".rds"))
  }

  readr::write_rds(gg, fname)
}

#' Plotmath labels in scientific notation (×10^k)
#'
#' Formats numeric values as expressions like 10^k suitable for axis labels.
#'
#' @param x Numeric vector.
#'
#' @return An expression vector for use in plot labels.
scientific_10 <- function(x) {
  parse(
    text = scales::scientific_format()(x) |>
      str_replace("e\\+*", " %*% 10^") |> 
      str_replace("^1 \\%\\*\\% 10\\^", "10^")
  )
}

#' Table of QC thresholds for reporting
#'
#' Creates a two-column human-readable table of per-cell and per-gene cutoffs.
#'
#' @param params A list/data frame with fields min_detected_genes_per_cell, min_counts_per_cell,
#'   max_mt_fraction, and min_detections_per_gene.
#'
#' @return A tibble with columns Group, Quantity, and Limit.
make_table_cuts <- function(params) {
  tribble(
    ~Group, ~Quantity, ~Limit,
    "Per cell", "Number of detected genes per cell", str_c(">", params$min_detected_genes_per_cell),
    "Per cell", "Number of counts per cell", str_c(">", params$min_counts_per_cell),
    "Per cell", "Mitochondrial count fraction", str_c("<", params$max_mt_fraction),
    "Per gene", "Number of cells a gene is detected in", str_c(">", params$min_detections_per_gene)
  )
}

#' Distributions of per-cell quantities
#'
#' Density plots for total counts, detected genes, and mitochondrial fraction.
#'
#' @param cell_data Tibble of cell-level metrics with columns total_count, n_genes, mito_fraction.
#' @param params Optional list with QC thresholds to draw reference lines.
#' @param with_vlines Logical; add vertical lines at thresholds.
#'
#' @return A patchwork of three ggplots.
plot_cell_marginals <- function(cell_data, params = NULL, with_vlines = TRUE) {
  pl_ <- function(what, xlab, log_scale = TRUE) {
    pl <- cell_data |> 
      mutate(val = get(what)) |> 
    ggplot(aes(x = val)) +
      th +
      theme(legend.position = "none") +
      #geom_density_ridges() +
      geom_density(fill = "grey70") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
      #scale_fill_manual(values = okabe_ito_palette) +
      labs(x = xlab, y = NULL)
    if(log_scale)
      pl <- pl + scale_x_log10(labels = scales::comma)
    pl
  }
  
  pl1 <- pl_("total_count", "Total count per cell")
  if(with_vlines) pl1 <- pl1 + geom_vline(xintercept = params$min_counts_per_cell, linetype = "dotted")
  pl2 <- pl_("n_genes", "Genes detected per cell")
  if(with_vlines) pl2 <-pl2 + geom_vline(xintercept = params$min_detected_genes_per_cell, linetype = "dotted")
  pl3 <- pl_("mito_fraction", "Mitochondrial fraction", log_scale = FALSE) +
  if(with_vlines) pl3 <- pl3 + geom_vline(xintercept = params$max_mt_fraction, linetype = "dotted")
  
  #cowplot::plot_grid(pl1, pl2, pl3, nrow = 1, align = "h", labels = c("A", "B", "C"))
  
  # cowplot throws cryptic errors while in targets; using patchwork instead
  pl1 + pl2 + pl3
  
}

#' Histogram of gene detection counts
#'
#' For each group, counts how many genes are detected in 0..k cells and plots per group.
#'
#' @param gene_detection Tibble with one row per gene and columns per group containing detection counts.
#' @param max_det Integer; values above are collapsed into >max_det.
#'
#' @return A faceted ggplot.
plot_gene_detections <- function(gene_detection, max_det = 8) {
  gene_detection |> 
    pivot_longer(-gene_id, names_to = "group", values_to = "detections") |> 
    group_by(group, detections) |> 
    tally() |> 
    mutate(
      det = as.character(detections),
      det = if_else(detections > max_det, sprintf(">%d", max_det), det) |> as_factor()
    ) |> 
    group_by(group, det) |>
    summarise(tot = sum(n)) |> 
  ggplot(aes(x = det, y = tot, fill = group)) +
    th +
    geom_col() +
    facet_wrap(~group) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of detections", y = "Number of genes") 
}

#' Genes detected vs total counts per cell
#'
#' Scatter of n_genes vs total_count coloured by mito_fraction with optional cut lines.
#'
#' @param cell_data Tibble of cell-level metrics.
#' @param params Optional list with min_detected_genes_per_cell and min_counts_per_cell.
#' @param with_lines Logical; add threshold lines.
#'
#' @return A ggplot object.
plot_count_det <- function(cell_data, params = NULL, with_lines = TRUE) {
  pl <- cell_data |> 
    ggplot(aes(x = total_count, y = n_genes, colour = mito_fraction)) +
    th +
    geom_point(size = 0.5) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_viridis_c(option = "cividis") +
    labs(x = "Total count per cell", y = "Genes detected per cell", colour = "MT fraction")
  if(with_lines) {
    pl <- pl +
      geom_hline(yintercept = params$min_detected_genes_per_cell, linetype = "dotted") +
      geom_vline(xintercept = params$min_counts_per_cell, linetype = "dotted")
  }
  return(pl)
}


#' Counts and sizes across samples
#'
#' Three side-by-side bars: number of cells, total reads (millions), and median genes per cell.
#'
#' @param cell_data_all Tibble of cell-level metrics with sample, total_count.
#' @param gene_detection_all Tibble with gene-by-sample detections.
#'
#' @return A cowplot grid of three ggplots.
plot_all_counts <- function(cell_data_all, gene_detection_all) {
  env <- new.env(parent = globalenv())

  env$d1 <- cell_data_all |> 
    count(sample)
  env$d3 <- gene_detection_all |> 
    pivot_longer(-c(gene_symbol, sample)) |> 
    summarise(tot = sum(value), .by = c("gene_symbol", "sample")) |> 
    summarise(n = median(tot), .by = "sample")
  env$d2 <- cell_data_all |> 
    summarise(n = sum(total_count) / 1e6, .by = "sample")

  with(env, {
    pl_ <- function(d, ylab) {
      ggplot(d, aes(x = sample, y = n)) +
      th +
      geom_col() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Sample", y = ylab)
    }

    pl1 <- pl_(d1, "Cell count")
    pl3 <- pl_(d3, "Median gene count per cell")
    pl2 <- pl_(d2, "Total read count (millions)")
  
    plot_grid(pl1, pl2, pl3, nrow = 1)
  })
}

#' Spatial plot of a metadata field
#'
#' Renders a Visium tile map filled by a metadata column, with optional divergent palette.
#'
#' @param seu Seurat object.
#' @param meta_column Metadata column name to show.
#' @param palette Either "viridis" (continuous) or "divergent" (RdBu centred).
#'
#' @return A ggplot of polygons per spot.
plot_spatial_meta <- function(seu, meta_column = "pseudotime", palette = "viridis") {
  xy <- get_spatial_coordinates(seu, swap_xy = FALSE) |> 
    mutate(val = seu@meta.data[[meta_column]])
  corners <- get_resolution_grid(xy)
  ex <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    mutate(val = get(meta_column)) |> 
    select(cell_id, val)
  d <- inner_join(ex, corners, by = "cell_id")
  mx <- max(abs(d$val))
  pl <-  d |> 
  ggplot(aes(x = x_corner, y = y_corner, fill = val)) +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = val), colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = meta_column)
  if(palette == "viridis") {
    pl <- pl + scale_fill_viridis_c(option = "cividis")
  } else if(palette == "Purples") {
    pl <- pl + scale_fill_distiller(palette = "Purples", direction = 1)
  } else if(palette == "divergent") {
    pl <- pl + scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx))
  }
  pl
}

#' Percentile bands ignoring zero expression
#'
#' Replaces zero values with 0 and assigns non-zeros to percentile bands [0, percentile)
#' and >= percentile.
#'
#' @param df Tibble with column expression.
#' @param percentile Numeric in (0,1); band split threshold.
#'
#' @return Tibble with added columns prc and factor val in {0,1,2}.
percentiles_no_zeroes <- function(df, percentile) {
  df |> 
    mutate(
      prc = if_else(expression == 0, NA, expression),
      prc = rank(prc, na.last = "keep") / length(na.omit(prc)),
      prc = replace_na(prc, 0), 
      val = case_when(
        prc == 0 ~ 0,
        prc < percentile ~ 1,
        prc >= percentile ~ 2
      ) |> as.factor()
    )
}

#' Spatial plots for selected genes
#'
#' Spatial plot for each gene's expression or discrete percentile bands; optional highlighted spots.
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols.
#' @param box Named list with x1, x2, y1, y2 to crop (defaults to full extent).
#' @param percentile Optional numeric in (0,1) for banded fill using \code{percentiles_no_zeroes()}.
#' @param spots Optional tibble with columns x,y to outline circular regions.
#' @param ncol Number of columns in facet layout.
#'
#' @return A faceted ggplot.
plot_spatial_selection_genes <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf),
                                         percentile = NULL, spots = NULL, ncol = 2) {
  ex <- get_gene_expression(seu, genes) |> 
    mutate(val = expression)
  if(!is.null(percentile)) {
    # percentiles, but ignore zeroes
    ex <- ex |> 
      group_by(gene_symbol) |>
      percentiles_no_zeroes(percentile) |> 
      ungroup()
  }

  xy <- get_spatial_coordinates(seu, swap_xy = FALSE) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, ex, by = "cell_id", relationship = "many-to-many")

  pl <- d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = val), colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    facet_wrap(~gene_symbol, ncol = ncol)

  if(is.null(percentile)) {
    pl <- pl + scale_fill_viridis_c(option = "cividis")
  } else {
    pl <- pl + scale_fill_manual(values = c("beige", "darkseagreen", "black")) +
      theme(legend.position = "none")
  }

  r <- 20
  a <- 2 * pi * seq(1, 30) / 29
  if(!is.null(spots)) {
    sp <- spots |> 
      mutate(idx = row_number()) |> 
      group_split(idx) |> 
      map(function(w) {
        tibble(
          idx = w$idx,
          fx = w$x + r * cos(a),
          fy = w$y + r * sin(a)
        )
      }) |> 
      list_rbind()
    pl <- pl +
      geom_path(data = sp, aes(x = fx, y = fy, group = idx), colour = "red", linewidth = 0.5, alpha = 0.7)
  }
  pl
}

#' Animated spatial gene map
#'
#' Creates a GIF cycling through genes, showing expression or percentile bands per frame.
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols (order defines animation).
#' @param box Named list with x1, x2, y1, y2 to crop.
#' @param percentile Optional numeric in (0,1) for banded fill.
#' @param transition_length Transition length passed to gganimate.
#'
#' @return Saves doc/images/spatial_<genes>.gif; invisibly returns the gganimate object.
animate_spatial_genes <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf),
                                  percentile = NULL, transition_length = 1) {
  ex <- get_gene_expression(seu, genes) |> 
    mutate(val = expression)
  if(!is.null(percentile)) {
    # percentiles, but ignore zeroes
    ex <- ex |> 
      group_by(gene_symbol) |>
      percentiles_no_zeroes(percentile) |> 
      ungroup()
  }

  xy <- get_spatial_coordinates(seu) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, ex, by = "cell_id", relationship = "many-to-many") |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes))

  pl <- d |> 
    ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = val)) +
    thm +
    geom_polygon(colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    theme(legend.position = "none") +
    labs(title = "{closest_state}") +
    gganimate::transition_states(gene_symbol, transition_length = transition_length)

    if(is.null(percentile)) {
    pl <- pl + scale_fill_viridis_c(option = "cividis")
  } else {
    pl <- pl + scale_fill_manual(values = c("beige", "darkseagreen", "black")) +
      theme(legend.position = "none")
  }

  sgenes <- str_c(genes, collapse = "_")
  fname <- str_glue("spatial_{sgenes}.gif")
  gganimate::anim_save(fname, pl, path = file.path("doc", "images"), renderer = gganimate::gifski_renderer())
}

#' Spatial plot of categorical metadata
#'
#' Fills Visium tiles by a discrete metadata column (e.g., cell_type) using a fixed palette. Designed
#' for cell type (that is, don't try it with a different meta_column)
#'
#' @param seu Seurat object.
#' @param meta_column Metadata column name to show (default "cell_type").
#' @param box Named list with x1, x2, y1, y2 to crop.
#'
#' @return A ggplot of polygons per spot.
plot_spatial_selection_meta <- function(seu, meta_column = "cell_type", box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf)) {
  mt <- get_cell_data(seu) |> 
    mutate(val = seu@meta.data[[meta_column]]) |> 
    select(cell_id, val)
  xy <- get_spatial_coordinates(seu) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, mt, by = "cell_id", relationship = "many-to-many") 

  pal <- c(okabe_ito_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))
  d |> 
    ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = val)) +
    thm +
    geom_polygon(colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    scale_fill_manual(values = pal)
}

#' Reduction plot coloured by metadata
#'
#' Scatter on reduction coordinates coloured by a continuous metadata variable.
#'
#' @param seu Seurat object.
#' @param reduction Reduction name ("pca", "umap", etc.).
#' @param meta_column Metadata column to colour by.
#' @param point_size Point size.
#'
#' @return A ggplot object.
plot_reduction_meta <- function(seu, reduction = "pca", meta_column = "pseudotime", point_size = 0.6) {
  Embeddings(seu, reduction = reduction) |>
    as_tibble(rownames = "cell_id") |>
    select(1:3) |> 
    set_names(c("cell_id", "x", "y")) |> 
    mutate(val = seu@meta.data[[meta_column]]) |> 
  ggplot(aes(x = x, y = y, colour = val)) +
    thm +
    geom_point(size = point_size, shape = 15, stroke = 0.1) +
    scale_colour_viridis_c(option = "cividis") +
    labs(x = NULL, y = NULL, colour = meta_column)
}



#' UMAP or spatial plot for one quantity
#'
#' Plots either a discrete metadata column or gene expression (one or multiple genes),
#' on UMAP or spatial tiles.
#'
#' @param seu Seurat object.
#' @param resolution Optional clustering resolution used to set Idents (e.g. "0.2").
#' @param point_size Point size for scatter/polygons.
#' @param what One of "umap" or "spatial".
#' @param meta_column Name of metadata column when plotting categories.
#' @param genes Optional gene symbol(s); if given, expression is plotted.
#' @param with_legend Logical; show legend.
#' @param label_size Size of cluster labels when plotting discrete groups.
#' @param combine_genes Logical; if multiple genes, average per cell before plotting.
#' @param discrete_palette Optional vector of colours for discrete groups.
#'
#' @return A ggplot object (optionally faceted when multiple genes).
plot_one_xy <- function(seu, resolution = NULL, point_size = 0.2, what = c("umap", "spatial"), meta_column = "cluster", suffix = "Spatial.008um",
                        genes = NULL, with_legend = FALSE, label_size = 6, combine_genes = FALSE, discrete_palette = NULL, swap_xy = FALSE)  {
  what <- match.arg(what)

  # Convert resolutions (e.g. "0.1") into metadata column name ("Spatial.008um_snn_res.0.1")
  if(!is.null(resolution))
    res <- get_resolution_columns(resolution, suffix = suffix)

  # Create x-y coordinates for umap or spatial
  if(what == "umap") {
    xy <- FetchData(seu, vars = c("umap_1", "umap_2")) |> 
      as_tibble(rownames = "cell_id") |> 
      rename(x = umap_1, y = umap_2)
    point_shape <- 16
  } else {
    xy <- get_spatial_coordinates(seu, swap_xy)
    corners <- get_resolution_grid(xy)
    point_shape <- 15
  }

  # if gene(s) provided, make an expression plot
  if(!is.null(genes)) {

    ex <- get_gene_expression(seu, genes, norm_per_gene = TRUE)
    if(combine_genes) {
      ex <- ex |> 
        group_by(cell_id) |> 
        summarise(expression = mean(expression)) |> 
        add_column(gene_symbol = "Combined")
    }

    if(what == "umap") {
      d <- inner_join(xy, ex, by = "cell_id") 
    # Plot zeroes first, then expressed genes on top
      d0 <- d |> filter(expression == 0)
      d1 <- d |> filter(expression > 0)
      pl <- d |> 
        ggplot(aes(x = x, y = y, colour = expression)) +
          thm +
          geom_point(data = d0, size = point_size, shape = point_shape, stroke = 0.1) +
          geom_point(data = d1, size = point_size, shape = point_shape, stroke = 0.1) +
          scale_colour_distiller(palette = "Purples", direction = 1) +
          labs(x = NULL, y = NULL)
    } else {
      d <- inner_join(corners, ex, by = "cell_id") 
      pl <- d |> 
        ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = expression)) +
          thm +
          geom_polygon(colour = NA, linewidth = 0) +
          scale_fill_distiller(palette = "Purples", direction = 1) +
          labs(x = NULL, y = NULL)  
    }


   if(length(genes) > 1 & !combine_genes)
      pl <- pl + facet_wrap(~ gene_symbol)
    
    
  # Otherwise make clustering or cell type plot
  } else {

    if(meta_column == "cluster") {
      col <- res
    } else {
      col <- meta_column
    }

    cl <- get_meta_column(seu, col)
    

    if(what == "umap") {
      xy_cl <- inner_join(xy, cl, by = "cell_id")
      pl <- xy_cl |>
      ggplot(aes(x = x, y = y)) +
      thm +
      geom_point(aes(colour = value), size = point_size, shape = point_shape, stroke = 0.1) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      labs(x = NULL, y = NULL)
      if(is.null(discrete_palette)) {
        pl <- pl + scale_colour_brewer(palette = "Paired")
      } else {
        pl <- pl + scale_colour_manual(values = discrete_palette)
      }
      xy_cl_mean <- xy_cl |>
        group_by(value) |> 
        summarise(mx = mean(x), my = mean(y))
      pl <- pl +  
        geom_text(data = xy_cl_mean, aes(label = value, x = mx, y = my), colour = "black", size = label_size)

    } else {
      xy_cl <- inner_join(corners, cl, by = "cell_id")
      pl <- xy_cl |>
        ggplot(aes(x = x_corner, y = y_corner, group = centre_id, fill = value)) +
        thm +
        geom_polygon(colour = NA, linewidth = 0) +
        guides(fill = guide_legend(override.aes = list(size = 3))) +
        labs(x = NULL, y = NULL)  
        if(is.null(discrete_palette)) {
          pl <- pl + scale_fill_brewer(palette = "Paired")
        } else {
          pl <- pl + scale_fill_manual(values = discrete_palette)
        }
  
    }
  }

  if(!with_legend)
      pl <- pl + theme(legend.position = "none")
  
  return(pl)
}


#' Compare multiple clustering resolutions
#'
#' For each resolution, shows UMAP (points) and spatial (tiles) side-by-side.
#'
#' @param seu Seurat object.
#' @param resolutions Numeric/character vector of resolutions.
#' @param umap_point_size Point size for UMAP.
#' @param label_size Label size for cluster text.
#'
#' @return A cowplot grid.
plot_resolutions <- function(seu, resolutions, umap_point_size = 0.2, label_size = 5) {
  pl <- map(resolutions, function(res) {

    pl1 <- plot_one_xy(seu, resolution = res, point_size = umap_point_size, what = "umap", label_size = label_size)
    pl2 <- plot_one_xy(seu, resolution = res, what = "spatial")
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  })

  plot_grid(plotlist = pl, ncol = 1, labels = resolutions)
}

#' Cell type maps (UMAP and spatial)
#'
#' Side-by-side UMAP and spatial plots coloured by cell_type.
#'
#' @param seu Seurat object.
#' @param umap_point_size Point size for UMAP.
#' @param label_size Label size for text labels.
#'
#' @return A cowplot grid.
plot_cell_type <- function(seu, umap_point_size = 0.3, label_size = 3) {
  pal <- c(okabe_ito_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))
  pl1 <- plot_one_xy(seu, point_size = umap_point_size, label_size = label_size, what = "umap", meta_column = "cell_type", discrete_palette = pal, with_legend = FALSE)
  pl2 <- plot_one_xy(seu, label_size = label_size, what = "spatial", meta_column = "cell_type", discrete_palette = pal, with_legend = TRUE) + guides(colour = guide_legend(override.aes = list(size = 3)))
  plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1.4))
}

#' Gene expression grids (UMAP + spatial)
#'
#' For each gene (or combined average), shows UMAP and spatial expression maps.
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols.
#' @param umap_point_size Point size for UMAP plots.
#' @param label_size Label size.
#' @param grid_label_size Facet label size for the grid.
#' @param combine_genes Logical; average genes per cell instead of faceting.
#'
#' @return A cowplot grid or single panel when combine_genes = TRUE.
plot_expression <- function(seu, genes, umap_point_size = 0.1, label_size = 5, grid_label_size = 10,
                            combine_genes = FALSE) {

  plot_side_by_side <- function(gn) {
    pl1 <- plot_one_xy(seu, genes = gn, point_size = umap_point_size, what = "umap", combine_genes = combine_genes, label_size = label_size)
    pl2 <- plot_one_xy(seu, genes = gn, what = "spatial", combine_genes = combine_genes, label_size = label_size)
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  }

  if(combine_genes) {
    plot_side_by_side(genes)
  } else {
    pl <- map(genes, function(gene) {
      plot_side_by_side(gene)
    })
    plot_grid(plotlist = pl, ncol = 1, labels = genes, label_size = grid_label_size, label_y = 0.95)
  }

}

#' Hybrid plot for gene expression by group
#'
#' Dot/quasirandom points per cell atop median line; bottom bars reflect % cells expressed with CI.
#'
#' @param seu Seurat object.
#' @param gids Character vector of gene symbols.
#' @param layer Data layer ("data" or "counts").
#' @param meta_column Metadata column defining groups (default "cluster").
#' @param width Quasirandom jitter width.
#' @param ncol Number of columns in facet layout (NULL for auto).
#' @param div Fraction of y-range used for percent bars below zero.
#' @param col_width Half-width of summary bars.
#'
#' @return A faceted ggplot.
plot_gene_expression <- function(seu, gids, layer = "data", meta_column = "cluster", width = 0.3, ncol = NULL,
                                 div = 0.3, col_width = 0.4) {
  
  # Necessary to avoid keeping large (~1 GB) environment within ggplot object
  env <- new.env(parent = globalenv())
  env$width <- width
  env$ncol <- ncol
  
  mt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id")
  
  env$d <- get_gene_expression(seu, gids, layer = layer) |> 
    left_join(mt, by = "cell_id") |> 
    mutate(
      gene_symbol = factor(gene_symbol, levels = gids),
      grp = get(meta_column),
      igrp = as.integer(grp)
    )
  env$groups <- levels(env$d$grp)
  env$m <- env$d |> 
    group_by(gene_symbol, igrp) |> 
    summarise(
      n_tot = n(),
      mn = mean(expression[expression > 0], na.rm = TRUE),
      md = median(expression[expression > 0], na.rm = TRUE),
      n_expr = length(which(expression > 0)),
      prop = n_expr / n_tot
    ) |> 
    ungroup() |> 
    nest(data = c(n_expr, n_tot)) |>
    mutate(res = purrr::map(data, ~PropCIs::exactci(.x$n_expr, .x$n_tot, conf.level = 0.95) |> broom::glance())) |>
    unnest(c(data, res)) |> 
    mutate(gene_symbol = factor(gene_symbol, levels = gids))
  
  env$max_y <- max(env$d$expression)
  env$min_y <- -div * env$max_y

  env$m <- env$m |>
    mutate(
      y = (1 - prop) * env$min_y,
      y_lo = (1 - conf.low) * env$min_y,
      y_up = (1 - conf.high) * env$min_y
    )
  
  env$d <- env$d |> 
    filter(expression > 0)
  
  env$wid <- col_width
  with(env, {
    ggplot(d, aes(x = igrp)) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      geom_quasirandom(aes(y = expression), size = 0.5, colour = "grey60", width = width) +
      geom_segment(data = m, aes(x = igrp - wid, xend = igrp + wid, y = md, yend = md), linewidth = 1) +
      geom_col(data = m, aes(y = min_y), fill = "brown", width = wid, alpha = 0.6) +
      geom_col(data = m, aes(y = y), fill = "grey80", width = wid) +
      geom_errorbar(data = m, aes(ymin = y_lo, ymax = y_up), width = wid / 2) +
      facet_wrap( ~ gene_symbol, ncol = ncol) +
      geom_hline(yintercept = 0, colour = "grey70") +
      scale_x_continuous(breaks = seq_along(groups), labels = groups) +
      scale_y_continuous(limits = c(min_y, max_y), expand = expansion(mult = c(0, 0.03))) +
      labs(x = NULL, y = "Value")
  })
}



#' Heatmap of average expression per cluster
#'
#' Computes AverageExpression() per cluster and shows gene-wise relative expression (z/mean scaling).
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param resolution Optional clustering resolution to set identities before averaging.
#'
#' @return A ggplot heatmap.
plot_gene_heatmap <- function(seu, gns, resolution = NULL) {
  env <- new.env(parent = globalenv())
  
  if(!is.null(resolution)) {
    res <- get_resolution_columns(resolution)
    Idents(seu) <- res  
  }

  d <- AverageExpression(seu, layer = "data", features = gns, group.by = "ident") |> 
    pluck("Spatial.008um") |> 
    as_tibble(rownames = "gene_symbol") |> 
    pivot_longer(-gene_symbol, names_to = "cluster") |> 
    group_by(gene_symbol) |> 
    mutate(
      z = value / mean(value)
    ) |> 
    ungroup()
  
  env$tab <- d |>
    pivot_wider(id_cols = gene_symbol, names_from = cluster, values_from = z) |>
    column_to_rownames("gene_symbol")

  with(env, {
    ggheatmap(tab, order.col = FALSE, with.y.text = TRUE, palette = "viridis",
              legend.name = "rel expr")  
  })
}

#' Dot plot of gene-by-group expression
#'
#' Wrapper around Seurat DotPlot() with cleaned aesthetics and optional legend suppression.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param meta_column Grouping metadata column.
#' @param with_legend Logical; include legend.
#' @param title Optional plot title.
#'
#' @return A ggplot dot plot.
plot_gene_dots <- function(seu, gns, meta_column = "cluster", with_legend = TRUE, title = NULL) {
  env <- new.env(parent = globalenv())
  env$meta_column <- meta_column
  env$with_legend <- with_legend
  env$title <- title

  Idents(seu) <- meta_column
  dp <- DotPlot(seu, features = gns)
  env$d <- dp$data |> 
    as_tibble(rownames = "gene_symbol")
   
  with(env, {
    pl <- ggplot(d, aes(x = id, y = fct_rev(features.plot), colour = avg.exp.scaled, size = pct.exp)) +
      th +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_point() +
      scale_colour_gradient(low = "grey90", high = "blue") +
      scale_radius(range = c(1, 6)) +
      labs(x = meta_column, y = "Gene", size = "Percent\nexpressed", colour = "Average\nexpresion", title = title)
    if(!with_legend)
      pl <- pl + theme(legend.position = "none")
    pl
  })
}

#' Dot plots of top markers per cell type
#'
#' For each cluster in cell_type_markers, shows dot plot of top n_top genes across cell_type.
#'
#' @param seu Seurat object.
#' @param cell_type_markers Tibble with cluster, gene_symbol, and avg_log2FC.
#' @param n_top Number of markers per cell type.
#'
#' @return A cowplot grid of dot plots.
plot_cell_type_marker_dots <- function(seu, cell_type_markers, n_top = 30) {
  cell_types <- cell_type_markers$cluster |> unique() |> as.character()
  pl <- map(cell_types, function(ct) {
    gns <- cell_type_markers |> 
      filter(cluster == ct) |> 
      slice_max(avg_log2FC, n = n_top) |> 
      arrange(desc(avg_log2FC)) |> 
      pull(gene_symbol)
    plot_gene_dots(seu, gns, meta_column = "cell_type", with_legend = FALSE, title = ct)
  })
  plot_grid(plotlist = pl, ncol = 2)
}


#' Save UMAP+spatial expression grids for a gene list
#'
#' Produces a PNG in doc/images for the given genes.
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols.
#' @param file_name Base output name.
#' @param width Plot width in inches.
#'
#' @return Invisibly, the file path written.
plot_png_expression_genes <- function(seu, genes, file_name, width = 8) {
  plot_expression(seu, genes, umap_point_size = 0.4, grid_label_size = 10) |> 
    gs(file_name, width = width, height = 4 * length(genes))
}

#' Save key marker panels per cell type
#'
#' For each cell type in key_markers, saves UMAP+spatial expression PNG for its genes.
#'
#' @param seu Seurat object.
#' @param key_markers Tibble with columns cell_type and gene_symbol.
#' @param file_prefix Prefix for output files.
#' @param width Plot width in inches.
#'
#' @return A list of file paths (invisibly).
plot_key_markers_umap_spatial <- function(seu, key_markers, file_prefix, width = 8) {
  all_feat <- Features(seu)

  km <- key_markers |> 
    filter(gene_symbol %in% all_feat)
  nums <- key_markers |> 
    count(cell_type)

  map(nums$cell_type, function(ct) {
    gname <- str_glue("{file_prefix}_{ct}")
    gns <- km |> 
      filter(cell_type == ct) |> 
      pull(gene_symbol) |> 
      as.character()
    plot_png_expression_genes(seu, gns, gname, width)
  }) 
}

#' Save top cluster marker panels
#'
#' For each cluster, saves UMAP+spatial expression PNG of its top markers.
#'
#' @param seu Seurat object.
#' @param top_cluster_markers Tibble with cluster and gene_symbol.
#' @param file_prefix Output file prefix.
#'
#' @return A list of file paths (invisibly).
plot_top_cluster_umap_spatial <- function(seu, top_cluster_markers, file_prefix = "top_clust_gene") {
  
  clusters <- top_cluster_markers$cluster |> 
    unique() |> 
    sort()
  
  map(clusters, function(cl) {
    gname <- sprintf("%s_cl%02d", file_prefix, cl)
    gns <- top_cluster_markers |> 
      filter(cluster == cl) |> 
      pull(gene_symbol)
    plot_png_expression_genes(seu, gns, gname)
  })
}

#' Correlation heatmap between genes
#'
#' Filters a correlation table by FDR and/or correlation, selects top genes, and draws a symmetric heatmap.
#'
#' @param tbl_cor Tibble with gene_1, gene_2, correlation, p_adj.
#' @param fdr_limit Adjusted p-value cutoff.
#' @param min_cor Minimum absolute correlation to include.
#' @param n_top Optional top-N genes to show.
#' @param cluster_rows Logical; cluster rows.
#' @param cluster_columns Logical; cluster columns.
#'
#' @return A ComplexHeatmap heatmap.
plot_cor_heat <- function(tbl_cor, fdr_limit = 0.01, min_cor = NULL, n_top = NULL, cluster_rows = TRUE, cluster_columns = TRUE) {
  env <- new.env(parent = globalenv())

  top_genes <- find_top_cor_genes(tbl_cor, fdr_limit, min_cor, n_top)

  # Make matrix symmetric for the heatmap plot
  tb <- symmetrise_table(tbl_cor)

  if(length(top_genes) == 0)
    return(NULL)
  tab <- tb |> 
    filter(gene_2 %in% top_genes & p_adj < fdr_limit) |> 
    pivot_wider(id_cols = gene_2, names_from = gene_1, values_from = correlation, values_fill = NA) |> 
    column_to_rownames("gene_2") |> 
    as.matrix()
  
  env$tab <- tab[, intersect(levels(tbl_cor$gene_1), colnames(tab))]
  env$palet <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#FFFFFF", "#B2182B"), space = "LAB")
  env$cluster_columns <- cluster_columns
  env$cluster_rows = cluster_rows

  with(env, {
    ComplexHeatmap::Heatmap(tab, col = palet, name = "Correlation", cluster_columns = cluster_columns, cluster_rows = cluster_rows)
  })
  #ggheatmap(tab, palette = "distiller", with.y.text = TRUE, legend.name = "Correlation")
}

#' Mean expression of gene sets vs gene sets
#'
#' Compares average expression of two sets (x,y) per cell and colours by a grouping variable.
#'
#' @param seu Seurat object.
#' @param genes1 Character vector for x-axis mean.
#' @param genes2 Character vector for y-axis mean.
#' @param group_var Metadata column for colour (default "cluster").
#' @param group_sel Optional subset of groups to highlight; others collapsed to "Other".
#' @param point_size Point size.
#'
#' @return A ggplot scatter.
plot_gene_vs_gene <- function(seu, genes1, genes2, group_var = "cluster", group_sel = NULL, point_size = 0.5) {

  mt <- get_cell_data(seu) |> 
    mutate(group = get(group_var)) |> 
    select(cell_id, group)

  palette <- okabe_ito_palette
  if(!is.null(group_sel)) {
    mt <- mt |> 
      mutate(
        group = as.character(group),
        group = if_else(group %in% group_sel, group, "Other"),
        group = fct_relevel(group, "Other")
      ) 
    palette <- c("grey90", okabe_ito_palette)
  }

  d1 <- get_gene_expression(seu, genes1, norm_per_gene = FALSE) |> 
    summarise(gr1 = mean(expression), .by = cell_id)
  d2 <- get_gene_expression(seu, genes2, norm_per_gene = FALSE) |> 
    summarise(gr2 = mean(expression), .by = cell_id)

  inner_join(d1, d2, by = "cell_id") |> 
    left_join(mt, by = "cell_id") |> 
  ggplot(aes(x = gr1, y = gr2, colour = group)) +
    th +
    geom_point(size = point_size) +
    scale_colour_manual(values = palette) +
    labs(x = str_c(genes1, collapse = ","), y = str_c(genes2, collapse = ","), colour = group_var)
}


#' PCA of cell-type marker matrix
#'
#' Runs PCA on group-by-gene average expression (rows: genes, columns: cell types).
#'
#' @param seu Seurat object.
#' @param select_markers Tibble with gene_symbol and grouping column "cell_type".
#'
#' @return A ggplot of PCA1 vs PCA2.
plot_cell_marker_pca <- function(seu, select_markers) {

  env <- new.env(parent = globalenv())

  tab <- get_group_gene_expression(seu, select_markers, group_var = "cell_type") |> 
    select(-cell_id) |> 
    as.matrix() |> 
    t()

  tab <- tab[apply(tab, 1, function(v) sum(v) > 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]

  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  
  env$pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  env$pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])

  env$xy <- pca$x |>
    as_tibble(rownames = "id") |>
    select(x = PC1, y = PC2, id)

  with(env, {
    xy |> 
      ggplot(aes(x = x, y = y)) +
      th +
      geom_point() +
      labs(x = pca1, y = pca2)
  })
}

#' List of sets for UpSet from expression matrix
#'
#' Converts a wide matrix of 0/1 (or probability) per cell to a list of sets crossing a threshold.
#'
#' @param mex Tibble with cell_id and columns for each set (numeric).
#'
#' @return A named list of character vectors (cell IDs per set).
make_list_for_upset <- function(mex) {
  d <- mex |> 
    pivot_longer(-cell_id) |> 
    filter(value > 0.5)
  groups <- unique(d$name)

  map(groups, function(gr) {
    d |> 
      filter(name == gr) |> 
      pull(cell_id)
  }) |> 
    set_names(groups)
}

#' Combination matrix for UpSet
#'
#' Creates a ComplexHeatmap::comb_mat from a wide 0/1 matrix of membership.
#'
#' @param mex Tibble with cell_id and set columns.
#' @param … Passed to ComplexHeatmap::make_comb_mat().
#'
#' @return A comb_mat object.
make_combmat_for_upset <- function(mex, ...) {
  mex |> 
    make_list_for_upset() |> 
    ComplexHeatmap::make_comb_mat(...)
}

# cnt_tab needs column what, "sample" and "n".
plot_counted_cell_perc <- function(cnt_tab, what) {
  cnt_tab |> 
    mutate(grp = get(what)) |> 
    group_by(sample) |> 
    mutate(perc = 100 *n / sum(n)) |> 
    ggplot(aes(x = grp, y = perc, fill = sample)) +
    th +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = what, y = "Percentage of cells")
}

#' Jaccard heatmap: clusters vs cell types
#'
#' Computes and plots Jaccard indices of overlap between clusters and annotated cell types.
#'
#' @param seu Seurat object.
#'
#' @return A ggplot tile heatmap.
plot_cluster_cell_type_jaccard <- function(seu) {
  env <- new.env(parent = globalenv())

  env$d <- seu |>
    jaccard_clusters_cell_types()
  with(env, {
    d |> 
      ggplot(aes(x = cluster, y = cell_type, fill = j)) +
      th +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis") +
      labs(fill = "Jaccard")
  })
}

#' Sankey diagram: clusters to cell types
#'
#' Alluvial flow from annotated cell_type to clusters at the given resolution.
#'
#' @param seu Seurat object.
#' @param resolution Resolution string used to fetch cluster column.
#'
#' @return A ggplot Sankey diagram.
plot_cluster_cell_type_sankey <- function(seu, resolution = "0.4") {
  env <- new.env(parent = globalenv())

  res <- get_resolution_columns(resolution)

  d <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, cell_type, cluster := !!res) |> 
    mutate(across(where(is.factor), as.character)) |> 
    make_long(cell_type, cluster) |> 
    mutate(
      node = factor(node, levels = c("PP", "Epsilon", "Delta", "Beta", "Alpha", "Ductal", "Acinar", "MPC", "Ventral MP", "Dorsal MP",
      "EP", "Immune", "Blood", "VSM", "Nervous", "Matrix", "Unknown", "1", "2", "6", "5", "3", "10", "9", "8", "4", "7"))
    )

  env <- new.env(parent = globalenv())
  env$d <- d

  with(env, {
    ggplot(d, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
    theme_sankey() +
    theme(
      legend.position = "none"
    ) +
    geom_sankey(flow.alpha = 0.6, flow.colour = "grey70", flow.fill = "grey80", node.alpha = 0.6) +
    geom_sankey_label() +
    #scale_fill_manual(values = c(okabe_ito_palette[1:4], cell_col$colour),drop = FALSE) +
    labs(x = NULL)
  })
}

#' Expression vs pseudo-time for genes
#'
#' Points of expression vs a pseudo-time metadata column, faceted by gene.
#'
#' @param seu Seurat object.
#' @param gids Character vector of genes.
#' @param pt_column Metadata column name containing pseudo-time.
#' @param point_size Point size.
#'
#' @return A faceted ggplot.
plot_expression_pseudotime <- function(seu, gids, pt_column = "pseudotime_beta", point_size = 0.8) {
  ex <- get_gene_expression(seu, gids) |> 
    filter(expression > 0)
  pt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, pseudo_time = all_of(pt_column)) |> 
    drop_na()
  inner_join(ex, pt, by = "cell_id") |> 
    ggplot(aes(x = pseudo_time, y = expression)) +
    th +
    geom_point(size = point_size) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    labs(x = "Pseudo-time", y = "Normalised expression")
}

#' Internal helper for MA/volcano plotting
#'
#' Draws group-faceted scatter with selected points highlighted.
#'
#' @param res Tibble with columns x, y, sel, and grouping column.
#' @param group_var Grouping column name to facet by.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A ggplot object.
plot_volma <- function(res, group_var, point_size, point_alpha) {
  r <- res |>
    mutate(
      group = get(group_var)
    ) |>
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)

  rm(res, r)  # Minimise environment for serialisation

  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group)

  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

#' MA plot by contrast (faceted)
#'
#' MA plot of average expression vs log2FC with significant points highlighted.
#'
#' @param res Differential results tibble.
#' @param a Column name for average expression (x-axis).
#' @param fc Column name for log2 fold change (y-axis).
#' @param p Column name for p-value.
#' @param fdr Column name for adjusted p-value.
#' @param group_var Grouping column for facets.
#' @param fdr_limit Significance cutoff.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A ggplot object.
plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group_var = "contrast",
                    fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(a),
      y = get(fc),
      sel = get(fdr) < fdr_limit
    ) |>
    plot_volma(group_var, point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

#' Volcano plot by contrast (faceted)
#'
#' Volcano plot of log2FC vs -log10(p) with significant points highlighted.
#'
#' @param res Differential results tibble.
#' @param fc Column name for log2 fold change (x-axis).
#' @param p Column name for p-value.
#' @param fdr Column name for adjusted p-value.
#' @param group_var Grouping column for facets.
#' @param fdr_limit Significance cutoff.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A ggplot object.
plot_volcano <- function(res, fc = "avg_log2FC", p = "p_val", fdr = "p_val_adj", group_var = "contrast",
                         fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      pval = get(p),
      pval = if_else(pval == 0, .Machine$double.xmin, pval),
      x = get(fc),
      y = -log10(pval),
      sel = get(fdr) < fdr_limit
    ) |>
    plot_volma(group_var, point_size, point_alpha) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

#' P-value distribution histograms
#'
#' Faceted histograms of p-values to assess uniformity/enrichment at 0.
#'
#' @param res Differential results tibble.
#' @param p Column name for p-value.
#' @param group_var Grouping column for facets.
#' @param n_bins Number of bins.
#'
#' @return A ggplot histogram.
plot_pdist <- function(res, p = "PValue", group_var = "contrast", n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |>
    mutate(p = get(p), grp = get(group_var)) |>
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}

#' VMP: Volcano, MA, and P-value panels
#'
#' Convenience function returning all three diagnostic plots.
#'
#' @param res Differential results tibble.
#' @param a Column name for average expression.
#' @param fc Column name for log2 fold change.
#' @param p Column name for p-value.
#' @param fdr Column name for adjusted p-value.
#' @param group Column name for grouping (facets).
#' @param fdr_limit Significance threshold.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param n_bins Number of histogram bins.
#'
#' @return A named list with volcano, ma, and pdist ggplots.
plot_vmp <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                     fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5, n_bins = 50){
  list(
    volcano = plot_volcano(res, fc, p, fdr, group, fdr_limit, point_size, point_alpha),
    ma = plot_ma(res, a, fc, p, fdr, group, fdr_limit, point_size, point_alpha),
    pdist = plot_pdist(res, p, group, n_bins)
  )
}

#' Spatial maps of cell-type scores
#'
#' Plots enrichment/score meta-columns (ES_<cell_type>) for each tissue/cell type.
#'
#' @param seu Seurat object.
#' @param tissue_markers Tibble with cell_type column.
#'
#' @return A cowplot grid.
plot_cell_type_scores <- function(seu, tissue_markers) {
  cell_types <- tissue_markers$cell_type |> unique()
  pls <- map(cell_types, function(ct) {
    ct <- str_glue("ES_{ct}")
    plot_spatial_meta(seu, ct)
  })
  plot_grid(plotlist = pls, ncol = 2)
}

#' Spatial overlap of two genes
#'
#' Classifies each spot into four states (gene1 only, gene2 only, both, none) and plots spatial tiles.
#'
#' @param seu Seurat object.
#' @param gene_1 First gene symbol.
#' @param gene_2 Second gene symbol.
#' @param section Optional tibble of a polyline to overlay with columns x,y.
#'
#' @return A ggplot polygon map.
plot_gene_overlap <- function(seu, gene_1, gene_2, section = NULL) {
  only_1 <- str_glue("{gene_1} only")
  only_2 <- str_glue("{gene_2} only")

  xy <- get_spatial_coordinates(seu, swap_xy = FALSE)
  corners <- get_resolution_grid(xy)
  ex <- get_gene_expression(seu, c(gene_1, gene_2), layer = "count") |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression) |> 
    rename(g1 := !!gene_1, g2 := !!gene_2) |> 
    mutate(
      state = case_when(
        g1 > 0 & g2 > 0 ~ "Both",
        g1 == 0 & g2 > 0 ~ only_2,
        g1 > 0 & g2 == 0 ~ only_1,
        g1 == 0 & g2 == 0 ~ "None"
      ) |> 
        factor(levels = c(only_1, only_2, "Both", "None"))
    )
    
  pl <- inner_join(corners, ex, by = "cell_id") |> 
    ggplot(aes(x = x_corner, y = y_corner)) +
    thm +
    geom_polygon(aes(group = cell_id, fill = state), colour = NA, linewidth = 0) +
    scale_fill_manual(values = c(okabe_ito_palette[c(2, 3, 1)], "grey90")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = NULL, y = NULL, colour = "Gene expression")

  if(!is.null(section)) {
    pl <- pl + 
      geom_path(data = section, aes(x = x, y = y), colour = "black", linewidth = 1)
  }

  pl
}

#' Distance of points to a line segment
#'
#' Computes projection parameter and perpendicular distance from (x,y) points to a segment.
#'
#' @param xy Tibble with columns x, y.
#' @param segment Tibble with two rows and columns x, y.
#'
#' @return Tibble with added t, proj_x, proj_y, and dst.
dist_to_segment <- function(xy, segment) {
  x1 <- segment$x[1]
  x2 <- segment$x[2]
  y1 <- segment$y[1]
  y2 <- segment$y[2]
  dx <- x2 - x1
  dy <- y2 - y1

  xy |> 
    mutate(
      # projection parameter t
      t = ((x - x1) * dx + (y - y1) * dy) / (dx^2 + dy^2),
      # projection point
      proj_x = x1 + t * dx,
      proj_y = y1 + t * dy,
      # Euclidean distance
      dst = sqrt((x - proj_x)^2 + (y - proj_y)^2)
    )
}


#' Expression along a spatial section
#'
#' Extracts spots near a line segment and plots expression vs segment parameter for two genes.
#'
#' @param seu Seurat object.
#' @param gene_1 First gene symbol.
#' @param gene_2 Second gene symbol.
#' @param section Tibble with two points (x,y) defining a segment.
#' @param layer Data layer ("data" or "count").
#' @param slice Image key used to infer spot radius.
#'
#' @return A ggplot line plot (two series).
plot_expression_section <- function(seu, gene_1, gene_2, section, layer = "data", slice = "slice1.008um") {
  
  rad <- seu@images[[slice]]@boundaries$centroids@radius
  ex <- get_gene_expression(seu, c(gene_1, gene_2), layer = layer) |> 
    mutate(gene_symbol = factor(gene_symbol, levels = c(gene_1, gene_2)))
  ylab <- ifelse(layer == "count", "UMI count", "log1p normalised count")

  pl <- seu |> 
    get_spatial_coordinates() |> 
    dist_to_segment(section) |> 
    filter(dst < rad / 2 & t >=0 & t <= 1) |> 
    arrange(t) |> 
    left_join(ex, by = "cell_id") |> 
  ggplot(aes(x = t, y = expression, colour = gene_symbol)) +
    th +
    geom_line(, linewidth = 1.3) +
    scale_colour_manual(values = okabe_ito_palette[c(2, 3)]) +
    labs(x = "t", y = ylab, colour = "Gene")

  if(layer == "count")
    pl <- pl + geom_hline(yintercept = 0, colour = "brown", alpha = 0.3)
  pl
}

#' Spatial overlap + section profile
#'
#' Two-panel layout: spatial overlap map and line-profile of expression along the section.
#'
#' @param seu Seurat object.
#' @param gene_1 First gene symbol.
#' @param gene_2 Second gene symbol.
#' @param section Tibble with two points (x,y) defining a segment.
#' @param layer Data layer passed to plot_expression_section().
#' @param slice Image key passed to plot_expression_section().
#'
#' @return A cowplot grid.
plot_spatial_section <- function(seu, gene_1, gene_2, section, layer = "count", slice = "slice1.008um") {
  pl1 <- plot_gene_overlap(seu, gene_1, gene_2, section)
  pl2 <- plot_expression_section(seu, gene_1, gene_2, section, layer = layer, slice = slice)

  plot_grid(pl1, pl2, ncol = 1, rel_heights = c(2, 1), align = "v")
}

#' Spatial difference between two meta scores
#'
#' Adds type_1 - type_2 as a new meta column and plots it with a divergent palette.
#'
#' @param seu Seurat object.
#' @param type_1 Name of first meta column.
#' @param type_2 Name of second meta column.
#'
#' @return A ggplot polygon map.
plot_spatial_difference <- function(seu, type_1, type_2) {
  newcol <- str_glue("{type_1} - {type_2}")
  seu[[newcol]] <- seu[[type_1]] - seu[[type_2]]
  plot_spatial_meta(seu, meta_column = newcol, palette = "divergent") +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.4, 'cm')
    )
}

# pairs is a list of vectors

#' Batch spatial differences for pairs
#'
#' For each pair of names, plots the difference of corresponding ES_ meta-columns.
#'
#' @param seu Seurat object.
#' @param pairs List of length-2 character vectors, e.g., list(c("Alpha","Beta")).
#'
#' @return A cowplot grid.
plot_spatial_differences <- function(seu, pairs) {
  pl <- map(pairs, function(pair) {
    ct1 <- paste0("ES_", pair[1])
    ct2 <- paste0("ES_", pair[2])
    plot_spatial_difference(seu, ct1, ct2)
  })

  plot_grid(plotlist = pl, ncol = 1)
}

#' UpSet plot of co-expressed cells
#'
#' Builds an UpSet diagram of cells expressing any of the given genes above a threshold.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param expr_limit Numeric threshold on expression to call expressed.
#'
#' @return A ComplexHeatmap UpSet object.
plot_upset_expressed_cells <- function(seu, gns, expr_limit = 0, layer = "data") {
  m <- get_gene_expression(seu, gns, layer = layer) |> 
    filter(expression > expr_limit) |> 
    mutate(expression = 1) |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression, values_fill = 0) |> 
    ComplexHeatmap::make_comb_mat()
  
  ComplexHeatmap::UpSet(
    m,
    top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE, numbers_rot = 0),
    right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE),
    comb_col = c("black", "brown", "darkblue", "darkgreen", "darkorchid3")[ComplexHeatmap::comb_degree(m)]
  )
}

#' Save comb_mat for UpSet to RDS
#'
#' Computes a combination matrix and writes it to doc/tables/<suffix>_<sample>.rds.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param suffix String identifier for output.
#' @param sample Sample identifier.
#' @param expr_limit Numeric threshold on expression to call expressed.
#'
#' @return Invisibly, the RDS path written.
combat_upset_expressed_cells <- function(seu, gns, suffix, sample, expr_limit = 0) {
  rfile <- file.path("doc", "tables", str_glue("{suffix}_{sample}.rds"))
  m <- get_gene_expression(seu, gns, layer = "data") |> 
    filter(expression > expr_limit) |> 
    mutate(expression = 1) |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression, values_fill = 0) |> 
    ComplexHeatmap::make_comb_mat()
  readr::write_rds(m, rfile)
}

#' Save PNG UpSet diagram
#'
#' Renders UpSet as a PNG in doc/images/<suffix>_<sample>.png.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param suffix String identifier for output.
#' @param sample Sample identifier.
#' @param width Image width in pixels.
#' @param height Image height in pixels.
#' @param expr_limit Numeric expression threshold.
#'
#' @return Invisibly, the PNG path written.
png_upset_expressed_cells <- function(seu, gns, suffix, sample, width, height, expr_limit = 0) {
  pfile <- file.path("doc", "images", str_glue("{suffix}_{sample}.png"))

  png(pfile, width = width, height = height)

  plot_upset_expressed_cells(seu, gns, expr_limit)

  dev.off()
}

#' Spatial scatter with selection box
#'
#' Spatial point map of a gene with a square box overlay at (x_cen, y_cen) with half-size size.
#'
#' @param seu Seurat object.
#' @param gene Gene symbol.
#' @param x_cen X centre of box.
#' @param y_cen Y centre of box.
#' @param size Half-width/height of box (in coordinate units).
#' @param point_size Point size.
#'
#' @return A ggplot object.
plot_spatial_box <- function(seu, gene, x_cen, y_cen, size, point_size = 0.6) {
  xy <- get_spatial_coordinates(seu, swap_xy = FALSE)
  corners <- get_resolution_grid(xy)
  ex <- get_gene_expression(seu, gene)

  bx <- tibble(
    x = c(x_cen - size, x_cen + size, x_cen + size, x_cen - size, x_cen - size),
    y = c(y_cen - size, y_cen - size, y_cen + size, y_cen + size, y_cen - size)
  )

  corners |> 
    inner_join(ex, by = "cell_id") |> 
    ggplot() +
    th +
    theme(legend.position = "none") +
    geom_polygon(aes(group = cell_id, x = x_corner, y = y_corner, fill = expression), colour = NA, linewidth = 0) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    geom_path(data = bx, aes(x = x, y = y)) +
    labs(x = "x", y = "y")
}

#' Compare 8µm vs 2µm resolution spatial tiling
#'
#' Plots polygons for 8µm and 2µm spot centres with their respective corner grids.
#'
#' @param res_comp Tibble of spatial centres with columns cell_id, x, y, and resolution in {“8um”,“2um”}.
#'
#' @return A ggplot showing both resolutions.
plot_resolution_comparison <- function(res_comp) {
  rg <- res_comp |> 
    filter(resolution == "8um") |> 
    get_resolution_grid()

  # Add fifth corner for completion
  rg5 <- rg |> 
    filter(corner_idx == 1) |> 
    mutate(corner_idx = 5)
  rg <- bind_rows(rg, rg5)

  rg <- rg |> 
    arrange(centre_id, corner_idx)

  res_comp |> 
    ggplot() +
    th +
    theme(legend.position = "none") +
    geom_text(aes(x = x, y = y, colour = resolution, label = expression)) +
    geom_path(data = rg, aes(x = x_corner, y = y_corner, group = centre_id)) +
    facet_wrap(~gene_symbol) +
    scale_colour_manual(values = c("grey70", "black"))
}


#' Compute square polygon corners for each spot
#'
#' Given centre coordinates on a regular grid, computes four corner coordinates per spot (plus indexing).
#'
#' @param res_comp Tibble with columns cell_id, x, y; may include centre_id.
#'
#' @return Tibble with cell_id, centre_id, corner_idx, x_corner, and y_corner.
get_resolution_grid <- function(res_comp) {
  df <- res_comp |> 
    mutate(centre_id = cell_id |> as.factor() |> as.integer()) |> 
    distinct(cell_id, centre_id, x, y)

  # This part was written by ChatGPT o4-mini-high and improved by MG:

  knn_res <- FNN::get.knn(df[, c("x","y")], k = 2)
  nearest_idx <- knn_res$nn.index[1, 2]   # the second‐closest to row 1
  nearest_dist <- knn_res$nn.dist[1, 2]   # the distance to that neighbour

  # Extract the vector from centre[1,] to its nearest neighbour:
  v_x <- df$x[nearest_idx] - df$x[1]
  v_y <- df$y[nearest_idx] - df$y[1]

  # (3) Side‐length = s = distance between two adjacent centres:
  s <- sqrt(v_x^2 + v_y^2)

  # Form a unit‐vector along that grid direction:
  u1 <- c(v_x, v_y) / s    # length‐1
  # A perpendicular unit‐vector (rotate u1 by +90° CCW):
  u2 <- c(-u1[2], u1[1])   # also length‐1

  # (4) Corner‐offsets = ±(s/2)*u1 ± (s/2)*u2.  Define the two half‐offsets:
  half1 <- (s/2) * u1      # half along axis‐1
  half2 <- (s/2) * u2      # half along axis‐2

  # Precompute the four combinations of ( ±half1  ±half2 ):
  corner_offsets <- rbind(
    half1 + half2,
    half1 - half2,
    -half1 - half2,
    -half1 + half2
  )
  # corner_offsets is a 4×2 matrix; row 1 = +half1 +half2, etc.

  map(1:4, function(k) {
    df |> 
      mutate(
        corner_idx = k,
        x_corner = x +  corner_offsets[k, 1],
        y_corner = y + corner_offsets[k, 2]
      )
  }) |> 
    list_rbind() |> 
    group_by(centre_id) |>     
    arrange(cell_id, corner_idx) |> 
    ungroup()
}

#' Cluster composition for a single cell type
#'
#' Bar plot of the proportion of one annotated cell_type within each clustering label.
#'
#' @param seu Seurat object.
#' @param ct Cell type of interest (string).
#'
#' @return A ggplot column chart.
plot_cluster_cell_type_overlap <- function(seu, ct = "Unknown") {
  ov <- cluster_cell_type_overlap(seu)

  ov |> 
    filter(cell_type == ct) |> 
    ggplot(aes(x = cluster, y = prop)) +
    th +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Cluster", y = "Proporion of spots in cluster", title = ct)
}

#' Spatial binning by a gene's expression
#'
#' Bins expression of gene_x into intervals; assigns a bin meta and plots spatial polygons by bin.
#'
#' @param seu Seurat object.
#' @param gene_x Gene symbol used for binning.
#' @param bins Number of equal-width intervals.
#'
#' @return A ggplot polygon map with legend.
plot_binned_spatial <- function(seu, gene_x, bins = 6) {
  gx <- get_gene_expression(seu, gene_x)
  g1 <- gx |> 
    filter(expression > 0) |> 
    mutate(bin = cut_interval(expression, n = bins) |> fct_expand("0", after = 0))
  g0 <- gx |> 
    filter(expression == 0) |> 
    mutate(bin = factor("0", levels = levels(g1$bin)))
  gx <- bind_rows(g0, g1)
  # reorder before inserting into meta
  m <- tibble(cell_id = colnames(seu)) |> 
    left_join(gx, by = "cell_id")
  seu$bin <- m$bin

  plot_one_xy(seu, what = "spatial", meta_column = "bin") + theme(legend.position = "right")
}


#' Expression by bins (dot/quasi plot)
#'
#' Computes bins by gene_x and shows expression of genes_y across bins via plot_gene_expression().
#'
#' @param seu Seurat object.
#' @param gene_x Gene symbol used for binning.
#' @param genes_y Character vector of genes to plot by bin.
#' @param bins Number of bins.
#'
#' @return A faceted ggplot.
plot_binned_expression <- function(seu, gene_x, genes_y, bins = 6) {
  b <- gene_binned_gradient(seu, gene_x, bins)
  seu$bin <- b$bin

  plot_gene_expression(seu, genes_y, meta_column = "bin") + ylab("Expression")
}

#' Summary statistics of expression across bins (obsolete)
#'
#' For each gene in genes_y, summarises mean and 95% CI across bins of gene_x.
#'
#' @param seu Seurat object.
#' @param gene_x Gene symbol used for binning.
#' @param genes_y Character vector of genes to summarise.
#' @param what One of "expression" (mean) or "pct" (percentage expressed).
#' @param bins Number of bins.
#' @param include_x_zeroes Logical; include zero-expression cells when binning.
#' @param include_y_zeroes Logical; include zero-expression cells in summaries.
#'
#' @return A faceted ggplot with points and error bars.
plot_binned_xy_expression <- function(seu, gene_x, genes_y, what = "expression", bins = 6,
                                      include_x_zeroes = FALSE, include_y_zeroes = FALSE) {
  env <- new.env(parent = globalenv())
  env$gene_x <- gene_x

  gx <- gene_binned_gradient(seu, gene_x, bins, include_x_zeroes)

  gy <- get_gene_expression(seu, genes_y)
  if(!include_y_zeroes & what == "expression") {
    gy <- gy |> 
      filter(expression > 0)
  }

  d <- gy |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes_y)) |> 
    inner_join(gx, by = "cell_id")  |> 
    group_by(gene_symbol, bin) |> 
    summarise(
      n_tot = n(),
      m = mean(expression),
      se = sd(expression) / sqrt(n()),
      ci = se * qt(0.975, n_tot - 1),
      n_expr = length(which(expression > 0)),
      prop = n_expr / n_tot
    ) |> 
    ungroup() |> 
    nest(data = c(n_expr, n_tot)) |>
    mutate(res = purrr::map(data, ~PropCIs::exactci(.x$n_expr, .x$n_tot, conf.level = 0.95) |> broom::glance())) |>
    unnest(c(data, res))

  if(what == "expression") {
    env$d <- d |> 
      mutate(ymin = m - ci, ymax = m + ci) |> 
      select(gene_symbol, bin, y = m, ymin, ymax)
    env$ylab <- "Mean expression"
  } else if(what == "proportion") {
    env$d <- d |> 
      select(gene_symbol, bin, y = prop, ymin = conf.low, ymax = conf.high)
    env$ylab <- "Proportion expressed"
  } else {
    stop("Huh?")
  }
  
  with(env, {
    d |> 
      ggplot(aes(x = bin, y = y, ymin = ymin, ymax = ymax)) +
      th +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6),
        legend.position = "none"
      ) +
      geom_errorbar(width = 0, colour = "grey70") +
      geom_point(aes(colour = bin == "0")) +
      scale_colour_manual(values = c("black", "orange")) +
      facet_wrap(~gene_symbol, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, NA)) +
      labs(x = str_glue("{gene_x} expression bin"), y = ylab)
  })
    
}


#' Module enrichment heatmap
#'
#' Displays enrichment ratios for terms (rows) vs module IDs (columns) for one module/ontology.
#'
#' @param me Tibble with columns module, ontology, term_name, ids, and enrichment.
#' @param mdl Module label to select.
#' @param ont Ontology label to select.
#' @param max_name_length Maximum characters of term names to display.
#'
#' @return A heatmap grob.
plot_module_enrichment_heatmap <- function(me, mdl = "red", ont = "re", max_name_length = 50) {
  M <- me |>
    filter(module == mdl & ontology == ont) |> 
    mutate(name = str_sub(term_name, start = 1, end = max_name_length)) |>
    separate_longer_delim(ids, delim = ", ") |> 
    pivot_wider(id_cols = name, names_from = ids, values_from = enrichment) |>
    column_to_rownames("name") |>
    as.matrix()
  M[is.infinite(M)] <- 0
  M[is.na(M)] <- 0

  ggheatmap(M, with.y.text = TRUE, palette = "Purples", legend.name = "Enrichment ratio")
}

#' Any gene expressed per cell
#'
#' For a set of genes, returns whether any have counts > 0 in each cell.
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols.
#'
#' @return Tibble with cell_id and logical any_expressed.
get_any_expressed <- function(seu, genes) {
  get_gene_expression(seu, genes, layer = "counts") |> 
    group_by(cell_id) |> 
    summarise(any_expressed = any(expression > 0))
}

#' Sankey: MPC → Final cell type → EP status
#'
#' Builds a three-stage flow: MPC/non-MPC, final assigned cell type (by highest ES with p<limit), 
#' and EP/non-EP.
#'
#' @param seu Seurat object.
#' @param p_val_limit P-value threshold for ES column selection.
#' @param EP_markers Marker genes for EP state.
#' @param MPC_markers Marker genes for MPC state.
#'
#' @return A ggplot Sankey diagram.
plot_sankey_cell_type <- function(seu, p_val_limit = 0.05, EP_markers = c("NEUROD1", "NEUROG3", "PAX4", "PAX6"),
                                  MPC_markers = c("PTF1A", "SOX9", "NKX6-1", "PDX1")) {
  cell_col <- tibble(
    cell_type = c("Epsilon", "Delta", "Beta", "Alpha", "Ductal", "Acinar", "Blood", "Immune", "VSM", "Nervous"),
    colour = c(rep("lightblue", 4), rep("bisque", 2), rep("lightgoldenrod1", 3), "lightgreen")
  ) |> 
    mutate(
      es_name = str_glue("ES_{cell_type}"),
      pval_name = str_glue("p_val_{cell_type}")
    )

  EP <- get_any_expressed(seu, EP_markers) |> 
    mutate(EP = if_else(any_expressed, "EP", "non-EP")) |> 
    select(-any_expressed)
  MPC <- get_any_expressed(seu, MPC_markers) |> 
    mutate(MPC = if_else(any_expressed, "MPC", "non-MPC")) |> 
    select(-any_expressed)
  epmpc <- full_join(EP, MPC, by = "cell_id")

  d <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, all_of(c(cell_col$es_name, cell_col$pval_name))) |> 
    rename_with(~ str_replace(., "p_val", "pval"))  |> 
    pivot_longer(-cell_id, names_sep = "_", names_to = c("quantity", "cell_type")) |> 
    pivot_wider(id_cols = c(cell_id, cell_type), names_from = quantity, values_from = value) |> 
    filter(pval < p_val_limit) |> 
    group_by(cell_id) |> 
    summarise(Final = cell_type[which.max(ES)]) |> 
    inner_join(epmpc, by = "cell_id") |> 
    make_long(MPC, Final, EP) |> 
    mutate(node = fct_relevel(node, c("MPC", "non-MPC", "EP", "non-EP", cell_col$cell_type))) 

  env <- new.env(parent = globalenv())
  env$d <- d
  env$cell_col <- cell_col

  with(env, {
    ggplot(d, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
    theme_sankey() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank()
    ) +
    geom_sankey(flow.alpha = 0.6, flow.colour = "grey70", flow.fill = "grey80", node.alpha = 0.6) +
    geom_sankey_label() +
    scale_fill_manual(values = c(okabe_ito_palette[1:4], cell_col$colour),drop = FALSE) +
    labs(x = NULL)
  })
}

#' Heatmap of enrichment gradient across bins
#'
#' Orders terms by centre of mass across bins and shows log10(odds ratio) excluding bin 0.
#'
#' @param enr Tibble with columns term_name, bin, and odds_ratio.
#'
#' @return A heatmap grob.
plot_enrichment_gradient <- function(enr) {
  mx <- enr |> 
    filter(!is.infinite(odds_ratio)) |> 
    pull(odds_ratio) |> 
    max()
  d <- enr |> 
    filter(bin != 0) |>
    mutate(
      odds_ratio = if_else(is.infinite(odds_ratio), mx * 1.1, odds_ratio) |> log(10),
      ibin = as.integer(bin)
    ) |>
    group_by(term_name) |> 
    mutate(cm = sum(ibin * odds_ratio) /  sum(odds_ratio)) |> 
    ungroup() |> 
    mutate(term_name = fct_reorder(term_name, cm)) |> 
    pivot_wider(id_cols = term_name, names_from = bin, values_from = odds_ratio, values_fill = 0) |>
    arrange(term_name) |> 
    column_to_rownames("term_name")

  ggheatmap(d, order.col = FALSE, order.row = FALSE, with.y.text = TRUE, palette = "viridis", legend.name = expression(log[10]~OR))
}


plot_cluster_match_samples <- function(cluster_markers_1, cluster_markers_2, n_top = 50) {
  jaccard_cluster_cluster(cluster_markers_1, cluster_markers_2, n_top) |> 
    ggplot(aes(x = cluster_1, y = cluster_2, fill = j)) +
    th +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis") +
    labs(fill = "Jaccard")
}


plot_cluster_cell_type_jaccard <- function(seu) {
  env <- new.env(parent = globalenv())

  env$d <- seu |>
    jaccard_clusters_cell_types()
  with(env, {
    d |> 
      ggplot(aes(x = cluster, y = cell_type, fill = j)) +
      th +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis") +
      labs(fill = "Jaccard")
  })
}




feature_plot <- function(seu, gids, genes, assay = "RNA", layer = "data", min_cutoff = 0.1, n_col = NULL,
                         point_size = 0.4, max_n = 12, reorder = FALSE) {
  
  env <- new.env(parent = globalenv())
  env$n_col = n_col
  env$point_size = point_size
  
  # Extract umap embeddings
  um <- seu[["umap"]]@cell.embeddings |>
    as_tibble(rownames = "cell_id")
  
  # Extract normalised expression for the selection of genes
  expr <- GetAssayData(seu, assay = assay, layer = layer)
  all_genes <- rownames(expr)
  gids <- intersect(gids, all_genes)
  stopifnot(length(gids) > 0)
  expr <- expr[gids, , drop = FALSE] |> 
    t() |> 
    as.data.frame() |> 
    as_tibble(rownames = "cell_id")
  
  if(reorder) {
    # Order genes from highest to lowest expression
    gene_order <- expr |> 
      pivot_longer(-cell_id, names_to = "gene_symbol") |> 
      group_by(gene_symbol) |> 
      summarise(m = median(value)) |> 
      arrange(desc(m)) |> 
      pull(gene_symbol)
  } else {
    gene_order <-  gids
  }
  
  env$d <- um |> 
    inner_join(expr, by = "cell_id") |> 
    pivot_longer(-c(cell_id, umap_1, umap_2), names_to = "gene_symbol", values_to = "Expression") |> 
    filter(Expression > quantile(Expression, probs = min_cutoff)) |> 
    mutate(gene_symbol = factor(gene_symbol, levels = gene_order)) |> 
    filter(gene_symbol %in% head(gene_order, max_n))

  with(env, {
    d |> 
      ggplot(aes(x = umap_1, y = umap_2, colour = Expression)) +
      th +
      facet_wrap(~gene_symbol, ncol = n_col) +
      geom_point(size = point_size) +
      scale_colour_distiller(palette = "Purples", direction = 1) +
      labs(x = "umap 1", y = "umap 2")
  })
}

top_conserved_genes <- function(cons, max_n = 12) {
  cons |> 
    group_by(cluster, gene_symbol) |> 
    summarise(
      mfc = mean(avg_log2FC)
    ) |>
    top_n(max_n, mfc) |> 
    arrange(desc(mfc), .by_group = TRUE) |> 
    ungroup()
}




plot_top_cluster_features <- function(seu_clust, conserved_markers, genes, max_n = 12, n_col = 4,
                                      x_size = 11, y_size = 10, file_prefix = "top_clust_gene") {
  top_genes <- top_conserved_genes(conserved_markers, max_n)
  clusters <- top_genes$cluster |> 
    as.character() |> 
    unique() |> 
    sort()
  
  map(clusters, function(cl) {
    gname <- str_glue("{file_prefix}_{cl}")
    tg <- top_genes |> 
      filter(cluster == cl)
    gids <- tg$gene_symbol
    feature_plot(seu_clust, gids, genes, n_col = n_col) |> 
      gs(gname, width = x_size, height = y_size)
  })
}




plot_gene_segmentation <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf),
                                   seu_grid = NULL, with_centroids = FALSE) {
  poly <- get_polygons(seu)
  ex <- get_gene_expression(seu, genes)

  if(with_centroids) {
    centroids <- get_spatial_coordinates(seu, swap_xy = FALSE)
  }

  if(!is.null(seu_grid)) {
    xy <- get_spatial_coordinates(seu_grid, swap_xy = FALSE)
    corners <- get_resolution_grid(xy)
  }

  d <- inner_join(poly, ex, by = "cell_id", relationship = "many-to-many") 
  pl <- d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x, y = y, group = cell_id, fill = expression), colour = NA, linewidth = 0) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    facet_wrap(~gene_symbol) +
    coord_cartesian(xlim = c(box$x1, box$x2), ylim = c(box$y1, box$y2)) +
    labs(x = NULL, y = NULL)  
  if(!is.null(seu_grid)) {
    pl <- pl +
      geom_path(data = corners, aes(x = x_corner, y = y_corner, group = cell_id), colour = "black", linewidth = 0.3)
  }
  if(with_centroids) {
    pl <- pl +
      geom_point(data = centroids, aes(x = x, y = y), colour = "black", size = 0.5)
  }
  pl
}


plot_cell_area_distribution <- function(seu, bins = 160) {
  seu |> 
    get_polygons() |>
    group_by(cell_id) |>
    summarise(ar = polygon_area(x, y)) |>
  ggplot(aes(x = ar)) +
    th +
    geom_histogram(bins = bins) +
    geom_vline(xintercept = 64) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Cell area (µm²)", y = "Count")
}