
simple_cpm <- function(data, library_size = NULL, size_factors = NULL) {
  if (is.null(library_size)) {
    # Use diferent functions for sparse vs dense matrices
    if (inherits(data, "Matrix")) {
      library_size <- Matrix::colSums(data)
    } else {
      library_size <- colSums(data)
    }
  }
  
  if (length(library_size) != ncol(data)) {
    stop("The length of 'library_size' doesn't match the number of columns in 'data'.")
  }
  
  if (!is.null(size_factors)) {
    if (length(size_factors) != length(library_size)) {
      stop("'library_size' and 'size_factors' are not the same length.")
    }
    library_size <- library_size * size_factors
  }
  
  # Preserve sparsity for sparse matrices
  if (inherits(data, "sparseMatrix")) {
    data@x <- as.numeric(data@x / rep(library_size / 1e6, diff(data@p)))
    return(data)
  } else {
    return(sweep(data, 2, library_size, "/") * 1e6)
  }
}

simple_log2norm <- function(data, library_size = NULL, size_factors = NULL, pseudocount = 0.5) {
  if (pseudocount == 1 & inherits(data, "sparseMatrix")) {
    # Preserves sparsity
    sparse_data <- simple_cpm(data, library_size = library_size,
                              size_factors = size_factors)
    sparse_data@x <- log2(sparse_data@x + pseudocount)
    return(sparse_data)
  } else {
    return(log2(simple_cpm(data, library_size = library_size, size_factors = size_factors)
                + pseudocount))
  }
}



find_sex_mismatches <- function(metadata, data,
                                sample_colname = "specimenID",
                                sex_colname = "sex",
                                y_expr_threshold = 2.0) {
  y_genes <- c(ENSG00000129824 = "RPS4Y1",
               ENSG00000198692 = "EIF1AY",
               ENSG00000067048 = "DDX3Y",
               ENSG00000012817 = "KDM5D")
  sex_genes <- c(ENSG00000229807 = "XIST", y_genes)
  
  # Strip any version numbers off of the gene names if they are Ensembl IDs
  # before checking that all necessary genes exist
  rownames(data) <- stringr::str_replace(rownames(data), "\\.[0-9]+", "")
  
  if (!all(names(sex_genes) %in% rownames(data)) &&
      !all(sex_genes %in% rownames(data))) {
    stop(paste0("Data is missing the required genes. Ensure that these genes ",
                "are present: ", paste(sex_genes, collapse = ", "), " (",
                paste(names(sex_genes), collapse = ", "), ")."))
  }
  
  # If genes are Ensembl IDs, convert to gene symbols
  if (any(grepl("ENSG00", rownames(data)))) {
    data <- data[names(sex_genes), ]
    rownames(data) <- sex_genes
  }
  
  if (!(sample_colname %in% colnames(metadata))) {
    stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
  }
  
  if (!(sex_colname %in% colnames(metadata))) {
    stop(paste0("\"", sex_colname, "\" is not a valid column in `metadata`."))
  }
  
  if (!all(colnames(data) %in% metadata[, sample_colname])) {
    stop("`metadata` is missing samples that are present in `data`.")
  }
  
  sex_check <- metadata |>
    dplyr::select(dplyr::all_of(c(sample_colname, sex_colname))) |>
    merge(
      t(data[sex_genes, ]),
      by.x = sample_colname, by.y = "row.names"
    ) |>
    dplyr::mutate(
      mean_Y = rowMeans(dplyr::across(dplyr::all_of(y_genes))),
      reported_sex = .data[[sex_colname]],
      estimated_sex = dplyr::case_when(mean_Y > y_expr_threshold ~ "male",
                                       .default = "female"),
      sex_valid = (reported_sex == estimated_sex)
    ) |>
    as.data.frame()
  
  mismatches <- sex_check[, sample_colname][!sex_check$sex_valid]
  
  return(list(
    sex_check_df = sex_check,
    y_expr_threshold = y_expr_threshold,
    mismatches = mismatches
  ))
}

plot_sex_mismatch_results <- function(sex_check_df, y_expr_threshold = 2.0,
                                      print_plot = TRUE) {
  # Put FALSE last so they are plotted on top of other dots
  sex_check_df <- dplyr::arrange(sex_check_df, dplyr::desc(sex_valid))
  
  plt1 <- ggplot2::ggplot(sex_check_df,
                          ggplot2::aes(x = XIST, y = mean_Y, color = reported_sex)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  status_colors <- c("TRUE" = "black", "FALSE" = "red")
  
  plt2 <- ggplot2::ggplot(sex_check_df,
                          ggplot2::aes(x = XIST, y = mean_Y, color = sex_valid)) +
    ggplot2::geom_point(size = ifelse(sex_check_df$sex_valid, 0.5, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::geom_hline(yintercept = y_expr_threshold, linetype = "dotdash",
                        color = "red", alpha = 0.8)
  
  if (print_plot) {
    print(plt1)
    print(plt2)
  }
  return(list(plt1, plt2))
}


find_pca_outliers <- function(data,
                              n_sds = 4,
                              n_pcs = 3,
                              metadata = NULL,
                              sample_colname = "specimenID",
                              cutoff_method = "default",
                              gene_info = NULL) {
  
  # # If gene biotype information is provided, use protein-coding autosomal genes
  # # only. Otherwise, use all genes in the data.
  # if (!is.null(gene_info) & "gene_biotype" %in% colnames(gene_info) &
  #     "chromosome_name" %in% colnames(gene_info)) {
  #   genes_use <- subset(gene_info, gene_biotype == "protein_coding" &
  #                         !grepl("(X|Y|M|K|G)", chromosome_name)) |>
  #     dplyr::pull(ensembl_gene_id)
  # } else {
  #   genes_use <- rownames(data)
  # }
  # 
  # data <- data[intersect(genes_use, rownames(data)), ]
  # 
  # # Remove genes that are mostly 0's, which may be 0 or a negative number in
  # # log2-scale. For PCA, restrict to genes expressed in >= 80% of samples
  # genes_keep <- rowSums(data > min(data)) >= 0.8 * ncol(data)
  # data <- data[genes_keep, ]
  # 
  # # Remove genes with very low variance
  # variance <- matrixStats::rowVars(data)
  # genes_keep <- names(variance)[variance > 0.001]
  # data <- data[genes_keep, ]
  
  pca_res <- stats::prcomp(t(data), center = TRUE, scale. = TRUE)
  
  if (!is.null(metadata)) {
    if (!(sample_colname %in% colnames(metadata))) {
      stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
    }
    
    if (!all(colnames(data) %in% metadata[, sample_colname])) {
      stop("`metadata` is missing samples that are present in `data`.")
    }
    
    pca_df <- merge(metadata, pca_res$x,
                    by.x = sample_colname, by.y = "row.names")
  } else {
    sample_colname <- "sample"
    pca_df <- as.data.frame(pca_res$x)
    pca_df$sample <- rownames(pca_df)
  }
  
  pc_names <- paste0("PC", 1:n_pcs)
  
  thresholds <- sapply(pc_names, function(pc) {
    stats::sd(pca_df[, pc]) * n_sds
  })
  
  if (cutoff_method == "default") {
    # Check that each PC value is within its threshold, for each PC separately
    in_bounds <- lapply(pc_names, function(pc) {
      abs(pca_df[, pc]) < thresholds[pc]
    }) |>
      as.data.frame() |>
      rowSums() == n_pcs # All values should be true across the row
    
  } else if (cutoff_method == "ellipse") {
    # Check that each point is inside an ellipse (or ellipsoid, if n_pcs > 2),
    # where each radius is one of the thresholds.
    # The formula for an ellipse/ellipsoid is (x^2 / a^2) + (y^2 / b^2) + ... = 1,
    # where (x, y, ...) are PCs and (a, b, ...) are radii defined by n_sds * sd(PC).
    in_bounds <- lapply(pc_names, function(pc) {
      pca_df[, pc]^2 / thresholds[pc]^2
    }) |>
      as.data.frame() |>
      rowSums() <= 1
  }
  
  pca_df$is_outlier <- !in_bounds
  outliers <- pca_df[, sample_colname][pca_df$is_outlier]
  
  return(list(
    pca_df = pca_df,
    thresholds = thresholds,
    outliers = outliers
  ))
}

find_pca_outliers_by_group <- function(data, pca_group,
                                       n_sds = 4,
                                       metadata = NULL,
                                       sample_colname = "specimenID",
                                       gene_info = NULL,
                                       min_group_size = 10) {
  if (length(pca_group) == 1) {
    if (is.null(metadata)) {
      stop("`pca_group` is a single value but no metadata was supplied.")
    }
    if (!(pca_group %in% colnames(metadata))) {
      stop(paste0("\"", pca_group, "\" is not a valid column in `metadata`."))
    }
    pca_group <- metadata[, pca_group]
  } else if (length(pca_group) != ncol(data)) {
    stop("The number of samples in `pca_group` does not match the number of samples in `data`.")
  }
  
  if (!is.null(metadata)) {
    if (!(sample_colname %in% colnames(metadata))) {
      stop(paste0("\"", sample_colname, "\" is not a valid column in `metadata`."))
    }
    if (!all(colnames(data) %in% metadata[, sample_colname])) {
      stop("`metadata` is missing samples that are present in `data`.")
    }
    # Make sure data and metadata are in the same sample order
    data <- data[, metadata[, sample_colname]]
  }
  
  results <- lapply(unique(pca_group), function(grp) {
    if (sum(pca_group == grp) < min_group_size) {
      message(paste0("Group '", grp, "' is too small to run outlier detection. Skipping..."))
      return(NULL)
    }
    
    data_group <- data[, pca_group == grp]
    
    if (is.null(metadata)) {
      meta_group <- NULL
    } else {
      meta_group <- metadata[pca_group == grp, ]
    }
    
    find_pca_outliers(data_group,
                      n_sds = n_sds,
                      metadata = meta_group,
                      sample_colname = sample_colname,
                      gene_info = gene_info)
  })
  
  names(results) <- unique(pca_group)
  
  # Remove any NULL entries from groups that were too small
  results <- results[lengths(results) > 0]
  
  return(list(
    group_results = results,
    outliers = lapply(results, "[[", "outliers") |> unlist() |> as.vector()
  ))
}


#' Plot PCA Outlier Detection Results

plot_pca_outliers_ellipse <- function(pca_df, pc1_threshold, pc2_threshold,
                                      print_plot = TRUE,
                                      color = "is_outlier", ...) {
  # TODO add ggrepel labels and make outlier points larger. Find a way to get
  # title information
  # Formula for an ellipse is (x^2 / a^2) + (y^2 / b^2) = 1, so
  # y = +/- sqrt((1 - x^2 / a^2) * b^2)
  ellipse_points <- function(axis_A, axis_B) {
    x <- seq(from = -axis_A, to = axis_A, length.out = 1000)
    y <- sqrt((1 - x^2 / axis_A^2) * axis_B^2)
    
    data.frame(x = c(x, rev(x)), y = c(y, rev(-y)))
  }
  
  ellipse_df <- ellipse_points(pc1_threshold, pc2_threshold)
  
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (X %in% colnames(pca_df)) {
      return(rlang::sym(X))
    } else {
      message(paste0("\"", X, "\" is not a valid column in `pca_df`. ",
                     "It will not be used in the plot."))
    }
  })
  
  if (!is.null(color)) {
    if (color %in% colnames(pca_df)) {
      aes_opts$color = ifelse(is.character(color), rlang::sym(color), color)
    } else {
      message(paste0("\"", color, "\" is not a valid column in `pca_df`. ",
                     "No color aes will be used."))
    }
  }
  
  plt <- ggplot2::ggplot(pca_df,
                         ggplot2::aes(x = PC1, y = PC2, !!!aes_opts)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_path(data = ellipse_df,
                       ggplot2::aes(x = x, y = y),
                       linetype = "dotdash",
                       color = "blue",
                       inherit.aes = FALSE) +
    ggplot2::theme_bw()
  
  if (print_plot) {
    print(plt)
  }
  
  return(plt)
}

plot_pca_outliers <- function(pca_df, thresholds,
                              print_plot = TRUE,
                              color = "is_outlier", ...) {
  # TODO add ggrepel labels and make outlier points larger. Find a way to get
  # title information
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (X %in% colnames(pca_df)) {
      return(rlang::sym(X))
    } else {
      message(paste0("\"", X, "\" is not a valid column in `pca_df`. ",
                     "It will not be used in the plot."))
    }
  })
  
  if (!is.null(color)) {
    if (color %in% colnames(pca_df)) {
      aes_opts$color = ifelse(is.character(color), rlang::sym(color), color)
    } else {
      message(paste0("\"", color, "\" is not a valid column in `pca_df`. ",
                     "No color aes will be used."))
    }
  }
  
  plt <- ggplot2::ggplot(pca_df,
                         ggplot2::aes(x = PC1, y = PC2, !!!aes_opts)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -PC1, xmax = PC1,
                                    ymin = -PC2, ymax = PC2),
                       data = as.data.frame(t(thresholds)),
                       linetype = "dotdash",
                       color = "blue",
                       fill = NA,
                       inherit.aes = FALSE) +
    ggplot2::theme_bw()
  
  if (print_plot) {
    print(plt)
  }
  
  return(plt)
}

plot_pca_outliers_13 <- function(pca_df, thresholds,
                              print_plot = TRUE,
                              color = "is_outlier", ...) {
  # TODO add ggrepel labels and make outlier points larger. Find a way to get
  # title information
  # The lines below will auto-inject "..." into the aes() statement
  aes_opts <- lapply(list(...), function(X) {
    if (X %in% colnames(pca_df)) {
      return(rlang::sym(X))
    } else {
      message(paste0("\"", X, "\" is not a valid column in `pca_df`. ",
                     "It will not be used in the plot."))
    }
  })
  
  if (!is.null(color)) {
    if (color %in% colnames(pca_df)) {
      aes_opts$color = ifelse(is.character(color), rlang::sym(color), color)
    } else {
      message(paste0("\"", color, "\" is not a valid column in `pca_df`. ",
                     "No color aes will be used."))
    }
  }
  
  plt <- ggplot2::ggplot(pca_df,
                         ggplot2::aes(x = PC1, y = PC3, !!!aes_opts)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -PC1, xmax = PC1,
                                    ymin = -PC3, ymax = PC3),
                       data = as.data.frame(t(thresholds)),
                       linetype = "dotdash",
                       color = "blue",
                       fill = NA,
                       inherit.aes = FALSE) +
    ggplot2::theme_bw()
  
  if (print_plot) {
    print(plt)
  }
  
  return(plt)
}


# Bin ages in 5-year increments, except for categories "under 80" and "90+"
bin_ages <- function(ageDeath) {
  age_labels <- c("missing or unknown", "Under 80",
                  "80-84", "85-89", "90+")
  
  ageDeath <- case_when(ageDeath == "90+" ~ 95, # 'cut' will convert back to 90+
                        ageDeath == "missing or unknown" ~ -1, # 'cut' will convert back to missing
                        .default = suppressWarnings(as.numeric(ageDeath)))
  
  cut(as.numeric(ageDeath),
      breaks = c(-1, 0, 80, 85, 90, 100),
      labels = age_labels,
      right = FALSE) |>
    droplevels()
}

remove_unusable_covariates <- function(data, always_keep = c(), verbose = TRUE) {
  data <- as.data.frame(data)
  
  # Columns that have only one unique value or are all NA
  same_value <- data |>
    select(-any_of(always_keep)) |>
    summarize(across(everything(), ~ length(na.omit(unique(.x))) <= 1))
  
  # Columns that have a unique value for each sample
  all_unique <- data |>
    select(-any_of(always_keep), -where(is.numeric)) |>
    summarize(across(everything(), ~ length(unique(.x)) == nrow(data)))
  
  to_remove <- c(colnames(same_value)[same_value == TRUE],
                 colnames(all_unique)[all_unique == TRUE]) |>
    unique()
  
  if (verbose) {
    print(paste("Removing", length(to_remove), "columns:",
                paste(to_remove, collapse = ", ")))
  }
  
  data |> select(-any_of(to_remove))
}

remove_correlated_covariates <- function(data,
                                         id_cols = c(),
                                         always_keep = c(),
                                         R2_threshold = 0.5,
                                         verbose = TRUE) {
  data <- as.data.frame(data)
  
  # Number of NA entries in each column
  na_vars <- data |>
    summarize(across(everything(), ~sum(is.na(.x))))
  
  cols_analyze <- setdiff(colnames(data), id_cols)
  
  form <- paste("~", paste(cols_analyze, collapse = " + "))
  
  cor_mat <- variancePartition::canCorPairs(form, data, showWarnings = FALSE)
  r2_mat <- cor_mat^2
  
  to_remove <- .get_removals(r2_mat, na_vars, R2_threshold, always_keep)
  
  if (verbose) {
    print(paste("Removing", length(to_remove), "columns:",
                paste(to_remove, collapse = ", ")))
  }
  
  data |> select(-any_of(to_remove))
}

#' Get Which Variables to Remove Based on Correlation
#'
#' This is a helper function for [remove_correlated_covariates()] to decide
#' which variable in each highly-correlated pair of covariates should get
#' removed.
#'
#' @param r2_mat a matrix of R^2 values (or other positive association-like
#'   values), which must have row and column names. This matrix must be square.
#' @param na_vars a one-row data.frame where the columns are the covariates and
#'   the values are the number of `NA` values in each column.
#' @inheritParams remove_correlated_covariates
#'
#' @return a character vector with the names of the columns that should be
#'   removed. May also be an empty vector.
.get_removals <- function(r2_mat, na_vars, R2_threshold = 0.5, always_keep = c()) {
  if (!(nrow(r2_mat) == ncol(r2_mat)) ||
      !(all(rownames(r2_mat) == colnames(r2_mat))) ||
      !isSymmetric(r2_mat)) {
    stop("`r2_mat` is not square or symmetrical")
  }
  
  # Remove self-correlation
  diag(r2_mat) <- NA
  
  r2_melt <- r2_mat
  
  # Avoids picking up both (a vs b) and (b vs a)
  r2_melt[upper.tri(r2_melt, diag = TRUE)] <- 0
  
  r2_melt <- r2_melt |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(cols = -var1,
                        names_to = "var2", values_to = "value") |>
    subset(value >= R2_threshold) |>  # pairs with R^2 > 0.5 (cor ~ 0.7) only
    dplyr::arrange(dplyr::desc(value))
  
  if (nrow(r2_melt) == 0) {
    return(c())
  }
  
  removed <- c()
  
  for (R in 1:nrow(r2_melt)) {
    vars <- as.character(c(r2_melt$var1[R], r2_melt$var2[R]))
    
    if (any(vars %in% removed)) {
      # No need to re-check if we're already removing one of these variables
      next
    } else if (any(vars %in% always_keep)) {
      # If either variable is in the always_keep list, only remove the variable
      # that isn't in always_keep. If both are in the list, neither one gets
      # removed.
      removed <- c(removed, setdiff(vars, always_keep))
    } else if (any(na_vars[1, vars] > 0) && (na_vars[1, vars[1]] != na_vars[1, vars[2]])) {
      # If either variable has any NA values, and they don't have the same number
      # of NAs, remove the one with the most NAs
      removed <- c(removed,
                   vars[which.max(na_vars[1, vars])])
    } else {
      # No NAs, no variables are in `always_keep` or `removed`
      cur_vars <- setdiff(rownames(r2_mat), removed)
      mean_r2 <- rowMeans(r2_mat[cur_vars, cur_vars], na.rm = TRUE)
      
      # Remove the variable with the largest mean R^2 with the other remaining variables
      removed <- c(removed,
                   vars[which.max(mean_r2[vars])])
    }
  }
  
  return(unique(removed))
}

plot_resid_qq <- function(resid_mat) {
  set.seed(100)
  rand_genes_100 <- sample(rownames(resid_mat), 100)
  rand_genes_25 <- sample(rand_genes_100, 25)
  
  resid_df <- as.data.frame(t(resid_mat[rand_genes_100, ])) |>
    tidyr::pivot_longer(everything(), names_to = "gene", values_to = "expr")
  
  plt1 <- ggplot(subset(resid_df, gene %in% rand_genes_25), aes(sample = expr)) +
    stat_qq(size = 0.1) + stat_qq_line(color = "red") +
    theme_bw() +
    facet_wrap(~gene, nrow = 5) +
    ggtitle("QQ plot for 25 random genes")
  
  print(plt1)
  
  plt2 <- ggplot(resid_df, aes(sample = expr)) +
    stat_qq(size = 0.1) + stat_qq_line(color = "red") +
    theme_bw() +
    ggtitle("QQ plot for 100 random genes together")
  
  print(plt2)
}
