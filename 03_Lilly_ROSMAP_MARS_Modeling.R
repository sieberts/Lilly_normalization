library(dplyr)
library(stringr)
library(matrixStats)
library(mvIC)
library(variancePartition)

load("QC_Analysis_Functions.R")

raw_data <- readRDS("Lilly_ROSMAP_MARS_data_qced.rds")
cqn_data <- readRDS("Lilly_ROSMAP_MARS_data_cqn.rds")

bio_vars <- c("ADoutcome", "ageDeath", "apoeGenotype", "RaceEthnicity", "sex")

# Clean up and filter technical variables --------------------------------------

numerical_vars <- raw_data$metadata |>
  dplyr::select(specimenID, RIN, DV200, PMI, RIN.imputed, DV200.imputed) |>
  
  # Ensure all data is numeric
  mutate(across(c(-specimenID), as.numeric))

## Set sample with missing RIN and DV200 back to NA for reimputation-Now fixed in upstream analysis
#numerical_vars$RIN.imputed[is.na(numerical_vars$RIN)]<-NA
#numerical_vars$DV200.imputed[is.na(numerical_vars$RIN)]<-NA

numerical_vars$DV200<-numerical_vars$DV200.imputed
numerical_vars <- numerical_vars |>
  dplyr::select(specimenID, RIN, DV200, PMI)


categorical_vars <- raw_data$metadata |>
  dplyr::select(
    specimenID, all_of(bio_vars),
    # Technical variables
    rnaBatch, libraryBatch, sequencingBatch, cohort
  ) |>
  mutate(ageDeath = bin_ages(ageDeath))

all_covariates <- merge(numerical_vars, categorical_vars)


  # Impute NA values and scale numeric covariates.
  clean_covariates <- all_covariates |>
    
    # DV200 and RIN can be imputed as the median of their rna batch. This is
    # done before removing unusable covariates in case RNA batch gets removed
    # (which happens with Columbia data).
    group_by(rnaBatch) |>
    mutate(across(any_of(c("DV200", "RIN")),
                  ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) |>
    ungroup() |>
    
    # Remove covariates that have all NAs or only 1 unique value in the column
    remove_unusable_covariates(always_keep = c("specimenID", "tissue")) |>
    
    mutate(
      # Other numeric variables with NAs (e.g. PMI) are set to the median of the whole column
      across(c(where(is.numeric), -any_of(c("DV200", "RIN"))),
             ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)),
      
      # Scale
      across(where(is.numeric), ~ as.numeric(scale(.x))),
      
      # Change NA values for factors to "missing or unknown"
      across(c(where(is.character), where(is.factor), -specimenID),
             ~ factor(ifelse(is.na(.x), "missing or unknown", as.character(.x))))
    ) |>
    
    # Convert from tibble -> data.frame
    as.data.frame()
  

  # Rush only:
  # - There are mostly 2 libraryBatches to 1 rnaBatch, with very few
  #   libraryBatches that have more than one rnaBatch. sequencingBatch is
  #   almost exactly 1:many with rnaBatch. 
  #   There is an extremely uneven distribution of race across
  #   some rna, library and sequencing Batches in the DLPFC, so I've decided to use
  #   sequencingBatch for that tissue instead even though there are only 3 batches.
  # - Remove cohort as a variable: the MARS cohort has only samples that are
  #   Black / African American, and some cohorts only have 1-2 samples total.
    clean_covariates <- clean_covariates |>
      select(-cohort, -rnaBatch, -libraryBatch)


saveRDS(clean_covariates,
        "RUSH_cleaned_covariates_for_regression.rds")


#var_part_list <- lapply(clean_covariates, function(covariates) {
  tissue <- "DLPFC"
  
  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(75000)
  
  data_tissue <- cqn_data$y + cqn_data$offset
  data_tissue <- data_tissue[, clean_covariates$specimenID]
  
  # Do this on the top 5000 variable genes only so it doesn't take too long to run
  var_genes <- rowVars(data_tissue) |> sort(decreasing = TRUE) |> names()
  
  data_tissue <- data_tissue[var_genes[1:5000], ]
  
  # Remove specimenID and tissue from the possible covariates
  covariates <- clean_covariates |>
    tibble::column_to_rownames("specimenID") 
  
  # Just to see a printout of what this function would remove. Doesn't actually
  # remove anything.
  tmp <- covariates |> remove_correlated_covariates()
  # Would remove DV200
  
  covariates_tech <- covariates |>
    select(-any_of(bio_vars))
  
  # All categorical values have to be mixed variables for fitExtractVarPartModel
  mixed_vars <- covariates |>
    select(where(is.character), where(is.factor)) |>
    colnames()
  
  mixed_vars_tech <- setdiff(mixed_vars, bio_vars)
  
  fixed_vars <- setdiff(colnames(covariates), mixed_vars)
  
  if (length(mixed_vars) > 0) {
    mixed_vars <- paste0("(1 | ", mixed_vars, ")")
  }
  
  if (length(mixed_vars_tech) > 0) {
    mixed_vars_tech <- paste0("(1 | ", mixed_vars_tech, ")")
  }
  
  mixed_vars_hybrid <- mixed_vars[!grepl("AD|apoe|Race", mixed_vars)]
  
  form_all <- paste("~", paste(c(fixed_vars, mixed_vars), collapse = " + "))
  form_tech <- paste("~", paste(c(fixed_vars, mixed_vars_tech), collapse = " + "))
  form_hybrid <- paste("~", paste(c(fixed_vars, mixed_vars_hybrid), collapse = " + "))
  
  var_part1 <- fitExtractVarPartModel(data_tissue, form_all, covariates) |> #,
   #                                   BPPARAM = MulticoreParam(n_cores)) |>
    sortCols()
  
  var_part2 <- fitExtractVarPartModel(data_tissue, form_tech, covariates) |> #,
  #                                    BPPARAM = MulticoreParam(n_cores)) |>
    sortCols()

  var_part3 <- fitExtractVarPartModel(data_tissue, form_hybrid, covariates) |> #,
    #                                    BPPARAM = MulticoreParam(n_cores)) |>
    sortCols()  
  
  
  plt1 <- plotVarPart(var_part1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "All variables, sorted by median")
  print(plt1)
  
  plt2 <- plotVarPart(sortCols(var_part1, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "All variables, sorted by mean")
  print(plt2)
  
  plt3 <- plotVarPart(var_part2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "Technical variables, sorted by median")
  print(plt3)
  
  plt4 <- plotVarPart(sortCols(var_part2, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "Technical variables, sorted by mean")
  print(plt4)
  
  plt5 <- plotVarPart(var_part3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "Variables, sorted by median")
  print(plt5)
  
  plt6 <- plotVarPart(sortCols(var_part3, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Sources of variance",
         subtitle = "Variables, sorted by mean")
  print(plt6)
  
  vp_summary <- function(var_part) {
    var_part |>
      as.data.frame() |>
      tibble::rownames_to_column("gene") |>
      select(-Residuals) |>
      tidyr::pivot_longer(-gene, names_to = "covariate", values_to = "value") |>
      group_by(covariate) |>
      summarize(mean_pct = mean(value),
                median_pct = median(value)) |>
      mutate(rank_mean = rank(-mean_pct),
             rank_median = rank(-median_pct),
             total_rank = (rank_mean + rank_median) / 2) |>
      arrange(desc(mean_pct))
  }
  
  

    vp <- vp_summary(var_part1)
    print(vp)
    
    covariates <- clean_covariates |>
      tibble::column_to_rownames("specimenID") 

    
    form_cc <- paste("~", paste0(colnames(covariates), collapse = " + "))
    cc <- canCorPairs(form_cc, covariates, showWarnings = FALSE)
    cc[upper.tri(cc, diag = TRUE)] <- NA
    
    print(cc)
    
#                       RIN      DV200        PMI ADoutcome  ageDeath apoeGenotype       sex sequencingBatch
#RaceEthnicity   0.42731255 0.42617669 0.17181550 0.1457006 0.2177274    0.2997285 0.2112233       0.5527521
    
    cor_pairs <- cc |>
      as.data.frame() |>
      tibble::rownames_to_column("var1") |>
      tidyr::pivot_longer(-var1, names_to = "var2", values_to = "cor_value") |>
      subset(cor_value > 0.7) |>
      arrange(desc(cor_value)) |>
      rowwise() |>
      mutate(
        var1_rank = vp$rank_mean[vp$covariate == var1],
        var2_rank = vp$rank_mean[vp$covariate == var2],
        lowest_pct = which.max(c(var1_rank, var2_rank)),
        var_remove = c(var1, var2)[lowest_pct]
      ) |>
      ungroup()
    
    print(cor_pairs)
    # Remove DV200
    
    to_remove <- c()
    for (N in 1:nrow(cor_pairs)) {
      pair <- c(cor_pairs$var1[N], cor_pairs$var2[N])
      
      # Only remove one of the two variables if neither has been removed yet
      if (!any(pair %in% to_remove)) {
        to_remove <- c(to_remove, cor_pairs$var_remove[N])
      }
    }
    
    print(paste("Removing", length(to_remove), "variables:",
                paste(to_remove, collapse = ", ")))
    
    
    keep_effects<-c("ADoutcome", "apoeGenotype", "RaceEthnicity")
    
    vars_keep <- vp |>
      arrange(desc(mean_pct)) |>
      subset(mean_pct > 0.01 & !(covariate %in% c(to_remove, keep_effects))) |>
      pull(covariate)
    
    print(paste("Recommended variables:",
                paste(vars_keep, collapse = ", ")))
#    "Recommended variables: RIN, sex, sequencingBatch"
    
    saveRDS(list(recommended_vars = vars_keep,
                 var_part = list(var_part_all = var_part1,
                                 var_part_tech = var_part2,
                                 var_part_model = var_part3),
                 vp_summary_all = vp,
                 vp_summary_tech = vp_summary(var_part2),
                 vp_summary_model = vp_summary(var_part3),
                 correlations = cc,
                 cor_pairs = cor_pairs),
            file.path("Lilly_ROSMAP_MARS_data_regression_var_part_results.rds"))
  
  
  
# Find a potential model for each tissue using mvIC ----------------------------


  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(16)
  
  data_sub <- cqn_data
  data_sub <- data_sub$y + data_sub$offset
  
  meta_sub <- clean_covariates
  data_sub <- data_sub[, meta_sub$specimenID]
  
  stopifnot(all(colnames(data_sub) == meta_sub$specimenID))
  
  variables <- setdiff(colnames(meta_sub), c("specimenID", "tissue"))
  baseFormula <- "~1"
  
  mixed_vars <- intersect(variables,
                          c("sequencingBatch",  "rnaBatch", "libraryBatch", "cohort"))
  
  fixed_vars <- setdiff(variables, mixed_vars)
  
  if (length(mixed_vars) > 0) {
    mixed_vars <- paste0("(1 | ", mixed_vars, ")")
  }
  
  results <- mvIC::mvForwardStepwise(data_sub,
                                     baseFormula = baseFormula,
                                     data = meta_sub,
                                     variables = c(fixed_vars, mixed_vars))
  print(results$formula)
# ~1 + DV200  
  saveRDS(results,
          file = "Lilly_ROSMAP_MARS_data_regression__mvIC_results.rds")




