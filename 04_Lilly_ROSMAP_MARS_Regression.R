library(dplyr)
library(stringr)
library(matrixStats)
library(variancePartition)
#library(reformulas)
library(synapser)

load("QC_Analysis_Functions.R")

synLogin()

project_pid<-"syn72386347"

cqn_data <- readRDS("Lilly_ROSMAP_MARS_data_cqn.rds")
clean_covariates <- readRDS("RUSH_cleaned_covariates_for_regression.rds")

biovars<-c("ADoutcome", "ageDeath", "apoeGenotype", "RaceEthnicity", "sex")

# Turn sex into binary coding to include as fixed effect
clean_covariates$sex<-as.numeric(factor(clean_covariates$sex))-1

#Model recommendations from `variancePartition` and `mvIC`
models_vp <- 
  readRDS("Lilly_ROSMAP_MARS_data_regression_var_part_results.rds")


models_mvic <- 
  readRDS("Lilly_ROSMAP_MARS_data_regression__mvIC_results.rds")


# variancePartition recommendations
variables <- models_vp$recommended_vars

mixed_effects <- clean_covariates |>
  select(where(is.character), where(is.factor)) |>
#  select(!all_of(biovars)) |>
  colnames() |>
  intersect(variables)

fixed_effects <- setdiff(variables, mixed_effects)

if (length(mixed_effects) > 0) {
  mixed_effects <- paste0("(1 | ", mixed_effects, ")")
}

formul<-paste("~", paste(c(fixed_effects, mixed_effects), collapse = " + "))


#`mvIC` recommendations
formsul.mvIC <- 
  paste(as.character(models_mvic$formula), collapse = " ")

print(formsul.mvIC)



# Percent variance explained by covariates

  vp <- models_vp$var_part$var_part_all
  
  plt <- plotVarPart(sortCols(vp)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "(Sorted by median)")
  print(plt)
  
  plt <- plotVarPart(sortCols(vp, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "(Sorted by mean)")
  print(plt)

  models_vp$vp_summary_all |>
    mutate(mean_pct = round(mean_pct, digits = 5),
           median_pct = round(median_pct, digits = 5)) |>
#    mutate(tissue = tissue, .before = covariate) |>
    select(-total_rank) |>
    as.data.frame() |>
    print()

  # Batch / race overlaps
  plot_batch_overlaps <- function(batch1, batch2) {
    overlap <- table(batch1, batch2)
    
    # Force 0 to be white
    colors <- c("#FFFFFF", corrplot::COL1("Reds", n = 199))
    corrplot::corrplot(overlap, method = "color", is.corr = FALSE,
                       col = colors)#, addCoef.col = "gray")
  }
    covariates <- clean_covariates
    plot_batch_overlaps(covariates$RaceEthnicity, covariates$sequencingBatch)

    
    
    # Regression
    

      # In case there is a random element anywhere in this process, set a seed for 
      # reproducibility
      set.seed(250000)
      
      data_tissue <- cqn_data
      data_tissue <- data_tissue$y + data_tissue$offset
      
      covariates <- clean_covariates
      data_tissue <- data_tissue[, covariates$specimenID]
      
      rownames(covariates) <- covariates$specimenID
      
      stopifnot(all(colnames(data_tissue) == covariates$specimenID))
      
      
      fit <- fitVarPartModel(data_tissue, formul, covariates) #,
#                             BPPARAM = MulticoreParam(n_cores))
      
      resid_tissue <- residuals(fit, data_tissue)
      rownames(resid_tissue) <- rownames(data_tissue)
      
      # Since this is lmm, use fixef to extract intercept
      intercep_tissue<-lapply(fit, function(x){
        lme4::fixef(x)[1]
      })
      intercep_tissue<-unlist(intercep_tissue)
      adj_tissue<-resid_tissue+intercep_tissue
      
      plot(apply(data_tissue, 1, mean), apply(adj_tissue, 1, mean),
           xlab="Mean Expression (Unadjusted)", ylab="Mean Expression (Adjusted)", 
           main ="~ RIN + sex + (1|sequencingBatch)")
      cor(apply(adj_tissue, 1, mean), apply(data_tissue, 1, mean))

      plot_resid_qq(resid_tissue)
      
      write.csv(resid_tissue,"Lilly_ROSMAP_MARS_regression_varPar_model_residuals.csv",
                quote = FALSE)
      
      write.csv(adj_tissue,"Lilly_ROSMAP_MARS_regression_varPar_model_adjusted.csv",
                quote = FALSE)
      
      mdfile<-File("Lilly_ROSMAP_MARS_regression_varPar_model_adjusted.csv", parent=project_pid, description="Expression adjusted for RIN, Sex, Sequencing Batch")
      mdfile<-synStore(mdfile)
      
      # Write metadata values
      model_vars <- all.vars(as.formula(formul))
      write.csv(covariates[, c("specimenID", model_vars)],
                "Lilly_ROSMAP_MARS_regression_varPar_model_covariates.csv",
                row.names = FALSE, quote = FALSE)
      
    
      fit.mvIC <- fitVarPartModel(data_tissue, formsul.mvIC, covariates) #,
      #                             BPPARAM = MulticoreParam(n_cores))
      
      resid_tissue.mvIC <- residuals(fit.mvIC)
      rownames(resid_tissue.mvIC) <- rownames(data_tissue)
      
      # Since this is lm vs lmm, use different method to extract intercept
      intercep_tissue.mvIC<-lapply(fit.mvIC, function(x){
        x$coefficients[1]
      })
      intercep_tissue.mvIC<-unlist(intercep_tissue.mvIC)
      adj_tissue.mvIC<-resid_tissue.mvIC+intercep_tissue.mvIC
      
      plot(apply(adj_tissue.mvIC, 1, mean), apply(data_tissue, 1, mean))
      plot(apply(data_tissue, 1, mean), apply(adj_tissue.mvIC, 1, mean),
           xlab="Mean Expression (Unadjusted)", ylab="Mean Expression (Adjusted)", 
           main ="~ DV200")
      cor(apply(adj_tissue.mvIC, 1, mean), apply(data_tissue, 1, mean))

      plot_resid_qq(resid_tissue.mvIC)
      
      write.csv(resid_tissue.mvIC,"Lilly_ROSMAP_MARS_regression_mvIC_model_residuals.csv",
                quote = FALSE)
      
      write.csv(adj_tissue.mvIC,"Lilly_ROSMAP_MARS_regression_mvIC_model_adjusted.csv",
                quote = FALSE)
      
      mdfile<-File("Lilly_ROSMAP_MARS_regression_mvIC_model_adjusted.csv", parent=project_pid, description="Expression adjusted for DV200")
      mdfile<-synStore(mdfile)
      
      # Write metadata values
      model_vars.mvIC <- all.vars(as.formula(formsul.mvIC))
      write.csv(covariates[, c("specimenID", model_vars.mvIC)],
                "Lilly_ROSMAP_MARS_regression_mvIC_model_covariates.csv",
                row.names = FALSE, quote = FALSE)
      
      
      
      
      
      #Some plots
      
      library(ggplot2)
      library(viridis)
      
      
      mvIC_plt_data <-   adj_tissue.mvIC |> # Un-adjusted log2-cpm data
        as.data.frame() |>
        tibble::rownames_to_column("ensembl_gene_id") |>
        tidyr::pivot_longer(cols = -ensembl_gene_id,
                            names_to = "specimenID",
                            values_to = "raw_expr")
      
      plt1 <- ggplot(mvIC_plt_data, aes(x = raw_expr, color = specimenID)) +
        geom_density() +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title = "Adjusted data (mvIC Model)", x = "Adjusted expression", y = "Density") +
        scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)
      
      varpar_plt_data <-   adj_tissue |> # Un-adjusted log2-cpm data
        as.data.frame() |>
        tibble::rownames_to_column("ensembl_gene_id") |>
        tidyr::pivot_longer(cols = -ensembl_gene_id,
                            names_to = "specimenID",
                            values_to = "raw_expr")
      
      plt2 <- ggplot(varpar_plt_data, aes(x = raw_expr, color = specimenID)) +
        geom_density() +
        theme_bw() +
        theme(legend.position = "none") +
        labs(title = "Adjusted data (variancePartition Model)", x = "Adjusted expression", y = "Density") +
        scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)
          