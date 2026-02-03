library(edgeR)
library(cqn)
library(stringr)
library(ggplot2)
library(viridis)
library(patchwork)
library(synapser)

load("QC_Analysis_Functions.R")

project_pid<-"syn72386347"

gene_info <- read.csv("GC_content_GRCh38_v49.csv")

data <- readRDS("Lilly_ROSMAP_MARS_data_qced.rds")
data <- DGEList(round(data$counts), samples = data$metadata)

data$samples$group <- data$samples$ADoutcome

# Filter to genes that are expressed in at least one group
genes_keep <- edgeR::filterByExpr(data, group = data$samples$group)
data_sub<-data[genes_keep, ]

rownames(gene_info) <- stringr::str_replace(gene_info$ensembl_gene_id, "\\.[0-9]+", "")


# ---- run-cqn ----


  # There seems to be a random element to cqn() so we set a seed for reproducibility
set.seed(500000)

gene_info_sub <- gene_info[rownames(data_sub$counts), ]
  
all(rownames(data_sub$counts) == rownames(gene_info_sub$ensembl_gene_id))

missinggenes<-rownames(data_sub$counts)[which(is.na(gene_info_sub$ensembl_gene_id))]

#
# Remove 2 genes missing from gene_info
#

genes_keep2<- !rownames(data_sub$counts) %in% missinggenes

data_sub<-data_sub[genes_keep2, ]
gene_info_sub <- gene_info[rownames(data_sub$counts), ]

  
cqn_data <- cqn(data_sub$counts,
                x = gene_info_sub$percent_gc_content,
                lengths = gene_info_sub$gene_length,
                lengthMethod = "smooth")

cqn_data<-cqn_data[c("y", "offset", "glm.offset")]


saveRDS(cqn_data, "Lilly_ROSMAP_MARS_data_cqn.rds")

# Write CQN data to a CSV file and optionally upload to Synapse

cqn_file <- "Lilly_ROSMAP_MARS_data_cqn.csv"
write.csv(cqn_data$y + cqn_data$offset,
            cqn_file,
            quote = FALSE)


syn_file <- File(cqn_file, parent = project_pid)
syn_file <- synStore(
  syn_file,
  forceVersion = FALSE
)
  



cqn_obj<-cqn_data
cqn_plt_data <-   (cqn_obj$y + cqn_obj$offset) |>
    as.data.frame() |>
    tibble::rownames_to_column("ensembl_gene_id") |>
    tidyr::pivot_longer(cols = -ensembl_gene_id,
                        names_to = "specimenID",
                        values_to = "cqn_expr")



raw_plt_data <-   cqn_obj$y |> # Un-adjusted log2-cpm data
    as.data.frame() |>
    tibble::rownames_to_column("ensembl_gene_id") |>
    tidyr::pivot_longer(cols = -ensembl_gene_id,
                        names_to = "specimenID",
                        values_to = "raw_expr")



  cqn_plt_tissue <- cqn_plt_data
  raw_plt_tissue <- raw_plt_data
  
  plt1 <- ggplot(raw_plt_tissue, aes(x = raw_expr, color = specimenID)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Raw log2-CPM data", x = "Raw log2 expression", y = "Density") +
    scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)
  
  plt2 <- ggplot(cqn_plt_tissue, aes(x = cqn_expr, color = specimenID)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "CQN", x = "CQN expression", y = "Density") +
    scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)
  
  print(plt1 + plt2 + plot_annotation(title = "ROSMAP/MARS DLPFC"))


