library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(furrr)
library(piano)
library(monocle3)

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)
options(stringsAsFactors = FALSE)

### Source files and functions including those from our previous sci-plex repository 
### (https://github.com/cole-trapnell-lab/sci-plex)
source("sci-plex/bin/cell_cycle.R")
source("sci-plex/bin/loadGSCSafe.R")
source("sci-plex/bin/GSA_helper_functions.R")
cc.genes = readRDS("sci-plex/bin/cc.genes.RDS")

#################################################################################
calculate_aggregate_expression_score <- function(cds, signature_genes, from_id = FALSE){

  if(from_id == TRUE){
    cds_subset = cds[rowData(cds)$id %in% signature_genes,]
  }
  else(cds_subset = cds[rowData(cds)$gene_short_name %in% signature_genes,])
  aggregate_signature_expression = exprs(cds_subset)
  #aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  #signature_score = log(aggregate_signature_expression+1)
  return(aggregate_signature_expression)
}

test_signature_kd_dependence <- function(signature_cds_subset, 
  signature_cds, 
  reference_cells=NULL, 
  min_fraction_cells=0.05, 
  pseudodose=0.001, 
  residualModelFormulaTerms=NULL){

  cell_ids_to_test = as.character(signature_cds_subset$cell)
  print(length(cell_ids_to_test))

  signature_cds_subset = signature_cds[,union(cell_ids_to_test, reference_cells)]

  target = unique(pData(signature_cds_subset)$gene_id)
  target = target[target != "NTC"]

  colData(signature_cds_subset)$gene_id <- factor(colData(signature_cds_subset)$gene_id, 
                                                  levels = c("NTC",target))

  message(paste("Testing cells perturbed for", target))
  message(paste("Fitting knockdown effect models to", nrow(signature_cds_subset), "genes"))

  #colData(signature_cds_subset)$pseudodose <- log(colData(signature_cds_subset)$dose + pseudodose)

  modelFormulaStr = "~gene_id + replicate + log(dose + 0.01)"

  if (is.null(residualModelFormulaTerms) == FALSE & length(residualModelFormulaTerms) > 0){
    for (i in (1:length(residualModelFormulaTerms))){
      if (length(unique(pData(signature_cds_subset)[,residualModelFormulaTerms[i]])) > 1){
        modelFormulaStr = paste(modelFormulaStr, "+", residualModelFormulaTerms[i])
      }
    }
  }
  
  model_tbl = fit_models(signature_cds_subset, 
                         model_formula_str = modelFormulaStr, 
                         reltol=1e-5, 
                         verbose=TRUE)
  model_coefs = coefficient_table(model_tbl) %>% dplyr::select(-model,-model_summary)
  model_coefs

}

#######################################################################
cds.list <- readRDS("cds.list.RDS")

MT_genes <- read.table("MT_genes.txt", header = FALSE)
MT_genes <- as.character(MT_genes$V1)

for(cell_type in names(cds.list)){

  print(dim(cds.list[[cell_type]]))

  cds.list[[cell_type]] <- cds.list[[cell_type]][!(rowData(cds.list[[cell_type]])$gene_short_name %in% MT_genes),
                                                 !is.na(colData(cds.list[[cell_type]])$gene_id) &
                                                 colData(cds.list[[cell_type]])$guide_number == 1]

  cds.list[[cell_type]] <- estimate_size_factors(cds.list[[cell_type]])
  cds.list[[cell_type]] <-detect_genes(cds.list[[cell_type]])

  print(dim(cds.list[[cell_type]]))

}

NTC_cds.list <- list()

for(cell_type in names(cds.list)){

  NTC_cds.list[[cell_type]] <- cds.list[[cell_type]][,colData(cds.list[[cell_type]])$gene_id %in% c("NTC","random") &
                                                      colData(cds.list[[cell_type]])$gRNA_maxCount >= 5]
  
  NTC_cds.list[[cell_type]] <- estimate_size_factors(NTC_cds.list[[cell_type]])

  print(dim(NTC_cds.list[[cell_type]]))

}

expressed_genes.list <- list()

for(cell_type in names(NTC_cds.list)){

genes_expressed_per_drug = 
  colData(NTC_cds.list[[cell_type]]) %>% 
  as.data.frame() %>%
  group_by(treatment) %>%
  nest() %>%
  mutate(fraction_genes = purrr::map(data, .f=function(pdata_subset, cds) {
    cds_subset = cds[,as.character(pdata_subset$cell)]
    cds_subset = detect_genes(cds_subset)
    tibble(id = rowData(cds_subset)$id,
           gene_short_name = rowData(cds_subset)$gene_short_name,
           num_cells_expressed = rowData(cds_subset)$num_cells_expressed,
           fraction_cells_expressed = rowData(cds_subset)$num_cells_expressed / ncol(cds_subset))
  }, NTC_cds.list[[cell_type]]))

genes_expressed_per_drug = 
  genes_expressed_per_drug %>% 
  unnest(fraction_genes) %>%
  dplyr::select(everything(),-data)

expressed_genes.list[[cell_type]] = 
  genes_expressed_per_drug %>% 
  filter(fraction_cells_expressed > 0.05) %>% 
  ungroup() %>% 
  pull(id) %>%
  unique()

  print(length(expressed_genes.list[[cell_type]]))

}

treatment_dose_dependent_diff_test.list <- list()

for(cell_type in names(NTC_cds.list)){

  treatment_dose_dependent_diff_test.list[[cell_type]] <- list()

  for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

    message(paste0("Determining drug/dose responsive gene expression in ",drug , " exposed ", cell_type, " cells"))

    cds_subset <- NTC_cds.list[[cell_type]][,colData(NTC_cds.list[[cell_type]])$treatment %in% c("vehicle",drug)]

    diff_test <- fit_models(cds_subset[expressed_genes.list[[cell_type]],], 
                            model_formula_str = "~log(dose + 0.01) + replicate",
                            cores =  15)

    treatment_dose_dependent_diff_test.list[[cell_type]][[drug]] <- coefficient_table(diff_test) %>% dplyr::select(-model,-model_summary)
    treatment_dose_dependent_diff_test.list[[cell_type]][[drug]]$cell_type <- rep(cell_type,dim(treatment_dose_dependent_diff_test.list[[cell_type]][[drug]])[1])
    treatment_dose_dependent_diff_test.list[[cell_type]][[drug]]$treatment <- rep(drug,dim(treatment_dose_dependent_diff_test.list[[cell_type]][[drug]])[1])

  }

  treatment_dose_dependent_diff_test.list[[cell_type]] <- do.call("rbind",treatment_dose_dependent_diff_test.list[[cell_type]])
  rm(cds_subset,diff_test)
}

# saveRDS(treatment_dose_dependent_diff_test.list, "treatment_dose_dependent_diff_test_updated_113020.list.rds")
treatment_dose_dependent_diff_test.list <- readRDS("treatment_dose_dependent_diff_test_updated_113020.list.rds")

treatment_dose_dependent_diff_test.list[["A172"]] %>%  
  filter(grepl("dose", term)) %>%
  ggplot() +
  geom_point(aes(x = normalized_effect,
                 y = -log10(q_value),
                 color = q_value < 0.01),
             size = 0.25, stroke = 0) +
  facet_wrap(~treatment, ncol = 2) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.01), color = "dimgrey", linetype = "dotted",size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("β coefficient") +
  ylab("-log(q value)") +
  scale_color_manual("Pass\nFilter\n(FDR < 1%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("QC_plots/Effect_size_distribution_NTC_all_treatment_dose_response_diff_test_A172_Supplementary_Figure_6A.png", 
         dpi = 600, width = 2, height = 1.75)

treatment_dose_dependent_diff_test.list[["T98G"]] %>%  
  filter(grepl("dose", term)) %>%
  ggplot() +
  geom_point(aes(x = normalized_effect,
                 y = -log10(q_value),
                 color = q_value < 0.01),
             size = 0.25, stroke = 0) +
  facet_wrap(~treatment, ncol = 2) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.01), color = "dimgrey", linetype = "dotted",size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("β coefficient") +
  ylab("-log10(q value)") +
  scale_color_manual("Pass\nFilter\n(FDR < 1%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("QC_plots/Effect_size_distribution_NTC_all_treatment_dose_response_diff_test_T98G_Supplementary_Figure_6B.png", 
         dpi = 600, width = 2, height = 1.75)

treatment_dose_dependent_diff_test.list[["U87MG"]] %>%  
  filter(grepl("dose", term)) %>%
  ggplot() +
  geom_point(aes(x = normalized_effect,
                 y = -log10(q_value),
                 color = q_value < 0.01),
             size = 0.25, stroke = 0) +
  facet_wrap(~treatment, ncol = 2) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.01), color = "dimgrey", linetype = "dotted",size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("β coefficient") +
  ylab("-log10(q value)") +
  scale_color_manual("Pass\nFilter\n(FDR < 1%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("QC_plots/Effect_size_distribution_NTC_all_treatment_dose_response_diff_test_U87MG_Supplementary_Figure_6C.png", 
         dpi = 600, width = 2, height = 1.75)

  drug_dose_response_sig_genes.list <- list()

for(cell_type in names(treatment_dose_dependent_diff_test.list)){

  drug_dose_response_sig_genes.list[[cell_type]] <- list()

  for(drug in unique(treatment_dose_dependent_diff_test.list[[cell_type]]$treatment)){

    drug_dose_response_sig_genes.list[[cell_type]][[drug]] <- treatment_dose_dependent_diff_test.list[[cell_type]] %>%
      filter(treatment == drug & grepl("dose",term) &  q_value < 0.01) %>%
      pull(id) %>%
      unique()

    print(length(drug_dose_response_sig_genes.list[[cell_type]][[drug]]))

  }
}

drug_dose_response_sig_genes_df.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_sig_genes_df.list[[cell_type]] <- list()

  for(drug in names(drug_dose_response_sig_genes.list[[cell_type]])){

  drug_dose_response_sig_genes_df.list[[cell_type]][[drug]] <- data.frame(id = drug_dose_response_sig_genes.list[[cell_type]][[drug]],
                                                                          treatment = rep(drug, length(drug_dose_response_sig_genes.list[[cell_type]][[drug]])))    

  }

}

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_sig_genes_df.list[[cell_type]] <- do.call("rbind",drug_dose_response_sig_genes_df.list[[cell_type]])

}

drug_dose_response_sig_genes_UpsetR_df.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_sig_genes_UpsetR_df.list[[cell_type]] <- drug_dose_response_sig_genes_df.list[[cell_type]] %>%
    reshape2::dcast(formula = id ~ treatment, fun.aggregate = length)

}

dir.create("DEG_testing")

png("DEG_testing/A172_DEG_upSetR_plot_Supplementary_Figure_6D.png", width = 7, height = 5, units = "in", res = 1200)
UpSetR::upset(drug_dose_response_sig_genes_UpsetR_df.list[["A172"]],
              sets = c("lapatinib","nintedanib","trametinib","zstk474"),
              order.by = "freq",
              text.scale = 1.75,
              point.size = 2,
              line.size = 0.5)
dev.off()

png("DEG_testing/T98G_DEG_upSetR_plot_Supplementary_Figure_6E.png", width = 7, height = 5, units = "in", res = 1200)
UpSetR::upset(drug_dose_response_sig_genes_UpsetR_df.list[["T98G"]],
              sets = c("lapatinib","nintedanib","trametinib","zstk474"),
              order.by = "freq",
              text.scale = 1.75,
              point.size = 2,
              line.size = .5)
dev.off()

png("DEG_testing/U87MG_DEG_upSetR_plot_Supplementary_Figure_6F.png", width = 7, height = 5, units = "in", res = 1200)
UpSetR::upset(drug_dose_response_sig_genes_UpsetR_df.list[["U87MG"]],
              sets = c("lapatinib","nintedanib","trametinib","zstk474"),
              order.by = "freq",
              text.scale = 1.75,
              point.size = 2,
              line.size = 0.5)
dev.off()

####
### Inspect the dynamics of kinome expression downstream of drug exposure

kinase_protein_to_gene_map <- readRDS("kinase_protein_to_gene_map.rds")

cell_type_kinase_deg.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  message(cell_type,":",(length(intersect(kinase_protein_to_gene_map$id %>% unique(), drug_dose_response_sig_genes_df.list[[cell_type]]$id %>% unique())))," kinase DEGs")

  cell_type_kinase_deg.list[[cell_type]] <- intersect(kinase_protein_to_gene_map$id %>% unique(), drug_dose_response_sig_genes_df.list[[cell_type]]$id %>% unique())
}

cell_type_kinase_degs <- do.call("c",cell_type_kinase_deg.list) %>% unique()

length(intersect(cell_type_kinase_deg.list[["A172"]],cell_type_kinase_deg.list[["T98G"]]))
length(intersect(cell_type_kinase_deg.list[["A172"]],cell_type_kinase_deg.list[["U87MG"]]))
length(intersect(cell_type_kinase_deg.list[["T98G"]],cell_type_kinase_deg.list[["U87MG"]]))

cell_type_perturbing_kinases <- Reduce("intersect",cell_type_kinase_deg.list)
####

drug_dose_response_sig_genes_by_cell_type <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_sig_genes_by_cell_type[[cell_type]] <- do.call("c",drug_dose_response_sig_genes.list[[cell_type]])
  names(drug_dose_response_sig_genes_by_cell_type[[cell_type]]) <- NULL
  drug_dose_response_sig_genes_by_cell_type[[cell_type]] <- unique(drug_dose_response_sig_genes_by_cell_type[[cell_type]])

  print(length(drug_dose_response_sig_genes_by_cell_type[[cell_type]]))
}

drug_dose_response_heatmap_matrices.list <- list()

for(cell_type in names(drug_dose_response_sig_genes_by_cell_type)){

  drug_dose_response_heatmap_matrices.list[[cell_type]] <- list()

  for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

    cds_subset <- NTC_cds.list[[cell_type]][drug_dose_response_sig_genes_by_cell_type[[cell_type]],colData(NTC_cds.list[[cell_type]])$treatment %in% c("vehicle",drug)]
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- t(t(cds_exprs)/colData(cds_subset)$Size_Factor)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    colnames(cds_exprs) <- c("id","cell","expression")

    cD <- colData(cds_subset) %>% as.data.frame()

    cds_exprs <- left_join(cds_exprs,cD,by="cell")
    cds_exprs <- cds_exprs %>%
      group_by(id,dose) %>%
      summarize(mean_expression = mean(expression)) %>%
      mutate(mean_expression = log(mean_expression + 1)) %>%
      tidyr::spread(key = dose, value = mean_expression) %>%
      as.data.frame()
    row.names(cds_exprs) <- cds_exprs$id
    cds_exprs$id <- NULL

    colnames(cds_exprs) <- paste0(drug,"_",colnames(cds_exprs))
    drug_dose_response_heatmap_matrices.list[[cell_type]][[drug]] <- cds_exprs

  }

}

identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["nintedanib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["trametinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["A172"]][["zstk474"]]))

identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["nintedanib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["trametinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["T98G"]][["zstk474"]]))

identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["nintedanib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["trametinib"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["nintedanib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["zstk474"]]))
identical(row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["trametinib"]]),row.names(drug_dose_response_heatmap_matrices.list[["U87MG"]][["zstk474"]]))

drug_dose_response_heatmap_matrices_joint.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_heatmap_matrices_joint.list[[cell_type]] <- do.call("cbind",drug_dose_response_heatmap_matrices.list[[cell_type]])
  colnames(drug_dose_response_heatmap_matrices_joint.list[[cell_type]]) <- sapply(colnames(drug_dose_response_heatmap_matrices_joint.list[[cell_type]]),function(x){stringr::str_split(x,pattern = "\\.")[[1]][2]})

}

length(row.names(drug_dose_response_heatmap_matrices_joint.list[["A172"]]))
length(unique(row.names(drug_dose_response_heatmap_matrices_joint.list[["A172"]])))

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_heatmap_matrices_joint.list[[cell_type]] <- scale(t(scale(t(drug_dose_response_heatmap_matrices_joint.list[[cell_type]]))))
  drug_dose_response_heatmap_matrices_joint.list[[cell_type]][drug_dose_response_heatmap_matrices_joint.list[[cell_type]] > 2] <- 2
  drug_dose_response_heatmap_matrices_joint.list[[cell_type]][drug_dose_response_heatmap_matrices_joint.list[[cell_type]] < -2] <- -2
  drug_dose_response_heatmap_matrices_joint.list[[cell_type]][is.na(drug_dose_response_heatmap_matrices_joint.list[[cell_type]])] <- 0

}

length(row.names(drug_dose_response_heatmap_matrices_joint.list[["A172"]]))
length(unique(row.names(drug_dose_response_heatmap_matrices_joint.list[["A172"]])))

ph.list <- list()

ph.list[["A172"]] <- pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["A172"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 1,
                   height = 3.5,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/A172_dose_dependent_deg_heatmap_Figure_4A.png")

ph.list[["T98G"]] <- pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["T98G"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 1,
                   height = 3.5,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/T98G_dose_dependent_deg_heatmap_Figure_4B.png")

ph.list[["U87MG"]] <- pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["U87MG"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 1,
                   height = 3.5,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/U87MG_dose_dependent_deg_heatmap_Figure_4C.png")

cutree_df.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  cutree_df.list[[cell_type]] <- cutree(ph.list[[cell_type]]$tree_row, 6)
  cutree_df.list[[cell_type]] <-  data.frame(id = names(cutree_df.list[[cell_type]]),
                                             cluster = cutree_df.list[[cell_type]])

}

# saveRDS(cutree_df.list, "NTC_drug_response_cutree_df.list.rds")
cutree_df.list <- readRDS("NTC_drug_response_cutree_df.list.rds")

### Recreate the above heatmaps for gene cluster assignment
pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["A172"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   annotation_row =  data.frame(row.names = cutree_df.list[["A172"]]$id,
                                                cluster = as.factor(cutree_df.list[["A172"]]$cluster)),
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 3,
                   height = 5,
                   treeheight_row = 10,
                   legend = TRUE,
                   file = "Heatmaps/A172_cutree_dose_dependent_deg_heatmap.png")

pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["T98G"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   annotation_row =  data.frame(row.names = cutree_df.list[["T98G"]]$id,
                                                cluster = as.factor(cutree_df.list[["T98G"]]$cluster)),
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 3,
                   height = 5,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/T98G_cutree_dose_dependent_deg_heatmap.png")

pheatmap::pheatmap(drug_dose_response_heatmap_matrices_joint.list[["U87MG"]],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   annotation_row =  data.frame(row.names = cutree_df.list[["U87MG"]]$id,
                                                cluster = as.factor(cutree_df.list[["U87MG"]]$cluster)),
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   width = 3,
                   height = 5,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/U87MG_cutree_dose_dependent_deg_heatmap.png")


### For each gene cluster in the heatmaps above aggregate the mean expression within that cluster
drug_dose_response_summary_df.list <- list()

for(cell_type in names(drug_dose_response_sig_genes.list)){

  drug_dose_response_summary_df.list[[cell_type]] <- as.data.frame(drug_dose_response_heatmap_matrices.list[[cell_type]])
  drug_dose_response_summary_df.list[[cell_type]]$id <- row.names(drug_dose_response_summary_df.list[[cell_type]])
  drug_dose_response_summary_df.list[[cell_type]] <- reshape2::melt(drug_dose_response_summary_df.list[[cell_type]])
  colnames(drug_dose_response_summary_df.list[[cell_type]]) <- c("id","condition","expression")
  
  drug_dose_response_summary_df.list[[cell_type]]$treatment <- sapply(drug_dose_response_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "\\.")[[1]][1]})
  drug_dose_response_summary_df.list[[cell_type]]$dose <- sapply(drug_dose_response_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "_")[[1]][2]})

  print(unique(drug_dose_response_summary_df.list[[cell_type]]$dose))
  drug_dose_response_summary_df.list[[cell_type]] <- left_join(drug_dose_response_summary_df.list[[cell_type]],cutree_df.list[[cell_type]],by="id")

  drug_dose_response_summary_df.list[[cell_type]] <- drug_dose_response_summary_df.list[[cell_type]] %>%
    group_by(treatment,dose,cluster) %>%
    mutate(aggregate_expression = mean(expression))
}

###########################################################################
#### Re-order cluster numbers so they go from top to bottom in Figure ####
###  A172
# order c("5","3","6","4","1","2")
###  T98G
# order c("6","2","3","5","4","1")
####  U87MG
# order c("6","3","2","5","1","4")

for(cell_type in names(drug_dose_response_summary_df.list)){drug_dose_response_summary_df.list[[cell_type]]$cell_type <-  rep(cell_type,dim(drug_dose_response_summary_df.list[[cell_type]])[1])}

### Rename clusters by plot numbers ###

drug_dose_response_summary_df.list[["A172"]]$plot_cluster <- sapply(drug_dose_response_summary_df.list[["A172"]]$cluster,function(x){

  if(x == 5)return("1")
  if(x == 3)return("2")
  if(x == 6)return("3")
  if(x == 4)return("4")
  if(x == 1)return("5")
  if(x == 2)return("6")
  return(NA)

  })

drug_dose_response_summary_df.list[["T98G"]]$plot_cluster <- sapply(drug_dose_response_summary_df.list[["T98G"]]$cluster,function(x){

  if(x == 6)return("1")
  if(x == 2)return("2")
  if(x == 3)return("3")
  if(x == 5)return("4")
  if(x == 4)return("5")
  if(x == 1)return("6")
  return(NA)

  })

drug_dose_response_summary_df.list[["U87MG"]]$plot_cluster <- sapply(drug_dose_response_summary_df.list[["U87MG"]]$cluster,function(x){

  if(x == 6)return("1")
  if(x == 3)return("2")
  if(x == 2)return("3")
  if(x == 5)return("4")
  if(x == 1)return("5")
  if(x == 4)return("6")
  return(NA)

  })

ggplot(drug_dose_response_summary_df.list[["A172"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 1, color = "black", shape  = 21, stroke = 0.05) +
  facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=drug_dose_response_summary_df.list[["A172"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/A172_dose_dependent_deg_heatmap_centroids_Figure_3A_right_panels.png",
    width = 0.6, height = 3, dpi = 600)

ggplot(drug_dose_response_summary_df.list[["T98G"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 1, color = "black", shape  = 21, stroke = 0.05) +
  facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=drug_dose_response_summary_df.list[["T98G"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/T98G_dose_dependent_deg_heatmap_centroids_Figure_3B_right_panels.png",
    width = 0.6, height = 3, dpi = 600)

ggplot(drug_dose_response_summary_df.list[["U87MG"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 1, color = "black", shape  = 21, stroke = 0.05) +
  facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=drug_dose_response_summary_df.list[["U87MG"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/U87MG_dose_dependent_deg_heatmap_centroids_Figure_3B_right_panels.png",
    width = 0.6, height = 3, dpi = 600)



#### Todo Re-do everything associated w marker plots
KinHub_kinome_metadata <- read.table("KinHub_kinome_metadata.txt", sep = "\t", header = TRUE)
kinase_protein_to_gene_map <- readRDS("kinase_protein_to_gene_map.rds")
cutree_df.list <- readRDS("NTC_drug_response_cutree_df.list.rds")

KinHub_kinome_joint_metadata <- left_join(KinHub_kinome_metadata,
                                    kinase_protein_to_gene_map,
                                    by="gene_id")

kinome_degs_cutree.list <- list()

for(cell_type in names(drug_dose_response_summary_df.list)){

  kinome_degs_cutree.list[[cell_type]] <- drug_dose_response_summary_df.list[[cell_type]] %>% 
    filter(id %in% KinHub_kinome_joint_metadata$id)

  kinome_degs_cutree.list[[cell_type]] <- left_join(kinome_degs_cutree.list[[cell_type]],KinHub_kinome_joint_metadata,by="id")

}

### Plot the exrpession of dose-responsive kinases of interest
dir.create("NTC_drug_responsive_kinases")
############
### A172 ###
############

ERBB4_cds <- NTC_cds.list[["A172"]][rowData(NTC_cds.list[["A172"]])$gene_short_name == "ERBB4",]

colData(ERBB4_cds)$treatment_dose <-  paste0(colData(ERBB4_cds)$treatment,"_",colData(ERBB4_cds)$dose)
colData(ERBB4_cds)$treatment_dose <- factor(colData(ERBB4_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(ERBB4_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2"))
  ggsave("NTC_drug_responsive_kinases/ERBB4_expression_by_treatment_A172_Supplementary_Figure_6G_left_panel_1.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

WEE1_cds <- NTC_cds.list[["A172"]][rowData(NTC_cds.list[["A172"]])$gene_short_name == "WEE1",]

colData(WEE1_cds)$treatment_dose <-  paste0(colData(WEE1_cds)$treatment,"_",colData(WEE1_cds)$dose)
colData(WEE1_cds)$treatment_dose <- factor(colData(WEE1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(WEE1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2"))
  ggsave("NTC_drug_responsive_kinases/WEE1_expression_by_treatment_A172_Supplementary_Figure_6G_left_panel_2.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

MAPK10_cds <- NTC_cds.list[["A172"]][rowData(NTC_cds.list[["A172"]])$gene_short_name == "MAPK10",]

colData(MAPK10_cds)$treatment_dose <-  paste0(colData(MAPK10_cds)$treatment,"_",colData(MAPK10_cds)$dose)
colData(MAPK10_cds)$treatment_dose <- factor(colData(MAPK10_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(MAPK10_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2"))
  ggsave("NTC_drug_responsive_kinases/MAPK10_expression_by_treatment_A172_Supplementary_Figure_6G_left_panel_3.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

CDK1_cds <- NTC_cds.list[["A172"]][rowData(NTC_cds.list[["A172"]])$gene_short_name == "CDK1",]

colData(CDK1_cds)$treatment_dose <-  paste0(colData(CDK1_cds)$treatment,"_",colData(CDK1_cds)$dose)
colData(CDK1_cds)$treatment_dose <- factor(colData(CDK1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(CDK1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2"))
  ggsave("NTC_drug_responsive_kinases/CDK1_expression_by_treatment_A172_Supplementary_Figure_6G_left_panel_4.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

IGF1R_cds <- NTC_cds.list[["A172"]][rowData(NTC_cds.list[["A172"]])$gene_short_name == "IGF1R",]

colData(IGF1R_cds)$treatment_dose <-  paste0(colData(IGF1R_cds)$treatment,"_",colData(IGF1R_cds)$dose)
colData(IGF1R_cds)$treatment_dose <- factor(colData(IGF1R_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(IGF1R_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2"))
  ggsave("NTC_drug_responsive_kinases/IGF1R_expression_by_treatment_A172_Supplementary_Figure_6G_left_panel_5.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

############
### T98G ###
############

EPHA5_cds <- NTC_cds.list[["T98G"]][rowData(NTC_cds.list[["T98G"]])$gene_short_name == "EPHA5",]

colData(EPHA5_cds)$treatment_dose <-  paste0(colData(EPHA5_cds)$treatment,"_",colData(EPHA5_cds)$dose)
colData(EPHA5_cds)$treatment_dose <- factor(colData(EPHA5_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(EPHA5_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/EPHA5_expression_by_treatment_T98G_Supplementary_Figure_6G_middle_panel_1.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)  

WEE1_cds <- NTC_cds.list[["T98G"]][rowData(NTC_cds.list[["T98G"]])$gene_short_name == "WEE1",]

colData(WEE1_cds)$treatment_dose <-  paste0(colData(WEE1_cds)$treatment,"_",colData(WEE1_cds)$dose)
colData(WEE1_cds)$treatment_dose <- factor(colData(WEE1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(WEE1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/WEE1_expression_by_treatment_T98G_Supplementary_Figure_6G_middle_panel_2.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)

CDK4_cds <- NTC_cds.list[["T98G"]][rowData(NTC_cds.list[["T98G"]])$gene_short_name == "CDK4",]

colData(CDK4_cds)$treatment_dose <-  paste0(colData(CDK4_cds)$treatment,"_",colData(CDK4_cds)$dose)
colData(CDK4_cds)$treatment_dose <- factor(colData(CDK4_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(CDK4_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/CDK4_expression_by_treatment_T98G_Supplementary_Figure_6G_middle_panel_3.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

JAK1_cds <- NTC_cds.list[["T98G"]][rowData(NTC_cds.list[["T98G"]])$gene_short_name == "JAK1",]

colData(JAK1_cds)$treatment_dose <-  paste0(colData(JAK1_cds)$treatment,"_",colData(JAK1_cds)$dose)
colData(JAK1_cds)$treatment_dose <- factor(colData(JAK1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(JAK1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/JAK1_expression_by_treatment_T98G_Supplementary_Figure_6G_middle_panel_4.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600)  

FGFR1_cds <- NTC_cds.list[["T98G"]][rowData(NTC_cds.list[["T98G"]])$gene_short_name == "FGFR1",]

colData(FGFR1_cds)$treatment_dose <-  paste0(colData(FGFR1_cds)$treatment,"_",colData(FGFR1_cds)$dose)
colData(FGFR1_cds)$treatment_dose <- factor(colData(FGFR1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(FGFR1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/FGFR1_expression_by_treatment_T98G_Supplementary_Figure_6G_middle_panel_5.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

#############
### U87MG ###
#############

ERBB4_cds <- NTC_cds.list[["U87MG"]][rowData(NTC_cds.list[["U87MG"]])$gene_short_name == "ERBB4",]

colData(ERBB4_cds)$treatment_dose <-  paste0(colData(ERBB4_cds)$treatment,"_",colData(ERBB4_cds)$dose)
colData(ERBB4_cds)$treatment_dose <- factor(colData(ERBB4_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(ERBB4_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/ERBB4_expression_by_treatment_U87MG_Supplementary_Figure_6G_right_panel_1.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

WEE1_cds <- NTC_cds.list[["U87MG"]][rowData(NTC_cds.list[["U87MG"]])$gene_short_name == "WEE1",]

colData(WEE1_cds)$treatment_dose <-  paste0(colData(WEE1_cds)$treatment,"_",colData(WEE1_cds)$dose)
colData(WEE1_cds)$treatment_dose <- factor(colData(WEE1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(WEE1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/WEE1_expression_by_treatment_U87MG_Supplementary_Figure_6G_right_panel_2.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

CDK6_cds <- NTC_cds.list[["U87MG"]][rowData(NTC_cds.list[["U87MG"]])$gene_short_name == "CDK6",]

colData(CDK6_cds)$treatment_dose <-  paste0(colData(CDK6_cds)$treatment,"_",colData(CDK6_cds)$dose)
colData(CDK6_cds)$treatment_dose <- factor(colData(CDK6_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(CDK6_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/CDK6_expression_by_treatment_U87MG_Supplementary_Figure_6G_right_panel_3.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

PLK1_cds <- NTC_cds.list[["U87MG"]][rowData(NTC_cds.list[["U87MG"]])$gene_short_name == "PLK1",]

colData(PLK1_cds)$treatment_dose <-  paste0(colData(PLK1_cds)$treatment,"_",colData(PLK1_cds)$dose)
colData(PLK1_cds)$treatment_dose <- factor(colData(PLK1_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(PLK1_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/PLK1_expression_by_treatment_U87MG_Supplementary_Figure_6G_right_panel_4.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

IGF1R_cds <- NTC_cds.list[["U87MG"]][rowData(NTC_cds.list[["U87MG"]])$gene_short_name == "IGF1R",]

colData(IGF1R_cds)$treatment_dose <-  paste0(colData(IGF1R_cds)$treatment,"_",colData(IGF1R_cds)$dose)
colData(IGF1R_cds)$treatment_dose <- factor(colData(IGF1R_cds)$treatment_dose,
  levels = c("vehicle_0","lapatinib_1","lapatinib_10","nintedanib_1","nintedanib_10",
             "trametinib_1","trametinib_10","zstk474_1","zstk474_10"))


plot_percent_cells_positive(IGF1R_cds,
                            group_cells_by = "treatment_dose",
                            bootstrap_samples = 100) +
  theme(text =  element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black",size = 0.25)) +
  scale_fill_manual(values = c("vehicle_0" =  "dimgrey", "lapatinib_1" = "navy", "lapatinib_10" = "navy", 
                               "nintedanib_1" = "darkolivegreen", "nintedanib_10" = "darkolivegreen",
                               "trametinib_1" = "firebrick2", "trametinib_10" = "firebrick2",
                               "zstk474_1" = "darkorange2", "zstk474_10" = "darkorange2")) 
  ggsave("NTC_drug_responsive_kinases/IGF1R_expression_by_treatment_U87MG_Supplementary_Figure_6G_right_panel_5.png",
         width = 0.6,
         height = 0.6,
         dpi  = 600) 

### Inspect the overlap of genes across dose-responsive gene clusters using the Jaccard coefficient

for(cell_type in names(cutree_df.list)){cutree_df.list[[cell_type]]$cell_type <-  rep(cell_type,dim(cutree_df.list[[cell_type]])[1])}

# Rename clusters by plot numbers ###

cutree_df.list[["A172"]]$plot_cluster <- sapply(cutree_df.list[["A172"]]$cluster,function(x){

  if(x == 5)return("1")
  if(x == 3)return("2")
  if(x == 6)return("3")
  if(x == 4)return("4")
  if(x == 1)return("5")
  if(x == 2)return("6")
  return(NA)

  })

cutree_df.list[["T98G"]]$plot_cluster <- sapply(cutree_df.list[["T98G"]]$cluster,function(x){

  if(x == 6)return("1")
  if(x == 2)return("2")
  if(x == 3)return("3")
  if(x == 5)return("4")
  if(x == 4)return("5")
  if(x == 1)return("6")
  return(NA)

  })

cutree_df.list[["U87MG"]]$plot_cluster <- sapply(cutree_df.list[["U87MG"]]$cluster,function(x){

  if(x == 6)return("1")
  if(x == 3)return("2")
  if(x == 2)return("3")
  if(x == 5)return("4")
  if(x == 1)return("5")
  if(x == 4)return("6")
  return(NA)

  })

for(cell_type in names(cutree_df.list)){cutree_df.list[[cell_type]]$cluster_id <-  paste0(cutree_df.list[[cell_type]]$cell_type,"_",cutree_df.list[[cell_type]]$plot_cluster)}

cutree_joint_df <- do.call("rbind",cutree_df.list) %>% distinct()

A172_to_T98G <- expand.grid(unique(cutree_df.list[["A172"]]$cluster_id),unique(cutree_df.list[["T98G"]]$cluster_id))
A172_to_U87MG <- expand.grid(unique(cutree_df.list[["A172"]]$cluster_id),unique(cutree_df.list[["U87MG"]]$cluster_id))
T98G_to_U87MG <- expand.grid(unique(cutree_df.list[["T98G"]]$cluster_id),unique(cutree_df.list[["U87MG"]]$cluster_id))

cluster_comparisons <- do.call("rbind",list(A172_to_T98G,A172_to_U87MG,
                                            T98G_to_U87MG))

colnames(cluster_comparisons) <- c("cluster_id_1","cluster_id_2")
cluster_comparisons$cluster_id_1 <- as.character(cluster_comparisons$cluster_id_1)
cluster_comparisons$cluster_id_2 <- as.character(cluster_comparisons$cluster_id_2)

cluster_comparisons$jaccard_index <- sapply(1:nrow(cluster_comparisons),function(x){

  id_1 <- cluster_comparisons[x,]$cluster_id_1
  id_2 <- cluster_comparisons[x,]$cluster_id_2

  set_1 <- cutree_joint_df %>%
             filter(cluster_id == id_1) %>%
             pull(id) %>%
             unique()

  set_2 <- cutree_joint_df %>%
             filter(cluster_id == id_2) %>%
             pull(id) %>%
             unique()

  jaccard_index <- bayesbio::jaccardSets(set_1,set_2)

  jaccard_index
})

cluster_comparisons$cell_type_1 <- sapply(cluster_comparisons$cluster_id_1,function(x)stringr::str_split(x,pattern="_")[[1]][1])
cluster_comparisons$cell_type_2 <- sapply(cluster_comparisons$cluster_id_2,function(x)stringr::str_split(x,pattern="_")[[1]][1])

cluster_comparisons %>%
  filter(cell_type_1 == "A172" & cell_type_2 ==  "T98G") %>%
  arrange(desc(jaccard_index))

cluster_comparisons %>%
  filter(cell_type_1 == "A172" & cell_type_2 ==  "T98G") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        axis.text = element_blank()) +
  xlab("\nA172 cluster") +
  ylab("T98G cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_A172_to_T98G_Supplementary_Figure_6H_top_panel.png",
         height = 1,
         width = 1.5,
         dpi = 600)

cluster_comparisons %>%
  filter(cell_type_1 == "A172" & cell_type_2 ==  "U87MG") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
    theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        axis.text = element_blank()) +
  xlab("\nA172 cluster") +
  ylab("U87MG cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_A172_to_U87MG_Supplementary_Figure_6H_middle_panel.png",
         height = 1,
         width = 1.5,
         dpi = 600)

cluster_comparisons %>%
  filter(cell_type_1 == "T98G" & cell_type_2 ==  "U87MG") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
    theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        axis.text = element_blank()) +
  xlab("\nT98G cluster") +
  ylab("U87MG cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_T98G_to_U87MG_Supplementary_Figure_6H_bottom_panel.png",
         height = 1,
         width = 1.5,
         dpi = 600)

cluster_comparisons %>%
  filter(cell_type_1 == "A172" & cell_type_2 ==  "T98G") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line")) +
  xlab("\nA172 cluster") +
  ylab("T98G cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_A172_to_T98G_for_axis.png",
         height = 1,
         width = 3,
         dpi = 600)

cluster_comparisons %>%
  filter(cell_type_1 == "A172" & cell_type_2 ==  "U87MG") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
    theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line")) +
  xlab("\nA172 cluster") +
  ylab("U87MG cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_A172_to_U87MG_for_axis.png",
         height = 1,
         width = 3,
         dpi = 600)

cluster_comparisons %>%
  filter(cell_type_1 == "T98G" & cell_type_2 ==  "U87MG") %>%
  mutate(jaccard_index = ifelse(jaccard_index >= 0.1,jaccard_index,0)) %>%
  ggplot() +
  geom_tile(aes(x = factor(cluster_id_1,levels = sort(unique(cluster_id_1))), 
                y = factor(cluster_id_2,levels = sort(unique(cluster_id_2))), 
                fill = jaccard_index)) +
  scale_fill_gradient2(low = "white", high = "red", na.value = "white",
                       name = "Jaccard\nindex") +
  monocle3:::monocle_theme_opts()  +
    theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line")) +
  xlab("\nT98G cluster") +
  ylab("U87MG cluster\n") 
  ggsave("Heatmaps/Jaccard_index_between_NTC_drug_response_clusters_T98G_to_U87MG_for_axis.png",
         height = 1,
         width = 3,
         dpi = 600)

### As a control, permute cluster labels and re-do jaccard calculation
set.seed(2016L)
permutted_cutree_joint_df <- cutree_joint_df %>%
  group_by(cell_type) %>%
  mutate(permuted_cluster_id = sample(cluster_id, replace = FALSE))

cluster_comparisons$permutted_jaccard_index <- sapply(1:nrow(cluster_comparisons),function(x){

  id_1 <- cluster_comparisons[x,]$cluster_id_1
  id_2 <- cluster_comparisons[x,]$cluster_id_2

  set_1 <- permutted_cutree_joint_df %>%
             filter(permuted_cluster_id == id_1) %>%
             pull(id) %>%
             unique()

  set_2 <- permutted_cutree_joint_df %>%
             filter(permuted_cluster_id == id_2) %>%
             pull(id) %>%
             unique()

  jaccard_index <- bayesbio::jaccardSets(set_1,set_2)

  jaccard_index
})

index_df <- cluster_comparisons %>%
  dplyr::select(jaccard_index, permutted_jaccard_index) %>%
  reshape2::melt() 
colnames(index_df) <- c("test", "index")

index_df <- index_df %>%
  mutate(comparison = ifelse(test == "permutted_jaccard_index", "permutted","real")) %>%
  mutate(comparison = factor(comparison, levels = c("real","permutted")))

ggplot(index_df, aes(x = index, fill = comparison)) + 
  geom_density(size = 0.25, alpha = 0.25) +
  scale_fill_manual("Pairwise\ncomparison", values = c("real" = "black", "permutted" = "firebrick3")) +
  geom_vline(xintercept = 0.1, color = "black", linetype = "dashed", size = 0.25) +
  monocle3:::monocle_theme_opts() +
  xlab("Jaccard index") +
  ylab("Density") +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.4,"line"), 
        legend.key.height = unit(0.4,"line"),
        axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("Heatmaps/Jaccard_coefficient_distribution_Supplementary_Figure_6K.png",
    width = 2,
    height = 1.5,
    dpi = 600)

### Collapse clusters with Jaccared coefficients over 0.1
cutree_joint_df$super_cluster_id <- sapply(cutree_joint_df$cluster_id, function(x){

  if(x %in% c("A172_4","T98G_4","U87MG_4"))return("upregulated_super_cluster_1")
  if(x %in% c("A172_3","T98G_3","U87MG_2"))return("downregulated_super_cluster_2")
  return(x)

})

upregulated_super_cluster_1_genes <- cutree_joint_df %>%
  filter(super_cluster_id == "upregulated_super_cluster_1") %>%
  pull(id) %>%
  unique()

downregulated_super_cluster_2_genes <- cutree_joint_df %>%
  filter(super_cluster_id == "downregulated_super_cluster_2") %>%
  pull(id) %>%
  unique()

### Plot the dynamics of converved signatures across drug and dose for each cell line
upregulated_dose_responsive_genes_matrix.list <- list()

for(cell_type in names(NTC_cds.list)){

  upregulated_dose_responsive_genes_matrix.list[[cell_type]] <- list()

  for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

    cds_subset <- NTC_cds.list[[cell_type]][upregulated_dose_responsive_genes,colData(NTC_cds.list[[cell_type]])$treatment %in% c("vehicle",drug)]
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- t(t(cds_exprs)/colData(cds_subset)$Size_Factor)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    colnames(cds_exprs) <- c("id","cell","expression")

    cD <- colData(cds_subset) %>% as.data.frame()

    cds_exprs <- left_join(cds_exprs,cD,by="cell")
    cds_exprs <- cds_exprs %>%
      group_by(id,dose) %>%
      summarize(mean_expression = mean(expression)) %>%
      mutate(mean_expression = log(mean_expression + 1)) %>%
      tidyr::spread(key = dose, value = mean_expression) %>%
      as.data.frame()
    row.names(cds_exprs) <- cds_exprs$id
    cds_exprs$id <- NULL

    colnames(cds_exprs) <- paste0(drug,"_",colnames(cds_exprs))
    upregulated_dose_responsive_genes_matrix.list[[cell_type]][[drug]] <- cds_exprs

  }

}

upregulated_dose_responsive_genes_summary_df.list <- list()

for(cell_type in names(upregulated_dose_responsive_genes_matrix.list)){

  upregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- as.data.frame(upregulated_dose_responsive_genes_matrix.list[[cell_type]])
  upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$id <- row.names(upregulated_dose_responsive_genes_summary_df.list[[cell_type]])
  upregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- reshape2::melt(upregulated_dose_responsive_genes_summary_df.list[[cell_type]])
  colnames(upregulated_dose_responsive_genes_summary_df.list[[cell_type]]) <- c("id","condition","expression")
  
  upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$treatment <- sapply(upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "\\.")[[1]][1]})
  upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$dose <- sapply(upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "_")[[1]][2]})

  print(unique(upregulated_dose_responsive_genes_summary_df.list[[cell_type]]$dose))
  upregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- upregulated_dose_responsive_genes_summary_df.list[[cell_type]] %>%
    group_by(treatment,dose) %>%
    mutate(aggregate_expression = mean(expression))
}

ggplot(upregulated_dose_responsive_genes_summary_df.list[["A172"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  #facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["A172"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/A172_mean_upregulated_dose_responsive_genes_Supplementary_Figure_6L_top_left_panel.png",
    width = 1, height = 0.65, dpi = 600)

ggplot(upregulated_dose_responsive_genes_summary_df.list[["T98G"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  #facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["T98G"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/T98G_mean_upregulated_dose_responsive_genes_Supplementary_Figure_6L_top_middle_panel.png",
    width = 1, height = 0.65, dpi = 600)

ggplot(upregulated_dose_responsive_genes_summary_df.list[["U87MG"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  #facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["U87MG"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/U87MG_mean_upregulated_dose_responsive_genes_Supplementary_Figure_6L_top_right_panel.png",
    width = 1, height = 0.65, dpi = 600)

downregulated_dose_responsive_genes_matrix.list <- list()

for(cell_type in names(NTC_cds.list)){

  downregulated_dose_responsive_genes_matrix.list[[cell_type]] <- list()

  for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

    cds_subset <- NTC_cds.list[[cell_type]][downregulated_dose_responsive_genes,colData(NTC_cds.list[[cell_type]])$treatment %in% c("vehicle",drug)]
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- t(t(cds_exprs)/colData(cds_subset)$Size_Factor)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    colnames(cds_exprs) <- c("id","cell","expression")

    cD <- colData(cds_subset) %>% as.data.frame()

    cds_exprs <- left_join(cds_exprs,cD,by="cell")
    cds_exprs <- cds_exprs %>%
      group_by(id,dose) %>%
      summarize(mean_expression = mean(expression)) %>%
      mutate(mean_expression = log(mean_expression + 1)) %>%
      tidyr::spread(key = dose, value = mean_expression) %>%
      as.data.frame()
    row.names(cds_exprs) <- cds_exprs$id
    cds_exprs$id <- NULL

    colnames(cds_exprs) <- paste0(drug,"_",colnames(cds_exprs))
    downregulated_dose_responsive_genes_matrix.list[[cell_type]][[drug]] <- cds_exprs

  }

}

downregulated_dose_responsive_genes_summary_df.list <- list()

for(cell_type in names(downregulated_dose_responsive_genes_matrix.list)){

  downregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- as.data.frame(downregulated_dose_responsive_genes_matrix.list[[cell_type]])
  downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$id <- row.names(downregulated_dose_responsive_genes_summary_df.list[[cell_type]])
  downregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- reshape2::melt(downregulated_dose_responsive_genes_summary_df.list[[cell_type]])
  colnames(downregulated_dose_responsive_genes_summary_df.list[[cell_type]]) <- c("id","condition","expression")
  
  downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$treatment <- sapply(downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "\\.")[[1]][1]})
  downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$dose <- sapply(downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$condition,
                                                                      function(x){stringr::str_split(x,pattern = "_")[[1]][2]})

  print(unique(downregulated_dose_responsive_genes_summary_df.list[[cell_type]]$dose))e_genes_summary_df.list[[cell_type]],cutree_df.list[[cell_type]],by="id")

  downregulated_dose_responsive_genes_summary_df.list[[cell_type]] <- downregulated_dose_responsive_genes_summary_df.list[[cell_type]] %>%
    group_by(treatment,dose) %>%
    mutate(aggregate_expression = mean(expression))
}

ggplot(downregulated_dose_responsive_genes_summary_df.list[["A172"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["A172"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/A172_mean_downregulated_dose_responsive_genes_Supplementary_Figure_6L_bottom_left_panel.png",
    width = 1, height = 0.65, dpi = 600)

ggplot(downregulated_dose_responsive_genes_summary_df.list[["T98G"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  #facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["T98G"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/T98G_mean_downregulated_dose_responsive_genes_Supplementary_Figure_6L_bottom_middle_panel.png",
    width = 1, height = 0.65, dpi = 600)

ggplot(downregulated_dose_responsive_genes_summary_df.list[["U87MG"]], aes(x = log10(as.numeric(dose) + 0.001), y = aggregate_expression, color = treatment,  fill = treatment)) + 
  geom_line(size = 0.5) +
  geom_point(size = 2, color = "black", shape  = 21, stroke = 0.05) +
  #facet_wrap(~plot_cluster,ncol = 1, scale = "free") +
  scale_x_discrete(labels=upregulated_dose_responsive_genes_summary_df.list[["U87MG"]]$dose) +
  scale_color_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  scale_fill_manual(values = c("lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Scaled aggregate expression") 
  ggsave("Heatmaps/NTC_drug_responses/U87MG_mean_downregulated_dose_responsive_genes_Supplementary_Figure_6L_bottom_right_panel.png",
    width = 1, height = 0.65, dpi = 600)

### Gene set enrichment anbalysis of conserved signatures
hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")
oncogenicSignaturesGSC <- loadGSCSafe(file="c6.all.v6.0.OncogenicSignatures.symbols.gmt")

replace_gene_names_vec <- function(input_vec, name_vec, retain_inds = c(-1,-2)) {
 temp <- merge(name_vec, input_vec, by="row.names")
 temp2 <- temp[,retain_inds]
 names(temp2) <- temp[,2]
 return(temp2)
}

collect_gsa_hyper_results_clusters <- function (genes_list, clusters, gsc) {
    gene_universe <- unique(as.character(genes_list))
    gsa_results <- list()
    cluster_ids <- unique(clusters)

    for (i in cluster_ids) {
        cluster_genes <- unique(names(clusters[clusters == i]))

        gsaRes <- runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe, 
            adjMethod = "BH")
        gsa_results[[length(gsa_results) + 1]] <- gsaRes
    }
    names(gsa_results) <- cluster_ids
    gsa_results
}

gsea_bar_plots <- function(GSAhyper_list, qval_cutoff, pattern, width, height, sample, gsc){
    
    for(cluster in names(GSAhyper_list)){
        
        print(cluster)
        
            GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
            GSAhyper_df$gene_set <- row.names(GSAhyper_df)
            colnames(GSAhyper_df) <- c("qval","gene_set")

            if(is.null(pattern) == FALSE){
                GSAhyper_df$gene_set <- stringr::str_replace(string = GSAhyper_df$gene_set, pattern = pattern, replace = "")
            }

            GSAhyper_df_cutoff <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
            mutate(gene_set = factor(gene_set, levels = gene_set))

            plot_title <- paste0(sample,"_",as.character(cluster),"_",gsc,".png")
            print(plot_title)
            
            ggplot(GSAhyper_df_cutoff, aes(x = gene_set, y = -log10(qval))) + 
            geom_bar(stat = "identity") + 
            coord_flip() +
            theme_classic(base_size = 8) +
            ggsave(plot_title, width = width, height = height)
        
        
        
    }
    
}

expressed_genes <-  Reduce(union,expressed_genes.list)

Ensembl_GSAlist <- as.matrix(rowData(cds.list[["A172"]][expressed_genes,])$gene_short_name)
rownames(Ensembl_GSAlist)<-row.names(rowData(cds.list[["A172"]][expressed_genes,]))
colnames(Ensembl_GSAlist) <- c("gene_short_name")
Ensembl_GSAlist<-Ensembl_GSAlist[,1]
Ensembl_GSAlist<-toupper(Ensembl_GSAlist)
length(Ensembl_GSAlist)

upregulated_cutree_joint_df <- cutree_joint_df %>%
  filter(super_cluster_id == "upregulated_super_cluster_1" &
         id %in% upregulated_super_cluster_1_genes) %>%
  dplyr::select(id,super_cluster_id) %>%
  distinct()

downregulated_cutree_joint_df <- cutree_joint_df %>%
  filter(super_cluster_id == "downregulated_super_cluster_2" &
         id %in% downregulated_super_cluster_2_genes) %>%
  dplyr::select(id,super_cluster_id) %>%
  distinct()

GSA_cutree_joint_df <-  rbind(upregulated_cutree_joint_df,downregulated_cutree_joint_df)

GSA_cutree_joint_df$gene_short_name <- rowData(cds.list[["A172"]][GSA_cutree_joint_df$id,])$gene_short_name
clusters_to_test <- GSA_cutree_joint_df$super_cluster_id
names(clusters_to_test) <- GSA_cutree_joint_df$gene_short_name

hallmarks_GSA <- collect_gsa_hyper_results_clusters(genes_list = Ensembl_GSAlist,
                                                    clusters = clusters_to_test,
                                                    gsc = hallmarksGSC)

gsea_bar_plots(hallmarks_GSA, 
                   qval_cutoff = 0.1, pattern = "HALLMARK_", 
                   width = 8, height = 10, 
                   sample = "GSA/NTC_drug_responses/Jaccard_collapsed_clusters", 
                   gsc = "Hallmarks")

oncogenicSignatures_GSA <- collect_gsa_hyper_results_clusters(genes_list = Ensembl_GSAlist,
                                                    clusters = clusters_to_test,
                                                    gsc = oncogenicSignaturesGSC)

gsea_bar_plots(oncogenicSignatures_GSA, 
                   qval_cutoff = 0.1, pattern = NULL, 
                   width = 8, height = 10, 
                   sample = "GSA/NTC_drug_responses/Jaccard_collapsed_clusters", 
                   gsc = "Oncogenic_Signatures")

### Create UMAP embeddings of NTC cells exposed to drugs and plot expression of signatures of interest
### using signature genes as features. Aim is to highlight relationship between dose, treatment, proliferation
### and other metrics with signature expression.

Jaccard_collapsed_clusters <- GSA_cutree_joint_df

upregulated_dose_responsive_genes <- Jaccard_collapsed_clusters %>%
  filter(super_cluster_id == "upregulated_super_cluster_1") %>%
  pull(id)

downregulated_dose_responsive_genes <- Jaccard_collapsed_clusters %>%
  filter(super_cluster_id == "downregulated_super_cluster_2") %>%
  pull(id)

signature_genes <- c(upregulated_dose_responsive_genes,downregulated_dose_responsive_genes)

for(cell_type in names(cds.list)){

  NTC_cds.list[[cell_type]] <- preprocess_cds(NTC_cds.list[[cell_type]],
                                      method = "PCA",
                                      num_dim = 20,
                                      norm_method = "log",
                                      use_genes = signature_genes)
  
  NTC_cds.list[[cell_type]] <- align_cds(NTC_cds.list[[cell_type]],
                                 residual_model_formula_str = "~replicate",
                                 alignment_group = "replicate")
  
  NTC_cds.list[[cell_type]] <- reduce_dimension(NTC_cds.list[[cell_type]],
                                        max_components = 2,
                                        preprocess_method = "Aligned",
                                        reduction_method = "UMAP",
                                        umap.metric = "cosine",
                                        umap.n_neighbors = 20,
                                        umap.min_dist = 0.1,
                                        umap.fast_sgd=FALSE,
                                        cores=1,
                                        verbose = T)
  
  colData(NTC_cds.list[[cell_type]])$UMAP1 <- reducedDims(NTC_cds.list[[cell_type]])[["UMAP"]][,1]
  colData(NTC_cds.list[[cell_type]])$UMAP2 <- reducedDims(NTC_cds.list[[cell_type]])[["UMAP"]][,2]

}

for(cell_type in names(cds.list)){

  NTC_cds.list[[cell_type]] <- estimate_cell_cycle(NTC_cds.list[[cell_type]],
                                                   g1s_markers = cc.genes$s.genes,
                                                   g2m_markers = cc.genes$g2m.genes)

  colData(NTC_cds.list[[cell_type]])$upregulated_signature_score <- calculate_aggregate_expression_score(NTC_cds.list[[cell_type]],upregulated_dose_responsive_genes, from_id = TRUE)
  colData(NTC_cds.list[[cell_type]])$downregulated_signature_score <- calculate_aggregate_expression_score(NTC_cds.list[[cell_type]],downregulated_dose_responsive_genes, from_id = TRUE)
}

for(cell_type in names(cds.list)){

  NTC_cds.list[[cell_type]] <- estimate_cell_cycle(NTC_cds.list[[cell_type]],
                                                   g1s_markers = cc.genes$s.genes,
                                                   g2m_markers = cc.genes$g2m.genes)

  colData(NTC_cds.list[[cell_type]])$upregulated_signature_score <- calculate_aggregate_expression_score(NTC_cds.list[[cell_type]],upregulated_dose_responsive_genes, from_id = TRUE)
  colData(NTC_cds.list[[cell_type]])$downregulated_signature_score <- calculate_aggregate_expression_score(NTC_cds.list[[cell_type]],downregulated_dose_responsive_genes, from_id = TRUE)
}

ggplot(colData(NTC_cds.list[["A172"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = treatment)) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Treatment", values = c("vehicle" = "grey70", "lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
  ggsave("UMAPs/NTC_UMAP_by_treatment_A172_Figure_4D_top_row.png", dpi = 600, width = 0.5, height = 0.5)


ggplot(colData(NTC_cds.list[["A172"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = as.character(dose))) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Dose", values = c("0" = "grey70", "1" = "#721F81FF", "10" = "#FEAF77FF")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
  ggsave("UMAPs/NTC_UMAP_by_dose_A172_Figure_4D_top_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["A172"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = proliferation_index)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
  ggsave("UMAPs/NTC_UMAP_by_proliferation_index_A172_Figure_4D_top_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["A172"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = upregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Up\nsignature\nscore", option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
  ggsave("UMAPs/NTC_UMAP_by_upregulated_signature_score_A172_Figure_4D_top_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["A172"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = downregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Down\nsignature\nscore") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
  ggsave("UMAPs/NTC_UMAP_by_downregulated_signature_score_A172_Figure_4D_top_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["T98G"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = treatment)) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Treatment", values = c("vehicle" = "grey70", "lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/NTC_UMAP_by_treatment_T98G_Figure_4D_middle_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["T98G"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = as.character(dose))) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Dose", values = c("0" = "grey70", "1" = "#721F81FF", "10" = "#FEAF77FF")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/NTC_UMAP_by_dose_T98G_Figure_4D_middle_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["T98G"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = proliferation_index)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_proliferation_index_T98G_Figure_4D_middle_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["T98G"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = upregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Up\nsignature\nscore", option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_upregulated_signature_score_T98G_Figure_4D_middle_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["T98G"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = downregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Down\nsignature\nscore") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_downregulated_signature_score_T98G_Figure_4D_middle_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = treatment)) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Treatment", values = c("vehicle" = "grey70", "lapatinib" = "navy", "nintedanib" = "darkolivegreen", "trametinib" = "firebrick2", "zstk474" = "darkorange2")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/NTC_UMAP_by_treatment_U87MG_Figure_4D_bottom_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = as.character(dose))) +
  geom_point(size = 0.2, stroke = 0.02) +
  scale_color_manual("Dose", values = c("0" = "grey70", "1" = "#721F81FF", "10" = "#FEAF77FF")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/NTC_UMAP_by_dose_U87MG_Figure_4D_bottom_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = proliferation_index)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_proliferation_index_U87MG_Figure_4D_bottom_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = upregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Up\nsignature\nscore", option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_upregulated_signature_score_U87MG_Figure_4D_bottom_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = downregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Down\nsignature\nscore") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_downregulated_signature_score_U87MG_Figure_4D_bottom_row.png", dpi = 600, width = 0.5, height = 0.5)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = proliferation_index)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Proliferation\nscore") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_proliferation_index_U87MG_for_legend.png", dpi = 600, width = 1, height = 1)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = upregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("UP\nsignature\nscore", option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_upregulated_signature_score_U87MG_for_legend.png", dpi = 600, width = 1, height = 1)

ggplot(colData(NTC_cds.list[["U87MG"]]) %>% as.data.frame(),aes(x = UMAP1, y  = UMAP2, color = downregulated_signature_score)) +
  geom_point(size = 0.2, stroke = 0.02) +
  viridis::scale_color_viridis("Down\nsignature\nscore") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.25,"line")) 
ggsave("UMAPs/NTC_UMAP_by_downregulated_signature_score_U87MG_for_legend.png", dpi = 600, width = 1, height = 1)

