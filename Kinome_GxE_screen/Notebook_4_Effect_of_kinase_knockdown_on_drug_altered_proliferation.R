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
cc.genes = readRDS("~/sci-plex/bin/cc.genes.RDS")

#################################################################################
calculate_aggregate_expression_score <- function(cds, signature_genes, from_id = FALSE){

  if(from_id == TRUE){
    cds_subset = cds[rowData(cds)$id %in% signature_genes,]
  }
  else(cds_subset = cds[rowData(cds)$gene_short_name %in% signature_genes,])
  aggregate_signature_expression = exprs(cds_subset)
  aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  aggregate_signature_expression = log(aggregate_signature_expression+1)
  return(aggregate_signature_expression)
}

fit_proliferation_index_models <- function(df_subset,df, drug = drug, pseudocount = 0.01){

  target = unique(df_subset$target)
  df_to_test = df %>% 
    filter(treatment %in% c("vehicle",drug) & gene_id %in% c(target,"NTC")) %>%
    mutate(gene_id = factor(gene_id, 
         levels = c("NTC",target)),
           dose = as.numeric(dose))

  proliferation_index_diff_test <- speedglm::speedglm(data = df_to_test, 
                           formula = proliferation_index ~ log(dose + pseudocount) + gene_id + gene_id:log(dose + pseudocount) + replicate)

  proliferation_index_diff_test <- broom::tidy(proliferation_index_diff_test)
  proliferation_index_diff_test$p.value <- as.numeric(proliferation_index_diff_test$p.value)
  proliferation_index_diff_test$treatment <- rep(drug,dim(proliferation_index_diff_test)[1])
  proliferation_index_diff_test$gene_id <- rep(target,dim(proliferation_index_diff_test)[1])


  return(proliferation_index_diff_test)

}
#################################################################################

### Load dataset. Corresponds to GSM7056149_sciPlexGxE_2_preprocessed_cds.list.RDS in the GEO repository.
cds.list <- readRDS("cds.list.RDS")

### Filter dataset for cells expressing 1 sgRNA and remove MT encoded genes from downstream analyses.
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

### Score cells for the aggregate expression of genes associaed with the cell cycle and proliferation
for(cell_type in names(cds.list)){

  cds.list[[cell_type]] <- estimate_cell_cycle(cds.list[[cell_type]],
                                               g1s_markers = cc.genes$s.genes, 
                                               g2m_markers = cc.genes$g2m.genes)
}

### Identify kinases whose loss leads to a signficant difference in the anti-proliferative response
### induced by targeting the RTK pathway.
options(future.globals.maxSize = 891289600 * 10)
plan(multiprocess, workers = 10)

A172_proliferation_index_test_res.list <- list()
for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){
A172_proliferation_index_test_res.list[[drug]] =
  colData(cds.list[["A172"]]) %>%
  as.data.frame() %>%
  filter(gene_id != "NTC") %>%
  mutate(target = gene_id) %>%
  group_by(gene_id) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = fit_proliferation_index_models,
      df = colData(cds.list[["A172"]]) %>% as.data.frame(),
      drug=drug,
      .progress = FALSE
    )
  )
  print(paste("finished",drug))
}

# saveRDS(A172_proliferation_index_test_res.list, "A172_proliferation_index_test_res.list.rds")

T98G_proliferation_index_test_res.list <- list()
for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){
T98G_proliferation_index_test_res.list[[drug]] =
  colData(cds.list[["T98G"]]) %>%
  as.data.frame() %>%
  filter(gene_id != "NTC") %>%
  mutate(target = gene_id) %>%
  group_by(gene_id) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = fit_proliferation_index_models,
      df = colData(cds.list[["T98G"]]) %>% as.data.frame(),
      drug=drug,
      .progress = FALSE
    )
  )
  print(paste("finished",drug))
}

# saveRDS(T98G_proliferation_index_test_res.list, "T98G_proliferation_index_test_res.list.rds")

U87MG_proliferation_index_test_res.list <- list()
for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){
U87MG_proliferation_index_test_res.list[[drug]] =
  colData(cds.list[["U87MG"]]) %>%
  as.data.frame() %>%
  filter(gene_id != "NTC") %>%
  mutate(target = gene_id) %>%
  group_by(gene_id) %>%
  nest() %>%
  mutate(
    test_res = furrr::future_map(
      data,
      .f = fit_proliferation_index_models,
      df = colData(cds.list[["U87MG"]]) %>% as.data.frame(),
      drug=drug,
      .progress = FALSE
    )
  )
  print(paste("finished",drug))
}

# saveRDS(U87MG_proliferation_index_test_res.list, "U87MG_proliferation_index_test_res.list.rds")

A172_proliferation_index_test_res.list <- readRDS("A172_proliferation_index_test_res.list.rds")
T98G_proliferation_index_test_res.list <- readRDS("T98G_proliferation_index_test_res.list.rds")
U87MG_proliferation_index_test_res.list <- readRDS("U87MG_proliferation_index_test_res.list.rds")

for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

  A172_proliferation_index_test_res.list[[drug]] <- A172_proliferation_index_test_res.list[[drug]] %>% ungroup() %>% dplyr::select(test_res)
  T98G_proliferation_index_test_res.list[[drug]] <- T98G_proliferation_index_test_res.list[[drug]] %>% ungroup() %>% dplyr::select(test_res)
  U87MG_proliferation_index_test_res.list[[drug]] <- U87MG_proliferation_index_test_res.list[[drug]] %>% ungroup() %>% dplyr::select(test_res)
  
}

A172_proliferation_index_test_res_df <- do.call("rbind",A172_proliferation_index_test_res.list)
T98G_proliferation_index_test_res_df <- do.call("rbind",T98G_proliferation_index_test_res.list)
U87MG_proliferation_index_test_res_df <- do.call("rbind",U87MG_proliferation_index_test_res.list)

A172_proliferation_index_test_res_df <- do.call("rbind",A172_proliferation_index_test_res_df$test_res)
T98G_proliferation_index_test_res_df <- do.call("rbind",T98G_proliferation_index_test_res_df$test_res)
U87MG_proliferation_index_test_res_df <- do.call("rbind",U87MG_proliferation_index_test_res_df$test_res)

A172_proliferation_index_test_res_df$cell_line <- rep("A712", dim(A172_proliferation_index_test_res_df)[1])
T98G_proliferation_index_test_res_df$cell_line <- rep("T98G", dim(T98G_proliferation_index_test_res_df)[1])
U87MG_proliferation_index_test_res_df$cell_line <- rep("U87MG", dim(U87MG_proliferation_index_test_res_df)[1])

A172_proliferation_index_test_res_df$q_value <- p.adjust(A172_proliferation_index_test_res_df$p.value, method = "BH")
T98G_proliferation_index_test_res_df$q_value <- p.adjust(T98G_proliferation_index_test_res_df$p.value, method = "BH")
U87MG_proliferation_index_test_res_df$q_value <- p.adjust(U87MG_proliferation_index_test_res_df$p.value, method = "BH")

A172_PI_perturbing_kinases <- A172_proliferation_index_test_res_df %>% 
	filter(grepl(":gene", term), q_value < 0.05) %>%
	pull(gene_id) %>%
	unique()

T98G_PI_perturbing_kinases <- T98G_proliferation_index_test_res_df %>% 
	filter(grepl(":gene", term), q_value < 0.05) %>%
	pull(gene_id) %>%
	unique()

U87MG_PI_perturbing_kinases <- U87MG_proliferation_index_test_res_df %>% 
	filter(grepl(":gene", term), q_value < 0.05) %>%
	pull(gene_id) %>%
	unique()

proliferation_index_test_res_df <- Reduce("rbind",list(A172_proliferation_index_test_res_df,T98G_proliferation_index_test_res_df,U87MG_proliferation_index_test_res_df))
# saveRDS(proliferation_index_test_res_df, "proliferation_index_test_res_df.rds")

PI_perturbing_kinases <- Reduce("union",list(A172_PI_perturbing_kinases,T98G_PI_perturbing_kinases,U87MG_PI_perturbing_kinases)) %>% unique()

### Exclude kinases with very low numbers of cells per genotype per condition (i.e., drug, dose combinations).
cell_number_summary.list <- list()

for(cell_type in names(cds.list)){

	cell_number_summary.list[[cell_type]] <- colData(cds.list[[cell_type]]) %>%
		as.data.frame() %>%
		group_by(treatment, dose, gene_id) %>%
		summarize(total_cells = n()) %>%
		arrange(total_cells)


  	}

cells_over_cutoff_A172 <- cell_number_summary.list[["A172"]] %>% 
	filter(total_cells > 5) %>% 
	group_by(gene_id) %>% 
	summarize(total_bins = n()) %>% 
	filter(total_bins == 9) %>% pull(gene_id)

cells_over_cutoff_T98G <- cell_number_summary.list[["T98G"]] %>% 
	filter(total_cells > 5) %>% 
	group_by(gene_id) %>% 
	summarize(total_bins = n()) %>% 
	filter(total_bins == 9) %>% pull(gene_id)

cells_over_cutoff_U87MG <- cell_number_summary.list[["U87MG"]] %>% 
	filter(total_cells > 5) %>% 
	group_by(gene_id) %>% 
	summarize(total_bins = n()) %>% 
	filter(total_bins == 9) %>% pull(gene_id)

genotypes_over_cutoff <- Reduce("intersect",list(cells_over_cutoff_A172,cells_over_cutoff_T98G,cells_over_cutoff_U87MG))
PI_perturbing_kinases_over_cutoff <- PI_perturbing_kinases[PI_perturbing_kinases %in% genotypes_over_cutoff]

### Generate pseudobulk heatmaps of kinases whose loss leads to a signficant difference in the anti-proliferative response
### induced by targeting the RTK pathway.
drug_dose_proliferation_index_heatmap_matrices.list.list <- list()

for(cell_type in names(cds.list)){

  drug_dose_proliferation_index_heatmap_matrices.list.list[[cell_type]] <- list()

  for(drug in c("lapatinib","nintedanib","trametinib","zstk474")){

    df_subset <- colData(cds.list[[cell_type]]) %>%
		as.data.frame() %>%
		filter(treatment %in% c("vehicle",drug)) %>%
		group_by(dose, gene_id) %>%
		summarize(mean_proliferation_index = mean(proliferation_index), total_cells = n()) %>%
		group_by(gene_id) %>%
		mutate(mean_proliferation_index = mean_proliferation_index-mean_proliferation_index[dose == 0]) %>%
		dplyr:::select(dose, gene_id, mean_proliferation_index) %>%
		tidyr::spread(key = dose, value = mean_proliferation_index) %>%
		as.data.frame()

    row.names(df_subset) <- df_subset$gene_id
    df_subset$gene_id <- NULL

    colnames(df_subset) <- paste0(cell_type,"_",drug,"_",colnames(df_subset))
    drug_dose_proliferation_index_heatmap_matrices.list.list[[cell_type]][[drug]] <- df_subset

  }

}

# Check that the order of the rows is correct
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["nintedanib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["trametinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["A172"]][["zstk474"]]))

identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["nintedanib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["trametinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["T98G"]][["zstk474"]]))

identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["nintedanib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["lapatinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["trametinib"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["nintedanib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["zstk474"]]))
identical(row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["trametinib"]]),row.names(drug_dose_proliferation_index_heatmap_matrices.list.list[["U87MG"]][["zstk474"]]))


drug_dose_proliferation_index_heatmap_matrices_joint.list <- list()

for(cell_type in names(drug_dose_proliferation_index_heatmap_matrices.list.list)){

  drug_dose_proliferation_index_heatmap_matrices_joint.list[[cell_type]] <- do.call("cbind",drug_dose_proliferation_index_heatmap_matrices.list.list[[cell_type]])
  colnames(drug_dose_proliferation_index_heatmap_matrices_joint.list[[cell_type]]) <- sapply(colnames(drug_dose_proliferation_index_heatmap_matrices_joint.list[[cell_type]]),function(x){stringr::str_split(x,pattern = "\\.")[[1]][2]})

}

length(row.names(drug_dose_proliferation_index_heatmap_matrices_joint.list[["A172"]]))
length(unique(row.names(drug_dose_proliferation_index_heatmap_matrices_joint.list[["A172"]])))

length(row.names(drug_dose_proliferation_index_heatmap_matrices_joint.list[["A172"]]))
length(unique(row.names(drug_dose_proliferation_index_heatmap_matrices_joint.list[["A172"]])))

drug_dose_proliferation_index_heatmap_matrices_joint_df <- do.call("cbind",drug_dose_proliferation_index_heatmap_matrices_joint.list)

paletteLength <- 35
hmcols<-colorRampPalette(c("navy","white","firebrick2"))(paletteLength)

PI_Breaks <- c(seq(min(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),][!is.na(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),])]), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),][!is.na(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),])])/paletteLength, max(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),][!is.na(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases),])]), length.out=floor(paletteLength/2)))

proliferation_index_ph.list <- list()

### Examine the grouping of drug-induced changes in proliferation perturbing kinases.
proliferation_index_ph.list[["joint"]] <- pheatmap::pheatmap(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases_over_cutoff),],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = TRUE,
                   color = viridis::viridis(35),
                   #breaks = PI_Breaks,
                   gaps_col = c(3,6,9,12,15,18,21,24,27,30,33),
                   cutree_rows = 6,
                   fontsize = 3,
                   fontsize_row = 3,
                   width = 2.5,
                   height = 4,
                   treeheight_row = 10,
                   legend = FALSE,
                   file = "Heatmaps/Genotype_proliferation_index_by_drug/Joint.png")

proliferation_index_ph.list[["joint_for_legend"]] <- pheatmap::pheatmap(drug_dose_proliferation_index_heatmap_matrices_joint_df[c("NTC","random",PI_perturbing_kinases_over_cutoff),],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = TRUE,
                   scale = "none",
                   color = viridis::viridis(35),
                   #breaks = PI_Breaks,
                   gaps_col = c(3,6,9,12,15,18,21,24,27,30,33),
                   cutree_rows = 6,
                   fontsize = 6,
                   fontsize_row = 6,
                   width = 5,
                   height = 8,
                   treeheight_row = 10,
                   legend = TRUE,
                   file = "Heatmaps/Genotype_proliferation_index_by_drug/Joint_for_legend.png")

# Hard to get good resolution on the joint heatmap, instead made each heatmaps individual for each line and ordered by the dendrogram of the joint analysis.
proliferation_index_ph.list[["A172"]] <- pheatmap::pheatmap(drug_dose_proliferation_index_heatmap_matrices_joint.list[["A172"]][c("NTC","random",PI_perturbing_kinases_over_cutoff),],
                   #clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   cluster_rows = proliferation_index_ph.list[["joint"]]$tree_row,
                   show_colnames = FALSE,
                   show_rownames = TRUE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   fontsize_row = 5,
                   width = 1.25,
                   height = 6,
                   treeheight_row = 5,
                   legend = FALSE,
                   file = "Heatmaps/Genotype_proliferation_index_by_drug/A172_PI_heatmap_Figure_3F_top_panel.png")

proliferation_index_ph.list[["T98G"]] <- pheatmap::pheatmap(drug_dose_proliferation_index_heatmap_matrices_joint.list[["T98G"]][c("NTC","random",PI_perturbing_kinases_over_cutoff),],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   cluster_rows = proliferation_index_ph.list[["joint"]]$tree_row,
                   show_colnames = FALSE,
                   show_rownames = TRUE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   fontsize_row = 5,
                   width = 1.25,
                   height = 6,
                   treeheight_row = 5,
                   legend = FALSE,
                   file = "Heatmaps/Genotype_proliferation_index_by_drug/T98G_PI_heatmap_Figure_3F_middle_panel.png")

proliferation_index_ph.list[["U87MG"]] <- pheatmap::pheatmap(drug_dose_proliferation_index_heatmap_matrices_joint.list[["U87MG"]][c("NTC","random",PI_perturbing_kinases_over_cutoff),],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   cluster_rows = proliferation_index_ph.list[["joint"]]$tree_row,
                   show_colnames = FALSE,
                   show_rownames = TRUE,
                   color = viridis::viridis(35),
                   gaps_col = c(3,6,9),
                   cutree_rows = 6,
                   fontsize = 6,
                   fontsize_row = 5,
                   width = 1.25,
                   height = 6,
                   treeheight_row = 5,
                   legend = FALSE,
                   file = "Heatmaps/Genotype_proliferation_index_by_drug/U87MG_PI_heatmap_Figure_3F_bottom_panel.png")


### Examine the effect size of interaction coefficients across all treatments for each line
A172_PI_perturbing_kinases <- proliferation_index_test_res_df %>%
  mutate(estimate = as.numeric(estimate)) %>%
  filter(grepl(":gene_id",term), cell_line == "A712", gene_id %in% PI_perturbing_kinases) %>%  
  dplyr::select(gene_id,treatment,estimate) %>%
  distinct()

A172_PI_perturbing_kinases_df <- A172_PI_perturbing_kinases %>%
	tidyr::spread(key = treatment, value = estimate) %>%
	as.data.frame() 

row.names(A172_PI_perturbing_kinases_df) <- A172_PI_perturbing_kinases_df$gene_id
A172_PI_perturbing_kinases_df$gene_id <- NULL

png("Heatmaps/Genotype_proliferation_index_by_drug/A172_PI_perturbing_kinases_beta_coefficients_heatmap_Figure_3E_top_left_panel.png", res = 600, pointsize = 6, units = "in", width = 1.5, height = 1.5)
circlize::circos.heatmap(A172_PI_perturbing_kinases_df[row.names(A172_PI_perturbing_kinases_df) %in% genotypes_over_cutoff,] %>% as.matrix(), 
			   col = circlize::colorRamp2(c(min(A172_PI_perturbing_kinases_df), 0, max(A172_PI_perturbing_kinases_df)), c("blue", "white", "red")), 
			   dend.side = "inside",
			   rownames.side = "outside",
			   dend.track.height = 0.2)

circlize::circos.clear()
dev.off()

T98G_PI_perturbing_kinases <- proliferation_index_test_res_df %>%
  mutate(estimate = as.numeric(estimate)) %>%
  filter(grepl(":gene_id",term), cell_line == "T98G", gene_id %in% PI_perturbing_kinases) %>%  
  dplyr::select(gene_id,treatment,estimate) %>%
  distinct()

T98G_PI_perturbing_kinases_df <- T98G_PI_perturbing_kinases %>%
  tidyr::spread(key = treatment, value = estimate) %>%
  as.data.frame() 

row.names(T98G_PI_perturbing_kinases_df) <- T98G_PI_perturbing_kinases_df$gene_id
T98G_PI_perturbing_kinases_df$gene_id <- NULL

png("Heatmaps/Genotype_proliferation_index_by_drug/T98G_PI_perturbing_kinases_beta_coefficients_heatmap_Figure_3E_top_right_panel.png", res = 600, pointsize = 6, units = "in", width = 1.5, height = 1.5)
circlize::circos.heatmap(T98G_PI_perturbing_kinases_df[row.names(T98G_PI_perturbing_kinases_df) %in% genotypes_over_cutoff,] %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(T98G_PI_perturbing_kinases_df), 0, max(T98G_PI_perturbing_kinases_df)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
			   dend.track.height = 0.2)
circlize::circos.clear()
dev.off()

U87MG_PI_perturbing_kinases <- proliferation_index_test_res_df %>%
  mutate(estimate = as.numeric(estimate)) %>%
  filter(grepl(":gene_id",term), cell_line == "U87MG", gene_id %in% PI_perturbing_kinases) %>%  
  dplyr::select(gene_id,treatment,estimate) %>%
  distinct()

U87MG_PI_perturbing_kinases_df <- U87MG_PI_perturbing_kinases %>%
  tidyr::spread(key = treatment, value = estimate) %>%
  as.data.frame() 

row.names(U87MG_PI_perturbing_kinases_df) <- U87MG_PI_perturbing_kinases_df$gene_id
U87MG_PI_perturbing_kinases_df$gene_id <- NULL

png("Heatmaps/Genotype_proliferation_index_by_drug/U87MG_PI_perturbing_kinases_beta_coefficients_heatmap_Figure_3E_bottom_middle_panel.png", res = 600, pointsize = 6, units = "in", width = 1.5, height = 1.5)
circlize::circos.heatmap(U87MG_PI_perturbing_kinases_df[row.names(U87MG_PI_perturbing_kinases_df) %in% genotypes_over_cutoff,] %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(U87MG_PI_perturbing_kinases_df), 0, max(U87MG_PI_perturbing_kinases_df)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
			   dend.track.height = 0.2)
circlize::circos.clear()
dev.off()

png("Heatmaps/Genotype_proliferation_index_by_drug/PI_perturbing_kinases_beta_coefficients_heatmap_nnotation_legend.png", res = 600, pointsize = 6, units = "in", width = 1.5, height = 1.5)
lgd = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
							 grid_height = unit(4, "mm"),grid_width = unit(4, "mm"),)
grid::grid.draw(lgd)
dev.off()

### Plot the expression on proliferation of a few kinases of interest
joint_colData_df <- do.call("rbind",list(colData(cds.list[["A172"]]),colData(cds.list[["T98G"]]),colData(cds.list[["U87MG"]]))) %>% as.data.frame()
joint_colData_df <- joint_colData_df %>%
	group_by(cell_line,gene_id) %>%
	mutate(proliferation_index = proliferation_index - mean(proliferation_index[dose == 0])) %>% 
	mutate(dose = factor(as.character(dose), levels = c("0","1","10")))

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","SMG1"), treatment %in% c("vehicle","nintedanib"), cell_line == "A172") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","SMG1"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "SMG1" = "deepskyblue4")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/A172_SMG1_nintedanib_Figure_3D_1st_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","SMG1"), treatment %in% c("vehicle","trametinib"), cell_line == "A172") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","SMG1"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "SMG1" = "deepskyblue4")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/A172_SMG1_trametinib_Figure_3D_2nd_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","STK36"), treatment %in% c("vehicle","nintedanib"), cell_line == "A172") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","STK36"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "STK36" = "deepskyblue4")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/A172_STK36_nintedanib_Figure_3D_3rd_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","STK36"), treatment %in% c("vehicle","trametinib"), cell_line == "A172") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","STK36"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "STK36" = "deepskyblue4")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/A172_STK36_trametinib_Figure_3D_4th_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","STK11"), treatment %in% c("vehicle","nintedanib"), cell_line == "U87MG") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","STK11"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "STK11" = "firebrick2")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/U87MG_STK11_nintedanib_Figure_3E_3rd_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","MAP2K4"), treatment %in% c("vehicle","nintedanib"), cell_line == "U87MG") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","MAP2K4"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "MAP2K4" = "firebrick2")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/U87MG_MAP2K4_nintedanib_Figure_3E_4th_panel.png", width = 1.25, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","ACVR1B","ACVR2A"), treatment %in% c("vehicle","trametinib"), cell_line == "T98G") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","ACVR1B","ACVR2A"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "ACVR1B" = "brown4", "ACVR2A" = "firebrick2")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/T98G_ACVR1B_ACVR2A_trametinib_Figure_3E_1st_panel.png", width = 1.5, height = 1.25, dpi = 600)

ggplot(joint_colData_df %>% 
	filter(gene_id %in% c("NTC","random","CDK18"), treatment %in% c("vehicle","zstk474"), cell_line == "T98G") %>%
	mutate(gene_id = factor(gene_id, levels = c("NTC","random","CDK18"))), aes(x = gene_id, y = proliferation_index, fill = gene_id)) +
	geom_violin(scale = "width", size = 0.1) + 
	facet_wrap(~dose, ncol = 3) +
	scale_fill_manual(values = c("NTC" = "grey70", "random" = "dimgrey", "CDK18" = "firebrick2")) +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	stat_summary(fun.data=mean_sdl, 
	             geom="point", color="black", size = 0.75, stroke = 0) +
	monocle3:::monocle_theme_opts() +
	theme(text = element_text(size = 6),
		  legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
	xlab("Genotype") +
	ylab("Proliferation index") +
	ggsave("Heatmaps/Genotype_proliferation_index_by_drug/T98G_CDK18_zstk474_Figure_3E_2nd_panel.png", width = 1.25, height = 1.25, dpi = 600)

