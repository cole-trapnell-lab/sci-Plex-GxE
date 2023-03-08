suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(furrr)
  library(tictoc)
  library(monocle3)})

suppressPackageStartupMessages({
  source("~/sci-plex/bin/cell_cycle.R")
  source("~/sci-plex/bin/loadGSCSafe.R")
  cc.genes = readRDS("~/sci-plex/bin/cc.genes.RDS")
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)
options(stringsAsFactors = FALSE)

#################################################################
#### Functions ####

calculate_mean_target_knockdown <- function(cds, gene_short_name_list, protein_to_gene_map){

  cds_subset <- cds[protein_to_gene_map$id,]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/monocle3:::size_factors(cds_subset))
  cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  
  colnames(cds_exprs) <- c("id","cell","expression")
  cds_exprs$id <- as.character(cds_exprs$id)
  cds_exprs$cell <- as.character(cds_exprs$cell)

  ##### protein_to_gene_map contains ensembl id, HGNC gene_short_name and the protein name 
  cds_exprs <- left_join(cds_exprs,protein_to_gene_map, by = "id")

  Mean_expression_by_target <- list()
  
  # Clean this up to a map function
  targets <- colData(cds) %>% 
    as.data.frame() %>%
    filter(!(gene_id %in% c("NTC","random"))) %>%
    pull(gene_id) %>%
    unique() %>%
    sort()

    for(target in targets){

      target_cells <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == target) %>%
        pull(cell)

      target_cells_2_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == target & gRNA_maxCount >= 2) %>%
        pull(cell)

      target_cells_5_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == target & gRNA_maxCount >= 5) %>%
        pull(cell)

      target_cells_10_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == target & gRNA_maxCount >= 10) %>%
        pull(cell)

      target_exprs <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_cells) %>%
        pull(expression) %>%
        mean()

      target_exprs_2_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_cells_2_cutoff) %>%
        pull(expression) %>%
        mean()

      target_exprs_5_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_cells_5_cutoff) %>%
        pull(expression) %>%
        mean()

      target_exprs_10_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_cells_10_cutoff) %>%
        pull(expression) %>%
        mean()

      NTC_cells <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == "NTC") %>%
        pull(cell)

      NTC_cells_2_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == "NTC" & gRNA_maxCount >= 2) %>%
        pull(cell)

      NTC_cells_5_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == "NTC" & gRNA_maxCount >= 5) %>%
        pull(cell)

      NTC_cells_10_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gene_id == "NTC" & gRNA_maxCount >= 10) %>%
        pull(cell)

      NTC_exprs <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_cells) %>%
        pull(expression) %>%
        mean()

      NTC_exprs_2_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_cells_2_cutoff) %>%
        pull(expression) %>%
        mean()

      NTC_exprs_5_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_cells_5_cutoff) %>%
        pull(expression) %>%
        mean()

      NTC_exprs_10_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_cells_10_cutoff) %>%
        pull(expression) %>%
        mean()

      knockdown <- target_exprs/NTC_exprs
      knockdown_2_cutoff <- target_exprs_2_cutoff/NTC_exprs_2_cutoff
      knockdown_5_cutoff <- target_exprs_5_cutoff/NTC_exprs_5_cutoff
      knockdown_10_cutoff <- target_exprs_10_cutoff/NTC_exprs_10_cutoff
      
      set.seed(2016L)
      target_random_cells <- colData(cds) %>%
        as.data.frame() %>%
        sample_n(length(target_cells)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      target_random_cells_2_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 2) %>%
        sample_n(length(target_cells_2_cutoff)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      target_random_cells_5_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 5) %>%
        sample_n(length(target_cells_5_cutoff)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      target_random_cells_10_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 10) %>%
        sample_n(length(target_cells_10_cutoff)) %>%
        pull(cell) %>%
        as.character()

      target_random_exprs <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_random_cells) %>%
        pull(expression) %>%
        mean()

      target_random_exprs_2_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_random_cells_2_cutoff) %>%
        pull(expression) %>%
        mean()

      target_random_exprs_5_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_random_cells_5_cutoff) %>%
        pull(expression) %>%
        mean()

      target_random_exprs_10_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% target_random_cells_10_cutoff) %>%
        pull(expression) %>%
        mean()

      set.seed(2016L)
      NTC_random_cells <- colData(cds) %>%
        as.data.frame() %>%
        sample_n(length(NTC_cells)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      NTC_random_cells_2_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 2) %>%
        sample_n(length(NTC_cells_2_cutoff)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      NTC_random_cells_5_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 5) %>%
        sample_n(length(NTC_cells_5_cutoff)) %>%
        pull(cell) %>%
        as.character()

      set.seed(2016L)
      NTC_random_cells_10_cutoff <- colData(cds) %>%
        as.data.frame() %>%
        filter(gRNA_maxCount >= 10) %>%
        sample_n(length(NTC_cells_10_cutoff)) %>%
        pull(cell) %>%
        as.character()

      NTC_random_exprs <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_random_cells) %>%
        pull(expression) %>%
        mean()

      NTC_random_exprs_2_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_random_cells_2_cutoff) %>%
        pull(expression) %>%
        mean()

      NTC_random_exprs_5_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_random_cells_5_cutoff) %>%
        pull(expression) %>%
        mean()

      NTC_random_exprs_10_cutoff <- cds_exprs %>% 
        filter(gene_id == target & 
               cell %in% NTC_random_cells_10_cutoff) %>%
        pull(expression) %>%
        mean()

      random_knockdown  <- target_random_exprs/NTC_random_exprs
      random_knockdown_2_cutoff  <- target_random_exprs_2_cutoff/NTC_random_exprs_2_cutoff
      random_knockdown_5_cutoff  <- target_random_exprs_5_cutoff/NTC_random_exprs_5_cutoff
      random_knockdown_10_cutoff  <- target_random_exprs_10_cutoff/NTC_random_exprs_10_cutoff
      
      Mean_expression_by_target[[target]] <- data.frame(gene_id = target,
                                                        mean_knockdown = knockdown,
                                                        mean_random_knockdown = random_knockdown,
                                                        mean_knockdown_2_cutoff = knockdown_2_cutoff,
                                                        mean_random_knockdown_2_cutoff = random_knockdown_2_cutoff,
                                                        mean_knockdown_5_cutoff = knockdown_5_cutoff,
                                                        mean_random_knockdown_5_cutoff = random_knockdown_5_cutoff,
                                                        mean_knockdown_10_cutoff = knockdown_10_cutoff,
                                                        mean_random_knockdown_10_cutoff = random_knockdown_10_cutoff)

      message("finished ",target)
      
    }
  
  Mean_expression <- do.call("rbind",Mean_expression_by_target)

  return(Mean_expression)
}

#################################################################
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

mean_exprs.list <- list()

for(cell_type in names(cds.list)){

  rep1_cells <- colData(cds.list[[cell_type]]) %>% 
    as.data.frame() %>%
    mutate(replicate =  as.character(replicate)) %>%
    filter(replicate ==  "replicate_1") %>%
    pull(cell)

  rep2_cells <- colData(cds.list[[cell_type]]) %>% 
    as.data.frame() %>%
    mutate(replicate =  as.character(replicate)) %>%
    filter(replicate ==  "replicate_2") %>%
    pull(cell)

  rep1_cds <- cds.list[[cell_type]][,rep1_cells]
  rep2_cds <- cds.list[[cell_type]][,rep2_cells]

  cds_exprs_rep1 <- exprs(rep1_cds)
  cds_exprs_rep1 <- Matrix::t(Matrix::t(cds_exprs_rep1)/colData(rep1_cds)$Size_Factor)
  cds_exprs_rep1 <- rowSums(cds_exprs_rep1)
  cds_exprs_rep1 <- data.frame(id = names(cds_exprs_rep1), 
                               rep1_expression = cds_exprs_rep1)

  cds_exprs_rep2 <- exprs(rep2_cds)
  cds_exprs_rep2 <- Matrix::t(Matrix::t(cds_exprs_rep2)/colData(rep2_cds)$Size_Factor)
  cds_exprs_rep2 <- rowSums(cds_exprs_rep2) 
  cds_exprs_rep2 <- data.frame(id = names(cds_exprs_rep2), 
                               rep2_expression = cds_exprs_rep2)

  joint_exprs = left_join(cds_exprs_rep1,cds_exprs_rep2,by="id")
  joint_exprs$cell_line <- rep(cell_type,dim(joint_exprs)[1])

  mean_exprs.list[[cell_type]] <- joint_exprs

  rm(cds_exprs_rep1,cds_exprs_rep2,rep1_cds,rep2_cds,rep1_cells,rep2_cells,joint_exprs)

}

mean_exprs <- do.call("rbind",mean_exprs.list)

ggplot(mean_exprs, aes(x = log10(rep1_expression + 1), y = log10(rep2_expression + 1))) +
  geom_point(size = 0.25, stroke = 0) +
  facet_wrap(~cell_line,  ncol = 3, scale =  "free") +
  geom_smooth(method = "lm", se=TRUE, color="red", formula = y ~ x, size = 0.25) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE, size = 0.8) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  xlab("log10(UMI count + 1)\nReplicate 1 pseudobulk") +
  ylab("log10(UMI count + 1)\nReplicate 2 pseudobulk")
  ggsave("QC_plots/Correlation_between_replicates_Supplementary_Figure_4B.png",
         width = 2.5,
         height = 1.25,
         dpi = 600)

#######################################################################################
#### Summarize effect on cell viability across drugs ####
#######################################################################################
#### Account for differing number of wells per hash plate
# 12 vehicle wells per hash plate 
# 12 10µM exposed wells per treatment per hash plate 
# 9 1µM exposed wells per treatment per hash plate

mean_count_df.list <- list()

for(cell_type in names(cds.list)){

  mean_count_df.list[[cell_type]] <- colData(cds.list[[cell_type]]) %>%
    as.data.frame() %>%
    group_by(gene_id, treatment, dose, top_oligo, replicate) %>%
    summarize(cells_per_well = n()) %>%
    group_by(gene_id, treatment, dose, replicate) %>%
    mutate(mean_cells_per_well = mean(cells_per_well), sd_cells_per_well = sd(cells_per_well)) %>%
    group_by(gene_id, replicate) %>%
    dplyr::select(gene_id, treatment, dose, replicate, mean_cells_per_well, sd_cells_per_well) %>%
    mutate(control_growth = mean_cells_per_well/mean_cells_per_well[treatment == "vehicle"]) %>%
    distinct()

  untreated_df <- mean_count_df.list[[cell_type]] %>% filter(dose == 0)
  lapatinib_untreated_df <- untreated_df %>% mutate(treatment = "lapatinib")
  nintedanib_untreated_df <- untreated_df %>% mutate(treatment = "nintedanib")
  trametinib_untreated_df <- untreated_df %>% mutate(treatment = "trametinib")
  zstk474_untreated_df <- untreated_df %>% mutate(treatment = "zstk474")

  mean_count_df.list[[cell_type]] <- rbind(mean_count_df.list[[cell_type]],
                       lapatinib_untreated_df,
                       nintedanib_untreated_df,
                       trametinib_untreated_df,
                       zstk474_untreated_df)

  mean_count_df.list[[cell_type]] <- mean_count_df.list[[cell_type]] %>% 
    filter(treatment != "vehicle") %>%
    mutate(treatment = factor(treatment, levels = c("lapatinib","nintedanib","trametinib","zstk474")))

  mean_count_df.list[[cell_type]]$cell_line <- rep(cell_type, dim(mean_count_df.list[[cell_type]])[1])  

}

mean_count_df <- do.call("rbind", mean_count_df.list)

max_value = mean_count_df %>%
  filter(gene_id %in% c("NTC","random")) %>%
  pull(control_growth) %>%
  max()

dir.create("Viability_estimates")

mean_count_df %>%
  filter(gene_id %in% c("NTC","random")) %>%
  mutate(treatment = ifelse(treatment == "lapatinib", "lap", ifelse(treatment == "nintedanib", "nin", ifelse(treatment == "trametinib", "tra", "zst"))),
         replicate = ifelse(replicate == "replicate_1", "rep1", "rep2"), 
         gene_id = paste0(gene_id,"_",replicate), 
         gene_id = factor(gene_id, levels = c("random_rep2","random_rep1","NTC_rep2", "NTC_rep1"))) %>%
  ggplot() +
  geom_tile(aes(x = as.character(dose), y = gene_id, fill = control_growth)) +
  facet_wrap(treatment~cell_line, ncol = 3) +
  viridis::scale_fill_viridis("% Control\ngrowth", option = "magma", limits = c(0,max_value)) +
  monocle3:::monocle_theme_opts() +
  xlab("Dose (µM)") +
  ylab("Genotype") +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.75,"line"),
        legend.key.height = unit(0.75,"line"),
        axis.text = element_text(size = 5)) 
  ggsave("Viability_estimates/Drug_associated_viability_estimates_Supplementary_Figure_5B.png", width = 2.5, height = 3, dpi = 600)


#######################################################################################
#### Investigate median knockdown levels across kinase perturbations ####
#######################################################################################
 vehicle_cds.list <- list()

for(cell_type in names(cds.list)){

  vehicle_cds.list[[cell_type]] <- cds.list[[cell_type]][,colData(cds.list[[cell_type]])$treatment == "vehicle" &
                                                          !is.na(colData(cds.list[[cell_type]])$gene_id) &
                                                          colData(cds.list[[cell_type]])$guide_number == 1]

  vehicle_cds.list[[cell_type]] <- estimate_size_factors(vehicle_cds.list[[cell_type]])

  vehicle_cds.list[[cell_type]] <- estimate_cell_cycle(vehicle_cds.list[[cell_type]],
                                                       g1s_markers = cc.genes$s.genes, 
                                                       g2m_markers = cc.genes$g2m.genes)

  ###### Remove cluster assignments from data QC steps
  colData(vehicle_cds.list[[cell_type]])$Cluster <- NULL
  colData(vehicle_cds.list[[cell_type]])$Partition <- NULL

  print(dim(vehicle_cds.list[[cell_type]]))

}

rm(cds.list)

median(mean_target_knockdown.list[["A172"]][!(is.na(mean_target_knockdown.list[["A172"]]$mean_knockdown) | is.infinite(mean_target_knockdown.list[["A172"]]$mean_knockdown)),]$mean_knockdown)
median(mean_target_knockdown.list[["A172"]][!(is.na(mean_target_knockdown.list[["A172"]]$mean_random_knockdown)| is.infinite(mean_target_knockdown.list[["A172"]]$mean_random_knockdown)),]$mean_random_knockdown)

median(mean_target_knockdown.list[["T98G"]][!(is.na(mean_target_knockdown.list[["T98G"]]$mean_knockdown) | is.infinite(mean_target_knockdown.list[["T98G"]]$mean_knockdown)),]$mean_knockdown)
median(mean_target_knockdown.list[["T98G"]][!(is.na(mean_target_knockdown.list[["T98G"]]$mean_random_knockdown) | is.infinite(mean_target_knockdown.list[["T98G"]]$mean_random_knockdown)),]$mean_random_knockdown)

median(mean_target_knockdown.list[["U87MG"]][!(is.na(mean_target_knockdown.list[["U87MG"]]$mean_knockdown) | is.infinite(mean_target_knockdown.list[["U87MG"]]$mean_knockdown)),]$mean_knockdown)
median(mean_target_knockdown.list[["U87MG"]][!(is.na(mean_target_knockdown.list[["U87MG"]]$mean_random_knockdown) | is.infinite(mean_target_knockdown.list[["U87MG"]]$mean_random_knockdown)),]$mean_random_knockdown)

# saveRDS(mean_target_knockdown.list,"mean_target_knockdown.list_updated_110120.rds")
# mean_target_knockdown.list <- readRDS("mean_target_knockdown.list_updated_110120.rds")

mean_target_knockdown_df <- reshape2::melt(mean_target_knockdown.list)
colnames(mean_target_knockdown_df) <- c("gene_id","type","mean_knockdown","cell_line")

mean_target_knockdown_df <- mean_target_knockdown_df %>%
  mutate(cutoff = ifelse(grepl("10",type),"10",ifelse(grepl("5",type),"5",ifelse(grepl("2",type),"2","none"))),
         type = ifelse(grepl("random",type), "random", "real")) %>%
  mutate(cutoff = factor(cutoff, levels = c("none","2","5","10")),
         cell_line =  factor(cell_line, levels = c("A172","T98G","U87MG")),
         type = factor(type, levels = c("real","random"))) %>%
  dplyr::select("cell_line","gene_id","mean_knockdown","cutoff","type")

ggplot(mean_target_knockdown_df %>%  filter(!is.na(mean_knockdown) & cutoff == "10"), 
       aes(x = type, y = mean_knockdown)) +
  geom_bar(stat = "summary", fun = "median", position = "dodge", color = "black", size = 0.1) +
  facet_wrap(~cell_line, scales = "free_y", ncol = 3) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ylab("Median knockdown level\nacross perturbations") +
  xlab("sgRNA assignment") +
  ggsave("QC_plots/Median_knockdown_levels_by_cell_line_Figure_3B.png",
    width = 2.5,  height = 1.75,dpi = 600)

ggpubr::compare_means(mean_knockdown~type, 
                      data = mean_target_knockdown_df %>%  
                              filter(!is.na(mean_knockdown) & cell_line == "A172" & cutoff == "10"))

ggpubr::compare_means(mean_knockdown~type, 
                      data = mean_target_knockdown_df %>%  
                              filter(!is.na(mean_knockdown) & cell_line == "T98G" & cutoff == "10"))

ggpubr::compare_means(mean_knockdown~type, 
                      data = mean_target_knockdown_df %>%  
                              filter(!is.na(mean_knockdown) & cell_line == "U87MG" & cutoff == "10"))

##########################################################
#### Determine the extent of enrichment and depletion of individual kinases in the experiment
##########################################################

plasmid_sgRNA_distribution <- read.table("CROPseq_kinome_sgRNA_plasmid_library_sequencing_results.txt")
colnames(plasmid_sgRNA_distribution) <- c("sequence")
plasmid_sgRNA_distribution$sequence <- as.character(plasmid_sgRNA_distribution$sequence)
plasmid_sgRNA_distribution$protospacer_sequence <- sapply(plasmid_sgRNA_distribution$sequence, function(x){protospacer_sequence <- substr(x,24,42)})

kinome_gRNA_metadata <- read.table("kinome_CRISPRi_guides_gRNASampleSheet.txt", header = FALSE)
colnames(kinome_gRNA_metadata) <- c("gRNA_id","protospacer_sequence")
kinome_gRNA_metadata$gRNA_id <- as.character(kinome_gRNA_metadata$gRNA_id)
kinome_gRNA_metadata$protospacer_sequence <- as.character(kinome_gRNA_metadata$protospacer_sequence)
kinome_gRNA_metadata$gene_id <- sapply(kinome_gRNA_metadata$gRNA_id, function(x){stringr::str_split(x, pattern = "_")[[1]][1]})

kinome_gRNA_metadata$gene_id <- sapply(kinome_gRNA_metadata$gene_id,function(x){gene_id <- ifelse(x %in% c("NTC","originalIDTorder"), "NTC",x)})

plasmid_sgRNA_distribution_annotated <- left_join(plasmid_sgRNA_distribution,kinome_gRNA_metadata, by = "protospacer_sequence")
plasmid_sgRNA_distribution_annotated <- plasmid_sgRNA_distribution_annotated[!is.na(plasmid_sgRNA_distribution_annotated$gene_id),]

plasmid_sgRNA_distribution_summary_df <- plasmid_sgRNA_distribution_annotated %>%
  group_by(gene_id) %>%
  summarize(plasmid_count = n()) %>%
  ungroup() %>%
  mutate(plasmid_frequency = (plasmid_count/sum(plasmid_count))*100)

A172_sgRNA_distribution_summary_df <- colData(vehicle_cds.list[["A172"]]) %>%
  as.data.frame() %>%
  filter(gRNA_maxCount >= 5) %>%
  dplyr::select(gene_id) %>%
  group_by(gene_id) %>%
  summarize(A172_count = n()) %>%
  ungroup() %>%
  mutate(A172_frequency = (A172_count/sum(A172_count))*100)

T98G_sgRNA_distribution_summary_df <- colData(vehicle_cds.list[["T98G"]]) %>%
  as.data.frame() %>%
  filter(gRNA_maxCount >= 5) %>%
  dplyr::select(gene_id) %>%
  group_by(gene_id) %>%
  summarize(T98G_count = n()) %>%
  ungroup() %>%
  mutate(T98G_frequency = (T98G_count/sum(T98G_count))*100)

U87MG_sgRNA_distribution_summary_df <- colData(vehicle_cds.list[["U87MG"]]) %>%
  as.data.frame() %>%
  filter(gRNA_maxCount >= 5) %>%
  dplyr::select(gene_id) %>%
  group_by(gene_id) %>%
  summarize(U87MG_count = n()) %>%
  ungroup() %>%
  mutate(U87MG_frequency = (U87MG_count/sum(U87MG_count))*100)

A172_sgRNA_proportions <- inner_join(A172_sgRNA_distribution_summary_df,plasmid_sgRNA_distribution_summary_df,by = "gene_id")
T98G_sgRNA_proportions <- inner_join(T98G_sgRNA_distribution_summary_df,plasmid_sgRNA_distribution_summary_df,by = "gene_id")
U87MG_sgRNA_proportions <- inner_join(U87MG_sgRNA_distribution_summary_df,plasmid_sgRNA_distribution_summary_df,by = "gene_id")

A172_sgRNA_proportions$relative_proportion <- log(A172_sgRNA_proportions$A172_frequency/A172_sgRNA_proportions$plasmid_frequency) %>% scale()
T98G_sgRNA_proportions$relative_proportion <- log(T98G_sgRNA_proportions$T98G_frequency/T98G_sgRNA_proportions$plasmid_frequency) %>% scale()
U87MG_sgRNA_proportions$relative_proportion <- log(U87MG_sgRNA_proportions$U87MG_frequency/U87MG_sgRNA_proportions$plasmid_frequency) %>% scale()

A172_sgRNA_proportions$cell_line <- rep("A172", dim(A172_sgRNA_proportions)[1])
T98G_sgRNA_proportions$cell_line <- rep("T98G", dim(T98G_sgRNA_proportions)[1])
U87MG_sgRNA_proportions$cell_line <- rep("U87MG", dim(U87MG_sgRNA_proportions)[1])

#sgRNA_proportions_df <- do.call("rbind", list(A172_sgRNA_proportions,T98G_sgRNA_proportions,U87MG_sgRNA_proportions))
#saveRDS(sgRNA_proportions_df, "vehicle_sgRNA_target_proportions_to_plasmid_df.rds")
A172_rank <- A172_sgRNA_proportions %>% arrange(relative_proportion) %>% pull(gene_id)
T98G_rank <- T98G_sgRNA_proportions %>% arrange(relative_proportion) %>% pull(gene_id)
U87MG_rank <- U87MG_sgRNA_proportions %>% arrange(relative_proportion) %>% pull(gene_id)

ggplot(A172_sgRNA_proportions, aes(x = factor(gene_id, levels = A172_rank), y = relative_proportion, label = gene_id, color =  abs(relative_proportion) > 2)) +
geom_point(size = 1, stroke = 0) +
ggrepel::geom_label_repel(data = subset(A172_sgRNA_proportions, abs(relative_proportion) > 2), 
                          size = 1, 
                          segment.size = 0.1,
                          label.size = 0.1,
                          box.padding = 0.1,
                          label.padding = 0.1,
                          color = "black") +
monocle3:::monocle_theme_opts() +
scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
theme(text = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none") +
ylab("Proportion relative to\n plasmid library") +
xlab("Perturbation") +
ggsave("QC_plots/gRNA_distribution_relative_to_plasmid_A172_Supplementary_Figure_5A_top_panel.png", 
  height = 1, width = 1.5, dpi = 600)

ggplot(T98G_sgRNA_proportions, aes(x = factor(gene_id, levels = T98G_rank), y = relative_proportion, label = gene_id, color =  abs(relative_proportion) > 2)) +
geom_point(size = 1, stroke = 0) +
ggrepel::geom_label_repel(data = subset(T98G_sgRNA_proportions, abs(relative_proportion) > 2), 
                          size = 1, 
                          segment.size = 0.1,
                          label.size = 0.1,
                          box.padding = 0.1,
                          label.padding = 0.1,
                          color = "black") +
monocle3:::monocle_theme_opts() +
scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
theme(text = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none") +
ylab("Proportion relative to\n plasmid library") +
xlab("Perturbation") +
ggsave("QC_plots/gRNA_distribution_relative_to_plasmid_T98G_Supplementary_Figure_5A_middle_panel.png", 
  height = 1, width = 1.5, dpi = 600)

ggplot(U87MG_sgRNA_proportions, aes(x = factor(gene_id, levels = U87MG_rank), y = relative_proportion, label = gene_id, color =  abs(relative_proportion) > 2)) +
geom_point(size = 1, stroke = 0) +
ggrepel::geom_label_repel(data = subset(U87MG_sgRNA_proportions, abs(relative_proportion) > 2), 
                          size = 1, 
                          segment.size = 0.1,
                          label.size = 0.1,
                          box.padding = 0.1,
                          label.padding = 0.1,
                          color = "black") +
monocle3:::monocle_theme_opts() +
scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
theme(text = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none") +
ylab("Proportion relative to\n plasmid library") +
xlab("Perturbation") +
ggsave("QC_plots/gRNA_distribution_relative_to_plasmid_U87MG_Supplementary_Figure_5A_bottom_panel.png", 
  height = 1, width = 1.5, dpi = 600)









