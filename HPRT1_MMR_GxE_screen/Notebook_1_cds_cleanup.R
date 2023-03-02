library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(devtools)
library(monocle3)

### Start from the hash annotated cds that is output by process_sci-Plex_sci-Plex-GxE_pipeline.sh
cds <- readRDS("path_to_cds")

### Remove cells with no hashes
cds <- cds[,!is.na(colData(cds)$top_oligo_W)]

### Clean up and add new metadata based on the assigned hash
names(colData(cds))[names(colData(cds)) == 'Cell'] <- 'cell'
colData(cds)$new_cell <- sapply(colData(cds)$cell,function(x){substr(x, start = 5, stop = nchar(x))})
  
CRISPR <- sapply(colData(cds)$top_oligo_W,function(x){stringr::str_split(x, pattern = "_")[[1]][4]})
colData(cds)$CRISPR_hash <- CRISPR

gRNA_library <- sapply(colData(cds)$top_oligo_W,function(x){stringr::str_split(x, pattern = "_")[[1]][5]})
colData(cds)$gRNA_library <- gRNA_library
  
dose <- sapply(colData(cds)$top_oligo_W,function(x){stringr::str_split(x, pattern = "_")[[1]][6]})
colData(cds)$dose <- as.numeric(dose)
  
treatment <- sapply(colData(cds)$top_oligo_W,function(x){stringr::str_split(x, pattern = "_")[[1]][7]})
colData(cds)$treatment <- treatment
  
time_point <-  sapply(colData(cds)$top_oligo_W,function(x){stringr::str_split(x, pattern = "_")[[1]][8]})
colData(cds)$time_point <- time_point

### Incorporate sgRNA enrichment results output by process_sci-Plex_sci-Plex-GxE_pipeline.sh, can be found on GEO repository
gRNA_enrichment_results <- readRDS("path_to_sciPlexGxE_1_sgRNATable_out.rds")

# Read in sgRNA metadata, can be found on GEO repository
CRISPR_library <- read.table("sciPlexGxE_1_sgRNA_sequences.txt.gz", header = T)
CRISPR_library$gene <- as.character(CRISPR_library$gene)
CRISPR_library$CRISPR <- as.character(CRISPR_library$CRISPR)
CRISPR_library$gRNA_sequence <- as.character(CRISPR_library$gRNA_sequence)
CRISPR_library$rank <- as.character(CRISPR_library$rank)
CRISPR_library$gRNA_sequence <- substr(CRISPR_library$gRNA_sequence, start = 2, stop = 20)

dim(gRNA_enrichment_results)
dim(gRNA_enrichment_results[gRNA_enrichment_results$gRNA %in% CRISPR_library$gRNA,])

CRISPR_library$gene <- sapply(CRISPR_library$gene, function(x){
  if(x == "negative_control")return("NTC")
  return(x)
  })

colnames(CRISPR_library) <- c("gene_id","CRISPR_gRNA","gRNA","gRNA_rank")

CRISPR_library$gRNA_id <- paste0(CRISPR_library$gene,"_",as.character(CRISPR_library$gRNA_rank))
CRISPR_library <- CRISPR_library %>% dplyr::select(gRNA_id,gene_id,gRNA,CRISPR_gRNA,gRNA_rank)

# Add a key to align sgRNA enrichment and gene expression data, then calculate sgRNA enrichment
gRNA_enrichment_results_subset <- gRNA_enrichment_results[gRNA_enrichment_results$new_cell %in% colData(cds)$new_cell,]
gRNA_enrichment_results_subset$duplication_rate <-  1-(gRNA_enrichment_results_subset$umi_count/gRNA_enrichment_results_subset$read_count)

median(gRNA_enrichment_results_subset[gRNA_enrichment_results_subset$read_count > 1,]$duplication_rate)

gRNA_enrichment_results_summary <- gRNA_enrichment_results_subset %>% 
  group_by(new_cell,gRNA) %>% 
  mutate(total_read_count = sum(read_count), total_umi_count = sum(umi_count))  %>%
  dplyr::select(new_cell,gRNA, total_read_count, total_umi_count) %>% 
  distinct() 

gRNA_proportions <- gRNA_enrichment_results_summary %>% ungroup() %>%
  dplyr::filter(gRNA %in% CRISPR_library$gRNA) %>%
  dplyr::select(new_cell,gRNA,total_read_count) %>%
  group_by(new_cell) %>%
  mutate(total = sum(total_read_count)) %>%
  arrange(new_cell, -total_read_count) %>%
  dplyr::slice(1:5) %>%
  mutate(rank = row_number(),
         ratio = total_read_count/total,
         ratio1 = max(total_read_count)/total) %>%
  mutate(ratio2 = ifelse(length(rank) > 1,total_read_count[rank == 2]/total,0),
         ratio3 = ifelse(length(rank) > 2,total_read_count[rank == 3]/total,0)) %>%
  dplyr::select(everything(), maxCount = total_read_count)

gRNA_proportions <- gRNA_proportions %>% 
mutate(top_to_second_best_ratio = ratio1/ratio2,
  top_to_third_best_ratio = ratio1/ratio3,
  second_to_third_best_ratio = ratio2/ratio3) %>%
ungroup()

dir.create("gRNA_QC_plots")

ggplot(gRNA_proportions[gRNA_proportions$ratio > 0,], aes(x = ratio, y= log10(total))) +
  geom_point(size = 0.1, stroke = 0) + 
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_ratios_per_total_Supplementary_Figure_1A.png", 
         width = 2, height = 1.5, dpi = 600)

ggplot(gRNA_proportions, aes(x = ratio1, y= ratio2, color = log(total))) +
  geom_point(size = 0.1, stroke = 0) + 
  viridis::scale_color_viridis() +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_proportions_per_cell_Supplementary_Figure_1B.png", 
         width = 2, height = 1.5, dpi = 600)
  
ggplot(gRNA_proportions, aes(x = log10(top_to_second_best_ratio))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_top_to_second_best_ratio_Supplementary_Figure_1C.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions %>% dplyr::select(new_cell,ratio1) %>% distinct(), aes(x = ratio1)) +
  geom_density() +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_ratio1.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions, aes(x = ratio)) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
    monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_ratio1.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions, aes(x = log10(total))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_by_total_read_counts.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions, aes(x = log10(maxCount))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_by_maxCount.png", 
         width = 2, height = 2, dpi = 600)


ggplot(gRNA_proportions[gRNA_proportions$ratio > 0,], aes(x = ratio, y= log10(maxCount))) +
  geom_point(size = 0.01, stroke = 0) + 
  geom_hline(yintercept = log10(10), color = "dimgrey", linetype = "dashed") +
  geom_vline(xintercept = 0.75, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dashed") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPR_MMR-HPRT1_gRNA_ratios_per_maxCount.png", 
         width = 2, height = 1.5, dpi = 600)

# Filter for sgRNAs that have 3 or less sgRNAs (1 sgRNA should be the most porevalent given thte low MOI)
gRNA_proportions %>% 
  filter(maxCount >= 10 & ratio > 0.3) %>% 
  pull(new_cell) %>% 
  unique() %>% 
  length()

gRNA_proportions <- left_join(gRNA_proportions,
                              CRISPR_library,
                              by = "gRNA")

gRNA_proportions_assignment_subset <- gRNA_proportions %>%
filter(maxCount >= 10 & ratio >= 0.3) %>%
dplyr::select(new_cell, gRNA,gRNA_id, gene_id, CRISPR_gRNA, rank)


gRNA_proportions_assignment_subset_wide <- gRNA_proportions_assignment_subset %>%
group_by(new_cell) %>%
mutate(protospacer_sequence = paste0(gRNA, collapse = ","),
  gRNA_id = paste0(gRNA_id, collapse = ","),
  gene_id = paste0(gene_id, collapse = ","),
  CRISPR_gRNA = paste0(CRISPR_gRNA, collapse = ","),
  guide_number = sum(rank)) %>%
mutate(guide_number = ifelse(guide_number == 3,2,ifelse(guide_number == 6, 3, 1))) %>%
dplyr::select(new_cell,protospacer_sequence,gRNA_id,gene_id,CRISPR_gRNA,guide_number) %>%
as.data.frame() %>%
group_by(new_cell) %>%
dplyr::slice(1:1) %>%
ungroup()

### Create a new object with update metadata
new_colData <- colData(cds) %>% as.data.frame()
new_colData <- left_join(new_colData,gRNA_proportions_assignment_subset_wide, by = "new_cell")
new_colData <- as(new_colData, "DataFrame")
row.names(new_colData) <- new_colData$cell

identical(row.names(new_colData), row.names(colData(cds)))
identical(row.names(new_colData), colnames(exprs(cds)))

colData(cds) <- new_colData

### Save preprocessed Monocle3 cds object

saveRDS(cds, "sciPlexGxE_1_preprocessed_cds.rds")



