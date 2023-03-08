library(ggplot2)
library(ggridges)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(devtools)
library(monocle3)

# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e7)

### Load updated cds object (notebook 0). Starting cell number at a 500 UMI cutoff -> 1052205
cds <- readRDS("updated_cds.RDS")

### Add background loading calculated using Jonathan's method to colData
background_projections_df <- readRDS("background_projections_df_JS_Packer_method.RDS")
background_projections_df$cell <- as.character(background_projections_df$cell)

identical(colData(cds)$cell,background_projections_df$cell)

colData(cds)$background_loading <- background_projections_df$background_projections
### Remove cells without a hash -> 1052117 (99.99164%)
cds <- cds[,!is.na(colData(cds)$top_oligo)]

### Cleaneup hash data and remove any cells with hashes mapped from a different replicate
### i.e., remove rare instances where a cell in replicate 1, RT wells 1 to 192, has a 
### hash assigned to replicate 2, RT wells 193 to 384, or vice-versa
rep1_sample_sheet <- read.table("GSM7056149_sciPlexGxE_2_hash_whitelist_rep1.txt", 
  sep = "\t", header = TRUE)
rep2_sample_sheet <- read.table("GSM7056149_sciPlexGxE_2_hash_whitelist_rep2.txt", 
  sep = "\t", header = TRUE)

rep1_colData <- colData(cds)[colData(cds)$RT_well %in%  c(1:192),]
rep1_colData$hash_oligo <- as.character(rep1_colData$top_oligo)

rep1_hash_collisions <- setdiff(rep1_colData$hash_oligo,rep1_sample_sheet$hash_oligo)
rep1_cell_collisions <- rep1_colData[rep1_colData$hash_oligo %in% rep1_hash_collisions,]$cell

rep1_colData <- rep1_colData[!(row.names(rep1_colData) %in% rep1_cell_collisions),]

rep2_colData <- colData(cds)[colData(cds)$RT_well %in%  c(193:384),]
rep2_colData$hash_oligo <- as.character(rep2_colData$top_oligo)

rep2_hash_collisions <- setdiff(rep2_colData$hash_oligo,rep2_sample_sheet$hash_oligo)
rep2_cell_collisions <- rep2_colData[rep2_colData$hash_oligo %in% rep2_hash_collisions,]$cell

rep2_colData <- rep2_colData[!(row.names(rep2_colData) %in% rep2_cell_collisions),]

rep1_colData <-merge(rep1_colData,rep1_sample_sheet, by = "hash_oligo")
rep2_colData <-merge(rep2_colData,rep2_sample_sheet, by = "hash_oligo")

hash_colData <- rbind(rep1_colData,rep2_colData)
row.names(hash_colData) <- hash_colData$cell
hash_colData$hash_oligo <- NULL

# Remove hash replicate collisions -> 1051983 (99.99789%)
cds <- cds[,!(row.names(colData(cds)) %in% c(rep1_cell_collisions,rep2_cell_collisions))]

### Replace metadata with metadata containig hash information
hash_colData <- hash_colData[row.names(colData(cds)),]
identical(row.names(hash_colData),row.names(colData(cds)))
identical(row.names(hash_colData),colnames(exprs(cds)))

colData(cds) <- hash_colData

### Plot some useful metrics per cell
dir.create("QC_plots")

cell_line_stats <- colData(cds) %>% 
  as.data.frame() %>%
  group_by(cell_line) %>%
  summarize(total_cells = n(), 
            median_umis = median(n.umi))

ggplot(cell_line_stats, aes(x = cell_line, y = total_cells)) +
  geom_bar(stat = "identity") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cell line") +
  ylab("Cell count") +
  ggsave("QC_plots/Total_cells_per_cell_line.png",
         width = 1, height = 1, dpi = 600)

colData(cds) %>% 
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = cell_line, y = log10(n.umi)), outlier.shape = NA, size = 0.1) +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cell line") +
  ylab("log(UMI count)") +
  ggsave("QC_plots/UMI_counts_per_cell_line.png",
         width = 1, height = 1, dpi = 600)

cell_count_df <- data.frame(subset = c("All", "Hash", "Hash > 10 umis"),
                            count = c(1052205,1052117,1047065))

ggplot(cell_count_df, aes(x = subset, y = count)) +
  geom_bar(stat = "identity") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cell subset") +
  ylab("Cell count") +
  ggsave("QC_plots/Total_cells_per_cell_subset.png",
         width = 1, height = 1, dpi = 600)

#### Append metadta with the results of sgRNA enrichment
kinome_gRNA_assignment <- read.table("GSM7056149_sciPlexGxE_2_gRNASampleSheet.txt", header = FALSE , sep = "\t")
colnames(kinome_gRNA_assignment) <- c("sgRNA", "protospacer_sequence")
kinome_gRNA_assignment$sgRNA <- as.character(kinome_gRNA_assignment$sgRNA)
kinome_gRNA_assignment$protospacer_sequence <- as.character(kinome_gRNA_assignment$protospacer_sequence)
kinome_gRNA_assignment$gene <- sapply(kinome_gRNA_assignment$sgRNA, 
                                function(x) stringr::str_split(x, pattern = "_")[[1]][1])

### Random region control sgRNAs have an extra base in the whitelist, this removes it
new_protospacer_sequence <- sapply(row.names(kinome_gRNA_assignment), function(x){

  ifelse(kinome_gRNA_assignment[x,]$gene %in% c("random","originalIDTorder"),
         substring(kinome_gRNA_assignment[x,]$protospacer_sequence,1,19),
         kinome_gRNA_assignment[x,]$protospacer_sequence)
  
  })

kinome_gRNA_assignment$protospacer_sequence <- new_protospacer_sequence

### Cleanup sgRNA metadata for consitent naming of control sgRNAs
kinome_gRNA_assignment[kinome_gRNA_assignment$sgRNA == "originalIDTorder_1",]$gene <- "random"
kinome_gRNA_assignment[kinome_gRNA_assignment$sgRNA == "originalIDTorder_1",]$sgRNA <- "random_region_35"

### Import kinome metadata
kinome_metadata <- read.table("kinome_metadata.txt", header = TRUE)
colnames(kinome_metadata) <- c("gene","protein","kinase_family")
kinome_metadata$gene <- as.character(kinome_metadata$gene)
kinome_metadata$protein <- as.character(kinome_metadata$protein)
kinome_metadata$kinase_family <- as.character(kinome_metadata$kinase_family)
kinome_metadata <- rbind(kinome_metadata, c("NTC","NTC","control"))
kinome_metadata <- rbind(kinome_metadata, c("random","random","control"))

kinome_gRNA_assignment <- left_join(kinome_gRNA_assignment, kinome_metadata, by = "gene")
colnames(kinome_gRNA_assignment) <- c("sgRNA","gRNA","gene","protein","kinase_family")
  
kinome_gRNA_assignment$gRNA <- toupper(kinome_gRNA_assignment$gRNA)

gRNA_enrichment_results <- data.table::fread("GSM7056149_sciPlexGxE_2_gRNATable_reads.out.txt")
gRNA_enrichment_results_gRNA_match <- gRNA_enrichment_results[gRNA_enrichment_results$gRNA %in% kinome_gRNA_assignment$gRNA,]

### Number of cells found in gRNA enrichment library -> 1050981 (99.88367%)
length(intersect(colData(cds)$new_cell, gRNA_enrichment_results$new_cell))

### Number of cells with an exact math to a gRNA in the whitelist -> 988276 (93.92428%)
length(intersect(colData(cds)$new_cell, gRNA_enrichment_results_gRNA_match$new_cell))

### Number of cells with at least 2 gRNA reads per cell -> 761062 (72.3302%)
length(intersect(colData(cds)$new_cell, gRNA_enrichment_results_gRNA_match[gRNA_enrichment_results_gRNA_match$read_count >= 2,]$new_cell))

### Calculate per cell enrichment of sgRNAs over background
gRNA_enrichment_results_subset <- gRNA_enrichment_results_gRNA_match[gRNA_enrichment_results_gRNA_match$new_cell %in% colData(cds)$new_cell,]
gRNA_enrichment_results_subset$duplication_rate <- 1-(gRNA_enrichment_results_subset$umi_count/gRNA_enrichment_results_subset$read_count)
median(gRNA_enrichment_results_subset$duplication_rate)
median(gRNA_enrichment_results_subset[gRNA_enrichment_results_subset$read_count > 1,]$duplication_rate)

gRNA_enrichment_results_summary <- gRNA_enrichment_results_subset %>% 
  group_by(new_cell,gRNA) %>% 
  mutate(total_read_count = sum(read_count), total_umi_count = sum(umi_count))  %>%
  dplyr::select(new_cell,gRNA, total_read_count, total_umi_count) %>% 
  distinct() 

gRNA_proportions <- gRNA_enrichment_results_summary %>% ungroup() %>%
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

length(intersect(colData(cds)$new_cell, gRNA_proportions$new_cell))

### Save intermediate gRNA proportions results
# saveRDS(gRNA_proportions, "gRNA_proportions.RDS")

### Merge gRNA proportions to kinome sgRNA whitelist
gRNA_proportions_assignment_subset <- left_join(gRNA_proportions,
                                                kinome_gRNA_assignment,
                                                by = "gRNA")

### Plot some useful QC metrics
dir.create("gRNA_QC_plots")

ggplot(gRNA_proportions, aes(x = ratio1, y= ratio2, color = log(total))) +
  geom_point(size = 0.1, stroke = 0) + 
  #geom_density2d(size = .1, color = "dimgrey") +
  scale_color_viridis() +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_proportions_per_cell.png", 
         width = 2, height = 1.5, dpi = 600)

ggplot(gRNA_proportions, aes(x = log10(top_to_second_best_ratio))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_top_to_second_best_ratio.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions, aes(x = log10(total))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_by_total_read_counts.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions, aes(x = log10(maxCount))) +
  geom_density() +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "dimgrey") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_by_maxCount.png", 
         width = 2, height = 2, dpi = 600)

ggplot(gRNA_proportions[gRNA_proportions$ratio > 0,], aes(x = ratio, y= log10(total))) +
  geom_point(size = 0.1, stroke = 0) + 
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_ratios_per_total.png", 
         width = 2, height = 1.5, dpi = 600)

ggplot(gRNA_proportions[gRNA_proportions$ratio > 0,], aes(x = ratio, y= log10(maxCount))) +
  geom_point(size = 0.01, stroke = 0) + 
  geom_hline(yintercept = log10(10), color = "dimgrey", linetype = "dashed") +
  geom_vline(xintercept = 0.75, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dashed") +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.1, "in"),
        legend.key.height = unit(0.1, "in")) +
  ggsave("gRNA_QC_plots/CRISPRi_kinome_gRNA_ratios_per_maxCount.png", 
         width = 2, height = 1.5, dpi = 600)

### Given the low MOI (~0.1) we do not expect to have many cells with multiple sgRNAs. Filter cells for those with 3 or less.
gRNA_proportions_assignment_subset_wide <- gRNA_proportions_assignment_subset %>%
filter(ratio > 0.3 | (top_to_second_best_ratio >= 2  & rank == "1") | (second_to_third_best_ratio >= 2  & rank == "1"))  %>%
dplyr::select(new_cell,gRNA,ratio,maxCount,rank,sgRNA,gene,protein,kinase_family)

gRNA_proportions_assignment_subset_wide <- gRNA_proportions_assignment_subset_wide %>%
group_by(new_cell) %>%
mutate(protospacer_sequence = paste0(gRNA, collapse = ","),
  gRNA_id = paste0(sgRNA, collapse = ","),
  gene_id = paste0(gene, collapse = ","),
  guide_number = sum(rank),
  gRNA_maxCount = max(maxCount),
  gRNA_topRatio = max(ratio)) %>%
mutate(guide_number = ifelse(guide_number == 3,2,ifelse(guide_number == 6, 3, 1))) %>%
dplyr::select(new_cell,protospacer_sequence,gRNA_id,gene_id,
  guide_number,gRNA_maxCount,gRNA_topRatio) %>%
as.data.frame() %>%
group_by(new_cell) %>%
dplyr::slice(1:1) %>%
ungroup()

### Update metadata to contain sgRNA information
new_colData <- colData(cds) %>% as.data.frame()
new_colData <- left_join(new_colData,gRNA_proportions_assignment_subset_wide, by = "new_cell")
new_colData <- as(new_colData, "DataFrame")
row.names(new_colData) <- new_colData$cell

identical(row.names(colData(cds)), colnames(exprs(cds)))
identical(row.names(new_colData), row.names(colData(cds)))

colData(cds) <- new_colData
 
# saveRDS(cds, "sgRNA_updated_cds.rds")


