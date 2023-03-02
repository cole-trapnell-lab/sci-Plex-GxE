library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(devtools)
library(monocle3)

### Source files and functions from our previous sci-plex repository (https://github.com/cole-trapnell-lab/sci-plex)
source("path_to_repo/sci-plex/bin/cell_cycle.R")
source("path_to_repo/sci-plex/sci-plex/bin/viability.R")
cc.genes = readRDS("path_to_repo/sci-plex/bin/cc.genes.RDS")

### Load expression data
cds <- readRDS("path_to_sciPlexGxE_1_preprocessed_cds.rds")
dim(cds)
dim(cds[,!is.na(colData(cds)$top_to_second_best_ratio_W)])

### Filter dataset by the confidence of hash assignment
cds <- cds[,colData(cds)$top_to_second_best_ratio_W >= 2.5  &
             colData(cds)$hash_umis_W >= 5]

dim(cds)
median(colData(cds)$n.umi)

dim(cds[,!is.na(colData(cds)$CRISPR_gRNA)])

dim(cds[,colData(cds)$gRNA_library == "HPRT1" &
          colData(cds)$gene_id %in% c("HPRT1","NTC")])

### Filter dataset by experiment (HPRT1) and for cells with 1 expressed sgRNA
HPRT1_cds <- cds[,colData(cds)$gRNA_library == "HPRT1" &
                  colData(cds)$guide_number == 1 &
                  colData(cds)$gene_id %in% c("HPRT1","NTC")]

HPRT1_cds <- estimate_size_factors(HPRT1_cds)

HPRT1_cds <- estimate_cell_cycle(HPRT1_cds, 
                           g1s_markers = cc.genes$s.genes, 
                           g2m_markers = cc.genes$g2m.genes)

### Inspect how often we miassign sgRNAs to cells known to expressed CRISPRi and CRISPRa sgRNAs (i.e., hashed separately)
ggplot(colData(HPRT1_cds) %>%as.data.frame(),
       aes(x = CRISPR_hash, fill = CRISPR_gRNA)) +
  geom_bar(color = "black",  size =  0.1) +
  facet_wrap(~factor(CRISPR_hash, levels = c("CRISPRi","CRISPRa")), scales = "free") +
  monocle_theme_opts() +
  scale_fill_manual("gRNA\nCRISPR\nassignment",
                    values = c("control" = "grey70", "CRISPRi" = "navy", "CRISPRa" = "brown4"),
                    labels = c("control" = "NTC", "CRISPRi" = "CRISPRi", "CRISPRa" = "CRISPRa")) +
  xlab("Hash CRISPR assignment") +
  ylab("Count") +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.4,"line"), 
        legend.key.height = unit(0.4,"line"),
        strip.text = element_blank()) +
  ggsave("gRNA_QC_plots/CRISPR_barnyard_results_Figure_1B.png", 
         dpi = 600, height = 1, width = 2)

### Inspect the expression of HPRT1 by CRISPR type 
HPRT1_plot_subset <- HPRT1_cds[rowData(HPRT1_cds)$gene_short_name == "HPRT1",]

colData(HPRT1_plot_subset) %>% 
  as.data.frame() %>%
  filter(dose == 0) %>% 
  group_by(CRISPR_hash, gene_id) %>% 
  summarise(n = n())

colData(HPRT1_plot_subset)$plot_label <- ifelse(colData(HPRT1_plot_subset)$gene_id == "HPRT1",
                                                paste0(colData(HPRT1_plot_subset)$gene, "_",
                                                colData(HPRT1_plot_subset)$CRISPR_gRNA),
                                                colData(HPRT1_plot_subset)$gene_id)

colData(HPRT1_plot_subset)$plot_label <- factor(colData(HPRT1_plot_subset)$plot_label,
                                                levels = c("NTC",
                                                           "HPRT1_CRISPRi",
                                                           "HPRT1_CRISPRa"))  
# plot_genes_violin(HPRT1_plot_subset,
#                   group_cells_by = "plot_label") +
#   scale_fill_manual(values = c("NTC" = "dimgrey",
#                                "HPRT1_CRISPRi" = "navy",
#                                "HPRT1_CRISPRa" = "brown4"),
#                     labels = c("NTC" = "dimgrey",
#                                "HPRT1_CRISPRi" = "navy",
#                                "HPRT1_CRISPRa" = "brown4")) +
#   theme(text = element_text(size = 6),
#         legend.key.width = unit(0.4,"line"), 
#         legend.key.height = unit(0.4,"line"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.x = element_blank()) +
#   ggsave("gRNA_QC_plots/HPRT1_levels_across_CRISPRs.png", 
#          dpi = 600, height = 1.5, width = 1)

colData(HPRT1_plot_subset)$gene_id <- factor(colData(HPRT1_plot_subset)$gene_id, levels = c("NTC","HPRT1"))

plot_percent_cells_positive(HPRT1_plot_subset[,colData(HPRT1_plot_subset)$CRISPR_hash == "CRISPRi"],
                            group_cells_by = "gene_id", panel_order = c("NTC","HPRT1")) +
  scale_fill_manual(values = c("NTC" = "grey70",
                               "HPRT1" = "navy")) +
  ggtitle("CRISPRi") +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.7,"line"), 
        legend.key.height = unit(0.7,"line"),
        title = element_text(hjust = 0, size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text.x = element_blank()) +
  ggsave("gRNA_QC_plots/HPRT1_levels_across_CRISPRi_Figuere_1C.png", 
         dpi = 600, height = 1.25, width = 0.8)

plot_percent_cells_positive(HPRT1_plot_subset[,colData(HPRT1_plot_subset)$CRISPR_hash == "CRISPRa"],
                            group_cells_by = "gene_id", panel_order = c("NTC","HPRT1")) +
  scale_fill_manual(values = c("NTC" = "grey70",
                               "HPRT1" = "brown4")) +
  ggtitle("CRISPRa") +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.7,"line"), 
        legend.key.height = unit(0.7,"line"),
        title = element_text(hjust = 0, size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text.x = element_blank()) +
  ggsave("gRNA_QC_plots/HPRT1_levels_across_CRISPRa_Figuere_1D.png", 
         dpi = 600, height = 1.25, width = 0.8)

### Inspect the expression of proliferation associated genes by genotype and dose of 6-thioguanine
dir.create("Heatmaps")

colData(HPRT1_cds) %>%
  as.data.frame() %>%
  filter(CRISPR_hash == "CRISPRi") %>%
  mutate(dose = factor(as.character(dose), levels = c("0","0.1","0.5","1",
                                                      "5","10","50","100"))) %>%
  group_by(CRISPR_hash,gene_id,dose) %>%
  mutate(mean_proliferation_index = mean(proliferation_index)) %>%
  dplyr::select(CRISPR_hash,gene_id,dose,mean_proliferation_index) %>%
  distinct() %>%
  ungroup() %>%
  ggplot() +
  geom_tile(aes(x = gene_id, y = dose, fill = mean_proliferation_index)) +
  viridis::scale_fill_viridis("Proliferation\nindex", option = "viridis") +
  ylab("[6-Thioguanine] (µM)") +
  xlab("Pertubation") +
  coord_flip() +
  monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.6,"line"), 
        legend.key.height = unit(0.5,"line")) +
  ggsave("Heatmaps/HPRT1_proliferation_heatmap_Figuer_1F.png",
         height = 1, width = 3.5, dpi = 600)



  
  
  

  