library(dplyr)
library(ggplot2)
library(tidyr)
library(furrr)
library(tictoc)
library(monocle3)})

### Source files and functions including those from our previous sci-plex repository 
### (https://github.com/cole-trapnell-lab/sci-plex)
source("sci-plex/bin/cell_cycle.R")
cc.genes = readRDS("sci-plex/bin/cc.genes.RDS")

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size = 1000e7)
options(stringsAsFactors = FALSE)

### Load in data from notebook 1 cleanup
cds <- readRDS("sgRNA_updated_cds.rds")

cds <- estimate_size_factors(cds)
cds <- estimate_cell_cycle(cds, 
                           g1s_markers = cc.genes$s.genes, 
                           g2m_markers = cc.genes$g2m.genes)

### Filter doublets using scrublet
colData(cds) %>% 
as.data.frame() %>%
ggplot() +
geom_density(aes(x = log10(scrublet_score))) +
geom_vline(xintercept = log10(0.25), linetype = "dashed") +
monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("QC_plots/Scrublet_scores_distribution.png", 
         width = 2, height = 2, dpi = 600)

colData(cds) %>% 
as.data.frame() %>%
mutate(scrublet_cutoff = ifelse(scrublet_score < 0.25, "singlet","multiplet")) %>%
ggplot() +
geom_boxplot(aes(x = factor(scrublet_cutoff,levels = c("singlet","multiplet")), y = log10(n.umi))) +
xlab("Scrublet prediction") +
monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("QC_plots/Scrublet_scores_cutoff_n.umi.png", 
         width = 2, height = 2, dpi = 600)

### Number of cells after doublet removal = 1042721 (0.8804325% inferred doublet rate)
cds <- cds[,colData(cds)$scrublet_score < 0.25]

ggplot(colData(cds) %>% as.data.frame(), aes(x = log10(top_to_second_best_ratio))) +
  geom_density() +
  geom_vline(xintercept = log10(2.5), linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("QC_plots/hash_top_to_second_best_ratio.png", 
         width = 2, height = 2, dpi = 600)

ggplot(colData(cds) %>% as.data.frame(), aes(x = log10(hash_umis))) +
  geom_density() +
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "dimgrey") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) +
  ggsave("QC_plots/hash_umis.png", 
         width = 2, height = 2, dpi = 600)

### Filter doublets using hashes. Number of cells after hash filtering 991940 (4.870047% empirical doublet rate)
cds <- cds[,colData(cds)$top_to_second_best_ratio >= 2.5 &
    colData(cds)$hash_umis >= 5]

### Given the size of the datatset. There are a small number of cells missassigned by the hashes across cell types.
### Below we use dimensionality reduction to filter for these instances, which although rare may have an outsize
### effect on downstream analyses.
genes_expressed_per_cell_line = colData(cds) %>%
                                as.data.frame() %>%
                                group_by(cell_line) %>%
                                nest() %>%
                                mutate(fraction_genes = purrr::map(data, .f = function(pdata_subset, cds) {
                                    cds_subset = cds[,as.character(pdata_subset$cell)]
                                    cds_subset = detect_genes(cds_subset)
                                    tibble(id = rowData(cds_subset)$id,
                                           gene_short_name = rowData(cds_subset)$gene_short_name,
                                           num_cells_expressed = rowData(cds_subset)$num_cells_expressed,
                                           fraction_cells_expressed = rowData(cds_subset)$num_cells_expressed / ncol(cds_subset))
                                }, cds))

genes_expressed_per_cell_line = genes_expressed_per_cell_line %>% 
  unnest(fraction_genes) %>%
  dplyr::select(everything(),-data)

expressed_genes = genes_expressed_per_cell_line %>% 
  filter(fraction_cells_expressed > 0.1) %>% 
  ungroup() %>% 
  dplyr::select(id) %>%
  distinct()

expressed_genes = as.character(expressed_genes$id)

length(expressed_genes)

cds <- preprocess_cds(cds,
                      method = "PCA",
                      num_dim = 25,
                      norm_method = "log",
                      use_genes = expressed_genes,
                      verbose =  TRUE)

cds <- reduce_dimension(cds,
                        max_components = 2,
                        reduction_method = "UMAP",
                        umap.metric = "cosine",
                        umap.n_neighbors = 20,
                        umap.min_dist = 0.1,
                        umap.fast_sgd=FALSE, 
                        cores=1,
                        verbose = T)

colData(cds)$UMAP1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 <- reducedDims(cds)[["UMAP"]][,2]

cds <- cluster_cells(cds, resolution = 1e-6)

colData(cds)$Cluster <- clusters(cds, reduction_method = "UMAP")
colData(cds)$Partition <- partitions(cds, reduction_method = "UMAP")

### Inspect the relationship between cluster and hash cell assignment
dir.create("UMAPs")

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = Cluster), 
             size = 0.01, stroke = 0) +
  facet_wrap(~cell_line)  +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.4,"line"))
  ggsave("UMAPs/UMAP_by_cluster_and_cell_line.png", dpi = 600, height = 1, width = 1.75)

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = Partition), 
             size = 0.01, stroke = 0) +
  facet_wrap(~cell_line)  +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "right",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
  ggsave("UMAPs/UMAP_by_partition_and_cell_line.png", dpi = 600, height = 1, width = 2.5)

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = treatment), 
             size = 0.01, stroke = 0) +
  facet_wrap(~cell_line)  +
  scale_color_manual("Treatment", values = c("lapatinib" = "navy","nintedanib" = "darkolivegreen",
                                "trametinib" = "brown4", "zstk474" = "darkcyan",
                                "vehicle" = "grey70")) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "right",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.4,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
  ggsave("UMAPs/UMAP_by_treatment_and_cell_line.png", dpi = 600, height = 1, width = 2.5)

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = proliferation_index), 
             size = 0.01, stroke = 0) +
  facet_wrap(~cell_line)  +
  scale_color_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "right",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.4,"line"))
  ggsave("UMAPs/UMAP_by_proliferation_index_and_cell_line.png", dpi = 600, height = 1, width = 2.5)

colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = log10(top_to_second_best_ratio)), 
             size = 0.01, stroke = 0) +
  facet_wrap(~cell_line)  +
  scale_color_viridis("Hash\npurity") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "right",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.4,"line"))
  ggsave("UMAPs/UMAP_by_hash_purity_and_cell_line.png", dpi = 600, height = 1, width = 2.5)


colData(cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = Cluster, y  = log10(top_to_second_best_ratio)), 
             size = 0.01, stroke = 0, outlier.shape = NA) +
  facet_wrap(~cell_line)  +
  #scale_color_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
    legend.position = "none",
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.4,"line"))
    ggsave("UMAPs/Cluster_hash_purity_and_cell_line.png", dpi = 600, height = 1, width = 2.5)

### Note: If you rerun this pipeline the cluster assignments below may change depending on your environment.
### You can identify the filtered cells from the cds object in GEO (GSM7056149_sciPlexGxE_2_preprocessed_cds.list.RDS) 
### which are already prefiltered by this code.
### Keep cells were hash assignment and cluster/partition distribution agrees across the majority of cells
A172_cell_list <- row.names(colData(cds)[colData(cds)$cell_line == "A172" &
                              colData(cds)$Partition %in% c("1","5"),])
# saveRDS(A172_cell_list,"A172_cell_list.RDS")

T98G_cell_list <- row.names(colData(cds)[colData(cds)$cell_line == "T98G" &
                              colData(cds)$Partition %in% c("2"),])
# saveRDS(T98G_cell_list,"T98G_cell_list.RDS")

U87MG_cell_list <- row.names(colData(cds)[colData(cds)$cell_line == "U87MG" &
                              colData(cds)$Partition %in% c("3","4"),])
# saveRDS(U87MG_cell_list,"U87MG_cell_list.RDS")

A172_cell_list <- readRDS("A172_cell_list.RDS")
T98G_cell_list <- readRDS("T98G_cell_list.RDS")
U87MG_cell_list <- readRDS("U87MG_cell_list.RDS")

length(A172_cell_list)
length(T98G_cell_list)
length(U87MG_cell_list)

### Divide cds objects by cell line. This list of cds objects is the main file used for downstream analyses.
cds.list <- list()

A172_cell_list_updated <- intersect(A172_cell_list,row.names(colData(cds)))
T98G_cell_list_updated <- intersect(T98G_cell_list,row.names(colData(cds)))
U87MG_cell_list_updated <- intersect(U87MG_cell_list,row.names(colData(cds)))

length(A172_cell_list_updated)
length(T98G_cell_list_updated)
length(U87MG_cell_list_updated)

cds.list[["A172"]] <- cds[,A172_cell_list_updated]
cds.list[["T98G"]] <- cds[,T98G_cell_list_updated]
cds.list[["U87MG"]] <- cds[,U87MG_cell_list_updated]

### Save final object. Corresponds to GSM7056149_sciPlexGxE_2_preprocessed_cds.list.RDS in the GEO repository.
# saveRDS(cds.list, "cds.list.RDS")










