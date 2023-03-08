library(ggplot2)
library(ggridges)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(devtools)
library(monocle3)

### Function to update from monocle 2 to monocle 3 (should not have to do this anymore) ###
update_cds <- function(cds){
  pd <- Biobase::pData(cds)
  fd <- Biobase::fData(cds)
  exprs <- Biobase::exprs(cds)
  new_cds <- new_cell_data_set(expression_data = exprs, 
                               cell_metadata = pd, 
                               gene_metadata = fd)
}

### Read in cds object output by the processing pipeline post hash assignment
cds <- readRDS("path_to_cds_2.RDS")

cds <- update_cds(cds)

### Rename metadata columns
colnames(colData(cds)) <- c("cell","sample","Size_Factor","n.umi",
                            "hash_umis", "pval","qval","top_to_second_best_ratio",
                            "top_oligo","RT_well","Lig_well","P7_index","P5_index")

### Add a metadata column with a new cell name based only on the P5 indix to align
### the result from sequencing sgRNA enrichment libraries
colData(cds)$PCR_plate <- sapply(colData(cds)$cell,function(x){PCR_plate <- substring(x, 3, 3)})
colData(cds)$P5_cell <- sapply(colData(cds)$cell, function(x){P5_cell <- substring(x, 5)})
colData(cds)$new_cell <- paste0(colData(cds)$PCR_plate,"_",colData(cds)$P5_cell)

# saveRDS(cds, "updated_cds.RDS")
