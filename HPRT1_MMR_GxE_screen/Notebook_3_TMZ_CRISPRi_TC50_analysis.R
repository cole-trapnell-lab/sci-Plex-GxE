library(devtools)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(piano)
library(UpSetR)
library(parallel)
library(monocle3)

### Source files and functions including those from our previous sci-plex repository
### (https://github.com/cole-trapnell-lab/sci-plex)
source("perturbation_dose_response.R")
source("sci-plex/bin/cell_cycle.R")
source("sci-plex/bin/GSA_helper_functions.R")
cc.genes = readRDS("sci-plex/bin/cc.genes.RDS")

### Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
options(stringsAsFactors = FALSE)

#### Functions ####
### Define Helper function for downstream Gene Set Enrichment Analysis
replace_gene_names_vec <- function(input_vec, name_vec, retain_inds = c(-1,-2)) {
  temp <- merge(name_vec, input_vec, by="row.names")
  temp2 <- temp[,retain_inds]
  names(temp2) <- temp[,2]
  return(temp2)
}

loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml")) 
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    if (addUserInfo == "skip") 
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con, 
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip") 
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip") 
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip") 
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }

  addInfo <- as.data.frame(addInfo)

  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2) 
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
                              ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}

collect_gsa_hyper_results_clusters <- function (genes_list, clusters, gsc) 
{
  gene_universe <- unique(as.character(genes_list))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (cluster in cluster_ids) {
    cluster_genes <- unique(names(clusters[clusters == cluster]))
    gsaRes <- runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe, 
                          adjMethod = "BH")
    gsa_results[[cluster]] <- gsaRes
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
            theme_classic(base_size = 8)
            ggsave(plot_title, width = width, height = height)   
        
    }
    
}

### Calculate aggregate expression scores
calculate_signature_score <- function(cds, signature_genes){
  cds_subset = cds[rowData(cds)$id %in% signature_genes,] 
  aggregate_signature_expression = exprs(cds_subset)
  aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  signature_score = log(aggregate_signature_expression+1)
  return(signature_score)
}

### Define functions for calculating angular distance -----------------------

bin.ang.distance = function(cell1, cell2, count.mat) {
  
  cell_1 = 
    count.mat[,cell1] %>% 
    unit_vector()
  cell_2 = 
    count.mat[,cell2] %>% 
    unit_vector()
  
  ang.dist = angular_distance(cell_1,cell_2 )
  
  return(ang.dist)
}

# Calculate the norm of a vector
norm_vec  =
  function(x) {
    sqrt(sum(x^2))
  }

# Compute the angular distance between two cells
angular_distance = 
  function(cell_1, cell_2){
    if(identical(cell_1,cell_2)){
      0
    }
    else{
      (2/pi * acos(cell_1 %*% cell_2 /(norm_vec(cell_1) * norm_vec(cell_2))))^2
    }
  }

# Re-scale so that the magnitude of the vector is unit
unit_vector = 
  function(vector){
    vector / norm_vec(vector)
  }

# Take a normalized expression count matrix
# Return the unit mean vector for that matrix
unit_mean =
  function(count_matrix){
    Matrix::rowSums(count_matrix) %>%
      unit_vector()
  }

# Calculate JSD
kld = function(prob.1, prob.2) {
  if (any(!(prob.2 > 0)))
    warning("Zero values in prob.2")
  
  tmp = ifelse(prob.1 > 0, log2(prob.1/prob.2), 0)
  return(sum(prob.1 * tmp))
}

bin.js.distance = function(cell1, cell2, count.mat) {
  count.mat.1 = count.mat[, cell1]
  count.mat.2 = count.mat[, cell2]
  
  prob.1 = count.mat.1 / sum(count.mat.1)
  prob.2 = count.mat.2 / sum(count.mat.2)
  
  genes.1 = which(prob.1 > 0)
  genes.2 = which(prob.2 > 0)
  common.genes = unique(c(genes.1, genes.2))
  
  prob.1 = prob.1[common.genes]
  prob.2 = prob.2[common.genes]
  center = (prob.1 + prob.2) / 2.0
  
  js.div = 0.5 * kld(prob.1, center) + 0.5 * kld(prob.2, center)
  return(sqrt(js.div))
}

# Calculate the pairwise angular transcriptome distance of every cell to the average profile 
# of NTC cells exposed highest dose of drug

calculate_pairwise_perturbed_angular_distances <- function(cds, feature_genes, from_proportions = FALSE, cores){
  
  treatment.list <- list()
  
  treatments = as.character(unique(colData(cds)$treatment))
  colData(cds)$dose_character <- as.character(colData(cds)$dose)
  
  cds <- cds[feature_genes,]
  colData(cds) <- colData(cds)[,c("cell","Size_Factor", "dose","treatment","gene_id")]
  NTC_cds <- cds[,colData(cds)$dose == "100" & colData(cds)$gene_id == "NTC"]
  NTC_mean_matrix <- monocle3:::normalize_expr_data(NTC_cds, norm_method = "size_only")
  NTC_mean_matrix <- data.frame(row.names = names(rowMeans(NTC_mean_matrix)),NTC_mean =  rowMeans(NTC_mean_matrix)) %>%
    as.matrix() %>%
    as("dgCMatrix")

  cds_exprs <-  monocle3:::normalize_expr_data(cds, norm_method = "size_only")
  message("Features of mean NTC cells agree with dataset ", identical(row.names(cds_exprs), row.names(NTC_mean_matrix)))

  cds_exprs <- cbind(cds_exprs,NTC_mean_matrix)

  cds_colData <- colData(cds)
  NTC_colData <- data.frame(matrix(ncol = ncol(cds_colData), nrow = 1))
  colnames(NTC_colData) <- colnames(cds_colData)
  row.names(NTC_colData) <- c("NTC_mean")
  NTC_colData$cell <- row.names(NTC_colData)
  NTC_colData$Size_Factor <- 1
  NTC_colData$dose <- 100 ## changed here was
  NTC_colData$treatment <- c("temozolomide")
  NTC_colData$gene_id <- c("NTC_mean")

  cds_colData <- rbind(cds_colData,NTC_colData)
  cds_colData$Size_Factor <- 1

  cds_rowData <- rowData(cds)

  cds <- monocle3::new_cell_data_set(expression_data = cds_exprs, cell_metadata = cds_colData, gene_metadata = cds_rowData)

  colData(cds)$dose_character <- as.character(colData(cds)$dose)
  NTC_cells <- "NTC_mean"

  print(length(NTC_cells))
  
  colData(cds)$dose_character <- as.character(colData(cds)$dose)
  
  for (drug in treatments){
    
    dose.list <- list()
    
    this.treatment = drug
    
    for(dose in sort(as.character(unique(colData(cds)$dose_character)))){
      
      message(paste0("Determining the pairwise angular distance  of ",dose, " µM ",drug," exposed cells to the mean expression profile of ",unique(NTC_colData$dose)," µM  exposed NTC cells"))
      
      this.dose = dose
      
      target.list = list()
      
      for(target in sort(unique(colData(cds)$gene_id))){
        
        this.target = target
        
        cells_of_interest =
          colData(cds) %>%
          as.data.frame() %>%
          dplyr::mutate(dose_character = as.character(dose_character)) %>%
          filter(dose_character == this.dose,
                 gene_id == this.target,
                 treatment == this.treatment) %>%
          pull(cell) %>%
          as.character()

        if(length(cells_of_interest) > 1){
          
          # Get a pairwise pairing of every set of cells 
          
          pairwise_cell_mat = expand.grid(cells_of_interest,NTC_cells)
          colnames(pairwise_cell_mat) = c("cell_1","cell_2")
          pairwise_cell_mat$cell_1 <- as.character(pairwise_cell_mat$cell_1)
          pairwise_cell_mat$cell_2 <- as.character(pairwise_cell_mat$cell_2)
          
          # Remove rows where the two cells in comparison are the same
          pairwise_cell_mat =
            pairwise_cell_mat %>%
            filter(cell_1 != cell_2)
          
          pairwise_cell_mat$treatment = this.treatment
          pairwise_cell_mat$dose = this.dose
          pairwise_cell_mat$gene_id = this.target
          pairwise_cell_mat$total_cells = dim(pairwise_cell_mat)[1]
          
          #### Create cds subset over feature genes (dose-dependent DEGs)
          temp_cds =
            cds[feature_genes,union(cells_of_interest,NTC_cells)]
          
          if(from_proportions == FALSE){
            temp_cds = estimate_size_factors(temp_cds)
            
            count_mat =
              temp_cds %>%
              monocle3:::normalize_expr_data(norm_method = "log",
                                             pseudo_count = 1)
          }
          
          else(count_mat = exprs(temp_cds))
          
          pairwise_cell_mat_subset = pairwise_cell_mat
          
          pairwise_cell_mat_subset$ang.distance =
            with(pairwise_cell_mat_subset, mcmapply(function(x, y) {
              bin.ang.distance(x, y, count_mat)
            }, cell_1, cell_2, mc.cores = cores))
          
          pairwise_cell_mat_subset$js.distance =
            with(pairwise_cell_mat_subset, mcmapply(function(x, y) {
              bin.js.distance(x, y, count_mat)
            }, cell_1, cell_2, mc.cores = cores))
          
          
          target.list[[this.target]] =
            pairwise_cell_mat_subset
          message("finished ", this.target)
        }
        
      }
      
      dose.list[[as.character(this.dose)]] = target.list
      dose.list[[dose]] <- do.call("rbind",dose.list[[dose]])
      
    }
    
    treatment.list[[drug]] <- do.call("rbind",dose.list)
    
  }
  
  treatment_df <- do.call("rbind",treatment.list)
  return(treatment_df)
}

### Load expression data
cds <- readRDS("GSM7056148_sciPlexGxE_1_preprocessed_cds.rds")

### Filter dataset by the confidence of hash assignment
cds <- cds[,colData(cds)$top_to_second_best_ratio_W >= 2.5  &
             colData(cds)$hash_umis_W >= 5]
cds <- cds[,!is.na(colData(cds)$gRNA_id)]

### Create some useful metadata columns
colData(cds)$cell_type <- rep("A172", dim(cds)[2])
colData(cds)$vehicle <- sapply(colData(cds)$dose,function(x){ifelse(x == 0,TRUE,FALSE)})
colData(cds)$hash_plate <- colData(cds)$gRNA_library
colData(cds)$technical_replicate <- colData(cds)$RT_well

### Filter dataset by experiment (MMR) and for cells with 1 expressed sgRNA
TMZ_cds <- cds[,colData(cds)$gRNA_library == "MMR" &
                 colData(cds)$guide_number == 1 &
                 colData(cds)$gene_id %in% 
                 c("NTC","MGMT","MSH2","MSH3","MSH6","PMS2","MLH1")]


colData(TMZ_cds)$treatment <- sapply(colData(TMZ_cds)$treatment,function(x){ifelse(x == "dmso","temozolomide",x)})
colData(TMZ_cds)$dose_character <- as.character(colData(TMZ_cds)$dose)

colData(TMZ_cds)$dose_character <- factor(colData(TMZ_cds)$dose,
                                          levels = c("0","0.1","0.5","1",
                                                     "5","10","50","100"))

TMZ_cds <- detect_genes(TMZ_cds)
TMZ_cds <- estimate_size_factors(TMZ_cds)
TMZ_cds <- estimate_cell_cycle(TMZ_cds, 
                               g1s_markers = cc.genes$s.genes, 
                               g2m_markers = cc.genes$g2m.genes)

expressed_genes <- row.names(fData(TMZ_cds)[Matrix::rowSums(exprs(TMZ_cds) > 0) > 
                                              dim(TMZ_cds)[2]*0.05 ,])

### Estimate kd efficiency across genotypes
untreated_subset_cds <- TMZ_cds[rowData(TMZ_cds)$gene_short_name %in% c("MGMT","MSH2","MSH3","MSH6","PMS2","MLH1"),
                                colData(TMZ_cds)$dose == 0]

colData(untreated_subset_cds)$gene_id <- factor(colData(untreated_subset_cds)$gene_id,
                                                levels = c("NTC","MGMT","MSH2","MSH3","MSH6","MLH1","PMS2"))

colData(untreated_subset_cds)$gRNA_id <- factor(colData(untreated_subset_cds)$gRNA_id,
                                                levels = c("NTC_1","NTC_2","NTC_3","NTC_4","NTC_5","NTC_6",
                                                           "MGMT_1","MGMT_2","MGMT_3",
                                                           "MSH2_1","MSH2_2","MSH2_3",
                                                           "MSH3_1","MSH3_2","MSH3_3",
                                                           "MSH6_1","MSH6_2","MSH6_3",
                                                           "MLH1_1","MLH1_2","MLH1_3",
                                                           "PMS2_1","PMS2_2","PMS2_3"))

targets <- c("MSH2","MSH3","MSH6","MLH1","PMS2","MGMT")

plot_subset.list <- list()
plot_target_subset.list <- list()

for(target in targets){
  
  plot_subset_cds <- untreated_subset_cds[rowData(untreated_subset_cds)$gene_short_name == target,
                                          colData(untreated_subset_cds)$gene_id %in% c("NTC",target)]
  colData(plot_subset_cds)$gene <- factor(colData(plot_subset_cds)$gene_id,
                                          levels = c("NTC",target))
  
  plot_subset.list[[target]] <- plot_percent_cells_positive(plot_subset_cds,
                                                            group_cells_by = "gRNA_id") +
    theme(text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 4),
          legend.position = "none") 
  
  plot_target_subset.list[[target]] <- plot_percent_cells_positive(plot_subset_cds,
                                                            group_cells_by = "gene_id") +
    theme(text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 4),
          legend.position = "none")
  
}

pdf("gRNA_QC_plots/MMR_array_Supplementary_Figure_2A.pdf", useDingbats = FALSE, width = 2.4, height = 3)
do.call("grid.arrange", c(plot_subset.list, ncol = 2))
dev.off()

pdf("gRNA_QC_plots/MMR_target_array_Supplementary_Figure_2B.pdf", useDingbats = FALSE, width = 1.5, height = 3)
do.call("grid.arrange", c(plot_target_subset.list, ncol = 2))
dev.off()

### Inspect the expression of genes associated with DNA damage response signaling and proliferation
marker <- "CDKN1A"
marker_cds <- TMZ_cds[rowData(cds)$gene_short_name == marker,]
marker_exprs <- Matrix::t(Matrix::t(exprs(marker_cds))/colData(marker_cds)$Size_Factor)
marker_exprs <- log(marker_exprs + 1)

colData(marker_cds)$marker_exprs <- marker_exprs

colData(marker_cds) %>%
  as.data.frame() %>%
  mutate(dose = factor(as.character(dose), levels = c("0","0.1","0.5","1",
                                                      "5","10","50","100"))) %>%
  group_by(gene_id,dose) %>%
  mutate(mean_expression = mean(marker_exprs)) %>%
  dplyr::select(gene_id,dose,mean_expression) %>%
  distinct() %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = c("MSH2","MSH6","MLH1","PMS2","MSH3","MGMT","NTC"))) %>% 
  ggplot() +
  geom_tile(aes(x = gene_id, y = dose, fill = mean_expression)) +
  #facet_wrap(~CRISPR, scales = "free") +
  viridis::scale_fill_viridis("CDKN1A\nexpression", option = "magma") +
  coord_flip() +
  ylab("[Temozolomide] (µM)") +
  xlab("Pertubation") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.y  = element_text(face = "italic"),
        legend.key.width = unit(0.4,"line"), 
        legend.key.height = unit(0.4,"line"))
ggsave("Heatmaps/CDKN1A_expression_heatmap_Figure_2B.png",
         height = 1, width = 2, dpi = 600)

colData(TMZ_cds) %>%
  as.data.frame() %>%
  mutate(dose = factor(as.character(dose), levels = c("0","0.1","0.5","1",
                                                      "5","10","50","100"))) %>%
  group_by(gene_id,dose) %>%
  mutate(mean_proliferation_index = mean(proliferation_index)) %>%
  dplyr::select(gene_id,dose,mean_proliferation_index) %>%
  distinct() %>%
  ungroup() %>%
  mutate(gene_id = factor(gene_id, levels = c("MSH2","MSH6","MLH1","PMS2","MSH3","MGMT","NTC"))) %>% 
  ggplot() +
  geom_tile(aes(x = gene_id, y = dose, fill = mean_proliferation_index)) +
  #facet_wrap(~CRISPR, scales = "free") +
  viridis::scale_fill_viridis("Proliferation\nindex", option = "viridis") +
  coord_flip() +
  ylab("[Temozolomide] (µM)") +
  xlab("Pertubation") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.y  = element_text(face = "italic"),
        legend.key.width = unit(0.4,"line"), 
        legend.key.height = unit(0.4,"line"))
ggsave("Heatmaps/MMR_proliferation_heatmap_Figure_2C.png",
         height = 1, width = 2, dpi = 600)

CDKN1A_cds_subset <- TMZ_cds[rowData(TMZ_cds)$gene_short_name == "CDKN1A",colData(TMZ_cds)$dose %in% c(100)]

colData(CDKN1A_cds_subset)$gene_id <- factor(colData(CDKN1A_cds_subset)$gene_id,
                                             levels = rev(c("MSH2","MSH6","MLH1","PMS2","MSH3","MGMT","NTC")))

plot_genes_violin(CDKN1A_cds_subset, min_expr = 0.1, group_cells_by = "gene_id") +
  theme(text = element_text(size = 6),
        strip.text.x = element_text(size = 9)) +
  scale_fill_manual(values = c("NTC" = "grey70","MGMT" = "grey70", "MSH3" = "grey70",
                               "MSH2" = "darkorange2","MSH6" = "darkorange2",
                               "MLH1" = "darkorange2", "PMS2" = "darkorange2")) +
  xlab("Perturbation") +
  ggsave("Marker_plots/CDKN1A_high_dose_violin.png",
         width = 2,
         height  = 1.5,
         dpi = 600)

# Identify differentially expressed genes as a function of genotype for every dose of drug tested
targets <- unique(colData(TMZ_cds)$gene_id)
targets <- targets[targets != "NTC"]

NTC_cells.list <- list()

for(dose in sort(unique(colData(TMZ_cds)$dose_character))){

  NTC_cells.list[[dose]] <- row.names(colData(TMZ_cds)[colData(TMZ_cds)$gene_id == "NTC" &
                                            colData(TMZ_cds)$dose_character == dose,])

  print(length(NTC_cells.list[[dose]]))

}

kd_diff_test.list <- list()

for(target in targets){

  message(target)
  kd_diff_test.list[[target]] <- list()
  cds_subset <- TMZ_cds[expressed_genes,colData(TMZ_cds)$gene_id %in% c("NTC",target)]

  for(dose in  sort(unique(colData(TMZ_cds)$dose_character))){
    message(dose)

    cds_dose_subset <- cds_subset[,NTC_cells.list[[dose]]]
    cds_dose_subset <- detect_genes(cds_dose_subset)
    genes_in_reference_set <- row.names(subset(rowData(cds_dose_subset),
                                               num_cells_expressed > ncol(cds_dose_subset)*0.05))

    dose_subset_cells <- row.names(colData(TMZ_cds)[colData(TMZ_cds)$gene_id == target &
                                                      colData(TMZ_cds)$dose_character == dose,])

    cds_dose_subset <- cds_subset[,dose_subset_cells]
    cds_dose_subset <- detect_genes(cds_dose_subset)
    genes_in_treated_set <- row.names(subset(rowData(cds_dose_subset),
                                             num_cells_expressed > ncol(cds_dose_subset)*0.05))


    cds_dose_subset <- cds_subset[union(genes_in_reference_set,genes_in_treated_set),
                                  union(NTC_cells.list[[dose]],dose_subset_cells)]

    colData(cds_dose_subset)$gene_id <- factor(colData(cds_dose_subset)$gene_id,
                                               levels = c("NTC",target))

    diff_test <- fit_models(cds_dose_subset,
                            model_formula_str = "~gene_id")

    kd_diff_test.list[[target]][[as.character(dose)]] <- coefficient_table(diff_test) %>%
      dplyr::select(-model,-model_summary)

    kd_diff_test.list[[target]][[as.character(dose)]]$dose <- rep(as.character(dose),
                                                                  dim(kd_diff_test.list[[target]][[as.character(dose)]])[1])

  }

  kd_diff_test.list[[target]] <- do.call("rbind",kd_diff_test.list[[target]])
  kd_diff_test.list[[target]]$target <- rep(target, dim(kd_diff_test.list[[target]])[1])
  message("Done")
}

# Save result of DEG test
# saveRDS(kd_diff_test.list,"genotype_dependent_degs_tests_by_dose.list.rds")
# write.table(kd_diff_test_results, "MMR_perturbation_diff_test_results.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# kd_diff_test.list <- readRDS("genotype_dependent_degs_tests_by_dose.list.rds")
kd_diff_test_results <- do.call("rbind", kd_diff_test.list)
kd_diff_test_results <- kd_diff_test_results[!is.na(kd_diff_test_results$p_value),]

kd_diff_test_results$q_value <- p.adjust(kd_diff_test_results$p_value,
                                         method = "BH")

kd_diff_test_results_filtered <- kd_diff_test_results %>%
  filter(grepl("gene_id", term))

dir.create("DEG_testing")

ggplot(kd_diff_test_results_filtered, aes(x = normalized_effect,
                                          y = -log(q_value),
                                          color = q_value < 0.005)) +
  geom_point(size = 0.2, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.005), color = "dimgrey", linetype = "dotted",size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab("Effect size") +
  ylab("-log(p adj)") +
  scale_color_manual("Pass\nFilter\n(FDR < 0.5%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
  ggsave("DEG_testing/Effect_size_distribution.png", 
         dpi = 600, width = 1.5, height = 1)

ggplot(kd_diff_test_results_filtered, aes(x = normalized_effect,
                                          y = -log(q_value),
                                          color = q_value < 0.005)) +
  geom_point(size = 0.2, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.005), color = "dimgrey", linetype = "dotted",size = 0.1) +
  facet_wrap(~factor(target, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), ncol = 7) +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab(expression(paste("Normalized ",beta," coefficient"))) +
  ylab("-log(p adj)") +
  scale_color_manual("Pass\nFilter\n(FDR < 0.5%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("DEG_testing/Effect_size_distribution_faceted_Supplementary_Figure_2C.png", 
         dpi = 600, width = 4, height = 1.3)

ggplot(kd_diff_test_results_filtered, aes(x = normalized_effect,
                                          y = -log(q_value),
                                          color = target == "MGMT" & q_value < 0.005)) +
  geom_point(size = 0.2, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.005), color = "dimgrey", linetype = "dotted",size = 0.1) +
  geom_hline(yintercept = -log(0.0001), color = "brown4", size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab("Effect size") +
  ylab("-log(p adj)") +
  scale_color_manual("MGMT\nDEG\n(FDR < 0.5%)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("DEG_testing/Effect_size_distribution_MGMT.png", 
         dpi = 600, width = 1.5, height = 1)

ggplot(kd_diff_test_results_filtered, aes(x = p_value)) +
  geom_histogram(bins = 100) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab("p value") +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("DEG_testing/p_value_distribution.png", 
         dpi = 600, width = 1.5, height = 1)

ggplot(kd_diff_test_results_filtered, aes(x = q_value)) +
  geom_histogram(bins = 100) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab("q value") +
  guides(guides(colour = guide_legend(override.aes = list(size=2)))) +
  ggsave("DEG_testing/q_value_distribution.png", 
         dpi = 600, width = 1.5, height = 1)

kd_sig_genes <- kd_diff_test_results_filtered %>%
  filter(q_value < 0.005) %>%
  pull(id) %>% unique()

MGMT_sig_genes <- kd_diff_test_results_filtered %>%
  filter(target == "MGMT",q_value < 0.005) %>%
  pull(id) %>% unique()

kd_sig_genes <- unique(kd_sig_genes[!(kd_sig_genes %in% MGMT_sig_genes)])

sig_genes_df <- kd_diff_test_results_filtered %>%
  filter(q_value < 0.005) %>%
  dplyr::select(target,id)%>%
  group_by(target) %>% 
  distinct() %>%
  reshape2::dcast(formula = id ~ target, fun.aggregate = length)

png("DEG_testing/DEG_upSetR_plot_Supplementary_Figure_2D.png", width = 3, height = 1.75, units = "in", res = 1000)
UpSetR::upset(sig_genes_df,
      sets = c("MSH2","MSH6","MLH1","PMS2","MSH3","MGMT"),
      order.by = "freq",
      text.scale = 0.5,
      point.size = .5,
      line.size = 0.2)
dev.off()

vehicle_sig_genes_df <- kd_diff_test_results_filtered %>%
  filter(q_value < 0.005  & dose == "0") %>%
  dplyr::select(target,id)%>%
  group_by(target) %>% 
  distinct() %>% 
  reshape2::dcast(formula = id ~ target, fun.aggregate = length)

png("DEG_testing/Vehicle_DEG_upSetR_plot_Supplementary_Figure_2E.png", width = 3, height = 1.75, units = "in", res = 1000)
UpSetR::upset(vehicle_sig_genes_df,
      sets = c("MSH2","MSH6","MLH1","PMS2","MSH3","MGMT"),
      order.by = "freq",
      text.scale = 0.5,
      point.size = .5,
      line.size = 0.2)
dev.off()

### Inspect the correlation structure across genotypes, drug and doses in the experiment
sig_genes_matrix <- kd_diff_test_results_filtered %>%
  mutate(target_dose = paste0(target,"_",dose)) %>%
  filter(id %in% kd_sig_genes) %>%
  dplyr::select(id, target_dose, normalized_effect) %>%
  tidyr::spread(key = target_dose, value = normalized_effect)

sig_genes_matrix<- as.data.frame(sig_genes_matrix)
row.names(sig_genes_matrix) <- sig_genes_matrix$id
sig_genes_matrix$id <- NULL

sig_genes_cor_matrix <- cor(sig_genes_matrix,
                            use = "complete", 
                            method = "pearson")

hmcols <- colorRampPalette(c("white","lightpink2","firebrick2","firebrick4"))(35)

dir.create("Heatmaps")

correlation_ph <- pheatmap::pheatmap(sig_genes_cor_matrix,
                                     file = "Heatmaps/Correlation_heatmap_Supplementary_Figure_2F.png",
                                     color = hmcols,
                                     treeheight_row = 15,
                                     treeheight_col = 15,
                                     cutree_rows = 5,
                                     cutree_cols = 5,
                                     height = 3,
                                     width = 3,
                                     fontsize = 6,
                                     fontsize_col = 3,
                                     fontsize_row = 3) 

### Inspect correlation coefficients at the top doses of drug
MMR_pertubed_sig_genes_cor_matrix_subset <- sig_genes_cor_matrix[c("MSH2_50","MSH6_50","MLH1_50","PMS2_50",
                                                                   "MSH2_100","MSH6_100","MLH1_100","PMS2_100"),
                                                                 c("MSH2_50","MSH6_50","MLH1_50","PMS2_50",
                                                                   "MSH2_100","MSH6_100","MLH1_100","PMS2_100")]

max(MMR_pertubed_sig_genes_cor_matrix_subset[MMR_pertubed_sig_genes_cor_matrix_subset != 1])
min(MMR_pertubed_sig_genes_cor_matrix_subset)

### Inspect the distribution of cell states as a function of genotype and treatment
TMZ_cds <- preprocess_cds(TMZ_cds,
                          method = "PCA",
                          num_dim = 25,
                          norm_method = "log",
                          use_genes = kd_sig_genes)

TMZ_cds <- reduce_dimension(TMZ_cds,
                            max_components = 2,
                            preprocess_method = "PCA",
                            reduction_method = "UMAP",
                            umap.metric = "cosine",
                            umap.n_neighbors = 20,
                            umap.min_dist = 0.1,
                            umap.fast_sgd=FALSE,
                            cores=1,
                            verbose = T)

colData(TMZ_cds)$UMAP1 <- reducedDims(TMZ_cds)[["UMAP"]][,1]
colData(TMZ_cds)$UMAP2 <- reducedDims(TMZ_cds)[["UMAP"]][,2]

TMZ_cds <- cluster_cells(TMZ_cds, 
                         resolution = 5e-3,
                        reduction_method = "PCA")

colData(TMZ_cds)$Cluster <- clusters(TMZ_cds, reduction_method = "PCA")
length(unique(colData(TMZ_cds)$Cluster))
colData(TMZ_cds)$Partition <- partitions(TMZ_cds, reduction_method = "PCA")
length(unique(colData(TMZ_cds)$Partition))

dir.create("UMAPs")

colData(TMZ_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, fill = Cluster), size = 0.5, stroke = 0.01, shape = 21, color = "black") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/UMAP_by_cluster_Supplementary_Figure_3A.png", dpi = 600, height = 1, width = 1.5)

colData(TMZ_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, fill = log(dose + 0.01)), size = 0.5, stroke = 0.01, shape = 21, color = "black") +
  viridis::scale_fill_viridis("log(dose) µM", option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.4,"line"),
        legend.key.height = unit(0.4,"line"))
ggsave("UMAPs/UMAP_by_dose_Supplementary_Figure_3B.png", dpi = 600, height = 1, width = 1.5)

colData(TMZ_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, fill = proliferation_index), size = 0.5, stroke = 0.01, shape = 21, color = "black") +
  viridis::scale_fill_viridis("Proliferation\nindex") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.4,"line"),
        legend.key.height = unit(0.4,"line"))
ggsave("UMAPs/UMAP_by_proliferation_index_Supplementary_Figure_3C.png", dpi = 600, height = 1, width = 1.5)

colData(TMZ_cds) %>%
  as.data.frame() %>%
  mutate(gene_id = factor(gene_id, levels = c("NTC","MGMT","MSH3",
                                              "MSH2","MSH6","MLH1","PMS2"))) %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = gene_id, fill = gene_id), size = 0.5, stroke = 0.01, shape = 21) +
  scale_color_manual("Genotype",values = c("NTC" = "black","MGMT" = "black", "MSH3" = "black",
                                           "MSH2" = "black","MSH6" = "black",
                                           "MLH1" = "black", "PMS2" = "black")) +
  scale_fill_manual("Genotype",values = c("NTC" = "grey70","MGMT" = "dimgrey", "MSH3" = "dimgrey",
                                          "MSH2" = "firebrick3","MSH6" = "brown4",
                                          "MLH1" = "navy", "PMS2" = "darkcyan")) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/UMAP_by_genotype_Supplementary_Figure_3D.png", dpi = 600, height = 1, width = 1.5)

ggplot(colData(TMZ_cds) %>%
           as.data.frame() %>%
           mutate(gene_id = factor(gene_id, levels = c("NTC","MGMT","MSH3",
                                                       "MSH2","MSH6","MLH1","PMS2"))),
       aes(x = UMAP1, y  = UMAP2)) +
  geom_density2d() +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  facet_wrap(~factor(gene_id, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), ncol = 7) +
  viridis::scale_fill_viridis("Cell\ndensity",option = "magma") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.4,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/UMAP_by_genotype_faceted_Supplementary_Figure_3E.png", dpi = 600, height = 1, width = 6)

### Inspect the relationship between cell state proportions and treatment at the sgRNA and target gene levels
gRNA_cluster_counts <- reshape2::acast(
  colData(TMZ_cds) %>%
    as.data.frame() %>%
    dplyr::mutate(dummy = 1) %>%
    dplyr::select(gRNA_id, Cluster, dummy),
  gRNA_id ~ Cluster, value.var = "dummy",
  fun.aggregate = sum, fill = 0)

gRNA_cluster_proportions <- t(gRNA_cluster_counts/rowSums(gRNA_cluster_counts))

gRNA_cluster_proportions_cds <- new_cell_data_set(expression_data = gRNA_cluster_proportions,
                                                  cell_metadata = data.frame(row.names = colnames(gRNA_cluster_proportions),
                                                                             gRNA_id = colnames(gRNA_cluster_proportions),
                                                                             gene_id = sapply(colnames(gRNA_cluster_proportions),function(x)stringr::str_split(x,pattern = "_")[[1]][1])),
                                                  gene_metadata = data.frame(row.names = row.names(gRNA_cluster_proportions),
                                                                             id = row.names(gRNA_cluster_proportions),
                                                                             gene_short_name = row.names(gRNA_cluster_proportions)))


reducedDims(gRNA_cluster_proportions_cds)[["PCA"]] <- t(as.matrix(exprs(gRNA_cluster_proportions_cds)))

gRNA_cluster_proportions_cds <- reduce_dimension(gRNA_cluster_proportions_cds,
                                                 max_components = 2,
                                                 reduction_method = "UMAP",
                                                 umap.metric = "cosine",
                                                 umap.n_neighbors = 3,
                                                 umap.min_dist = 0.2,
                                                 umap.fast_sgd=FALSE,
                                                 cores=1,
                                                 verbose = T)

colData(gRNA_cluster_proportions_cds)$UMAP1 <- reducedDims(gRNA_cluster_proportions_cds)[["UMAP"]][,1]
colData(gRNA_cluster_proportions_cds)$UMAP2 <- reducedDims(gRNA_cluster_proportions_cds)[["UMAP"]][,2]

colData(gRNA_cluster_proportions_cds)$predicted_sentivity <- sapply(colData(gRNA_cluster_proportions_cds)$gene_id,
                                                                    function(x){ifelse(x %in% c("MSH2","MSH6","MLH1","PMS2"),
                                                                                       "Resistant",
                                                                                       "Sensitive")})
colData(gRNA_cluster_proportions_cds) %>%
  as.data.frame() %>%
  mutate(gene_id = factor(gene_id, levels = c("NTC","MGMT","MSH3",
                                              "MSH2","MSH6","MLH1","PMS2"))) %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = predicted_sentivity, fill = predicted_sentivity), size = 2, stroke = 0.1, shape = 21) +
  scale_fill_manual("Predicted\nsensitivity", values = c("Resistant" = "darkorange2", "Sensitive" =  "grey70")) +
  scale_color_manual("Predicted\nsensitivity", values = c("Resistant" = "black", "Sensitive" =  "black")) +
  ggrepel::geom_text_repel(data = colData(gRNA_cluster_proportions_cds) %>%
                             as.data.frame(),
                           aes(x = UMAP1, y  = UMAP2, label =  gRNA_id), 
                           size = 1.25,
                           force = 50,
                           segment.size = 0.1,
                           max.iter = 10000,
                           #label.size = 0.1,
                           box.padding = 0.05,
                           #label.padding = 0.1,
                           color = "black") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/UMAP_gRNA_cluster_proportions_Figure_2E.png", dpi = 600, height = 1, width = 1)

genotype_cluster_counts <- reshape2::acast(
  colData(TMZ_cds) %>%
    as.data.frame() %>%
    dplyr::mutate(dummy = 1) %>%
    dplyr::select(gene_id, Cluster, dummy),
  gene_id ~ Cluster, value.var = "dummy",
  fun.aggregate = sum, fill = 0)

genotype_cluster_proportions <- t(genotype_cluster_counts/rowSums(genotype_cluster_counts))

genotype_cluster_proportions_cds <- new_cell_data_set(expression_data = genotype_cluster_proportions,
                                                      cell_metadata = data.frame(row.names = colnames(genotype_cluster_proportions),
                                                                                 gene_id = colnames(genotype_cluster_proportions)),
                                                      gene_metadata = data.frame(row.names = row.names(genotype_cluster_proportions),
                                                                                 id = row.names(genotype_cluster_proportions),
                                                                                 gene_short_name = row.names(genotype_cluster_proportions)))


reducedDims(genotype_cluster_proportions_cds)[["PCA"]] <- t(as.matrix(exprs(genotype_cluster_proportions_cds)))

genotype_cluster_proportions_cds <- reduce_dimension(genotype_cluster_proportions_cds,
                                                     max_components = 2,
                                                     reduction_method = "UMAP",
                                                     umap.metric = "cosine",
                                                     umap.n_neighbors = 3,
                                                     umap.min_dist = 0.2,
                                                     umap.fast_sgd=FALSE,
                                                     cores=1,
                                                     verbose = T)

colData(genotype_cluster_proportions_cds)$UMAP1 <- reducedDims(genotype_cluster_proportions_cds)[["UMAP"]][,1]
colData(genotype_cluster_proportions_cds)$UMAP2 <- reducedDims(genotype_cluster_proportions_cds)[["UMAP"]][,2]

colData(genotype_cluster_proportions_cds)$predicted_sentivity <- sapply(colData(genotype_cluster_proportions_cds)$gene_id,
                                                                        function(x){ifelse(x %in% c("MSH2","MSH6","MLH1","PMS2"),
                                                                                           "Resistant",
                                                                                           "Sensitive")})

colData(genotype_cluster_proportions_cds) %>%
  as.data.frame() %>%
  mutate(gene_id = factor(gene_id, levels = c("NTC","MGMT","MSH3",
                                              "MSH2","MSH6","MLH1","PMS2"))) %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y  = UMAP2, color = predicted_sentivity, fill = predicted_sentivity), size = 2.5, stroke = 0.1, shape = 21) +
  scale_fill_manual("Predicted\nsensitivity", values = c("Resistant" = "darkorange2", "Sensitive" =  "grey70")) +
  scale_color_manual("Predicted\nsensitivity", values = c("Resistant" = "black", "Sensitive" =  "black")) +
  ggrepel::geom_label_repel(data = colData(genotype_cluster_proportions_cds) %>%
                              as.data.frame(),
                            aes(x = UMAP1, y  = UMAP2, label =  gene_id),
                            size = 1.25,
                            force = 50,
                            segment.size = 0.1,
                            label.size = 0.1,
                            box.padding = 0.1,
                            label.padding = 0.1,
                            color = "black") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.25,"line")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("UMAPs/UMAP_genotype_cluster_proportions_Figure_2F.png", dpi = 600, height = 1, width = 1)

### Define gene modules across the set of genotype dependent DE genes
gene_module_df <- find_gene_modules(TMZ_cds[kd_sig_genes,],
                                    resolution = 0.002,
                                    umap.fast_sgd = FALSE,
                                    random_seed = 2016L,
                                    cores = 1)

# saveRDS(gene_module_df, "gene_modules_across_genotype_dependent_degs_df.rds")
#gene_module_df <- readRDS("gene_modules_across_genotype_dependent_degs_df.rds")
#write.table(gene_module_df, "Gene_modules_across_MMR_pertrubed_DEGs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cell_group_df <- tibble::tibble(cell=row.names(colData(TMZ_cds[,colData(TMZ_cds)$dose %in% c(10,50,100)])), 
                                cell_group=colData(TMZ_cds[,colData(TMZ_cds)$dose %in% c(10,50,100)])$gene_id)
agg_mat <- aggregate_gene_expression(TMZ_cds[,colData(TMZ_cds)$dose %in% c(10,50,100)], gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))
colnames(agg_mat) <- c("MGMT","MLH1","MSH2","MSH3","MSH6","NTC","PMS2")

paletteLength <- 35
myBreaks <- c(seq(min(agg_mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(agg_mat)/paletteLength, max(agg_mat), length.out=floor(paletteLength/2)))

pheatmap(agg_mat, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         clustering_method="ward.D2",
         fontsize=6, 
         #color = viridis::magma(35),
         color = colorRampPalette(c("blue", "white", "red"))(paletteLength),
         breaks = myBreaks,
         width = 2.5,
         height = 2.5,
         treeheight_col = 15,
         treeheight_row = 15,
         fontsize_row = 6,
         fontsize_col = 6,
         file = "Heatmaps/Gene_modules_by_perturbation_forLegend.png")

aggregate_expression_ph <- pheatmap(agg_mat, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         clustering_method="ward.D2",
         fontsize=6, 
         #color = viridis::magma(35),
         color = colorRampPalette(c("blue", "white", "red"))(paletteLength),
         legend = FALSE,
         breaks = myBreaks,
         width = 1.5,
         height = 1.25,
         show_colnames = TRUE,
         #annotation_col = data.frame(row.names =  row.names(agg_mat_annotation),
         #                             Response = agg_mat_annotation$expected_response),
         treeheight_col = 5,
         treeheight_row = 5,
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize_row = 6,
         fontsize_col = 6,
         file = "Heatmaps/Gene_modules_by_perturbation_Figure_2D.png")

## Load Gene Set Collections
hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")

Ensembl_GSAlist <- as.matrix(rowData(TMZ_cds[expressed_genes,])$gene_short_name)
rownames(Ensembl_GSAlist)<-row.names(rowData(TMZ_cds[expressed_genes,]))
colnames(Ensembl_GSAlist) <- c("gene_short_name")
Ensembl_GSAlist<-Ensembl_GSAlist[,1]
Ensembl_GSAlist<-toupper(Ensembl_GSAlist)
length(Ensembl_GSAlist)

gene_module_GSA <- as.character(gene_module_df$module)
names(gene_module_GSA) <- gene_module_df$id

hallmarks_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                         replace_gene_names_vec(gene_module_GSA,
                                                                                Ensembl_GSAlist),
                                                         hallmarksGSC)

dir.create("GSA/DEG_modules")
gsea_bar_plots(hallmarks_GSAhyper, qval_cutoff = 0.1, pattern = "HALLMARK_", 
               width = 8, height = 10, sample = "GSA/DEG_modules/", gsc = "Hallmarks")

#### Identifying magnitude of perturbation effects in the absence of genotype DEG testing
NTC_cds <- TMZ_cds[,colData(TMZ_cds)$gene_id == "NTC"]
NTC_cds <- estimate_size_factors(NTC_cds)

TMZ_dose_response_diff_test <- fit_models(NTC_cds[expressed_genes,],
                                          model_formula_str = "~log(dose + 0.01)",
                                          cores = 1)

TMZ_dose_response_diff_test <- coefficient_table(TMZ_dose_response_diff_test) %>%
  dplyr::select(-model,-model_summary)


# saveRDS(TMZ_dose_response_diff_test, "NTC_TMZ_dose_response_diff_test.rds")
#TMZ_dose_response_diff_test <- readRDS("NTC_TMZ_dose_response_diff_test.rds")
#write.table(TMZ_dose_response_diff_test, "NTC_TMZ_dose_diff_test.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)

TMZ_dose_response_diff_test <- TMZ_dose_response_diff_test %>%
  filter(!is.na(p_value)) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

TMZ_dose_response_diff_test %>%  
  filter(grepl("dose", term)) %>%
  ggplot() +
  geom_point(aes(x = normalized_effect,
                 y = -log(q_value),
                 color = q_value < 0.01 & abs(normalized_effect) > 0.1),
             size = 0.1, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  geom_hline(yintercept = -log(0.01), color = "black", linetype = "dotted",size = 0.1) +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"),
        legend.key.height = unit(0.2,"line")) +
  xlab("Effect size") +
  ylab("-log(p adj)") +
  scale_color_manual("Pass\nFilter\n(FDR < 1% &\n|Beta| > 0.1)", values = c("TRUE" = "red", "FALSE" = "black")) +
  guides(guides(colour = guide_legend(override.aes = list(size=2))))
ggsave("DEG_testing/Effect_size_distribution_NTC_TMZ_dose_response_diff_test_Supplementary_Figure_3F.png", 
         dpi = 600, width = 2, height = 1.25)

TMZ_dose_response_sig_genes <- TMZ_dose_response_diff_test %>%
  filter(grepl("dose",term) &  q_value < 0.01 & abs(normalized_effect) > 0.1) %>%
  pull(id) %>%
  unique()

### Calculate the pairwise transcriptome distance of every cell to NTC cells treated with 100 mM TMZ
### across dose-responsive genes
TMZ_dose_response_angular_distance_df <- calculate_pairwise_perturbed_angular_distances(TMZ_cds,
                                                                                        TMZ_dose_response_sig_genes,
                                                                                        cores = detectCores()-2)

# saveRDS(TMZ_dose_response_angular_distance_df,"TMZ_dose_response_angular_distance_to_top_dose_df.rds")
# TMZ_dose_response_angular_distance_df <- readRDS("TMZ_dose_response_angular_distance_to_top_dose_df.rds")

dir.create("Angular_distance")

ggplot(TMZ_dose_response_angular_distance_df, aes(x = log10(js.distance), y = log10(ang.distance))) +
  geom_point(size =  0.2, stroke = 0) +
  geom_smooth(method = "lm", formula = y~x, size = 0.25, color = "brown4") +
  ggpubr::stat_cor(label.x = -0.32, label.y = -1,  size = 1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 4)) +
  xlab("log10(Jensen-Shannon\ndistance)") +
  ylab("log10(Angular distance)") 
ggsave("Angular_distance/Relationship_ang_and_JS_distance_Supplementary_Figure_3G.png",
         width = 1.15,
         height = 1.2,
         dpi = 600)

median_NTC_0 <- TMZ_dose_response_angular_distance_df %>% 
  filter(dose ==  "0" & gene_id == "NTC") %>%
  summarize(median_dist  = median(log(ang.distance)))

median_NTC_10 <- TMZ_dose_response_angular_distance_df %>% 
  filter(dose ==  "10" & gene_id == "NTC") %>%
  summarize(median_dist  = median(log(ang.distance)))

median_NTC_100 <- TMZ_dose_response_angular_distance_df %>% 
  filter(dose ==  "100" & gene_id == "NTC") %>%
  summarize(median_dist  = median(log(ang.distance)))

ggplot(TMZ_dose_response_angular_distance_df %>% filter(dose %in%  c("0","10","100")),
       aes(x = log(ang.distance), fill = factor(gene_id, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")))) +
  geom_density() +
  geom_vline(data = TMZ_dose_response_angular_distance_df %>% filter(dose %in%  c("0")),
             aes(xintercept = median_NTC_0$median_dist), linetype = "dashed", color = "grey70", size = 0.5) +
  geom_vline(data = TMZ_dose_response_angular_distance_df %>% filter(dose %in%  c("10")),
             aes(xintercept = median_NTC_10$median_dist), linetype = "dashed", color = "grey70", size = 0.5) +
  geom_vline(data = TMZ_dose_response_angular_distance_df %>% filter(dose %in%  c("100")),
             aes(xintercept = median_NTC_100$median_dist), linetype = "dashed", color = "grey70", size = 0.5) +
  facet_wrap(factor(gene_id, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")) ~ dose, ncol = 3, scales = "free_y") +
  scale_fill_viridis_d() +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  xlab("log(angular distance)") +
  ylab("Cell density")
ggsave("Angular_distance/Distribution_of_angular_distance_0_10_100uM_TMZ_Supplementary_Figure_3I.png",
         width = 3,
         height = 4,
         dpi = 600)

angular_distance_summary_df <- TMZ_dose_response_angular_distance_df %>%
  group_by(dose, gene_id) %>%
  summarize(median_ang_distance = median(ang.distance))

ggplot(TMZ_dose_response_angular_distance_df, 
       aes(x = factor(dose, levels = c("0","0.1","0.5","1","5","10","50","100")), 
           y = log(ang.distance), 
           fill = factor(dose, levels = c("0","0.1","0.5","1","5","10","50","100")))) +
  geom_violin(color = "black", size = 0.1) +
  facet_wrap(~factor(gene_id,levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), nrow = 1) +
  stat_summary(fun=mean, geom="point", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual("Dose (µM)",
                    values = c("0"="gray", "0.1"="#1D1147FF", 
                               "0.5"="#51127CFF", "1"="#822681FF",
                               "5"="#B63679FF","10"="#E65164FF", 
                               "50" = "#FB8861FF", "100"="#FEC287FF")) +
  theme(text = element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_text(size = 4, angle = 45, hjust = 1)) +
  xlab("[Temozolomide] (µM)") +
  ylab("Pairwise angular\ndistance to 100 µM NTC")
  #guides(guides(fill = guide_legend(override.aes = list(size=2))))
ggsave("Angular_distance/Angular_distance_to_mean_100uM_NTC_TMZ_sig_genes_Figure_2H.png", dpi = 600, height = 1.25, width = 3.5)

ggplot(TMZ_dose_response_angular_distance_df, 
       aes(x = factor(dose, levels = c("0","0.1","0.5","1","5","10","50","100")), 
           y = log(js.distance), 
           fill = factor(dose, levels = c("0","0.1","0.5","1","5","10","50","100")))) +
  geom_violin(color = "black", size = 0.1) +
  facet_wrap(~factor(gene_id,levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), nrow = 1) +
  stat_summary(fun=mean, geom="point", size = 0.1) +
  monocle3:::monocle_themxe_opts() +
  scale_fill_manual("Dose (µM)",
                    values = c("0"="gray", "0.1"="#1D1147FF", 
                               "0.5"="#51127CFF", "1"="#822681FF",
                               "5"="#B63679FF","10"="#E65164FF", 
                               "50" = "#FB8861FF", "100"="#FEC287FF")) +
  theme(text = element_text(size = 6),
        legend.position = "none",
        axis.text.x = element_text(size = 4, angle = 45, hjust = 1)) +
    xlab("[Temozolomide] (µM)") +
  ylab("Pairwise Jensen-Shannon\ndistance to 100 µM NTC")
  #guides(guides(fill = guide_legend(override.aes = list(size=2))))
ggsave("Angular_distance/Jensen-Shannon_distance_to_mean_100uM_NTC_TMZ_sig_genes_Supplementary_Figure_2J.png", dpi = 600, height = 1.25, width = 3.5)

TMZ_dose_response_angular_distance_matrix <- TMZ_dose_response_angular_distance_df %>%
  group_by(gene_id,dose) %>%
  dplyr::summarise(mean_angular_distance =  mean(log(ang.distance))) %>%
  tidyr::spread(key = gene_id,value = mean_angular_distance)

TMZ_dose_response_angular_distance_matrix <- as.data.frame(TMZ_dose_response_angular_distance_matrix)
row.names(TMZ_dose_response_angular_distance_matrix) <- TMZ_dose_response_angular_distance_matrix$dose
TMZ_dose_response_angular_distance_matrix$dose <- NULL

TMZ_dose_response_angular_distance_matrix <- apply(TMZ_dose_response_angular_distance_matrix,2,
                                                   function(x){x-x[1]})
TMZ_dose_response_angular_distance_matrix <- t(TMZ_dose_response_angular_distance_matrix)

pheatmap(TMZ_dose_response_angular_distance_matrix[c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2"),
                                                   c("0","0.1","0.5","1","5","10","50","100")],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = viridis::magma(35),
         clustering_method = "ward.D2",
         treeheight_row = 10,
         gaps_row = c(3,5),
         scale = "none",
         width = 2.2,
         height = 1.25,
         fontsize = 6,
         file = "Angular_distance/Angular_distance_to_NTC_heatmap_Supplementary_Figure_3H.png")

### Summarize the effect of perturbation on angular distance to 100 µM NTC and obtain fit information for 
### each genotype
angular_distance_summary <- TMZ_dose_response_angular_distance_df %>%
  mutate(ang.distance = as.numeric(ang.distance), dose = as.numeric(dose)) %>%
  group_by(dose, gene_id) %>%
  dplyr::summarise(mean_angular_distance = mean(ang.distance),
                   sd_angular_distance = sd(ang.distance))  %>%
  group_by(gene_id) %>%
  dplyr::mutate(norm_mean_angular_distance = mean_angular_distance/mean_angular_distance[dose == "0"]) %>%
  ungroup() %>%
  mutate(dose = as.numeric(dose))

y_min <- min(log10(angular_distance_summary$norm_mean_angular_distance))

ggplot(angular_distance_summary, aes(x = log10(dose + 0.01), y =  log10(norm_mean_angular_distance))) +
  geom_point(size = 0.5) +
  facet_wrap(~factor(gene_id, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), ncol = 3) +
  geom_smooth(method = "lm", formula = y~x, color = "brown4", size = 0.5) +
  ggpubr::stat_regline_equation(size = 1.5, label.y = -0.25) +
  ylim(y_min,NA) +
  monocle3:::monocle_theme_opts() +
  xlab("log10(Dose + 0.01)") +
  ylab("log10(Normalized mean angular distance)") +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.4,"line"),
        legend.key.height = unit(0.4,"line"))
ggsave("Angular_distance/Mean_angular_distances_by_pertubation_Supplementary_Figure_3K.png",
         width =  2.5, height = 3, dpi =  600)

### Repeat the above fits downsampling to 25 cells per genotype of interest
set.seed(2016L)
TMZ_dose_response_angular_distance_df_downsampled <- TMZ_dose_response_angular_distance_df %>%
  group_by(gene_id,dose) %>%
  filter(cell_1 %in% sample(cell_1, 25, replace = FALSE))

angular_distance_summary_downsampled <- TMZ_dose_response_angular_distance_df_downsampled %>%
  mutate(ang.distance = as.numeric(ang.distance), dose = as.numeric(dose)) %>%
  group_by(dose, gene_id) %>%
  dplyr::summarise(mean_angular_distance = mean(ang.distance),
                   sd_angular_distance = sd(ang.distance))  %>%
  group_by(gene_id) %>%
  dplyr::mutate(norm_mean_angular_distance = mean_angular_distance/mean_angular_distance[dose == "0"]) %>%
  #dplyr::mutate(norm_mean_angular_distance = ifelse(norm_mean_angular_distance < 1,norm_mean_angular_distance,1)) %>%
  ungroup() %>%
  mutate(dose = as.numeric(dose))

y_min_downsampled <- min(log10(angular_distance_summary_downsampled$norm_mean_angular_distance))

ggplot(angular_distance_summary_downsampled, aes(x = log10(dose + 0.01), y =  log10(norm_mean_angular_distance))) +
  geom_point(size = 0.5, stroke = 0) +
  facet_wrap(~factor(gene_id, levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")), ncol = 3) +
  geom_smooth(method = "lm", formula = y~x, color = "brown4", size = 0.25) +
  ggpubr::stat_regline_equation(size = 0.7, label.y = -0.25) +
  #ggpubr::stat_cor(label.y = -0.12, size = 1)+ 
  ylim(y_min_downsampled,NA) +
  monocle3:::monocle_theme_opts() +
  xlab("log10(Dose + 0.01)") +
  ylab("log10(Normalized mean angular distance)") +
  theme(text = element_text(size = 6),
        axis.text = element_blank(),
        legend.key.width = unit(0.4,"line"),
        legend.key.height = unit(0.4,"line"))
ggsave("Angular_distance/Mean_angular_distances_by_pertubation_downsampled_Supplementary_Figure_3L.png",
       width = 1.25, height = 2, dpi =  600)

#### Extract the transcriptional effective concentration 50 (TC50) for NTC
NTC_summary <- angular_distance_summary[angular_distance_summary$gene_id == "NTC",]
NTC_drm_model <- drc::drm(mean_angular_distance ~ dose, 
    data = NTC_summary, 
    fct=drc::LL.4(), 
    type = "Poisson")

NTC_drm_model_coef <- coef(NTC_drm_model)
NTC_drm_model_ed <- drc::ED(NTC_drm_model,50, interval="delta", display=FALSE)

NTC_EC50 <- NTC_drm_model_ed[1]

#### Extract the transcriptional effective concentration 50 (TC50) for NTC after downsampling
NTC_summary_downsampled <- angular_distance_summary_downsampled[angular_distance_summary_downsampled$gene_id == "NTC",]
NTC_drm_model_downsampled <- drc::drm(mean_angular_distance ~ dose, 
                          data = NTC_summary_downsampled, 
                          fct=drc::LL.4(), 
                          type = "Poisson")

NTC_drm_model_coef_downsampled <- coef(NTC_drm_model_downsampled)
NTC_drm_model_ed_downsampled <- drc::ED(NTC_drm_model_downsampled,50, interval="delta", display=FALSE)

NTC_EC50_downsampled <- NTC_drm_model_ed_downsampled[1]

### Extrapolate the amount of TMZ necessary for each genotype to reach the pairwise transcriptome distance to control at the NTC TC50
### This was calculated from the fits above
Dq_df <- data.frame(row.names = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2"),
                    gene_id = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2"),
                    Inferred_EC50 = c(7.260921,11.94835032,15.61118201,20130.32185,3050842.216,9.27833e14,9.3116e49))

Dq_df <- Dq_df %>% mutate(gene_id = factor(gene_id,levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")))

max_molarity <- 55.49*1e6

ggplot(Dq_df,aes(x = gene_id, y = log10(Inferred_EC50),  fill = gene_id)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_hline(yintercept = log10(max_molarity), color = "black", size = 0.125, linetype = "dashed") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle  = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual("Perturbation",values = c("NTC" = "grey70","MGMT" = "grey70", "MSH3" = "grey70",
                                          "MSH2" = "darkorange2","MSH6" = "darkorange2",
                                          "MLH1" = "darkorange2", "PMS2" = "darkorange2")) +
  scale_y_log10() +
  coord_flip() +
  xlab("Perturbation") +
  ylab("Inferred Transcriptional\nEC50 (µM)")
ggsave("Angular_distance/Relative_EC50_of_angular_distance_fits_Figure_2I.png",
         width = 1.5, height = 1.25, dpi = 600)

# Recalculate inferred TC50s for the downsampled data
Dq_df_downsampled <- data.frame(row.names = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2"),
                    gene_id = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2"),
                    Inferred_EC50 = c(8.579367,6.33920968,33.71360704,778.9984155,8529.461533,1.00286e51,1.52927e27))

Dq_df_downsampled <- Dq_df_downsampled %>% mutate(gene_id = factor(gene_id,levels = c("NTC","MGMT","MSH3","PMS2","MLH1","MSH6","MSH2")))

ggplot(Dq_df_downsampled,aes(x = gene_id, y = log10(Inferred_EC50),  fill = gene_id)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_hline(yintercept = log10(max_molarity), color = "black", size = 0.125, linetype = "dashed") +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle  = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual("Perturbation",values = c("NTC" = "grey70","MGMT" = "grey70", "MSH3" = "grey70",
                                              "MSH2" = "darkorange2","MSH6" = "darkorange2",
                                              "MLH1" = "darkorange2", "PMS2" = "darkorange2")) +
  scale_y_log10() +
  coord_flip() +
  xlab("Perturbation") +
  ylab("Inferred Transcriptional\nEC50 (µM)")
ggsave("Angular_distance/Relative_EC50_of_angular_distance_fits_downsampled_Supplementary_Figure_3M.png",
       width = 1.5, height = 1.25, dpi = 600)
  