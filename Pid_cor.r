true_persister_gene_set_rel <- "true_persister_gene_set_rel.Rds"
true_persister_gene_set_d15 <- "true_persister_gene_set_d15.Rds"
library(SC3)
library(Seurat)
library(tidyverse)
library(pheatmap)
load("ALL_data_Feb22.Rda")

pid <- "ALL3"

dat.sub <- dat[,dat$dna_cell_type=="blasts"  & 
                 dat$rna_cell_type=="blasts" & 
                 dat$patient_id==pid & 
                 dat$cell_phase == "G1" &
                 dat$timepoint_short== "d0"]

# Split ALL3 cells at diagnosis by clone ("dna_best_class") and generate correlation matrix of persister 
# genes per clone

clones <- unique(dat.sub$dna_best_class)

counts <- lapply(clones, function(x){
  c <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@counts)
}) 

names(counts) <- lapply(clones, function(x) x)

counts_per <- list()

for (i in 1:length(clones)){
  counts_per[[i]] <- counts[[i]][true_persister_gene_set_d15,]
}

for (i in 1:length(clones)){
  counts_per[[i]] <- (counts_per[[i]][rowSums(counts_per[[i]])>0,])
}

gene_list_for_scaled <- lapply(counts_per, function(x){
  rownames(x)
})

counts_per_scale <- lapply(clones, function(x){
  df <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@scale.data)
})

names(counts_per_scale) <- lapply(clones, function(x) x)

# Subset count data-frame to only include genes that are expressed for each patient
for (i in 1:length(clones)){
  counts_per_scale[[i]] <- counts_per_scale[[i]][gene_list_for_scaled[[i]],]
}


res_per <- lapply(counts_per_scale, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

filtered_res_per <- list()

for (i in 1:length(clones)){
  filtered_res_per[[i]] <- res_per[[i]] %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2 & abs(Freq) > 0.35) %>%
    filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2), pmin(as.character(Var1), as.character(Var2)))))) %>%
    arrange(desc(Freq))
}

names(filtered_res_per) <- lapply(clones, function(x) x)

######################################################################################################################
# Correlation heatmaps

paletteLength <- 100
plots_per <- lapply(res_per, function(x){
  pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T, 
           color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
           breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
})

######################################################################################################################
# Same but for relapse genes 
counts_rel <- list()

for (i in 1:length(clones)){
  counts_rel[[i]] <- counts[[i]][true_persister_gene_set_rel,]
}

for (i in 1:length(clones)){
  counts_rel[[i]] <- (counts_rel[[i]][rowSums(counts_rel[[i]])>0,])
}

gene_list_for_scaled_rel <- lapply(counts_rel, function(x){
  rownames(x)
})

counts_rel_scale <- lapply(clones, function(x){
  df <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@scale.data)
})

names(counts_rel_scale) <- lapply(clones, function(x) x)

# Subset count data-frame to only include genes that are expressed for each patient
for (i in 1:length(clones)){
  counts_rel_scale[[i]] <- counts_rel_scale[[i]][gene_list_for_scaled[[i]],]
}


res_rel <- lapply(counts_rel_scale, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

filtered_res_rel <- list()

for (i in 1:length(clones)){
  filtered_res_rel[[i]] <- res_rel[[i]] %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2 & abs(Freq) > 0.35) %>%
    filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2), pmin(as.character(Var1), as.character(Var2)))))) %>%
    arrange(desc(Freq))
}

names(filtered_res_rel) <- lapply(clones, function(x) x)

######################################################################################################################

paletteLength <- 100
plots_rel <- lapply(res_rel, function(x){
  pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T, 
           color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
           breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
})

