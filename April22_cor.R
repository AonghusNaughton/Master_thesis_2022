library(tidyverse)
library(scran)
library(amap)
library(beeswarm)

dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
persister_genes <- readRDS("persister_genes.Rds")
names <- unique(dat.sub$patient_id)

df <- lapply(names, function(x){
  if (x=="ALL3"){
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@counts)
    df2 <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@scale.data)
  } else {
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)
    df2 <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@scale.data)
  }
  df <- df[rowSums(df)>0,]
  genes <- rownames(df)
  df2 <- df2[genes,]
  return(df2)
}) 

names(df) <- lapply(names, function(x){
  x
})

clust <- lapply(df, function(x){
  Kmeans(x, centers=25, iter.max = 100, method = "pearson")
})

names(clust) <- lapply(names, function(x) x)

cluster_plots <- lapply(clust, function(x){
  x$cluster %>%
    as.table() %>% as.data.frame() %>%
    filter(Var1 %in% persister_genes) %>%
    mutate(x="") %>%
    ggplot(aes(x=x, y=Freq)) +
    geom_beeswarm() +
    theme_classic() +
    labs(x="Proportion of persister genes per cluster",
         y="Cluster")
})

saveRDS(clust, "kmeans_clusters_allpids_all_genes.Rds")

clusters <- lapply(clust, function(x){
  x$cluster %>%
    as.table() %>% as.data.frame() %>%
    filter(Var1 %in% persister_genes)
})

clusters$ALL1$Var1[clusters$ALL1$Freq==15]
clusters$ALL1$Var1[clusters$ALL1$Freq==1]
