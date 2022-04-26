# Correlation in ALL3
library(tidyverse)
library(scran)
library(amap)

dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
dat.sub <- dat.sub[,dat.sub$patient_id=="ALL3" & !is.na(dat.sub$dna_cell_type)]
clones <- unique(dat.sub$dna_best_class)

persister_genes <- readRDS("persister_genes.Rds")

df1 <- lapply(clones, function(x){
 df <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@counts)
 df2 <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@scale.data)
 df <- df[persister_genes,]
 df <- df[rowSums(df)>0,]
 genes <- rownames(df)
 df2 <- df2[genes,]
 return(df2)
}) 

names(df1) <- lapply(clones, function(x) x)

correlated_pairs <- lapply(df1, function(x){
  res <- correlatePairs(x)
  p.vals <- res$p.value
  res$hochberg <- p.adjust(p.vals, method = "hochberg")
  res$holm <- p.adjust(p.vals, method = "holm")
  res$bonferroni <- p.adjust(p.vals, method = "bonferroni")
  res <- as.data.frame(res[res$FDR <= 0.05,])
  return(res)
})


# On average, how well do persister genes correlate with each other? 

df1 <- lapply(clones, function(x){
  df <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@counts)
  df2 <- as.data.frame(dat.sub[,dat.sub$dna_best_class==x]@assays$RNA@scale.data)
  df <- df[rowSums(df)>0,]
  genes <- rownames(df)
  df2 <- df2[genes,]
  return(df2)
}) 

# kmeans clustering

clust <- lapply(df1, function(x){
  Kmeans(x, centers=25, iter.max = 100, method = "pearson")
})

names(clust) <- lapply(clones, function(x) x)

table(clust$ALL3_2$cluster[names(clust$ALL3_2$cluster) %in% persister_genes])
table(clust$ALL3_4$cluster[names(clust$ALL3_4$cluster) %in% persister_genes])

names(df1) <- lapply(clones, function(x) x)

cor <- lapply(df1, function(x){
  res <- correlatePairs(x)
  p.vals <- res$p.value
  res <- as.data.frame(res[res$FDR <= 0.05,])
  return(res)
})

cor <- lapply(df1, function(x){
  cor <- corr.test(x = as.matrix(t(x[rownames(x) %in% persister_genes,])),
                   y = as.matrix(t(x[!rownames(x) %in% persister_genes,])),
                   method = "spearman")
})



# How do other genes correlate with persister genes at diagnosis? 


