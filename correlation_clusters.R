dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
true_persister_gene_set_rel <- readRDS("true_persister_gene_set_rel_new1.Rds")
true_persister_gene_set_d15 <- readRDS("true_persister_gene_set_d15_new1.Rds")
B_cell_states <- readRDS("B_cell_states.rds")
geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel, B_cell_states$`Immature-B`, B_cell_states$`Mature-B`))
dat.sub <- dat.sub[,dat.sub$cell_phase=="G1" & dat.sub$timepoint=="diagnosis"]

names <- unique(dat.sub$patient_id)

df1 <- lapply(names, function(x){
  if (x=="ALL3"){
    c <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@counts)
  } else {
    c <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)
  }
  return(c)
}) 

names(df1) <- lapply(names, function(x){
  x
})

for (i in names){
  df1[[i]] <- df1[[i]][true_persister_gene_set_d15,]
}

for (i in names){
  df1[[i]] <- (df1[[i]][rowSums(df1[[i]])>0,])
}

gene_list_for_scaled <- lapply(df1, function(x){
  rownames(x)
})

df2 <- lapply(names, function(x){
  if (x=="ALL3"){
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@scale.data)
  } else {
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@scale.data)
  }
  return(df)
})

names(df2) <- lapply(names, function(x){
  x
})

for (i in names){
  df2[[i]] <- df2[[i]][gene_list_for_scaled[[i]],]
}


res <- lapply(df2, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})


clusters <- lapply(res, function(x){
  clust <- hclust(as.dist(1-abs(x)))
})

lapply(clusters, function(x){
  plot(x)
})

cmb <- combn(unique(dat.sub$patient_id), 2)

apply(cmb, 2, function(x){
  dendlist(as.dendrogram(clusters[[x[1]]], as.dendrogram(clusters[[x[2]]]))) %>%
    untangle(method = "step1side") %>%
    tanglegram()
})

dendlist(as.dendrogram(clusters$MRD_ALL71), as.dendrogram(clusters$MRD_ALL67)) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()                       # Draw the two dendrograms

a <- hclust(as.dist(1-abs(res$MRD_ALL71)))
plot(a)

plot(a, hang = -1, cex = 0.6)


nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(a,  xlab = "Height",
     nodePar = nodePar, horiz = TRUE)


heatmap.2(res$MRD_ALL71, main = "Hierarchical Cluster", dendrogram="column",trace="none",col=greenred(10))


