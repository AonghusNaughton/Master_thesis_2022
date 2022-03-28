# correlations
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
dat.sub <- dat.sub[,dat.sub$cell_phase=="G1" & dat.sub$timepoint=="diagnosis"]
geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel))

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
  df1[[i]] <- df1[[i]][geneSets,]
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

saveRDS(res, "cor_res_combined.Rds")

######################################################################################################################


# Three nmodified correlation matrices -- One for persister genes, one for relapse genes and one with persister genes on rows and 
# relapse genes on columns

res_per <- lapply(res, function(x){
  res <- x[intersect(rownames(x), true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_d15)]
})

res_rel <- lapply(res, function(x){
  res <- x[intersect(rownames(x), true_persister_gene_set_rel), intersect(colnames(x), true_persister_gene_set_rel)]
})

res_combined<- lapply(res, function(x){
  res <- x[intersect(rownames(x),true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_rel)]
})

######################################################################################################################
######################################################################################################################
# Average correlation of each persister gene to each relapse gene 

persister_genes_x_coord <- lapply(res_modified, function(x){
  rowMeans(x)
})

relapse_genes_y_coord <- lapply(res_modified, function(x){
  colMeans(x)
})

######################################################################################################################

# Average correlation of each perister/relapse gene to all other genes in same set

persister_genes_y_coord <- lapply(res_per, function(x){
  rowMeans(x)
})

relapse_genes_x_coord <- lapply(res_rel, function(x){
  rowMeans(x)
})

######################################################################################################################
# Create data structure that holds x and y coordinates for each gene and plot as points on a coordinate plane.

persister_coordinates <- list()

for (i in names){
  persister_coordinates[[i]] <- data.frame(persister_genes_x_coord[[i]], persister_genes_y_coord[[i]], "persister")
}

names(persister_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(persister_coordinates[[i]]) <- c("x", "y", "geneset")
}

relapse_coordinates <- list()

for (i in names){
  relapse_coordinates[[i]] <- data.frame(relapse_genes_x_coord[[i]], relapse_genes_y_coord[[i]], "relapse")
}
names(relapse_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(relapse_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(persister_coordinates[[i]], (relapse_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab("Average correlation with relapse genes") +
    ylab("Average correlation with persister genes") +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1)
})

# lapply(names(plots), function(x){
#   ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/ALL3_clonal_cor/", x, "_avg_cor_d0.pdf", sep = ""), 
#          plot = plots[[x]]) 
# })



