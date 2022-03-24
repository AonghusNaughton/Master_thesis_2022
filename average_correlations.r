
# Average correlation of each gene to the genes in the same set (persister/relapse) or to the other set.

dat.sub <- dat[,(dat$dna_cell_type=="blasts" | is.na(dat$dna_cell_type))  & 
                            dat$rna_cell_type=="blasts" &
                            dat$cell_phase == "G1" &
                            dat$timepoint_short== "d0"]

true_persister_gene_set_rel <- readRDS("true_persister_gene_set_rel_new1.Rds")
true_persister_gene_set_d15 <- readRDS("true_persister_gene_set_d15_new1.Rds")

# Matrix with all genes from both sets

pids <- unique(dat.sub$patient_id)

counts <- lapply(pids, function(x){
  if (x=="ALL3"){
    c <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@counts)
  } else {
    c <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)
  }
  return(c)
})

names(counts) <- lapply(pids, function(x) x)

geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel))

counts_combined <- list()

for (i in pids){
  counts_combined[[i]] <- counts[[i]][geneSets,]
}

for (i in pids){
  counts_combined[[i]] <- (counts_combined[[i]][rowSums(counts_combined[[i]])>0,])
}

gene_list_for_scaled <- lapply(counts_combined, function(x){
  rownames(x)
})

counts_scale <- lapply(pids, function(x){
  if (x=="ALL3"){
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts" & dat.sub$rna_cell_type=="blasts"]@assays$RNA@scale.data)
  } else {
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@scale.data)
  }
  return(df)
})

names(counts_scale) <- lapply(pids, function(x) x)

counts_combined_scale <- list()

for (i in pids){
  counts_combined_scale[[i]] <- counts_scale[[i]][gene_list_for_scaled[[i]],]
}


res <- lapply(counts_combined_scale, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

######################################################################################################################


# Modify correlation matrix so that persister genes are on rows and relapse genes are on columns. 

res_modified <- lapply(res, function(x){
  res <- x[intersect(rownames(x),true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_rel)]
})

# Average correlation of each persister gene to each relapse gene 

persister_genes_x_coord <- lapply(res_modified, function(x){
  rowMeans(x)
})

relapse_genes_y_coord <- lapply(res_modified, function(x){
  colMeans(x)
})

######################################################################################################################

# Average correlation of each perister/relapse gene to all other genes in same set

res_relapse <- readRDS("res_relapse.Rds")
res_persister <- readRDS("res_persister.Rds")

persister_genes_y_coord <- lapply(res_persister, function(x){
  rowMeans(x)
})

relapse_genes_x_coord <- lapply(res_relapse, function(x){
  rowMeans(x)
})

######################################################################################################################
# Create data structure that holds x and y coordinates for each gene and plot as points on a coordinate plane.

persister_coordinates <- list()

for (i in pids){
  persister_coordinates[[i]] <- data.frame(persister_genes_x_coord[[i]], persister_genes_y_coord[[i]], "persister")
}

names(persister_coordinates) <- lapply(pids, function(x) x)

for (i in pids){
  colnames(persister_coordinates[[i]]) <- c("x", "y", "geneset")
}

relapse_coordinates <- list()

for (i in pids){
  relapse_coordinates[[i]] <- data.frame(relapse_genes_x_coord[[i]], relapse_genes_y_coord[[i]], "relapse")
}
names(relapse_coordinates) <- lapply(pids, function(x) x)

for (i in pids){
  colnames(relapse_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in pids){
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

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/new_3_pids/", x, "3pids.pdf", sep = ""), 
         plot = plots[[x]]) 
})

######################################################################################################################
# Make a plot for all patients 




  