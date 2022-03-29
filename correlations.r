# correlations
library(tidyverse)
load("ALL_data_Feb22.Rda")
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
true_persister_gene_set_rel <- readRDS("true_persister_gene_set_rel_new1.Rds")
true_persister_gene_set_d15 <- readRDS("true_persister_gene_set_d15_new1.Rds")
geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel))
dat.sub <- dat.sub[,dat.sub$cell_phase=="G1" & dat.sub$timepoint=="diagnosis"]
geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel, pseudo_persister, pseudo_relapse))

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

# Pseudo matrices -- One with pseudo_persister genes only (1), one with pseudo_relapse genes (2) only and one with pseudo_persister genes
# as rows with persister genes on columns (3) and another with pseudo_relapse genes on rows with relapse_genes on columns (4) 

res_per <- lapply(res, function(x){
  res <- x[intersect(rownames(x), true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_d15)]
})

res_rel <- lapply(res, function(x){
  res <- x[intersect(rownames(x), true_persister_gene_set_rel), intersect(colnames(x), true_persister_gene_set_rel)]
})

res_combined<- lapply(res, function(x){
  res <- x[intersect(rownames(x),true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_rel)]
})

pseudo_persister_only <- lapply(res, function(x){
  res <- x[intersect(rownames(x), pseudo_persister), intersect(colnames(x), pseudo_persister)]
})

pseudo_rel_only <- lapply(res, function(x){
  res <- x[intersect(rownames(x), pseudo_relapse), intersect(colnames(x), pseudo_relapse)]
})

pseudo_persister_mix <- lapply(res, function(x){
  res <- x[intersect(rownames(x), pseudo_persister), intersect(colnames(x), true_persister_gene_set_d15)]
})

pseudo_rel_mix <- lapply(res, function(x){
  res <- x[intersect(rownames(x), pseudo_relapse), intersect(colnames(x), true_persister_gene_set_rel)]
})


######################################################################################################################
# Average correlation of each persister gene to each relapse gene 

persister_genes_x_coord <- lapply(res_combined, function(x){
  rowMeans(x)
})

relapse_genes_y_coord <- lapply(res_combined, function(x){
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

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/new_3_names", x, "_avg_cor_d0.pdf", sep = ""),
         plot = plots[[x]])
})

######################################################################################################################
# pseudo 

pseudo_rel_x_coord <- lapply(pseudo_rel_only, function(x){
  rowMeans(x)
})

pseudo_rel_y_coord <- lapply(pseudo_rel_mix, function(x){
  rowMeans(x)
})

relapse_x_coord <- lapply(pseudo_rel_mix, function(x){
  colMeans(x)
})

relapse_y_coord <- lapply(res_rel, function(x){
  colMeans(x)
})
 
pseudo_persister_x_coord <- lapply(pseudo_persister_only, function(x){
  rowMeans(x)
})

pseudo_persister_y_coord <- lapply(pseudo_persister_mix, function(x){
  rowMeans(x)
})

persister_x_coord <- lapply(pseudo_persister_mix, function(x){
  colMeans(x)
})

persister_y_coord <- lapply(res_per, function(x){
  colMeans(x)
})

######################################################################################################################
persister_coordinates <- list ()
for (i in names){
  persister_coordinates[[i]] <- data.frame(persister_x_coord[[i]], persister_y_coord[[i]], "persister")
}

names(persister_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(persister_coordinates[[i]]) <- c("x", "y", "geneset")
}

pseudo_persister_coordinates <- list()

for (i in names){
  pseudo_persister_coordinates[[i]] <- data.frame(pseudo_persister_x_coord[[i]], pseudo_persister_y_coord[[i]], "pseudo_persister")
}
names(pseudo_persister_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(pseudo_persister_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(persister_coordinates[[i]], (pseudo_persister_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab("Average correlation with pseudo persister genes") +
    ylab("Average correlation with persister genes") +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1)
})

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/pseudo/persister/", i, ".pdf", sep = ""),
         plot = plots[[i]])
}

######################################################################################################################
relapse_coordinates <- list ()
for (i in names){
  relapse_coordinates[[i]] <- data.frame(relapse_x_coord[[i]], relapse_y_coord[[i]], "relapse")
}

names(relapse_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(relapse_coordinates[[i]]) <- c("x", "y", "geneset")
}

pseudo_relapse_coordinates <- list()

for (i in names){
  pseudo_relapse_coordinates[[i]] <- data.frame(pseudo_rel_x_coord[[i]], pseudo_rel_y_coord[[i]], "pseudo_relapse")
}
names(pseudo_relapse_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(pseudo_relapse_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(relapse_coordinates[[i]], (pseudo_relapse_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab("Average correlation with pseudo relapse genes") +
    ylab("Average correlation with relapse genes") +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1)
})

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/pseudo/relapse/", i, ".pdf", sep = ""),
         plot = plots[[i]])
}
######################################################################################################################

# To address internal correlation structures 
# Peristers 

all.combos <- list()
for (i in names){
  all.combos[[i]] <- res_per[[i]] %>%
         as.table() %>% as.data.frame() %>%
         subset(Var1 != Var2)  %>%
    mutate(combined=paste(Var1, Var2, sep = "."))
}
names(all.combos) <- lapply(names, function(x) x)

cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)


df_for_plot <- apply(cmb, 2, function(i){
  x <- all.combos[[i[1]]]
  y <- all.combos[[i[2]]]
  common <- intersect(x$combined, y$combined)
  df <- data.frame(x[x$combined %in% common,], y[y$combined %in% common,])
  df <- df[, c("Freq", "Freq.1", "combined")]
  return(df)
})

names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(names[x[1]], names[x[2]], sep = ".")
})
  
plots <- lapply(names(df_for_plot), function(x){
  a <- gsub("\\..*", "", x)
  b <- gsub(".*\\.", "", x)
  ggplot(df_for_plot[[x]], aes(Freq, Freq.1, color=combined)) +
    geom_point() +
    geom_smooth(method = lm, se=FALSE, col='red', size=2) +
    xlab(paste("Gene pair correlations in ", a, sep = "")) +
    ylab(paste("Gene pair correlations in ", b, sep = "")) +
    theme(legend.position = "none")
})

names(plots) <- lapply(names(df_for_plot), function(x) x)

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/correlation_structure/persister/", i, ".pdf", sep = ""),
         plot = plots[[i]])
}


######################################################################################################################
# Relapse
rel_genes.combos <- list()
for (i in names){
  rel_genes.combos[[i]] <- res_rel[[i]] %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2)  %>%
    mutate(combined=paste(Var1, Var2, sep = "."))
}
names(rel_genes.combos) <- lapply(names, function(x) x)

cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)


df_for_plot <- apply(cmb, 2, function(i){
  x <- rel_genes.combos[[i[1]]]
  y <- rel_genes.combos[[i[2]]]
  common <- intersect(x$combined, y$combined)
  df <- data.frame(x[x$combined %in% common,], y[y$combined %in% common,])
  df <- df[, c("Freq", "Freq.1", "combined")]
  return(df)
})

names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(names[x[1]], names[x[2]], sep = ".")
})

plots <- lapply(names(df_for_plot), function(x){
  a <- gsub("\\..*", "", x)
  b <- gsub(".*\\.", "", x)
  ggplot(df_for_plot[[x]], aes(Freq, Freq.1, color=combined)) +
    geom_point() +
    geom_smooth(method = lm, se=FALSE, col='red', size=2) +
    xlab(paste("Gene pair correlations in ", a, sep = "")) +
    ylab(paste("Gene pair correlations in ", b, sep = "")) +
    theme(legend.position = "none")
})

names(plots) <- lapply(names(df_for_plot), function(x) x)

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/correlation_structure/relapse/", i, ".pdf", sep = ""),
         plot = plots[[i]])
}
