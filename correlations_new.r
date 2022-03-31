dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
true_persister_gene_set_rel <- readRDS("true_persister_gene_set_rel_new1.Rds")
true_persister_gene_set_d15 <- readRDS("true_persister_gene_set_d15_new1.Rds")
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


# Correlation matrices for all combinations of gene sets

res_list <- list(
  res_Immature.Mature <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Immature-B`), intersect(colnames(x), B_cell_states$`Mature-B`)]
  }),
  res_Immature.Persister <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Immature-B`), intersect(colnames(x), true_persister_gene_set_d15)]
  }),
  res_Immature.Relapse <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Immature-B`), intersect(colnames(x), true_persister_gene_set_rel)]
  }),
  res_Mature.Persister <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Mature-B`), intersect(colnames(x), true_persister_gene_set_d15)]
  }),
  res_Mature.Relapse <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Mature-B`), intersect(colnames(x), true_persister_gene_set_rel)]
  }),
  res_Immature <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Immature-B`), intersect(colnames(x), B_cell_states$`Immature-B`)]
  }),
  res_Mature <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$`Mature-B`), intersect(colnames(x), B_cell_states$`Mature-B`)]
  }),
  res_Persister <- lapply(res, function(x){
    df <- x[intersect(rownames(x), true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_d15)]
  }),
  res_Relapse <- lapply(res, function(x){
    df <- x[intersect(rownames(x), true_persister_gene_set_rel), intersect(colnames(x), true_persister_gene_set_rel)]
  }),
  res_Persister.Relapse <- lapply(res, function(x){
    df <- x[intersect(rownames(x), true_persister_gene_set_d15), intersect(colnames(x), true_persister_gene_set_rel)]
  })
)

names(res_list) <- c("res_Immature.Mature", "res_Immature.Persister",
                     "res_Immature.Relapse", "res_Mature.Persister", 
                     "res_Mature.Relapse", "res_Immature",
                     "res_Mature", "res_Persister",
                     "res_Relapse", "res_Persister.Relapse")

plots <- list()
for (i in names(res_list)){
  plots[[i]] <- lapply(res_list[[i]], function(x){
    pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T,
             color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
             breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1),
                        seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
  })
}

lapply(names(plots), function(x){
  dir.create(paste("/Users/aonghusnaughton/Proj_eng/March22/", x, "_cor_heatmaps", sep = ""))
})

lapply(names(plots), function(x){
  lapply(names(plots[[x]]), function(y){
    pdf(paste("/Users/aonghusnaughton/Proj_eng/March22/", x, "_cor_heatmaps/", y, ".pdf", sep = ""), height = 10, width = 10)
    grid::grid.newpage()
    grid::grid.draw(plots[[x]][[y]]$gtable)
    dev.off()
  })
})

######################################################################################################################

# Average correlation of itself to persister genes 
Immature_genes_x_coord <- lapply(res_Immature.Relapse, function(x){
  rowMeans(x)
})

# Average correlation of itself to genes in same set 
Immature_genes_y_coord <- lapply(res_Immature, function(x){
  rowMeans(x)
})

# Average correlation of itself to genes in same set 
Relapse_genes_x_coord <- lapply(res_Relapse, function(x){
  rowMeans(x)
})

# Average correlation of itself to ImImmature b-cell genes 
Relapse_genes_y_coord <- lapply(res_Immature.Relapse, function(x){
  colMeans(x)
})

######################################################################################################################
# Create data structure that holds x and y coordinates for each gene and plot as points on a coordinate plane.

Immature_coordinates <- list()

for (i in names){
  Immature_coordinates[[i]] <- data.frame(Immature_genes_x_coord[[i]], Immature_genes_y_coord[[i]], "Immature")
}

names(Immature_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(Immature_coordinates[[i]]) <- c("x", "y", "geneset")
}

Relapse_coordinates <- list()

for (i in names){
  Relapse_coordinates[[i]] <- data.frame(Relapse_genes_x_coord[[i]], Relapse_genes_y_coord[[i]], "Relapse")
}
names(Relapse_coordinates) <- lapply(names, function(x) x)

for (i in names){
  colnames(Relapse_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(Immature_coordinates[[i]], (Relapse_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab("Average correlation with Relapse genes") +
    ylab("Average correlation with Immature genes") +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1)
})

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/Immature_vs_Relapse/", x, ".pdf", sep = ""),
         plot = plots[[x]],
         width = 8,
         height = 5)
})
