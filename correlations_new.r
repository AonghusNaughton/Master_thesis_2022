library(tidyverse)
library(pheatmap)

dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
true_persister_gene_set_rel <- readRDS("true_persister_gene_set_rel_new1.Rds")
true_persister_gene_set_d15 <- readRDS("true_persister_gene_set_d15_new1.Rds")
B_cell_states <- readRDS("B_cell_states.rds")
geneSets <- unique(c(true_persister_gene_set_d15, true_persister_gene_set_rel, B_cell_states$`Immature-B`, B_cell_states$`Mature-B`, B_cell_states$HSPC))
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


# saveRDS(res, "correlation_matrix.Rds")

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
  }),
  res_HSPC <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$HSPC), intersect(colnames(x), B_cell_states$HSPC)]
  }),
  res_HSPC.Mature <- lapply(res, function(x){
    df <- x[intersect(rownames(x), B_cell_states$HSPC), intersect(colnames(x), B_cell_states$`Mature-B`)]
  })
)

names(res_list) <- c("res_Immature.Mature", "res_Immature.Persister",
                     "res_Immature.Relapse", "res_Mature.Persister", 
                     "res_Mature.Relapse", "res_Immature",
                     "res_Mature", "res_Persister",
                     "res_Relapse", "res_Persister.Relapse",
                     "res_HSPC", "res_HSPC.Mature")

saveRDS(res_list, "res_list.Rds")

plots <- list()
paletteLength <- 100
for (i in names(res_list)){
  plots[[i]] <- lapply(res_list[[i]], function(x){
    pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T,
             color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
             breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1),
                        seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
  })
}

# lapply(names(plots), function(x){
#   dir.create(paste("/Users/aonghusnaughton/Proj_eng/March22/", x, "_cor_heatmaps", sep = ""))
# })

lapply(names(plots), function(x){
  lapply(names(plots[[x]]), function(y){
    pdf(paste("/Users/aonghusnaughton/Proj_eng/March22/", x, "_cor_heatmaps/", y, ".pdf", sep = ""), height = 12, width = 12)
    grid::grid.newpage()
    grid::grid.draw(plots[[x]][[y]]$gtable)
    dev.off()
  })
})

######################################################################################################################

y <- "HSPC"
x <- "Mature"

# Average correlation of y-axis genes to x-axis genes 
assign(paste0(y, "_genes_x_coord"), lapply(eval(parse(text=paste0("res_list[['res_", y, ".", x, "']]"))), function(i){
  rowMeans(i)
}))

# Average correlation of y-axis genes to y-axis genes 
assign(paste0(y, "_genes_y_coord"), lapply(eval(parse(text=paste0("res_list[['res_", y,"']]"))), function(i){
  rowMeans(i)
}))


# Average correlation of x-axis genes x-axis genes 
assign(paste0(x, "_genes_x_coord"), lapply(eval(parse(text = paste0("res_list[['res_", x, "']]"))), function(i){
  rowMeans(i)
}))

# Average correlation of x-axis genes to y-axis genes 
assign(paste0(x, "_genes_y_coord"), lapply(eval(parse(text = paste0("res_list[['res_", y, ".", x, "']]"))), function(i){
  colMeans(i)
}))


######################################################################################################################

# Create data structure that holds x and y coordinates for each gene and plot as points on a coordinate plane.


y_gene_coordinates <- vector(mode = "list", length = length(names))
names(y_gene_coordinates) <- lapply(names, function(x) x)

for (i in names){
  y_gene_coordinates[[i]] <- data.frame(eval(parse(text = paste0(y, "_genes_x_coord")))[[i]],
                                        eval(parse(text = paste0(y, "_genes_y_coord")))[[i]], 
                                        y)
} 

for (i in names){
  colnames(y_gene_coordinates[[i]]) <- c("x", "y", "geneset")
}

x_gene_coordinates <- vector(mode = "list", length = length(names))

for (i in names){
  x_gene_coordinates[[i]] <- data.frame(eval(parse(text = paste0(x, "_genes_x_coord")))[[i]], 
                                        eval(parse(text = paste0(x, "_genes_y_coord")))[[i]], 
                                        x)
} 

for (i in names){
  colnames(x_gene_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(x_gene_coordinates[[i]], (y_gene_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab(paste0(x)) +
    ylab(paste0(y)) +
    ylim(-0.15,0.15) +
    xlim(-0.15,0.15) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_classic() +
    theme(legend.position = "right")
})

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/Average_correlations/HSPC_vs_Mature/", x, ".pdf", sep = ""),
         plot = plots[[x]],
         width = 8,
         height = 5)
})

######################################################################################################################

# compute similarity between persister matrices 
# Overall correlation

# 11 correlation matrices 
cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
names_cmb <- combn(unique(dat.sub$patient_id), 2)

matched_res <- apply(cmb, 2, function(x){
  df <- list(m1=res_list$res_Persister[[x[[1]]]][intersect(rownames(res_list$res_Persister[[x[[1]]]]), 
                                                   rownames(res_list$res_Persister[[x[[2]]]])), 
                                                 intersect(colnames(res_list$res_Persister[[x[[1]]]]),
                                                 colnames(res_list$res_Persister[[x[[2]]]]))],
             m2=res_list$res_Persister[[x[[2]]]][intersect(rownames(res_list$res_Persister[[x[[1]]]]), 
                                                        rownames(res_list$res_Persister[[x[[2]]]])),
                                                 intersect(colnames(res_list$res_Persister[[x[[1]]]]),
                                                 colnames(res_list$res_Persister[[x[[2]]]]))])
})

names(matched_res) <- apply(names_cmb, 2, function(x){
  paste(x[[1]], x[[2]], sep = ".")
})

lapply(matched_res, function(x){
  print(paste(length(rownames(x[["m1"]])), length(rownames(x[["m2"]])), sep = " "))
  print(paste(length(colnames(x[["m1"]])), length(colnames(x[["m2"]])), sep = " "))
})

extracted_tri.all.combos <- lapply(matched_res, function(x){
  list <- list(t(x[["m1"]])[lower.tri(t(x[["m1"]]))],
               t(x[["m2"]])[lower.tri(t(x[["m2"]]))])
  names(list) <- c("m1", "m2")
  return(list)
})

lapply(extracted_tri.all.combos, function(x){
  print(paste(length(x[["m1"]]), length(x[["m2"]]), sep = " "))
})

spearman_correlations <- lapply(extracted_tri.all.combos, function(x){
  cor.test(x[["m1"]], x[["m2"]], method = "spearman", exact=F)
})


# Clusters 
