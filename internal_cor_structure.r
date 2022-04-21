# Pairwise gene correlations in each patient 
library(tidyverse)
library(Hmisc)
dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
names <- unique(dat.sub$patient_id)

persister_genes <- readRDS("persister_genes.Rds")
relapse_genes <- readRDS("relapse_genes.Rds")
genes <- unique(c(persister_genes, relapse_genes))

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

df1 <- lapply(df1, function(x){
  x <- x[genes,]
  x <- x[rowSums(x) > 0,]
  return(x)
})

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


var.corr <- lapply(df2, function(x){
  res <- correlatePairs(x)
})

sig.corr <- lapply(var.corr, function(x){
  x %>%
    as.data.frame() %>%
    filter(FDR <= 0.05) %>%
    filter(FDR != 0)
  # res <- as.data.frame(x[x$FDR <= 0.05,])
})

sig.corr <- lapply(sig.corr, function(x){
  x %>%
    mutate(gene_set=case_when(gene1 %in% persister_genes & gene2 %in% persister_genes ~ "Persister genes",
                              gene1 %in% persister_genes & gene2 %in% relapse_genes ~ "Mixed pair",
                              gene2 %in% persister_genes & gene1 %in% relapse_genes ~ "Mixed pair",
                              gene1 %in% relapse_genes & gene2 %in% relapse_genes ~ "Relapse genes")) %>%
    mutate(combined=paste(gene1, gene2, sep = "."))
})

cmb <- combn(c(unique(dat.sub$patient_id)), 2)

df_for_plot <- apply(cmb, 2, function(i){
  x <- sig.corr[[i[1]]]
  y <- sig.corr[[i[2]]]
  common <- intersect(x$combined, y$combined)
  dfx <- data.frame(x[x$combined %in% common,])
  order <- dfx$combined
  dfy <- data.frame(y[y$combined %in% common,])
  dfy <- dfy %>%
    arrange(factor(combined, levels=order))
  df <- data.frame(dfx, dfy)
  df <- df[, c("rho", "rho.1", "combined", "gene_set")]
  return(df)
})

names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(x[1], x[2], sep = ".")
})

plots <- lapply(names(df_for_plot), function(x){
    a <- gsub("\\..*", "", x)
    b <- gsub(".*\\.", "", x)
    ggplot(df_for_plot[[x]], aes(rho, rho.1, color=gene_set)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, col='red', size=1) +
      xlab(paste("Gene pair correlations in ", a, sep = "")) +
      ylab(paste("Gene pair correlations in ", b, sep = "")) +
      theme_classic() +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0)
  })

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/int_cor_structure/persister_and_relapse/p_val_filtered/", i, ".pdf", sep = ""),
         plot = plots[[i]],
         width = 8,
         height = 5)
}

# corr <- lapply(df2, function(x){
#   res <- corr.test(as.matrix(t(x)), method = "spearman")
# })
# 
# res <- lapply(corr, function(x){
#   x$p %>%
#     {(function(x){x[lower.tri(x)] <- NA; x})(.)} %>%
#     as.table() %>% as.data.frame() %>%
#     subset(Var1 != Var2) %>%
#     subset(!is.na(Freq)) %>%
#     subset(Freq < 0.05) %>%
#     mutate(combined=paste(Var1, Var2, sep = "."),
#            combined_backwards=paste(Var2, Var1, sep = "."))
#   
#   x$r %>%
#     as.table() %>% as.data.frame() %>%
#     mutate(combined=paste(Var1, Var2, sep = ".")) %>%
#     subset((combined  %in% a$combined | combined %in% a$combined_backwards)) %>%
#     mutate(geneSet=case_when(Var1 %in% persister_genes & Var2 %in% persister_genes ~ "Persister",
#                              Var1 %in% persister_genes & Var2 %in% relapse_genes ~ "Relapse-Persister",
#                              Var1 %in% relapse_genes & Var2 %in% persister_genes ~ "Persister-Relapse",
#                              Var1 %in% relapse_genes & Var2 %in% relapse_genes ~ "Relapse"))
#     
# })
# 
# cmb <- combn(c(unique(dat.sub$patient_id)), 2)
# 
# df_for_plot <- apply(cmb, 2, function(i){
#   x <- res[[i[1]]]
#   y <- res[[i[2]]]
#   common <- intersect(x$combined, y$combined)
#   df <- data.frame(x[x$combined %in% common,], y[y$combined %in% common,])
#   df <- df[, c("Freq", "Freq.1", "combined", "genes")]
#   return(df)
# })
# 
# 
# names(df_for_plot) <- apply(cmb, 2, function(x){
#   paste(x[1], x[2], sep = ".")
# })
# 
# 
# filtered <- lapply(df_for_plot, function(x){
#   x %>%
#     filter((Freq > 0.1 & Freq.1 > 0.1) | (Freq < -0.1 & Freq.1 < -0.1))
# })
# 
# plots <- lapply(names(filtered), function(x){
#   a <- gsub("\\..*", "", x)
#   b <- gsub(".*\\.", "", x)
#   ggplot(filtered[[x]], aes(Freq, Freq.1, color=genes)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col='red', size=1) +
#     xlab(paste("Gene pair correlations in ", a, sep = "")) +
#     ylab(paste("Gene pair correlations in ", b, sep = "")) +
#     theme_classic() +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0)
# })
# 
names(plots) <- lapply(names(df_for_plot), function(x) x)
# 
# for (i in names(plots)){
#   ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/int_cor_structure/persister_and_relapse/filtered/", i, ".pdf", sep = ""),
#          plot = plots[[i]],
#          width = 8,
#          height = 5)
# }


