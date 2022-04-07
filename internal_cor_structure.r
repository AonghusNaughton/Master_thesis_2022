# Pairwise gene correlations in each patient 
library(tidyverse)
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
dat.sub <- dat.sub[, dat.sub$timepoint_short=="d0"]
names <- unique(dat.sub$patient_id)

persister_genes <- readRDS("persister_genes.Rds")
relapse_genes <- readRDS("relapse_genes.Rds")
genes <- c(persister_genes, relapse_genes)

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

res <- lapply(df2, function(x){
  t(x) %>%
    as.matrix %>%
    cor %>% {(function(x){x[upper.tri(x)] <- NA; x})(.)} %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2) %>%
    subset(!is.na(Freq)) %>%
    mutate(combined=paste(Var1, Var2, sep = ".")) %>% 
    mutate(genes=case_when(Var1 %in% persister_genes & Var2 %in% persister_genes ~ "Persister pair",
                             (Var1 %in% persister_genes & Var2 %in% relapse_genes | 
                                Var1 %in% relapse_genes & Var2 %in% persister_genes) ~ "Mix pair",
                             Var1 %in% relapse_genes & Var2 %in% relapse_genes ~ "Relapse pair")) %>%
    select(-c(Var1, Var2)) %>%
    arrange(combined) 
})

cmb <- combn(c(unique(dat.sub$patient_id)), 2)

df_for_plot <- apply(cmb, 2, function(i){
  x <- res[[i[1]]]
  y <- res[[i[2]]]
  common <- intersect(x$combined, y$combined)
  df <- data.frame(x[x$combined %in% common,], y[y$combined %in% common,])
  df <- df[, c("Freq", "Freq.1", "combined", "genes")]
  return(df)
})


names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(x[1], x[2], sep = ".")
})


filtered <- lapply(df_for_plot, function(x){
  x %>%
    filter((Freq > 0.1 & Freq.1 > 0.1) | (Freq < -0.1 & Freq.1 < -0.1))
})

plots <- lapply(names(filtered), function(x){
  a <- gsub("\\..*", "", x)
  b <- gsub(".*\\.", "", x)
  ggplot(filtered[[x]], aes(Freq, Freq.1, color=genes)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, col='red', size=1) +
    xlab(paste("Gene pair correlations in ", a, sep = "")) +
    ylab(paste("Gene pair correlations in ", b, sep = "")) +
    theme_classic() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
})

names(plots) <- lapply(names(test), function(x) x)

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/int_cor_structure/persister_and_relapse/filtered/", i, ".pdf", sep = ""),
         plot = plots[[i]],
         width = 8,
         height = 5)
}


