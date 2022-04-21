library(tidyverse)
highly_exp_genes <- readRDS("highly_exp_genes.Rds")
dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
counts <- readRDS("counts_per_pid_10.Rds")
scaled_counts <- readRDS("scale_data_counts_10.Rds")
persister_genes <- readRDS("persister_genes.Rds")
names <- unique(dat.sub$patient_id)

counts <- lapply(counts, function(x){
  x <- x[rowSums(x)>0,]
})

gene_list_for_scaled <- lapply(counts, function(x){
  rownames(x)
})

for (i in names){
  scaled_counts[[i]] <- scaled_counts[[i]][gene_list_for_scaled[[i]],]
}

correlations <- lapply(scaled_counts, function(x){
  res <- correlatePairs(x, pairings = list(rownames(x)[!rownames(x) %in% persister_genes],
                                           rownames(x)[rownames(x) %in% persister_genes]))
})

sig.corr <- lapply(correlations, function(x){
  x %>%
    as.data.frame() %>%
    filter(FDR <= 0.05) %>%
    filter(FDR != 0) %>%
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
  df <- df[, c("rho", "rho.1", "combined", "gene1", "gene2")]
  return(df)
})

names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(x[1], x[2], sep = ".")
})

plots <- lapply(names(df_for_plot), function(x){
  a <- gsub("\\..*", "", x)
  b <- gsub(".*\\.", "", x)
  ggplot(df_for_plot[[x]], aes(rho, rho.1, label=c(gene1, gene2))) +
    geom_point() +
    geom_text(size=3, check_overlap = TRUE, aes(color=)) +
    geom_smooth(method = "lm", se = FALSE, col='red', size=1) +
    xlab(paste("Gene pair correlations in ", a, sep = "")) +
    ylab(paste("Gene pair correlations in ", b, sep = "")) +
    theme_classic() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
})
names(plots) <- lapply(names(df_for_plot), function(x) x)

for (i in names(plots)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/int_cor_structure/persister_and_relapse/AllvsPersister/", i, ".pdf", sep = ""),
         plot = plots[[i]],
         width = 8,
         height = 5)
}

# correlations <- lapply(scaled_counts, function(x){
#   res <- cor(x=as.matrix(t(x[rownames(x)[!rownames(x) %in% persister_genes],])),
#              y=as.matrix(t(x[rownames(x)[rownames(x) %in% persister_genes],])),
#              method = "spearman")
# })
# 
# saveRDS(correlations, "correlations_allVsPersister.Rds")
# 
# average_cors <- lapply(names, function(x){
#   y <- data.frame(rowMeans(correlations[[x]]))
#   colnames(y) <- paste0(x)
#   return(y)
# })
# 
# names(average_cors) <- lapply(names, function(x) x)
# 
# average_cors <- lapply(names(average_cors), function(x){
#   name <- paste0(x)
#   absent <- setdiff(highly_exp_genes, rownames(average_cors[[x]]))
#   new_rows <- lapply(absent, function(y){
#     df <- data.frame(y, 0)
#     rownames(df) <- df[[1]]
#     df[[1]] <- NULL
#     colnames(df) <- c(name)
#     return(df)
#   })
#   new_rows <- do.call(rbind, new_rows)
#   cors <- rbind(average_cors[[x]], new_rows)
#   return(cors)
# })
# 
# names(average_cors) <- lapply(names, function(x) x)
# 
# average_cors_df <- do.call(cbind, average_cors)
# 
# MRD_cases <- readRDS("persister_cases.Rds")
# Relapse_cases <- readRDS("relapsed_cases.Rds")
# 
# tidy <- average_cors_df %>%
#   as.matrix %>% as.table() %>% as.data.frame() %>%
#   rename(Gene=Var1, Patient=Var2, Average_correlation=Freq) %>%
#   mutate(MRD_or_Relapse=case_when(Patient %in% MRD_cases ~ "MRD",
#                                   Patient %in% Relapse_cases ~ "Relapse")) %>%
#   group_by(Gene) %>%
#   summarize(Mean_average_exp = mean(Average_correlation, na.rm=TRUE)) %>%
#   arrange(desc(Mean_average_exp))
#   
# ranks <- tidy$Mean_average_exp
# names(ranks) <- tidy$Gene
# ranks <- sort(ranks, decreasing = T)
# 
# gse <- gseGO(geneList=ranks, 
#              ont ="BP", 
#              keyType = "SYMBOL", 
#              nPermSimple = 100000, 
#              minGSSize = 3, 
#              maxGSSize = 800, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = org.Hs.eg.db, 
#              pAdjustMethod = "none")
# 
# pdf("test.pdf",height = 12, width = 10)
# dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=55))
# dev.off()
# 
# 
# 
# # Gene pairs 
# gene_pairs <- lapply(correlations, function(x){
#   x %>%
#     as.table() %>% as.data.frame() %>%
#     mutate(combined=paste(Var1, Var2, sep = "."))
# })




  