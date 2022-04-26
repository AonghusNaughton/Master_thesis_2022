# Pairwise gene correlations in each patient 
library(tidyverse)
library(Hmisc)
library(ggrepel)
dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
names <- unique(dat.sub$patient_id)

persister_genes <- readRDS("persister_genes.Rds")

df <- lapply(names, function(x){
  if (x=="ALL3"){
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts"]@assays$RNA@counts)
    df2 <- as.data.frame(dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts"]@assays$RNA@scale.data)
  } else {
    df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)
    df2 <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@scale.data)
  }
  df <- df[persister_genes,]
  df <- df[rowSums(df)>0,]
  genes <- rownames(df)
  df2 <- df2[genes,]
  return(df2)
})

names(df) <- lapply(names, function(x) x)


var.corr <- lapply(df, function(x){
  res <- correlatePairs(x)
})

sig.corr <- lapply(var.corr, function(x){
  x %>%
    as.data.frame() %>%
    filter(p.value <= 0.05)
})

sig.corr <- lapply(sig.corr, function(x){
  x %>%
    mutate(combined=paste(gene1, gene2, sep = " & "))
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
  df <- df[, c("rho", "rho.1", "combined")]
  df <- df %>%
    mutate(high_or_low=case_when(rho > 0.15 & rho.1 > 0.15 ~ "high",
                                 rho < -0.15 & rho.1 < -0.15 ~ "low",
                                 ((rho < 0.15 & rho.1 > -0.15) | (rho > 0.15 & rho.1 > -0.15) | (rho < 0.15 & rho.1 < -0.15))  ~ "neither"))
  
  return(df)
})

names(df_for_plot) <- apply(cmb, 2, function(x){
  paste(x[1], x[2], sep = ".")
})

combinations_of_combos <- combn(names(df_for_plot), 2)


temp_df <- apply(combinations_of_combos, 2, function(i){
  x <- df_for_plot[[i[1]]]
  y <- df_for_plot[[i[2]]]
  high_x <- x %>%
    filter(high_or_low == "high") %>%
    pull(combined)
  low_x <- x %>%
    filter(high_or_low == "low") %>%
    pull(combined)
  high_y <- y %>%
    filter(high_or_low == "high") %>%
    pull(combined)
  low_y <- y %>%
    filter(high_or_low == "low") %>%
    pull(combined)
  common_high <- intersect(high_x, high_y)
  common_low <- intersect(low_x, low_y)
  return(list(common_high,
              common_low))
})

names(temp_df) <- apply(combinations_of_combos, 2, function(x){
  paste(x[1], x[2], sep = ".")
})

common_high <- c()
common_low <- c()
for(i in temp_df){
  common_high <- unique(c(common_high, i[[1]])) 
  common_low <- unique(c(common_low, i[[2]]))
}

df_for_plot <- lapply(df_for_plot, function(x){
  x %>%
    mutate(commonality=case_when(combined %in% common_high ~ "Commonly correlated",
                                 combined %in% common_low ~ "Commonly anti-correlated",
                                 (!combined %in% common_high & !combined %in% common_low) ~ "Not"))
})

cols <- c("Commonly correlated" = "darkgreen", "Commonly anti-correlated" = "red", "Not" = "blue")

plots <- lapply(names(df_for_plot), function(x){
    a <- gsub("\\..*", "", x)
    b <- gsub(".*\\.", "", x)
    ggplot(df_for_plot[[x]], aes(rho, rho.1, color=commonality, label=combined)) +
      geom_point() +
      # geom_smooth(method = "lm", se = FALSE, col='red', size=1) +
      xlab(paste("Gene pair correlations in ", a, sep = "")) +
      ylab(paste("Gene pair correlations in ", b, sep = "")) +
      theme_classic() +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      theme(legend.position = "right",
            legend.title=element_blank()) +
      geom_label_repel(aes(label = ifelse((((commonality=="Commonly correlated" | commonality=="Commonly anti-correlated") & ((rho > 0.2 & rho.1 > 0.2) | (rho < -0.2 & rho.1 < -0.2))) | ((rho > 0.25 & rho.1 > 0.25) | (rho < -0.25 & rho.1 < -0.25))), as.character(combined), "")),
                       box.padding   = 0.1,
                       segment.color = 'grey50',
                       max.overlaps = Inf,
                       size=1,
                       show.legend = F) +
      scale_colour_manual(values = cols)
  })

names(plots) <- lapply(names(df_for_plot), function(x) x)



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

# 
# for (i in names(plots)){
#   ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/int_cor_structure/persister_and_relapse/filtered/", i, ".pdf", sep = ""),
#          plot = plots[[i]],
#          width = 8,
#          height = 5)
# }


