library(tidyverse)

dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
relapse_genes <- readRDS("relapse_genes.Rds")
persister_genes <- readRDS("persister_genes.Rds")

geneSets <- unique(c(persister_genes, relapse_genes))
# geneSets <- unique(c(relapse_genes, random_relapse))
timepoint <- "d0"
dat.sub <- dat.sub[,dat.sub$cell_phase=="G1" & dat.sub$timepoint_short==timepoint]

names <- unique(dat.sub$patient_id)
if (timepoint == "d15"){
  names <- names[!names %in% "MRD_ALL47"]
}

minCount <- 10

names <- unique(dat.sub$patient_id)
df1 <- lapply(names, function(x) as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)) 
names(df1) <- lapply(names, function(x) x)

df_na <- lapply(df1, function(x){
  na_if(x, 0)
})

genes_above_threshold <- lapply(df_na, function(x){
  genes <- apply(x, 1, function(row) sum(!is.na(row)))
  genes <- names(genes[genes > minCount])
  return(genes)
})

for (i in names){
  df1[[i]] <- df1[[i]][genes_above_threshold[[i]],]
}

cmb1 <- combn(unique(dat.sub$patient_id), 2)

common_to_two <- apply(cmb1, 2, function(x){
  intersect(rownames(df1[[x[1]]]), rownames(df1[[x[2]]]))
})

names(common_to_two) <- apply(cmb1, 2, function(x){
  paste(x[1], "&", x[2], sep = "")
})

genes_in_at_least_two <- c()
for (i in 1:length(names(common_to_two))){
  genes_in_at_least_two <- c(genes_in_at_least_two, common_to_two[[i]])
}

genes_in_at_least_two <- unique(genes_in_at_least_two)
# saveRDS(genes_in_at_least_two, "highly_exp_genes.Rds")

dat.sub <- subset(dat.sub, features= genes_in_at_least_two)
saveRDS(dat.sub, "dat.sub_filtered_10_cells_d0.Rds")

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

# saveRDS(df1, "counts_per_pid_10.Rds")

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

# saveRDS(df2, "scale_data_counts_10.Rds")

for (i in names){
  df2[[i]] <- df2[[i]][gene_list_for_scaled[[i]],]
}

corr <- lapply(df2, function(x){
  res <- corr.test(as.matrix(t(x)), method = "spearman")
})

corr_filtered <- lapply(corr, function(matrix){
  # a <- matrix$p %>%
  #   {(function(x){x[lower.tri(x)] <- NA; x})(.)} %>%
  #   as.table() %>% as.data.frame() %>%
  #   subset(Var1 != Var2) %>%
  #   subset(!is.na(Freq)) %>%
  #   subset(Freq < 0.1) %>%
  #   mutate(combined=paste(Var1, Var2, sep = "."),
  #          combined_backwards=paste(Var2, Var1, sep = "."))
  
  b <- matrix$r %>%
    as.table() %>% as.data.frame() %>%
    mutate(combined=paste(Var1, Var2, sep = ".")) %>%
    # subset((combined  %in% a$combined | combined %in% a$combined_backwards)) %>%
    mutate(geneSet=case_when(Var1 %in% persister_genes & Var2 %in% persister_genes ~ "Persister",
                             Var1 %in% persister_genes & Var2 %in% relapse_genes ~ "Relapse-Persister",
                             Var1 %in% relapse_genes & Var2 %in% persister_genes ~ "Persister-Relapse",
                             Var1 %in% relapse_genes & Var2 %in% relapse_genes ~ "Relapse")) %>%
  group_by(geneSet) %>%
  group_by(Var2, .add=T) %>%
  summarise(avg=mean(Freq)) %>%
    mutate(x_axis=case_when(geneSet %in% "Relapse" ~ avg,  geneSet %in% "Persister-Relapse" ~ avg)) %>%
    mutate(y_axis=case_when(geneSet %in% "Persister" ~ avg, geneSet %in% "Relapse-Persister" ~ avg)) %>%
    mutate(gene_set=case_when(Var2 %in% persister_genes ~ "Persister",
                              Var2 %in% relapse_genes ~ "Relapse")) %>%
    ungroup(Var2) %>%
    ungroup(geneSet) %>%
    dplyr::select(-c(geneSet, avg)) %>%
    replace_na(list(x_axis=0, y_axis=0)) %>%
    group_by(Var2, gene_set) %>%
    summarise(x_axis=sum(x_axis),
              y_axis=sum(y_axis))
  return(b) 
})

plots <- lapply(corr_filtered, function(x){
  x %>%
    ggplot(aes(x_axis,y_axis, color=gene_set)) +
    geom_point() + 
    xlab(paste0("Relapse")) +
    ylab(paste0("Persister")) +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_classic()
})

lapply(names(plots), function(i){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/Average_correlations_new/filtered_p_vals/", i, ".pdf", sep = ""),
         plot = plots[[i]],
         width = 8,
         height = 5)
})


# saveRDS(res, "correlation_matrix.Rds")



######################################################################################################################

# compute similarity between persister matrices 

# 11 correlation matrices 
# # cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
# names_cmb <- combn(unique(dat.sub$patient_id), 2)
# 
# matched_res <- apply(cmb, 2, function(x){
#   df <- list(m1=res_list$res_Persister[[x[[1]]]][intersect(rownames(res_list$res_Persister[[x[[1]]]]), 
#                                                    rownames(res_list$res_Persister[[x[[2]]]])), 
#                                                  intersect(colnames(res_list$res_Persister[[x[[1]]]]),
#                                                  colnames(res_list$res_Persister[[x[[2]]]]))],
#              m2=res_list$res_Persister[[x[[2]]]][intersect(rownames(res_list$res_Persister[[x[[1]]]]), 
#                                                         rownames(res_list$res_Persister[[x[[2]]]])),
#                                                  intersect(colnames(res_list$res_Persister[[x[[1]]]]),
#                                                  colnames(res_list$res_Persister[[x[[2]]]]))])
# })
# 
# names(matched_res) <- apply(names_cmb, 2, function(x){
#   paste(x[[1]], x[[2]], sep = ".")
# })
# 
# lapply(matched_res, function(x){
#   print(paste(length(rownames(x[["m1"]])), length(rownames(x[["m2"]])), sep = " "))
#   print(paste(length(colnames(x[["m1"]])), length(colnames(x[["m2"]])), sep = " "))
# })
# 
# extracted_tri.all.combos <- lapply(matched_res, function(x){
#   list <- list(t(x[["m1"]])[lower.tri(t(x[["m1"]]))],
#                t(x[["m2"]])[lower.tri(t(x[["m2"]]))])
#   names(list) <- c("m1", "m2")
#   return(list)
# })
# 
# lapply(extracted_tri.all.combos, function(x){
#   print(paste(length(x[["m1"]]), length(x[["m2"]]), sep = " "))
# })
# 
# spearman_correlations <- lapply(extracted_tri.all.combos, function(x){
#   cor.test(x[["m1"]], x[["m2"]], method = "spearman", exact=F)
# })

