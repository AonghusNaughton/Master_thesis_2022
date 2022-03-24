load("~/Proj_eng/ALL_data_Feb22.Rda")
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
source("general.r")
library(pheatmap)
library(Seurat)
library(vegan)
library(scran)
library(SC3)
library(tidyverse)
library(UCell)
library(ggplotify)
library(reshape2)
library(gplots)
library(dichromat)
library(plyr)
# Pseudo bulk approach

# Subet to only include cases that have diagnostic and d15 samples
dat.sub <- SetIdent(dat.sub, value = "pid_time")
names <- unique(dat.sub$pid_time)
names.d15 <- names[grepl(".d15", names)]
names.d0 <- gsub(".d15", ".d0", names.d15)
names.all <- c(names.d0, names.d15)
names.all <- names.all[!names.all %in% c("MRD_ALL47.d0", "MRD_ALL47.d15")] #only one cell at d15
dat.sub <- subset(dat.sub, idents = names.all[1:length(names.all)])

table(dat.sub$pid_time)

dim(dat.sub)

dat.sub <- SetIdent(dat.sub, value = "timepoint_short")

m <- FindAllMarkers(dat.sub)
m_sig_up_d15 <- m[m$p_val_adj<0.05 & m$avg_log2FC>0.25 & m$cluster=="d15",]
rownames(m_sig_up_d15) <- m_sig_up_d15$gene

d15_genes <- m_sig_up_d15$gene

dat.sub <- dat.sub[,dat.sub$timepoint_short=="d0"]

names <- unique(dat.sub$patient_id)

# Need to do first for raw counts to enable only genes that are expressed to be analysed on the scaled data.  

df1 <- lapply(names, function(x){
  df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@counts)
}) 

names(df1) <- lapply(names, function(x){
  x
})

# Remove MALAT1 as not relevant 
for (i in 1:length(names)){
  df1[[i]] <- df1[[i]][d15_genes[2:length(d15_genes)],]
}

for (i in 1:length(names)){
  df1[[i]] <- (df1[[i]][rowSums(df1[[i]])>0,])
}

gene_list_for_scaled <- lapply(df1, function(x){
  rownames(x)
})

df2 <- lapply(names, function(x){
  df <- as.data.frame(dat.sub[,dat.sub$patient_id==x]@assays$RNA@scale.data)
})

names(df2) <- lapply(names, function(x){
  x
})

# Subset count data-frame to only include genes that are expressed for each patient
for (i in 1:length(names)){
  df2[[i]] <- df2[[i]][gene_list_for_scaled[[i]],]
}


res <- lapply(df2, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

filtered_res <- list()

for (i in 1:length(names)){
  filtered_res[[i]] <- res[[i]] %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2 & abs(Freq) > 0.35) %>%
    filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2), pmin(as.character(Var1), as.character(Var2)))))) %>%
    arrange(desc(Freq))
}

names(filtered_res) <- lapply(names, function(x){
  x
})

saveRDS(filtered_res, "correlated_d15persister_genes_d0.Rds")

pdf("cor_pids_scaled.pdf", height = 30, width = 35)
lapply(res, function(x){
  pheatmap(as.matrix(x), cluster_rows = TRUE, cluster_cols = TRUE)
})
dev.off()

pdf("cor_ALL5_scaled.pdf", height = 30, width = 35)
pheatmap(as.matrix(res$ALL5), cluster_rows = T, cluster_cols = T)
dev.off()

#################################################################################################################################

dat.sub <- SetIdent(dat.sub, value = "timepoint_short")
# Find differentially expressed genes for each case with very low log-fold cutoff
marks <- lapply(unique(dat.sub$patient_id), function(x){
  m <- FindMarkers(dat.sub[,dat.sub$patient_id==x & (dat.sub$timepoint_short=="d0" | dat.sub$timepoint_short=="d15")], 
                   ident.1 = "d15", 
                   ident.2 = "d0", 
                   logfc.threshold = 0.05,
                   min.pct = 0.01)
  return(m)
})

names(marks) <- lapply(unique(dat.sub$patient_id), function(x){
  x
})

genenamesMarks <- lapply(marks, function(x){
  rownames(x)
})

# Subset markers to only include those that are common to all cases  
common_genes <-genenamesMarks %>%
  reduce(intersect)

filt_marks <- lapply(marks, function(x){
  m <- x[rownames(x) %in% common_genes,]
  return(m)
})
names(filt_marks) <- lapply(names(marks), function(x){
  x
})

# all combinations of patients and their indices
cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
cmb1 <- combn(unique(dat.sub$patient_id), 2)

# create dataframes that have stricter log-fold cut-offs
m_up_d15_high <- list()
for (i in 1:length(unique(dat.sub$patient_id))){
  m_up_d15_high[[i]] <- filt_marks[[i]][filt_marks[[i]]$avg_log2FC>0.5 & (filt_marks[[i]]$pct.1>0.25 | filt_marks[[i]]$pct.2>0.25),]
}

names(m_up_d15_high) <- lapply(unique(dat.sub$patient_id), function(x){
  x
})

#################################################################################################################################

# List of all commonly up-regulated genes at d15 between all combinations of patients 
combinations_high <- apply(cmb, 2, function(x){
  intersect(rownames(m_up_d15_high[[x[1]]]),  rownames(m_up_d15_high[[x[2]]]))
})

names(combinations_high) <- apply(cmb1, 2, function(x){
  paste(x[1], "&", x[2], sep = "")
})

#################################################################################################################################

# Must have a high logfold change in at least 3 cases 

persister_genes_high <- c()
for (i in 1:length(names(combinations_high))){
  for (j in 1:length(combinations_high[[i]])){
    for (k in 1:length(names(combinations_high))){
      if (k != i) {
        if (combinations_high[[i]][j] %in% combinations_high[[k]]){
          persister_genes_high <- c(persister_genes_high, combinations_high[[i]][j])
        }
      }
    }
  }
}

persister_genes_high <- unique(persister_genes_high)


# Must have a high log fold change in at least two cases. 
# Append to list if gene is present in any pair of patients. Option to increase this threshold to 3 patients above (commented code)
# persister_genes_high <- c()
# 
# for (i in 1:length(names(combinations_high))){
#   persister_genes_high <- c(persister_genes_high, combinations_high[[i]])
# }
# 
# persister_genes_high <- unique(persister_genes_high)

lapply(filt_marks, function(x){
  table(persister_genes_high %in% rownames(x))
})
#################################################################################################################################

# Arrange dataframes and extract vectors of p-value. 

filt_marks <- lapply(filt_marks, function(x){
  x[persister_genes_high,]
})

filt_marks <- lapply(filt_marks, function(x){
  x[order(rownames(x)),]
})

p_val_vecs <- lapply(filt_marks, function(x){
  x$p_val_adj
})


combined_pvals <- metapod::combineParallelPValues(list(p_val_vecs$MRD_ALL71, 
                                 p_val_vecs$MRD_ALL64,
                                 p_val_vecs$MRD_ALL67,
                                 p_val_vecs$MRD_ALL68,
                                 p_val_vecs$ALL5), 
                                 method = "fisher", 
                                 log.p = F)

table(combined_pvals$p.value < 0.05)

combined_pvals_bool <- combined_pvals$p.value<0.05
indices_sig <- which(combined_pvals_bool %in% 1)

filt_sig_marks <- lapply(filt_marks, function(x){
  x[indices_sig,]
})

true_persister_genes_d15_new <- lapply(filt_sig_marks, function(x){
  rownames(x)
})
true_persister_gene_set_d15_new <- true_persister_genes_d15_new[[1]]
saveRDS(true_persister_gene_set_d15_new, "true_persister_gene_set_d15_new1.Rds")
#################################################################################################################################
# Expression heatmaps 

DoHeatmap(dat.sub[, dat.sub$timepoint_short=="d15"], features = true_persister_gene_set_d15, group.by = "patient_id")
exp <- DotPlot(dat.sub, features = true_persister_gene_set_d15, group.by = "patient_id", split.by = "timepoint")
View(exp$data)



#################################################################################################################################
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
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
  df1[[i]] <- df1[[i]][true_persister_gene_set_d15,]
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


res_persister <- lapply(df2, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

saveRDS(res_persister, "res_persister.Rds")

filtered_res <- list()

for (i in names){
  filtered_res[[i]] <- res_persister[[i]] %>%
    as.table() %>% as.data.frame() %>%
    subset(Var1 != Var2 & abs(Freq) > 0.35) %>%
    filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2), pmin(as.character(Var1), as.character(Var2)))))) %>%
    arrange(desc(Freq))
}

names(filtered_res) <- lapply(names, function(x){
  x
})
#################################################################################################################################

# Create pseudo gene set that represents persister gene set (similar pct expression)

exp <- FetchData(dat.sub[,dat.sub$timepoint_short=="d15"], true_persister_gene_set_d15)
average_pct.exp_persister <- colMeans(as.matrix(colMeans(exp  > 0))*100)

exp.vals_for_pseudo <- FetchData(dat.sub[,dat.sub$timepoint_short=="d15"], 
                                 rownames(dat.sub[,dat.sub$timepoint_short=="d15"])[rownames(dat.sub[,dat.sub$timepoint_short=="d15"]) %!in% 
                                                                                      true_persister_gene_set_d15])
pct.exp <- as.matrix(colMeans(exp.vals_for_pseudo  > 0))*100
pct.exp_filt <- pct.exp[pct.exp>51,]
random_indices <- floor(runif(length(true_persister_gene_set_d15), min = 0, max = length(pct.exp_filt) + 1))
pseudo_persister <- pct.exp_filt[random_indices]

pseudo_persister <- names(pseudo_persister)

saveRDS(pseudo_persister, "pseudo_persister_genes.Rds")



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
  df1[[i]] <- df1[[i]][pseudo_persister,]
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


res_pseudo_persister <- lapply(df2, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

saveRDS(res_pseudo_persister, "res_pseudo_persister.Rds")

#################################################################################################################################
paletteLength <- 100
plots <- lapply(res_persister, function(x){
  pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T, 
           color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
           breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
})

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/correlations_3pids//persister_genes/persister_cor_new_", x, ".pdf", sep = ""), 
         plot = plots[[x]],
         width = 15,
         height = 12) 
})


#################################################################################################################################
# Which persister genes are co-expressed at diagnosis on a case-by-case basis?

co.expressed_persisters <- lapply(filtered_res, function(x){
  geneset <- c(unique(x$Var1), unique(x$Var2))
}) 

#################################################################################################################################

d0.rowSums <- lapply(res, function(x){
  sort(rowSums(x), decreasing = T)
})

df <- lapply(d0.rowSums, function(x){
  dat.frame <- data.frame(x)
})

rownames <- lapply(df, function(x){
  rownames(x)
})

common <- Reduce(intersect, rownames)


for (i in 1:length(df)){
  df[[i]] <- df[[i]][common, ]
}

df1 <- lapply(df, function(x){
  

})

df2 <- as.data.frame(do.call(cbind, df))


#################################################################################################################################

# Module scores for diagnostic samples (to each patient individually and then as cohort) - UCell 

dat.sub <- dat[, (dat$dna_cell_type=="blasts" | is.na(dat$dna_cell_type)) & 
                 dat$rna_cell_type=="blasts" & 
                 dat$patient_id=="ALL3" & 
                 dat$cell_phase=="G1"]

split1 <- split(sc3_split$ALL3.d0$cell_id, sc3_split$ALL3.d0$sc3_2_clusters)
# split2 <- split(sc3_ALL71_d0_split_wNa$cell_id, sc3_ALL71_d0_split_wNa$sc3_4_clusters)
split1.rel <- split(sc3_split_wNa$ALL3.rel$cell_id, sc3_split_wNa$ALL3.rel$sc3_5_clusters)
split1.rel2 <- split(sc3_split$ALL3.rel2$cell_id, sc3_split$ALL3.rel2$sc3_5_clusters)

annot <- dat.sub@meta.data %>%
  mutate(sc3_clusters_to_compare=case_when(cell_id %in% split1$`1` ~ 1,
                                           cell_id %in% split1$`2` ~ 2,
                                           cell_id %in% split1.rel$`1` ~ 3,
                                           cell_id %in% split1.rel$`2` ~ 4,
                                           cell_id %in% split1.rel$`3` ~ 5,
                                           cell_id %in% split1.rel$`4` ~ 6,
                                           cell_id %in% split1.rel$`5` ~ 7,
                                           cell_id %in% split1.rel2$`1` ~ 8,
                                           cell_id %in% split1.rel2$`2`~ 9,
                                           cell_id %in% split1.rel2$`3`~ 10,
                                           cell_id %in% split1.rel2$`4`~ 11,
                                           cell_id %in% split1.rel2$`5`~ 12))
dat.sub <- AddMetaData(dat.sub, annot)
dat.sub <- AddModuleScore_UCell(dat.sub, features = list(true_persister_gene_set_d15))

colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(list(true_persister_gene_set_d15))+1):ncol(dat.sub@meta.data)] <- "d15_persister_gene_score"

dat.sub <- dat.sub[, dat.sub$sc3_clusters_to_compare %in% c(1:12)]

d15.score <- tibble(dat.sub$d15_persister_gene_score, dat.sub$sc3_clusters_to_compare, dat.sub$timepoint, dat.sub$dna_best_class)
bcell.score <- tibble(dat.sub$bcell_score, dat.sub$sc3_clusters_to_compare, dat.sub$timepoint, dat.sub$dna_best_class)

plt_d15.scores.sc3 <- d15.score %>%
  ggplot(aes(x=dat.sub$sc3_clusters_to_compare, 
             y=dat.sub$d15_persister_gene_score, 
             fill=as.factor(dat.sub$sc3_clusters_to_compare),
             colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "SC3 cluster") +
  scale_color_discrete(name = "Timepoint") +
  ylab("Persister score") + xlab(NULL) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0.1, 0.28)

plt_d15.scores.clone <- d15.score %>%
  ggplot(aes(x=dat.sub$dna_best_class, 
             y=dat.sub$d15_persister_gene_score, 
             fill=as.factor(dat.sub$dna_best_class),
             colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Clone") +
  scale_color_discrete(name = "Timepoint") +
  ylab("Persister score") + xlab(NULL) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  ylim(0.1, 0.28)

plt_bcell.scores.sc3 <- bcell.score %>%
  ggplot(aes(x=dat.sub$sc3_clusters_to_compare, 
             y=dat.sub$bcell_score, 
             fill=as.factor(dat.sub$sc3_clusters_to_compare),
             colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "SC3 cluster") +
  scale_color_discrete(name = "Timepoint") +
  ylab("B-cell score") + xlab(NULL) +
  theme(legend.position = "right")

plt_bcell.scores.clone <- bcell.score %>%
  ggplot(aes(x=dat.sub$dna_best_class, 
             y=dat.sub$bcell_score, 
             fill=as.factor(dat.sub$dna_best_class),
             colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Clone") +
  scale_color_discrete(name = "Timepoint") +
  ylab("B-cell score") + xlab(NULL) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# this is a test 
