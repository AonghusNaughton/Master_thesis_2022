load("ALL_data_Feb22.Rda")
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
# Pseudo bulk approach

# Subet to only include cases that have diagnostic and d15 samples
dat.sub <- SetIdent(dat.sub, value = "pid_time")
names <- unique(dat.sub$pid_time)
names.rel <- names[grepl(".rel", names)]
names.d0 <- gsub(".rel", ".d0", names.rel)
names.d0 <- unique(gsub("2", "", names.d0))
names.all <- c(names.d0, names.rel)
# names.all <- names.all[!names.all %in% c("MRD_ALL47.d0", "MRD_ALL47.d15")] #only one cell at d15
dat.sub <- subset(dat.sub, idents = names.all[1:length(names.all)])
to_remove <- dat.sub[,dat.sub$patient_id=="ALL3" & is.na(dat.sub$dna_cell_type)]$cell_id
dat.sub <- dat.sub[, !colnames(dat.sub) %in% to_remove]
table(dat.sub$pid_time)

dim(dat.sub)


dat.sub <- SetIdent(dat.sub, value = "timepoint_short")

table(dat.sub[,(dat.sub$timepoint_short=="d0" | dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")]$pid_time)

# Find differentially expressed genes for each case with very low log-fold cutoff
marks_rel <- lapply(unique(dat.sub$patient_id), function(x){
  m <- FindMarkers(dat.sub[,dat.sub$patient_id==x & 
                             (dat.sub$timepoint_short=="d0" | dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")], 
                   if (x == "ALL3"){
                     ident.1 = c("rel", "rel2") } else {
                       ident.1 = "rel"
                     }, 
                   ident.2 = "d0", 
                   logfc.threshold = 0.05,
                   min.pct = 0.01)
  return(m)
})

names(marks_rel) <- lapply(unique(dat.sub$patient_id), function(x){
  x
})

genenamesMarks <- lapply(marks_rel, function(x){
  rownames(x)
})

# Subset markers to only include those that are common to all cases  
common_genes <-genenamesMarks %>%
  reduce(intersect)

filt_marks <- lapply(marks_rel, function(x){
  m <- x[rownames(x) %in% common_genes,]
  return(m)
})
names(filt_marks) <- lapply(names(marks_rel), function(x){
  x
})

# all combinations of patients and their indices
cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
cmb1 <- combn(unique(dat.sub$patient_id), 2)

# create dataframes that have stricter log-fold cut-offs
m_up_rel_high <- list()
for (i in 1:length(unique(dat.sub$patient_id))){
  m_up_rel_high[[i]] <- filt_marks[[i]][filt_marks[[i]]$avg_log2FC>0.5 & (filt_marks[[i]]$pct.1>0.25 | filt_marks[[i]]$pct.2>0.25),]
}

names(m_up_rel_high) <- lapply(unique(dat.sub$patient_id), function(x){
  x
})

#################################################################################################################################

# List of all commonly up-regulated genes at d15 between all combinations of patients 
combinations_high_rel <- apply(cmb, 2, function(x){
  intersect(rownames(m_up_rel_high[[x[1]]]),  rownames(m_up_rel_high[[x[2]]]))
})

names(combinations_high_rel) <- apply(cmb1, 2, function(x){
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
# persister_genes_high_rel <- c()
# 
# for (i in 1:length(names(combinations_high_rel))){
#   persister_genes_high_rel <- c(persister_genes_high_rel, combinations_high_rel[[i]])
# }
# 
# persister_genes_high_rel <- unique(persister_genes_high_rel)

lapply(filt_marks_rel, function(x){
  table(persister_genes_high_rel %in% rownames(x))
})
#################################################################################################################################

# Arrange dataframes and extract vectors of p-value. 

filt_marks_rel <- lapply(filt_marks_rel, function(x){
  x[persister_genes_high_rel,]
})

filt_marks_rel <- lapply(filt_marks_rel, function(x){
  x[order(rownames(x)),]
})

p_val_vecs_rel <- lapply(filt_marks_rel, function(x){
  x$p_val_adj
})


combined_pvals_rel <- metapod::combineParallelPValues(list(p_val_vecs_rel[[1]], 
                                                       p_val_vecs_rel[[2]],
                                                       p_val_vecs_rel[[3]],
                                                       p_val_vecs_rel[[4]]), 
                                                  method = "fisher", 
                                                  log.p = F)

table(combined_pvals_rel$p.value < 0.05)

combined_pvals_bool <- combined_pvals_rel$p.value<0.05
indices_sig <- which(combined_pvals_bool %in% 1)

filt_sig_marks_rel <- lapply(filt_marks_rel, function(x){
  x[indices_sig,]
})

true_persister_genes_rel <- lapply(filt_sig_marks_rel, function(x){
  rownames(x)
})
true_persister_gene_set_rel_new <- true_persister_genes_rel[[1]]
saveRDS(true_persister_gene_set_rel_new, "true_persister_gene_set_rel_new1.Rds")

#################################################################################################################################
# Expression heatmaps 

DoHeatmap(dat.sub[, (dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")], features = true_persister_gene_set_rel_new, group.by = "patient_id")

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
  df1[[i]] <- df1[[i]][true_persister_gene_set_rel,]
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


res_relapse <- lapply(df2, function(x){
  res <- cor(as.matrix(t(x)), method = "spearman")
})

saveRDS(res_relapse, "res_relapse.Rds")


#################################################################################################################################
paletteLength <- 100
plots <- lapply(res_relapse, function(x){
  pheatmap(as.matrix(x), cluster_rows = T, cluster_cols = T, 
           color = colorRampPalette(c("yellow", "white", "blue"))(paletteLength),
           breaks = c(seq(min(as.matrix(x)), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(as.matrix(x))/paletteLength, max(as.matrix(x)), length.out=floor(paletteLength/2))))
})

lapply(names(plots), function(x){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/March22/correlations_3pids/relapse_genes/relapse_cor_new_", x, ".pdf", sep = ""), 
         plot = plots[[x]],
         width = 12,
         height = 10) 
})

#################################################################################################################################

# 
# 
# dat.sub <- dat[, (dat$dna_cell_type=="blasts" | is.na(dat$dna_cell_type)) & 
#                  dat$rna_cell_type=="blasts" & 
#                  dat$patient_id=="ALL3" & 
#                  dat$cell_phase=="G1"]
# 
# split1 <- split(sc3_split$ALL3.d0$cell_id, sc3_split$ALL3.d0$sc3_2_clusters)
# # split2 <- split(sc3_ALL71_d0_split_wNa$cell_id, sc3_ALL71_d0_split_wNa$sc3_4_clusters)
# split1.rel <- split(sc3_split_wNa$ALL3.rel$cell_id, sc3_split_wNa$ALL3.rel$sc3_5_clusters)
# split1.rel2 <- split(sc3_split$ALL3.rel2$cell_id, sc3_split$ALL3.rel2$sc3_5_clusters)
# 
# annot <- dat.sub@meta.data %>%
#   mutate(sc3_clusters_to_compare=case_when(cell_id %in% split1$`1` ~ 1,
#                                            cell_id %in% split1$`2` ~ 2,
#                                            cell_id %in% split1.rel$`1` ~ 3,
#                                            cell_id %in% split1.rel$`2` ~ 4,
#                                            cell_id %in% split1.rel$`3` ~ 5,
#                                            cell_id %in% split1.rel$`4` ~ 6,
#                                            cell_id %in% split1.rel$`5` ~ 7,
#                                            cell_id %in% split1.rel2$`1` ~ 8,
#                                            cell_id %in% split1.rel2$`2`~ 9,
#                                            cell_id %in% split1.rel2$`3`~ 10,
#                                            cell_id %in% split1.rel2$`4`~ 11,
#                                            cell_id %in% split1.rel2$`5`~ 12))
# dat.sub <- AddMetaData(dat.sub, annot)
# dat.sub <- AddModuleScore_UCell(dat.sub, features = list(true_persister_gene_set_rel))
# 
# colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(list(true_persister_gene_set_rel))+1):ncol(dat.sub@meta.data)] <- "rel_score"
# 
# dat.sub <- dat.sub[, dat.sub$sc3_clusters_to_compare %in% c(1:12)]
# 
# rel.score <- tibble(dat.sub$rel_score, dat.sub$sc3_clusters_to_compare, dat.sub$timepoint, dat.sub$dna_best_class)
# bcell.score <- tibble(dat.sub$bcell_score, dat.sub$sc3_clusters_to_compare, dat.sub$timepoint, dat.sub$dna_best_class)
# 
# plt_rel.scores.sc3 <- rel.score %>%
#   ggplot(aes(x=dat.sub$sc3_clusters_to_compare, 
#              y=dat.sub$rel_score, 
#              fill=as.factor(dat.sub$sc3_clusters_to_compare),
#              colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "SC3 cluster") +
#   scale_color_discrete(name = "Timepoint") +
#   ylab("Relapse score") + xlab(NULL) +
#   theme(legend.position = "right",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) + 
#   ylim(0.2, 0.35)
# 
# plt_rel.scores.clone <- rel.score %>%
#   ggplot(aes(x=dat.sub$dna_best_class, 
#              y=dat.sub$rel_score, 
#              fill=as.factor(dat.sub$dna_best_class),
#              colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "Clone") +
#   scale_color_discrete(name = "Timepoint") +
#   ylab("Relapse score") + xlab(NULL) +
#   theme(legend.position = "right",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) + 
#   ylim(0.2, 0.35)
# 
# plt_bcell.scores.sc3 <- bcell.score %>%
#   ggplot(aes(x=dat.sub$sc3_clusters_to_compare, 
#              y=dat.sub$bcell_score, 
#              fill=as.factor(dat.sub$sc3_clusters_to_compare),
#              colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "SC3 cluster") +
#   scale_color_discrete(name = "Timepoint") +
#   ylab("B-cell score") + xlab(NULL) +
#   theme(legend.position = "right",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# plt_bcell.scores.clone <- bcell.score %>%
#   ggplot(aes(x=dat.sub$dna_best_class, 
#              y=dat.sub$bcell_score, 
#              fill=as.factor(dat.sub$dna_best_class),
#              colour=factor(dat.sub$timepoint, levels = c("diagnosis", "relapse 1", "relapse 2")))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "Clone") +
#   scale_color_discrete(name = "Timepoint") +
#   ylab("B-cell score") + xlab(NULL) +
#   theme(legend.position = "right",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
