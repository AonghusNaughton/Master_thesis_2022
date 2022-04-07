# New relapse geneset 

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

# Subet to only include cases that have diagnostic and relapse samples
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


# To include only highly expressed genes: 
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

dat.sub <- subset(dat.sub, features= genes_in_at_least_two)

###################################################################################################
# Find differentially expressed genes 

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

# Subset markers to only include those that are at least expressed in all cases  
# common_genes <-genenamesMarks %>%
#   reduce(intersect)
# 
# filt_marks_rel <- lapply(marks_rel, function(x){
#   m <- x[rownames(x) %in% common_genes,]
#   return(m)
# })
# names(filt_marks_rel) <- lapply(names(marks_rel), function(x){
#   x
# })

# all combinations of patients and their indices
cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
cmb1 <- combn(unique(dat.sub$patient_id), 2)

# create dataframes that have stricter log-fold cut-offs
m_up_rel_high <- list()
for (i in 1:length(unique(dat.sub$patient_id))){
  m_up_rel_high[[i]] <- marks_rel[[i]][marks_rel[[i]]$avg_log2FC>0.5 &
                                         (marks_rel[[i]]$pct.1>0.25 | marks_rel[[i]]$pct.2>0.25),]
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
# 
# relapse_genes_high <- c()
# for (i in 1:length(names(combinations_high_rel))){
#   for (j in 1:length(combinations_high_rel[[i]])){
#     for (k in 1:length(names(combinations_high_rel))){
#       if (k != i) {
#         if (combinations_high_rel[[i]][j] %in% combinations_high_rel[[k]]){
#           relapse_genes_high <- c(relapse_genes_high, combinations_high_rel[[i]][j])
#         }
#       }
#     }
#   }
# }
# 
# relapse_genes_high <- unique(relapse_genes_high)


# Must have a high log fold change in at least two cases. 
# Append to list if gene is present in any pair of patients. Option to increase this threshold to 3 patients above (commented code)

relapse_genes_high <- c()

for (i in 1:length(names(combinations_high_rel))){
  relapse_genes_high <- c(relapse_genes_high, combinations_high_rel[[i]])
}

relapse_genes_high <- unique(relapse_genes_high)

lapply(marks_rel, function(x){
  table(relapse_genes_high %in% rownames(x))
})

lapply(marks_rel, function(x){
  which(!relapse_genes_high %in% rownames(x))
})
#################################################################################################################################

# Arrange dataframes and extract vectors of p-value. 

filt_marks_rel <- lapply(marks_rel, function(x){
  x <- x[relapse_genes_high,]
  x <- na.omit(x)
  absent <- setdiff(relapse_genes_high, rownames(x))
  new_rows <- lapply(absent, function(y){
    data.frame(y, "p_val"=1, "avg_log2FC"=0, "pct.1"=0, "pct.2"=0, "p_val_adj"=1)
  })
  new_rows <- do.call(rbind, new_rows)
  rownames(new_rows) <- new_rows[[1]]
  new_rows <- new_rows[, 2:length(colnames(new_rows))]
  x <- rbind(x, new_rows)
  x <- x[(rownames(x) %in% relapse_genes_high),]
  return(x)
})

lapply(filt_marks_rel, function(x) length(rownames(x)))

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

lapply(filt_sig_marks_rel, function(x) length(rownames(x)))

relapse_genes <- lapply(filt_sig_marks_rel, function(x){
  rownames(x)
})

relapse_genes <- relapse_genes[[1]]
saveRDS(relapse_genes, "relapse_genes.Rds")

pdf("/Users/aonghusnaughton/Proj_eng/April22/exp_heatmaps/relapse.pdf" , width = 12, height = 12)
DoHeatmap(dat.sub, features = relapse_genes, group.by = "timepoint")
dev.off()
