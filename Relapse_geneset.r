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
# Create pseudo gene set that represents relapse gene set (similar pct expression)

# exp <- FetchData(dat.sub[,(dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")], true_persister_gene_set_rel)
# average_pct.exp_rel <- colMeans(as.matrix(colMeans(exp  > 0))*100)
# st_dev.exp_rel <- sd(as.matrix(colMeans(exp >0)*100)[,1])
# range <- (average_pct.exp_rel - st_dev.exp_rel) : (average_pct.exp_rel + st_dev.exp_rel)
# 
# exp.vals_for_pseudo <- FetchData(dat.sub[,(dat.sub$timepoint_short=="rel" | 
#                                              dat.sub$timepoint_short=="rel2")], 
#                                  rownames(dat.sub[,(dat.sub$timepoint_short=="rel" | 
#                                                       dat.sub$timepoint_short=="rel2")])[rownames(dat.sub[,(dat.sub$timepoint_short=="rel" | 
#                                                                                                               dat.sub$timepoint_short=="rel2")]) %!in% true_persister_gene_set_rel])
# pct.exp <- as.matrix(colMeans(exp.vals_for_pseudo  > 0))*100
# pct.exp_filt <- pct.exp[pct.exp>range[1] & pct.exp<range[length(range)],]
# random_indices <- sample(length(pct.exp_filt) + 1, length(true_persister_gene_set_rel) + 1, replace = F)
# pseudo_relapse<- pct.exp_filt[random_indices]
# mean(pseudo_relapse)
# 
# pseudo_relapse <- names(pseudo_relapse)


#################################################################################################################################

# Pseudo gene set that represents a relapse-like gene set
# More accuarate way -- create n bins of all genes expressed where n=length of persister genes. Bins are based on average expression
# across cells across patients that persister gene set was derived from. 

object <- dat.sub[,(dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")]
nbin <- length(true_persister_gene_set_rel)
ctrl <- 10
pool <- NULL %||% rownames(x = object)
data.avg <- Matrix::rowMeans(x = GetAssayData(object = object)[pool, ])
data.avg <- data.avg[order(data.avg)]
data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
names(x = data.cut) <- names(x = data.avg)
features.use <- c()
cluster.length <- nbin
features.use <- true_persister_gene_set_rel
ctrl.use <- vector(mode = "list", length = cluster.length)
for (i in 1:nbin) {
  for (j in 1:length(x = features.use)) {
    ctrl.use[[i]] <- c(
      ctrl.use[[i]],
      names(x = sample(
        x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
        size = ctrl,
        replace = FALSE
      ))
    )
  }
}

pseudo_relapse <- c()

for (i in 1:length(true_persister_gene_set_rel)){
  pseudo_relapse <- c(pseudo_relapse, sample(ctrl.use[[i]], 1))
}

test <- FetchData(object = object, true_persister_gene_set_rel)
sort(as.matrix(colMeans(test  > 0))*100)

test1 <- FetchData(object, pseudo_relapse)
sort(as.matrix(colMeans(test1  > 0))*100)

saveRDS(pseudo_relapse, "pseudo_relapse.Rds")
