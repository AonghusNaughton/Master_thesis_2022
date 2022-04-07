# New persister geneset

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

#################################################################################################################################

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

saveRDS(dat.sub, "dat.sub_filtered_10_cells.Rds")

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

# genenamesMarks <- lapply(marks, function(x){
#   rownames(x)
# })
# 
# # Subset markers to only include those that are common to all cases  
# common_genes <-genenamesMarks %>%
#   reduce(intersect)
# 
# filt_marks <- lapply(marks, function(x){
#   m <- x[rownames(x) %in% common_genes,]
#   return(m)
# })
# names(filt_marks) <- lapply(names(marks), function(x){
#   x
# })

# all combinations of patients and their indices
cmb <- combn(c(1:length(unique(dat.sub$patient_id))), 2)
cmb1 <- combn(unique(dat.sub$patient_id), 2)

# create dataframes that have stricter log-fold cut-offs
m_up_d15_high <- list()
for (i in 1:length(unique(dat.sub$patient_id))){
  m_up_d15_high[[i]] <- marks[[i]][marks[[i]]$avg_log2FC>0.5 & (marks[[i]]$pct.1>0.25 | marks[[i]]$pct.2>0.25),]
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

lapply(marks, function(x){
  table(persister_genes_high %in% rownames(x))
})
#################################################################################################################################

# Arrange dataframes and extract vectors of p-value. 

filt_marks <- lapply(marks, function(x){
  x <- x[persister_genes_high,]
  x <- na.omit(x)
  absent <- setdiff(persister_genes_high, rownames(x))
  new_rows <- lapply(absent, function(y){
    data.frame(y, "p_val"=1, "avg_log2FC"=0, "pct.1"=0, "pct.2"=0, "p_val_adj"=1)
  })
  new_rows <- do.call(rbind, new_rows)
  rownames(new_rows) <- new_rows[[1]]
  new_rows <- new_rows[, 2:length(colnames(new_rows))]
  x <- rbind(x, new_rows)
  x <- x[(rownames(x) %in% persister_genes_high),]
  return(x)
})

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

persister_genes <- lapply(filt_sig_marks, function(x){
  rownames(x)
})
persister_genes <- persister_genes[[1]]
saveRDS(persister_genes, "persister_genes.Rds")

pdf("/Users/aonghusnaughton/Proj_eng/April22/exp_heatmaps/persisters.pdf" , width = 12, height = 12)
DoHeatmap(dat.sub, features = persister_genes, group.by = "timepoint")
dev.off()
############################################################################################################

# dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
dat.sub <- readRDS("dat.sub_filtered_10_cells.Rds")
object <- dat.sub[,dat.sub$timepoint_short=="d15"]
nbin <- length(persister_genes)
ctrl = 10
pool <- NULL %||% rownames(x = object)
data.avg <- Matrix::rowMeans(x = GetAssayData(object = object)[pool, ])
data.avg <- data.avg[order(data.avg)]
data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = nbin)
features.use <- persister_genes
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

random_persister <- c()

for (i in 1:length(persister_genes)){
  random_persister <- c(random_persister, sample(ctrl.use[[i]], 1))
}

test <- FetchData(object = object, persister_genes)
sort(as.matrix(colMeans(test  > 0))*100)
test1 <- FetchData(object, random_persister)
sort(as.matrix(colMeans(test1  > 0))*100)

saveRDS(random_persister, "random_persister_220407.Rds")

