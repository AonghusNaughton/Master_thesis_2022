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

# exp <- FetchData(dat.sub[,dat.sub$timepoint_short=="d15"], true_persister_gene_set_d15)
# average_pct.exp_persister <- colMeans(as.matrix(colMeans(exp  > 0))*100)
# st_dev.exp_persister <- sd(as.matrix(colMeans(exp >0)*100)[,1])
# range <- (average_pct.exp_persister - st_dev.exp_persister) : (average_pct.exp_persister + st_dev.exp_persister)
# 
# exp.vals_for_pseudo <- FetchData(dat.sub[,dat.sub$timepoint_short=="d15"], 
#                                  rownames(dat.sub[,dat.sub$timepoint_short=="d15"])[rownames(dat.sub[,dat.sub$timepoint_short=="d15"]) %!in% 
#                                                                                       true_persister_gene_set_d15])
# pct.exp <- as.matrix(colMeans(exp.vals_for_pseudo  > 0))*100
# pct.exp_filt <- pct.exp[pct.exp>range[1] & pct.exp<range[length(range)],]
# random_indices <- sample(length(pct.exp_filt) + 1, length(true_persister_gene_set_d15) + 1, replace = F)
# pseudo_persister <- pct.exp_filt[random_indices]
# 
# pseudo_persister <- names(pseudo_persister)


#################################################################################################################################

# Pseudo gene set that represents a persister-like gene set. 
# More accuarate way -- create n bins of all genes expressed where n=length of persister genes. Bins are based on average expression
# across cells across patients that persister gene set was derived from. 

object <- dat.sub[,dat.sub$timepoint_short=="d15"]
nbin <- length(true_persister_gene_set_d15)
ctrl = 10
pool <- NULL %||% rownames(x = object)
data.avg <- Matrix::rowMeans(x = GetAssayData(object = object)[pool, ])
data.avg <- data.avg[order(data.avg)]
data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = nbin)
features.use <- true_persister_gene_set_rel
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

pseudo_persister <- c()

for (i in 1:length(true_persister_gene_set_d15)){
  pseudo_persister <- c(pseudo_persister, sample(ctrl.use[[i]], 1))
}

test <- FetchData(object = object, true_persister_gene_set_d15)
sort(as.matrix(colMeans(test  > 0))*100)
test1 <- FetchData(object, pseudo_persister)
sort(as.matrix(colMeans(test1  > 0))*100)

saveRDS(pseudo_persister, "pseudo_persister.Rds")


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

