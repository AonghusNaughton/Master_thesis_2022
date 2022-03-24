library(Seurat)
library(tidyverse)
library(SC3)
library(ggplotify)
library(pheatmap)
library(ggplotify)
load("~/Proj_eng/ALL_data_Feb22.Rda")
sc3_split <- readRDS("~/Proj_eng/sc3_split.Rds")
sc3_split_wNa <- readRDS("~/Proj_eng/sc3_split_wNA.Rds")
sc3_ALL71_d5_split_wNa <- readRDS("sc3_ALL71_d5_split_wNa.Rds")

dat.sub <- dat[,dat$dna_cell_type=="blasts"  & 
                 dat$rna_cell_type=="blasts" & 
                 dat$patient_id=="ALL3" & 
                 dat$cell_phase == "G1"]

split <- split(sc3_split$ALL3.d0$cell_id, sc3_split$ALL3.d0$sc3_2_clusters)


annot <- dat.sub@meta.data %>%
  mutate(sc3_clusters_to_compare=case_when(cell_id %in% split$`1` ~ 1,
                                           cell_id %in% split$`2` ~ 2))
                                           # cell_id %in% split$`3` ~ 3,
                                           # cell_id %in% split$`4` ~ 4,
                                           # cell_id %in% split$`5` ~ 5))
                                           # cell_id %in% split$`4` ~ 6))
                                           # cell_id %in% split1$`3` ~ 3))
dat.sub <- AddMetaData(dat.sub, annot)
dat.sub <- SetIdent(dat.sub, value = "sc3_clusters_to_compare")
m <- FindAllMarkers(dat.sub[,dat.sub$timepoint_short=="d0" & dat.sub$sc3_clusters_to_compare %in% c(1:2)])#& dat.sub$dna_cell_type=="blasts"])
m_sig <- m[m$p_val_adj<0.05 & m$avg_log2FC>0.5,]


DoHeatmap(dat.sub[,dat.sub$timepoint_short=="d0" & dat.sub$sc3_clusters_to_compare %in% c(1:2)], features = m_sig$gene, group.by = "sc3_clusters_to_compare") + 
  theme(axis.text.y = element_text(size = 10))

#####################################################################################################################

# clones.per.cluster takes data from sc3_split_wNa (NA values present but not assigned to DNA class)

mat <- sapply(clones.per.cluster_wNa$MRD_ALL67.d0$sc3_5_clusters, function(x){
  as.matrix(table(x))
})
mat_by_clone <- apply(mat, 2, function(x){
  x/sum(x)
})
mat_by_state <- apply(mat, 1, function(x){
  x/sum(x)
})
rownames(mat_by_clone) <- c(1:length(mat_by_clone[,1]))
colnames(mat_by_state) <- c(1:length(mat_by_state[1,]))

lapply(list(mat_by_clone, mat_by_state), function(x){
  as.ggplot(pheatmap(x, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none"))
})

####################################################################################################################
clusters <- sapply(2:4, function(x){
  c(paste("sc3_",x,"_clusters",sep = ""))
})

clones.per.cluster_ALL71_split <- lapply(clusters, function(y){
  split(sc3_ALL71_split1.d5[[y]], sc3_ALL71_split1.d5$dna_best_class)
})

names(clones.per.cluster_ALL71_split) <- lapply(clusters, function(x){
    x
})

mat <- sapply(clones.per.cluster$MRD_ALL71.d15$sc3_3_clusters, function(x){
  as.matrix(table(x))
})

mat_1 <- sapply(clones.per.cluster_ALL71_split$sc3_3_clusters, function(x){
  as.matrix(table(x))
})

merged_mat <- rbind(mat_1, mat[2,])

mat_by_clone <- apply(mat, 2, function(x){
  x/sum(x)
})
mat_by_state <- apply(mat, 1, function(x){
  x/sum(x)
})
rownames(mat_by_clone) <- c(1:length(mat_by_clone[,1]))
colnames(mat_by_state) <- c(1:length(mat_by_state[1,]))

lapply(list(mat_by_clone, mat_by_state), function(x){
  as.ggplot(pheatmap(x, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none"))
})

#####################################################################################################################

as.ggplot(sc3_plot_consensus(sc3_ALL71_d0_split_wNa, k=3,
                   show_pdata = c("nFeature_RNA",
                                  "nCount_RNA",
                                  "bcell_score",
                                  "dna_best_class")))
#####################################################################################################################

split <- split(sc3_split_wNa$MRD_ALL47.d2$cell_id, sc3_split_wNa$MRD_ALL47.d2$sc3_2_clusters)



annot <- dat.sub@meta.data %>%
  mutate(sc3_clusters_to_compare=case_when(cell_id %in% split$`1` ~ 1,
                                           cell_id %in% split$`2` ~ 2))
                                           # cell_id %in% split$`3` ~ 3))
dat.sub <- AddMetaData(dat.sub, annot)

#####################################################################################################################

sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(dat.sub[,dat.sub$sc3_clusters_to_compare ==2 ]@assays$RNA@counts),
    logcounts = log2(as.matrix(dat.sub[,dat.sub$sc3_clusters_to_compare == 2]@assays$RNA@counts) + 1)
  ),
  colData = dat.sub[,dat.sub$sc3_clusters_to_compare == 2]@meta.data, 
)
rowData(sce)$feature_symbol <- rownames(sce)
sc3_ALL5_split2.d15 <- sc3(sce, ks= 2:4, biology = T)

saveRDS(sc3_ALL5_split.d0, "sc3_ALL5_d0_split_wNa.Rds")

#####################################################################################################################

sce <- lapply(c(1:2), function(x){
  r <- SingleCellExperiment(
    assays = list(
      counts =as.matrix(dat.sub[,dat.sub$sc3_clusters_to_compare==x]@assays$RNA@counts),
      logcounts = log2(as.matrix(dat.sub[,dat.sub$sc3_clusters_to_compare==x]@assays$RNA@counts) + 1)
    ),
    colData = dat.sub[,dat.sub$sc3_clusters_to_compare==x]@meta.data,
  )
})

names(sce) <- lapply(c(1:2), function(x){
  x
})

for (i in c(1:2)){
  rowData(sce[[i]])$feature_symbol <- rownames(sce[[i]])
}

sc3_ALL47_split.d2 <- lapply(sce, function(x){
  c <- sc3(x, ks = 2:4, biology = T,)
  return(c)
})

saveRDS(sc3_ALL47_split.d2, "sc3_ALL47_split.d2.Rds")

#####################################################################################################################

b <- data.frame(dat.sub$bcell_score, dat.sub$clone_time)
b <- b[b$dat.sub.clone_time %in% c("ALL3_2.d0", "ALL3_4.d0", "ALL3_9.rel2"),]
b1 <- data.frame(dat.sub$bcell_score, dat.sub$dna_best_class)
b1 <- b1[b1$dat.sub.dna_best_class == "ALL3_9",]

b2 <- rbind(b, b1) 

b.m <- melt(b2, id.var = "dat.sub.dna_best_class")

p2 <- ggplot(data = b.m, aes(x=dat.sub.dna_best_class, y = value, fill=as.factor(dat.sub.dna_best_class))) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "DNA clone x time \n 2/4 = Diagnostic \n 9 = rel/rel2 combined") +
  ylab("b-cell score") + xlab(NULL) +
  theme(legend.position = "right")

beeswarm(b$dat.sub.bcell_score ~ b$dat.sub.clone_time,
         col = c("#3FA0FF", "#136b39", "#e53935"))


#####################################################################################################################

# Fractions of cells not in G1 split by clone

df <- as.data.frame.matrix((table(dat.sub$cell_phase, dat.sub$clone_time)))

df$base <- rownames(df)
df <- melt(df)
barplot(as.matrix(df))
ggplot(df, aes(fill = base, y = value, x = variable)) + geom_bar(stat = "identity") + theme_classic()

#####################################################################################################################



        