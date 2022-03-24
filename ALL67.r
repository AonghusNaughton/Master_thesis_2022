library(tidyverse)
library(cowplot)
library(scales)
library(Seurat)
library(ggrepel)
library(ggsci)
library(EnhancedVolcano)
library(gridExtra)
library(ggplot2)
library(SC3)
library(reshape2)
library(cmapR)
source("general.R")
load("rna_seurat_wAnnot211101.Rda")

sc3.d0 <- sc3_plot_consensus(sc3_split_wNa$MRD_ALL67.d0, k = 5,
                             show_pdata = c("nFeature_RNA",
                                            "nCount_RNA",
                                            "bcell_score",
                                            "dna_best_class"))
sc3.15 <- sc3_plot_consensus(sc3_split_wNa$MRD_ALL67.d15, k = 2,
                             show_pdata = c("nFeature_RNA",
                                            "nCount_RNA",
                                            "bcell_score",
                                            "dna_best_class"))


m_phase_genes <- parse_grp("geneset_m.grp")
m_phase_genes <- m_phase_genes[3:length(m_phase_genes)]
s_phase_genes <- parse_grp("geneset_S.grp")
s_phase_genes <- s_phase_genes[3:length(s_phase_genes)]
apoptosis_genes <- parse_grp("apoptosis_geneset.grp")
apoptosis_genes <- apoptosis_genes[3:length(apoptosis_genes)]

dat.sub <- dat[,dat$dna_cell_type=="blasts" & dat$rna_cell_type=="blasts" & dat$patient_id=="MRD_ALL67" & dat$cell_phase=="G1"]

phase_sets <- list(s_phase_genes, m_phase_genes, apoptosis_genes)
names(phase_sets) <- c("s", "m", "apoptosis")
dat.sub <- AddModuleScore(dat.sub, features = phase_sets)
colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(phase_sets)+1):ncol(dat.sub@meta.data)] <- paste(names(phase_sets), "score_new", sep = "_")


split.d0 <- split(sc3_split_wNa$MRD_ALL67.d0$cell_id, sc3_split_wNa$MRD_ALL67.d0$sc3_5_clusters)
split.d15 <- split(sc3_split_wNa$MRD_ALL67.d15$cell_id, sc3_split_wNa$MRD_ALL67.d15$sc3_2_clusters)

annot <- dat.sub@meta.data %>%
  mutate(sc3_clusters_to_compare=case_when(cell_id %in% split.d0$`1` ~ 1,
                                           cell_id %in% split.d0$`2` ~ 2,
                                           cell_id %in% split.d0$`3` ~ 3,
                                           cell_id %in% split.d0$`4` ~ 4,
                                           cell_id %in% split.d0$`5` ~ 5,
                                           cell_id %in% split.d15$`1` ~ 6,
                                           cell_id %in% split.d15$`2`~ 7))
dat.sub <- AddMetaData(dat.sub, annot)

boxplot(split(dat.sub[,dat.sub$timepoint=="diagnosis"]$bcell_score, dat.sub[,dat.sub$timepoint=="diagnosis"]$sc3_clusters_to_compare))
boxplot(split(dat.sub[,dat.sub$timepoint=="diagnosis"]$apoptosis_score_new, dat.sub[,dat.sub$timepoint=="diagnosis"]$sc3_clusters_to_compare))

dat.sub <- SetIdent(dat.sub, value = "sc3_clusters_to_compare")

m <- FindAllMarkers(dat.sub[,dat.sub$sc3_clusters_to_compare%in%c(1,2,3,5)])

m_sig_up <- m[m$p_val_adj<0.05 & m$avg_log2FC>0.5 ,]
m_sig_up_split <- split(m_sig_up$gene, m_sig_up$cluster)


# Cluster 3 genes dying -- RPL genes/TP53 

dat.sub <- AddModuleScore(dat.sub, features = m_sig_up_split)

colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)] <- paste("d0_", names(m_sig_up_split), sep = "")

FeaturePlot(dat.sub, features = c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)]), split.by = "timepoint")  

clusters <- c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)])

mat <- lapply(c(6:7), function(x){
  m <- matrix(c(dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_5",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_3",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_2",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_1",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$sc3_clusters_to_compare,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$s_score_new, 
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$m_score_new,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$apoptosis_score_new),
              ncol = 8)
})

names(mat) <- lapply(c(6:7), function(x){
  x
})

for (i in 1:2){
  colnames(mat[[i]]) <- c(c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)]), "sc3_cluster", "s_score", "m_score", "apoptosis_score")
}


rownames(mat$`6`) <- colnames(dat.sub[,dat.sub$sc3_clusters_to_compare==6])
rownames(mat$`7`) <- colnames(dat.sub[,dat.sub$sc3_clusters_to_compare==7])

df <- as.data.frame(do.call(rbind, mat))

df.m <- melt(df, id.var = "sc3_cluster")

p <- ggplot(data = df.m, aes(x=variable, y=value, fill=as.factor(sc3_cluster))) + 
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Day 15 cluster") +
  ylab("score") +
  theme(legend.position = "right")

b <- data.frame(dat.sub$bcell_score, dat.sub$sc3_clusters_to_compare)
b.m <- melt(b, id.var = "dat.sub.sc3_clusters_to_compare")

p1 <- ggplot(data = b.m, aes(x=dat.sub.sc3_clusters_to_compare, y = value, fill=as.factor(dat.sub.sc3_clusters_to_compare)))+
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "SC3 cluster") +
  ylab("b-cell score") + xlab(NULL)
theme(legend.position = "right")


# Differential gene expression analysis between day 15 cluster

m_15 <- FindAllMarkers(dat.sub[,dat.sub$timepoint_short=="d15"])

cl6_genes <- m_15[m_15$cluster==6,]
rownames(cl6_genes) <- cl6_genes$gene
cl6_sig <- cl6_genes[cl6_genes$p_val<0.05 & (cl6_genes$avg_log2FC > 0.5 | cl6_genes$avg_log2FC < -0.5),]


EnhancedVolcano(cl6_genes,
                lab = rownames(cl6_genes),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Differential expression between cluster d15 (specific to cluster 6)",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                legendLabels=c('Not sig.','Log (base 2) FC','adj p-value',
                               'adj p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 5,
                legendIconSize = 2.0)
