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
library(ggplotify)
source("general.R")
load("rna_seurat_wAnnot211101.Rda")


# Differential expression between diagnostic clusters in MRD_ALL64 and assigning marker gene sets to clusters at day 15 

dat.sub <- dat[,dat$dna_cell_type=="blasts" & dat$rna_cell_type=="blasts" & dat$patient_id=="MRD_ALL64" & dat$cell_phase=="G1"]

m_phase_genes <- parse_grp("geneset_m.grp")
m_phase_genes <- m_phase_genes[3:length(m_phase_genes)]
s_phase_genes <- parse_grp("geneset_S.grp")
s_phase_genes <- s_phase_genes[3:length(s_phase_genes)]
apoptosis_genes <- parse_grp("apoptosis_geneset.grp")
apoptosis_genes <- apoptosis_genes[3:length(apoptosis_genes)]


phase_sets <- list(s_phase_genes, m_phase_genes, apoptosis_genes)
names(phase_sets) <- c("s", "m", "apoptosis")
dat.sub <- AddModuleScore(dat.sub, features = phase_sets)
colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(phase_sets)+1):ncol(dat.sub@meta.data)] <- paste(names(phase_sets), "score_new", sep = "_")


split.d0 <- split(sc3_split_wNa$MRD_ALL64.d0$cell_id, sc3_split_wNa$MRD_ALL64.d0$sc3_4_clusters)
split.d15 <- split(sc3_split_wNa$MRD_ALL64.d15$cell_id, sc3_split_wNa$MRD_ALL64.d15$sc3_2_clusters)


annot <- dat.sub@meta.data %>%
  mutate(sc3_clusters_to_compare=case_when(cell_id %in% split.d0$`1` ~ 1,
                                           cell_id %in% split.d0$`2` ~ 2,
                                           cell_id %in% split.d0$`3` ~ 3,
                                           cell_id %in% split.d0$`4` ~ 4,
                                           cell_id %in% split.d15$`1` ~ 5,
                                           cell_id %in% split.d15$`2`~ 6))
dat.sub <- AddMetaData(dat.sub, annot) 

boxplot(split(dat.sub[,dat.sub$timepoint=="diagnosis"]$bcell_score, dat.sub[,dat.sub$timepoint=="diagnosis"]$sc3_clusters_to_compare))
boxplot(split(dat.sub[,dat.sub$timepoint=="day 15"]$bcell_score, dat.sub[,dat.sub$timepoint=="day 15"]$sc3_clusters_to_compare))
boxplot(split(dat.sub$bcell_score, dat.sub$sc3_clusters_to_compare))
boxplot(split(dat.sub[,dat.sub$timepoint=="diagnosis"]$apoptosis_score_new, dat.sub[,dat.sub$timepoint=="diagnosis"]$sc3_clusters_to_compare))

dat.sub <- SetIdent(dat.sub, value = "sc3_clusters_to_compare")

m <- FindAllMarkers(dat.sub[,(dat.sub$sc3_clusters_to_compare==1 | dat.sub$sc3_clusters_to_compare==2 | dat.sub$sc3_clusters_to_compare==3)])

m_sig <- m[(m$p_val_adj<0.05 & m$avg_log2FC>0.5) | (m$p_val_adj<0.05 & m$avg_log2FC < -0.5),]
m_sig_up <- m[m$p_val_adj<0.05 & m$avg_log2FC>0.5 ,]

m_sig_up_split <- split(m_sig_up$gene, m_sig_up$cluster)

dat.sub <- AddModuleScore(dat.sub, features = m_sig_up_split)

colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)] <- paste("d0_", names(m_sig_up_split), sep = "")

FeaturePlot(dat.sub, features = c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)]), split.by = "timepoint")  

clusters <- c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)])

mat <- lapply(c(5:6), function(x){
  m <- matrix(c(dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_2",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_3",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_1",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$sc3_clusters_to_compare,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$s_score_new, 
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$m_score_new,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$apoptosis_score_new),
              ncol = 7)
})


names(mat) <- lapply(c(5:6), function(x){
  x
})

for (i in 1:2){
  colnames(mat[[i]]) <- c(c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)]), "sc3_cluster", "s_score", "m_score", "apoptosis_score")
}


rownames(mat$`5`) <- colnames(dat.sub[,dat.sub$sc3_clusters_to_compare==5])
rownames(mat$`6`) <- colnames(dat.sub[,dat.sub$sc3_clusters_to_compare==6])

df <- as.data.frame(do.call(rbind, mat))
sorted_cols <- sort(colnames(df))
df <- df[,sorted_cols]

df.m <- melt(df, id.var = "sc3_cluster")

p <- ggplot(data = df.m, aes(x=variable, y=value, fill=as.factor(sc3_cluster))) + 
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Day 15 cluster") +
  ylab("score") +
  theme(legend.position = "right") + 
  ylim(-1, 3)

# Self values

matr <- lapply(c(1:4), function(x){
  m <- matrix(c(dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_2",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_3",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$"d0_1",
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$sc3_clusters_to_compare,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$s_score_new, 
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$m_score_new,
                dat.sub[,dat.sub$sc3_clusters_to_compare==x]$apoptosis_score_new),
              ncol = 7)
}) 


names(matr) <- lapply(c(1:4), function(x){
  x
})



for (i in 1:4){
  colnames(matr[[i]]) <- c(c(colnames(dat.sub@meta.data)[(ncol(dat.sub@meta.data)-length(m_sig_up_split)+1):ncol(dat.sub@meta.data)]), "sc3_cluster", "s_score", "m_score", "apoptosis_score")
}


for (i in 1:4){
  rownames(matr[[i]]) <- colnames(dat.sub[,dat.sub$sc3_clusters_to_compare==i])
}


df_self <- as.data.frame(do.call(rbind, matr))
sorted_self_cols <- sort(colnames(df_self))
df_self <- df_self[,sorted_self_cols]
df_self.m <- melt(df_self, id.var = "sc3_cluster")

p1 <- ggplot(data = df_self.m, aes(x=variable, y=value, fill=as.factor(sc3_cluster))) + 
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Diagnosis cluster") +
  ylab("diagnostic cluster score") +
  theme(legend.position = "right")

b <- data.frame(dat.sub$bcell_score, dat.sub$sc3_clusters_to_compare)
b.m <- melt(b, id.var = "dat.sub.sc3_clusters_to_compare")

p2 <- ggplot(data = b.m, aes(x=dat.sub.sc3_clusters_to_compare, y = value, fill=as.factor(dat.sub.sc3_clusters_to_compare)))+
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "SC3 cluster") +
  ylab("b-cell score") + xlab(NULL)
  theme(legend.position = "right")

sc3.d0 <- as.ggplot(sc3_plot_consensus(sc3_split_wNa$MRD_ALL64.d0, k = 4,
                                    show_pdata = c("nCount_RNA",
                                    "nFeature_RNA",
                                    "bcell_score",
                                    "dna_best_class")))
sc3.d15 <- as.ggplot(sc3_plot_consensus(sc3_split_wNa$MRD_ALL64.d15, k = 2,
                                        show_pdata = c("nCount_RNA",
                                        "nFeature_RNA",
                                        "bcell_score",
                                        "dna_best_class")))

figure1 <- plot_grid(plot_grid(sc3.d0, NULL, sc3.d15, labels = c("d0", "", "d15"), nrow = 1, rel_widths = c(0.48, 0.01, 0.48)), NULL, 
          plot_grid(NULL, p2, NULL, labels = c("b_cell scores", "", ""), nrow = 1, rel_widths = c(0.3, 0.4, 0.3)), nrow = 3, rel_heights = c(0.7, 0.1, 0.3)) 

figure2 <- plot_grid(p1, p, labels = c("'self' scores", "d0 scores at d15"))

heat <- DoHeatmap(dat.sub[, dat.sub$timepoint=="diagnosis"], features = m_sig$gene, group.by = "sc3_clusters_to_compare", label = F) + 
  theme(text = element_text(size = 6))


# linear model


