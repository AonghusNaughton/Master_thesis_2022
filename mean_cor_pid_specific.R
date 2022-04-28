dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
dat.sub <- dat.sub[,(dat.sub$timepoint_short=="d0" | 
                       dat.sub$timepoint_short=="rel") & 
                     dat.sub$patient_id=="ALL4"]
counts <- as.data.frame(dat.sub@assays$RNA@counts)
counts_na <- na_if(counts, 0)
minCount <- 10
genes <- apply(counts_na, 1, function(row) sum(!is.na(row)))
genes <- names(genes[genes > minCount])
counts <- counts[genes,]

clones <- unique(dat.sub$clone_time[!dat.sub$clone_time %in% c("NA", "ALL4_6.rel")])

data <- lapply(clones, function(x){
  counts <-  as.data.frame(dat.sub[,dat.sub$clone_time==x]@assays$RNA@counts)
  scaled <- as.data.frame(dat.sub[,dat.sub$clone_time==x]@assays$RNA@scale.data)
  
  counts <- counts[persister_genes[persister_genes %in% rownames(counts)],]
  counts <- counts[rowSums(counts) > 0,]
  genes <- rownames(counts)
  scaled <- scaled[genes,]
  return(scaled)
})
names(data) <- lapply(clones, function(x) x)

cor <- lapply(data, function(x){
  cor(as.matrix(t(x)), method="spearman")
})

timepoints <- c()
cln <- c()
for (i in clones){
  timepoints <- c(timepoints, unique(as.character(dat.sub[,dat.sub$clone_time==i]$timepoint_short)))
  cln <- c(cln, unique(dat.sub[,dat.sub$clone_time==i]$dna_best_class))
}
timepoints <- as.list(timepoints)
names(timepoints) <- lapply(clones, function(x) x)
cln <- as.list(cln)
names(cln) <- lapply(clones, function(x) x)

cor_mod <- lapply(clones, function(x){
  cor[[x]] %>%
    as.table() %>% as.data.frame() %>%
    mutate(timepoint=timepoints[[x]],
           clone=cln[[x]])
})
names(cor_mod) <- lapply(clones, function(x) x)
combined_cor <- do.call("rbind", cor_mod)

filtered.cor <- combined_cor %>%
  filter(Var1 != Var2) %>%
  group_by(clone, .add = T) %>%
  group_by(timepoint, .add = T) %>%
  group_by(Var1, .add = T) %>%
  summarise(avg=mean(Freq)) %>%
  ungroup(timepoint, clone) %>%
  group_by(clone, .add = T)

p1 <- filtered.cor %>%
  ggplot(aes(x=clone, y=avg)) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  # scale_fill_discrete(name = "Clone") +
  # scale_color_discrete(name = "Cell type") +
  ylab("Mean internal correlation") + xlab("Clone") +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  facet_grid(.~timepoint, scales = "free_x")

p1 <- filtered.cor %>%
  ggplot(aes(x=x_axis, y=y_axis)) +
  geom_point() + 
  xlab(paste0("Relapse")) +
  ylab(paste0("Diagnosis")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  # geom_label(data= filtered.cor %>%
  #   filter((x_axis==max(x_axis) & y_axis==max(y_axis))),
  # aes(label=Var1), hjust=-0.5) +
  facet_grid(.~ clone)

a <- filtered.cor %>%
  filter(clone=="ALL4_3")
table(a$y_axis > 0)

clones <- unique(dat.sub$dna_best_class)

lapply(clones, function(x){
  a <- filtered.cor %>%
    filter(clone==x)
  return(table(a$y_axis>0))
})

lapply(clones, function(x){
  a <- filtered.cor %>%
    filter(clone==x) %>%
    group_by(clone) %>%
    summarise(mean_of_means=mean(y_axis))
})

#######################################################################################
# ALL3

dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
dat.sub <- dat.sub[,(dat.sub$timepoint_short=="d0" | 
                       dat.sub$timepoint_short=="rel" | 
                       dat.sub$timepoint_short=="rel2") & 
                     dat.sub$patient_id=="ALL3" &
                     !is.na(dat.sub$dna_cell_type)]
counts <- as.data.frame(dat.sub@assays$RNA@counts)
counts_na <- na_if(counts, 0)
minCount <- 10
genes <- apply(counts_na, 1, function(row) sum(!is.na(row)))
genes <- names(genes[genes > minCount])
counts <- counts[genes,]

clones <- unique(dat.sub$clone_time[!dat.sub$clone_time %in% c("ALL3_2.rel2", "ALL3_7.rel2", 
                                                               "ALL3_8.rel2", "ALL3_6.rel2")])

data <- lapply(clones, function(x){
  counts <-  as.data.frame(dat.sub[,dat.sub$clone_time==x]@assays$RNA@counts)
  scaled <- as.data.frame(dat.sub[,dat.sub$clone_time==x]@assays$RNA@scale.data)
  
  counts <- counts[persister_genes[persister_genes %in% rownames(counts)],]
  counts <- counts[rowSums(counts) > 0,]
  genes <- rownames(counts)
  scaled <- scaled[genes,]
  return(scaled)
})
names(data) <- lapply(clones, function(x) x)

cor <- lapply(data, function(x){
  cor(as.matrix(t(x)), method="spearman")
})

timepoints <- c()
cln <- c()
for (i in clones){
  timepoints <- c(timepoints, unique(as.character(dat.sub[,dat.sub$clone_time==i]$timepoint_short)))
  cln <- c(cln, unique(dat.sub[,dat.sub$clone_time==i]$dna_best_class))
}
timepoints <- as.list(timepoints)
names(timepoints) <- lapply(clones, function(x) x)
cln <- as.list(cln)
names(cln) <- lapply(clones, function(x) x)

cor_mod <- lapply(clones, function(x){
  cor[[x]] %>%
    as.table() %>% as.data.frame() %>%
    mutate(timepoint=timepoints[[x]],
           clone=cln[[x]])
})
names(cor_mod) <- lapply(clones, function(x) x)

combined_cor <- do.call("rbind", cor_mod)

filtered.cor <- combined_cor %>%
  filter(Var1 != Var2) %>%
  group_by(clone, .add = T) %>%
  group_by(timepoint, .add = T) %>%
  group_by(Var1, .add = T) %>%
  summarise(avg=mean(Freq)) %>%
  ungroup(timepoint, clone) %>%
  group_by(clone, .add = T)


p1 <- filtered.cor %>%
  ggplot(aes(x=clone, y=avg)) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  # scale_fill_discrete(name = "Clone") +
  # scale_color_discrete(name = "Cell type") +
  ylab("Mean internal correlation") + xlab("Clone") +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  facet_grid(.~timepoint, scales = "free_x")


