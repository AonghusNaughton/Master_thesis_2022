dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")

names <- unique(dat.sub[,(dat.sub$timepoint_short=="d15" | dat.sub$timepoint_short=="rel")]$patient_id)

dat.sub <- dat.sub[,(dat.sub$timepoint_short=="d0" | 
                       dat.sub$timepoint_short=="rel" | 
                       dat.sub$timepoint_short=="rel2" | 
                       dat.sub$timepoint_short=="d15") & 
                     dat.sub$patient_id %in% names[names != "MRD_ALL47"]]

to_remove <- dat.sub[,dat.sub$patient_id=="ALL3" & is.na(dat.sub$dna_cell_type)]$cell_id
dat.sub <- dat.sub[, !colnames(dat.sub) %in% to_remove]

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

names <- unique(dat.sub$pid_time)
data <- lapply(names, function(x){
  counts <-  as.data.frame(dat.sub[,dat.sub$pid_time==x]@assays$RNA@counts)
  scaled <- as.data.frame(dat.sub[,dat.sub$pid_time==x]@assays$RNA@scale.data)
  
  counts <- counts[persister_genes[persister_genes %in% rownames(counts)],]
  counts <- counts[rowSums(counts) > 0,]
  genes <- rownames(counts)
  scaled <- scaled[genes,]
  return(scaled)
})

names(data) <- lapply(names, function(x) x)

cor <- lapply(data, function(x){
  cor(as.matrix(t(x)), method="spearman")
})

timepoints <- c()
pids <- c()
for (i in names){
  timepoints <- c(timepoints, unique(as.character(dat.sub[,dat.sub$pid_time==i]$timepoint_short)))
  pids <- c(pids, unique(dat.sub[,dat.sub$pid_time==i]$patient_id))
}
timepoints <- as.list(timepoints)
pids <- as.list(pids)
names(timepoints) <- lapply(names, function(x) x)
names(pids) <- lapply(names, function(x) x)

MRD_pids <- unique(dat.sub[,dat.sub$timepoint_short=="d15"]$patient_id)
rel_pids <- unique(dat.sub[,dat.sub$timepoint_short=="rel"]$patient_id)

cor_mod <- lapply(names, function(x){
  cor[[x]] %>%
    as.table() %>% as.data.frame() %>%
    mutate(timepoint=timepoints[[x]],
           patient=pids[[x]])
})

names(cor_mod) <- lapply(names, function(x) x)

combined_cor <- do.call("rbind", cor_mod)

filtered.cor <- combined_cor %>%
  filter(Var1 != Var2) %>%
  group_by(patient, .add = T) %>%
  group_by(timepoint, .add = T) %>%
  group_by(Var1, .add = T) %>%
  summarise(avg=mean(Freq)) %>%
  ungroup(timepoint) %>%
  group_by(patient, .add = T) %>%
  mutate(status=case_when(patient %in% MRD_pids ~ "MRD",
                          patient %in% rel_pids ~ "Relapse"))

p1 <- filtered.cor %>%
  ggplot(aes(x=patient, y=avg, fill=status, color=patient)) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(NULL) +
  scale_color_discrete(name = "patient") +
  ylab("Mean internal correlation") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  facet_grid(.~timepoint, scales = "free_x")


