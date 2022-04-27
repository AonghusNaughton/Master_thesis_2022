library(tidyverse)
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
persister_genes <- readRDS("persister_genes.Rds")

# geneSets <- unique(c(relapse_genes, random_relapse))
timepoint <- c("d0", "rel")
pid <- "ALL4"
dat.sub <- dat.sub[,dat.sub$cell_phase=="G1" & 
                     (dat.sub$timepoint_short==timepoint[1] |dat.sub$timepoint_short==timepoint[2]) &
                     dat.sub$patient_id==pid]

minCount <- 10

counts <- as.data.frame(dat.sub@assays$RNA@counts)

df_na <- na_if(counts, 0)

genes_above_threshold <- apply(df_na, 1, function(row) sum(!is.na(row)))
genes_above_threshold <- names(genes_above_threshold[genes_above_threshold > minCount])
dat.sub <- subset(dat.sub, features= genes_above_threshold)


data <- list(d0 <- list(counts <- as.data.frame(dat.sub[,dat.sub$timepoint_short=="d0"]@assays$RNA@counts),
                        scaled <- as.data.frame(dat.sub[,dat.sub$timepoint_short=="d0"]@assays$RNA@scale.data)),
             rel <- list(counts <- as.data.frame(dat.sub[,dat.sub$timepoint_short=="rel"]@assays$RNA@counts),
                         scaled <- as.data.frame(dat.sub[,dat.sub$timepoint_short=="rel"]@assays$RNA@scale.data)))
names(data) <- c("d0", "rel")

data <- lapply(data, function(x){
  df <- x[[1]][persister_genes[persister_genes %in% rownames(x[[1]])],]
  df <- df[rowSums(df) > 0,]
  genes <- rownames(df)
  df2 <- x[[2]][genes,]
  return(df2)
})

cor <- lapply(data, function(x){
  cor(as.matrix(t(x)), method="spearman")
})


cor$d0 <- cor$d0 %>%
  as.table() %>% as.data.frame() %>%
  mutate(timepoint="d0")

cor$rel <- cor$rel %>%
  as.table() %>% as.data.frame() %>%
  mutate(timepoint="rel")

cor_for_filter <- do.call("rbind", cor)
filtered.cor  <- cor_for_filter %>%
    filter(Var1 != Var2) %>%
    group_by(timepoint, .add=T) %>%
    group_by(Var1, .add = T) %>%
    summarise(avg=mean(Freq)) %>%
  mutate(x_axis=case_when(timepoint=="rel"~ avg),
         y_axis=case_when(timepoint=="d0"~avg)) %>%
  ungroup(Var1) %>%
  ungroup(timepoint) %>%
  dplyr::select(-c(timepoint, avg)) %>%
  replace_na(list(x_axis=0, y_axis=0)) %>%
  group_by(Var1) %>%
  summarise(x_axis=sum(x_axis),
            y_axis=sum(y_axis))
  # ggplot(aes(x=x_axis, y=y_axis)) +
  # geom_point() + 
  # xlab(paste0("Relapse")) +
  # ylab(paste0("Diagnosis")) +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  # theme_classic()

filtered.cor %>%
  ggplot(aes(x=x_axis, y=y_axis)) +
  geom_point() + 
  xlab(paste0("Relapse")) +
  ylab(paste0("Diagnosis")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, col='red', size=1)
