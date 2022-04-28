library(UCell)
library(Seurat)

# Persister module scores at relapse 

persister_genes <- readRDS("persister_genes.Rds")
persister_genes <- list(persister_genes)

dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
names <- unique(dat.sub$patient_id)

dat.pid_UCell <- lapply(names, function(x){
  if(x=="ALL3"){
    sub <- dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts"]
  } else {
    sub <- dat.sub[,dat.sub$patient_id==x]
  }
  sub <- AddModuleScore_UCell(sub, features = persister_genes)
  colnames(sub@meta.data)[(ncol(sub@meta.data)-length(persister_genes)+1):ncol(sub@meta.data)] <- "persister_score"
  return(sub)
})

dat.pid_seurat <- lapply(names, function(x){
  if(x=="ALL3"){
    sub <- dat.sub[,dat.sub$patient_id==x & dat.sub$dna_cell_type=="blasts"]
  } else {
    sub <- dat.sub[,dat.sub$patient_id==x]
  }
  sub <- AddModuleScore(sub, features = persister_genes)
  colnames(sub@meta.data)[(ncol(sub@meta.data)-length(persister_genes)+1):ncol(sub@meta.data)] <- "persister_score"
  return(sub)
})

names(dat.pid_UCell) <- lapply(names,function(x) x)
names(dat.pid_seurat) <- lapply(names, function(x) x)

relapse_cases <- unique(dat.sub[,(dat.sub$timepoint_short=="rel" | dat.sub$timepoint_short=="rel2")]$patient_id)
MRD_cases <- unique(dat.sub$patient_id[!dat.sub$patient_id %in% relapse_cases])

score_data_UCell <- lapply(dat.pid_UCell, function(x){
  tibble(x$patient_id, x$persister_score, x$dna_best_class, x$timepoint, x$bcell_score) %>%
    mutate(case=case_when(`x$patient_id` %in% relapse_cases ~ "Relapse",
                          `x$patient_id` %in% MRD_cases ~ "MRD"))
})

score_data_seurat <- lapply(dat.pid_seurat, function(x){
  tibble(x$patient_id, x$persister_score, x$dna_best_class, x$timepoint, x$bcell_score) %>%
    mutate(case=case_when(`x$patient_id` %in% relapse_cases ~ "Relapse",
                          `x$patient_id` %in% MRD_cases ~ "MRD"))
})

# persister_plot <- ggplot(bind_rows(score_data, .id = "df"), 
#                           aes(x=`x$patient_id`,
#                               y=`x$persister_score`, 
#                               fill=as.factor(`x$timepoint`),
#                               colour= as.factor(`x$dna_best_class`))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "Clone") +
#   scale_color_discrete(name = "Timepoint") +
#   ylab("Persister score") + xlab(NULL) +
#   theme(legend.position = "right",
#         legend.title = element_text(size = 5), 
#         legend.text = element_text(size = 5)) + 
#   guides(shape = guide_legend(override.aes = list(size = 0.25)),
#          color = guide_legend(override.aes = list(size = 0.25))) +
#   facet_grid(.~`case`, scales = "free_x")
# ggsave("/Users/aonghusnaughton/Proj_eng/April22/Persister_scores/persister_score.pdf", plot = persister_plot, width = 15, height = 8)


# timepoint_plot <- ggplot(bind_rows(score_data, .id = "df"), 
#                          aes(x=`x$patient_id`,
#                              y=`x$persister_score`, 
#                              fill=as.factor(`x$timepoint`))) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "Timepoint") +
#   # scale_color_discrete(name = "Timepoint") +
#   ylab("Persister score") + xlab(NULL) +
#   theme(legend.position = "right",
#         legend.title = element_text(size = 5), 
#         legend.text = element_text(size = 5)) + 
#   guides(shape = guide_legend(override.aes = list(size = 0.25)),
#          color = guide_legend(override.aes = list(size = 0.25)))  +
#   ylim(-0.5, 1.2) +
#   facet_grid(.~`case`, scales = "free_x")

bcell_plot_seurat <- ggplot(bind_rows(score_data_seurat, .id = "df"), 
                         aes(x=`x$patient_id`,
                             y=`x$bcell_score`, 
                             fill=as.factor(`x$timepoint`))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Timepoint") +
  # scale_color_discrete(name = "Timepoint") +
  ylab("B-cell score_seurat") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.25)),
         color = guide_legend(override.aes = list(size = 0.25))) +
  ylim(-0.5, 1.2) +
  facet_grid(.~`case`, scales = "free_x")

bcell_plot_UCell <- ggplot(bind_rows(score_data_seurat, .id = "df"),
                            aes(x=`x$patient_id`,
                                y=`x$bcell_score`,
                                fill=as.factor(`x$timepoint`))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) +
  theme_classic() +
  scale_fill_discrete(name = "Timepoint") +
  # scale_color_discrete(name = "Timepoint") +
  ylab("B-cell score_UCell") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  guides(shape = guide_legend(override.aes = list(size = 0.25)),
         color = guide_legend(override.aes = list(size = 0.25))) +
  ylim(-0.5, 1.2) +
  facet_grid(.~`case`, scales = "free_x")


seurat_comp <- bcell_plot_UCell + bcell_plot_seurat

ggsave("/Users/aonghusnaughton/Proj_eng/April22/Persister_scores/bcell_score_seurat.pdf",plot = seurat_comp, width = 25, height = 10)





