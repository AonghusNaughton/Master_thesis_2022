library(ggpubr)

# b cell score comparisons
dat.sub <- dat[,dat$patient_id=="ALL3" & 
                 dat$cell_phase=="G1" & !is.na(dat$dna_cell_type)]

bcell.score <- na.omit(tibble(dat.sub$bcell_score, dat.sub$timepoint, dat.sub$dna_best_class, dat.sub$rna_cell_subtype))

plt_bcell.scores.clone <- bcell.score %>%
  ggplot(aes(x=as.factor(`dat.sub$dna_best_class`), 
             y=`dat.sub$bcell_score`, 
             fill=as.factor(`dat.sub$dna_best_class`),
             colour= as.factor(`dat.sub$rna_cell_subtype`))) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Clone") +
  scale_color_discrete(name = "Cell type") +
  ylab("B-cell score") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  ylim(-0.5, 0.5)

ggsave("bcell_ALL3_clones.pdf")

dat.sub1 <- dat[,(dat$dna_cell_type=="blasts" | is.na(dat$dna_cell_type)) & dat$rna_cell_type=="blasts" & dat$cell_phase=="G1"]
dat.sub1 <- SetIdent(dat.sub1, value = "pid_time")
names <- unique(dat.sub1$pid_time)
names.d15 <- names[grepl(".d15", names)]
names.d0 <- gsub(".d15", ".d0", names.d15)
names.all <- c(names.d0, names.d15)
names.all <- names.all[!names.all %in% c("MRD_ALL47.d0", "MRD_ALL47.d15")] #only one cell at d15
dat.sub1 <- subset(dat.sub1, idents = names.all[1:length(names.all)])
table(dat.sub1$pid_time)
dim(dat.sub1)

bcell.score1 <- tibble(dat.sub1$bcell_score, dat.sub1$timepoint, dat.sub1$patient_id)

plt <- bcell.score1 %>%
  ggplot(aes(x=factor(`dat.sub1$timepoint`, levels = c("diagnosis", "day 15")), 
             y=`dat.sub1$bcell_score`, 
             fill=factor(`dat.sub1$timepoint`, levels = c("diagnosis", "day 15")),
             colour=`dat.sub1$patient_id`)) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Timepoint") +
  scale_color_discrete(name = "Patient ID") +
  ylab("B-cell score") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  ylim(-0.5, 0.5)

ggsave("../bcell_d0vsd15.pdf")

plt2 <- bcell.score1 %>%
  ggplot(aes(x=`dat.sub1$patient_id`,
             y=`dat.sub1$bcell_score`, 
             fill=factor(`dat.sub1$timepoint`, levels = c("diagnosis", "day 15")),
             colour=`dat.sub1$patient_id`)) +
  stat_boxplot(geom = 'errorbar', width= 0.6)+
  geom_boxplot(width = 0.6) +
  expand_limits(y = 0) + 
  theme_classic() +
  scale_fill_discrete(name = "Timepoint") +
  scale_color_discrete(name = "Patient ID") +
  ylab("B-cell score") + xlab(NULL) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5)) + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  ylim(-0.5, 0.5)
ggsave("../bcell_d0vsd15_2.pdf")

ggarrange(plt_bcell.scores.clone, plt2, nrow = 1)
ggsave("../bcell_ALL3vsMRD.pdf", width = 12, height = 5)

# bcell.score <- tibble(dat.sub$bcell_score, dat.sub$timepoint)
# 
# plt3 <- bcell.score %>%
#   ggplot(aes(x=`dat.sub$timepoint`,
#              y=`dat.sub$bcell_score`,
#              fill=`dat.sub$timepoint`)) +
#   stat_boxplot(geom = 'errorbar', width= 0.6)+
#   geom_boxplot(width = 0.6) +
#   expand_limits(y = 0) + 
#   theme_classic() +
#   scale_fill_discrete(name = "Timepoint") +
#   ylab("B-cell score") + xlab(NULL) +
#   theme(legend.position = "right")
# ggsave("../bcell_d0vsd15_3.pdf")


