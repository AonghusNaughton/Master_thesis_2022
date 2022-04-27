# Stacked barplots with numbers

dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")
to_remove <- c(dat.sub[,dat.sub$patient_id=="ALL3" & is.na(dat.sub$dna_cell_type)]$cell_id , dat.sub[,dat.sub$patient_id=="ALL2"]$cell_id)
dat.sub <- dat.sub[, !colnames(dat.sub) %in% to_remove]

df <- tibble(dat.sub$patient_id, dat.sub$cell_id, dat.sub$dna_class_short, dat.sub$timepoint) %>%
  group_by(`dat.sub$patient_id`, .add = T) %>%
  group_by(`dat.sub$timepoint`, .add = T) %>%
  group_by(`dat.sub$dna_class_short`, .add = T) %>%
  tally(name = "n_clone") %>%
  ungroup(`dat.sub$timepoint`)
totals <- df %>%
  group_by(`dat.sub$patient_id`, .add = T) %>%
  group_by(`dat.sub$timepoint`, .add=T) %>%
  summarise(n_timepoint = sum(n_clone))
p <- ggplot(df, aes(x=factor(`dat.sub$timepoint`, levels = c("diagnosis", "day 2", "day 3", "day 5", 
                                                             "day 15", "day 29", "relapse", "relapse 1", 
                                                             "relapse 2")),
                    y=n_clone,
                    fill=`dat.sub$dna_class_short`)) +
  geom_col(position = "fill") +
  theme_classic() +
  facet_wrap(.~`dat.sub$patient_id`, scales = "free_x", nrow = 2) +
  labs(y="Fraction of cells per clone", x="Timepoint", fill="DNA class") +
  geom_text(data = totals, aes(x=`dat.sub$timepoint`,
                               y=n_timepoint,
                               label=n_timepoint,
                               fill=NULL), nudge_y = 10)
  
  
ggsave("/Users/aonghusnaughton/Proj_eng/April22/clonal_fractions.pdf", height = 7, width = 11)
  