# Relapse gene sets 11 April 
library(Seurat)
library(tidyverse)
library(psych)


dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")

# Subet to only include cases that have diagnostic and relapse samples
dat.sub <- SetIdent(dat.sub, value = "pid_time")
names <- unique(dat.sub$pid_time)
names.rel <- names[grepl(".rel", names)]
names.d0 <- gsub(".rel", ".d0", names.rel)
names.d0 <- unique(gsub("2", "", names.d0))
names.all <- c(names.d0, names.rel)
# names.all <- names.all[!names.all %in% c("MRD_ALL47.d0", "MRD_ALL47.d15")] #only one cell at d15
dat.sub <- subset(dat.sub, idents = names.all[1:length(names.all)])
to_remove <- dat.sub[,dat.sub$patient_id=="ALL3" & is.na(dat.sub$dna_cell_type)]$cell_id
dat.sub <- dat.sub[, !colnames(dat.sub) %in% to_remove]
table(dat.sub$pid_time)
dim(dat.sub)


# To include only highly expressed genes: 
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

dat.sub <- SetIdent(dat.sub, value = "timepoint_short")
markers <- lapply(names, function(x){
  m <- FindMarkers(dat.sub[, dat.sub$patient_id==x], if (x == "ALL3"){
    ident.1 = c("rel", "rel2") } else {
      ident.1 = "rel"
    }, 
    ident.2 = "d0",
    logfc.threshold = 0.5,
    min.pct = 0.25)
  m <- m[m$avg_log2FC > 0 & m$p_val_adj < 0.05,]
  return(m)
})
names(markers) <- lapply(names, function(x) x)

markers$ALL4 <- markers$ALL4[markers$ALL4$avg_log2FC > 1,]

# How do ALL3 and ALL4 relapse genes correlate (seperately)

dat.sub <- readRDS("dat.sub_filtered_10_cells_d0.Rds")
names <- unique(dat.sub$patient_id)
counts <- readRDS("counts_per_pid_10.Rds")
scaled_counts <- readRDS("scale_data_counts_10.Rds")
persister_genes <- readRDS("persister_genes.Rds")
genes <-  unique(c(rownames(markers$ALL4), rownames(markers$ALL3), persister_genes))


counts <- lapply(counts, function(x){
  x <- x[genes,]
  x <- x[rowSums(x)>0,]
  return(x)
})

gene_list_for_scaled <- lapply(counts, function(x){
  rownames(x)
})

for (i in names){
  scaled_counts[[i]] <- scaled_counts[[i]][gene_list_for_scaled[[i]],]
}

correlations <- lapply(scaled_counts, function(x){
  corr.test(as.matrix(t(x)), method = "spearman")
})

correlations_filtered
















correlations <- lapply(scaled_counts, function(x){
  res <- cor(x=as.matrix(t(x[rownames(x)[rownames(x) %in% unique(c(rownames(markers$ALL4), rownames(markers$ALL3), persister_genes))],])),
             method = "spearman")
})

res_Relapse_ALL3<- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), rownames(markers$ALL3)), intersect(colnames(x), rownames(markers$ALL3))]
})
res_Relapse_ALL4 <- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), rownames(markers$ALL4)), intersect(colnames(x), rownames(markers$ALL4))]
})
res_Relapse_ALL3.Relapse_ALL4 <- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), rownames(markers$ALL3)), intersect(colnames(x), rownames(markers$ALL4))]
})
res_Persister <- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), persister_genes), intersect(colnames(x), persister_genes)]
})
res_Persister.Relapse_ALL3 <- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), persister_genes), intersect(colnames(x), rownames(markers$ALL3))]
})
res_Persister.Relapse_ALL4 <- lapply(correlations, function(x){
  df <- x[intersect(rownames(x), persister_genes), intersect(colnames(x), rownames(markers$ALL4))]
})


y <- "Persister"
x <- "Relapse_ALL4"

# Average correlation of y-axis genes to x-axis genes 
assign(paste0(y, "_genes_x_coord"), lapply(eval(parse(text=paste0("res_", y, ".", x))), function(i){
  rowMeans(i)
}))

# Average correlation of y-axis genes to y-axis genes 
assign(paste0(y, "_genes_y_coord"), lapply(eval(parse(text=paste0("res_", y))), function(i){
  rowMeans(i)
}))


# Average correlation of x-axis genes x-axis genes 
assign(paste0(x, "_genes_x_coord"), lapply(eval(parse(text = paste0("res_", x))), function(i){
  rowMeans(i)
}))

# Average correlation of x-axis genes to y-axis genes 
assign(paste0(x, "_genes_y_coord"), lapply(eval(parse(text = paste0("res_", y, ".", x))), function(i){
  colMeans(i)
}))


######################################################################################################################

# Create data structure that holds x and y coordinates for each gene and plot as points on a coordinate plane.


y_gene_coordinates <- vector(mode = "list", length = length(names))
names(y_gene_coordinates) <- lapply(names, function(x) x)

for (i in names){
  y_gene_coordinates[[i]] <- data.frame(eval(parse(text = paste0(y, "_genes_x_coord")))[[i]], 
                                        eval(parse(text = paste0(y, "_genes_y_coord")))[[i]], 
                                        y)
} 

for (i in names){
  colnames(y_gene_coordinates[[i]]) <- c("x", "y", "geneset")
}

x_gene_coordinates <- vector(mode = "list", length = length(names))

for (i in names){
  x_gene_coordinates[[i]] <- data.frame(eval(parse(text = paste0(x, "_genes_x_coord")))[[i]], 
                                        eval(parse(text = paste0(x, "_genes_y_coord")))[[i]], 
                                        x)
} 

for (i in names){
  colnames(x_gene_coordinates[[i]]) <- c("x", "y", "geneset")
}

combined <- list()

for (i in names){
  combined[[i]] <- rbind(x_gene_coordinates[[i]], (y_gene_coordinates[[i]]))
}

plots <- lapply(combined, function(i){
  ggplot(i, aes(x,y, color=geneset)) +
    geom_point() + 
    xlab(paste0("Average correlation with ", x, " genes")) +
    ylab(paste0("Average correlation with ",y, " genes")) +
    ylim(-0.1,0.1) +
    xlim(-0.1,0.1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_classic()
})

lapply(names(plots), function(i){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/Average_correlations_new/two_sets_relapse/Persister_v_RelALL4/", i, ".pdf", sep = ""),
         plot = plots[[i]],
         width = 8,
         height = 5)
})


# correlations <- lapply(correlations, function(x){
#   diag(x) <- NA
# })
# 
# for (i in names){
#   diag(correlations[[i]]) <- NA
# }

