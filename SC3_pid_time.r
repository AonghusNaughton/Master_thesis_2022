load("ALL_data_Feb22.Rda")
library(Seurat)
library(SC3)
library(SingleCellExperiment)

dat.sub_wNa <- dat[, (dat$dna_cell_type=="blasts" | is.na(dat$dna_cell_type)) & 
                     dat$rna_cell_type=="blasts" & 
                     dat$cell_phase=="G1"]

dat.sub <-dat[, dat$dna_cell_type=="blasts" & 
                dat$rna_cell_type=="blasts" &
                dat$cell_phase=="G1"]

# dat.sub_wNa_annot@meta.data <- within(dat.sub@meta.data, dna_best_class[is.na(dna_best_class)] <- paste(patient_id[is.na(dna_best_class)], "NA", sep = "_"))

names <- unique(dat.sub_wNa$pid_time)

sce <- lapply(names, function(x){
  r <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(dat.sub_wNa[, dat.sub_wNa$pid_time==x]@assays$RNA@counts),
      logcounts = log2(as.matrix(dat.sub_wNa[, dat.sub_wNa$pid_time==x]@assays$RNA@counts) + 1)
    ),
    colData = dat.sub_wNa[, dat.sub_wNa$pid_time==x]@meta.data, 
  )
})

names(sce) <- lapply(names, function(x){
  x
})

# Adds rowData$feature_symbol that is required in SC3

for (i in names){
  rowData(sce[[i]])$feature_symbol <- rownames(sce[[i]])
}

# Remove MRD_ALL47.d15 as too few cells for clustering
sce[["MRD_ALL47.d15"]] = NULL


# Estimate optimal k for all and append to object
# sce <- lapply(sce, function(x){
#   sc3_estimate_k(x)
# })
# 
# # View value estimated fort each
# lapply(sce, function(x){
#   str(metadata(x)$sc3)
# })

# run sc3 for each patient x time and calculate 2-12 clusters for each
sc3_split_wNa <- lapply(sce, function(x){
  c <- sc3(x, ks = 2:7, biology = T,)
  return(c)
})

