# UMAPs
dat.sub <- readRDS("dat.sub_wNa_feb22.Rds")

dat.sub <- SetIdent(dat.sub, value = "pid_time")
names <- unique(dat.sub$pid_time)
names <- names[!names %in% c("ALL5.d29","MRD_ALL47.d15")]
dat.sub <- subset(dat.sub, idents = names[1:length(names)])

dat.s <- lapply(names, function(x){
  if("ALL3" %in% x){
    d <- list(as.data.frame(dat.sub[,dat.sub$pid_time==x & dat.sub$dna_cell_type=="blasts"]@assays$RNA@counts),
              as.data.frame(dat.sub[,dat.sub$pid_time==x & dat.sub$dna_cell_type=="blasts"]@meta.data))
    names(d) <- c("counts", "meta.data")
  } else {
    d <- list(as.data.frame(dat.sub[,dat.sub$pid_time==x]@assays$RNA@counts),
              as.data.frame(dat.sub[,dat.sub$pid_time==x]@meta.data))
    names(d) <- c("counts", "meta.data")
  }
  return(d)
}) 

names(dat.s) <- lapply(names, function(x) x)


dat.split <- lapply(dat.s, function(x){
  dat <- CreateSeuratObject(x$counts, names.field=1, names.delim="//.", assay="RNA")
  dat <- AddMetaData(dat, x$meta)
  dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
  dat <- FindVariableFeatures(dat, selection.method="dispersion", nfeatures = 2000, verbose=F) # Seurat feature finder (vst, dispersion)
  dat <- ScaleData(dat, verbose=F, rownames(dat))
  dat <- RunPCA(dat, npcs=10, verbose=F)
  dat <- FindNeighbors(dat, dims=1:10, verbose=F)
  dat <- FindClusters(dat, resolution = 0.5, algorithm=1, n.iter=1000, n.start=100, verbose=F)
  dat <- RunUMAP(dat, dims = 1:10, verbose=F)
  return(dat)
})

names(dat.split) <- lapply(names, function(x) x)

UMAPs <- lapply(dat.split, function(x){
  DimPlot(x, group.by = "dna_best_class")
})

for (i in names(UMAPs)){
  ggsave(filename = paste("/Users/aonghusnaughton/Proj_eng/April22/UMAPs/", i, ".pdf", sep = ""),
         plot = UMAPs[[i]],
         width = 5,
         height = 5)
}

