counts_list <- list(reads_SKA002 <- read.table("/Users/aonghusnaughton/Proj_eng/data_Solrun/SKA002_counts.tsv",  sep="\t", stringsAsFactors=F, header=T),
                    reads_SKA_5.SKA_7 <- read.table("/Users/aonghusnaughton/Proj_eng/data_Solrun/SKA_5_and_7.tsv",  sep="\t", stringsAsFactors=F, header=T),
                    reads_SKA_metseq6 <- read.table("/Users/aonghusnaughton/Proj_eng/data_Solrun/SKA_metseq6.tsv",  sep="\t", stringsAsFactors=F, header=T),
                    reads_PDAC_RNA_counts_combined <- read.csv("/Users/aonghusnaughton/Proj_eng/data_Solrun/PDAC_RNA_counts_combined.csv",stringsAsFactors=F, header=T),
                    reads_13.14.15 <- read.table("/Users/aonghusnaughton/Proj_eng/data_Solrun/rna_counts.tsv", sep="\t", stringsAsFactors=F, header=T))

counts_list2 <- lapply(counts_list, function(x){
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  colnames(x) <- sub("-","\\.",colnames(x)) # Convert to cell_id
  return(x)
})
counts <- do.call(cbind, counts_list2)

ercc.idx <- grepl("^ERCC-", rownames(counts)) # ERCC spike-ins
counts.f <- counts[!grepl("^__|_",rownames(counts)) & !ercc.idx,] # Throws HTseq rows, ERCC and underscored genes (2)
counts.ercc <- counts[ercc.idx,]
ercc.frac <- colSums(counts.ercc)/(colSums(counts.f)+colSums(counts.ercc))
dim(counts.f)
duplicated_names <- duplicated(colnames(counts.f))
counts.f <- counts.f[!duplicated_names]
counts.f[is.na(counts.f)] <- 0

# Filter
mincount <- 5000
mingenes <- 500

filtered.counts <- counts.f %>% 
  subset(rowSums(.)>0) %>% # Drop genes (rows) with 0 counts in all samples
  subset(select=colSums(.)>mincount) %>% # Drop cells (cols) with total read count < minimum count
  subset(select=colSums(. > 0)>mingenes) %>% # Drop cells with genes below minimum gene (Default 500)
  as.matrix()

ALL_cells <- grep("VZA.*", colnames(filtered.counts))
SKB_cells <- grep("SKB.*", colnames(filtered.counts))
EBA_cells <- grep("EBA.*", colnames(filtered.counts))
filtered.counts <- subset(filtered.counts, select = -c(ALL_cells, SKB_cells, EBA_cells))

dim(filtered.counts)

meta_seurat <- data.frame(cell_id=colnames(filtered.counts))
meta_seurat$plate_id <- substr(meta_seurat$cell_id, 1, 7)
rownames(meta_seurat) <- meta_seurat$cell_id


dat <- CreateSeuratObject(filtered.counts, names.field=1, names.delim="//.", assay="RNA")
dat <- AddMetaData(dat, meta_seurat)
dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
dat <- FindVariableFeatures(dat, selection.method="dispersion", nfeatures = 2000, verbose=F) # Seurat feature finder (vst, dispersion)
dat <- ScaleData(dat, verbose=F, rownames(dat))
dat <- RunPCA(dat, npcs=30, verbose=F)
dat <- FindNeighbors(dat, dims=1:30, verbose=F)
dat <- FindClusters(dat, resolution = 1.5, algorithm=1, n.iter=1000, n.start=100, verbose=F)
dat <- RunUMAP(dat, n.components=2, dims=1:30,verbose=F)

UMAP <- DimPlot(dat, group.by = "plate_id") 

pdf("UMAP_martin_legend.pdf", height = 5, width = 6)
UMAP
dev.off()

