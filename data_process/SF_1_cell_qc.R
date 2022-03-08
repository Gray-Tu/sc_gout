
setwd("SYN_data")
# SYN data

library(hdf5r)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)

Convert("syn146.h5ad", dest = "h5seurat", overwrite = TRUE)
Convert("syn151.h5ad", dest = "h5seurat", overwrite = TRUE)
Convert("syn183.h5ad", dest = "h5seurat", overwrite = TRUE)

syn146 <- LoadH5Seurat("syn146.h5seurat")
syn151 <- LoadH5Seurat("syn151.h5seurat")
syn183 <- LoadH5Seurat("syn183.h5seurat")

syn146$CASE <- "Gout146"
syn146$UA <- 12.1
syn146$GROUP <- "" # NO IRB
syn146$GENDER <- "female"
syn146$AGE <- 83

syn151$CASE <- "Gout151"
syn151$UA <- 8.2
syn151$GROUP <- "" # NO IRB
syn151$GENDER <- "male"
syn151$AGE <- 33

syn183$CASE <- "Gout183"
syn183$UA <- 7.8
syn183$GROUP <- "acute"
syn183$GENDER <- "male"
syn183$AGE <- 41

samplelist <- list(
  syn146, syn151, syn183
)

get_count_feature_mt <- function(d){
  temp <- d
  nCount = colSums(x = d, slot = "counts")  # nCount_RNA
  nFeature = colSums(x = GetAssayData(object = d, slot = "counts") > 0)
  temp[['nCount_RNA']] <- nCount
  temp[['nFeature_RNA']] <- nFeature
  temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern="^MT-")
  temp2 <- GetAssayData(d, slot="count")
  ncells <- rowSums(x = temp2 > 0)
  temp <- temp[which(ncells >= 5),]
  return(temp)
}

# subset under 500-4500 genes/ 15 mt% / 30000 reads and NormalizeData for cell cycle score
subset_normal_VST <- function(data) {
  data <- subset(data,
                 subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & nCount_RNA < 30000)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
  return (data)
}

Idents(syn146) <- "syn146"
Idents(syn151) <- "syn151"
Idents(syn183) <- "syn183"

VlnPlot(get_count_feature_mt(syn146), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(get_count_feature_mt(syn151), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(get_count_feature_mt(syn183), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#1 remove predicted doublets
samplelist <- lapply(samplelist, function(x) subset(x, predicted_doublets==FALSE))
#2 minic Read10X return
samplelist <- lapply(samplelist, function(x) get_count_feature_mt(x))
#3 subset under 500-4500 genes/ 7.5 mt% / 30000 reads and NormalizeData for cell cycle score
samplelist <- lapply(samplelist, function(x) subset_normal_VST(x))

#4 Cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
samplelist <- lapply(samplelist,
                     function(x) CellCycleScoring(object = x,
                                                  s.features = s.genes,
                                                  g2m.features = g2m.genes,
                                                  set.ident = TRUE))
#5 SCT normalize and regress out of cell cycle stage, mt%, genes, reads
samplelist <- lapply(samplelist,
                     function(x)
                       SCTransform(x,
                                   vars.to.regress = c("S.Score",
                                                       "G2M.Score",
                                                       "percent.mt",
                                                       "nCount_RNA",
                                                       "nFeature_RNA"),
                                   verbose = FALSE)
)
save(samplelist, file="SYNsamplelist.rda")

# For SCT and rPCA
features <- SelectIntegrationFeatures(object.list = samplelist, nfeatures = 3000)
samplelist <- PrepSCTIntegration(object.list = samplelist,
                                 anchor.features = features)
samplelist <- lapply(X = samplelist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# rPCA
anchors <- FindIntegrationAnchors(object.list = samplelist,
                                  anchor.features = features,
                                  normalization.method = "SCT",
                                  reduction = "rpca")
syn.combined <- IntegrateData(anchorset = anchors,
                               normalization.method = "SCT")
save(syn.combined, file = "qcm2_pri_integrated_rpca_syn_combine.Rda")

# dim 50
syn.combined <- RunPCA(syn.combined, npcs = 50, verbose = FALSE)
syn.combined <- RunUMAP(syn.combined, reduction = "pca", dims = 1:50)
syn.combined <- FindNeighbors(syn.combined, reduction = "pca", dims = 1:50)
syn.combined <- FindClusters(syn.combined, resolution = 0.8)

DimPlot(syn.combined, split.by = "CASE")
save(syn.combined, file = "pc50qcm2_ready_rpca_syn_combine.Rda")

# dim 30
syn.combined.30 <- RunPCA(syn.combined, npcs = 30, verbose = FALSE)
syn.combined.30 <- RunUMAP(syn.combined.30, reduction = "pca", dims = 1:30)
syn.combined.30 <- FindNeighbors(syn.combined.30, reduction = "pca", dims = 1:30)
syn.combined.30 <- FindClusters(syn.combined.30, resolution = 0.8)

DimPlot(syn.combined.30, split.by = "CASE")

# dim 20
syn.combined.20 <- RunPCA(syn.combined, npcs = 20, verbose = FALSE)
syn.combined.20 <- RunUMAP(syn.combined.20, reduction = "pca", dims = 1:20)
syn.combined.20 <- FindNeighbors(syn.combined.20, reduction = "pca", dims = 1:20)
syn.combined.20 <- FindClusters(syn.combined.20, resolution = 0.8)

DimPlot(syn.combined.20, split.by = "CASE")
