# 1115 update MH/RG
library(hdf5r)
library(Seurat)
library(SeuratDisk)
library(patchwork)

#setwd("CellReport_scrublet")
#setwd("/home2/gray/workspace/20210705_GoutSingleCellPaper/20210819_plan/SINGLE_CELL/CellReport_scrublet")
setwd("/home2/gray/workspace/20210705_GoutSingleCellPaper/20210819_plan/SINGLE_CELL/CellReport_scrublet/Test_Parameters")
#Convert("RJ.h5ad", dest = "h5seurat", overwrite = TRUE)
#Convert("YB.h5ad", dest = "h5seurat", overwrite = TRUE)
#Convert("MH.h5ad", dest = "h5seurat", overwrite = TRUE)
#Convert("RG.h5ad", dest = "h5seurat", overwrite = TRUE)
#RJ <- LoadH5Seurat("../RJ.h5seurat")
YB <- LoadH5Seurat("../YB.h5seurat")
MH <- LoadH5Seurat("../MH.h5seurat")
RG <- LoadH5Seurat("../RG.h5seurat")
N01F <- LoadH5Seurat("/home2/gray/workspace/20211109_AHUS/20211109_SINGLECELL/N01F.h5seurat")
N05S <- LoadH5Seurat("/home2/gray/workspace/20211109_AHUS/20211109_SINGLECELL/N05S.h5seurat")

N01F$CASE <- "N01F"
N05S$CASE <- "N05S"

N01F$GROUP <- "Health"
N05S$GROUP <- "Health"


#for (case in c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')){
#  for (stage in c('a', 'b', 'c')){
#    name = paste(paste(case, stage, sep="_"),"h5ad",sep='.')
#    Convert(name, dest = "h5seurat", overwrite = TRUE)
#  }
#}

c1_a <- LoadH5Seurat("../c1_a.h5seurat")
c1_b <- LoadH5Seurat("../c1_b.h5seurat")
c1_c <- LoadH5Seurat("../c1_c.h5seurat")
c2_a <- LoadH5Seurat("../c2_a.h5seurat")
c2_b <- LoadH5Seurat("../c2_b.h5seurat")
c2_c <- LoadH5Seurat("../c2_c.h5seurat")
c3_a <- LoadH5Seurat("../c3_a.h5seurat")
c3_b <- LoadH5Seurat("../c3_b.h5seurat")
c3_c <- LoadH5Seurat("../c3_c.h5seurat")
c4_a <- LoadH5Seurat("../c4_a.h5seurat")
c4_b <- LoadH5Seurat("../c4_b.h5seurat")
c4_c <- LoadH5Seurat("../c4_c.h5seurat")
c5_a <- LoadH5Seurat("../c5_a.h5seurat")
c5_b <- LoadH5Seurat("../c5_b.h5seurat")
c5_c <- LoadH5Seurat("../c5_c.h5seurat")
c6_a <- LoadH5Seurat("../c6_a.h5seurat")
c6_b <- LoadH5Seurat("../c6_b.h5seurat")
c6_c <- LoadH5Seurat("../c6_c.h5seurat")
c7_a <- LoadH5Seurat("../c7_a.h5seurat")
c7_b <- LoadH5Seurat("../c7_b.h5seurat")
c7_c <- LoadH5Seurat("../c7_c.h5seurat")

c1_a$CASE <- "Case01"
c1_b$CASE <- "Case01"
c1_c$CASE <- "Case01"
c2_a$CASE <- "Case02"
c2_b$CASE <- "Case02"
c2_c$CASE <- "Case02"
c3_a$CASE <- "Case03"
c3_b$CASE <- "Case03"
c3_c$CASE <- "Case03"
c4_a$CASE <- "Case04"
c4_b$CASE <- "Case04"
c4_c$CASE <- "Case04"
c5_a$CASE <- "Case05"
c5_b$CASE <- "Case05"
c5_c$CASE <- "Case05"
c6_a$CASE <- "Case06"
c6_b$CASE <- "Case06"
c6_c$CASE <- "Case06"
c7_a$CASE <- "Case07"
c7_b$CASE <- "Case07"
c7_c$CASE <- "Case07"
#RJ$CASE <- "RJ"
YB$CASE <- "YB"
MH$CASE <- "MH"
RG$CASE <- "RG"

c1_a$GROUP <- "A_acute"
c2_a$GROUP <- "A_acute"
c3_a$GROUP <- "A_acute"
c4_a$GROUP <- "A_acute"
c5_a$GROUP <- "A_acute"
c6_a$GROUP <- "A_acute"
c7_a$GROUP <- "A_acute"
c1_b$GROUP <- "B_week"
c2_b$GROUP <- "B_week"
c3_b$GROUP <- "B_week"
c4_b$GROUP <- "B_week"
c5_b$GROUP <- "B_week"
c6_b$GROUP <- "B_week"
c7_b$GROUP <- "B_week"
c1_c$GROUP <- "C_months"
c2_c$GROUP <- "C_months"
c3_c$GROUP <- "C_months"
c4_c$GROUP <- "C_months"
c5_c$GROUP <- "C_months"
c6_c$GROUP <- "C_months"
c7_c$GROUP <- "C_months"
YB$GROUP <- "Health"
#RJ$GROUP <- "Health"
MH$GROUP <- "Health"
RG$GROUP <- "Health"

samplelist <- list(
    #RJ, 
    MH, RG, N01F, N05S, YB, 
    c1_a, c2_a, c3_a, c4_a, c5_a, c6_a, c7_a,
    c1_b, c2_b, c3_b, c4_b, c5_b, c6_b, c7_b,
    c1_c, c2_c, c3_c, c4_c, c5_c, c6_c, c7_c
  )

# Function of minic Read10X return
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

# subset under 500-3500 genes/ 7.5 mt% / 20000 reads and NormalizeData for cell cycle score
#              500-4500        15  mt% / 30000 reads
subset_normal_VST <- function(data) {
  data <- subset(data,
                 subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & nCount_RNA < 30000)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
  return (data)
}

#1 remove predicted doublets
samplelist <- lapply(samplelist, function(x) subset(x, predicted_doublets==FALSE))
#2 minic Read10X return
samplelist <- lapply(samplelist, function(x) get_count_feature_mt(x))
#3 subset under  genes/  mt% /  reads and NormalizeData for cell cycle score
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
save(samplelist, file="m2h5_samplelist.rda")

# For SCT and rPCA
features <- SelectIntegrationFeatures(object.list = samplelist, nfeatures = 3000)
samplelist <- PrepSCTIntegration(object.list = samplelist,
                                 anchor.features = features)
samplelist <- lapply(X = samplelist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# rPCA and select HC as reference set
anchors <- FindIntegrationAnchors(object.list = samplelist,
                                  anchor.features = features,
                                  normalization.method = "SCT",
                                  reference = c(1, 2, 3, 4, 5),
                                  reduction = "rpca")
save(anchors, file = "qcm2h5_ref_rpca_Gout_anchors.Rda")
Gout.combined <- IntegrateData(anchorset = anchors,
                          normalization.method = "SCT")
save(Gout.combined, file = "qcm2h5_pri_integrated_rpca_Gout_combinePBMC.Rda")

# Note SCT no run ScaleData, dim 30
Gout.combined <- RunPCA(Gout.combined, npcs = 30, verbose = FALSE)
Gout.combined <- RunUMAP(Gout.combined, reduction = "pca", dims = 1:30)
Gout.combined <- FindNeighbors(Gout.combined, reduction = "pca", dims = 1:30)
Gout.combined <- FindClusters(Gout.combined, resolution = 0.8)

save(Gout.combined, file = "pc30qcm2h5_ready_rpca_Gout_combinePBMC.Rda")

                     

