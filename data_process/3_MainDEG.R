# DEG of main cell tyep cross SYN/Acute

setwd("With_SYN")

library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(Seurat)
library(dplyr)
# load data ----
load(file="addbpe_pc50m2h5_InfScore_rpca_Gout_combinePBMC.Rda")
Acute_HC_set <- subset(Gout.combined.50, GROUP %in% c("A_acute", "Health"))

#syn.combined
load(file="pc50_syn.rda")

# integrate dataset
DefaultAssay(Acute_HC_set) <- "RNA"
DefaultAssay(syn.combined) <- "RNA"

samplelist = list(Acute_HC_set, syn.combined)
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
features <- SelectIntegrationFeatures(object.list = samplelist, nfeatures = 3000)
samplelist <- PrepSCTIntegration(object.list = samplelist,
                                 anchor.features = features)
samplelist <- lapply(X = samplelist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = samplelist,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = "rpca")
pbmc.syn.combined <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT")
save(pbmc.syn.combined, file = "integrated_rpca_syn_PBMCAcuteHealth.Rda")


pbmc.syn.combined <- RunPCA(pbmc.syn.combined, npcs = 50, verbose = FALSE)
pbmc.syn.combined <- RunUMAP(pbmc.syn.combined, reduction = "pca", dims = 1:50)
pbmc.syn.combined <- FindNeighbors(pbmc.syn.combined, reduction = "pca", dims = 1:50)
pbmc.syn.combined <- FindClusters(pbmc.syn.combined, resolution = 0.8)

my_DF <- pbmc.syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(GROUP = ifelse(GROUP %in% c("acute", ""),
                            "Synovial fluid", 
                            GROUP))

pbmc.syn.combined$GROUP <- my_DF$GROUP

pbmc.syn.combined$cell_type <- factor(pbmc.syn.combined$cell_type, 
                     levels=c("Naive CD4 T cells",
                              "Follicular helper T cells",
                              "T regulatory cells",
                              "Terminal effector CD4 T cells",
                              "Th1 cells",
                              "Th1/Th17 cells",
                              "Th17 cells",
                              "Th2 cells",
                              "Naive CD8 T cells",
                              "Central memory CD8 T cells",
                              "Effector memory CD8 T cells",
                              "MAIT cells",
                              "Terminal effector CD8 T cells",
                              "Non-Vd2 gd T cells",
                              "Vd2 gd T cells",
                              "Natural killer cells",
                              "Classical monocytes",
                              "Intermediate monocytes",
                              "Non classical monocytes",
                              "Macrophages",
                              "Macrophages M1",
                              "INF-Macro",
                              "lung Macro",
                              "Myeloid dendritic cells",
                              "Plasmacytoid dendritic cells",
                              "Naive B cells",
                              "Non-switched memory B cells",
                              "Switched memory B cells",
                              "Plasmablasts",
                              "Exhausted B cells",
                              "Low-density basophils",
                              "Low-density neutrophils",
                              "Megakaryoctes",
                              "Progenitor cells"))

save(pbmc.syn.combined, file = "integrated_rpca_syn_PBMCAcuteHealth.Rda")

p1 <- DimPlot(pbmc.syn.combined, label=T, repel = T) + NoLegend()
p2 <- DimPlot(pbmc.syn.combined, group.by = "cell_type", 
              cols=DiscretePalette(35, "polychrome"),
              label.box = T,
              label=T, repel=T)
p3 <- DimPlot(pbmc.syn.combined, split.by = "GROUP", group.by = "cell_type",
              label=T, repel=T, cols=DiscretePalette(36, "polychrome"))
ggsave("temp2.pdf", plot=p1 / p2, width = 16, height = 16)
p3
ggsave("temp.pdf", 
       plot=p3+theme(legend.position="bottom"), 
       width = 14, 
       height = 10)

# MonoMaDC----
# Classical monocytes
CD14Mono <- subset(pbmc.syn.combined, cell_type=="Classical monocytes")
DimPlot(CD14Mono, split.by="GROUP")
CD14Mono <- RunPCA(CD14Mono, npcs = 50, verbose = FALSE)
CD14Mono <- RunUMAP(CD14Mono, reduction = "pca", dims = 1:50)
CD14Mono <- FindNeighbors(CD14Mono, reduction = "pca", dims = 1:50)
CD14Mono <- FindClusters(CD14Mono, resolution = 0.2)
DimPlot(CD14Mono, split.by="GROUP")
DimPlot(CD14Mono, group.by="GROUP")
#DEG_CD14Mono <- FindAllMarkers(CD14Mono)
#write.csv(DEG_CD14Mono, file="DEG_CD14Mono.csv")

# Intermediate monocytes
InterMono <- subset(pbmc.syn.combined, cell_type=="Intermediate monocytes")
DimPlot(InterMono, split.by="GROUP")
InterMono <- RunPCA(InterMono, npcs = 50, verbose = FALSE)
InterMono <- RunUMAP(InterMono, reduction = "pca", dims = 1:50)
InterMono <- FindNeighbors(InterMono, reduction = "pca", dims = 1:50)
InterMono <- FindClusters(InterMono, resolution = 0.6)
DimPlot(InterMono, split.by="GROUP")
DimPlot(InterMono, group.by="GROUP")
#DEG_InterMono <- FindAllMarkers(InterMono)
#write.csv(DEG_InterMono, file="DEG_InterMono.csv")


# Plasmacytoid dendritic cells
pDC <- subset(pbmc.syn.combined, cell_type=="Plasmacytoid dendritic cells")
DimPlot(pDC, split.by="GROUP")
pDC <- RunPCA(pDC, npcs = 50, verbose = FALSE)
pDC <- RunUMAP(pDC, reduction = "pca", dims = 1:50)
pDC <- FindNeighbors(pDC, reduction = "pca", dims = 1:50)
pDC <- FindClusters(pDC, resolution = 0.4)
DimPlot(pDC, split.by="GROUP")
DimPlot(pDC, group.by="GROUP")
#DEG_pDC <- FindAllMarkers(pDC)
#write.csv(DEG_pDC, file="DEG_pDC.csv")

# Myeloid dendritic cells
mDC <- subset(pbmc.syn.combined, cell_type=="Myeloid dendritic cells")
DimPlot(mDC, split.by="GROUP")
mDC <- RunPCA(mDC, npcs = 50, verbose = FALSE)
mDC <- RunUMAP(mDC, reduction = "pca", dims = 1:50)
mDC <- FindNeighbors(mDC, reduction = "pca", dims = 1:50)
mDC <- FindClusters(mDC, resolution = 0.2)
DimPlot(mDC, split.by="GROUP")
DimPlot(mDC, group.by="GROUP")
#DEG_mDC <- FindAllMarkers(mDC)
#write.csv(DEG_mDC, file="DEG_mDC.csv")

# T----
# Effector memory CD8 T cells
CD8Tem <- subset(pbmc.syn.combined, cell_type=="Effector memory CD8 T cells")
DimPlot(CD8Tem, split.by="GROUP")
CD8Tem <- RunPCA(CD8Tem, npcs = 50, verbose = FALSE)
CD8Tem <- RunUMAP(CD8Tem, reduction = "pca", dims = 1:50)
CD8Tem <- FindNeighbors(CD8Tem, reduction = "pca", dims = 1:50)
CD8Tem <- FindClusters(CD8Tem, resolution = 0.8)
DimPlot(CD8Tem, split.by="GROUP")
DimPlot(CD8Tem, group.by="GROUP")
#DEG_CD8Tem <- FindAllMarkers(CD8Tem)
#write.csv(DEG_CD8Tem, file="DEG_CD8Tem.csv")

# MAIT cells
MAIT <- subset(pbmc.syn.combined, cell_type=="MAIT cells")
DimPlot(MAIT, split.by="GROUP")
MAIT <- RunPCA(MAIT, npcs = 50, verbose = FALSE)
MAIT <- RunUMAP(MAIT, reduction = "pca", dims = 1:50)
MAIT <- FindNeighbors(MAIT, reduction = "pca", dims = 1:50)
MAIT <- FindClusters(MAIT, resolution = 0.5)
DimPlot(MAIT, split.by="GROUP")
DimPlot(MAIT, group.by="GROUP")
#DEG_MAIT <- FindAllMarkers(MAIT)
#write.csv(DEG_MAIT, file="DEG_MAIT.csv")

# Naive CD4 T cells
CELL = "Naive CD4 T cells"
NaiveCD4T <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(NaiveCD4T, split.by="GROUP")
NaiveCD4T <- RunPCA(NaiveCD4T, npcs = 50, verbose = FALSE)
NaiveCD4T <- RunUMAP(NaiveCD4T, reduction = "pca", dims = 1:50)
NaiveCD4T <- FindNeighbors(NaiveCD4T, reduction = "pca", dims = 1:50)
NaiveCD4T <- FindClusters(NaiveCD4T, resolution = 0.2)
DimPlot(NaiveCD4T, split.by="GROUP")
DimPlot(NaiveCD4T, group.by="GROUP")
#DEG_NaiveCD4T <- FindAllMarkers(NaiveCD4T)
#write.csv(DEG_NaiveCD4T, file="DEG_NaiveCD4T.csv")

# T regulatory cells
CELL = "T regulatory cells"
Treg <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Treg, split.by="GROUP")
Treg <- RunPCA(Treg, npcs = 50, verbose = FALSE)
Treg <- RunUMAP(Treg, reduction = "pca", dims = 1:50)
Treg <- FindNeighbors(Treg, reduction = "pca", dims = 1:50)
Treg <- FindClusters(Treg, resolution = 0.4)
DimPlot(Treg, split.by="GROUP")
DimPlot(Treg, group.by="GROUP")
#DEG_Treg <- FindAllMarkers(Treg)
#write.csv(DEG_Treg, file="DEG_Treg.csv")

# Th1/Th17 cells
CELL = "Th1/Th17 cells"
Th1_Th17 <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Th1_Th17, split.by="GROUP")
Th1_Th17 <- RunPCA(Th1_Th17, npcs = 50, verbose = FALSE)
Th1_Th17 <- RunUMAP(Th1_Th17, reduction = "pca", dims = 1:50)
Th1_Th17 <- FindNeighbors(Th1_Th17, reduction = "pca", dims = 1:50)
Th1_Th17 <- FindClusters(Th1_Th17, resolution = 0.5)
DimPlot(Th1_Th17, split.by="GROUP")
DimPlot(Th1_Th17, group.by="GROUP")
#DEG_Th1_Th17 <- FindAllMarkers(Th1_Th17)
#write.csv(DEG_Th1_Th17, file="DEG_Th1_Th17.csv")

# Vd2 gd T cells
CELL = "Vd2 gd T cells"
vd2gdT <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(vd2gdT, split.by="GROUP")
vd2gdT <- RunPCA(vd2gdT, npcs = 50, verbose = FALSE)
vd2gdT <- RunUMAP(vd2gdT, reduction = "pca", dims = 1:50)
vd2gdT <- FindNeighbors(vd2gdT, reduction = "pca", dims = 1:50)
vd2gdT <- FindClusters(vd2gdT, resolution = 0.1)
DimPlot(vd2gdT, split.by="GROUP")
DimPlot(vd2gdT, group.by="GROUP")
#DEG_vd2gdT <- FindAllMarkers(vd2gdT)
#write.csv(DEG_vd2gdT, file="DEG_vd2gdT.csv")


# RNA assay DEG ----
# In DE analysis, must set assay to "RNA"
DefaultAssay(CD14Mono) <- "RNA"
RNADEG_CD14Mono <- FindAllMarkers(CD14Mono)
write.csv(RNADEG_CD14Mono, file="RNADEG_CD14Mono.csv")

DefaultAssay(InterMono) <- "RNA"
RNADEG_InterMono <- FindAllMarkers(InterMono)
write.csv(RNADEG_InterMono, file="RNADEG_InterMono.csv")

DefaultAssay(pDC) <- "RNA"
RNADEG_pDC <- FindAllMarkers(pDC)
write.csv(RNADEG_pDC, file="RNADEG_pDC.csv")

DefaultAssay(pDC) <- "RNA"
RNADEG_mDC <- FindAllMarkers(mDC)
write.csv(RNADEG_mDC, file="RNADEG_mDC.csv")

DefaultAssay(CD8Tem) <- "RNA"
RNADEG_CD8Tem <- FindAllMarkers(CD8Tem)
write.csv(RNADEG_CD8Tem, file="RNADEG_CD8Tem.csv")

DefaultAssay(MAIT) <- "RNA"
RNADEG_MAIT <- FindAllMarkers(MAIT)
write.csv(RNADEG_MAIT, file="RNADEG_MAIT.csv")

DefaultAssay(NaiveCD4T) <- "RNA"
RNADEG_NaiveCD4T <- FindAllMarkers(NaiveCD4T)
write.csv(RNADEG_NaiveCD4T, file="RNADEG_NaiveCD4T.csv")

DefaultAssay(Treg) <- "RNA"
RNADEG_Treg <- FindAllMarkers(Treg)
write.csv(RNADEG_Treg, file="RNADEG_Treg.csv")

DefaultAssay(Th1_Th17) <- "RNA"
RNADEG_Th1_Th17 <- FindAllMarkers(Th1_Th17)
write.csv(RNADEG_Th1_Th17, file="RNADEG_Th1_Th17.csv")

DefaultAssay(vd2gdT) <- "RNA"
RNADEG_vd2gdT <- FindAllMarkers(vd2gdT)
write.csv(RNADEG_vd2gdT, file="RNADEG_vd2gdT.csv")


save(CD14Mono,
     InterMono,
     pDC,
     mDC,
     CD8Tem,
     MAIT,
     NaiveCD4T,
     Treg,
     Th1_Th17,
     vd2gdT,
     file = "MainDEGusage.rda"
     )

#---
load(file="MainDEGusage.rda")
library(pheatmap)
get_all_marker_topN <- function(markers, topN){
  temp <-  markers %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC)
  return(temp)
}


# GO BP
library(DOSE)
library(pathview)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)


# get all deg ----
DEG_CD14Mono <- read.csv("RNADEG_CD14Mono.csv") %>% 
  mutate(cell_type="Classical monocytes")
DEG_mDC      <- read.csv("RNADEG_mDC.csv") %>% 
  mutate(cell_type="Myeloid dendritic cells")
DEG_pDC      <- read.csv("RNADEG_pDC.csv") %>% 
  mutate(cell_type="Plasmacytoid dendritic cells")
DEG_InterMono<- read.csv("RNADEG_InterMono.csv") %>% 
  mutate(cell_type="Intermediate monocytes")
DEG_vd2gdT   <- read.csv("RNADEG_vd2gdT.csv") %>% 
  mutate(cell_type="Vd2 gd T cells")
DEG_Th1_Th17 <- read.csv("RNADEG_Th1_Th17.csv") %>% 
  mutate(cell_type="Th1/Th17 cells")
DEG_Treg     <- read.csv("RNADEG_Treg.csv") %>% 
  mutate(cell_type="T regulatory cells")
DEG_MAIT     <- read.csv("RNADEG_MAIT.csv") %>% 
  mutate(cell_type="MAIT cells")
DEG_CD8Tem   <- read.csv("RNADEG_CD8Tem.csv") %>% 
  mutate(cell_type="Effector memory CD8 T cells")

DEG_total <- rbind(DEG_CD14Mono,
      DEG_mDC,
      DEG_pDC,
      DEG_InterMono,
      DEG_vd2gdT,
      DEG_Th1_Th17,
      DEG_Treg,
      DEG_MAIT,
      DEG_CD8Tem
      )

DEG_total$UPDN <- ifelse(DEG_total$p_val_adj < 0.05 & DEG_total$avg_log2FC > 0.25, 'UP', 
       ifelse(DEG_total$p_val_adj < 0.05 & DEG_total$avg_log2FC < -0.25, "DN", ""))


write.csv(unique(DEG_total$gene), file="DEgene_names.csv")
# to python, convert ENSGID
ENSGID_df <- read.csv("mapping_ENSGID.csv")
# get entriz id
mapping_id <- bitr(ENSGID_df$gene_id, 
                   fromType="ENSEMBL", 
                   toType=c("SYMBOL", "ENTREZID"),
                   OrgDb=org.Hs.eg.db,
                   drop = FALSE)
# match entriz id
ENSGID_df$ENTREZID <- mapping_id[match(ENSGID_df$gene_id, 
                                       mapping_id$ENSEMBL),]$ENTREZID
final_map_df <- ENSGID_df[!is.na(ENSGID_df$ENTREZID),]
# map to DEG_total
DEG_total$ENTREZID <- final_map_df[match(DEG_total$gene, final_map_df$x), ]$ENTREZID

# GOBP
GOBP_bigdat <- data.frame()
for (cellbp in unique(DEG_total$cell_type)){
  rm(cell_set)
  cell_set <- subset(DEG_total, cell_type==cellbp)
  for (cls in unique(cell_set$cluster)){
    cat(cellbp, cls)
    cls_set <- subset(cell_set, cluster==cls)
     # up
    rm(up_set)
    up_set <- subset(cls_set, UPDN=="UP")
    bp_up <- try(enrichGO(gene = up_set$ENTREZID,
                          OrgDb = "org.Hs.eg.db", 
                          ont = "BP",
                          pvalueCutoff = 0.05))
    if ('try-error' %in% class(bp_up) | 'NULL' %in% class(bp_up)){
      print(paste("Error", cellbp, cls, "up", sep=" "))
      } 
    else {
      bp_up <- setReadable(bp_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      bp_up.result <- bp_up@result
      bp_up.result$group <- paste(cellbp, cls,"up", sep="_")
      bp_up.result$cell <- cellbp
      bp_up.result$cluster <- cls
      GOBP_bigdat <- rbind(GOBP_bigdat, bp_up.result)
      }
     # dn
    rm(dn_set)
    dn_set <- subset(cls_set, UPDN=="DN")
    bp_dn <- try(enrichGO(gene = dn_set$ENTREZID,
                          OrgDb = "org.Hs.eg.db",
                          ont = "BP",
                          pvalueCutoff = 0.05))
    if ('try-error' %in% class(bp_dn) | 'NULL' %in% class(bp_dn)){
      print(paste("Error", cellbp, cls, "dn", sep=" "))
      }
    else {
      bp_dn <- setReadable(bp_dn, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      bp_dn.result <- bp_dn@result
      bp_dn.result$group <- paste(cellbp, cls, "dn", sep="_")
      bp_dn.result$cell <- cellbp
      bp_dn.result$cluster <- cls
      GOBP_bigdat <- rbind(GOBP_bigdat, bp_dn.result)
      }
  }
}
write.csv(GOBP_bigdat, file="GOBP_deg_subcluster.csv")


# kegg
kegg_bigdat <- data.frame()
for (cellkk in unique(DEG_total$cell_type)){
  rm(cell_set)
  cell_set <- subset(DEG_total, cell_type==cellkk)
  for (cls in unique(cell_set$cluster)){
    cat(cellbp, cls)
    cls_set <- subset(cell_set, cluster==cls)
    # up
    rm(up_set)
    up_set <- subset(cls_set, UPDN=="UP")
    kk_up <- try(enrichKEGG(gene         = up_set$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.05))
    if ('try-error' %in% class(kk_up) | 'NULL' %in% class(kk_up)){
      print(paste("Error", cellkk, cls, "up", sep=" "))
    } 
    else {
      kk_up <- setReadable(kk_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      kk_up.result <- kk_up@result
      kk_up.result$group <- "up"
      kk_up.result$cell <- cellkk
      kk_up.result$cluster <- cls
      kegg_bigdat <- rbind(kegg_bigdat, kk_up.result)
    }
    # dn
    rm(dn_set)
    dn_set <- subset(cls_set, UPDN=="DN")
    kk_dn <- try(enrichKEGG(gene         = dn_set$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.05))
    if ('try-error' %in% class(kk_dn) | 'NULL' %in% class(kk_dn)){
      print(paste("Error", cellkk, cls, "dn", sep=" "))
    }
    else {
      kk_dn <- setReadable(kk_dn, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      kk_dn.result <- kk_dn@result
      kk_dn.result$group <- "dn"
      kk_dn.result$cell <- cellkk
      kk_dn.result$cluster <- cls
      kegg_bigdat <- rbind(kegg_bigdat, kk_dn.result)
    }
  }
}
write.csv(kegg_bigdat, file="kegg_deg_subcluster.csv")

# clusterProfiler ----

# CD14 Mono----

make_kegg_gobp_func <- function(cell, deg_df){
  up_gene <- subset(deg_df, cell_type==cell & UPDN=="UP")
  dn_gene <- subset(DEG_total, cell_type==cell & UPDN=="DN")
  cell <- str_replace(cell, "/", ".")
  
  up_glist <- list()
  for (cls in unique(up_gene$cluster)){
    cat(cls)
    cls_set <- subset(up_gene, cluster==cls)
    temp <- cls_set[!is.na(cls_set$ENTREZID),]$ENTREZID
    up_glist[as.character(cls)] <- list(temp)
  }
  dn_glist <- list()
  for (cls in unique(dn_gene$cluster)){
    cat(cls)
    cls_set <- subset(dn_gene, cluster==cls)
    temp <- cls_set[!is.na(cls_set$ENTREZID),]$ENTREZID
    dn_glist[as.character(cls)] <- list(temp)
  }
  cat(cell)
  print(" end gene list")
  up_ck <- compareCluster(geneCluster = up_glist, fun = enrichKEGG)
  up_ck <- setReadable(up_ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(up_ck@compareClusterResult,
            file=paste(cell, "up_kegg_compareCluster.csv"))
  
  dn_ck <- compareCluster(geneCluster = dn_glist, fun = enrichKEGG)
  dn_ck <- setReadable(dn_ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(dn_ck@compareClusterResult,
            file=paste(cell, "dn_kegg_compareCluster.csv"))
  cat(cell)
  print(" end kegg csv")
  up_bp <- compareCluster(geneCluster = up_glist,
                          OrgDb = org.Hs.eg.db,
                          fun = enrichGO, 
                          ont = "BP",
                          readable = TRUE)
  write.csv(up_bp@compareClusterResult,
            file=paste(cell, "up_GOBP_compareCluster.csv"))
  
  dn_bp <- compareCluster(geneCluster = dn_glist,
                          OrgDb = org.Hs.eg.db,
                          fun = enrichGO, 
                          ont = "BP",
                          readable = TRUE)
  write.csv(dn_bp@compareClusterResult,
            file=paste(cell, "dn_GOBP_compareCluster.csv"))
  cat(cell)
  print(" end GOBP csv")
  p_upck <- dotplot(up_ck, showCategory=10)
  p_dnck <- dotplot(dn_ck, showCategory=10)
  p_upbp <- dotplot(up_bp, showCategory=10)
  p_dnbp <- dotplot(dn_bp, showCategory=10)
  p_cnet_upck <- cnetplot(up_ck, showCategory=10)
  p_cnet_dnck <- cnetplot(dn_ck, showCategory=10)
  p_cnet_upbp <- cnetplot(up_bp, showCategory=10)
  p_cnet_dnbp <- cnetplot(dn_bp, showCategory=10)
  
  ggsave(paste("UP_kegg",cell, ".png"), 
         plot_grid(p_upck, p_cnet_upck, ncol=1, rel_heights = c(1.5, 1))+
           theme(plot.background = element_rect(fill = "white", colour = "white")), 
         width = 12, 
         height = 16)
  ggsave(paste("DN_kegg",cell, ".png"), 
         plot_grid(p_dnck, p_cnet_dnck, ncol=1, rel_heights = c(1.5, 1))+
           theme(plot.background = element_rect(fill = "white", colour = "white")), 
         width = 12, 
         height = 16)
  ggsave(paste("UP_GOBP",cell, ".png"), 
         plot_grid(p_upbp, p_cnet_upbp, ncol=1, rel_heights = c(1.5, 1))+
           theme(plot.background = element_rect(fill = "white", colour = "white")), 
         width = 12, 
         height = 16)
  ggsave(paste("DN_GOBP",cell, ".png"), 
         plot_grid(p_dnbp, p_cnet_dnbp, ncol=1, rel_heights = c(1.5, 1))+
           theme(plot.background = element_rect(fill = "white", colour = "white")), 
         width = 12, 
         height = 16)
  cat(cell)
  print(" end PLOT")
  returnlist <- list("upck" = up_ck,
                "dnck" = dn_ck,
                "upbp" = up_bp,
                "dnbp" = dn_bp)
  return(returnlist)
}

# [1] "Classical monocytes"          "Myeloid dendritic cells"     
# [3] "Plasmacytoid dendritic cells" "Intermediate monocytes"      
# [5] "Vd2 gd T cells"               "Th1/Th17 cells"              
# [7] "T regulatory cells"           "MAIT cells"                  
# [9] "Effector memory CD8 T cells" 

cell <- "Classical monocytes"
CD14Mono_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Myeloid dendritic cells"
mDC_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Plasmacytoid dendritic cells"
pDC_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Intermediate monocytes"
InterMono_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Vd2 gd T cells"
vd2gdT_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Th1/Th17 cells"
Th1Th17_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "T regulatory cells"
Treg_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "MAIT cells"                  
MAIT_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Effector memory CD8 T cells" 
CD8Tem_funList <- make_kegg_gobp_func(cell, DEG_total)

#cell <- "Classical monocytes"
#up_gene <- subset(DEG_total, cell_type==cell & UPDN=="UP")
#dn_gene <- subset(DEG_total, cell_type==cell & UPDN=="DN")
# up_glist <- list()
# for (cls in unique(up_gene$cluster)){
#   cat(cls)
#   cls_set <- subset(up_gene, cluster==cls)
#   temp <- cls_set[!is.na(cls_set$ENTREZID),]$ENTREZID
#   up_glist[as.character(cls)] <- list(temp)
# }
# dn_glist <- list()
# for (cls in unique(dn_gene$cluster)){
#   cat(cls)
#   cls_set <- subset(dn_gene, cluster==cls)
#   temp <- cls_set[!is.na(cls_set$ENTREZID),]$ENTREZID
#   dn_glist[as.character(cls)] <- list(temp)
# }
## kegg----
# up_ck <- compareCluster(geneCluster = up_glist, fun = enrichKEGG)
# up_ck <- setReadable(up_ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# write.csv(up_ck@compareClusterResult,
#           file=paste(cell, "up_kegg_compareCluster.csv"))
# 
# dn_ck <- compareCluster(geneCluster = dn_glist, fun = enrichKEGG)
# dn_ck <- setReadable(dn_ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# write.csv(dn_ck@compareClusterResult,
#           file=paste(cell, "dn_kegg_compareCluster.csv"))
## GOBP----
# up_bp <- compareCluster(geneCluster = up_glist,
#                         OrgDb = org.Hs.eg.db,
#                         fun = enrichGO, 
#                         ont = "BP",
#                         readable = TRUE)
# write.csv(up_bp@compareClusterResult,
#           file=paste(cell, "up_GOBP_compareCluster.csv"))
# 
# dn_bp <- compareCluster(geneCluster = dn_glist,
#                         OrgDb = org.Hs.eg.db,
#                         fun = enrichGO, 
#                         ont = "BP",
#                         readable = TRUE)
# write.csv(dn_bp@compareClusterResult,
#           file=paste(cell, "dn_GOBP_compareCluster.csv"))
# 
# p_upck_CD14Mono <- dotplot(up_ck, showCategory=10)
# p_dnck_CD14Mono <- dotplot(dn_ck, showCategory=10)
# p_upbp_CD14Mono <- dotplot(up_bp, showCategory=10)
# p_dnbp_CD14Mono <- dotplot(dn_bp, showCategory=10)
# p_cnet_upck_CD14Mono <- cnetplot(up_ck, showCategory=10)
# p_cnet_dnck_CD14Mono <- cnetplot(dn_ck, showCategory=10)
# p_cnet_upbp_CD14Mono <- cnetplot(up_bp, showCategory=10)
# p_cnet_dnbp_CD14Mono <- cnetplot(dn_bp, showCategory=10)
# 
# ggsave(paste("UP_kegg",cell, ".png"), 
#        plot_grid(p_upck_CD14Mono, p_cnet_upck_CD14Mono, ncol=1, rel_heights = c(1.5, 1))+
#          theme(plot.background = element_rect(fill = "white")), 
#        width = 12, 
#        height = 16)
# ggsave(paste("DN_kegg",cell, ".png"), 
#        plot_grid(p_dnck_CD14Mono, p_cnet_dnck_CD14Mono, ncol=1, rel_heights = c(1.5, 1))+
#          theme(plot.background = element_rect(fill = "white")), 
#        width = 12, 
#        height = 16)
# ggsave(paste("UP_GOBP",cell, ".png"), 
#        plot_grid(p_upbp_CD14Mono, p_cnet_upbp_CD14Mono, ncol=1, rel_heights = c(1.5, 1))+
#          theme(plot.background = element_rect(fill = "white")), 
#        width = 12, 
#        height = 16)
# ggsave(paste("DN_GOBP",cell, ".png"), 
#        plot_grid(p_dnbp_CD14Mono, p_cnet_dnbp_CD14Mono, ncol=1, rel_heights = c(1.5, 1))+
#          theme(plot.background = element_rect(fill = "white")), 
#        width = 12, 
#        height = 16)

# ggsave(paste("temp.png"),
#        plot_grid(p_upck_CD14Mono, p_cnet_upck_CD14Mono, ncol=1, rel_heights = c(1.5, 1))+
#                 theme(#plot.background = element_rect(fill = "white"),
#                       plot.background = element_rect(colour = "white", fill="white")),#, size=5)),
#        width = 12,
#        height = 16)


# DimPlot/BarPlot/Heatmap ----

library(ggplotify)

make_subcls_plot <- function(cell, DEG_df, seurat_obj){
  pos_DEG <- subset(DEG_df, cell_type==cell & UPDN=="UP")
  cell <- str_replace(cell, "/", ".")
  seurat_obj$GROUP <- factor(seurat_obj$GROUP, levels=c("Synovial fluid", 
                                                        "A_acute",
                                                        "Health"))
  DefaultAssay(seurat_obj) <- "RNA"
  DF <- seurat_obj[[]] %>% 
    group_by(CASE, GROUP, seurat_clusters) %>% 
    summarise(cell_number = n()) %>%
    mutate(total_number = sum(cell_number)) %>%
    mutate(percent = cell_number/total_number*100)
  t10 <- get_all_marker_topN(pos_DEG, 10)
  exp.avg <-  AverageExpression(seurat_obj,
                                assays = "RNA", 
                                return.seurat = TRUE,
                                slot = "data")
  mat <- GetAssayData(exp.avg, slot = "scale.data")
  mat <- as.matrix(mat[unique(t10$gene), ])
  
  p1_dim <- DimPlot(seurat_obj, split.by = "GROUP") +
    theme(plot.title = element_text(size = 5),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10))
  p2_rel_bar <- ggplot(DF, 
                       aes(x = CASE, y = percent, fill = seurat_clusters))+
    geom_bar(stat = "identity")+
    facet_grid(~GROUP, scales="free", space="free_x")+
    theme(axis.text.x=element_text(angle=90, hjust=1),
          legend.justification = "left")+
    labs(fill='subclusters', x="case")
  p3_heatmap <- pheatmap(t(mat),
                         cellheight = 12,
                         cellwidth = 12,
                         cluster_cols = T,
                         cluster_rows = T,
                         show_colnames = T,
                         fontsize = 8,
                         angle_col = 45)
  g3 = as.ggplot(p3_heatmap)
  p <- plot_grid(p1_dim, p2_rel_bar, g3, ncol=1, 
            labels=c("A", "B", "C"), 
            align = "v",
            axis="lr")+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  ggsave(paste("subcluster_", cell, ".pdf", sep=""), width = 12, height = 10)
  ggsave(paste("subcluster_", cell, ".png", sep=""), width = 12, height = 10)
  returnlist <- list("mat"=mat, "DF"=DF, "p"=p)
  return(returnlist)
}

cell <- "Classical monocytes"
CD14Mono_subList <- make_subcls_plot(cell, DEG_total, CD14Mono)
cell <- "Myeloid dendritic cells"
mDC_subList <- make_subcls_plot(cell, DEG_total, mDC)
cell <- "Plasmacytoid dendritic cells"
pDC_subList <- make_subcls_plot(cell, DEG_total, pDC)
cell <- "Intermediate monocytes"
InterMono_subList <- make_subcls_plot(cell, DEG_total, InterMono)
cell <- "Vd2 gd T cells"
vd2gdT_subList <- make_subcls_plot(cell, DEG_total, vd2gdT)
cell <- "Th1/Th17 cells"
Th1Th17_subList <- make_subcls_plot(cell, DEG_total, Th1_Th17)
cell <- "T regulatory cells"
Treg_subList <- make_subcls_plot(cell, DEG_total, Treg)
cell <- "MAIT cells"                  
MAIT_subList <- make_subcls_plot(cell, DEG_total, MAIT)
cell <- "Effector memory CD8 T cells" 
CD8Tem_subList <- make_subcls_plot(cell, DEG_total, CD8Tem)


## Classical monocytes
#DEG_CD14Mono_pos <- subset(DEG_CD14Mono, avg_log2FC > 0.25 & p_val_adj <0.05 )
#DEG_CD14Mono_neg <- subset(DEG_CD14Mono, avg_log2FC < -0.25 & p_val_adj <0.05 )
#CD14Mono$GROUP <- factor(CD14Mono$GROUP, levels=c("Synovial fluid", 
#                                                  "A_acute",
#                                                  "Health"))
# p1_dim <- DimPlot(CD14Mono, split.by = "GROUP") +
#   theme(plot.title = element_text(size = 5),
#         axis.title.x = element_text(size=10),
#         axis.title.y = element_text(size=10)
#   )

# DF <- CD14Mono[[]] %>%
#   group_by(CASE, GROUP, seurat_clusters) %>%
#   summarise(cell_number = n()) %>%
#   mutate(total_number = sum(cell_number)) %>%
#   mutate(percent = cell_number/total_number*100)
# 
# p2_rel_bar <- ggplot(DF, 
#                      aes(x = CASE, y = percent, fill = seurat_clusters))+
#   geom_bar(stat = "identity")+
#   facet_grid(~GROUP, scales="free", space="free_x")+
#   theme(axis.text.x=element_text(angle=90, hjust=1),
#         legend.justification = "left")+
#   labs(fill='subclusters', x="case")
# 
# t10 <- get_all_marker_topN(DEG_CD14Mono_pos, 10)
# # DefaultAssay(CD14Mono) <- "RNA"
# exp.avg <-  AverageExpression(CD14Mono,
#                               assays = "RNA", 
#                               return.seurat = TRUE,
#                               slot = "data")
# mat <- GetAssayData(exp.avg, slot = "scale.data")
# mat <- as.matrix(mat[unique(t10$gene), ])
# p3_heatmap <- pheatmap(t(mat),
#                        cellheight = 12,
#                        cellwidth = 12,
#                        cluster_cols = T,
#                        cluster_rows = T,
#                        show_colnames = T,
#                        fontsize = 8,
#                        angle_col = 45)
# g3 = as.ggplot(p3_heatmap)
# plot_grid(p1_dim, p2_rel_bar, g3, ncol=1, 
#           labels=c("A", "B", "C"), 
#           align = "v",
#           axis="lr")
# ggsave("subcluster_Classical monocytes.pdf", width = 12, height = 10)
# ggsave("subcluster_Classical monocytes.png", width = 12, height = 10)

# SAVE ----
save(
  CD14Mono_funList,
  mDC_funList,
  pDC_funList,
  InterMono_funList,
  vd2gdT_funList,
  Th1Th17_funList,
  Treg_funList,
  MAIT_funList,
  CD8Tem_funList,
  CD14Mono_subList,
  mDC_subList,
  pDC_subList,
  InterMono_subList,
  vd2gdT_subList,
  Th1Th17_subList,
  Treg_subList,
  MAIT_subList,
  CD8Tem_subList,
  file = "DEG_plot.rda"
)

