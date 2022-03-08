setwd("With_SYN")

library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplotify)
library(pheatmap)
# load data ----
load(file = "MainDEGusage_ScoreAdd.rda") # from 4-1_EnergyScore

CD14Mono
InterMono
mDC
pDC
NaiveCD4T
Treg
Th1_Th17
CD8Tem
vd2gdT
MAIT

load(file="integrated_rpca_syn_PBMCAcuteHealth.Rda")
# pbmc.syn.combined

# make subclusters of ----
# 1. CD16Mono
# 2. CD8Te
# 3. CD4Te
# 4. Th1
# 5. Th2
# 6. Th17
# 7. B (?)

# 1. CD16Mono
CELL = "Non classical monocytes"
CD16Mono <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(CD16Mono, split.by="GROUP")
CD16Mono <- RunPCA(CD16Mono, npcs = 50, verbose = FALSE)
CD16Mono <- RunUMAP(CD16Mono, reduction = "pca", dims = 1:50)
CD16Mono <- FindNeighbors(CD16Mono, reduction = "pca", dims = 1:50)
CD16Mono <- FindClusters(CD16Mono, resolution = 0.4)
DimPlot(CD16Mono, split.by="GROUP")
DimPlot(CD16Mono, group.by="GROUP")

# 2. CD8Te
CELL = "Terminal effector CD8 T cells"
CD8Te <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(CD8Te, split.by="GROUP")
CD8Te <- RunPCA(CD8Te, npcs = 50, verbose = FALSE)
CD8Te <- RunUMAP(CD8Te, reduction = "pca", dims = 1:50)
CD8Te <- FindNeighbors(CD8Te, reduction = "pca", dims = 1:50)
CD8Te <- FindClusters(CD8Te, resolution = 0.5)
DimPlot(CD8Te, split.by="GROUP")
DimPlot(CD8Te, group.by="GROUP")

# 3. CD4Te
CELL = "Terminal effector CD4 T cells"
CD4Te <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(CD4Te, split.by="GROUP")
CD4Te <- RunPCA(CD4Te, npcs = 50, verbose = FALSE)
CD4Te <- RunUMAP(CD4Te, reduction = "pca", dims = 1:50)
CD4Te <- FindNeighbors(CD4Te, reduction = "pca", dims = 1:50)
CD4Te <- FindClusters(CD4Te, resolution = 0.5)
DimPlot(CD4Te, split.by="GROUP")
DimPlot(CD4Te, group.by="GROUP")

# 4. Th1
CELL = "Th1 cells"
Th1 <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Th1, split.by="GROUP")
Th1 <- RunPCA(Th1, npcs = 50, verbose = FALSE)
Th1 <- RunUMAP(Th1, reduction = "pca", dims = 1:50)
Th1 <- FindNeighbors(Th1, reduction = "pca", dims = 1:50)
Th1 <- FindClusters(Th1, resolution = 0.4)
DimPlot(Th1, split.by="GROUP")
DimPlot(Th1, group.by="GROUP")

# 5. Th2
CELL = "Th2 cells"
Th2 <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Th2, split.by="GROUP")
Th2 <- RunPCA(Th2, npcs = 50, verbose = FALSE)
Th2 <- RunUMAP(Th2, reduction = "pca", dims = 1:50)
Th2 <- FindNeighbors(Th2, reduction = "pca", dims = 1:50)
Th2 <- FindClusters(Th2, resolution = 0.4)
DimPlot(Th2, split.by="GROUP")
DimPlot(Th2, group.by="GROUP")

# 6. Th17
CELL = "Th17 cells"
Th17 <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Th17, split.by="GROUP")
Th17 <- RunPCA(Th17, npcs = 50, verbose = FALSE)
Th17 <- RunUMAP(Th17, reduction = "pca", dims = 1:50)
Th17 <- FindNeighbors(Th17, reduction = "pca", dims = 1:50)
Th17 <- FindClusters(Th17, resolution = 0.6)
DimPlot(Th17, split.by="GROUP")
DimPlot(Th17, group.by="GROUP")

# 7. 
CELL = "Non-Vd2 gd T cells"
nonvd2gdT <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(nonvd2gdT, split.by="GROUP")
nonvd2gdT <- RunPCA(nonvd2gdT, npcs = 50, verbose = FALSE)
nonvd2gdT <- RunUMAP(nonvd2gdT, reduction = "pca", dims = 1:50)
nonvd2gdT <- FindNeighbors(nonvd2gdT, reduction = "pca", dims = 1:50)
nonvd2gdT <- FindClusters(nonvd2gdT, resolution = 0.3)
DimPlot(nonvd2gdT, split.by="GROUP")
DimPlot(nonvd2gdT, group.by="GROUP")

# make DEG
DefaultAssay(CD16Mono) <- "RNA"
RNADEG_CD16Mono <- FindAllMarkers(CD16Mono)
write.csv(RNADEG_CD16Mono, file="RNADEG_CD16Mono.csv")

DefaultAssay(CD8Te) <- "RNA"
RNADEG_CD8Te <- FindAllMarkers(CD8Te)
write.csv(RNADEG_CD8Te, file="RNADEG_CD8Te.csv")

DefaultAssay(CD4Te) <- "RNA"
RNADEG_CD4Te <- FindAllMarkers(CD4Te)
write.csv(RNADEG_CD4Te, file="RNADEG_CD4Te.csv")

DefaultAssay(Th1) <- "RNA"
RNADEG_Th1 <- FindAllMarkers(Th1)
write.csv(RNADEG_Th1, file="RNADEG_Th1.csv")

DefaultAssay(Th2) <- "RNA"
RNADEG_Th2 <- FindAllMarkers(Th2)
write.csv(RNADEG_Th2, file="RNADEG_Th2.csv")

DefaultAssay(Th17) <- "RNA"
RNADEG_Th17 <- FindAllMarkers(Th17)
write.csv(RNADEG_Th17, file="RNADEG_Th17.csv")

DefaultAssay(nonvd2gdT) <- "RNA"
RNADEG_nonvd2gdT <- FindAllMarkers(nonvd2gdT)
write.csv(RNADEG_nonvd2gdT, file="RNADEG_nonvd2gdT.csv")

# add score----
add_INFscore <- function(seurat_obj, gene_df){
  DefaultAssay(seurat_obj) <- "RNA"
  obj_gene <- row.names(seurat_obj)
  in_gene <- gene_df[gene_df$GENE %in% obj_gene, ]
  in_gene_list <- list(in_gene)
  return_seurat_obj <- AddModuleScore(seurat_obj, 
                                      features=in_gene_list,
                                      name="inflammatory_score")
  return(return_seurat_obj)
}
# Read genes 
INFgene_df <- read.csv("HALLMARK_INFLAMMATORY_RESPONSE_geneset.txt")
CD16Mono <- add_INFscore(CD16Mono, INFgene_df)
CD8Te <- add_INFscore(CD8Te, INFgene_df)
CD4Te <- add_INFscore(CD4Te, INFgene_df)
Th1 <- add_INFscore(Th1, INFgene_df)
Th2 <- add_INFscore(Th2, INFgene_df)
Th17 <- add_INFscore(Th17, INFgene_df)
nonvd2gdT <- add_INFscore(nonvd2gdT, INFgene_df)

total_genes <- data.frame()
gene_file_list <- Sys.glob("Energy_score_genes/Gene_list/*txt", dirmark = FALSE)
for (f in  gene_file_list){
  df <- read.csv(f, skip=1)
  total_genes <- union(df[,1], total_genes)
}
#total_genes # 595 genes

add_ENGscore <- function(seurat_obj, gene_list){
  DefaultAssay(seurat_obj) <- "RNA"
  obj_gene <- row.names(seurat_obj)
  in_gene <- gene_list[gene_list %in% obj_gene]
  in_gene_list <- list(in_gene)
  return_seurat_obj <- AddModuleScore(seurat_obj, 
                                      features=in_gene_list,
                                      name="energy_score")
  return(return_seurat_obj)
}
CD16Mono <- add_ENGscore(CD16Mono, total_genes)
CD8Te <- add_ENGscore(CD8Te, total_genes)
CD4Te <- add_ENGscore(CD4Te, total_genes)
Th1 <- add_ENGscore(Th1, total_genes)
Th2 <- add_ENGscore(Th2, total_genes)
Th17 <- add_ENGscore(Th17, total_genes)
nonvd2gdT <- add_ENGscore(nonvd2gdT, total_genes)
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
     CD16Mono,
     CD8Te,
     CD4Te,
     Th1,
     Th2,
     Th17,
     nonvd2gdT,
     file = "MainDEGusage_ScoreAdd_20220204.rda"
)

# make Func----
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
DEG_NaiveCD4T   <- read.csv("RNADEG_NaiveCD4T.csv") %>% 
  mutate(cell_type="Naive CD4 T cells")

DEG_CD16Mono   <- read.csv("RNADEG_CD16Mono.csv") %>% 
  mutate(cell_type="Non classical monocytes")
DEG_CD8Te   <- read.csv("RNADEG_CD8Te.csv") %>% 
  mutate(cell_type="Terminal effector CD8 T cells")
DEG_CD4Te   <- read.csv("RNADEG_CD8Tem.csv") %>% 
  mutate(cell_type="Terminal effector CD4 T cells")
DEG_Th1   <- read.csv("RNADEG_Th1.csv") %>% 
  mutate(cell_type="Th1 cells")
DEG_Th2   <- read.csv("RNADEG_Th2.csv") %>% 
  mutate(cell_type="Th2 cells")
DEG_Th17 <- read.csv("RNADEG_Th17.csv") %>% 
  mutate(cell_type="Th17 cells")
DEG_nonvd2gdT <- read.csv("RNADEG_nonvd2gdT.csv") %>%
  mutate(cell_type="Non-Vd2 gd T cells")

DEG_total <- rbind(DEG_CD14Mono,
                   DEG_mDC,
                   DEG_pDC,
                   DEG_InterMono,
                   DEG_vd2gdT,
                   DEG_Th1_Th17,
                   DEG_Treg,
                   DEG_MAIT,
                   DEG_CD8Tem,
                   DEG_NaiveCD4T,
                   DEG_CD16Mono,
                   DEG_CD8Te,
                   DEG_CD4Te,
                   DEG_Th1,
                   DEG_Th2,
                   DEG_Th17,
                   DEG_nonvd2gdT
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

# [1] "Classical monocytes"           "Myeloid dendritic cells"      
# [3] "Plasmacytoid dendritic cells"  "Intermediate monocytes"       
# [5] "Vd2 gd T cells"                "Th1/Th17 cells"               
# [7] "T regulatory cells"            "MAIT cells"                   
# [9] "Effector memory CD8 T cells"   "Naive CD4 T cells"            
# [11] "Non classical monocytes"       "Terminal effector CD8 T cells"
# [13] "Terminal effector CD4 T cells" "Th1 cells"                    
# [15] "Th2 cells"                     "Th17 cells"  
# "Non-Vd2 gd T cells"

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
cell <- "Naive CD4 T cells"
NaiveCD4T_funList <- make_kegg_gobp_func(cell, DEG_total)

cell <- "Non classical monocytes"
CD16Mono_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Terminal effector CD8 T cells"
CD8Te_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Terminal effector CD4 T cells"
CD4Te_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Th1 cells"
Th1_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Th2 cells"
Th2_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Th17 cells"
Th17_funList <- make_kegg_gobp_func(cell, DEG_total)
cell <- "Non-Vd2 gd T cells"
nonvd2gdT_funList <- make_kegg_gobp_func(cell, DEG_total)



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
cell <- "Naive CD4 T cells"
NaiveCD4T_subList <- make_subcls_plot(cell, DEG_total, NaiveCD4T)

cell <- "Non classical monocytes"
CD16Mono_subList <- make_subcls_plot(cell, DEG_total, CD16Mono)
cell <- "Terminal effector CD8 T cells"
CD8Te_subList <- make_subcls_plot(cell, DEG_total, CD8Te)
cell <- "Terminal effector CD4 T cells"
CD4Te_subList <- make_subcls_plot(cell, DEG_total, CD4Te)
cell <- "Th1 cells"
Th1_subList <- make_subcls_plot(cell, DEG_total, Th1)
cell <- "Th2 cells"
Th2_subList <- make_subcls_plot(cell, DEG_total, Th2)
cell <- "Th17 cells"
Th17_subList <- make_subcls_plot(cell, DEG_total, Th17)
cell <- "Non-Vd2 gd T cells"
nonvd2gdT_subList <- make_subcls_plot(cell, DEG_total, nonvd2gdT)




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
  NaiveCD4T_funList,
  CD16Mono_funList,
  CD8Te_funList,
  CD4Te_funList,
  Th1_funList,
  Th2_funList,
  Th17_funList,
  nonvd2gdT_funList,
  
  CD14Mono_subList,
  mDC_subList,
  pDC_subList,
  InterMono_subList,
  vd2gdT_subList,
  Th1Th17_subList,
  Treg_subList,
  MAIT_subList,
  CD8Tem_subList,
  NaiveCD4T_subList,
  CD16Mono_subList,
  CD8Te_subList,
  CD4Te_subList,
  Th1_subList,
  Th2_subList,
  Th17_subList,
  nonvd2gdT_subList,
  
  file = "DEG_plot.rda"
)

# 20220211---- 
load("MainDEGusage_ScoreAdd_20220204.rda")
# add 
# Plasmablasts
# Macrophages
# Naive CD8 T cells
# Central memory CD8 T cells
# Natural killer cells

## subclustering---- 
CELL = "Plasmablasts"
Plasmablasts <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Plasmablasts, split.by="GROUP")
Plasmablasts <- RunPCA(Plasmablasts, npcs = 50, verbose = FALSE)
Plasmablasts <- RunUMAP(Plasmablasts, reduction = "pca", dims = 1:50)
Plasmablasts <- FindNeighbors(Plasmablasts, reduction = "pca", dims = 1:50)
Plasmablasts <- FindClusters(Plasmablasts, resolution = 0.5)
DimPlot(Plasmablasts, split.by="GROUP")
DimPlot(Plasmablasts, group.by="GROUP")
DimPlot(Plasmablasts)

CELL = "Macrophages"
Macrophages <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(Macrophages, split.by="GROUP")
Macrophages <- RunPCA(Macrophages, npcs = 50, verbose = FALSE)
Macrophages <- RunUMAP(Macrophages, reduction = "pca", dims = 1:50)
Macrophages <- FindNeighbors(Macrophages, reduction = "pca", dims = 1:50)
Macrophages <- FindClusters(Macrophages, resolution = 0.6)
DimPlot(Macrophages, split.by="GROUP")
DimPlot(Macrophages, group.by="GROUP")
DimPlot(Macrophages)

CELL = "Naive CD8 T cells"
NaiveCD8T <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(NaiveCD8T, split.by="GROUP")
NaiveCD8T <- RunPCA(NaiveCD8T, npcs = 50, verbose = FALSE)
NaiveCD8T <- RunUMAP(NaiveCD8T, reduction = "pca", dims = 1:50)
NaiveCD8T <- FindNeighbors(NaiveCD8T, reduction = "pca", dims = 1:50)
NaiveCD8T <- FindClusters(NaiveCD8T, resolution = 0.6)
DimPlot(NaiveCD8T, split.by="GROUP")
DimPlot(NaiveCD8T, group.by="GROUP")
DimPlot(NaiveCD8T)

CELL = "Central memory CD8 T cells"
CD8Tcm <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(CD8Tcm, split.by="GROUP")
CD8Tcm <- RunPCA(CD8Tcm, npcs = 50, verbose = FALSE)
CD8Tcm <- RunUMAP(CD8Tcm, reduction = "pca", dims = 1:50)
CD8Tcm <- FindNeighbors(CD8Tcm, reduction = "pca", dims = 1:50)
CD8Tcm <- FindClusters(CD8Tcm, resolution = 0.6)
DimPlot(CD8Tcm, split.by="GROUP")
DimPlot(CD8Tcm, group.by="GROUP")
DimPlot(CD8Tcm)

#Natural killer cells
CELL = "Natural killer cells"
NK <- subset(pbmc.syn.combined, cell_type==CELL)
DimPlot(NK, split.by="GROUP")
NK <- RunPCA(NK, npcs = 50, verbose = FALSE)
NK <- RunUMAP(NK, reduction = "pca", dims = 1:50)
NK <- FindNeighbors(NK, reduction = "pca", dims = 1:50)
NK <- FindClusters(NK, resolution = 0.6)
DimPlot(NK, split.by="GROUP")
DimPlot(NK, group.by="GROUP")
DimPlot(NK)

## add score----
Plasmablasts <- add_INFscore(Plasmablasts, INFgene_df)
Macrophages <- add_INFscore(Macrophages, INFgene_df)
NaiveCD8T <- add_INFscore(NaiveCD8T, INFgene_df)
CD8Tcm <- add_INFscore(CD8Tcm, INFgene_df)
NK <- add_INFscore(NK, INFgene_df)

Plasmablasts <- add_ENGscore(Plasmablasts, total_genes)
Macrophages <- add_ENGscore(Macrophages, total_genes)
NaiveCD8T <- add_ENGscore(NaiveCD8T, total_genes)
CD8Tcm <- add_ENGscore(CD8Tcm, total_genes)
NK <- add_ENGscore(NK, total_genes)

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
     CD16Mono,
     CD8Te,
     CD4Te,
     Th1,
     Th2,
     Th17,
     nonvd2gdT,
     Plasmablasts,
     Macrophages,
     NaiveCD8T,
     CD8Tcm,
     NK,
     file = "MainDEGusage_ScoreAdd_20220211.rda"
)

## FindMarkers
DefaultAssay(Plasmablasts) <- "RNA"
RNADEG_Plasmablasts <- FindAllMarkers(Plasmablasts)
write.csv(RNADEG_Plasmablasts, file="RNADEG_Plasmablasts.csv")

DefaultAssay(Macrophages) <- "RNA"
RNADEG_Macrophages <- FindAllMarkers(Macrophages)
write.csv(RNADEG_Macrophages, file="RNADEG_Macrophages.csv")

DefaultAssay(NaiveCD8T) <- "RNA"
RNADEG_NaiveCD8T <- FindAllMarkers(NaiveCD8T)
write.csv(RNADEG_NaiveCD8T, file="RNADEG_NaiveCD8T.csv")

DefaultAssay(CD8Tcm) <- "RNA"
RNADEG_CD8Tcm <- FindAllMarkers(CD8Tcm)
write.csv(RNADEG_CD8Tcm, file="RNADEG_CD8Tcm.csv")

DefaultAssay(NK) <- "RNA"
RNADEG_NK <- FindAllMarkers(NK)
write.csv(RNADEG_NK, file="RNADEG_NK.csv")
