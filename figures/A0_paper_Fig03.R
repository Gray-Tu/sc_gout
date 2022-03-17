# subclusters abundance diff.
setwd("With_SYN")

library(Seurat)
library(dplyr)
library(ggpubr)
library(cowplot)

# get subclustered objects
load("MainDEGusage_ScoreAdd_20220211.rda")
# 
# CD14Mono,
# InterMono,
# CD16Mono,
# pDC,
# mDC,
# CD8Tem,
# CD8Te,
# CD4Te,
# MAIT,
# NaiveCD4T,
# Treg,
# Th1_Th17,
# Th1,
# Th2,
# Th17,
# vd2gdT,
# nonvd2gdT

# Plasmablasts
# Macrophages
# NaiveCD8T
# CD8Tcm
# NK

# make DF 

case_to_group <- as.data.frame(list(
  "CASE" = c("MH", "RG", "N01F", "N05S", "YB", 
             "Case01", "Case02", "Case03", "Case04",
             "Case05", "Case06", "Case07",
             "Gout146", "Gout151", "Gout183"),
  "GROUP" = c("Health", "Health", "Health", "Health", "Health",
              "A_acute", "A_acute", "A_acute", "A_acute", 
              "A_acute", "A_acute", "A_acute",
              "Synovial fluid", "Synovial fluid", "Synovial fluid")
))

get_abun_DF <- function(object){
  df <- object[[]]
  
  DF <- df %>%
    group_by(CASE, GROUP, cell_type, seurat_clusters) %>%
    summarise(cell_number = n()) %>%
    mutate(total_number = sum(cell_number)) %>%
  mutate(percent = cell_number/total_number*100)
  
  # add zero if no cell
  total_cls <- unique(DF$seurat_clusters)
  total_case <- c("MH", "RG", "N01F", "N05S", "YB", 
  "Case04", "Case05", "Case06", "Case01",
  "Case02", "Case03", "Case07", 
  "Gout146", "Gout151", "Gout183")
  addzero_df <- data.frame()
  cell_type <- unique(DF$cell_type)
  for (C in total_case){
    
    case_df <- subset(DF, CASE==C)
    G <- subset(case_to_group, CASE==C)$GROUP
    for (cls in total_cls){
      print(paste(G, C, cls))
      if (length(row.names(subset(case_df, seurat_clusters==cls))) == 0){
        cell_number = 0
        total_number = 0
        percent = 0.0
        GROUP = G
        CASE = C
        seurat_clusters=cls
        temp_df <- data.frame(CASE, GROUP, cell_type, seurat_clusters, cell_number,total_number,percent)
        addzero_df <- rbind(addzero_df, temp_df)
      }
    }
  }
  outDF <- rbind(DF, addzero_df)
  outDF$GROUP <- factor(outDF$GROUP,
                     levels=c("Synovial fluid", "A_acute", "Health"))
  outDF$CASE <- factor(outDF$CASE,
                    levels=c("MH", "RG", "N01F", "N05S", "YB", 
                             "Case04", "Case05", "Case06", "Case01",
                             "Case02", "Case03", "Case07",
                             "Gout146", "Gout151", "Gout183"))
  outDF$seurat_clusters.cell <- paste(outDF$cell_type, outDF$seurat_clusters, sep="_")
  
  return(outDF)
}

CD14Mono.DF <- get_abun_DF(CD14Mono)
InterMono.DF <- get_abun_DF(InterMono)
CD16Mono.DF <- get_abun_DF(CD16Mono)
pDC.DF <- get_abun_DF(pDC)
mDC.DF <- get_abun_DF(mDC)
CD8Tem.DF <- get_abun_DF(CD8Tem)
CD8Te.DF <- get_abun_DF(CD8Te)
CD4Te.DF <- get_abun_DF(CD4Te)
MAIT.DF <- get_abun_DF(MAIT)
NaiveCD4T.DF <- get_abun_DF(NaiveCD4T)
Treg.DF <- get_abun_DF(Treg)
Th1_Th17.DF <- get_abun_DF(Th1_Th17)
Th1.DF <- get_abun_DF(Th1)
Th2.DF <- get_abun_DF(Th2)
Th17.DF <- get_abun_DF(Th17)
vd2gdT.DF <- get_abun_DF(vd2gdT)
nonvd2gdT.DF <- get_abun_DF(nonvd2gdT)
Plasmablasts.DF <- get_abun_DF(Plasmablasts)
Macrophages.DF <- get_abun_DF(Macrophages)
NaiveCD8T.DF <- get_abun_DF(NaiveCD8T)
CD8Tcm.DF <- get_abun_DF(CD8Tcm)
NK.DF <- get_abun_DF(NK)

total.DF <- rbind(
  CD14Mono.DF,
  InterMono.DF,
  CD16Mono.DF,
  pDC.DF,
  mDC.DF,
  CD8Tem.DF,
  CD8Te.DF,
  CD4Te.DF,
  MAIT.DF, 
  NaiveCD4T.DF, 
  Treg.DF, 
  Th1_Th17.DF,
  Th1.DF, 
  Th2.DF, 
  Th17.DF, 
  vd2gdT.DF, 
  nonvd2gdT.DF, 
  Plasmablasts.DF, 
  Macrophages.DF, 
  NaiveCD8T.DF, 
  CD8Tcm.DF, 
  NK.DF
)

write.csv(total.DF, file="A0_addZero_fig3_Total_AbundanceTable.csv")





# make Wlicoxn test
make_wlicox <- function(DF, cell){
  wlicox_out <- compare_means(percent~GROUP, 
                              DF, 
                              method = "wilcox.test", 
                              paired = FALSE,
                              group.by = "seurat_clusters")
  
  wlicox_out$cell <- cell
  
  return(wlicox_out)
}
CD14Mono.wlicox <- make_wlicox(CD14Mono.DF, paste( unique(CD14Mono.DF$cell_type)))
InterMono.wlicox <- make_wlicox(InterMono.DF, paste( unique(InterMono.DF$cell_type)))
CD16Mono.wlicox <- make_wlicox(CD16Mono.DF, paste( unique(CD16Mono.DF$cell_type)))
pDC.wlicox <- make_wlicox(pDC.DF, paste( unique(pDC.DF$cell_type)))
mDC.wlicox <- make_wlicox(mDC.DF, paste( unique(mDC.DF$cell_type)))
CD8Tem.wlicox <- make_wlicox(CD8Tem.DF, paste( unique(CD8Tem.DF$cell_type)))
CD8Te.wlicox <- make_wlicox(CD8Te.DF, paste( unique(CD8Te.DF$cell_type)))
CD4Te.wlicox <- make_wlicox(CD4Te.DF, paste( unique(CD4Te.DF$cell_type)))
MAIT.wlicox <- make_wlicox(MAIT.DF, paste( unique(MAIT.DF$cell_type)))
NaiveCD4T.wlicox <- make_wlicox(NaiveCD4T.DF, paste( unique(NaiveCD4T.DF$cell_type)))
Treg.wlicox <- make_wlicox(Treg.DF, paste( unique(Treg.DF$cell_type)))
Th1_Th17.wlicox <- make_wlicox(Th1_Th17.DF, paste( unique(Th1_Th17.DF$cell_type)))
Th1.wlicox <- make_wlicox(Th1.DF, paste( unique(Th1.DF$cell_type)))
Th2.wlicox <- make_wlicox(Th2.DF, paste( unique(Th2.DF$cell_type)))
Th17.wlicox <- make_wlicox(Th17.DF, paste( unique(Th17.DF$cell_type)))
vd2gdT.wlicox <- make_wlicox(vd2gdT.DF, paste( unique(vd2gdT.DF$cell_type)))
nonvd2gdT.wlicox <- make_wlicox(nonvd2gdT.DF, paste( unique(nonvd2gdT.DF$cell_type)))
Plasmablasts.wlicox <- make_wlicox(Plasmablasts.DF, paste( unique(Plasmablasts.DF$cell_type)))
Macrophages.wlicox <- make_wlicox(Macrophages.DF, paste( unique(Macrophages.DF$cell_type)))
NaiveCD8T.wlicox <- make_wlicox(NaiveCD8T.DF, paste( unique(NaiveCD8T.DF$cell_type)))
CD8Tcm.wlicox <- make_wlicox(CD8Tcm.DF, paste( unique(CD8Tcm.DF$cell_type)))
NK.wlicox <- make_wlicox(NK.DF, paste( unique(NK.DF$cell_type)))


total.wlicox <- rbind(
  CD14Mono.wlicox,
  InterMono.wlicox,
  pDC.wlicox,
  NaiveCD4T.wlicox,
  Plasmablasts.wlicox,
  CD16Mono.wlicox,
  Macrophages.wlicox,
  Th1_Th17.wlicox,
  NaiveCD8T.wlicox,
  CD8Tcm.wlicox,
  CD8Tem.wlicox,
  MAIT.wlicox,
  NK.wlicox,
  mDC.wlicox,
  Treg.wlicox,
  Th1.wlicox,
  Th2.wlicox,
  Th17.wlicox,
  vd2gdT.wlicox,
  nonvd2gdT.wlicox,
  CD8Te.wlicox,
  CD4Te.wlicox
)

select_cell.wlicox <- subset(total.wlicox, p < 0.05)
unique(select_cell.wlicox$cell)
write.csv(select_cell.wlicox, file="A0_addZero_fig3_wilcoxTable.csv")
write.csv(total.wlicox, file="A0_addZero_fig3_Total_wilcoxTable.csv")

# boxplot
make_boxplot <- function(DF, ylabstr, compar_list=list(c("Synovial fluid", "Health"),
                                                    c("A_acute", "Synovial fluid"),
                                                    c("A_acute", "Health"))){
  
  stat.test <- compare_means(percent~GROUP, 
                             DF, 
                             method = "wilcox.test",
                             paired = FALSE,
                             group.by = "seurat_clusters")
  #print(stat.test)
  cls_y <- DF %>% 
    group_by(seurat_clusters) %>%
    mutate(ymax = max(percent)) %>% 
    select(seurat_clusters, ymax) %>%
    unique()
  
  y_list <- list()
  for (c in seq(1:length(rownames(stat.test)))){
    cls = as.character(stat.test[c,"seurat_clusters"])
    t2 <- cls_y %>% subset(seurat_clusters == cls)
    y_list[[c]] = as.numeric(t2["ymax"] * (1 + 0.1 *(c%%3 +1)))
  }
  stat.test$y.position <- as.numeric(y_list)
  #print(stat.test)
  stat.test <- stat.test %>% subset(p < 0.05)
  p <- ggboxplot(DF,           
                 x="GROUP", 
                 y="percent",
                 fill="seurat_clusters")+
    labs(y=ylabstr, fill="subclusters")+
    theme_classic()+
    facet_wrap(~seurat_clusters, scales = "free")+
    geom_point(aes(colour = CASE))
  
  boxplot_out <- p + stat_pvalue_manual(stat.test, 
                                        label = "p.signif", 
                                        tip.length=0, vjust=0.4)
  # boxplot_out <- p +  stat_compare_means(comparisons = compar_list,
  #                                             label = "p.signif",
  #                                             tip.length=0,
  #                                             vjust=0.4,
  #                                             method = "wilcox.test") 
  return(boxplot_out)
}

boxplot.CD14Mono <- make_boxplot(CD14Mono.DF, paste('% of', unique(CD14Mono.DF$cell_type)))
boxplot.InterMono <- make_boxplot(InterMono.DF, paste('% of', unique(InterMono.DF$cell_type)))
boxplot.CD16Mono <- make_boxplot(CD16Mono.DF, paste('% of', unique(CD16Mono.DF$cell_type)))
boxplot.pDC <- make_boxplot(pDC.DF, paste('% of', unique(pDC.DF$cell_type)))
boxplot.mDC <- make_boxplot(mDC.DF, paste('% of', unique(mDC.DF$cell_type)))
boxplot.CD8Tem <- make_boxplot(CD8Tem.DF, paste('% of', unique(CD8Tem.DF$cell_type)))
boxplot.CD8Te <- make_boxplot(CD8Te.DF, paste('% of', unique(CD8Te.DF$cell_type)))
boxplot.CD4Te <- make_boxplot(CD4Te.DF, paste('% of', unique(CD4Te.DF$cell_type)))
boxplot.MAIT <- make_boxplot(MAIT.DF, paste('% of', unique(MAIT.DF$cell_type)))
boxplot.NaiveCD4T <- make_boxplot(NaiveCD4T.DF, paste('% of', unique(NaiveCD4T.DF$cell_type)))
boxplot.Treg <- make_boxplot(Treg.DF, paste('% of', unique(Treg.DF$cell_type)))
boxplot.Th1_Th17 <- make_boxplot(Th1_Th17.DF, paste('% of', unique(Th1_Th17.DF$cell_type)))
boxplot.Th1 <- make_boxplot(Th1.DF, paste('% of', unique(Th1.DF$cell_type)))
boxplot.Th2 <- make_boxplot(Th2.DF, paste('% of', unique(Th2.DF$cell_type)))
boxplot.Th17 <- make_boxplot(Th17.DF, paste('% of', unique(Th17.DF$cell_type)))
boxplot.vd2gdT <- make_boxplot(vd2gdT.DF, paste('% of', unique(vd2gdT.DF$cell_type)))
boxplot.nonvd2gdT <- make_boxplot(nonvd2gdT.DF, paste('% of', unique(nonvd2gdT.DF$cell_type)))
boxplot.Plasmablasts <- make_boxplot(Plasmablasts.DF, paste('% of', unique(Plasmablasts.DF$cell_type)))
boxplot.Macrophages <- make_boxplot(Macrophages.DF, paste('% of', unique(Macrophages.DF$cell_type)))
boxplot.NaiveCD8T <- make_boxplot(NaiveCD8T.DF, paste('% of', unique(NaiveCD8T.DF$cell_type)))
boxplot.CD8Tcm <- make_boxplot(CD8Tcm.DF, paste('% of', unique(CD8Tcm.DF$cell_type)))
boxplot.NK <- make_boxplot(NK.DF, paste('% of', unique(NK.DF$cell_type)))


library(clusterProfiler)
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
DEG_Plasmablasts <- read.csv('RNADEG_Plasmablasts.csv') %>%
  mutate(cell_type="Plasmablasts")
DEG_Macrophages <- read.csv('RNADEG_Macrophages.csv') %>%
  mutate(cell_type="Macrophages")
DEG_NaiveCD8T <- read.csv('RNADEG_NaiveCD8T.csv') %>%
  mutate(cell_type="Naive CD8 T cells")
DEG_CD8Tcm <- read.csv('RNADEG_CD8Tcm.csv') %>%
  mutate(cell_type="Central memory CD8 T cells")
DEG_NK <- read.csv('RNADEG_NK.csv') %>%
  mutate(cell_type="Natural killer cells")
  

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
                   DEG_nonvd2gdT,
                   DEG_Plasmablasts,
                   DEG_Macrophages,
                   DEG_NaiveCD8T,
                   DEG_CD8Tcm,
                   DEG_NK
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

write.csv(DEG_total, file="ADV_03_subcluster_DEG.csv")

make_subcls_dimbardot_plot <- function(cell, DEG_df, seurat_obj, DF){
  pos_DEG <- subset(DEG_df, cell_type==cell & UPDN=="UP")
  cell <- str_replace(cell, "/", ".")
  seurat_obj$GROUP <- factor(seurat_obj$GROUP, levels=c("Synovial fluid", 
                                                        "A_acute",
                                                        "Health"))
  DefaultAssay(seurat_obj) <- "RNA"
  t10 <- get_all_marker_topN(pos_DEG, 10)
  write.csv(t10, paste("ADV_Fig03_top10DEG", cell, ".csv"))
  select_genes <- unique(t10$gene)
  p_genedot <- DotPlot(seurat_obj, 
                       features=rev(select_genes))+
    theme(axis.text.y = element_text(size=8))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))+
    labs(y="subclusters")
  
  p_dimplot <- DimPlot(seurat_obj, split.by = "GROUP", ncol=1) +
    theme(plot.title = element_text(size = 5),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10))+NoLegend()
  p_subbar <- ggplot(DF, 
                     aes(y = CASE, x = percent, fill = seurat_clusters))+
    geom_bar(stat = "identity")+
    facet_grid(rows="GROUP", scales="free", space="free")+
    theme(axis.text.x=element_text(angle=90, hjust=1),
          legend.justification = "left")+
    theme_classic2()+
    labs(fill='subclusters', y="", x="frequency")
  
  returnlist <- list("p_genedot"=p_genedot, "p_subbar"=p_subbar, "p_dimplot"=p_dimplot)
  return(returnlist)
}
get_all_marker_topN <- function(markers, topN){
  temp <-  markers %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC)
  return(temp)
}

make_advFig03 <- function(cell, DEG_df, seurat_obj, DF, boxplot_p){
  plist <- make_subcls_dimbardot_plot(cell, 
                                      DEG_df,
                                      seurat_obj,
                                      DF)
  d <- plist$p_genedot
  b <- plist$p_subbar
  a <- plist$p_dimplot
  c <- boxplot_p+facet_wrap(~seurat_clusters, scales = "free", ncol=2)+scale_x_discrete(guide = guide_axis(n.dodge=2))
  pabc <- plot_grid(
    a,
    b,
    c,
    ncol=3, 
    labels=c("a", "b", "c"),
    rel_widths = c(1, 0.8, 1.2)
  )+theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  pabcd <- plot_grid(
    pabc,
    d,
    ncol=1,
    labels=c("", "d"),
    rel_heights = c(1.8, 1)
  )+theme(plot.background = element_rect(fill = "white", colour = "white"))
  cell = str_replace(cell, "/", ".")
  ggsave(paste("ADV_Fig03_", cell, ".png", sep=""),
         plot=pabcd, width = 12, height = 12)
  ggsave(paste("ADV_Fig03_", cell, ".pdf", sep=""),
         plot=pabcd, width = 12, height = 12)
  ggsave(paste("ADV_Fig03_", cell, ".jpg", sep=""),
         plot=pabcd, width = 12, height = 12, dpi=600)
}

make_advFig03("Classical monocytes", 
              DEG_total,
              CD14Mono,
              CD14Mono.DF,
              boxplot.CD14Mono)
make_advFig03("Myeloid dendritic cells", 
              DEG_total,
              mDC,
              mDC.DF,
              boxplot.mDC)
make_advFig03("Plasmacytoid dendritic cells", 
              DEG_total,
              pDC,
              pDC.DF,
              boxplot.pDC)
make_advFig03("Intermediate monocytes", 
              DEG_total,
              InterMono,
              InterMono.DF,
              boxplot.InterMono)
make_advFig03("Vd2 gd T cells", 
              DEG_total,
              vd2gdT,
              vd2gdT.DF,
              boxplot.vd2gdT)
make_advFig03("Th1/Th17 cells", 
              DEG_total,
              Th1_Th17,
              Th1_Th17.DF,
              boxplot.Th1_Th17)
make_advFig03("T regulatory cells", 
              DEG_total,
              Treg,
              Treg.DF,
              boxplot.Treg)
make_advFig03("MAIT cells", 
              DEG_total,
              MAIT,
              MAIT.DF,
              boxplot.MAIT)
make_advFig03("Effector memory CD8 T cells", 
              DEG_total,
              CD8Tem,
              CD8Tem.DF,
              boxplot.CD8Tem)
make_advFig03("Naive CD4 T cells", 
              DEG_total,
              NaiveCD4T,
              NaiveCD4T.DF,
              boxplot.NaiveCD4T)
make_advFig03("Non classical monocytes", 
              DEG_total,
              CD16Mono,
              CD16Mono.DF,
              boxplot.CD16Mono)
make_advFig03("Terminal effector CD8 T cells", 
              DEG_total,
              CD8Te,
              CD8Te.DF,
              boxplot.CD8Te)
make_advFig03("Terminal effector CD4 T cells", 
              DEG_total,
              CD4Te,
              CD4Te.DF,
              boxplot.CD4Te)
make_advFig03("Terminal effector CD4 T cells", 
              DEG_total,
              CD4Te,
              CD4Te.DF,
              boxplot.CD4Te)
make_advFig03("Th1 cells", 
              DEG_total,
              Th1,
              Th1.DF,
              boxplot.Th1)
make_advFig03("Th2 cells", 
              DEG_total,
              Th2,
              Th2.DF,
              boxplot.Th2)
make_advFig03("Th17 cells", 
              DEG_total,
              Th17,
              Th17.DF,
              boxplot.Th17)
make_advFig03("Non-Vd2 gd T cells", 
              DEG_total,
              nonvd2gdT,
              nonvd2gdT.DF,
              boxplot.nonvd2gdT)

make_advFig03("Plasmablasts", 
              DEG_total,
              Plasmablasts,
              Plasmablasts.DF,
              boxplot.Plasmablasts)
make_advFig03("Macrophages", 
              DEG_total,
              Macrophages,
              Macrophages.DF,
              boxplot.Macrophages)
make_advFig03("Naive CD8 T cells", 
              DEG_total,
              NaiveCD8T,
              NaiveCD8T.DF,
              boxplot.NaiveCD8T)
make_advFig03("Central memory CD8 T cells", 
              DEG_total,
              CD8Tcm,
              CD8Tcm.DF,
              boxplot.CD8Tcm)
# make_advFig03("Natural killer cells", 
#               DEG_total,
#               NK,
#               NK.DF,
#               boxplot.NK)


# adjust for NK part
make_subcls_dimbardot_plot_NK <- function(cell, DEG_df, seurat_obj, DF){
  pos_DEG <- subset(DEG_df, cell_type==cell & UPDN=="UP")
  cell <- str_replace(cell, "/", ".")
  seurat_obj$GROUP <- factor(seurat_obj$GROUP, levels=c("Synovial fluid", 
                                                        "A_acute",
                                                        "Health"))
  DefaultAssay(seurat_obj) <- "RNA"
  t10 <- get_all_marker_topN(pos_DEG, 10)
  select_genes <- unique(t10$gene)
  p_genedot <- DotPlot(seurat_obj, 
                       features=rev(select_genes))+
    theme(axis.text.y = element_text(size=8))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8))+
    labs(y="subclusters")
  
  p_dimplot <- DimPlot(seurat_obj, split.by = "GROUP", ncol=1) +
    theme(plot.title = element_text(size = 5),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10))+NoLegend()
  p_subbar <- ggplot(DF, 
                     aes(y = CASE, x = percent, fill = seurat_clusters))+
    geom_bar(stat = "identity")+
    facet_grid(rows="GROUP", scales="free", space="free")+
    theme(axis.text.x=element_text(angle=90, hjust=1),
          legend.justification = "left")+
    theme_classic2()+
    labs(fill='subclusters', y="", x="frequency")
  
  returnlist <- list("p_genedot"=p_genedot, "p_subbar"=p_subbar, "p_dimplot"=p_dimplot)
  return(returnlist)
}
make_advFig03_NK <- function(cell, DEG_df, seurat_obj, DF, boxplot_p){
  plist <- make_subcls_dimbardot_plot_NK(cell, 
                                      DEG_df,
                                      seurat_obj,
                                      DF)
  d <- plist$p_genedot
  b <- plist$p_subbar
  a <- plist$p_dimplot
  c <- boxplot_p+facet_wrap(~seurat_clusters, scales = "free_y", ncol=3)+
    scale_x_discrete(guide = guide_axis(n.dodge=3))
  pabc <- plot_grid(
    a,
    b,
    c,
    ncol=3, 
    labels=c("a", "b", "c"),
    rel_widths = c(1, 0.8, 1.2)
  )+theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  pabcd <- plot_grid(
    pabc,
    d,
    ncol=1,
    labels=c("", "d"),
    rel_heights = c(1.8, 1)
  )+theme(plot.background = element_rect(fill = "white", colour = "white"))
  cell = str_replace(cell, "/", ".")
  ggsave(paste("ADV_Fig03_", cell, ".png", sep=""),
         plot=pabcd, width = 12, height = 12)
  ggsave(paste("ADV_Fig03_", cell, ".pdf", sep=""),
         plot=pabcd, width = 12, height = 12)
  ggsave(paste("ADV_Fig03_", cell, ".jpg", sep=""),
         plot=pabcd, width = 12, height = 12, dpi=600)
}
make_advFig03_NK("Natural killer cells", 
                 DEG_total,
                 NK,
                 NK.DF,
                 boxplot.NK)


