# use singleR results
library(SingleR)
library(celldex)
library(Seurat)
library(ggplot2)
setwd("With_SYN")
# PBMC----
load("../plan_1124/pc50m2h5_InfScore_rpca_Gout_combinePBMC.Rda")
# Gout.combined.50

# add bpe.fine.calls if M1 M2
bpe <- BlueprintEncodeData()
singleR_warp_fine <- function(data, ref){
  pred.fine <- SingleR(test=GetAssayData(data, slot = "count"),
                       assay.type.test=1,
                       ref=ref,
                       labels=ref$label.fine)
  return(pred.fine$labels)
}
Gout.combined.50$bpe.fine.calls <- singleR_warp_fine(Gout.combined.50, bpe)

# 
Macro_Mega <- c("lung Macro", "macrophage", "INF-Macro", "Megakaryocte")
my_DF <- Gout.combined.50[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(singleR.clShort %in% Macro_Mega, singleR.clShort, monaco.fine.calls)) 

Gout.combined.50$cell_type <- my_DF$cell_type
#
bpe_marco <- c("Macrophages M1", "Macrophages M2")
my_DF <- Gout.combined.50[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(bpe.fine.calls %in% bpe_marco,
                            bpe.fine.calls, 
                            cell_type)) 

Gout.combined.50$cell_type <- my_DF$cell_type
# Upper macrophage
my_DF <- Gout.combined.50[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="macrophage",
                            "Macrophages", 
                            cell_type))
Gout.combined.50$cell_type <- my_DF$cell_type
# add -s
my_DF <- Gout.combined.50[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="Megakaryocte",
                            "Megakaryoctes", 
                            cell_type))
Gout.combined.50$cell_type <- my_DF$cell_type

Gout.combined.50$cell_type <- factor(Gout.combined.50$cell_type, 
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


Acute_HC_set <- subset(Gout.combined.50, GROUP %in% c("A_acute", "Health"))
save(Gout.combined.50, file="addbpe_pc50m2h5_InfScore_rpca_Gout_combinePBMC.Rda")

p1 <- DimPlot(Acute_HC_set, 
        group.by = "cell_type", 
        cols=DiscretePalette(32, "polychrome"),
        split.by = "GROUP"
        )+NoLegend()

p2 <- DimPlot(Acute_HC_set, 
        group.by = "cell_type", 
        cols=DiscretePalette(32, "polychrome"),
        label.box=T,
        label.color = "white",
        label=T,
        repel=T,
        )

ggsave("AcuteHCPBMC_annotation.pdf", 
       plot=plot_grid(p1, p2, ncol = 1),
       width = 14, height = 14)

# SYN----
load("../SYN_data/pc50m2_autolabel_rpca_syn_combined.Rda")
# syn.combined
#
my_DF <- syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(singleR.clShort %in% Macro, singleR.clShort, monaco.fine.calls)) 

syn.combined$cell_type <- my_DF$cell_type
#
my_DF <- syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(bpe.fine.calls %in% bpe_marco,
                            bpe.fine.calls, 
                            cell_type)) 

syn.combined$cell_type <- my_DF$cell_type
# Upper macrophage
my_DF <- syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="macrophage",
                            "Macrophages", 
                            cell_type))
syn.combined$cell_type <- my_DF$cell_type
# add -s
my_DF <- syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="Megakaryocte",
                            "Megakaryoctes", 
                            cell_type))
syn.combined$cell_type <- my_DF$cell_type

syn.combined$cell_type <- factor(syn.combined$cell_type, 
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

p3 <- DimPlot(syn.combined, 
              group.by = "cell_type", 
              cols=DiscretePalette(32, "polychrome"),
              split.by = "CASE"
)+NoLegend()

p4 <- DimPlot(syn.combined, 
        group.by = "cell_type", 
        cols=DiscretePalette(32, "polychrome"),
        label.box=T,
        label.color = "white",
        label=T,
        repel=T,
)#+theme(legend.position = "bottom")
ggsave("SYN_annotation.pdf", 
       plot=plot_grid(p3, p4, ncol = 1, rel_heights=c(1, 1.2)),
       width = 14, height = 7)
save(syn.combined, file="pc50_syn.rda")
