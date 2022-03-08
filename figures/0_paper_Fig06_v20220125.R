# PBMC cell cell interaction

setwd("With_SYN")
library(Seurat)
library(CellChat)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

# https://github.com/sqjin/CellChat/issues/149
setIdent <- function(object, ident.use = NULL, levels = NULL){
  object@idents <- as.factor(object@meta[[ident.use]])
  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      object@net$prob <- object@net$prob[idx, idx, ]
      object@net$pval <- object@net$pval[idx, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    
  }
  return(object)
}

# load data ----
load(file="integrated_rpca_syn_PBMCAcuteHealth.Rda")
# pbmc.syn.combined
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
future::plan("multiprocess", workers = 16) 

ccCase01_A <- readRDS(file="cellchat.Case01.A_acute.rds")
ccCase02_A <- readRDS(file="cellchat.Case02.A_acute.rds")
ccCase03_A <- readRDS(file="cellchat.Case03.A_acute.rds")
ccCase04_A <- readRDS(file="cellchat.Case04.A_acute.rds")
ccCase05_A <- readRDS(file="cellchat.Case05.A_acute.rds")
ccCase06_A <- readRDS(file="cellchat.Case06.A_acute.rds")
ccCase07_A <- readRDS(file="cellchat.Case07.A_acute.rds")
#saveRDS(cellchat, file = paste("cellchat", case, group, "rds", sep="."))

case <- "Case01"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase01_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase01_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()


case <- "Case02"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase02_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase02_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "Case03"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase03_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase03_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "Case04"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase04_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase04_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "Case05"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase05_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase05_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "Case06"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase06_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase06_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "Case07"
group <- "A_acute"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccCase07_A, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccCase07_A, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

ccMH_H <- readRDS(file="cellchat.MH.Health.rds")
ccRG_H <- readRDS(file="cellchat.RG.Health.rds")
ccYB_H <- readRDS(file="cellchat.YB.Health.rds")
ccN05S_H <- readRDS(file="cellchat.N05S.Health.rds")
ccN01F_H <- readRDS(file="cellchat.N01F.Health.rds")

case <- "MH"
group <- "Health"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccMH_H, pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccMH_H, pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()


case <- "RG"
group <- "Health"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccRG_H, 
                                         pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccRG_H, 
                                         pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "YB"
group <- "Health"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccYB_H, 
                                         pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccYB_H, 
                                         pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "N05S"
group <- "Health"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccN05S_H, 
                                         pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccN05S_H, 
                                         pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

case <- "N01F"
group <- "Health"
png(paste("netAnalysis_signalingRole_heatmap", case, group, "png", sep="."),
    width = 3600, height = 2800, res=300)
ht1 <- netAnalysis_signalingRole_heatmap(ccN01F_H, 
                                         pattern = "outgoing",
                                         width = 12,
                                         height =15,)
ht2 <- netAnalysis_signalingRole_heatmap(ccN01F_H, 
                                         pattern = "incoming",
                                         width = 12,
                                         height =15,)

draw(ht1+ht2, column_title=paste(case, group, sep=" "))
dev.off()

# select cells to chrod plot



order_of_cells <- levels(ccCase01_A@idents)
names(order_of_cells) <- seq(1:length(order_of_cells))
select_cell <- c(
  "Naive CD4 T cells",
  "Classical monocytes",
  "Intermediate monocytes",
  "Plasmacytoid dendritic cells"
)
names(order_of_cells[order_of_cells %in% select_cell])
vertex.receiver = c()
netVisual_bubble(ccCase01_A, 
                 sources.use = c(1, 17, 18, 21), 
                 targets.use = c(1, 17, 18, 21), 
                 remove.isolate = FALSE)
netVisual_bubble(ccCase01_A, 
                 sources.use = c(1), 
                 targets.use = c(1, 17, 18, 21), 
                 remove.isolate = FALSE)
ccCase01_A@netP$pathways
pathways.show = "FLT3"
vertex.receiver = c(1, 17, 18, 21)
netVisual_individual(ccCase01_A, 
                     signaling = pathways.show, 
                     #pairLR.use = LR.show, 
                     #layout = "chord",
                     vertex.receiver=vertex.receiver)
pdf(file ="temp.pdf", width = 20, height =16)
#netVisual_aggregate()
netVisual_aggregate(ccCase01_A, 
                    signaling = pathways.show, 
                    layout = "chord",
                    vertex.receiver=vertex.receiver
)
dev.off()



make_source_dot <- function(ccobj, pathway_list, case, cell, id){
  netVisual_bubble(ccobj, 
                   sources.use = id, 
                   targets.use = c(1, 17, 18, 21), 
                   signaling = pathway_list, 
                   remove.isolate = FALSE)
  ggsave(paste("dotPlot_source",cell, case, "png", sep="."), 
         width = 4, height = 6)
}
make_target_dot <- function(ccobj, pathway_list, case, cell, id){
  netVisual_bubble(ccobj, 
                   sources.use = c(1, 17, 18, 21),
                   targets.use = id, 
                   signaling = pathway_list, 
                   remove.isolate = FALSE)
  ggsave(paste("dotPlot_target",cell, case, "png", sep="."), 
         width = 4, height = 6)
}
NaiveCD4T_L <- c(
  "IL16",  "MIF", "SEMA4", "MHC-I", "SELPLG", "ANNEXIN",
  "CD99", "ICAM", "CD40", "CD6", "FLT3"
)
cell="NaiveCD4T"
make_source_dot(ccCase01_A, NaiveCD4T_L, "Case01", cell, 1)
make_source_dot(ccCase02_A, NaiveCD4T_L, "Case02", cell, 1)
make_source_dot(ccCase03_A, NaiveCD4T_L, "Case03", cell, 1)
make_source_dot(ccCase04_A, NaiveCD4T_L, "Case04", cell, 1)
make_source_dot(ccCase05_A, NaiveCD4T_L, "Case05", cell, 1)
make_source_dot(ccCase06_A, NaiveCD4T_L, "Case06", cell, 1)
make_source_dot(ccCase07_A, NaiveCD4T_L, "Case07", cell, 1)
make_source_dot(ccMH_H, NaiveCD4T_L, "MH", cell, 1)
make_source_dot(ccRG_H, NaiveCD4T_L, "RG", cell, 1)
make_source_dot(ccYB_H, NaiveCD4T_L, "YB", cell, 1)
make_source_dot(ccN01F_H, NaiveCD4T_L, "N01F", cell, 1)
make_source_dot(ccN05S_H, NaiveCD4T_L, "N05S", cell, 1)

# CD14Mono
CD14Mono_L <- c(
  "GALECTIN",  "APP", "MHC-II", "ADGRE5", "RESISTIN",
  "IL16", "PECAM1", "SELPLG", "SEMA4", "MHC-I",
  "CD99", "ITGB2", "ANNEXIN", "VISFATIN", "GRN"
)
cell = "CD14Mono"
make_source_dot(ccCase01_A, CD14Mono_L, "Case01", cell, 17)
make_source_dot(ccCase02_A, CD14Mono_L, "Case02", cell, 17)
make_source_dot(ccCase03_A, CD14Mono_L, "Case03", cell, 17)
make_source_dot(ccCase04_A, CD14Mono_L, "Case04", cell, 17)
make_source_dot(ccCase05_A, CD14Mono_L, "Case05", cell, 17)
make_source_dot(ccCase06_A, CD14Mono_L, "Case06", cell, 17)
make_source_dot(ccCase07_A, CD14Mono_L, "Case07", cell, 17)
make_source_dot(ccMH_H, CD14Mono_L, "MH", cell, 17)
make_source_dot(ccRG_H, CD14Mono_L, "RG", cell, 17)
make_source_dot(ccYB_H, CD14Mono_L, "YB", cell, 17)
make_source_dot(ccN01F_H, CD14Mono_L, "N01F", cell, 17)
make_source_dot(ccN05S_H, CD14Mono_L, "N05S", cell, 17)

# InterMono
InterMono_L <- c(
  "ADGRE5", "MHC-II", "GALECTIN", "IL16", "SELPLG",
  "SEMA4", "PECAM1", "MHC-I", "MIF", "ICAM",
  "RESISTIN", "ITGB2", "APP", "ANNEXIN", "CD99",
  "VISFATIN", "GRN"
)
cell = "InterMono"
make_source_dot(ccCase01_A, InterMono_L, "Case01", cell, 18)
make_source_dot(ccCase02_A, InterMono_L, "Case02", cell, 18)
make_source_dot(ccCase03_A, InterMono_L, "Case03", cell, 18)
make_source_dot(ccCase04_A, InterMono_L, "Case04", cell, 18)
make_source_dot(ccCase05_A, InterMono_L, "Case05", cell, 18)
make_source_dot(ccCase06_A, InterMono_L, "Case06", cell, 18)
make_source_dot(ccCase07_A, InterMono_L, "Case07", cell, 18)
make_source_dot(ccMH_H, InterMono_L, "MH", cell, 18)
make_source_dot(ccRG_H, InterMono_L, "RG", cell, 18)
make_source_dot(ccYB_H, InterMono_L, "YB", cell, 18)
make_source_dot(ccN01F_H, InterMono_L, "N01F", cell, 18)
make_source_dot(ccN05S_H, InterMono_L, "N05S", cell, 18)

pDC_L <- c(
  "APP", "MHC-II", "IL16", "SELPLG", "SEMA4",
  "PECAM1", "MHC-I", "CD99", "BTLA", "MIF",
  "GRN", "GALECTIN", "VISFATIN", "ALCAM"
)
cell = "pDC"
make_source_dot(ccCase01_A, pDC_L, "Case01", cell, 21)
make_source_dot(ccCase02_A, pDC_L, "Case02", cell, 21)
make_source_dot(ccCase03_A, pDC_L, "Case03", cell, 21)
make_source_dot(ccCase04_A, pDC_L, "Case04", cell, 21)
make_source_dot(ccCase05_A, pDC_L, "Case05", cell, 21)
make_source_dot(ccCase06_A, pDC_L, "Case06", cell, 21)
make_source_dot(ccCase07_A, pDC_L, "Case07", cell, 21)
make_source_dot(ccMH_H, pDC_L, "MH", cell, 21)
make_source_dot(ccRG_H, pDC_L, "RG", cell, 21)
make_source_dot(ccYB_H, pDC_L, "YB", cell, 21)
make_source_dot(ccN01F_H, pDC_L, "N01F", cell, 21)
make_source_dot(ccN05S_H, pDC_L, "N05S", cell, 21)

# target
## NaiveCD4T
NaiveCD4T_R <- c(
  "SELPLG", "MHC-II", "APP", "IL16", "GALECTIN",
  "ADGRE5", "ITGB2", "RESISTIN", "ALCAM"
)
cell = "NaiveCD4T"
make_target_dot(ccCase01_A, NaiveCD4T_R, "Case01", cell, 1)
make_target_dot(ccCase02_A, NaiveCD4T_R, "Case02", cell, 1)
make_target_dot(ccCase03_A, NaiveCD4T_R, "Case03", cell, 1)
make_target_dot(ccCase04_A, NaiveCD4T_R, "Case04", cell, 1)
make_target_dot(ccCase05_A, NaiveCD4T_R, "Case05", cell, 1)
make_target_dot(ccCase06_A, NaiveCD4T_R, "Case06", cell, 1)
make_target_dot(ccCase07_A, NaiveCD4T_R, "Case07", cell, 1)
make_target_dot(ccMH_H, NaiveCD4T_R, "MH", cell, 1)
make_target_dot(ccRG_H, NaiveCD4T_R, "RG", cell, 1)
make_target_dot(ccYB_H, NaiveCD4T_R, "YB", cell, 1)
make_target_dot(ccN01F_H, NaiveCD4T_R, "N01F", cell, 1)
make_target_dot(ccN05S_H, NaiveCD4T_R, "N05S", cell, 1)

## CD14Mono
CD14Mono_R <- c(
  "CD99", "MHC-I", "SEMA4", "IL16", "SELPLG",
  "ANNEXIN", "MHC-II", "MIF", "PECAM1", "GALECTIN",
  "VISFATIN", "APP", "ADGRE5", "ICAM", "GRN", 
  "RESISTIN"
)
cell = "CD14Mono"
make_target_dot(ccCase01_A, CD14Mono_R, "Case01", cell, 17)
make_target_dot(ccCase02_A, CD14Mono_R, "Case02", cell, 17)
make_target_dot(ccCase03_A, CD14Mono_R, "Case03", cell, 17)
make_target_dot(ccCase04_A, CD14Mono_R, "Case04", cell, 17)
make_target_dot(ccCase05_A, CD14Mono_R, "Case05", cell, 17)
make_target_dot(ccCase06_A, CD14Mono_R, "Case06", cell, 17)
make_target_dot(ccCase07_A, CD14Mono_R, "Case07", cell, 17)
make_target_dot(ccMH_H, CD14Mono_R, "MH", cell, 17)
make_target_dot(ccRG_H, CD14Mono_R, "RG", cell, 17)
make_target_dot(ccYB_H, CD14Mono_R, "YB", cell, 17)
make_target_dot(ccN01F_H, CD14Mono_R, "N01F", cell, 17)
make_target_dot(ccN05S_H, CD14Mono_R, "N05S", cell, 17)

# InterMono
InterMono_R <- c(
  "CD99", "MHC-I", "SEMA4", "IL16", "ANNEXIN",
  "MHC-II", "PECAM1", "ITGB2", "MIF", "GRN",
  "GALECTIN", "APP", "ICAM", "ADGRE5", "RESISTIN",
  "CD48", "VISFATIN"
)
cell = "InterMono"
make_target_dot(ccCase01_A, InterMono_R, "Case01", cell, 18)
make_target_dot(ccCase02_A, InterMono_R, "Case02", cell, 18)
make_target_dot(ccCase03_A, InterMono_R, "Case03", cell, 18)
make_target_dot(ccCase04_A, InterMono_R, "Case04", cell, 18)
make_target_dot(ccCase05_A, InterMono_R, "Case05", cell, 18)
make_target_dot(ccCase06_A, InterMono_R, "Case06", cell, 18)
make_target_dot(ccCase07_A, InterMono_R, "Case07", cell, 18)
make_target_dot(ccMH_H, InterMono_R, "MH", cell, 18)
make_target_dot(ccRG_H, InterMono_R, "RG", cell, 18)
make_target_dot(ccYB_H, InterMono_R, "YB", cell, 18)
make_target_dot(ccN01F_H, InterMono_R, "N01F", cell, 18)
make_target_dot(ccN05S_H, InterMono_R, "N05S", cell, 18)

# pDC

pDC_R <- c(
  "IL16", "SELPLG", "SEMA4", "MHC-II", "PECAM1",
  "APP", "MIF", "GALECTIN", "MHC-I", "ADGRE5",
  "VISFATIN", "RESISTIN", "FLT3", "CD6"
)
cell="pDC"
make_target_dot(ccCase01_A, pDC_R, "Case01", cell, 21)
make_target_dot(ccCase02_A, pDC_R, "Case02", cell, 21)
make_target_dot(ccCase03_A, pDC_R, "Case03", cell, 21)
make_target_dot(ccCase04_A, pDC_R, "Case04", cell, 21)
make_target_dot(ccCase05_A, pDC_R, "Case05", cell, 21)
make_target_dot(ccCase06_A, pDC_R, "Case06", cell, 21)
make_target_dot(ccCase07_A, pDC_R, "Case07", cell, 21)
make_target_dot(ccMH_H, pDC_R, "MH", cell, 21)
make_target_dot(ccRG_H, pDC_R, "RG", cell, 21)
make_target_dot(ccYB_H, pDC_R, "YB", cell, 21)
make_target_dot(ccN01F_H, pDC_R, "N01F", cell, 21)
make_target_dot(ccN05S_H, pDC_R, "N05S", cell, 21)
