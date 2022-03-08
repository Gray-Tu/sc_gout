# PBMC cell cell interaction

setwd("With_SYN")
library(Seurat)
library(CellChat)
library(patchwork)
library(cowplot)
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

# GROUP only use A_acute, Health, 
# SF will break in CCI
#for (group in unique(pbmc.syn.combined$GROUP)){
for (group in ("A_acute", "Health")){
  g_set <- subset(pbmc.syn.combined, GROUP==group)
  for (case in unique(g_set$CASE)){
    c_set <- subset(g_set, CASE==case)
    Idents(c_set) <- c_set$cell_type
    labels <- Idents(c_set)
    data.input <- GetAssayData(c_set, assay = "RNA", slot = "data")
    meta <- data.frame(group = labels, row.names = names(labels))
    cellchat <- createCellChat(object = data.input,
                               meta = meta,
                               group.by = "group")
    # correct order!
    cell.levels <- unique(c_set$cell_type)[order(unique(c_set$cell_type))]
    
    cellchat <- setIdent(cellchat, 
                         ident.use = "group", 
                         levels = cell.levels)
    
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    df.net <- subsetCommunication(cellchat) # tables
    df.netP <- subsetCommunication(cellchat, slot.name = "netP") # tables
    write.csv(df.net, file=paste("CellChat_net", case, group, "csv", sep="."))
    write.csv(df.netP, file=paste("CellChat_netP", case, group, "csv", sep="."))
    mat <- cellchat@net$weight
    groupSize <- as.numeric(table(cellchat@idents))
    # P1 each cell type Plot Weight
    pdf(paste("CCI_weight",case, group, "pdf", sep="."), 
        height = 14, width = 14)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    # per cellchat@netP$pathways
    for (pathways.show in cellchat@netP$pathways){
      # Hierarchy plot
      vertex.receiver = seq(1,20) 
      pdf(file=paste("Hierarchy", pathways.show, case, group, "pdf", sep="."),
          height = 10, width = 20)
      netVisual_aggregate(cellchat, 
                          signaling = pathways.show,  
                          vertex.receiver = vertex.receiver,
                          layout="hierarchy")
      dev.off()
      
      # Circle plot
      pdf(file=paste("Circle", pathways.show, case, group, "pdf", sep="."),
          height = 14, width = 14)
      netVisual_aggregate(cellchat, 
                          signaling = pathways.show,  
                          vertex.receiver = vertex.receiver,
                          layout="circle")
      dev.off()
      
      # Chord plot
      pdf(file=paste("Chord", pathways.show, case, group, "pdf", sep="."),
          height = 10, width = 10)
      netVisual_aggregate(cellchat, 
                          signaling = pathways.show,  
                          vertex.receiver = vertex.receiver,
                          layout="chord")
      dev.off()
      # Heatmap
      pdf(file=paste("Heatmap_CellChat", pathways.show, case, group, "pdf", sep="."),
          height = 8, width = 8)
      p <- netVisual_heatmap(cellchat, 
                             signaling = pathways.show,
                             color.heatmap = "Reds")
      print(p)
      dev.off()
      
      # vlnPlot
      pdf(file=paste("VlnPlot_CellChat", pathways.show, case, group, "pdf", sep="."),
          height = 8, width = 8)
      p <- plotGeneExpression(cellchat, signaling = pathways.show)
      print(p)
      dev.off()
      # centrality scores
      pdf(file=paste("centrality_heatmap", pathways.show, case, group, "pdf", sep="."),
          width = 10, height = 3)
      netAnalysis_signalingRole_network(cellchat, 
                                        signaling = pathways.show, 
                                        width = 20, height = 3, font.size = 10)
      dev.off()
      
    }
    
    pdf(paste("netAnalysis_signalingRole_heatmap", case, group, "pdf", sep="."),
        width = 12, height = 10)
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                             width = 12,
                                             height =15,)
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                             width = 12,
                                             height =15,)
    print(ht1 + ht2)
    dev.off()
    
    saveRDS(cellchat, file = paste("cellchat", case, group, "rds", sep="."))
  }
}


