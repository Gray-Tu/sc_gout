setwd("With_SYN")

#paste("cellchat", case, group, "rds", sep=".")
CCIc1 <- readRDS("cellchat.Case01.A_acute.rds")
CCIc4 <- readRDS("cellchat.Case04.A_acute.rds")
CCIc2 <- readRDS("cellchat.Case02.A_acute.rds")
CCIc5 <- readRDS("cellchat.Case05.A_acute.rds")
CCIc7 <- readRDS("cellchat.Case07.A_acute.rds")
CCIN01F <- readRDS("cellchat.N01F.Health.rds")
CCIYB <- readRDS("cellchat.YB.Health.rds")
CCIRG <- readRDS("cellchat.RG.Health.rds")

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
# fix typo
# Megakaryoctes to Megakaryocytes
fix_typo <- function(cci){
  my_DF <- cci@meta %>% 
    mutate_all(as.character) %>% 
    mutate(cell_type = ifelse(group=="Megakaryoctes",
                              "Megakaryocytes",
                              group))
  cci@meta$cell_type <- my_DF$cell_type
  levels(cci@meta$group)[levels(cci@meta$group) =='Megakaryoctes' ] <- 'Megakaryocytes'
  cci@meta$cell_type <- factor(my_DF$cell_type,
                               levels=levels(cci@meta$group))
  cci_fix <- setIdent(cci, ident.use = "cell_type")
  return(cci_fix)
}
CCIc1 <- fix_typo(CCIc1)
CCIc2 <- fix_typo(CCIc2)
CCIc4 <- fix_typo(CCIc4)
CCIc5 <- fix_typo(CCIc5)
CCIc7 <- fix_typo(CCIc7)
CCIN01F <- fix_typo(CCIN01F)
CCIYB <- fix_typo(CCIYB)
CCIRG <- fix_typo(CCIRG)
# relabel need re-run those step, https://github.com/sqjin/CellChat/issues/35
# CellChatDB <- CellChatDB.human
# CellChatDB.use <- CellChatDB
re_run_CCI <- function(cellchat, case, group){
  #cellchat@DB <- CellChatDB.use
  #cellchat <- subsetData(cellchat)
  #cellchat <- identifyOverExpressedGenes(cellchat)
  #cellchat <- identifyOverExpressedInteractions(cellchat)
  #cellchat <- projectData(cellchat, PPI.human)
  #cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  df.net <- subsetCommunication(cellchat) # tables
  df.netP <- subsetCommunication(cellchat, slot.name = "netP") # tables
  write.csv(df.net, file=paste("fixCellChat_net", case, group, "csv", sep="."))
  write.csv(df.netP, file=paste("fixCellChat_netP", case, group, "csv", sep="."))
  
  return(cellchat)
}
CCIc1 <- re_run_CCI(CCIc1, "Case01", "A_acute")
CCIc2 <- re_run_CCI(CCIc2, "Case02", "A_acute")
CCIc4 <- re_run_CCI(CCIc4, "Case04", "A_acute")
CCIc5 <- re_run_CCI(CCIc5, "Case05", "A_acute")
CCIc7 <- re_run_CCI(CCIc7, "Case07", "A_acute")
CCIN01F <- re_run_CCI(CCIN01F, "N01F", "Health")
CCIYB <- re_run_CCI(CCIYB, "YB", "Health")
CCIRG <- re_run_CCI(CCIRG, "RG", "Health")

#pathways.show = "IL1" ----
pdf("ADV_FigCCI_IL1_C1C4_TNF_C2C5C7YB.pdf", height = 20, width = 30)
par(mfrow = c(2, 3), xpd=NA)

{ # IL1
  pathways.show = "IL1"
  netVisual_chord_cell(CCIc1, 
                    signaling = pathways.show,  
                    title.name ="",
                    lab.cex=1.5

                    
                    )
  fig_label("a", cex=4, font=2, region="plot")
  fig_label("Case01-IL1", cex=4, region="plot", pos="bottom")
  netVisual_chord_cell(CCIc4, 
                    signaling = pathways.show,  
                    title.name ="",
                    lab.cex=1.5

                    
  )
  fig_label("b", cex=4, font=2, region="plot")
  fig_label("Case04-IL1", cex=4, region="plot", pos="bottom")
  
  pathways.show = "TNF"
  netVisual_chord_cell(CCIc2, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("c", cex=4, font=2, region="plot")
  fig_label("Case02-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc5, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("d", cex=4, font=2, region="plot")
  fig_label("Case05-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc7, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("e", cex=4, font=2, region="plot")
  fig_label("Case07-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIYB, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("f", cex=4, region="plot")
  fig_label("YB-TNF", cex=4, region="plot", pos="bottom")
  }
dev.off()
# # png ----
# png("ADV_FigCCI_IL1_C1C4_TNF_C2C5C7YB.png", unit="in", res=300, height = 20, width = 30)
# par(mfrow = c(2, 3), xpd=NA)
# 
# { # IL1
#   pathways.show = "IL1"
#   netVisual_chord_cell(CCIc1, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("a", cex=4, region="plot")
#   fig_label("Case01-IL1", cex=4, region="plot", pos="bottom")
#   netVisual_chord_cell(CCIc4, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("b", cex=4, region="plot")
#   fig_label("Case04-IL1", cex=4, region="plot", pos="bottom")
#   
#   pathways.show = "TNF"
#   netVisual_chord_cell(CCIc2, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("c", cex=4, region="plot")
#   fig_label("Case02-TNF", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc5, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("d", cex=4, region="plot")
#   fig_label("Case05-TNF", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc7, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("e", cex=4, region="plot")
#   fig_label("Case07-TNF", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIYB, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        
#                        
#   )
#   fig_label("f", cex=4, region="plot")
#   fig_label("YB-TNF", cex=4, region="plot", pos="bottom")
# }
# dev.off()

# emf ----
require(devEMF)
emf("ADV_FigCCI_IL1_C1C4_TNF_C2C5C7YB.emf", height = 20, width = 30)
par(mfrow = c(2, 3), xpd=NA)

{ # IL1
  pathways.show = "IL1"
  netVisual_chord_cell(CCIc1, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("a", cex=4, font=2,  region="plot")
  fig_label("Case01-IL1", cex=4, region="plot", pos="bottom")
  netVisual_chord_cell(CCIc4, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("b", cex=4, font=2, region="plot")
  fig_label("Case04-IL1", cex=4, region="plot", pos="bottom")
  
  pathways.show = "TNF"
  netVisual_chord_cell(CCIc2, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("c", cex=4, font=2, region="plot")
  fig_label("Case02-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc5, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("d", cex=4, font=2, region="plot")
  fig_label("Case05-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc7, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("e", cex=4, font=2, region="plot")
  fig_label("Case07-TNF", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIYB, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
                       
                       
  )
  fig_label("f", cex=4, font=2, region="plot")
  fig_label("YB-TNF", cex=4, region="plot", pos="bottom")
}
dev.off()

# pdf("ADV_FigCCI_TGFB1_C1C2C4C5C7N01F.pdf", height = 20, width = 30)
# par(mfrow = c(2, 3), xpd=NA)
# 
# { 
#   pathways.show = "TGFb"
#   netVisual_chord_cell(CCIc1, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("A", cex=4, region="plot")
#   fig_label("Case01-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc2, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("B", cex=4, region="plot")
#   fig_label("Case02-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc4, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("C", cex=4, region="plot")
#   fig_label("Case04-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc5, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("D", cex=4, region="plot")
#   fig_label("Case05-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc7, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("E", cex=4, region="plot")
#   fig_label("Case07-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIN01F, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#   )
#   fig_label("F", cex=4, region="plot")
#   fig_label("N01F-TGFb", cex=4, region="plot", pos="bottom")
# }
# dev.off()
# 
# png("ADV_FigCCI_TGFB1_C1C2C4C5C7N01F.png", unit="in", res=300, height = 20, width = 30)
# par(mfrow = c(2, 3), xpd=NA)
# 
# { 
#   pathways.show = "TGFb"
#   netVisual_chord_cell(CCIc1, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("A", cex=4, region="plot")
#   fig_label("Case01-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc2, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("B", cex=4, region="plot")
#   fig_label("Case02-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc4, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("C", cex=4, region="plot")
#   fig_label("Case04-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc5, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("D", cex=4, region="plot")
#   fig_label("Case05-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIc7, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("E", cex=4, region="plot")
#   fig_label("Case07-TGFb", cex=4, region="plot", pos="bottom")
#   
#   netVisual_chord_cell(CCIN01F, 
#                        signaling = pathways.show,  
#                        title.name ="",
#                        lab.cex=1.3
#                        )
#   fig_label("F", cex=4, region="plot")
#   fig_label("N01F-TGFb", cex=4, region="plot", pos="bottom")
# }
# dev.off()


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  # https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

# TGFb ----
pdf("ADV_FigCCI_TGFB1_C1C2C4C5C7N01FRG.pdf", height = 20, width = 40)
par(mfrow = c(2, 4), xpd=NA)

{ 
  pathways.show = "TGFb"
  netVisual_chord_cell(CCIc1, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("a", cex=4, font=2, region="plot")
  fig_label("Case01-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc2, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("b", cex=4, font=2, region="plot")
  fig_label("Case02-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc4, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("c", cex=4, font=2, region="plot")
  fig_label("Case04-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc5, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("d", cex=4, font=2, region="plot")
  fig_label("Case05-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc7, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("e", cex=4, font=2, region="plot")
  fig_label("Case07-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIN01F, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("f", cex=4, font=2, region="plot")
  fig_label("N01F-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIRG, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5 
  )
  fig_label("g", cex=4, region="plot")
  fig_label("RG-TGFb", cex=4, region="plot", pos="bottom")
  
}
dev.off()
# emf ----
emf("ADV_FigCCI_TGFB1_C1C2C4C5C7N01FRG.emf", height = 20, width = 30)
par(mfrow = c(2, 4), xpd=NA)

{ 
  pathways.show = "TGFb"
  netVisual_chord_cell(CCIc1, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("a", cex=4, font=2, region="plot")
  fig_label("Case01-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc2, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("b", cex=4, font=2, region="plot")
  fig_label("Case02-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc4, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("c", cex=4, font=2, region="plot")
  fig_label("Case04-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc5, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("d", cex=4, font=2, region="plot")
  fig_label("Case05-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIc7, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("e", cex=4, font=2, region="plot")
  fig_label("Case07-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIN01F, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("f", cex=4, font=2, region="plot")
  fig_label("N01F-TGFb", cex=4, region="plot", pos="bottom")
  
  netVisual_chord_cell(CCIRG, 
                       signaling = pathways.show,  
                       title.name ="",
                       lab.cex=1.5
  )
  fig_label("g", cex=4, font=2, region="plot")
  fig_label("RG-TGFb", cex=4, region="plot", pos="bottom")
  
}
dev.off()
