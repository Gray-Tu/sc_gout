# INF-score and ENG-score in 
# 1. domain_cls-Acute vs domain_cls-SF
# 2. cls-1-Acute vs cls-1-SF

setwd("With_SYN")

library(Seurat)
library(dplyr)
library(ggpubr)
library(cowplot)

load("MainDEGusage_ScoreAdd_20220211.rda")

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

boxplot_main_INF <- function(DF){
  main_stat.test <- compare_means(inflammatory_score1~GROUP, 
                                  DF, 
                                  method = "wilcox.test",
                                  paired = FALSE)
  
  group_y <- DF %>% 
    mutate(ymax = max(inflammatory_score1)) %>% 
    dplyr::select(ymax) %>%
    unique()
  
  groupy_list <- list()
  for (g in seq(1:length(rownames(main_stat.test)))){
    groupy_list[[g]] = as.numeric(group_y["ymax"] * (1 + 0.1 *(g%%3 +1)))
  }
  
  main_stat.test$y.position <- as.numeric(groupy_list)
  main_stat.test <- main_stat.test %>% subset(p < 0.05)
  p <- ggboxplot(DF,           
                 x="GROUP", 
                 y="inflammatory_score1",
                 fill="GROUP")+
    labs(y="inflammatory score", fill="GROUP")+
    theme_classic()
  
  INF_main_p <- p + stat_pvalue_manual(main_stat.test, 
                                       label = "p.signif", 
                                       tip.length=0, vjust=0.4)
  
  return(INF_main_p)
}
# main ENG
boxplot_main_ENG <- function(DF){
  main_stat.test <- compare_means(energy_score1~GROUP, 
                                  DF, 
                                  method = "wilcox.test",
                                  paired = FALSE)
  
  group_y <- DF %>% 
    mutate(ymax = max(energy_score1)) %>% 
    dplyr::select(ymax) %>%
    unique()
  
  groupy_list <- list()
  for (g in seq(1:length(rownames(main_stat.test)))){
    groupy_list[[g]] = as.numeric(group_y["ymax"] * (1 + 0.1 *(g%%3 +1)))
  }
  
  main_stat.test$y.position <- as.numeric(groupy_list)
  main_stat.test <- main_stat.test %>% subset(p < 0.05)
  p <- ggboxplot(DF,           
                 x="GROUP", 
                 y="energy_score1",
                 fill="GROUP")+
    labs(y="energy score", fill="GROUP")+
    theme_classic()
  
  ENG_main_p <- p + stat_pvalue_manual(main_stat.test, 
                                       label = "p.signif", 
                                       tip.length=0, vjust=0.4)
  
  return(ENG_main_p)
}

boxplot_cls_INF <- function(DF){
  INF_cls_wlicon <- compare_means(inflammatory_score1~seurat_clusters, 
                                  DF, 
                                  method = "wilcox.test", 
                                  paired = FALSE) %>%
    subset(p<0.05)
  
  mylist <- list()
  for (i in c(1:length(rownames(INF_cls_wlicon)))){
    g1 <- as.character(INF_cls_wlicon[i,]$group1)
    g2 <- as.character(INF_cls_wlicon[i,]$group2)
    mylist[[i]] <-c(g1, g2)
  }
  
  INF_cls_p <- DF %>% 
    ggboxplot(x="seurat_clusters", 
              y="inflammatory_score1", 
              fill="seurat_clusters"
    )+
    labs(y="inflammatory score", x="", fill="subclusters")+
    theme_classic2()+
    stat_compare_means(comparisons = mylist,
                       label = "p.signif",
                       tip.length=0,
                       vjust=0.4,
                       method = "wilcox.test",
    )
  
  return(INF_cls_p)
}
boxplot_cls_ENG <- function(DF){
  ENG_cls_wlicon <- compare_means(energy_score1~seurat_clusters, 
                                  DF, 
                                  method = "wilcox.test", 
                                  paired = FALSE) %>%
    subset(p<0.05)
  
  mylist <- list()
  for (i in c(1:length(rownames(ENG_cls_wlicon)))){
    g1 <- as.character(ENG_cls_wlicon[i,]$group1)
    g2 <- as.character(ENG_cls_wlicon[i,]$group2)
    mylist[[i]] <-c(g1, g2)
  }
  
  ENG_cls_p <- DF %>% 
    ggboxplot(x="seurat_clusters", 
              y="energy_score1", 
              fill="seurat_clusters"
    )+
    labs(y="energy score", x="", fill="subclusters")+
    theme_classic2()+
    stat_compare_means(comparisons = mylist,
                       label = "p.signif",
                       tip.length=0,
                       vjust=0.4,
                       method = "wilcox.test",
    )
  
  return(ENG_cls_p)
}

get_test_table <- function(DF, cell_name){
  # main test 
  INF_main <- compare_means(inflammatory_score1~GROUP,
                            DF,
                            method = "wilcox.test",
                            paired = FALSE)
  ENG_main <- compare_means(energy_score1~GROUP,
                            DF,
                            method = "wilcox.test",
                            paired = FALSE)
  # cls test
  INF_cls <- compare_means(inflammatory_score1~seurat_clusters,
                           DF,
                           method = "wilcox.test",
                           paired = FALSE)
  ENG_cls <- compare_means(energy_score1~seurat_clusters,
                           DF,
                           method = "wilcox.test",
                           paired = FALSE)
  # cls-group test
  INF_clsgroup <- compare_means(inflammatory_score1~GROUP,
                                DF,
                                method = "wilcox.test",
                                paired = FALSE,
                                group.by = "seurat_clusters")
  
  ENG_clsgroup <- compare_means(energy_score1~GROUP,
                                DF,
                                method = "wilcox.test",
                                paired = FALSE,
                                group.by = "seurat_clusters")
  
  write.csv(INF_main, file=paste("ADV_fig05_INF_group", cell_name, "csv", sep="."))
  write.csv(ENG_main, file=paste("ADV_fig05_ENG_group", cell_name, "csv", sep="."))
  write.csv(INF_cls, file=paste("ADV_fig05_INF_cls", cell_name, "csv", sep="."))
  write.csv(ENG_cls, file=paste("ADV_fig05_ENG_cls", cell_name, "csv", sep="."))
  write.csv(INF_clsgroup, file=paste("ADV_fig05_INF_clsgroup", cell_name, "csv", sep="."))
  write.csv(ENG_clsgroup, file=paste("ADV_fig05_ENG_clsgroup", cell_name, "csv", sep="."))
  
}
#1. Intermediate monocytes (from Fig2B)
InterMono <- add_ENGscore(InterMono, total_genes)
InterMono <- add_INFscore(InterMono, INFgene_df)
InterMono$GROUP <- factor(InterMono$GROUP, levels=c("Synovial fluid", "A_acute", "Health"))
INF_dim <- FeaturePlot(InterMono, features="inflammatory_score1", label=T)+labs(title="Inflammatory score")
ENG_dim <- FeaturePlot(InterMono, features="energy_score1", label=T)+labs(title="Energy score")
DF <- InterMono[[]]
INF_main_p <- boxplot_main_INF(InterMono[[]])
ENG_main_p <- boxplot_main_ENG(InterMono[[]])
INF_cls_p <- boxplot_cls_INF(InterMono[[]])
ENG_cls_p <- boxplot_cls_ENG(InterMono[[]])

### INF sub_group----
stat.test <- compare_means(
  inflammatory_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.25, .3, .35, .4, .45, .35, .35))

INF_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
           y="inflammatory_score1", 
           fill="GROUP"
           
  )+
  labs(y="inflammatory score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.1, 0.5)

### ENG sub_group----
stat.test <- compare_means(
  energy_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.12, .12, .15, .12, .15, .12, .15))

ENG_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
           y="energy_score1", 
           fill="GROUP"
  )+
  labs(y="energy score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.05,0.18)

group_le <- get_legend(INF_main_p+theme(legend.text.align = 0))
cls_le <- get_legend(INF_cls_p+theme(legend.text.align = 0))
legend_p <- plot_grid(cls_le, group_le, ncol=1)

l1 <- plot_grid(INF_dim, INF_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)), 
                ENG_dim, ENG_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)),
                group_le, rel_widths = c(1.2,0.8,1.2,0.8,0.6),
                align = "h", axis="bt", ncol=5, labels=c("A", "B", "C", "D"))

b2 <- plot_grid(INF_cls_p+NoLegend(), 
                INF_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("E", "G"))

b3 <- plot_grid(ENG_cls_p+NoLegend(), ENG_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("F", "H"))

l23 <- plot_grid(b2, b3, legend_p, ncol=3, rel_widths = c(1,1,0.3))
plot_grid(l1, l23, ncol=1, rel_heights = c(1,2))
ggsave("ADV_Fig05_INF_ENG_InterMono.png", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_InterMono.pdf", width = 12, height = 8)
get_test_table(DF, "InterMono")


#2. MAIT cells (from Fig2N)
#MAIT
MAIT <- add_ENGscore(MAIT, total_genes)
MAIT <- add_INFscore(MAIT, INFgene_df)
MAIT$GROUP <- factor(MAIT$GROUP, levels=c("Synovial fluid", "A_acute", "Health"))
INF_dim <- FeaturePlot(MAIT, features="inflammatory_score1", label=T)+labs(title="Inflammatory score")
ENG_dim <- FeaturePlot(MAIT, features="energy_score1", label=T)+labs(title="Energy score")
DF <- MAIT[[]]
INF_main_p <- boxplot_main_INF(DF)
ENG_main_p <- boxplot_main_ENG(DF)
INF_cls_p <- boxplot_cls_INF(DF)
ENG_cls_p <- boxplot_cls_ENG(DF)
### INF sub_group----
stat.test <- compare_means(
  inflammatory_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.15, .18, .15, .18))

INF_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="inflammatory_score1", 
            fill="GROUP"
            
  )+
  labs(y="inflammatory score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.05, 0.2)

### ENG sub_group----
stat.test <- compare_means(
  energy_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.07, .07))

ENG_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="energy_score1", 
            fill="GROUP"
  )+
  labs(y="energy score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.03,0.08)

group_le <- get_legend(INF_main_p+theme(legend.text.align = 0))
cls_le <- get_legend(INF_cls_p+theme(legend.text.align = 0))
legend_p <- plot_grid(cls_le, group_le, ncol=1)

l1 <- plot_grid(INF_dim, INF_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)), 
                ENG_dim, ENG_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)),
                group_le, rel_widths = c(1.2,0.8,1.2,0.8,0.6),
                align = "h", axis="bt", ncol=5, labels=c("A", "B", "C", "D"))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

b2 <- plot_grid(INF_cls_p+NoLegend(), 
                INF_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("E", "G"))

b3 <- plot_grid(ENG_cls_p+NoLegend(), ENG_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("F", "H"))

l23 <- plot_grid(b2, b3, legend_p, ncol=3, rel_widths = c(1,1,0.3))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
plot_grid(l1, l23, ncol=1, rel_heights = c(1,2))
ggsave("ADV_Fig05_INF_ENG_MAIT.png", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_MAIT.pdf", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_MAIT.jpg", width = 12, height = 8, dpi=600)
get_test_table(DF, "MAIT")

#3. T regulatory cells
Treg <- add_ENGscore(Treg, total_genes)
Treg <- add_INFscore(Treg, INFgene_df)
Treg$GROUP <- factor(Treg$GROUP, levels=c("Synovial fluid", "A_acute", "Health"))
INF_dim <- FeaturePlot(Treg, features="inflammatory_score1", label=T)+labs(title="Inflammatory score")
ENG_dim <- FeaturePlot(Treg, features="energy_score1", label=T)+labs(title="Energy score")
DF <- Treg[[]]
INF_main_p <- boxplot_main_INF(DF)
ENG_main_p <- boxplot_main_ENG(DF)
INF_cls_p <- boxplot_cls_INF(DF)
ENG_cls_p <- boxplot_cls_ENG(DF)
### INF sub_group----
stat.test <- compare_means(
  inflammatory_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.12, .14, .12))

INF_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="inflammatory_score1", 
            fill="GROUP"
            
  )+
  labs(y="inflammatory score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.07, 0.15)

### ENG sub_group----
stat.test <- compare_means(
  energy_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.09))

ENG_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="energy_score1", 
            fill="GROUP"
  )+
  labs(y="energy score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.03,0.1)

group_le <- get_legend(INF_main_p+theme(legend.text.align = 0))
cls_le <- get_legend(INF_cls_p+theme(legend.text.align = 0))
legend_p <- plot_grid(cls_le, group_le, ncol=1)

l1 <- plot_grid(INF_dim, INF_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)), 
                ENG_dim, ENG_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)),
                group_le, rel_widths = c(1.2,0.8,1.2,0.8,0.6),
                align = "h", axis="bt", ncol=5, labels=c("A", "B", "C", "D"))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

b2 <- plot_grid(INF_cls_p+NoLegend(), 
                INF_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("E", "G"))

b3 <- plot_grid(ENG_cls_p+NoLegend(), ENG_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("F", "H"))

l23 <- plot_grid(b2, b3, legend_p, ncol=3, rel_widths = c(1,1,0.3))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
plot_grid(l1, l23, ncol=1, rel_heights = c(1,2))
ggsave("ADV_Fig05_INF_ENG_Treg.png", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_Treg.pdf", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_Treg.jpg", width = 12, height = 8, dpi=600)
get_test_table(DF, "Treg")

# other cell 
#4. CD14Mono
CD14Mono <- add_ENGscore(CD14Mono, total_genes)
CD14Mono <- add_INFscore(CD14Mono, INFgene_df)
CD14Mono$GROUP <- factor(CD14Mono$GROUP, levels=c("Synovial fluid", "A_acute", "Health"))
INF_dim <- FeaturePlot(CD14Mono, features="inflammatory_score1", label=T)+labs(title="Inflammatory score")
ENG_dim <- FeaturePlot(CD14Mono, features="energy_score1", label=T)+labs(title="Energy score")
DF <- CD14Mono[[]]
INF_main_p <- boxplot_main_INF(DF)
ENG_main_p <- boxplot_main_ENG(DF)
INF_cls_p <- boxplot_cls_INF(DF)
ENG_cls_p <- boxplot_cls_ENG(DF)
### INF sub_group----
stat.test <- compare_means(
  inflammatory_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.3, .3, .34, .38, .3, .34,
                        .38, .42, .46, .42, .46))

INF_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="inflammatory_score1", 
            fill="GROUP"
            
  )+
  labs(y="inflammatory score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.1, .48)

### ENG sub_group----
stat.test <- compare_means(
  energy_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.12, .12, .12, .14, .12, .14, .12))

ENG_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="energy_score1", 
            fill="GROUP"
  )+
  labs(y="energy score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.06,0.15)

group_le <- get_legend(INF_main_p+theme(legend.text.align = 0))
cls_le <- get_legend(INF_cls_p+theme(legend.text.align = 0))
legend_p <- plot_grid(cls_le, group_le, ncol=1)

l1 <- plot_grid(INF_dim, INF_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)), 
                ENG_dim, ENG_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)),
                group_le, rel_widths = c(1.2,0.8,1.2,0.8,0.6),
                align = "h", axis="bt", ncol=5, labels=c("A", "B", "C", "D"))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

b2 <- plot_grid(INF_cls_p+NoLegend(), 
                INF_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("E", "G"))

b3 <- plot_grid(ENG_cls_p+NoLegend(), ENG_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("F", "H"))

l23 <- plot_grid(b2, b3, legend_p, ncol=3, rel_widths = c(1,1,0.3))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
plot_grid(l1, l23, ncol=1, rel_heights = c(1,2))
ggsave("ADV_Fig05_INF_ENG_CD14Mono.png", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_CD14Mono.pdf", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_CD14Mono.jpg", width = 12, height = 8, dpi=600)
get_test_table(DF, "CD14Mono")

#5. mDC
mDC <- add_ENGscore(mDC, total_genes)
mDC <- add_INFscore(mDC, INFgene_df)
mDC$GROUP <- factor(mDC$GROUP, levels=c("Synovial fluid", "A_acute", "Health"))
INF_dim <- FeaturePlot(mDC, features="inflammatory_score1", label=T)+labs(title="Inflammatory score")
ENG_dim <- FeaturePlot(mDC, features="energy_score1", label=T)+labs(title="Energy score")
DF <- mDC[[]]
INF_main_p <- boxplot_main_INF(DF)
ENG_main_p <- boxplot_main_ENG(DF)
INF_cls_p <- boxplot_cls_INF(DF)
ENG_cls_p <- boxplot_cls_ENG(DF)
### INF sub_group----
stat.test <- compare_means(
  inflammatory_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.30, .34, .30, .34, .34))
                        

INF_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="inflammatory_score1", 
            fill="GROUP"
  )+
  labs(y="inflammatory score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.1, .36)

### ENG sub_group----
stat.test <- compare_means(
  energy_score1 ~ GROUP, 
  data = DF, 
  group.by = "seurat_clusters",
  method = "wilcox.test") %>% 
  subset(p<0.05) %>%
  mutate(y.position = c(.12, .14, .12, .14))

ENG_clsgroup_p <- DF %>% 
  ggboxplot(x="GROUP", 
            y="energy_score1", 
            fill="GROUP"
  )+
  labs(y="energy score", x="")+
  theme_classic2()+
  facet_wrap(~seurat_clusters, scales = "free_x", ncol=6)+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length=0, vjust=0.4)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylim(-0.05,0.15)

group_le <- get_legend(INF_main_p+theme(legend.text.align = 0))
cls_le <- get_legend(INF_cls_p+theme(legend.text.align = 0))
legend_p <- plot_grid(cls_le, group_le, ncol=1)

l1 <- plot_grid(INF_dim, INF_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)), 
                ENG_dim, ENG_main_p+NoLegend()+scale_x_discrete(guide = guide_axis(n.dodge=2)),
                group_le, rel_widths = c(1.2,0.8,1.2,0.8,0.6),
                align = "h", axis="bt", ncol=5, labels=c("A", "B", "C", "D"))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

b2 <- plot_grid(INF_cls_p+NoLegend(), 
                INF_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("E", "G"))

b3 <- plot_grid(ENG_cls_p+NoLegend(), ENG_clsgroup_p+NoLegend(), 
                ncol=1, align = "hv", axis="btlr",
                labels=c("F", "H"))

l23 <- plot_grid(b2, b3, legend_p, ncol=3, rel_widths = c(1,1,0.3))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
plot_grid(l1, l23, ncol=1, rel_heights = c(1,2))
ggsave("ADV_Fig05_INF_ENG_mDC.png", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_mDC.pdf", width = 12, height = 8)
ggsave("ADV_Fig05_INF_ENG_mDC.jpg", width = 12, height = 8, dpi=600)
get_test_table(DF, "mDC")

# ----
# paper cells
# 1. InterMono
# 2. Treg
# 3. MAIT
# 4. CD14Mono
# 5. mDC

# supp cells
# 6. ?






