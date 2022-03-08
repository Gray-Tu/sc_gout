# Fig2 combine (PBMC Acute vs Health) and (Acute PBMC vs SF)
# ALT vis
setwd("With_SYN")

library(Seurat)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

# load data ----
load(file="integrated_rpca_syn_PBMCAcuteHealth.Rda")
# pbmc.syn.combined
load(file="addbpe_pc50m2h5_InfScore_rpca_Gout_combinePBMC.Rda")
# Gout.combined.50
Acute_HC_set <- subset(Gout.combined.50, GROUP %in% c("A_acute", "Health"))

# Fix Typo
"Megakaryocytes"
my_DF <- Gout.combined.50[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="Megakaryoctes",
                            "Megakaryocytes",
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
                                              "Megakaryocytes",
                                              "Progenitor cells"))

my_DF <- pbmc.syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="Megakaryoctes",
                            "Megakaryocytes",
                            cell_type))
pbmc.syn.combined$cell_type <- my_DF$cell_type
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
                                               "Megakaryocytes",
                                               "Progenitor cells"))


T_cell <- c(
  "Naive CD4 T cells",
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
  "Vd2 gd T cells"
)

NK_cell <- c("Natural killer cells")
MonoMa <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "Non classical monocytes",
  "Macrophages",
  "Macrophages M1",
  "INF-Macro",
  "lung Macro"
)
DC <- c(
  "Myeloid dendritic cells",
  "Plasmacytoid dendritic cells"
)
B_cell <- c(
  "Naive B cells",
  "Non-switched memory B cells",
  "Switched memory B cells",
  "Plasmablasts",
  "Exhausted B cells"
)
OTHER <- c(
  "Low-density basophils",
  "Low-density neutrophils",
  "Megakaryocytes",
  "Progenitor cells"
)

df <- pbmc.syn.combined[[]]
df$GROUP <- factor(df$GROUP, 
                   levels=c("Synovial fluid", "A_acute", "Health"))

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

get_MainTypeAbun_DF <- function(df){ # subset B_cell df
  
  DF <- df %>%
    group_by(CASE, GROUP, cell_type) %>%
    summarise(cell_number = n()) %>%
    mutate(total_number = sum(cell_number)) %>%
    mutate(percent = cell_number/total_number*100)
  
  # add zero if no cell
  total_cell <- unique(DF$cell_type)
  total_case <- c("MH", "RG", "N01F", "N05S", "YB", 
                  "Case04", "Case05", "Case06", "Case01",
                  "Case02", "Case03", "Case07", 
                  "Gout146", "Gout151", "Gout183")
  addzero_df <- data.frame()
  for (C in total_case){
    case_df <- subset(DF, CASE==C)
    G <- subset(case_to_group, CASE==C)$GROUP
    for (cell in total_cell){
      print(paste(G, C, cell))
      if (length(row.names(subset(case_df, cell_type==cell))) == 0){
        cell_number = 0
        total_number = 0
        percent = 0.0
        GROUP = G
        CASE = C
        cell_type=cell
        temp_df <- data.frame(CASE, GROUP, cell_type, cell_number, total_number, percent)
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

  return(outDF)
}

# Test base B ----
B_DF <- df %>% 
  subset(cell_type%in%B_cell)%>%
  get_MainTypeAbun_DF()

wilcox_B <- compare_means(percent~GROUP, 
                          B_DF, 
                          method = "wilcox.test", 
                          paired = FALSE,
                          group.by = "cell_type")
subset(wilcox_B, p < 0.05)


# Test base T and NK ----
T_NK_DF <- df %>% 
  subset(cell_type %in% union(T_cell, NK_cell))%>%
  get_MainTypeAbun_DF()

wilcox_T_NK <- compare_means(percent~GROUP, 
                             T_NK_DF, 
                          method = "wilcox.test", 
                          paired = FALSE,
                          group.by = "cell_type")
subset(wilcox_T_NK, p < 0.05)

# Test base MonoMaDC ----
MonoMaDC_DF <- df %>% 
  subset(cell_type %in% union(MonoMa, DC))%>%
  get_MainTypeAbun_DF()

wilcox_MonoMaDC <- compare_means(percent~GROUP, 
                                 MonoMaDC_DF, 
                             method = "wilcox.test", 
                             paired = FALSE,
                             group.by = "cell_type")
subset(wilcox_MonoMaDC, p < 0.05)

# Test Other base All ---
OTHER_DF_all <- df %>% 
  get_MainTypeAbun_DF() %>%
  subset(cell_type %in% OTHER)
wilcox_OTHER <- compare_means(percent~GROUP, 
                              OTHER_DF_all, 
                              method = "wilcox.test", 
                              paired = FALSE,
                              group.by = "cell_type")
subset(wilcox_OTHER, p < 0.05)
  

wilcox_B$base <- "B"
wilcox_MonoMaDC$base <- "MonoMaDC"
wilcox_T_NK$base <- "T_NK"
wilcox_OTHER$base <- "all"
wilcox_tab <- rbind(wilcox_B, wilcox_MonoMaDC, wilcox_T_NK, wilcox_OTHER)
write.csv(wilcox_tab, file="ADV_Fig02_wilcox_table.csv")

select_res <- subset(wilcox_tab, p<0.05)

#1 get Acute vs Health PBMC
# Classical monocytes
AH_CD14Mono_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                 cell_type=="Classical monocytes" & GROUP %in% c("A_acute", "Health"))
p <- ggboxplot(AH_CD14Mono_baseMonoMaDC,           
               x="GROUP", 
               y="percent",
               fill="GROUP")+
  ylab("% of Mono. & Ma. & DC")+
  theme_classic()+
  facet_wrap(~cell_type, scales = "free")+
  scale_fill_manual(values = c("#00BA38", "#619CFF"))
  
AH_boxplot_CD14Mono_baseMonoMaDC <- p +  stat_compare_means(comparisons = list(c("A_acute", "Health")), 
                                                       label = "p.signif",
                                                       tip.length=0,
                                                       vjust=0.4,
                                                       method = "wilcox.test")
# Intermediate monocytes
AH_InterMono_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                    cell_type=="Intermediate monocytes" & GROUP %in% c("A_acute", "Health"))
p <- ggboxplot(AH_InterMono_baseMonoMaDC,           
               x="GROUP", 
               y="percent",
               fill="GROUP")+
  ylab("% of Mono. & Ma. & DC")+
  theme_classic()+
  facet_wrap(~cell_type, scales = "free")+
  scale_fill_manual(values = c("#00BA38", "#619CFF"))
AH_boxplot_InterMono_baseMonoMaDC <- p +  stat_compare_means(comparisons = list(c("A_acute", "Health")), 
                                                        label = "p.signif",
                                                        tip.length=0,
                                                        vjust=0.4,
                                                        method = "wilcox.test")
# Plasmacytoid dendritic cells
AH_pDC_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                    cell_type=="Plasmacytoid dendritic cells" & GROUP %in% c("A_acute", "Health"))
p <- ggboxplot(AH_pDC_baseMonoMaDC,           
               x="GROUP", 
               y="percent",
               fill="GROUP")+
  ylab("% of Mono. & Ma. & DC")+
  theme_classic()+
  facet_wrap(~cell_type, scales = "free")+
  scale_fill_manual(values = c("#00BA38", "#619CFF"))
AH_boxplot_pDC_baseMonoMaDC <- p +  stat_compare_means(comparisons = list(c("A_acute", "Health")), 
                                                             label = "p.signif",
                                                             tip.length=0,
                                                             vjust=0.4,
                                                             method = "wilcox.test")

# Naive CD4 T cells
AH_NavieCD4T_baseT_NK <- subset(T_NK_DF, 
                             cell_type=="Naive CD4 T cells" & GROUP %in% c("A_acute", "Health")) 
p <- ggboxplot(AH_NavieCD4T_baseT_NK,           
               x="GROUP", 
               y="percent",
               fill="GROUP")+
  ylab("% of T & NK cells")+
  theme_classic()+
  facet_wrap(~cell_type, scales = "free")+
  scale_fill_manual(values = c("#00BA38", "#619CFF"))
AH_boxplot_NavieCD4T_baseT_NK <- p +  stat_compare_means(comparisons = list(c("A_acute", "Health")), 
                                                      label = "p.signif",
                                                      tip.length=0,
                                                      vjust=0.4,
                                                      method = "wilcox.test") 


#2 get SF vs Acute
make_ADV_Fig2_boxplot <- function(DF, ylab_str, comparisons){
  p <- ggboxplot(DF,           
                 x="GROUP", 
                 y="percent",
                 fill="GROUP")+
    ylab(ylab_str)+
    theme_classic()+
    facet_wrap(~cell_type, scales = "free")+
    scale_fill_manual(values = c("#F8766D", "#00BA38"))
  
  DF_boxplot <- p +  stat_compare_means(comparisons = comparisons,
                                        label = "p.signif",
                                        tip.length=0,
                                        vjust=0.4,
                                        method = "wilcox.test") 
  return(DF_boxplot)
} 

# Plasmablasts
SA_plasmablaste_baseB <- subset(B_DF, 
                                cell_type=="Plasmablasts" & GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_plasmablaste_baseB <- make_ADV_Fig2_boxplot(SA_plasmablaste_baseB, "% of B cells", list(c("A_acute", "Synovial fluid")))

# Classical monocytes
SA_CD14Mono_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                cell_type=="Classical monocytes" & GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_CD14Mono_baseMonoMaDC <- make_ADV_Fig2_boxplot(SA_CD14Mono_baseMonoMaDC, 
                                                          "% of Mono. & Ma. & DC",
                                                          list(c("A_acute", "Synovial fluid")))

# Non classical monocytes
SA_CD16Mono_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                   cell_type=="Non classical monocytes" & 
                                     GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_CD16Mono_baseMonoMaDC <- make_ADV_Fig2_boxplot(SA_CD16Mono_baseMonoMaDC, 
                                                          "% of Mono. & Ma. & DC",
                                                          list(c("A_acute", "Synovial fluid")))

# Macrophages
SA_Ma_baseMonoMaDC <- subset(MonoMaDC_DF, 
                                   cell_type=="Macrophages" & 
                                     GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_Ma_baseMonoMaDC <- make_ADV_Fig2_boxplot(SA_Ma_baseMonoMaDC, 
                                                          "% of Mono. & Ma. & DC",
                                                          list(c("A_acute", "Synovial fluid")))

# Naive CD4 T cells
SA_NavieCD4T_baseT_NK <- subset(T_NK_DF, 
                             cell_type=="Naive CD4 T cells" & 
                               GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_NavieCD4T_baseT_NK <- make_ADV_Fig2_boxplot(SA_NavieCD4T_baseT_NK, 
                                                       "% of T & NK cells",
                                                    list(c("A_acute", "Synovial fluid")))

# Th1/Th17 cells
SA_Th1Th17_baseT_NK <- subset(T_NK_DF, 
                                cell_type=="Th1/Th17 cells" & 
                                  GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_Th1Th17_baseT_NK <- make_ADV_Fig2_boxplot(SA_Th1Th17_baseT_NK, 
                                                       "% of T & NK cells",
                                                       list(c("A_acute", "Synovial fluid")))

# Naive CD8 T cells
SA_NavieCD8T_baseT_NK <- subset(T_NK_DF, 
                                cell_type=="Naive CD8 T cells" & 
                                  GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_NavieCD8T_baseT_NK <- make_ADV_Fig2_boxplot(SA_NavieCD8T_baseT_NK, 
                                                       "% of T & NK cells",
                                                       list(c("A_acute", "Synovial fluid")))

# Central memory CD8 T cells
SA_CD8Tcm_baseT_NK <- subset(T_NK_DF, 
                                cell_type=="Central memory CD8 T cells" & 
                                  GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_CD8Tcm_baseT_NK <- make_ADV_Fig2_boxplot(SA_CD8Tcm_baseT_NK, 
                                                       "% of T & NK cells",
                                                       list(c("A_acute", "Synovial fluid")))

# Effector memory CD8 T cells
SA_CD8Tem_baseT_NK <- subset(T_NK_DF, 
                             cell_type=="Effector memory CD8 T cells" & 
                               GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_CD8Tem_baseT_NK <- make_ADV_Fig2_boxplot(SA_CD8Tem_baseT_NK, 
                                                    "% of T & NK cells",
                                                    list(c("A_acute", "Synovial fluid")))

# MAIT cells
SA_MAIT_baseT_NK <- subset(T_NK_DF, 
                             cell_type=="MAIT cells" & 
                               GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_MAIT_baseT_NK <- make_ADV_Fig2_boxplot(SA_MAIT_baseT_NK, 
                                                    "% of T & NK cells",
                                                    list(c("A_acute", "Synovial fluid")))

# Natural killer cells
SA_NK_baseT_NK <- subset(T_NK_DF, 
                           cell_type=="Natural killer cells" & 
                             GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_NK_baseT_NK <- make_ADV_Fig2_boxplot(SA_NK_baseT_NK, 
                                                  "% of T & NK cells",
                                                  list(c("A_acute", "Synovial fluid")))

# Low-density basophils
SA_LDbasophils_base_all <- subset(OTHER_DF_all, 
                         cell_type=="Low-density basophils" & 
                           GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_LDbasophils_base_all  <- make_ADV_Fig2_boxplot(SA_LDbasophils_base_all, 
                                                "% of all cells",
                                                list(c("A_acute", "Synovial fluid")))

# Megakaryoctes
SA_Megakaryoctes_base_all <- subset(OTHER_DF_all, 
                                  cell_type=="Megakaryoctes" & 
                                    GROUP %in% c("A_acute", "Synovial fluid")) 
SA_boxplot_Megakaryoctes_base_all  <- make_ADV_Fig2_boxplot(SA_Megakaryoctes_base_all, 
                                                          "% of all cells",
                                                          list(c("A_acute", "Synovial fluid")))


SA_boxplot_plasmablaste_baseB
SA_boxplot_CD14Mono_baseMonoMaDC
SA_boxplot_CD16Mono_baseMonoMaDC
SA_boxplot_Ma_baseMonoMaDC
SA_boxplot_NavieCD4T_baseT_NK
SA_boxplot_Th1Th17_baseT_NK
SA_boxplot_NavieCD8T_baseT_NK
SA_boxplot_CD8Tcm_baseT_NK
SA_boxplot_CD8Tem_baseT_NK
SA_boxplot_MAIT_baseT_NK
SA_boxplot_NK_baseT_NK
SA_boxplot_LDbasophils_base_all
SA_boxplot_Megakaryoctes_base_all

# get legend
legend <- get_legend(
  ggboxplot(T_NK_DF,           
            x="GROUP", 
            y="percent",
            fill="GROUP")+
    theme_classic2()#+
    #geom_point(aes(colour = CASE))+
    #guides(col = guide_legend(ncol = 2))#, title.position = "right"))
)
p_SA_13cells <- plot_grid(
  SA_boxplot_plasmablaste_baseB + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD14Mono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD16Mono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  SA_boxplot_Ma_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  SA_boxplot_NavieCD4T_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_Th1Th17_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_NavieCD8T_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD8Tcm_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD8Tem_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_MAIT_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_NK_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_LDbasophils_base_all + theme(legend.position="none") + labs(x=""),
  SA_boxplot_Megakaryoctes_base_all + theme(legend.position="none") + labs(x=""),
  legend,
  ncol=4, 
  align = 'vh',
  labels = c("E", "F", "G", "H", "I", "J",
             "K", "L", "M", "N", "O", "P", "Q"),
  hjust = -1
) +  theme(plot.background = element_rect(fill = "white", colour = "white"))

p_SA_13cells #+ prow + legend + plot_layout(widths = c(2, 0.3))
plot_grid(p_AH_4cells, p_SA_13cells, rel_heights = c(1,3), ncol=1)

# show
# AvsH
# CD14Mono, InterMono, pDC, NaiveCD4T
# SF vs A
# CD14Mono, Macrophage, Th1_Th17, CD8Tcm,
# CD8Tem, MAIT, NK

ADV_fig02 <- plot_grid(
  AH_boxplot_CD14Mono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  AH_boxplot_InterMono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  AH_boxplot_pDC_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  AH_boxplot_NavieCD4T_baseT_NK + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_plasmablaste_baseB + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD14Mono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_CD16Mono_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  SA_boxplot_Ma_baseMonoMaDC + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_NavieCD4T_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_Th1Th17_baseT_NK + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_NavieCD8T_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD8Tcm_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_CD8Tem_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_MAIT_baseT_NK + theme(legend.position="none") + labs(x=""),
  SA_boxplot_NK_baseT_NK + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_LDbasophils_base_all + theme(legend.position="none") + labs(x=""),
  #SA_boxplot_Megakaryoctes_base_all + theme(legend.position="none") + labs(x=""),
  legend,
  ncol=4, 
  align = 'vh',
  axis = "lr",
  labels = c("A", "B", "C", "D", "E",
             "F", "G", "H", "I", "J", "K"),
             #"K", "L", "M", "N", "O", "P", "Q"),
  hjust = -1
) +  theme(plot.background = element_rect(fill = "white", colour = "white"))

ADV_fig02
ggsave("ADValt_Fig02.pdf",  width = 12, height = 10)
ggsave("ADValt_Fig02.png",  width = 12, height = 10)
