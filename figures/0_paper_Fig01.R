#paper plot

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
load(file="pc50_syn.rda")
# syn.combined


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

my_DF <- syn.combined[[]] %>% 
  mutate_all(as.character) %>% 
  mutate(cell_type = ifelse(cell_type=="Megakaryoctes",
                            "Megakaryocytes",
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
                                               "Megakaryocytes",
                                               "Progenitor cells"))

Acute_HC_set <- subset(Gout.combined.50, GROUP %in% c("A_acute", "Health"))


# fix color order ----
cell_to_cols <- read.csv("mod01_col.txt", sep="\t")
rownames(cell_to_cols) <- cell_to_cols$cell

# case to day
case_to_day <- as.data.frame(list(
  "CASE" = c("MH", "RG", "N01F", "N05S", "YB", 
             "Case01", "Case02", "Case03", "Case04",
             "Case05", "Case06", "Case07",
             "Gout146", "Gout151", "Gout183"),
  "PainDay" = c(0, 0, 0, 0, 0,
                3, 3, 14, 1, 
                1, 2, 19,
                0, 0, 0)
))
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

rownames(case_to_day) <- case_to_day$CASE
d <- case_to_day$PainDay
names(d) <- case_to_day$CASE
pbmc.syn.combined$pain_day <- d[match(pbmc.syn.combined$CASE, names(d))]

pbmc.syn.combined$CASE <- factor(pbmc.syn.combined$CASE,
                                 levels=c("MH", "RG", "N01F", "N05S", "YB", 
                                          "Case04", "Case05", "Case06", "Case01",
                                          "Case02", "Case03", "Case07",
                                          "Gout146", "Gout151", "Gout183"
                                          ))
pbmc.syn.combined$GROUP <- factor(pbmc.syn.combined$GROUP,
                                  levels=c("Synovial fluid",
                                           "A_acute",
                                           "Health"))
case_to_day <- pbmc.syn.combined[[]]%>%
  dplyr::select(CASE, GROUP, pain_day) %>%
  unique()
         
# Fig01----
# show dimplot of PBMC 
# paper usage Fig01, abundance----
levels_select_1 = levels(Acute_HC_set$cell_type) %in% unique(Acute_HC_set$cell_type)
order_select_2 = levels(Acute_HC_set$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]
# PBMC Acute-Health
b_label <- c("white", "white", "black", "black",
             "white", "black", "white", "white",
             "black", "black", "black", "white",
             "white", "white", "white", "black",
             "black", "black", "black", "white",
             "white", "black", "white", "white",
             "black", "black", "black", "white",
             "white", "black", "black")
Fig01_b <- DimPlot(Acute_HC_set, group.by = "cell_type", 
              label=T, repel=T, label.box=T, label.color=b_label,
              #label.color="white",
              #cols=colors_select_3$cols
              cols=colors_select_3$mod1cols
)+labs(title="PBMC")#, tag="B")
#unique(ggplot_build(Fig01_b)$data[[1]]$colour)

levels_select_1 = levels(syn.combined$cell_type) %in% unique(syn.combined$cell_type)
order_select_2 = levels(syn.combined$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]
# SF threes samples

c_label <- c("black", "black", "white", "white",
             "black", "white", "white", "black",
             "white", "black", "white", "black",
             "black", "white", "black", "white",
             "black", "black", "white", "black",
             "black", "black", "white", "black",
             "white", "black", "black", "white")
Fig01_c <- DimPlot(syn.combined, 
              group.by = "cell_type", 
              label=T, repel=T, label.box=T, label.color = c_label,
              #label.color="white",
              #cols=colors_select_3$cols
              cols=colors_select_3$mod1cols
)+labs(title="SF")#, tag="D")
#Fig01_c

library(magick)
library(ggplotify)
library(grImport2)
Fig01_a <- image_read_svg('Fig01_A_Design_View.svg')# magick
Fig01_a <- image_crop(Fig01_a, "3631x2380+10+10")
Fig01_a <- as.grob(Fig01_a) 

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



get_MainTypeAbun_DF <- function(df){ 
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
  outDF$cell_type <- factor(outDF$cell_type,
                            levels=levels(DF$cell_type))
  
  return(outDF)
}

df <- pbmc.syn.combined[[]]
ALL_DF <- df %>% get_MainTypeAbun_DF()

pain_day_bar <- case_to_day %>% 
  ggplot(aes(x=CASE, y=pain_day))+
  geom_bar(stat="identity")+
  facet_grid(~GROUP, scales = "free", space="free_x")+
  theme_classic2()+
  labs(x="")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

B_baseALL_DF <- subset(ALL_DF, cell_type %in% B_cell)
T_NK_baseALL_DF <- subset(ALL_DF, cell_type %in% union(T_cell, NK_cell))
MonoDCMa_baseALL_DF <- subset(ALL_DF, cell_type %in% union(MonoMa, DC))
OTHER_baseALL_DF <- subset(ALL_DF, cell_type %in% OTHER)

subset_B <- subset(df, cell_type %in% B_cell)
subset_T_NK <- subset(df, cell_type %in% union(T_cell, NK_cell))
subset_MonoDCMa <- subset(df, cell_type %in% union(MonoMa, DC))
subset_OTHER <- subset(df, cell_type %in% OTHER)

# B ----
B_DF <- subset_B %>% get_MainTypeAbun_DF()
levels_select_1 = levels(B_DF$cell_type) %in% unique(B_DF$cell_type)
order_select_2 = levels(B_DF$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]

B_bar_xCASE <- ggplot(B_DF, 
                      aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme_classic2()+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='B cells', y="% of B cells", x="") 

B_bar_xCASE_all <- ggplot(B_baseALL_DF, 
       aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='B cells', y="% of all cells", x="") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+NoLegend()

B_bar_xCASE_2 <- B_bar_xCASE_all / B_bar_xCASE + plot_layout(guides = 'collect')

# T NK----
T_NK_DF <- subset_T_NK %>% get_MainTypeAbun_DF()
levels_select_1 = levels(T_NK_DF$cell_type) %in% unique(T_NK_DF$cell_type)
order_select_2 = levels(T_NK_DF$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]
T_NK_bar_xCASE <- ggplot(T_NK_DF, 
                      aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme_classic2()+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='T & NK cells', y="% of T & NK cells", x="") 

T_NK_bar_xCASE_all <- ggplot(T_NK_baseALL_DF, 
                          aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='T & NK cells', y="% of all cells", x="") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+NoLegend()

T_NK_bar_xCASE_2 <- T_NK_bar_xCASE_all / T_NK_bar_xCASE + plot_layout(guides = 'collect')

# MonoDCMa ----
MonoDCMa_DF <- subset_MonoDCMa %>% get_MainTypeAbun_DF()
levels_select_1 = levels(MonoDCMa_DF$cell_type) %in% unique(MonoDCMa_DF$cell_type)
order_select_2 = levels(MonoDCMa_DF$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]

MonoDCMa_bar_xCASE <- ggplot(MonoDCMa_DF, 
                         aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme_classic2()+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='Mono. & Ma. & DC', y="% of Mono. & Ma. & DC", x="") 

MonoDCMa_bar_xCASE_all <- ggplot(MonoDCMa_baseALL_DF,
                             aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+
  scale_fill_manual(values = colors_select_3$mod1cols)+
  labs(fill='Mono. & Ma. & DC', y="% of all cells", x="") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+NoLegend()

MonoDCMa_bar_xCASE_2 <- plot_grid(MonoDCMa_bar_xCASE_all,
                              MonoDCMa_bar_xCASE+NoLegend(),
                              ncol=1, align="v", axis="lr")


OTHER_DF <- subset_OTHER %>% get_MainTypeAbun_DF()
levels_select_1 = levels(OTHER_DF$cell_type) %in% unique(OTHER_DF$cell_type)
order_select_2 = levels(OTHER_DF$cell_type)[levels_select_1]
colors_select_3 = cell_to_cols[order_select_2,]

OTHER_bar_xCASE <- ggplot(OTHER_DF, 
                             aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+ 
  scale_fill_manual(values = colors_select_3$mod1cols)+ 
  labs(fill='Cells', y="% of other cells", x="")

OTHER_bar_xCASE_all <- ggplot(OTHER_baseALL_DF,
                                 aes(x = CASE, y = percent, fill = cell_type))+
  geom_bar(stat = "identity")+
  facet_grid(~GROUP, scales="free", space="free_x")+
  theme_classic2()+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.justification = "left")+
  #scale_fill_manual(values = colors_select_3$cols)+
  scale_fill_manual(values = colors_select_3$mod1cols)+
  labs(fill="Cells",
    y="% of all cells", x="") #+
  #theme(axis.ticks.x = element_blank(),
  #      axis.text.x = element_blank())+NoLegend()

OTHER_bar_xCASE_2 <- plot_grid(OTHER_bar_xCASE_all,
                               OTHER_bar_xCASE+NoLegend(),
                               ncol=1, align="v", axis="lr")


d_B <- B_bar_xCASE + theme(axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())
d_T_NK <- T_NK_bar_xCASE + theme(axis.ticks.x = element_blank(),
                                 axis.text.x = element_blank())
d_MonoDCMa <- MonoDCMa_bar_xCASE + theme(axis.ticks.x = element_blank(),
                                         axis.text.x = element_blank())
Fig01_d <- plot_grid(
  pain_day_bar, 
  d_B, 
  d_T_NK,
  #MonoDCMa_bar_xCASE,
  d_MonoDCMa,
  OTHER_bar_xCASE_all,
  rel_heights= c(0.25, 1, 1, 1, 0.4),
  align="v",
  axis="lr",
  labels=c("D", "E", "F", "G", "H"),
  ncol=1
)

ggsave("Fig01.pdf", plot = plot_grid(plot_grid(Fig01_a, Fig01_b, Fig01_c, ncol=1, labels=c("A", "B", "C")),
                                     Fig01_d, labels=c("", "")),
       width = 24, height = 18)
ggsave("Fig01.png", plot = plot_grid(plot_grid(Fig01_a, Fig01_b, Fig01_c, ncol=1, labels=c("A", "B", "C"))+theme(plot.background = element_rect(fill = "white", colour = "white")),
                                     Fig01_d, labels=c("", "")),
       width = 24, height = 18)

# table of DF
B_DF$base <- "B cells"
T_NK_DF$base <- "T & NK cells"
MonoDCMa_DF$base <- "Mono. & Ma. & DC"
OTHER_baseALL_DF$base <- "all cells"

write.csv(
  rbind(B_DF, T_NK_DF, MonoDCMa_DF, OTHER_baseALL_DF),
  "Table_Fig01abundance.csv")

# Fig01_D_B
# pain_day_bar / B_bar_xCASE + plot_layout(heights =  c(0.2, 1))
# ggsave("Fig01_D_Bcell.png", width = 12, height = 6)
# 
# # Fig01_D_T_NK
# pain_day_bar / T_NK_bar_xCASE + plot_layout(heights =  c(0.2, 1))
# ggsave("Fig01_D_T_NKcell.png", width = 12, height = 6)
# 
# # Fig01_D_MonoMaDC
# pain_day_bar / MonoDCMa_bar_xCASE + plot_layout(heights =  c(0.2, 1))
# ggsave("Fig01_D_MonoDCMacell.png", width = 12, height = 6)
# 


