#Figure 5####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
library(Seurat)
library(sctransform)
library(ggplot2)
library(WGCNA)
library(tidyverse)
library(ggpubr)
library(ggsci)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))
load(paste0(data_directory, "CrossTissueAnalysis_Zenodo.Rdata"))
load(paste0(data_directory, "MxIHC_TMAdata_Zenodo.Rdata"))

#Figure 5 A-C####
Sample_UMAP <- 
  Merged_MetaData %>%
  filter(Group %in% c("Pancreas", "Oral", "Colon")) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = Sample.type)) +
  geom_point(size = 0.1) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_color_npg(name = "Sample type") +
  theme(legend.position = "right",legend.key.size = unit(10, "pt"))+
  guides(colour = guide_legend(override.aes = list(size = 2)))

Class_UMAP <- 
  Merged_MetaData %>%
  filter(Group %in% c("Pancreas", "Oral", "Colon")) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = predicted.id)) +
  geom_point(size = 0.1) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_colour_manual(values = Fibs_col.palette, name = "Predicted class") +
  theme(legend.position = "right",legend.key.size = unit(10, "pt")) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

Prob_VlnPlot <- 
  Merged_MetaData %>%
  filter(Group %in% c("Pancreas", "Oral", "Colon")) %>%
  ggplot(aes(x = predicted.id, y = prediction.score.max, fill = predicted.id)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  #geom_jitter(size = 0.1, width = 0.1, alpha = 0.1, show.legend = F) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_fill_manual(values = Fibs_col.palette, name = "Predicted class") +
  rotate_x_text(angle = 45) +
  theme(legend.position = "right", axis.title.x = element_blank(),
        legend.key.size = unit(10, "pt")) +
  ylab("Classification Probability") +
  ylim(c(0,1))

Fig_5ABC <- ggarrange(Sample_UMAP, Class_UMAP, Prob_VlnPlot, nrow = 3,
                   align = "v")
Fig_5ABC


SupplFig_5E <- 
  Merged_MetaData %>%
  filter(Group %in% c("Pancreas", "Oral", "Colon", "Lung")) %>%
  filter(predicted.id == "Myo") %>%
  ggplot(aes(x = Group, y = prediction.score.max, fill = predicted.id)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  #geom_jitter(size = 0.1, width = 0.1, alpha = 0.1, show.legend = F) +
  #facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_fill_manual(values = Fibs_col.palette, name = "Predicted class") +
  rotate_x_text(angle = 45) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        legend.key.size = unit(10, "pt")) +
  ylab("Classification Probability") +
  scale_x_discrete(labels = c("IPF", "PDAC", "HNSCC", "CRC")) +
  ylim(c(0,1))
SupplFig_5E


#MxIHC data####
#Figure 7EF####
All_TMA.data.df.Fibroblasts <- All_TMA.data.df %>%
  filter(Cell.type2 %in% c("Alveolar", "Adventitial", "Myo"))
table(All_TMA.data.df.Fibroblasts$Cell.type2)
All_TMA.data.df.Fibroblasts$Cell.type2 <- factor(
  as.character(All_TMA.data.df.Fibroblasts$Cell.type2),
  levels = c("Adventitial", "Alveolar",  "Myo")
)
dt <- as.table(as.matrix(table(All_TMA.data.df.Fibroblasts$Core_ID,
                               All_TMA.data.df.Fibroblasts$Cell.type2)))
Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
CoreData_long <- merge(MxIHC_TMA_metaData, Sample.pct_long,
                       by.x = "Core_ID", by.y = "Var1")
names(CoreData_long)[names(CoreData_long) == "Freq"] <- "Core.pct"
names(CoreData_long)[names(CoreData_long) == "Var2"] <- "Fibs_SubPop"

CoreData_long$Group <- factor(CoreData_long$Structure,
                              levels = unique(CoreData_long$Structure)[c(3,5,2,7,6,4,1)],
                              labels = c("Pancreas", "Oral", "Colon", "Lung", "Skin", "Breast", "Kidney"))
names(CoreData_long)
Fig_5E <- 
  CoreData_long[] %>%
  drop_na(Structure.filtered) %>%
  filter(Fibs_SubPop == "Adventitial") %>%
  filter(Structure %in% c("COLON", "PANCREAS", "HNSCC")) %>%
  ggplot(aes(x = TvN, y = Core.pct)) +
  theme_pubr(base_size = 7) +
  facet_wrap(~Group) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,125)) +
  stat_compare_means(comparisons = list(c("Normal", "Tumour")), size = 3,
                     #label = "p.signif", method = "wilcox",
                     label.y = 110, size = 2.5) +
  ylab("% of all fibroblast per core\n(MxIHC)") +
  theme(axis.title.x = element_blank(), legend.position = "none")
Fig_5E

Fig_5F <- 
  CoreData_long[] %>%
  drop_na(Structure.filtered) %>%
  filter(Fibs_SubPop == "Myo") %>%
  filter(Structure %in% c("COLON", "PANCREAS", "HNSCC")) %>%
  ggplot(aes(x = TvN, y = Core.pct)) +
  theme_pubr(base_size = 7) +
  facet_wrap(~Group) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,125)) +
  stat_compare_means(comparisons = list(c("Normal", "Tumour")), size = 3,
                     #label = "p.signif", method = "wilcox",
                     label.y = 110, size = 2.5) +
  ylab("% of all fibroblast per core\n(MxIHC)") +
  theme(axis.title.x = element_blank(), legend.position = "none")
Fig_5F

#Figure 5D - Representative MxIHC images####
s = "PANCREAS"
PDAC_CTR_C01 <- 
  All_TMA.data.df %>%
  filter(TvN == "Normal" & Core == "C01") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")
PDAC_Tumour_CO2 <- 
  All_TMA.data.df %>%
  filter(TvN == "Tumour" & Core == "C02") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")

HNSCC_CTR_E07 <- 
  All_TMA.data.df %>%
  filter(TvN == "Normal" & Core == "E07") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")

HNSCC_Tumour_E07 <- 
  All_TMA.data.df %>%
  filter(TvN == "Tumour" & Core == "E07") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")

COLON_CTR_B03 <- 
  All_TMA.data.df %>%
  filter(TvN == "Normal" & Core == "B03") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")

COLON_Tumour_B05 <- 
  All_TMA.data.df %>%
  filter(TvN == "Tumour" & Core == "B05") %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl.,
             colour = Cell.type2)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_colour_manual(values = Fibs_col.palette,
                      na.value = "grey80")

Fig_7D <- ggarrange(PDAC_CTR_C01, PDAC_Tumour_CO2,
                 HNSCC_CTR_E07, HNSCC_Tumour_E07,
                 COLON_CTR_B03, COLON_Tumour_B05,
                 nrow = 6, ncol = 1)
Fig_7D

#Supplementary Figure 5 B-D####
Sample_UMAP.IPF <- 
  Merged_MetaData %>%
  filter(Group %in% c("Lung")) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = Sample.type)) +
  geom_point(size = 0.1) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_color_npg(name = "Sample type") +
  theme(legend.position = "right",
        legend.key.size = unit(10, "pt"),
        legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 2), title.position = "top"))

Class_UMAP.IPF <- 
  Merged_MetaData %>%
  filter(Group %in% c("Lung")) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = predicted.id)) +
  geom_point(size = 0.1) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_colour_manual(values = Fibs_col.palette, name = "Predicted class") +
  theme(legend.position = "right",
        legend.key.size = unit(10, "pt"),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2), title.position = "top"))

Prob_VlnPlot.IPF <- 
  Merged_MetaData %>%
  filter(Group %in% c("Lung")) %>%
  ggplot(aes(x = predicted.id, y = prediction.score.max, fill = predicted.id)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  #geom_jitter(size = 0.1, width = 0.1, alpha = 0.1, show.legend = F) +
  facet_wrap(~Group, scales = "free", nrow = 1) +
  theme_pubr(base_size = 7) +
  scale_fill_manual(values = Fibs_col.palette, name = "Predicted class") +
  rotate_x_text(angle = 45) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        legend.key.size = unit(10, "pt")) +
  ylab("Classification Probability") +
  guides(fill = guide_legend(title.position = "top"))+
  ylim(c(0,1))

SupplFig_5BCD <- cowplot::plot_grid(Sample_UMAP.IPF, Class_UMAP.IPF, Prob_VlnPlot.IPF,
                             nrow = 3,
                             align = "v", axis = "l")
SupplFig_5BCD
