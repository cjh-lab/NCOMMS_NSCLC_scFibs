#Figure 1####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
library(Seurat)
library(WGCNA)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "TLDS_AllCells_RPCAintegrated.Rdata"))
load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))
load(paste0(data_directory, "HCLA.SS_muralSig_data.Rdata"))
load(paste0(data_directory, "Qian_AllCells_RPCAintegrated.Rdata"))


#Figure 1B####
#Gender Data redacted for GDPR concerns
Stage_barplot <- 
  Sample_MetaData %>%
  filter(Dataset == "TLDS") %>%
  filter(!Sample.Subtype == "Control") %>%
  ggplot(aes(y = Dataset, fill = Stage)) +
  geom_bar(position = "fill") +
  facet_wrap(~Sample.Subtype) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", legend.key.size = unit(3, "pt"), legend.direction = "horizontal",
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  guides(fill = guide_legend(title.position = "top")) +
  xlab("Sample Fraction") +
  scale_fill_brewer(palette = "Set2", name = "Stage")
Smoking_barplot <- 
  Sample_MetaData %>%
  filter(Dataset == "TLDS") %>%
  filter(!Sample.Subtype == "Control") %>%
  drop_na(Smoking.status) %>%
  ggplot(aes(y = Dataset, fill = Smoking.status)) +
  geom_bar(position = "fill") +
  facet_wrap(~Sample.Subtype) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", legend.key.size = unit(3, "pt"), legend.direction = "horizontal",
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  guides(fill = guide_legend(title.position = "top")) +
  xlab("Sample Fraction") +
  scale_fill_brewer(palette = "Set2", name = "Smoking Status")
Grade_barplot <- 
  Sample_MetaData %>%
  filter(Dataset == "TLDS") %>%
  filter(!Sample.Subtype == "Control") %>%
  drop_na(Differentiation) %>%
  ggplot(aes(y = Dataset, fill = Differentiation)) +
  geom_bar(position = "fill") +
  facet_wrap(~Sample.Subtype) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", legend.key.size = unit(3, "pt"), legend.direction = "horizontal",
        axis.title.y = element_blank(), axis.text.y = element_blank()) +
  guides(fill = guide_legend(title.position = "top")) +
  xlab("Sample Fraction") +
  scale_fill_brewer(palette = "Set2", name = "Grade")

Fig_1B <- ggarrange(NULL,
  ggarrange(Smoking_barplot + theme(legend.justification = 0,
                                    axis.text.x = element_blank(), axis.title.x = element_blank(),
                                    axis.ticks.y = element_blank()),
           Stage_barplot+ theme(legend.justification = 0,
                               axis.text.x = element_blank(), axis.title.x = element_blank(),
                               axis.ticks.y = element_blank(),strip.background = element_blank(),
                               strip.text.x = element_blank()),
          Grade_barplot+ theme(legend.justification = 0,
                               axis.ticks.y = element_blank(),strip.background = element_blank(),
                               strip.text.x = element_blank()) +
            rotate_x_text(angle = 45),
          ncol = 1, nrow = 3, align = "v", heights = c(1.3,1,1.7)),
  ncol = 2, widths = c(0.05,0.95))
Fig_1B

#Supplementary Data 1####
TLDS_metaData <- Sample_MetaData %>%
  filter(Dataset == "TLDS")

cols2keep <- colnames(TLDS_metaData)[(!colSums(apply(TLDS_metaData,2,is.na)) == nrow(TLDS_metaData))]

Table_S1 <- TLDS_metaData[, cols2keep]
Control.samples <- as.character(Table_S1$PatientID)[Table_S1$Sample.type == "Control"]
Table_S1$Control.tissue.Sampled <- Table_S1$PatientID %in% Control.samples
Table_S1.final <- Table_S1[!Table_S1$Sample.type == "Control", ]
Table_S1.final <- Table_S1.final[, -c(20:21)]

#Mesenchymal cell processing####
TLDS.combined <- SetIdent(TLDS.combined, value = "Cell.type")

TLDS.combined$Cell.type3 <- factor(
  TLDS.combined$Cell.type2,
  levels = levels(TLDS.combined$Cell.type2),
  labels = c(levels(TLDS.combined$Cell.type2)[1],
             "Stroma",
             levels(TLDS.combined$Cell.type2)[3:8])
)
#Figure 1C####
Fig_1C <- 
  DimPlot(TLDS.combined, label = F, group.by = "Cell.type3")  &
  theme_pubr(base_size = 5) &
  #NoLegend() &
  theme(plot.title = element_blank(), 
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "right", legend.key.size = unit(1.5, "pt"),
        legend.justification = c(0,0.5)) &
  scale_color_brewer(palette = "Paired", 
                     labels = c("T cells", "Stroma", "Mono/Mac", "Epithelial", "B cells",
                                "Mast cells", "Endothelial", "AT1"), na.translate = F)
Fig_1C

table(TLDS.combined$orig.ident %in% Sample_MetaData$SampleID)
TLDS.combined$SampleID = TLDS.combined$orig.ident
Cell_MetaData <- TLDS.combined@meta.data[,"SampleID", drop = F]
matchSamples <- match(Cell_MetaData$SampleID, Sample_MetaData$SampleID)
table(is.na(matchSamples))
Cell_MetaData <- cbind(Cell_MetaData, Sample_MetaData[matchSamples, ])
TLDS.combined <- AddMetaData(TLDS.combined, Cell_MetaData)

#Supplementary Figure 1A####
#Calculate nearest neighbour z-scores on raw data
DefaultAssay(TLDS.combined) <- "RNA"
TLDS.combined <- NormalizeData(TLDS.combined)
TLDS.combined <- FindVariableFeatures(TLDS.combined)
TLDS.combined <- ScaleData(TLDS.combined, verbose = FALSE)
TLDS.combined <- RunPCA(TLDS.combined, npcs = 30, verbose = FALSE)
TLDS.combined <- FindNeighbors(TLDS.combined, reduction = "pca", dims = 1:30)
Batchtest_raw <- Batch.Quant(seurat_obj = TLDS.combined,
                             vars = c("PatientID"),
                             nPerm = 1000, assay_nn = "RNA")

#Calculate nn z-scores on integrated data
Batchtest_rpca <- Batch.Quant(seurat_obj = TLDS.combined,
                              vars = c("PatientID"),
                              nPerm = 1000, assay_nn = "integrated")

#Collate and plot results
Batch.test_res <- data.frame(
  raw = Batchtest_raw$PatientID,
  rpca  = Batchtest_rpca$PatientID
)

SupplFig_1A <- 
  reshape2::melt(Batch.test_res) %>%
  ggplot(aes(x = variable, y = value)) +
  theme_pubr(base_size = 5) +
  #geom_violin
  geom_boxplot(outlier.size = 0.1) +
  #scale_y_log10() + 
  geom_hline(yintercept = 1.96, linetype = "dashed") +
  rotate_x_text(angle = 45) +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("Raw Data", "rPCA Integrated", "Clustering")) +
  ylab("Intra-Patient nn overlap\n(Z-score)") +
  ggforce::facet_zoom(ylim = c(-1, 7))
SupplFig_1A

#Supplementary Figure 1B####
SupplFig_1B <- 
  TLDS.combined@meta.data %>%
  drop_na(Cell.type3) %>%
  ggplot(aes(x = orig.ident, fill = Cell.type3)) +
  theme_pubr(base_size = 5) +
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,"pt"),
        axis.text.x =  element_blank(), legend.position = "right") +
  facet_grid(~Sample.Subtype, scales = "free_x", space ="free_x") +
  geom_bar(position = "fill", alpha = 0.9) +
  scale_fill_brewer(palette = "Paired", labels = c("T cells", "Stroma", "Mono/Mac", "Epithelial", "B cells",
                                          "Mast cells", "Endothelial", "AT1")) +
  xlab("Sample") + ylab("Fraction") 
SupplFig_1B

#Supplementary Figure 1C&D####
SupplFig_1CD <-
    ggarrange(
  TLDS.combined@meta.data %>%
    drop_na(Cell.type3) %>%
    ggplot(aes(x = Cell.type3, y = nFeature_RNA)) +
    theme_pubr(base_size = 5) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    geom_violin(aes(fill = Cell.type3), alpha = 0.9, scale = "width") +
    geom_jitter(alpha = 0.1, size = 0.1, width = 0.2) +
    scale_y_log10() +
    scale_fill_brewer(palette = "Paired") +
    rotate_x_text(angle = 45) +
    ylab("nGenes"),
  TLDS.combined@meta.data %>%
    drop_na(Cell.type3) %>%
    ggplot(aes(x = Cell.type3, y = nCount_RNA)) +
    theme_pubr(base_size = 5) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    geom_violin(aes(fill = Cell.type3), alpha = 0.9, scale = "width") +
    geom_jitter(alpha = 0.1, size = 0.1, width = 0.2) +
    scale_y_log10() +
    scale_fill_npg() +
    rotate_x_text(angle = 45) +
    ylab("nCounts"),
  ncol = 2
)
SupplFig_1CD


#Supplementary Figure 1E####
#TLDS.combined.Markers <- FindAllMarkers(TLDS.combined, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
DefaultAssay(TLDS.combined) <- "RNA"
Markers <- c("TRBC2", "DCN", "LYZ", "KRT19", "CD79A", "KIT", "VWF", "AGER")
SupplFig_1E <- 
  FeaturePlot(TLDS.combined, features = Markers, ncol = 4, reduction = "umap") &
  theme_pubr(base_size = 5) &
  theme(plot.title = element_text(face = "italic"), legend.position = "none")
SupplFig_1E


#Figure 1D####
TLDS_Mesenchymal_seurat <- subset(TLDS.combined, idents = "Fibroblasts")
Fig_1D <- 
  FeaturePlot(TLDS_Mesenchymal_seurat, c("MCAM", "RGS5","ACTA2", "DPT"),
            pt.size = 0.1, reduction = "umap") &
  theme_pubr(base_size = 5) &
  ylim(c(2.5,10)) &
  theme(legend.position = "none",
        legend.justification = 0,
        legend.key.width = unit(2,"pt"),
        legend.key.height = unit(5, "pt"),
        legend.background = element_blank(),
        plot.title = element_text(face = "italic")) 
Fig_1D

#Separating Mural cells and fibroblasts
DefaultAssay(TLDS_Mesenchymal_seurat) <- "RNA"
TLDS_Mesenchymal_seurat <- AddModuleScore(
  TLDS_Mesenchymal_seurat, features = Consensus_FibsvMural_sigs
)
names(TLDS_Mesenchymal_seurat@meta.data)
TLDS_Mesenchymal_seurat$Fib_or_Mural <- "Undetermined"
TLDS_Mesenchymal_seurat$Fib_or_Mural[TLDS_Mesenchymal_seurat$Cluster1 > 0 & 
                                       TLDS_Mesenchymal_seurat$Cluster1 - TLDS_Mesenchymal_seurat$Cluster2 > 0.1] <- "Fibroblast"
TLDS_Mesenchymal_seurat$Fib_or_Mural[TLDS_Mesenchymal_seurat$Cluster2 > 0 & 
                                       TLDS_Mesenchymal_seurat$Cluster2 - TLDS_Mesenchymal_seurat$Cluster1 > 0.1] <- "Mural.cells"

TLDS_Mesenchymal_seurat$Fibrolast.Module.Score <- TLDS_Mesenchymal_seurat$Cluster1
TLDS_Mesenchymal_seurat$Mural.Module.Score <- TLDS_Mesenchymal_seurat$Cluster2

#Figure 1E-G: HCLA Mural cell score optimisation####
Fig_1E <-
  Mural.v.Fibs_markers %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), label = label, colour = p_val_adj < 0.01 & abs(avg_log2FC) > 0.5)) +
  geom_point(size = 0.1) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(size = 1.75, colour = "black", fontface = "italic", min.segment.length = 0, max.overlaps = 30) +
  scale_color_manual(values = c("grey70", "red")) + 
  ylab("-log10(adj.P)") +
  xlab("Mural vs Fibroblast (log2[FC])")
Fig_1E

HCLA.SS_muralSig_data$Fib_or_Mural[is.na(HCLA.SS_muralSig_data$Fib_or_Mural)] <- "Undetermined"
Fig_1F <- 
  HCLA.SS_muralSig_data %>%
  ggplot(aes(x = Mural.Cells.Module.Score,
             y = Fibroblast.Module.Score,
             colour = Fib_or_Mural)) +
  theme_pubr(base_size = 5) + 
  geom_point(size = 0.1) + 
  xlab("Mural Cells Module Score") + 
  ylab("Fibroblast Module Score") +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "top", plot.title = element_blank(),
        legend.key.height = unit(2, "pt"), legend.direction = "vertical",
        legend.justification = 0.25,
        legend.background = element_blank(), legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size = 2)))
Fig_1F

HCLA.SS_muralSig_data$free_annotation_merged2 <- factor(HCLA.SS_muralSig_data$free_annotation_merged,
                                                        labels = c("Mural cells", "Fibroblasts", "Other"))
Fig_1G <- 
  HCLA.SS_muralSig_data %>%
  ggplot(aes(x = free_annotation, fill = Fib_or_Mural)) +
  theme_pubr(base_size = 5) +
  geom_bar(position = "fill") +
  rotate_x_text(angle = 65) +
  theme(legend.key.size = unit(4,"pt")) +
  facet_grid(~free_annotation_merged2,  scale = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Dark2", name = "Cell Classification") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  ylab("Fraction") +
  xlab("Travaglini et al. Annotaion")
Fig_1G

#Figure 1H-J####
Fig_1HIJ <- 
  ggarrange(
  FeaturePlot(TLDS_Mesenchymal_seurat, c("Fibrolast.Module.Score"), pt.size = 0.1) &
    scale_color_viridis_c() &
    ylim(c(2.5,10)) &
    theme_pubr(base_size =7) &
    theme(legend.position = c(0, 0.9),
          legend.justification = 0,
          legend.key.width = unit(2,"pt"),
          legend.key.height = unit(5, "pt"),
          legend.background = element_blank()),
  FeaturePlot(TLDS_Mesenchymal_seurat, c("Mural.Module.Score"), pt.size = 0.1) &
    scale_color_viridis_c() &
    ylim(c(2.5,10)) &
    theme_pubr(base_size =7) &
    theme(legend.position = c(0, 0.9),
          legend.justification = 0,
          legend.key.width = unit(2,"pt"),
          legend.key.height = unit(5, "pt"),
          legend.background = element_blank()),
      DimPlot(TLDS_Mesenchymal_seurat, reduction = "umap",
            group.by = "Fib_or_Mural",
            label = F, repel = TRUE, pt.size = 0.1) &
      ylim(c(2.5,10)) &
      theme_pubr(base_size = 5) &
      theme(legend.position = c(0,0.9), plot.title = element_blank(),
            legend.key.height = unit(2, "pt"), legend.justification = c(0,0.5),
            legend.background = element_blank()) &
    scale_colour_brewer(palette = "Dark2") &
      guides(colour = guide_legend(override.aes = list(size = 2))),
  widths = c(1,1,1), ncol = 3, align = "h", labels = c("h","i","j"))
Fig_1HIJ



#Supplementary Figure 1F-G - Qian et al Mural cell score validation####
Qian.combined <- SetIdent(Qian.combined, value = "seurat_clusters")
DimPlot(Qian.combined)

test <- aggregate(Qian.combined$Cluster1 ~ Qian.combined$seurat_clusters,
                  FUN = median)
Fibs_cluster <- as.character(test$`Qian.combined$seurat_clusters`)[
  which.max(test$`Qian.combined$Cluster1`)]
Qian_Mesenchymal_seurat <- subset(Qian.combined, idents = Fibs_cluster)
Qian.combined$Stroma <- Qian.combined$seurat_clusters == Fibs_cluster

DefaultAssay(Qian_Mesenchymal_seurat) <- "RNA"
Qian_Mesenchymal_seurat <- AddModuleScore(
  Qian_Mesenchymal_seurat, features = Consensus_FibsvMural_sigs
)
names(Qian_Mesenchymal_seurat@meta.data)
Qian_Mesenchymal_seurat$Fib_or_Mural <- "Undetermined"
Qian_Mesenchymal_seurat$Fib_or_Mural[Qian_Mesenchymal_seurat$Cluster1 > 0 & 
                                       Qian_Mesenchymal_seurat$Cluster1 - Qian_Mesenchymal_seurat$Cluster2 > 0.1] <- "Fibroblast"
Qian_Mesenchymal_seurat$Fib_or_Mural[Qian_Mesenchymal_seurat$Cluster2 > 0 & 
                                       Qian_Mesenchymal_seurat$Cluster2 - Qian_Mesenchymal_seurat$Cluster1 > 0.1] <- "Mural.cells"

Qian_Mesenchymal_seurat$Fibrolast.Module.Score <- Qian_Mesenchymal_seurat$Cluster1
Qian_Mesenchymal_seurat$Mural.Module.Score <- Qian_Mesenchymal_seurat$Cluster2

SupplFig_1FG <- 
  ggarrange(
  DimPlot(Qian.combined, pt.size = 0.1, group.by = "Stroma") & theme_pubr(base_size = 5) &
    theme(plot.title = element_blank()) &
    scale_colour_d3(name = "Stroma") &
    theme(legend.position = c(0,0.1),
          legend.justification = 0,
          legend.key.size = unit(2,"pt"),
          legend.background = element_blank()),
  ggarrange(
    DimPlot(Qian_Mesenchymal_seurat, reduction = "umap",
            group.by = "Fib_or_Mural",
            label = F, repel = TRUE, pt.size = 0.1) &
      xlim(c(7.5,15)) & ylim(c(-7.5,0)) &
      theme_pubr(base_size = 5) &
      theme(legend.position = c(0,0.9), plot.title = element_blank(),
            legend.key.height = unit(2, "pt"), legend.justification = 0.5,
            legend.background = element_blank()) & 
      scale_color_brewer(palette = "Dark2"),
    Qian_Mesenchymal_seurat@meta.data %>%
      ggplot(aes(x = Cluster2, y = Cluster1, colour = Fib_or_Mural)) +
      theme_pubr(base_size = 5) + 
      geom_point(size = 0.1) + 
      theme(legend.position = "bottom", legend.title = element_blank()) +
      xlab("Mural Cells Module Score") + 
      ylab("Fibroblast Module Score")+
      theme(legend.position = c(1,0.9), plot.title = element_blank(),
            legend.key.height = unit(2, "pt"), legend.justification = 0.5,
            legend.background = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size = 2))) + 
      scale_color_brewer(palette = "Dark2"),
    ncol = 2, common.legend = T, legend = "top"),
  widths = c(1,1.5), ncol = 2, align = "h")
SupplFig_1FG

Figure_Panels <- c("Fig_1B", "Fig_1C", "Fig_1D", "Fig_1E", "Fig_1F", "Fig_1G", "Fig_1HIJ",
                  "SupplFig_1A", "SupplFig_1B", "SupplFig_1CD", "SupplFig_1E", "SupplFig_1FG")
rm(list = ls()[!ls() %in% Figure_Panels])
