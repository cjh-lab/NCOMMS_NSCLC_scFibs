#Figure 3####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
library(tidyverse)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(Seurat)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "MxIHC_Zenodo.Rdata"))
load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))
load(paste0(data_directory, "CIBERSORT_Zenodo.Rdata"))
load(paste0(data_directory, "BulkData_Zenodo.Rdata"))
fibs_names <- c("Adventitial", "Alveolar", "Myo")

#Figure 3A IHC Marker Identification####
names(All.Markers)
All.Markers$cluster <- factor(
  All.Markers$cluster,
  levels = levels(All.Markers$cluster)[3:1]
)
All.Markers$label <- All.Markers$gene
All.Markers$label[All.Markers$p_val_adj.Sample > 0.0001 | All.Markers$log2FC.Sample < 0.75 |
                    !All.Markers$HPA_reliability == "Enhanced" |
                    !All.Markers$Detected.Intracellular == T |
                    is.na(All.Markers$HPA_reliability) |
                    is.na(All.Markers$Detected.Intracellular)] <- NA
All.Markers$label2 <- All.Markers$gene
All.Markers$label2[!All.Markers$gene %in% c("AOC3", "CD34", "POSTN", "ACTA2")] <- NA
All.Markers$label[All.Markers$label %in% c("AOC3", "CD34", "POSTN", "ACTA2")] <- NA

Fig_3A <- 
  All.Markers[All.Markers$p_val_adj.Sample < 0.01, ] %>%
  drop_na(p_val_adj.Sample) %>%
  ggplot(aes(x = log2FC.Sample, y = -log10(p_val_adj.Sample) ,
             size = pct.1, 
             color = HPA_reliability == "Enhanced" & Detected.Intracellular == T &
               !is.na(HPA_reliability) & !is.na(Detected.Intracellular))) +
  theme_pubr(base_size = 7) + 
  geom_point() +
  scale_size(range = c(0.1,2), name = "Detection rate") +
  scale_color_manual(name = "IHC QC passed", values = c("grey70", "red")) +
  facet_wrap(~cluster, ncol = 3) +
  ggrepel::geom_text_repel(aes(label = label), colour = "black", min.segment.length = 0, size = 2, nudge_x = 0.2, fontface = "italic") + 
  ggrepel::geom_label_repel(aes(label = label2), colour = "black", label.padding = 0.15,
                            min.segment.length = 0, size = 3, nudge_x = 0.2, fontface = "italic") +
  theme(legend.position = "right") +
  ylab("Sample level -log10(adj.P)") +
  xlab("Sample level Average log2(FC)")
Fig_3A


#Supplementary Figure 3A - scRNA-seq marker expression####
levels(Fibs.integrated$Fibs_MajorClusters)
Fibs.integrated <- SetIdent(Fibs.integrated, value = "Fibs_MajorClusters")
Fibs.integrated_samples <- AverageExpression(Fibs.integrated, assays = "RNA", return.seurat = T, 
                                             group.by = c("Fibs_MajorClusters", "SampleID2"),
                                             slot = "data")
Fibs.integrated_samples <- NormalizeData(Fibs.integrated_samples)
#extract meta data from sample names
Fibs.integrated_samples@meta.data <- 
  Fibs.integrated_samples@meta.data %>%
  mutate(Fibs_MajorClusters = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,1],
         Dataset = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,2],
         SampleID = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,3])

SupplFig_3A <- 
  ggarrange(
    Plot_SCandSAMPLE_exprs(gene = "CD34", sc_seurat = Fibs.integrated, sample_seurat = Fibs.integrated_samples,
                       ident = "Fibs_MajorClusters"),
    Plot_SCandSAMPLE_exprs(gene = "AOC3", sc_seurat = Fibs.integrated, sample_seurat = Fibs.integrated_samples,
                       ident = "Fibs_MajorClusters"),
    Plot_SCandSAMPLE_exprs(gene = "POSTN", sc_seurat = Fibs.integrated, sample_seurat = Fibs.integrated_samples,
                       ident = "Fibs_MajorClusters"),
    Plot_SCandSAMPLE_exprs(gene = "ACTA2", sc_seurat = Fibs.integrated, sample_seurat = Fibs.integrated_samples,
                       ident = "Fibs_MajorClusters"),
    ncol = 4)
SupplFig_3A


#Supplementary Figure 3 B&C - MxIHC feature plots####
names(MxIHC_data)
MxIHC_Fibroblasts <- MxIHC_data[
  log(MxIHC_data$`Area.of.sub.objects.Nucleus..2...µm².`) > 3 &
    log(MxIHC_data$`Area.of.sub.objects.Nucleus..2...µm².`) < 4.5 &
    MxIHC_data$PanCK < 0.05 &
    MxIHC_data$MCAM < 0.05 &
    apply(MxIHC_data[,c("CD34", "POSTN", "ACTA2", "AOC3")],1,max) > 0.1, ]
MxIHC_Fibroblasts$Cell.type <- factor(MxIHC_Fibroblasts$Cell.type, levels = fibs_names)

#Not run - Result pre-saved in Rdata file 
#MxIHC_Fibroblasts.ds <- Fibroblasts[sample(1:nrow(Fibroblasts), 0.01*nrow(Fibroblasts)), ]
#umap_dr <- umap::umap(MxIHC_Fibroblasts.ds[,c("CD34", "AOC3", "POSTN", "ACTA2")])
#MxIHC_Fibroblasts.ds$UMAP1 <- umap_dr$layout[,1]
#MxIHC_Fibroblasts.ds$UMAP2 <- umap_dr$layout[,2]

SupplFig_3B <- 
  MxIHC_Fibroblasts.ds %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = Cell.type)) +
  geom_point(size = 0.1)+
  scale_colour_manual(values = Fibs_col.palette) +
  theme_pubr(base_size = 7) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1, "pt"))+
  guides(colour = guide_legend(override.aes = list(size = 2)))
SupplFig_3B

ACTA2.Fplot <- 
  MxIHC_Fibroblasts.ds %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = ACTA2 > 0.1)) +
  geom_point(size = 0.1)+
  #scale_color_viridis_c(option = "B", breaks = c(0,0.25,0.5,0.75)) +
  scale_color_viridis_d()+
  theme_pubr(base_size = 6) +
  guides(colour = guide_legend(title = "Area%>10", override.aes = list(size = 2))) +
  ggtitle("ACTA2")

POSTN.Fplot <- 
  MxIHC_Fibroblasts.ds %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = POSTN > 0.1)) +
  geom_point(size = 0.1)+
  #scale_color_viridis_c(option = "B", breaks = c(0,0.25,0.5,0.75)) +
  scale_color_viridis_d()+
  theme_pubr(base_size = 6) +
  guides(colour = guide_legend(title = "Area%>10", override.aes = list(size = 2))) +
  ggtitle("POSTN")

Cd34.Fplot <- 
  MxIHC_Fibroblasts.ds %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = CD34 > 0.1)) +
  geom_point(size = 0.1)+
  #scale_color_viridis_c(option = "B", breaks = c(0,0.25,0.5,0.75)) +
  scale_color_viridis_d()+
  theme_pubr(base_size = 6) +
  guides(colour = guide_legend(title = "Area%>10", override.aes = list(size = 2))) +
  ggtitle("CD34")

AOC3.Fplot <- 
  MxIHC_Fibroblasts.ds %>%
  ggplot(aes(x = UMAP1, y = UMAP2, colour = AOC3 > 0.1)) +
  geom_point(size = 0.1)+
  #scale_color_viridis_c(option = "B", breaks = c(0,0.25,0.5,0.75)) +
  scale_color_viridis_d()+
  theme_pubr(base_size = 6) +
  guides(colour = guide_legend(title = "Area%>10", override.aes = list(size = 2))) +
  ggtitle("AOC3")

SupplFig_3C <- ggarrange(Cd34.Fplot,AOC3.Fplot, POSTN.Fplot, ACTA2.Fplot,
                         ncol = 4, nrow = 1, common.legend = T, legend = "right")
SupplFig_3C

#CIBERSORTx Validation####
#Supplementary Figure 3D####
names(CSx_SignatureMtx)
CSx_SignatureMtx <- CSx_SignatureMtx[, names(CSx_SignatureMtx)[c(1,2,4,7,6,3,5,8)]]
CSx_SignatureMtx$max <- apply(CSx_SignatureMtx[,2:8], 1, which.max)
SupplFig_3D <- pheatmap::pheatmap(CSx_SignatureMtx[order(CSx_SignatureMtx$max),2:8], scale = "row",
                   cluster_rows = F, cluster_cols = F, show_rownames = F,
                   color = viridis::viridis(101),
                   fontsize = 7, 
                  )

SupplFig_3D

#Supplementary Figure 3E####
names(PB.Sim_CSx.res) <- paste("CSx_", names(PB.Sim_CSx.res), sep = "")
rownames(PB.Sim_CSx.res) <- PB.Sim_CSx.res$CSx_Mixture
cell_names <- colnames(PseudoBulk_1_metaData)

Compare_res <- merge(PseudoBulk_1_metaData,
                     PB.Sim_CSx.res[,],
                     by = 0)
Cell.level_cor <- list()
Cell.level_p <- list()
for (i in cell_names) {
  Cell.level_cor[[i]] <- cor(Compare_res[[i]], Compare_res[[paste("CSx_", i, sep = "")]])
  Cell.level_p[[i]] <- cor.test(Compare_res[[i]], Compare_res[[paste("CSx_", i, sep = "")]])$p.value
}  
Cell.level_cor_df <- as.data.frame(do.call(rbind, Cell.level_cor))
Cell.level_p_df <- do.call(rbind, Cell.level_p)
Cell.level_cor_df$cell_type <- names(Cell.level_cor)
Cell.level_cor_df$p <- Cell.level_p_df

cell_names
Cor_barplot <- 
  Cell.level_cor_df %>% #[Cell.level_cor_df$cell_type %in% fibs_names, ] %>%
  ggplot(aes(y = reorder(cell_type, V1), x = V1, label = paste("p =", signif(p,3))))  +
  geom_bar(stat = "identity", aes(fill = cell_type)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  xlab("CIBERSORTx:Ground truth (r)") +
  ylab("Cell type") +
  xlim(c(0,1)) +
  geom_text(hjust = 1.1, size = 2) +
  #ggtitle("Simulated PseudoBulk (n=200)") +
  scale_fill_manual(values = Fibs_col.palette, na.value = "grey90") 
Cor_barplot

names(Compare_res)
Myo_scatter <- 
  Compare_res %>%
  ggplot(aes(y = Myo, x = CSx_Myo)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm", formula = y ~ x) +
  theme_pubr(base_size = 7) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label =  paste(stat(rr.label))),
                        parse = TRUE, size = 3)    +
  xlab("Estimated Myo abundance\n(CSx Abs score)") +
  ylab("True Myo abundance\n(#cells per sample)") 


Alveolar_scatter <- 
  Compare_res %>%
  ggplot(aes(y = Alveolar, x = CSx_Alveolar)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm", formula = y ~ x) +
  theme_pubr(base_size = 7) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label =  paste(stat(rr.label))),
                        parse = TRUE, size = 3)   +
  xlab("Estimated Alveolar abundance\n(CSx Abs score)") +
  ylab("True Alveolar abundance\n(#cells per sample)")  

Adventitial_scatter <- 
  Compare_res %>%
  ggplot(aes(y = Adventitial, x = CSx_Adventitial)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm", formula = y ~ x) +
  theme_pubr(base_size = 7) +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label =  paste(stat(rr.label))),
                        parse = TRUE, size = 3)  +
  xlab("Estimated Adventitial abundance\n(CSx Abs score)") +
  ylab("True Adventitial abundance\n(#cells per sample)")

SupplFig_3E <- ggarrange(Cor_barplot, Adventitial_scatter, Alveolar_scatter, Myo_scatter, ncol = 2, nrow = 2, align = "hv")
SupplFig_3E

#Supplementary Figure 3F-H ####
MxIHC_Fibroblasts$TvN.Grouping.ID <- paste(MxIHC_Fibroblasts$Slide.name, MxIHC_Fibroblasts$Sample.Region.TvN, sep = "_")
dt <- as.table(as.matrix(table(MxIHC_Fibroblasts$TvN.Grouping.ID,
                               MxIHC_Fibroblasts$Cell.type)))
Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
Sample.pct_long$Slide.name <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,1]
Sample.pct_long$TvN <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,2]
TvN.Data_long <- merge(MxIHC_Sample.metaData, Sample.pct_long,
                       by = "Slide.name")
names(TvN.Data_long)[names(TvN.Data_long) == "Freq"] <- "Region.pct"
names(TvN.Data_long)[names(TvN.Data_long) == "Var2"] <- "Fibs_SubPop"

TvN.Data_long$SubtypeTvN <- TvN.Data_long$Subtype
TvN.Data_long$SubtypeTvN[TvN.Data_long$TvN == "Control"] <- "Control"
TvN.Data_long <- TvN.Data_long[!TvN.Data_long$TvN == "NA", ]

SupplFig_3F <- 
  TvN.Data_long[] %>%
  filter(SubtypeTvN == "Control" ) %>%
  filter(!Slide.name == "Drop-Seq Cohort 1") %>%
  ggplot(aes(x = Subtype, y = Region.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none")+
  facet_grid(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  #ggrepel::geom_text_repel(aes(label = Slide.name)) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,115)) +
  stat_compare_means(comparisons = list(c("LUAD", "LUSC")), size = 2.5,
                     #label = "p.signif",
                     method = "wilcox.test", label.y = 100) +
  scale_fill_manual(values = Fibs_col.palette) +
  xlab("Adjacent Tumour Subtype")
SupplFig_3F

names(TvN.Data_long)
SupplFig_3H <- 
  TvN.Data_long[] %>%
  filter(SubtypeTvN == "Control" ) %>%
  filter(!Slide.name == "Drop-Seq Cohort 1") %>%
  ggplot(aes(x = Control...Fibrosis.Inflammation == 1, y = Region.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_wrap(~Fibs_SubPop, ncol = 3) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  #ggrepel::geom_text_repel(aes(label = Slide.name)) +
  #geom_text(aes(label = Normal.tissue.comment), size = 2) +
  rotate_x_text(angle = 45) +
  scale_fill_manual(values = Fibs_col.palette) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,115)) +
  stat_compare_means(comparisons = list(c("FALSE", "TRUE")), size = 2.5,
                     #label = "p.signif",
                     label.y = 100,
                     method = "wilcox.test") +
  xlab("Histological evidence of fibrosis/inflammation")
SupplFig_3H


SupplFig_3G <- 
  MxIHC_Sample.metaData %>%
  #drop_na(Control...Fibrosis.Inflammation) %>%
  filter(!Slide.name == "Drop-Seq Cohort 1") %>%
  filter(!Control...Fibrosis.Inflammation == "NA") %>%
  ggplot(aes(x = Subtype, fill = factor(Control...Fibrosis.Inflammation, levels = c(0,1), labels = c(F, T)))) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(name = "Fibrosis or\nInflammation\nevident") +
  theme_pubr(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(5, "pt")) +
  #ggtitle("Fibrosis or Inflammation") +
  rotate_x_text(angle = 45) + ylab("Fraction")
SupplFig_3G

#Figure 3E####
Fibs.integrated$SampleID2 <- paste(Fibs.integrated$Dataset, Fibs.integrated$SampleID, sep = "_")
Sample_MetaData$SampleID2 <- paste(Sample_MetaData$Dataset, Sample_MetaData$SampleID, sep = "_")
Sample_MetaData$Sample.type[Sample_MetaData$Sample.type == "Tumor"] <- "Tumour"
Sample_MetaData$Sample.Subtype2 <- droplevels(Sample_MetaData$Sample.Subtype, exclude = c("Carcinoid" , "Pleiomorphic.carcinoma", "LCLC") )
dt <- as.table(as.matrix(table(Fibs.integrated$SampleID2,
                               Fibs.integrated$Fibs_MajorClusters)))

Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
Sample.pct <- reshape2::dcast(as.data.frame(dt), formula = Var1 ~ Var2)
rownames(Sample.pct) <- Sample.pct$Var1
Sample.pct <- Sample.pct[,2:4]/rowSums(Sample.pct[,2:4])*100


scRNAseq_MetaData2 <- merge(Sample_MetaData, Sample.pct_long,
                          by.x = "SampleID2", by.y = "Var1")
names(scRNAseq_MetaData2)[23] <- "Residual"
names(scRNAseq_MetaData2)[22] <- "Fibs_SubPop"
scRNAseq_MetaData2$Fibs_SubPop <- factor(
  as.character(scRNAseq_MetaData2$Fibs_SubPop),
  levels = c("Adventitial", "Alveolar", "Myo")
)

my.comparisons <- list(c("Control", "LUAD"), c("Control", "LUSC"), c("LUAD", "LUSC"))
Fig_3E <- 
  scRNAseq_MetaData2 %>%
  drop_na(Sample.Subtype2) %>%
  ggplot(aes(x = Sample.Subtype2, y = Residual)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  facet_grid(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  scale_fill_manual(values = Fibs_col.palette) +
  geom_jitter(alpha = 0.5, width = 0.15, size = 0.1) +
  rotate_x_text(angle = 45) +
  ylab("% of Fibroblasts per sample\n(scRNA-seq)") +
  #geom_text(aes(label = Slide.name), size = 2) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,135)) +
  stat_compare_means(comparisons = my.comparisons,
                     size = 2.5,
                     #label = "p.signif",
                     method = "wilcox.test", tip.length = 0.01) 
Fig_3E

#Figure 3F####
MxIHC_Fibroblasts$TvN.Grouping.ID <- paste(MxIHC_Fibroblasts$Slide.name, MxIHC_Fibroblasts$Sample.Region.TvN, sep = "_")
dt <- as.table(as.matrix(table(MxIHC_Fibroblasts$TvN.Grouping.ID,
                               MxIHC_Fibroblasts$Cell.type)))
Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
Sample.pct_long$Slide.name <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,1]
Sample.pct_long$TvN <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,2]

TvN.Data_long <- merge(MxIHC_Sample.metaData, Sample.pct_long,
                       by = "Slide.name")
names(TvN.Data_long)[names(TvN.Data_long) == "Freq"] <- "Region.pct"
names(TvN.Data_long)[names(TvN.Data_long) == "Var2"] <- "Fibs_SubPop"

TvN.Data_long$SubtypeTvN <- TvN.Data_long$Subtype
TvN.Data_long$SubtypeTvN[TvN.Data_long$TvN == "Control"] <- "Control"
TvN.Data_long <- TvN.Data_long[!TvN.Data_long$TvN == "NA", ]

names(TvN.Data_long)
Fig_3F <- 
  TvN.Data_long[] %>%
  ggplot(aes(x = SubtypeTvN, y = Region.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  facet_grid(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  scale_fill_manual(values = Fibs_col.palette) +
  geom_jitter(alpha = 0.5, width = 0.15, size = 0.1) +
  rotate_x_text(angle = 45) +
  ylab("% of Fibroblasts per region\n(MxIHC)") +
  #geom_text(aes(label = Slide.name), size = 2) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,135)) +
  stat_compare_means(comparisons = my.comparisons, size = 2.5,
                     #label = "p.signif",
                     method = "wilcox.test", tip.length = 0.01) 
Fig_3F

#Figure 3G####
filteredTraits_long <- reshape2::melt(filteredTraits, id.vars = c("Sample.Subtype", "Grade", "LUAD.MolSubtype_recalculated", "TP53.Status"),
                                      measure.vars = c("Adventitial_Fibs.pct", "Alveolar_Fibs.pct", "Myo_Fibs.pct"))
filteredTraits_long$Fibs_SubPop <- factor(
  filteredTraits_long$variable,
  levels = levels(filteredTraits_long$variable),
  labels = c("Adventitial", "Alveolar", "Myo")
)
Fig_3G <- 
  filteredTraits_long %>%
  drop_na(Sample.Subtype) %>%
  ggplot(aes(x = Sample.Subtype, y = value)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  facet_grid(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  scale_fill_manual(values = Fibs_col.palette) +
  geom_jitter(alpha = 0.5, width = 0.15, size = 0.1) +
  rotate_x_text(angle = 45) +
  ylab("% of Fibroblasts per sample\n(CIBERSORTx)") +
  #geom_text(aes(label = Slide.name), size = 2) +
  scale_y_continuous(breaks = c(0,25,50, 75, 100),
                     limits = c(0,135)) +
  stat_compare_means(comparisons = my.comparisons, size = 3,
                     label = "p.signif",
                     method = "wilcox.test", tip.length = 0.01) 
Fig_3G

# Figure 3C - Subtype MxIHC examples####
#LUAD example
Slide = "Lung Adeno 1"
sample_data <- MxIHC_Fibroblasts %>% filter(Sample == Slide)
XY_cols <- names(sample_data)[1:2]
convex.hull.idx <- chull(as.matrix(sample_data[sample_data$Sample.Region.TvN == "Tumour", XY_cols]))
convex.hull.pts <- rownames(sample_data[sample_data$Sample.Region.TvN == "Tumour",])[convex.hull.idx]
LUAD_eg <- 
  sample_data %>%
  slice_sample(prop = 0.25) %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl., colour = Cell.type)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "none") +
  scale_color_manual(values = Fibs_col.palette) +
  geom_polygon(data = MxIHC_Fibroblasts[convex.hull.pts, XY_cols],
               colour = "black", alpha = 0, linetype = "dotted", size = 2) 
LUAD_eg

#LUSC example
Slide = "Drop-Seq Cohort 10"
sample_data <- MxIHC_Fibroblasts %>% filter(Sample == Slide)
XY_cols <- names(sample_data)[1:2]
convex.hull.idx <- chull(as.matrix(sample_data[sample_data$Sample.Region.TvN == "Tumour", XY_cols]))
convex.hull.pts <- rownames(sample_data[sample_data$Sample.Region.TvN == "Tumour",])[convex.hull.idx]
LUSC_eg <- 
  sample_data %>%
  slice_sample(prop = 0.25) %>%
  ggplot(aes(x = X.Center..Pxl., y = Y.Center..Pxl., colour = Cell.type)) +
  geom_point(size = 0.1) +
  theme_void(base_size = 7) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "none") +
  scale_color_manual(values = Fibs_col.palette) +
  geom_polygon(data = MxIHC_Fibroblasts[convex.hull.pts, XY_cols],
               colour = "black", alpha = 0, linetype = "dotted", size = 2) 
LUSC_eg

layer_scales(LUSC_eg)$x$range$range
xmax.diff <- max(layer_scales(LUSC_eg)$x$range$range) - max(layer_scales(LUAD_eg)$x$range$range)
LUAD_xMax <- max(layer_scales(LUAD_eg)$x$range$range)
ymax.diff <- max(layer_scales(LUSC_eg)$y$range$range) - max(layer_scales(LUAD_eg)$y$range$range)
LUAD_yMax <- max(layer_scales(LUAD_eg)$y$range$range)

Fig_3C <- 
  ggarrange(LUAD_eg + xlim(c(-0.5*xmax.diff, LUAD_xMax + 0.5*xmax.diff)) +
            ylim(c(-0.5*ymax.diff, LUAD_yMax + 0.5*ymax.diff)),
            LUSC_eg , ncol = 1, nrow = 2)


Fig_3C
