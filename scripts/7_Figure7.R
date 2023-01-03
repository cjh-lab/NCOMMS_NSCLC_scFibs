#Figure 7####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
options(stringsAsFactors = F)
pardefault <- par()
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(survival)
library(survminer)
library(forestmodel)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "BulkData_Zenodo.Rdata"))
load(paste0(data_directory, "MxIHC_Zenodo.Rdata"))
load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))

LM22_names <- names(Merged.LUADtraits)[32:53]
fibs_names <- c("Adventitial_Fibs.pct", "Alveolar_Fibs.pct", "Myo_Fibs.pct")
#Figure 7A####
#G1 example
Slide = "Drop-Seq Cohort 8"
sample_data <- MxIHC_Fibroblasts %>% filter(Sample == Slide)
XY_cols <- names(sample_data)[1:2]
convex.hull.idx <- chull(as.matrix(sample_data[sample_data$Sample.Region.TvN == "Tumour", XY_cols]))
convex.hull.pts <- rownames(sample_data[sample_data$Sample.Region.TvN == "Tumour",])[convex.hull.idx]
G1_eg <- 
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
G1_eg

#G2 example
Slide = "Lung Adeno 2"
sample_data <- MxIHC_Fibroblasts %>% filter(Sample == Slide)
XY_cols <- names(sample_data)[1:2]
convex.hull.idx <- chull(as.matrix(sample_data[sample_data$Sample.Region.TvN == "Tumour", XY_cols]))
convex.hull.pts <- rownames(sample_data[sample_data$Sample.Region.TvN == "Tumour",])[convex.hull.idx]
G2_eg <- 
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
G2_eg

#G3 example
Slide = "Lung Adeno 5"
sample_data <- MxIHC_Fibroblasts %>% filter(Sample == Slide)
XY_cols <- names(sample_data)[1:2]
convex.hull.idx <- chull(as.matrix(sample_data[sample_data$Sample.Region.TvN == "Tumour", XY_cols]))
convex.hull.pts <- rownames(sample_data[sample_data$Sample.Region.TvN == "Tumour",])[convex.hull.idx]
G3_eg <- 
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
G3_eg


which.max(c(max(layer_scales(G1_eg)$x$range$range),
            max(layer_scales(G2_eg)$x$range$range),
            max(layer_scales(G3_eg)$x$range$range)))
which.max(c(max(layer_scales(G1_eg)$y$range$range),
            max(layer_scales(G2_eg)$y$range$range),
            max(layer_scales(G3_eg)$y$range$range)))
xmax <- max(layer_scales(G1_eg)$x$range$range)
ymax <- max(layer_scales(G2_eg)$y$range$range)


Fig_7A <- 
  ggarrange(G1_eg + expand_limits(y = c(median(layer_scales(G1_eg)$y$range$range) - 0.5*ymax,
                                        median(layer_scales(G1_eg)$y$range$range) + 0.5*ymax)),
            G2_eg + expand_limits(x = c(median(layer_scales(G2_eg)$x$range$range) - 0.5*xmax,
                                        median(layer_scales(G2_eg)$x$range$range) + 0.5*xmax)),
            G3_eg + expand_limits(x = c(median(layer_scales(G3_eg)$x$range$range) - 0.5*xmax,
                                        median(layer_scales(G3_eg)$x$range$range) + 0.5*xmax),
                                  y = c(median(layer_scales(G3_eg)$y$range$range) - 0.5*ymax,
                                        median(layer_scales(G3_eg)$y$range$range) + 0.5*ymax)),
            nrow = 1, ncol = 3)
Fig_7A


#Figure 7B####
Bulk_Fibs.pct <- reshape2::melt(Merged.LUADtraits, id.vars = c("Grade2", "Dataset.factor", "Stage", "Stage_4cat", "TP53.Status", "EGFR.Status", "KRAS.Status","LUAD.MolSubtype_factor"),
                                measure.vars = c("Myo_Fibs.pct", "Alveolar_Fibs.pct", "Adventitial_Fibs.pct") 
)
Bulk_Fibs.pct$Fibs_SubPop <- factor(
  Bulk_Fibs.pct$variable,
  levels = levels(Bulk_Fibs.pct$variable)[c(3,2,1)],
  labels = c("Adventitial", "Alveolar", "Myo")
)
names(Bulk_Fibs.pct)[names(Bulk_Fibs.pct) == "value"] <- "Fibs.pct"
my.comparisons <- list(c("G1/G2", "G3"))
Fig_7B <- 
  Bulk_Fibs.pct %>%
  drop_na(Grade2) %>%
  ggplot(aes(x = Grade2, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  facet_wrap(~Fibs_SubPop) +
  theme(legend.position = "none") +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = my.comparisons, size = 2.5) +
  xlab("Tumour Grade") + 
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
Fig_7B

#Supplementary Figure 7A####
SupplFig_7A <- 
  Bulk_Fibs.pct %>%
  drop_na(Grade2) %>%
  ggplot(aes(x = Grade2, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_grid(Dataset.factor~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = my.comparisons, size = 2.5) +
  xlab("Tumour Grade") + 
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
SupplFig_7A

#Figure 7C####
Fibs.integrated$SampleID2 <- paste(Fibs.integrated$Dataset, Fibs.integrated$SampleID, sep = "_")
Sample_MetaData$SampleID2 <- paste(Sample_MetaData$Dataset, Sample_MetaData$SampleID, sep = "_")
Sample_MetaData$Sample.type[Sample_MetaData$Sample.type == "Tumor"] <- "Tumour"
Sample_MetaData$Sample.Subtype2 <- droplevels(Sample_MetaData$Sample.Subtype, exclude = c("Carcinoid" , "Pleiomorphic.carcinoma", "LCLC") )
dt <- as.table(as.matrix(table(Fibs.integrated$SampleID2,
                               Fibs.integrated$Fibs_MajorClusters)))

Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
scRNAseq_Fibs.pct <- merge(Sample_MetaData, Sample.pct_long,
                          by.x = "SampleID2", by.y = "Var1")
names(scRNAseq_Fibs.pct)[23] <- "Fibs.pct"
names(scRNAseq_Fibs.pct)[22] <- "Fibs_SubPop"
scRNAseq_Fibs.pct$Fibs_SubPop <- factor(
  scRNAseq_Fibs.pct$Fibs_SubPop,
  levels = c("Adventitial", "Alveolar", "Myo")
)


scRNAseq_Fibs.pct$Differentiation2 <- factor(scRNAseq_Fibs.pct$Differentiation,
                                             levels = c("G1", "G2", "G3"),
                                             labels = c("G1/G2", "G1/G2", "G3"))
my.comparisons <- list(c("G1/G2", "G3"))
Fig_7C <-
  scRNAseq_Fibs.pct %>%
  filter(Sample.Subtype2 == "LUAD") %>%
  drop_na(Differentiation2) %>%
  ggplot(aes(x = Differentiation2, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_wrap(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = my.comparisons, size = 2.5) +
  xlab("Tumour Grade") + 
  ylab("% of fibroblasts per sample\n(scRNA-seq)")
Fig_7C


#Figure 7D####
MxIHC_Fibroblasts$TvN.Grouping.ID <- paste(MxIHC_Fibroblasts$Slide.name, MxIHC_Fibroblasts$Sample.Region.TvN, sep = "_")
dt <- as.table(as.matrix(table(MxIHC_Fibroblasts$TvN.Grouping.ID,
                               MxIHC_Fibroblasts$Cell.type)))
Sample.pct_long <- as.data.frame(dt/rowSums(dt)*100)
Sample.pct_long$Slide.name <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,1]
Sample.pct_long$TvN <- do.call(rbind, strsplit(as.character(Sample.pct_long$Var1), "_", fixed = T))[,2]
TvN.Data_long <- merge(MxIHC_Sample.metaData, Sample.pct_long,
                       by = "Slide.name")
names(TvN.Data_long)[names(TvN.Data_long) == "Freq"] <- "Fibs.pct"
names(TvN.Data_long)[names(TvN.Data_long) == "Var2"] <- "Fibs_SubPop"

TvN.Data_long$SubtypeTvN <- TvN.Data_long$Subtype
TvN.Data_long$SubtypeTvN[TvN.Data_long$TvN == "Control"] <- "Control"

MxIHC_Fibs.pct <- TvN.Data_long[!TvN.Data_long$TvN == "NA", ]
names(MxIHC_Fibs.pct)
MxIHC_Fibs.pct$Differentiation2 <- factor(MxIHC_Fibs.pct$Grade,
                                          levels = c("G1", "G2", "G3"),
                                          labels = c("G1/G2", "G1/G2", "G3"))
Fig_7D <- 
  MxIHC_Fibs.pct %>%
  filter(SubtypeTvN == "LUAD") %>%
  drop_na(Differentiation2) %>%
  ggplot(aes(x = Differentiation2, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_wrap(~Fibs_SubPop) +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = my.comparisons, size = 2.5, label.y = 100) +
  xlab("Tumour Grade") + 
  ylab("% of fibroblasts per sample\n(MxIHC)")
Fig_7D



#Figure 7E####
Fig_7E <- 
  Bulk_Fibs.pct %>%
  drop_na(TP53.Status) %>%
  ggplot(aes(x = TP53.Status, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_grid(~Fibs_SubPop, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = list(c("0","1")), size = 2.5,
                     label.y = 110) +
  xlab("TP53 Status") + 
  scale_x_discrete(labels = c("wt", "mut")) +
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
Fig_7E

#Supplementary Figure 7E####
SupplFig_7E <- 
  Bulk_Fibs.pct %>%
  drop_na(TP53.Status) %>%
  ggplot(aes(x = TP53.Status, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_grid(Dataset.factor~Fibs_SubPop, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  stat_compare_means(comparisons = list(c("0","1")), size = 2.5,
                     label.y = 110) +
  xlab("TP53 Status") + 
  scale_x_discrete(labels = c("wt", "mut")) +
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
SupplFig_7E

#Figure 7F####
Fig_7F <- 
  Bulk_Fibs.pct %>%
  drop_na(LUAD.MolSubtype_factor) %>%
  ggplot(aes(x = LUAD.MolSubtype_factor, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_grid(~Fibs_SubPop, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  xlab("Molecular Subtype") + 
  scale_x_discrete(labels = c("TRU", "PP", "PI")) +
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
Fig_7F

#SUpplementary Figure 7D####
SupplFig_7D <- 
  Bulk_Fibs.pct %>%
  drop_na(LUAD.MolSubtype_factor) %>%
  ggplot(aes(x = LUAD.MolSubtype_factor, y = Fibs.pct)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  facet_grid(Dataset.factor~Fibs_SubPop, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, aes(fill = Fibs_SubPop)) +
  geom_jitter(alpha = 0.5, width = 0.2, size = 0.1) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  scale_y_continuous(limits = c(0,125), breaks = c(0,25,50,75,100)) +
  xlab("Molecular Subtype") + 
  scale_x_discrete(labels = c("TRU", "PP", "PI")) +
  ylab("% of fibroblasts per sample\n(CIBERSORTx)")
SupplFig_7D

#Figure 7G####
LM22_names <- names(Merged.LUADtraits)[32:53]
fibs_names <- c("Adventitial_Fibs.pct", "Alveolar_Fibs.pct", "Myo_Fibs.pct")

Immune_cor <- WGCNA::corAndPvalue(Merged.LUADtraits[, LM22_names],
                         Merged.LUADtraits[, fibs_names], 
                         use = "pairwise.complete.obs",
                         method = "spearman")
Immune_cor$sig_cor <- Immune_cor$cor
Immune_cor$sig_cor[Immune_cor$p > 0.001] <- NA
Immune_cor$cor_df <- data.frame(
  Immune_cor$cor,
  Immune.cell = rownames(Immune_cor$cor)
)
Immune_cor$cor_df$Immune.cell.label <- Immune_cor$cor_df$Immune.cell
sig_cors <- rownames(na.omit(Immune_cor$sig_cor[,2:3]))
Immune_cor$cor_df$Immune.cell.label[!Immune_cor$cor_df$Immune.cell %in% sig_cors] <- NA

Fig_7G <- 
  Immune_cor$cor_df %>%
  ggplot(aes(y = Alveolar_Fibs.pct, x = Myo_Fibs.pct,
             colour = Alveolar_Fibs.pct > Myo_Fibs.pct)) +
  theme_pubr(base_size = 7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(size = abs(Myo_Fibs.pct))) +
  scale_size(range = c(0.1, 2)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c(Fibs_col.palette[[2]],Fibs_col.palette[[3]])) +
  ggrepel::geom_text_repel(aes(label = Immune.cell.label),
                           size = 2.5, min.segment.length = 0) +
  xlab("Cor to Myo (rho)") +
  ylab("Cor to Alveolar (rho)")
Fig_7G

#Supplementary Figure 7F####
Immune_cor.list <- list()
Immune_corplots.list <- list()
Merged.LUADtraits.list <- split(Merged.LUADtraits, f = Merged.LUADtraits$Dataset.factor)
for(i in names(Merged.LUADtraits.list)){
  Immune_cor.list[[i]] <- WGCNA::corAndPvalue(Merged.LUADtraits.list[[i]][, LM22_names],
                                              Merged.LUADtraits.list[[i]][, fibs_names], 
                                              use = "pairwise.complete.obs",
                                              method = "spearman")
  Immune_cor.list[[i]]$sig_cor <- Immune_cor$cor
  Immune_cor.list[[i]]$sig_cor[Immune_cor$p > 0.001] <- NA
  Immune_cor.list[[i]]$cor_df <- data.frame(
    Immune_cor.list[[i]]$cor,
    Immune.cell = rownames(Immune_cor.list[[i]]$cor)
  )
  sig_cors <- rownames(na.omit(Immune_cor.list[[i]]$sig_cor[,2:3]))
  Immune_cor.list[[i]]$cor_df$Immune.cell.label <- Immune_cor.list[[i]]$cor_df$Immune.cell
  Immune_cor.list[[i]]$cor_df$Immune.cell.label[!Immune_cor.list[[i]]$cor_df$Immune.cell %in% sig_cors] <- NA
}

cor_df.list <- list()
for(i in names(Immune_cor.list)){
  cor_df.list[[i]] <- Immune_cor.list[[i]]$cor_df
  cor_df.list[[i]]$DS <- i
}
cor_df_long <- do.call(rbind, cor_df.list)

Adv_immune.cor <- 
  cor_df_long %>%
  ggplot(aes(x = reorder(Immune.cell, Adventitial_Fibs.pct, na.rm = T),
             y = Adventitial_Fibs.pct)) +
  theme_pubr(base_size = 7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, width = 0.2) +
  coord_flip() +
  xlab("LM22 Immune cell subtype") +
  ylab("Correlation to Adventitial (rho)")

Alv_immune.cor <- 
  cor_df_long %>%
  ggplot(aes(x = reorder(Immune.cell, Alveolar_Fibs.pct, na.rm = T),
             y = Alveolar_Fibs.pct)) +
  theme_pubr(base_size = 7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, width = 0.2) +
  coord_flip() +
  xlab("LM22 Immune cell subtype") +
  ylab("Correlation to Alveolar (rho)")

Myo_immune.cor <- 
  cor_df_long %>%
  ggplot(aes(x = reorder(Immune.cell, Myo_Fibs.pct, na.rm = T),
             y = Myo_Fibs.pct)) +
  theme_pubr(base_size = 7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, width = 0.2) +
  coord_flip() +
  xlab("LM22 Immune cell subtype") +
  ylab("Correlation to Myo (rho)")

SupplFig_7F <- ggarrange(Adv_immune.cor, Alv_immune.cor, Myo_immune.cor,
          nrow = 3)
SupplFig_7F

