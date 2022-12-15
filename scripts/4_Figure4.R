#Figure 4####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
library(Seurat)
library(ggplot2)
library(WGCNA)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(destiny)
library(metap)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "Trajectory Analysis results.Rdata"))
load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))
load(paste0(data_directory, "MxIHC_Zenodo.Rdata"))
Fibs.integrated <- AddMetaData(Fibs.integrated, DPT_output)

A.M_pT_res_stats2$Gene <- rownames(A.M_pT_res_stats2)
T.M_pT_res_stats2$Gene <- rownames(T.M_pT_res_stats2)

#Figure 4A-C####
pA <- ggplot(DM_dim_data, aes(x = DC1, y = DC2, colour = Fibs_subpopulation)) +
  theme_pubr(base_size = 7) +
  theme(axis.text = element_blank(), legend.title = element_blank(), legend.position = "right") +
  geom_point(size = 0.1) +
  scale_color_manual(values = Fibs_col.palette) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
pB <- ggplot(DM_dim_data, aes(x = DC1, y = DC2, colour = factor(dpt_branch))) +
  theme_pubr(base_size = 7) +
  theme(axis.text = element_blank(), legend.position = "right") +
  geom_point(size = 0.1) +
  scale_color_brewer(name = "DPT branch") +
  guides(colour = guide_legend(override.aes = list(size = 2), title.position = "top", title.hjust =  0.5))
pC <- ggplot(DM_dim_data, aes(x = DC1, y = DC2, colour = All_pT)) +
  theme_pubr(base_size = 7) +
  theme(axis.text = element_blank(), legend.key.width = unit(5, "pt"),
        legend.position = "right", legend.key.height = unit(10, "pt")) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(name = "Pseudotime") +
  guides(colour = guide_colourbar(title.position = "top", title.hjust =  0.5))

Fig_4ABC <- ggarrange(pA,NULL,pB,NULL,pC, ncol = 5, widths = c(1,0.1,1,0.1), align = "h",
                  labels = c("a", "", "b", "", "c"))

Fig_4ABC



#Supplementary Figure 4A####
SupplFig_4A <- 
  ggplot(DM_dim_data, aes(x = DC1, y = DC2, colour = Fibs_subpopulation)) +
  theme_pubr(base_size = 7) +
  theme(axis.text = element_blank(), legend.title = element_blank()) +
  geom_point(size = 0.1) +
  scale_color_manual(values = Fibs_col.palette) +
  facet_wrap(~Dataset, nrow = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
SupplFig_4A

#pt GENE CLUSTERING####
#A.M_pT (Alveolar to Myo)####
names(A.M_pT_res)
datExpr_Full <- list()
for(i in names(A.M_pT_res)){
  datExpr_Full[[i]]<- t(A.M_pT_res[[i]]$pT_fit)
}

sig_genes <- as.character(na.omit(rownames(A.M_pT_res_stats2)[A.M_pT_res_stats2$meta_p_adj.BH < 1e-10]))
length(sig_genes)

datExpr_CTR.TUMOUR <- list()
for(i in names(datExpr_Full)[c(1:4,7,8)]){
  datExpr_CTR.TUMOUR[[i]] <- datExpr_Full[[i]][, sig_genes]
}
datExpr_CTR.TUMOUR_median <- 
  Reduce(pmedian, datExpr_CTR.TUMOUR)

fit <- hclust(as.dist(1-cor(datExpr_CTR.TUMOUR_median)), method="ward.D2")
plot(fit) # display dendogram
groups <- cutree(fit, k=4)

pheatmap::pheatmap(t(datExpr_CTR.TUMOUR_median[, order(groups)]), scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,
                   clustering_distance_rows = "correlation",
                   clustering_method = "ward.D2",
                   #cutree_rows = 4,
                   show_rownames = F,
                   annotation_row = data.frame(
                     row.names = sig_genes,
                     cluster = as.factor(groups)
                   ))
rowMatch <- match(rownames(A.M_pT_res_stats2) ,sig_genes)
A.M_pT_res_stats2$Cluster <- as.factor(groups)[rowMatch]

A.M_pT_res_stats2$Module <- factor(
  A.M_pT_res_stats2$Cluster,
  levels = c(1:4),
  labels = c("Progenitor","Differentiation","Transient","Proto.Differentiation")
)
levels(A.M_pT_res_stats2$Module)
A.M_pT_res_stats2$Module <- factor(A.M_pT_res_stats2$Module, levels = c("Progenitor", "Transient", "Proto.Differentiation", "Differentiation"))

gene_order <- factor(
  groups,
  levels = c(1:4),
  labels = c("Progenitor","Differentiation","Transient","Proto.Differentiation")
)
gene_order <- factor(gene_order, levels = c("Progenitor", "Transient", "Proto.Differentiation", "Differentiation"))

Mod_cols <- RColorBrewer::brewer.pal(5, "Set1")
names(Mod_cols) <- levels(A.M_pT_res_stats2$Module)

#Figure 4D####
AM.PT_TUMOUR_Heatmap <- 
  pheatmap::pheatmap(t(datExpr_CTR.TUMOUR_median)[order(gene_order), ], scale = "row",
                     cluster_cols = F, cluster_rows = F,
                     #labels_row = rownames_filtered,
                     #clustering_distance_rows = "euclidean",
                     #clustering_method = "ward.D2",
                     #cutree_rows = 4,
                     breaks = seq(-3,3,length=101),
                     show_rownames = F, show_colnames = F,
                     legend = F,
                     annotation_legend = F,
                     fontsize = 7,
                     annotation_names_row = F,
                     annotation_names_col = F,
                     gaps_row = c(
                       sum(as.numeric(gene_order) == 1, na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:2), na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:3), na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:4), na.rm = T)
                     ),
                     annotation_row = data.frame(
                       row.names = A.M_pT_res_stats2$Gene[A.M_pT_res_stats2$Gene %in% sig_genes],
                       Module = A.M_pT_res_stats2$Module[A.M_pT_res_stats2$Gene %in% sig_genes]
                     ),
                     annotation_col = data.frame(
                       row.names = rownames(datExpr_CTR.TUMOUR_median),
                       Pseudotime = as.numeric(rownames(datExpr_CTR.TUMOUR_median))
                     ),
                     annotation_colors = list(
                       Pseudotime = viridis::viridis(3),
                       Module = Mod_cols
                     )
  )

AM.PT_TUMOUR_Heatmap

DC_plotting.df <-   data.frame(
  DM_dim_data,
  Fibs.integrated@meta.data
)
names(DC_plotting.df)
AM_pT_plot <-  
  DC_plotting.df %>%
  ggplot(aes(x = DC1, y = DC2, colour = A.M_pT)) +
  geom_point(size = 0.1) + 
  theme_void(base_size = 7) +
  theme(legend.position = c(0.1,1), legend.justification = c("left", "top"),
        legend.direction = "horizontal",
        legend.key.width = unit(8, "pt"), legend.key.height = unit(2, "pt")) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_color_viridis_c(name = "DPT", na.value = "grey80") 
AM_pT_plot

Fig_4D <- ggarrange(
  ggarrange(NULL, AM_pT_plot, ncol = 2),
  AM.PT_TUMOUR_Heatmap$gtable,
  nrow = 2, heights = c(2.5,6.5)
)

#T.M_pT (Adventitial to Myo)####
names(T.M_pT_res)
datExpr_Full <- list()
for(i in names(T.M_pT_res)){
  datExpr_Full[[i]]<- t(T.M_pT_res[[i]]$pT_fit)
}
sig_genes <- as.character(na.omit(rownames(T.M_pT_res_stats2)[T.M_pT_res_stats2$meta_p_adj.BH < 1e-10]))
length(sig_genes)

datExpr_CTR.TUMOUR <- list()
for(i in names(datExpr_Full)[c(1:4,7,8)]){
  datExpr_CTR.TUMOUR[[i]] <- datExpr_Full[[i]][, sig_genes]
}
datExpr_CTR.TUMOUR_median <- 
  Reduce(pmedian, datExpr_CTR.TUMOUR)

fit <- hclust(as.dist(1-cor(datExpr_CTR.TUMOUR_median)), method="ward.D2")
plot(fit) # display dendogram
groups <- cutree(fit, k=4)

pheatmap::pheatmap(t(datExpr_CTR.TUMOUR_median[, order(groups)]), scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,
                   clustering_distance_rows = "correlation",
                   clustering_method = "ward.D2",
                   #cutree_rows = 4,
                   show_rownames = F,
                   annotation_row = data.frame(
                     row.names = sig_genes,
                     cluster = as.factor(groups)
                   ))
rowMatch <- match(rownames(T.M_pT_res_stats2) ,sig_genes)
T.M_pT_res_stats2$Cluster <- as.factor(groups)[rowMatch]

T.M_pT_res_stats2$Module <- factor(
  T.M_pT_res_stats2$Cluster,
  levels = c(1:4),
  labels = c("Proto.Differentiation","Progenitor","Transient","Differentiation")
)
levels(T.M_pT_res_stats2$Module)
T.M_pT_res_stats2$Module <- factor(T.M_pT_res_stats2$Module, levels = c("Progenitor", "Transient", "Proto.Differentiation", "Differentiation"))

gene_order <- factor(
  groups, levels = c(1:4),
  labels = c("Proto.Differentiation","Progenitor","Transient","Differentiation")
)
gene_order <- factor(gene_order, levels = c("Progenitor", "Transient", "Proto.Differentiation", "Differentiation"))

Mod_cols <- RColorBrewer::brewer.pal(5, "Set1")
names(Mod_cols) <- levels(T.M_pT_res_stats2$Module)

Heatmap_df <- t(datExpr_CTR.TUMOUR_median)[order(gene_order), ]
Heatmap_df <- Heatmap_df[!is.na(gene_order[order(gene_order)]), ]
dim(Heatmap_df)
#Figure 4E####
TM.PT_TUMOUR_Heatmap <- 
  pheatmap::pheatmap(Heatmap_df, scale = "row",
                     cluster_cols = F, cluster_rows = F,
                     #labels_row = rownames_filtered,
                     #clustering_distance_rows = "euclidean",
                     #clustering_method = "ward.D2",
                     #cutree_rows = 4,
                     breaks = seq(-3,3,length=101),
                     show_rownames = F, show_colnames = F,
                     legend = F,
                     annotation_legend = F,
                     fontsize = 7,
                     annotation_names_row = F,
                     annotation_names_col = F,
                     gaps_row = c(
                       sum(as.numeric(gene_order) == 1, na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:2), na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:3), na.rm = T),
                       sum(as.numeric(gene_order) %in% c(1:4), na.rm = T)
                     ),
                     annotation_row = data.frame(
                       row.names = T.M_pT_res_stats2$Gene[T.M_pT_res_stats2$Gene %in% sig_genes],
                       Module = T.M_pT_res_stats2$Module[T.M_pT_res_stats2$Gene %in% sig_genes]
                     ),
                     annotation_col = data.frame(
                       row.names = rownames(datExpr_CTR.TUMOUR_median),
                       Pseudotime = as.numeric(rownames(datExpr_CTR.TUMOUR_median))
                     ),
                     annotation_colors = list(
                       Pseudotime = viridis::viridis(3),
                       Module = Mod_cols
                     )
  )
TM_pT_plot <-  
  DC_plotting.df %>%
  ggplot(aes(x = DC1, y = DC2, colour = T.M_pT)) +
  geom_point(size = 0.1) + 
  theme_void(base_size = 7) +
  theme(legend.position = c(0.1,1), legend.justification = c("left", "top"),
        legend.direction = "horizontal",
        legend.key.width = unit(8, "pt"), legend.key.height = unit(2, "pt")) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  scale_color_viridis_c(name = "DPT", na.value = "grey80") 
TM_pT_plot

Fig_4E <- ggarrange(
  ggarrange(NULL, TM_pT_plot, ncol = 2),
  TM.PT_TUMOUR_Heatmap$gtable,
  nrow = 2, heights = c(2.5,6.5)
)
Fig_4E

#Core Module signatures in All_pT####
Mod.Overlap_df <-  cbind(A.M_pT_res_stats2[,c("Gene", "Module", "meta_p_adj.BH")],
                         T.M_pT_res_stats2[,c("Gene", "Module", "meta_p_adj.BH")])
names(Mod.Overlap_df) <- c("GeneA", "Alv2M_Module", "Alv2M_adjP", "GeneB", "Adv2M_Module", "Adv2M_adjP")

Core.Mods <- Mod.Overlap_df[Mod.Overlap_df$Alv2M_Module == Mod.Overlap_df$Adv2M_Module, ]
Module.Genes_list <- list()
for(i in levels(Core.Mods$Alv2M_Module)) {
  Module.Genes_list[[i]] <- as.character(na.omit(Core.Mods$GeneA[Core.Mods$Alv2M_Module == i]))
}

names(Fibs.integrated@meta.data)
DefaultAssay(Fibs.integrated) <- "RNA"
Fibs.integrated <- 
  AddModuleScore(Fibs.integrated,
                 features = Module.Genes_list)
names(Fibs.integrated@meta.data)
names(Fibs.integrated@meta.data)[grep("^Cluster", names(Fibs.integrated@meta.data))] <-
  paste0("Core_", levels(T.M_pT_res_stats2$Module))

pT_df <- data.frame(
  Fibs.integrated@meta.data[, c("Sample.type", "All_pT", "dpt_branch", "Fibs_MajorClusters", "Dataset", "Sample.Subtype",
                                paste("Core", levels(T.M_pT_res_stats2$Module), sep = "_"))]
  
)
variables = paste("Core", levels(T.M_pT_res_stats2$Module), sep = "_")
pT_df_long <- reshape2::melt(pT_df[, c("All_pT", "Sample.type", "Dataset",variables)],
                             id.vars = c("All_pT","Sample.type", "Dataset"), na.rm = T)
levels(pT_df_long$variable)
pT_df_long$Dataset_ordered <- 
  factor(pT_df_long$Dataset, levels = levels(pT_df_long$Dataset) )
pT_df_long$Dataset_type <- 
  factor(pT_df_long$Dataset_ordered, levels = levels(pT_df_long$Dataset_ordered),
         labels = c("CTR & Tumour", "CTR Only", "CTR Only", rep("CTR & Tumour", 6)))

names(pT_df_long)
levels(pT_df_long$variable)
Mod_cols2 <- Mod_cols
names(Mod_cols2) <- levels(pT_df_long$variable)
pT_df_long$Mod.label <- factor(pT_df_long$variable,
                               labels = c("Progenitor", "Early-activation",
                                          "Proto-differentiation", "Differentiation"))
#Supplementary Figure 4B - Faceted Core Module Loess plot####
SupplFig_4B <- 
  pT_df_long[!is.na(pT_df_long$variable), ] %>%
  ggplot(aes(x = All_pT, y = value, colour = variable)) +
  geom_point(size = 0.1, aes(colour = variable), alpha = 0.1) + 
  facet_grid(Mod.label~Dataset_type+Dataset, scales = "free_y") + 
  geom_smooth(method = "loess") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank(),
        #axis.text = element_blank(), axis.ticks = element_blank(),
        #strip.text.y = element_blank(),
        strip.text = element_text(size = 6)) +
  xlab("Pseudotime") + ylab("Module Score") +
  scale_color_manual(values = Mod_cols2)
SupplFig_4B


#Figure 4F####
Fig_4F <- 
  pT_df_long[!is.na(pT_df_long$variable), ] %>%
  ggplot(aes(x = All_pT, y = value, colour = variable)) +
  #geom_point(size = 0.1, aes(colour = variable), alpha = 0.1) + 
  facet_grid(Mod.label~Dataset_type, scales = "free_y") + 
  geom_smooth(method = "loess") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank(),
        #axis.text = element_blank(), axis.ticks = element_blank(),
        #strip.text.y = element_blank(),
        strip.text = element_text(size = 6)) +
  xlab("Pseudotime") + ylab("Module Score") +
  scale_color_manual(values = Mod_cols2)
Fig_4F


#Figure 4G####
Fibs.integrated$All_pT_5cat <- factor(Hmisc::cut2(Fibs.integrated$All_pT, g = 5),
                                      labels = paste0("pT q", 1:5))

test <- aggregate(formula = Core_Transient ~ All_pT_5cat + SampleID + Sample.type,
                  data = Fibs.integrated@meta.data,
                  FUN = median)
Transient.boxplot <- 
  test %>%
  drop_na(All_pT_5cat) %>%
  ggplot(aes(x = Sample.type, y = Core_Transient, fill = Sample.type)) +
  theme_pubr(base_size = 7) +
  theme(axis.text.x = element_blank()) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_grid(~All_pT_5cat) +
  ylab("Early-activation Module\n(Normalised exprs)") +
  stat_compare_means(comparisons = list(c("Control", "Tumour")),
                     #label = "p.signif",
                     size = 2.5,
                     bracket.size = 0.1, tip.length = 0.01,
                     label.y = max(test$Core_Transient)+0.1,
                     #vjust = 0.5,
                     hide.ns = T)+
  expand_limits(y = 2)

test <- aggregate(formula = Core_Proto.Differentiation ~ All_pT_5cat + SampleID + Sample.type,
                  data = Fibs.integrated@meta.data,
                  FUN = median)
PD.boxplot <- 
  test %>%
  drop_na(All_pT_5cat) %>%
  ggplot(aes(x = Sample.type, y = Core_Proto.Differentiation, fill = Sample.type)) +
  theme_pubr(base_size = 7) +
  theme(axis.text.x = element_blank()) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_grid(~All_pT_5cat) +
  ylab("Proto-Differentation Module\n(Normalised exprs)") +
  stat_compare_means(comparisons = list(c("Control", "Tumour")),
                     #label = "p.signif",
                     size = 2.5,
                     bracket.size = 0.1, tip.length = 0.01,
                     label.y = max(test$Core_Proto.Differentiation)+0.1,
                     #vjust = 0.5,
                     hide.ns = T)+
  expand_limits(y = 2)

test <- aggregate(formula = Core_Differentiation ~ All_pT_5cat + SampleID + Sample.type,
                  data = Fibs.integrated@meta.data,
                  FUN = median)
Differentiation.boxplot <- 
  test %>%
  drop_na(All_pT_5cat) %>%
  ggplot(aes(x = Sample.type, y = Core_Differentiation, fill = Sample.type)) +
  theme_pubr(base_size = 7) +
  theme(axis.text.x = element_blank()) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_grid(~All_pT_5cat) +
  ylab("Differentation Module\n(Normalised exprs)") +
  stat_compare_means(comparisons = list(c("Control", "Tumour")),
                     #label = "p.signif",
                     size = 2.5,
                     bracket.size = 0.1, tip.length = 0.01,
                     label.y = max(test$Core_Differentiation)+0.1,
                     #vjust = 0.5,
                     hide.ns = T)+
  expand_limits(y = 1.5)

Fig_4G <- 
  ggarrange(
    Transient.boxplot,
    PD.boxplot,
    Differentiation.boxplot,
    nrow = 3, common.legend = T, legend = "right")
Fig_4G

#Supplementary Figure 4C####
nMods = length(levels(Mod.Overlap_df$Alv2M_Module))
Mods <- levels(Mod.Overlap_df$Alv2M_Module)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nMods, ncol = nMods);
CountTbl = matrix(0, nrow = nMods, ncol = nMods);
# Execute all pairwaise comparisons
for (Alv.mod in 1:nMods)
  for (Adv.mod in 1:nMods)
  {
    AlvMembers = (Mod.Overlap_df$Alv2M_Module == Mods[Alv.mod])
    AdvMembers = (Mod.Overlap_df$Adv2M_Module == Mods[Adv.mod])
    pTable[Alv.mod, Adv.mod] = fisher.test(AlvMembers, AdvMembers, alternative = "greater")$p.value
    CountTbl[Alv.mod, Adv.mod] = sum(Mod.Overlap_df$Alv2M_Module == Mods[Alv.mod] &
                                       Mod.Overlap_df$Adv2M_Module == Mods[Adv.mod], na.rm = T)
  }
diag(pTable)
max(diag(pTable))

pTable.stars <- gtools::stars.pval(pTable)
pTable.stars[pTable.stars == " "] <- "ns"

text.mtx <- paste0(CountTbl, "\n(p=", signif(pTable,3), ")")

AlvModTotals = apply(CountTbl, 1, sum)
AdvModTotals = apply(CountTbl, 2, sum)
Mod.labels <- c("Progenitor", "Early-activation", "Proto-differentiation", "Differentiation")
colnames(CountTbl) <- Mods
rownames(CountTbl) <- Mods
Fisher.hm <- reshape2::melt(CountTbl)
names(Fisher.hm) <- c("Alv.to.Myo", "Adv.to.Myo", "Count")
Fisher.hm$p.value <- reshape2::melt(pTable)[,3]
Fisher.hm$p.stars <- gtools::stars.pval(Fisher.hm$p.value)
Fisher.hm$label <- paste0(Fisher.hm$Count, "\n", Fisher.hm$p.stars)

SupplFig_4C <- 
  Fisher.hm %>%
  ggplot(aes(x = factor(Alv.to.Myo, labels = c("Progenitor", "Early-activation", "Proto-differentiation", "Differentiation")),
             y = factor(Adv.to.Myo, labels = c("Progenitor", "Early-activation", "Proto-differentiation", "Differentiation")),
             label = label, fill = -log10(p.value))) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(7,"pt"), plot.title = element_text(hjust = 0.5)) +
  geom_tile() +
  geom_text(size = 2) +
  rotate_x_text(angle = 45) +
  rotate_y_text(angle = 45) +
  scale_fill_gradient(low = "white", high = "red", name = "-log10(p)") +
  ggtitle("Module overlap between trajectories") +
  ylab("Adventitial to Myo") +
  xlab("Alveolar to Myo")
SupplFig_4C

#Supplementary Figure 4D####
SupplFig_4D <- 
  DM_dim_data %>%
  drop_na(All_pT) %>%
  ggplot(aes(x = DC1, y = DC2, colour = factor(Hmisc::cut2(All_pT, g = 5),
                                               labels = paste0("pT q", 1:5)))) +
  theme_pubr(base_size = 7) +
  theme(legend.key.size = unit(1,"pt"),
        legend.position = c(0.01,0.01),
        legend.justification = c("left","bottom")) +
  geom_point(size = 0.1) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_d3(name = "Pseudotime bin")
SupplFig_4D

#Figure 4I - MxIHC HSPA1A####
Fibs_HSPA1A.agg <- aggregate(
  MxIHC_Fibroblasts$HSPA1A,
  by = list(Sample = MxIHC_Fibroblasts$Slide.name,
            Subpop = MxIHC_Fibroblasts$Cell.type,
            TvN = MxIHC_Fibroblasts$Sample.Region.TvN,
            Subtype = MxIHC_Fibroblasts$Subtype,
            Inflam.Fibro = MxIHC_Fibroblasts$Control...Fibrosis.Inflammation,
            exclusion = MxIHC_Fibroblasts$Include.Normal.95pct),
  median)

Fig_4I <- 
  Fibs_HSPA1A.agg %>%
  drop_na(TvN) %>%
  filter(!Sample == "Drop-Seq Cohort 1") %>%
  ggplot(aes(x = TvN, y = x)) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  facet_grid(~Subpop) +
  #geom_violin(scale = "width", aes(fill = Subpop)) +
  #geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(size = 0.1, aes(colour = Subpop)) +
  geom_line(aes(group = Sample, colour = Subpop)) +
  scale_colour_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  ylab("HSPA1A\n(MxIHC Median per cell exprs)") +
  stat_compare_means(comparisons = list(c("Control", "Tumour")), paired = T,
                     #label = "p.signif",
                     size = 2.5,
                     method = "wilcox", label.y = 0.18, tip.length = 0.01)+
  expand_limits(y = 0.2)
Fig_4I

