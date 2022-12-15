#Figure 2####
#Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
#Hanley et al
library(Seurat)
library(WGCNA)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(msigdbr)
library(GSVA)
library(limma)
library(readxl)

data_directory <- " " #Specify directory where Zenodo repo is saved 
source(paste0(data_directory, "0_NewFunctions.R"))

load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))
load(paste0(data_directory, "HCLA_annotation.Rdata"))
Fibroblast_GeneSets_v2_2021_0623 <- read_excel(paste0(data_directory, "Fibroblast GeneSets v2 2021_0623.xlsx"))

#Fibroblast clustering####
DefaultAssay(Fibs.integrated) <- "integrated"
Fibs.integrated <- ScaleData(Fibs.integrated)#, features = Non.ribo.genes)
Fibs.integrated <- RunPCA(Fibs.integrated)#, features = Non.ribo.genes)
#Jackstraw calculation #'d out to save time...
#Fibs.integrated <- JackStraw(Fibs.integrated, dims = 50)
#Fibs.integrated <- ScoreJackStraw(Fibs.integrated, dims = 1:50)
#table(Fibs.integrated@reductions$pca@jackstraw$overall.p.values[,2] < 0.01)
ElbowPlot(Fibs.integrated, ndims = 50)

PCs2use <- c(1:30)
Fibs.integrated <- RunUMAP(Fibs.integrated, dims = PCs2use)
Fibs.integrated <- FindNeighbors(Fibs.integrated, dims = PCs2use, verbose = FALSE)
Fibs.integrated <- AddMetaData(Fibs.integrated, as.data.frame(Fibs.integrated@reductions$umap@cell.embeddings))
Fibs.integrated <-  FindClusters(Fibs.integrated, resolution = 0.1, verbose = FALSE)
Fibs.integrated$Fibs_MajorClusters <- factor(Fibs.integrated$seurat_clusters, levels = c(2,1,0), labels = c("Adventitial", "Alveolar", "Myo"))

#Figure 2B - UMAP plot####
Fig_2B <- 
  DimPlot(Fibs.integrated, group.by = "Fibs_MajorClusters", pt.size = 0.1) &
  theme_pubr(base_size = 5) &
  scale_colour_manual(values = Fibs_col.palette) &
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        plot.title = element_blank(),
        legend.key.size = unit(2,"pt")
        )
Fig_2B


#Fibroblast marker analysis (single-cell)#### 
Fibs.integrated <- SetIdent(Fibs.integrated, value = "Fibs_MajorClusters")
All.Markers <- FindAllMarkers(Fibs.integrated, test.use = "bimod", assay = "RNA", only.pos = T)
library(org.Hs.eg.db)
library(hpar)
All.Markers$ENSEMBL <- mapIds(org.Hs.eg.db, keys = All.Markers$gene, keytype = "SYMBOL", column="ENSEMBL")
test <- hpar::getHpa(All.Markers$ENSEMBL)
table(test$Reliability)
test <- test[!duplicated(test$Gene.name), ]
rownames(test) <- test$Gene.name
All.Markers$HPA_reliability <- test[All.Markers$gene, "Reliability"]
test <- hpar::getHpa(All.Markers$ENSEMBL, hpadata = "hpaSubcellularLoc")
test <- test[!duplicated(test$Gene.name), ]
rownames(test) <- test$Gene.name
test$membrane <- 1:nrow(test) %in% grep("membrane", test$GO.id, fixed = T)
test$cytosol <- 1:nrow(test) %in% grep("Cytosol", test$GO.id, fixed = T)
test$filaments <- 1:nrow(test) %in% grep("filaments", test$GO.id, fixed = T)
test$nuclear <- 1:nrow(test) %in% grep("Nucle", test$GO.id, fixed = T)
All.Markers$subcell.loc <- test[All.Markers$gene,
                                c("membrane","cytosol","filaments","nuclear")]
All.Markers$pct.diff <- All.Markers$pct.1 - All.Markers$pct.2
All.Markers$Detected.Intracellular <- rowSums(All.Markers$subcell.loc) > 0

#Fibroblast Marker analysis (Sample-level)####
Fibs.integrated <- SetIdent(Fibs.integrated, value = "Fibs_MajorClusters")
Fibs.integrated_samples <- AverageExpression(Fibs.integrated, assays = "RNA", return.seurat = T, 
                                             group.by = c("Fibs_MajorClusters", "SampleID2"),
                                             slot = "data")
Fibs.integrated_samples <- NormalizeData(Fibs.integrated_samples)
dim(Fibs.integrated_samples)

#extract meta data from sample names
Fibs.integrated_samples@meta.data <- 
  Fibs.integrated_samples@meta.data %>%
  mutate(Fibs_MajorClusters = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,1],
         Dataset = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,2],
         SampleID = str_split_fixed(colnames(Fibs.integrated_samples), "_", 3)[,3])

#Collate meta data from sc dataset
Fibs.Sample_MetaData <- Fibs.integrated_samples@meta.data[,"SampleID", drop = F]
matchSamples <- match(Fibs.Sample_MetaData$SampleID, Sample_MetaData$SampleID)
table(is.na(matchSamples))
Fibs.Sample_MetaData <- cbind(Fibs.Sample_MetaData, Sample_MetaData[matchSamples, ])
keepCols_metaData <- c("SampleID", "nCount_RNA", "nFeature_RNA",
                       "Fibs_MajorClusters", "Dataset")
Fibs.integrated_samples@meta.data <- Fibs.integrated_samples@meta.data[, keepCols_metaData]
Fibs.integrated_samples <- AddMetaData(Fibs.integrated_samples, Fibs.Sample_MetaData)
Fibs.integrated_samples$Sample.type[Fibs.integrated_samples$Sample.type == "Tumor"] <- "Tumour"
Fibs.integrated_samples$pN <- droplevels(Fibs.integrated_samples$pN, exclude = c(NA, "NA"))
levels(Fibs.integrated_samples$pN)

#Run DE analysis
Fibs.integrated_samples <- SetIdent(Fibs.integrated_samples, value = "Fibs_MajorClusters")
Sample.level_Markers <- FindAllMarkers(Fibs.integrated_samples, test.use = "wilcox", assay = "RNA",
                                       only.pos = F, return.thresh = 1, min.pct = 0.25, logfc.threshold = 0)

#Collate sc and sample DE analysis results
Marker.match <- match(paste(All.Markers$cluster, All.Markers$gene),
                      paste(Sample.level_Markers$cluster, Sample.level_Markers$gene))
All.Markers$p_val_adj.Sample <- Sample.level_Markers$p_val_adj[Marker.match]
All.Markers$log2FC.Sample <- Sample.level_Markers$avg_log2FC[Marker.match]

#N.B this All.Markers object is pre-saved in Fibs.integrated Rdata file for use in later figures

#Figure 2C Heatmap####
All.Markers$cluster <- factor(as.character(All.Markers$cluster), levels = c("Adventitial", "Alveolar", "Myo"))
All.Markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0) %>%
  top_n(w = -p_val_adj.Sample, n = 10) -> top10_markers
DefaultAssay(Fibs.integrated_samples) <- "RNA"
Fibs.integrated_samples <- ScaleData(Fibs.integrated_samples, features = top10_markers$gene, assay = "RNA")
levels(Fibs.integrated_samples$Fibs_MajorClusters)
Fibs.integrated_samples$Fibs_MajorClusters <- factor(Fibs.integrated_samples$Fibs_MajorClusters,
  levels = c("Adventitial", "Alveolar", "Myo")
)
Fig_2C <-
  DoHeatmap(Fibs.integrated_samples, features = rev(top10_markers$gene),
            assay = "RNA",
            group.by = c("Fibs_MajorClusters") ,
            group.colors = Fibs_col.palette,
            label = F, disp.max = 2) &
  theme(axis.text.y = element_text(face = "italic", size = 5),legend.position = "none") 
Fig_2C


#Supplementary Data 2 - export####
View(All.Markers)
Table_S2 <- All.Markers[,-10]
Table_S2 <- Table_S2[,-c(14:20)]
Table_S2$Gene.Symbol <- paste0("'", Table_S2$gene)
Table_S2 <- Table_S2[!is.na(Table_S2$ENSEMBL), ]

  
collections <- msigdbr_collections()
#REACTOME Pathway analysis####
h_gene_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
#prepare in list format for fgsea
Reactome_pathwaysH = split(x = h_gene_sets$human_gene_symbol, f = h_gene_sets$gs_name)

REACTOME_gsva_res <- gsva(as.matrix(Fibs.integrated_samples@assays$RNA@data), gset.idx.list = Reactome_pathwaysH, method="gsva",
                          min.sz=10, max.sz=500)

f <- factor(Fibs.integrated_samples$Fibs_MajorClusters, levels = c("Myo","Alveolar", "Adventitial"))
design <- model.matrix(~0+f)
colnames(design) <- c("Myo", "Alveolar", "Adventitial")

REACTOME_gsva.fit <- lmFit(REACTOME_gsva_res, design)
REACTOME_gsva.fit <- eBayes(REACTOME_gsva.fit)
contrast.matrix <- makeContrasts(Myo.v.ALL = Myo-(Alveolar+Adventitial)/2,
                                 Alv.v.ALL = Alveolar-(Myo+Adventitial)/2,
                                 Adv.v.ALL = Adventitial-(Myo+Alveolar)/2,
                                 Myo.v.Alv = Myo - Alveolar,
                                 Myo.v.Adv = Myo-Adventitial,
                                 Alv.v.Adv = Alveolar-Adventitial,
                                 levels=design)
REACTOME_gsva.fit2 <- contrasts.fit(REACTOME_gsva.fit, contrast.matrix)
REACTOME_gsva.fit2 <- eBayes(REACTOME_gsva.fit2)
res2 <- decideTests(REACTOME_gsva.fit2, p.value=0.01, method = "global")
summary(res2)

REACTOME_tt.res <- list(
  Myo = topTable(REACTOME_gsva.fit2, coef = "Myo.v.ALL", n=1270),
  Alveolar = topTable(REACTOME_gsva.fit2, coef = "Alv.v.ALL", n=1270),
  Adventitial = topTable(REACTOME_gsva.fit2, coef = "Adv.v.ALL", n=1270)
)

REACTOME_tt.res_df <- do.call(rbind, REACTOME_tt.res)
REACTOME_tt.res_df$Cluster <- do.call(rbind, strsplit(rownames(REACTOME_tt.res_df), ".", 2))[,1]
REACTOME_tt.res_df$Pathway <- do.call(rbind, strsplit(rownames(REACTOME_tt.res_df), ".", 2))[,2]
REACTOME_tt.res_df$Pathway.label <- gsub("REACTOME_", "", REACTOME_tt.res_df$Pathway)
REACTOME_tt.res_df$Pathway.label <- gsub("_", " ", REACTOME_tt.res_df$Pathway.label)

REACTOME_tt.res_subset <- REACTOME_tt.res_df[REACTOME_tt.res_df$logFC > 0, ] %>%
  group_by(Cluster) %>%
  top_n(wt = -adj.P.Val, n = 5)
REACTOME_tt.res_subset <- REACTOME_tt.res_subset[-10, ]
REACTOME_tt.res_subset <- rbind(REACTOME_tt.res_subset, REACTOME_tt.res_df["Alveolar.REACTOME_TRP_CHANNELS", ])

Fig_2D <- 
  REACTOME_tt.res_subset %>%
  ggplot(aes(x = reorder(Pathway.label, logFC), y = logFC, label =  paste0("adj.P=", signif(adj.P.Val, 3)),
             fill = Cluster)) +
  theme_pubr(base_size = 5)+
  theme(axis.title.y = element_blank(),legend.title = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.position = "top", strip.text = element_blank(), legend.key.size = unit(10, "pt"))+
  geom_bar(stat = "identity", alpha = 0.75)+
  facet_grid(Cluster~., scales = "free", space = "free") +
  geom_text(hjust = 1.1, colour = "white", size = 2) +
  scale_fill_manual(values = Fibs_col.palette) +
  coord_flip()
Fig_2D

#Supplementary Data 3 export#####
Table_S3 <- REACTOME_tt.res_df

#Matrisome Gene analysis####
CGP_gene_sets = msigdbr(species = "human", category = "C2")
CGP_pathwaysH = split(x = CGP_gene_sets$human_gene_symbol, f = CGP_gene_sets$gs_name)
All.Markers$NABA_Matrisome <- All.Markers$gene %in% CGP_pathwaysH$NABA_MATRISOME
All.Markers$NABA_BASEMENT_MEMBRANES <- All.Markers$gene %in% CGP_pathwaysH$NABA_BASEMENT_MEMBRANES
All.Markers$NABA_COLLAGENS <- All.Markers$gene %in% CGP_pathwaysH$NABA_COLLAGENS
All.Markers$NABA_ECM_GLYCOPROTEINS <- All.Markers$gene %in% CGP_pathwaysH$NABA_ECM_GLYCOPROTEINS
All.Markers$NABA_PROTEOGLYCANS <- All.Markers$gene %in% CGP_pathwaysH$NABA_PROTEOGLYCANS
names(All.Markers)
All.Markers$type <- factor(
  apply(All.Markers[,c("NABA_BASEMENT_MEMBRANES","NABA_COLLAGENS","NABA_ECM_GLYCOPROTEINS","NABA_PROTEOGLYCANS")],1, which.max),
                           levels = 1:5,
                           labels = c("Basement.Membranes", "Interstitial.Collagens",
                                      "ECM.Glycoproteins", "Proteoglycans", "Other")
)
All.Markers$type[All.Markers$type == "Basement.Membranes" &
                   All.Markers$NABA_BASEMENT_MEMBRANES == F] <- "Other"
#Figure 2E####
Fig_2E <- 
  All.Markers[!All.Markers$type == "Other" &
                All.Markers$p_val_adj.Sample < 0.01, ] %>%
  ggplot(aes(x = type, fill = cluster)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "top", legend.key.size = unit(10,"pt"), legend.title = element_blank()) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 65) +
  ylab("Proportion of DEGs") +
  xlab("Matrisome Category")
Fig_2E


Fibs.integrated_samples <- AddModuleScore(
  Fibs.integrated_samples,
  features = list(
    CGP_pathwaysH$NABA_BASEMENT_MEMBRANES,
    CGP_pathwaysH$NABA_COLLAGENS,
    CGP_pathwaysH$NABA_ECM_GLYCOPROTEINS,
    CGP_pathwaysH$NABA_PROTEOGLYCANS,
    CGP_pathwaysH$NABA_COLLAGENS[!CGP_pathwaysH$NABA_COLLAGENS %in% CGP_pathwaysH$NABA_BASEMENT_MEMBRANES]
  ),
  assay = "RNA"
)
names(Fibs.integrated_samples@meta.data)[grep("^Cluster", names(Fibs.integrated_samples@meta.data))] <-
  c("Basement.Membranes", "Collagens", "ECM.Glycoproteins", "Proteoglycans", "Interstitial.Collagens")

#Figure 2FG####
comparisons <- list(c("Myo", "Adventitial"), c("Myo", "Alveolar"), c("Alveolar", "Adventitial"))
pA <- 
  Fibs.integrated_samples@meta.data %>%
  ggplot(aes(x = Fibs_MajorClusters, y = Basement.Membranes,
             fill = Fibs_MajorClusters)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = list(c("Myo", "Adventitial"), c("Alveolar", "Adventitial")),
                     #label = "p.signif",
                     size = 2.5,
                     hide.ns = T,label.y = c(0.65,0.8)) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  ylab("Basement Membranes") +
  guides(fill = "none") +  expand_limits(y = 0.9)

pB <- 
  Fibs.integrated_samples@meta.data %>%
  ggplot(aes(x = Fibs_MajorClusters, y = Interstitial.Collagens,
             fill = Fibs_MajorClusters)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", axis.title.x = element_blank()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = list(c("Myo", "Adventitial"), c("Myo", "Alveolar")),
                     #label = "p.signif",
                     size = 2.5,
                     hide.ns = T, label.y = c(0.9,1.1)) +
  scale_fill_manual(values = Fibs_col.palette) +
  ylab("Interstitial Collagens") +
  rotate_x_text(angle = 45) +
  guides(fill = "none") + expand_limits(y = 1.25)

Fig_2FG <- ggarrange(pA,pB, ncol = 2, labels = c("f", "g"))
Fig_2FG

#Figure 2H####
comparisons <- list(c("Control", "Tumour"))
Fig_2H <- 
  Fibs.integrated_samples@meta.data %>%
  drop_na(Sample.type) %>%
  ggplot(aes(x = Sample.type, y = Interstitial.Collagens)) +
  theme_pubr(base_size = 5) +
  #theme(legend.position = "right", axis.title.x = element_blank()) +
  facet_wrap(~Fibs_MajorClusters) +
  geom_boxplot(aes(fill = Fibs_MajorClusters), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons,
                     #label = "p.signif",
                     size = 2.5,
                     label.y = 0.9) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  ylab("Interstitial Collagens") +
  guides(fill = "none") +
  #ylim(c(-0.25,1.1)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") + expand_limits(y = 1.1)
Fig_2H

#Fibroblast Gene Signature analysis####
Fibs_GeneSets_filtered <- Fibroblast_GeneSets_v2_2021_0623[-1, ]
Fibs_GeneSets_list <- list()
for(i in names(Fibs_GeneSets_filtered)){
  Fibs_GeneSets_list[[i]] <- as.character(na.omit(Fibs_GeneSets_filtered[[i]]))
}

FibsSigs_gsva_res <- gsva(as.matrix(Fibs.integrated_samples@assays$RNA@data), gset.idx.list = Fibs_GeneSets_list, method="gsva")
Fibs.integrated_samples <- AddMetaData(Fibs.integrated_samples,
                                       metadata = as.data.frame(t(FibsSigs_gsva_res)))
f <- factor(Fibs.integrated_samples$Fibs_MajorClusters, levels = c("Myo","Alveolar", "Adventitial"))
design <- model.matrix(~0+f)
colnames(design) <- c("Myo", "Alveolar", "Adventitial")

FibsSigs_gsva.fit <- lmFit(FibsSigs_gsva_res, design)
FibsSigs_gsva.fit <- eBayes(FibsSigs_gsva.fit)
contrast.matrix <- makeContrasts(Myo.v.ALL = Myo-(Alveolar+Adventitial)/2,
                                 Alv.v.ALL = Alveolar-(Myo+Adventitial)/2,
                                 Adv.v.ALL = Adventitial-(Myo+Alveolar)/2,
                                 Myo.v.Alv = Myo - Alveolar,
                                 Myo.v.Adv = Myo-Adventitial,
                                 Alv.v.Adv = Alveolar-Adventitial,
                                 levels=design)
FibsSigs_gsva.fit2 <- contrasts.fit(FibsSigs_gsva.fit, contrast.matrix)
FibsSigs_gsva.fit2 <- eBayes(FibsSigs_gsva.fit2)
res2 <- decideTests(FibsSigs_gsva.fit2, p.value=0.01, method = "global")
summary(res2)

FibsSigs_tt.res <- list(
  Myo = topTable(FibsSigs_gsva.fit2, coef = "Myo.v.ALL", n=74),
  Alveolar = topTable(FibsSigs_gsva.fit2, coef = "Alv.v.ALL", n=74),
  Adventitial = topTable(FibsSigs_gsva.fit2, coef = "Adv.v.ALL", n=74)
)

FibsSigs_tt.res_df <- do.call(rbind, FibsSigs_tt.res)
FibsSigs_tt.res_df$Cluster <- do.call(rbind, strsplit(rownames(FibsSigs_tt.res_df), ".", 2))[,1]
FibsSigs_tt.res_df$Pathway <- gsub("Myo.", "", rownames(FibsSigs_tt.res_df), fixed = T)
FibsSigs_tt.res_df$Pathway <- gsub("Alveolar.", "", FibsSigs_tt.res_df$Pathway, fixed = T)
FibsSigs_tt.res_df$Pathway <- gsub("Adventitial.", "", FibsSigs_tt.res_df$Pathway, fixed = T)
FibsSigs_tt.res_df$Pathway.label <- gsub("_", " ", FibsSigs_tt.res_df$Pathway)

#Supplementary Data 4####
Table_S4 <- FibsSigs_tt.res_df

#Figure 2I####
Fig_2I <- 
  Fibs.integrated_samples@meta.data %>%
  drop_na(Sample.type) %>%
  ggplot(aes(x = Sample.type, y = Elyada_MyoCAFs_PDAC)) +
  theme_pubr(base_size = 5) +
  #theme(legend.position = "right", axis.title.x = element_blank()) +
  facet_wrap(~Fibs_MajorClusters) +
  geom_boxplot(aes(fill = Fibs_MajorClusters), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons,
                     #label = "p.signif",
                     size = 2.5,
                     label.y = 0.9) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  ylab("myoCAF Signature") +
  guides(fill = "none") +
  #ylim(c(-0.25,1.1)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") + expand_limits(y = 1.1)
Fig_2I

#Figure 2J####
comparisons <- list(c("Control", "Tumour"))
Fig_2J <- 
  Fibs.integrated_samples@meta.data %>%
  drop_na(Sample.type) %>%
  ggplot(aes(x = Sample.type, y = Elyada_iCAFs_PDAC)) +
  theme_pubr(base_size = 5) +
  #theme(legend.position = "right", axis.title.x = element_blank()) +
  facet_wrap(~Fibs_MajorClusters) +
  geom_boxplot(aes(fill = Fibs_MajorClusters), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons,
                     #label = "p.signif",
                     size = 2.5,
                     label.y = 0.8) +
  scale_fill_manual(values = Fibs_col.palette) +
  rotate_x_text(angle = 45) +
  ylab("iCAF Signature") +
  guides(fill = "none") +
  #ylim(c(-0.25,1.1)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") + expand_limits(y = 1)
Fig_2J


#Supplementary Figure 2A####
DefaultAssay(Fibs.integrated) <- "RNA"
Fibs.integrated <- NormalizeData(Fibs.integrated)
Fibs.integrated <- FindVariableFeatures(Fibs.integrated)
Fibs.integrated <- ScaleData(Fibs.integrated, verbose = FALSE)
Fibs.integrated <- RunPCA(Fibs.integrated, npcs = 30, verbose = FALSE)
Fibs.integrated <- FindNeighbors(Fibs.integrated, reduction = "pca", dims = 1:30)

Batchtest_raw <- Batch.Quant(seurat_obj = Fibs.integrated,
                             vars = c("Dataset"),
                             nPerm = 1000, assay_nn = "RNA")

Batchtest_cca <- Batch.Quant(seurat_obj = Fibs.integrated,
                              vars = c("Dataset"),
                              nPerm = 1000, assay_nn = "integrated")

Batch.test_res <- data.frame(
  raw = Batchtest_raw$Dataset,
  CCA  = Batchtest_cca$Dataset
)

SupplFig_2A <- 
  reshape2::melt(Batch.test_res) %>%
  ggplot(aes(x = variable, y = value)) +
  theme_pubr(base_size = 5) +
  #geom_violin
  geom_boxplot(outlier.size = 0.1) +
  #scale_y_log10() + 
  geom_hline(yintercept = 1.96, linetype = "dashed") +
  rotate_x_text(angle = 45) +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("Raw Data", "CCA Integrated")) +
  ylab("Intra-Dataset nn overlap\n(Z-score)") +
  ggforce::facet_zoom(ylim = c(-1, 12))
SupplFig_2A



#Supplementary Figure 2B####
DefaultAssay(Fibs.integrated) <- "integrated"
Fibs.integrated <-
  FindClusters(Fibs.integrated,
               resolution = c(0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5), verbose = FALSE)
Strict.Markers.samples <- list()
for(i in c("integrated_snn_res.0.05",
           "integrated_snn_res.0.1","integrated_snn_res.0.15",
           "integrated_snn_res.0.2","integrated_snn_res.0.25",
           "integrated_snn_res.0.3", "integrated_snn_res.0.4",
           "integrated_snn_res.0.5")){
  Fibs.integrated_samples.test <- AverageExpression(Fibs.integrated, assays = "RNA",
                                                    return.seurat = T, 
                                                    group.by = c(i,"SampleID2"),
                                                    slot = "data")
  Fibs.integrated_samples.test <- NormalizeData(Fibs.integrated_samples.test)
  Fibs.integrated_samples.test$Fibs_MajorClusters <- do.call(rbind, strsplit(colnames(Fibs.integrated_samples.test), "_", fixed = T))[,1]
  
  Fibs.integrated_samples.test <- SetIdent(Fibs.integrated_samples.test, value = "Fibs_MajorClusters")
  Strict.Markers.samples[[i]] <- FindAllMarkers(Fibs.integrated_samples.test, test.use = "wilcox", assay = "RNA",
                                                logfc.threshold = 1,
                                                min.pct = 0.5,
                                                only.pos = T)
  
  
}
Marker.counts <- list()
for(i in names(Strict.Markers.samples)){
  Marker.counts[[i]] <- as.data.frame(table(Strict.Markers.samples[[i]]$cluster))
  Marker.counts[[i]]$resolution.ID = i
  Marker.counts[[i]]$resolution = as.numeric(gsub("integrated_snn_res.", "", i))
}
Marker.counts <- do.call(rbind, Marker.counts)
names(Marker.counts) <- c("ClusterID", "nMarkers", "resolutionID", "resolution")

min.markers <- list()
for(i in c("integrated_snn_res.0.05",
           "integrated_snn_res.0.1", "integrated_snn_res.0.2",
           "integrated_snn_res.0.3", "integrated_snn_res.0.4",
           "integrated_snn_res.0.5")){
  min.markers[[i]] <- min(table(Strict.Markers.samples[[i]]$cluster))
}
Marker.counts$Min.n <- unlist(min.markers)[match(Marker.counts$resolutionID, names(min.markers))]  

nClusters <- aggregate(Marker.counts[, "ClusterID"], FUN = length, by = list(Marker.counts$resolution))
minMarkers <- aggregate(Marker.counts[,"nMarkers"], FUN = min, by = list(Marker.counts$resolution))
ClusteringQC_df <- minMarkers
names(ClusteringQC_df) <- c("Resolution", "Min.Markers")
ClusteringQC_df$nClusters <- nClusters$x
plot.df <- reshape2::melt(ClusteringQC_df, id.vars = "Resolution")

SupplFig_2B <- 
  plot.df %>%
  ggplot(aes(x = as.factor(Resolution))) +
  #geom_bar(aes(y = Min.Markers), hjust = 1, stat = "identity") +
  geom_bar(aes(y = value), stat = "identity") +
  facet_wrap(~variable, ncol = 1, scales = "free") +
  theme_pubr(base_size = 5) +
  xlab("Clustering resolution")
SupplFig_2B

#Supplementary Figure 2C####
Fibs.integrated <- AddMetaData(Fibs.integrated, metadata = HCLA.annotation)
SupplFig_2C <- 
  Fibs.integrated@meta.data %>%
  drop_na(HCLA.free_annotation.filtered) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = HCLA.free_annotation.filtered))+
  geom_point(data = Fibs.integrated@meta.data[, c("UMAP_1", "UMAP_2")], 
             aes(x = UMAP_1, y = UMAP_2), colour = "grey80", size = 0.1) +
  geom_point(size = 0.1) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right") +
  scale_color_npg(name = "Travaglini et al. Annotation") +
  guides(colour = guide_legend(override.aes = list(size = 2)))
SupplFig_2C



#Supplementary Figure 2D####
SupplFig_2D <- 
  Fibs.integrated@meta.data %>%
  #drop_na(HCLA.free_annotation.filtered) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = Fibs_MajorClusters))+
  facet_wrap(~Dataset, ncol = 4) +
  geom_point(data = Fibs.integrated@meta.data[, c("UMAP_1", "UMAP_2")], 
             aes(x = UMAP_1, y = UMAP_2), colour = "grey80", size = 0.1) +
  geom_point(size = 0.1) +
  theme_pubr(base_size = 5) +
  scale_colour_manual(values = Fibs_col.palette, name = "Fibroblast Subpopulations") +
  guides(colour = guide_legend(override.aes = list(size = 2)))
SupplFig_2D


#Supplementary Figure 2E####
levels(Fibs.integrated$Sample.Subtype)
Fibs.integrated$Sample.Subtype3 <-factor(
  as.character(Fibs.integrated$Sample.Subtype),
  levels = c("Control", "LUAD", "LUSC", "Carcinoid", "LCLC", "Pleiomorphic.carcinoma"),
  labels = c("Control", "LUAD", "LUSC", "Other", "Other", "Other")
)
SupplFig_2E <- 
  Fibs.integrated@meta.data %>%
  ggplot(aes(x = SampleID2, fill = Fibs_MajorClusters))+
  facet_grid(~Sample.Subtype3, scales = "free_x", space = "free_x")  +
  geom_bar(position = "fill") +
  theme_pubr(base_size = 5) +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = Fibs_col.palette, name = "Fibroblast Subpopulations") +
  ylab("Fraction") + xlab("Sample")
SupplFig_2E


#Supplementary Figure 2F####
SupplFig_2F <- 
  FibsSigs_tt.res_df[FibsSigs_tt.res_df$logFC > 0, ] %>%
  group_by(Cluster) %>%
  top_n(wt = -adj.P.Val, n = 5) %>%
  ggplot(aes(x = reorder(Pathway.label, logFC), y = logFC, label = paste0("adj.P=", signif(adj.P.Val, 3)),
             fill = Cluster)) +
  theme_pubr(base_size = 5)+
  theme(axis.title.y = element_blank(),legend.title = element_blank(),
        legend.position = "none", strip.text = element_blank(), legend.key.size = unit(10, "pt"))+
  geom_bar(stat = "identity", alpha = 0.75)+
  facet_grid(Cluster~., scales = "free", space = "free") +
  geom_text(hjust = 1.1, colour = "white", size = 2.5) +
  scale_fill_manual(values = Fibs_col.palette) +
  coord_flip()
SupplFig_2F


#Supplementary Figure 2G&H####
Fibs.integrated_samples <- SetIdent(Fibs.integrated_samples, value = "Fibs_MajorClusters")
levels(Fibs.integrated_samples@active.ident)
Fibs.Subpop.list <- SplitObject(subset(Fibs.integrated_samples, idents = c("Alveolar", "Adventitial", "Myo")),
                                split.by = "Fibs_MajorClusters")

TvN_markers_list <- list()
for(i in names(Fibs.Subpop.list)){
  Fibs.Subpop.list[[i]] <- SetIdent(Fibs.Subpop.list[[i]], value = "Sample.type")
  TvN_markers_list[[i]] <- FindMarkers(
    Fibs.Subpop.list[[i]], ident.1 = "Tumour", ident.2 = "Control", assay = "RNA", only.pos = F
  )
  TvN_markers_list[[i]]$Gene <- rownames(TvN_markers_list[[i]])
}
table(TvN_markers_list$Myo$p_val_adj < 0.01)
TvN_markers_list$Myo$label <- TvN_markers_list$Myo$Gene
TvN_markers_list$Myo$label[TvN_markers_list$Myo$p_val_adj > 0.01] <- NA

SupplFig_2H <- 
  TvN_markers_list$Myo %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = p_val_adj < 0.01)) + 
  theme_pubr(base_size = 5) +
  theme(legend.position = "none") +
  geom_point(aes(size = -log10(p_val_adj))) +
  scale_size(range = c(0.1,2)) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("Myofibroblasts") +
  xlab("Tumour v Control (logFC)")+
  ggrepel::geom_text_repel(aes(label = label), fontface = "italic", size = 2, min.segment.length = 0)
SupplFig_2H

TvN_markers_list$Adventitial$label <- TvN_markers_list$Adventitial$Gene
TvN_markers_list$Adventitial$label[TvN_markers_list$Adventitial$p_val_adj > 0.01] <- NA
SupplFig_2G <- 
  TvN_markers_list$Adventitial %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = p_val_adj < 0.01)) + 
  theme_pubr(base_size = 5) +
  theme(legend.position = "none") +
  geom_point(aes(size = -log10(p_val_adj))) +
  scale_size(range = c(0.1,2)) +
  ggtitle("Adventitial fibroblasts") +
  xlab("Tumour v Control (logFC)")+
  scale_color_manual(values = c("black", "red")) +
  xlim(c(-2.5, 2.5))+
  ggrepel::geom_text_repel(aes(label = label), fontface = "italic", size = 2, min.segment.length = 0)
SupplFig_2G


#Supplementary Data 5####
for(i in names(TvN_markers_list)){
  TvN_markers_list[[i]]$Cluster <- i
  TvN_markers_list[[i]]$label <- TvN_markers_list[[i]]$Gene
}
TvN_markers.df <- do.call(rbind, TvN_markers_list)
Table_S5 <- TvN_markers.df


