#Figure 6####
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
load(paste0(data_directory, "IntegratedFibs_Zenodo.Rdata"))

#Univariate Cox regression####
#LUAD
#Calculate period for which CoxPH assumptions hold
ys2test <- rev(seq(1,10, by = 1))
cox.zph_OK <- list()
for(DS in levels(Merged.LUADtraits$Dataset.factor)[]){
  surv_data <- Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                               Merged.LUADtraits$Dataset.factor == DS, ]
  
  cox.zph_p = 0
  for(i in ys2test){
    if (cox.zph_p < 0.05) {
      split_data <- survSplit(Surv(OS_YEARS, OS) ~ 
                                Myo_Fibs.pct,# + Stage_4cat + Age,
                              zero = -0.1,
                              data = surv_data[!is.na(surv_data$OS_YEARS), ],
                              cut = i, episode = "tgroup", id = "id")
      model.coxph2 <- coxph(Surv(OS_YEARS, OS) ~ 
                              Myo_Fibs.pct,
                            data = split_data[split_data$tgroup == 1,])
      cox.zph_res = cox.zph(model.coxph2)
      cox.zph_res.sum = c(i, cox.zph_res$table[1,3])
      cox.zph_p <- cox.zph_res$table[1,3]
      names(cox.zph_res.sum) <- c("OS_YEARS", "cox.zph_p")
    }
    cox.zph_OK[[DS]] <- cox.zph_res.sum
  }
}
do.call(rbind, cox.zph_OK)
max_year <- min(do.call(rbind, cox.zph_OK)[,1])

#Perform Cox regression analysis on each dataset individually
split_data <- list()
for(DS in levels(Merged.LUADtraits$Dataset.factor)[]){
  surv_data <- Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                               Merged.LUADtraits$Dataset.factor == DS, ]
  split_data[[DS]] <- survSplit(Surv(OS_YEARS, OS) ~ 
                                  Myo_Fibs.pct + Alveolar_Fibs.pct + Adventitial_Fibs.pct,
                                zero = -0.1,
                                data = surv_data[!is.na(surv_data$OS_YEARS), ],
                                cut = max_year, episode = "tgroup", id = "id")
}
Myo_cox.res <- list()
Myo_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Myo_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Myo_cox.res[["GSE68465"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Myo_Fibs.pct),
                                   data =  split_data[["GSE68465"]][split_data[["GSE68465"]]$tgroup == 1,])
Myo_cox.res[["GSE31210"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Myo_Fibs.pct),
                                   data =  split_data[["GSE31210"]][split_data[["GSE31210"]]$tgroup == 1,])
Myo_cox.res[["GSE72094"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      scale(Myo_Fibs.pct),
                                    data =  split_data[["GSE72094"]][split_data[["GSE72094"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Myo_cox.res, panels = panels)
Myo_LUAD_Zscores <- list()
for(DS in names(Myo_cox.res)){
  Myo_LUAD_Zscores[[DS]] <- summary(Myo_cox.res[[DS]])$coefficients[4:5]
}
Myo_LUAD_pvals <- list()
for(DS in names(Myo_cox.res)){
  Myo_LUAD_pvals[[DS]] <- summary(Myo_cox.res[[DS]])$coefficients[5]
}
Myo_LUAD_pvals
metap::sumz(unlist(Myo_LUAD_pvals))


Alveolar_cox.res <- list()
Alveolar_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Alveolar_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Alveolar_cox.res[["GSE68465"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Alveolar_Fibs.pct),
                                   data =  split_data[["GSE68465"]][split_data[["GSE68465"]]$tgroup == 1,])
Alveolar_cox.res[["GSE31210"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Alveolar_Fibs.pct),
                                   data =  split_data[["GSE31210"]][split_data[["GSE31210"]]$tgroup == 1,])
Alveolar_cox.res[["GSE72094"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      scale(Alveolar_Fibs.pct),
                                    data =  split_data[["GSE72094"]][split_data[["GSE72094"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Alveolar_cox.res, panels = panels)
Alveolar_LUAD_Zscores <- list()
for(DS in names(Alveolar_cox.res)){
  Alveolar_LUAD_Zscores[[DS]] <- summary(Alveolar_cox.res[[DS]])$coefficients[4:5]
}
Alveolar_LUAD_pvals <- list()
for(DS in names(Alveolar_cox.res)){
  Alveolar_LUAD_pvals[[DS]] <- summary(Alveolar_cox.res[[DS]])$coefficients[5]
}
Alveolar_LUAD_pvals
metap::sumz(unlist(Alveolar_LUAD_pvals))

Adventitial_cox.res <- list()
Adventitial_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Adventitial_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Adventitial_cox.res[["GSE68465"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Adventitial_Fibs.pct),
                                   data =  split_data[["GSE68465"]][split_data[["GSE68465"]]$tgroup == 1,])
Adventitial_cox.res[["GSE31210"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Adventitial_Fibs.pct),
                                   data =  split_data[["GSE31210"]][split_data[["GSE31210"]]$tgroup == 1,])
Adventitial_cox.res[["GSE72094"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      scale(Adventitial_Fibs.pct),
                                    data =  split_data[["GSE72094"]][split_data[["GSE72094"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Adventitial_cox.res, panels = panels)
Adventitial_LUAD_Zscores <- list()
for(DS in names(Adventitial_cox.res)){
  Adventitial_LUAD_Zscores[[DS]] <- summary(Adventitial_cox.res[[DS]])$coefficients[4:5]
}

LUAD_zScores.df <- data.frame(Z = do.call(rbind, Myo_LUAD_Zscores)[,1],
                              p = do.call(rbind, Myo_LUAD_Zscores)[,2],
                              SubPop = "Myo",
                              DS = rownames(do.call(rbind, Myo_LUAD_Zscores)))
LUAD_zScores.df <- rbind(LUAD_zScores.df,
                         data.frame(Z = do.call(rbind, Alveolar_LUAD_Zscores)[,1],
                                    p = do.call(rbind, Alveolar_LUAD_Zscores)[,2],
                                    SubPop = "Alveolar",
                                    DS = rownames(do.call(rbind, Alveolar_LUAD_Zscores))))
LUAD_zScores.df <- rbind(LUAD_zScores.df,
                         data.frame(Z = do.call(rbind, Adventitial_LUAD_Zscores)[,1],
                                    p = do.call(rbind, Adventitial_LUAD_Zscores)[,2],
                                    SubPop = "Adventitial",
                                    DS = rownames(do.call(rbind, Adventitial_LUAD_Zscores))))
LUAD_zScores.df$Subtype <- "LUAD"

reshape2::melt(LUAD_zScores.df, id.vars = c("SubPop", "DS"), measure.vars = "Z") %>%
  ggplot(aes(y = value, x = SubPop, fill = SubPop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1) + 
  theme_pubr(base_size = 7) +
  scale_fill_manual(values = Fibs_col.palette) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Good Survival <- CoxPH Zscore -> Poor Survival") +
  theme(axis.title.y = element_blank()) +
  coord_flip()
  


#Univariate Cox regression LUSC####
#LUSC
#Calculate period for which CoxPH assumptions hold
ys2test <- rev(seq(1,10, by = 1))
cox.zph_OK <- list()
for(DS in unique(Merged.LUSCtraits$Dataset)){
  surv_data <- Merged.LUSCtraits[Merged.LUSCtraits$Sample.Subtype == "LUSC" &
                                   Merged.LUSCtraits$Dataset == DS, ]
  
  cox.zph_p = 0
  for(i in ys2test){
    if (cox.zph_p < 0.05) {
      split_data <- survSplit(Surv(OS_YEARS, OS) ~ 
                                Myo_Fibs.pct,# + Stage_4cat + Age,
                              zero = -0.1,
                              data = surv_data[!is.na(surv_data$OS_YEARS), ],
                              cut = i, episode = "tgroup", id = "id")
      model.coxph2 <- coxph(Surv(OS_YEARS, OS) ~ 
                              Myo_Fibs.pct,
                            data = split_data[split_data$tgroup == 1,])
      cox.zph_res = cox.zph(model.coxph2)
      cox.zph_res.sum = c(i, cox.zph_res$table[1,3])
      cox.zph_p <- cox.zph_res$table[1,3]
      names(cox.zph_res.sum) <- c("OS_YEARS", "cox.zph_p")
    }
    cox.zph_OK[[DS]] <- cox.zph_res.sum
  }
}
do.call(rbind, cox.zph_OK)
max_year <- 4

#Perform Cox regression analysis on each dataset individually
split_data <- list()
for(DS in unique(Merged.LUSCtraits$Dataset)[]){
  surv_data <- Merged.LUSCtraits[Merged.LUSCtraits$Sample.Subtype == "LUSC" &
                                   Merged.LUSCtraits$Dataset == DS, ]
  split_data[[DS]] <- survSplit(Surv(OS_YEARS, OS) ~ 
                                  Myo_Fibs.pct + Alveolar_Fibs.pct + Adventitial_Fibs.pct,
                                zero = -0.1,
                                data = surv_data[!is.na(surv_data$OS_YEARS), ],
                                cut = max_year, episode = "tgroup", id = "id")
}
Myo_cox.res <- list()
Myo_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Myo_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Myo_cox.res[["GSE157009"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Myo_Fibs.pct),
                                   data =  split_data[["GSE157009"]][split_data[["GSE157009"]]$tgroup == 1,])
Myo_cox.res[["GSE157010"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Myo_Fibs.pct),
                                   data =  split_data[["GSE157010"]][split_data[["GSE157010"]]$tgroup == 1,])
Myo_cox.res[["GSE4573"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      scale(Myo_Fibs.pct),
                                    data =  split_data[["GSE4573"]][split_data[["GSE4573"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Myo_cox.res, panels = panels)
Myo_LUSC_Zscores <- list()
for(DS in names(Myo_cox.res)){
  Myo_LUSC_Zscores[[DS]] <- summary(Myo_cox.res[[DS]])$coefficients[4:5]
}

Alveolar_cox.res <- list()
Alveolar_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Alveolar_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Alveolar_cox.res[["GSE157009"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                       scale(Alveolar_Fibs.pct),
                                     data =  split_data[["GSE157009"]][split_data[["GSE157009"]]$tgroup == 1,])
Alveolar_cox.res[["GSE157010"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                       scale(Alveolar_Fibs.pct),
                                     data =  split_data[["GSE157010"]][split_data[["GSE157010"]]$tgroup == 1,])
Alveolar_cox.res[["GSE4573"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Alveolar_Fibs.pct),
                                   data =  split_data[["GSE4573"]][split_data[["GSE4573"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Alveolar_cox.res, panels = panels)
Alveolar_LUSC_Zscores <- list()
for(DS in names(Alveolar_cox.res)){
  Alveolar_LUSC_Zscores[[DS]] <- summary(Alveolar_cox.res[[DS]])$coefficients[4:5]
}

Adventitial_cox.res <- list()
Adventitial_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  scale(Adventitial_Fibs.pct),
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Adventitial_cox.res[["GSE157009"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                       scale(Adventitial_Fibs.pct),
                                     data =  split_data[["GSE157009"]][split_data[["GSE157009"]]$tgroup == 1,])
Adventitial_cox.res[["GSE157010"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                       scale(Adventitial_Fibs.pct),
                                     data =  split_data[["GSE157010"]][split_data[["GSE157010"]]$tgroup == 1,])
Adventitial_cox.res[["GSE4573"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                     scale(Adventitial_Fibs.pct),
                                   data =  split_data[["GSE4573"]][split_data[["GSE4573"]]$tgroup == 1,])

forestmodel::forest_model(model_list = Adventitial_cox.res, panels = panels)
Adventitial_LUSC_Zscores <- list()
for(DS in names(Adventitial_cox.res)){
  Adventitial_LUSC_Zscores[[DS]] <- summary(Adventitial_cox.res[[DS]])$coefficients[4:5]
}

LUSC_zScores.df <- data.frame(Z = do.call(rbind, Myo_LUSC_Zscores)[,1],
                              p = do.call(rbind, Myo_LUSC_Zscores)[,2],
                              SubPop = "Myo",
                              DS = rownames(do.call(rbind, Myo_LUSC_Zscores)))
LUSC_zScores.df <- rbind(LUSC_zScores.df,
                         data.frame(Z = do.call(rbind, Alveolar_LUSC_Zscores)[,1],
                                    p = do.call(rbind, Alveolar_LUSC_Zscores)[,2],
                                    SubPop = "Alveolar",
                                    DS = rownames(do.call(rbind, Alveolar_LUSC_Zscores))))
LUSC_zScores.df <- rbind(LUSC_zScores.df,
                         data.frame(Z = do.call(rbind, Adventitial_LUSC_Zscores)[,1],
                                    p = do.call(rbind, Adventitial_LUSC_Zscores)[,2],
                                    SubPop = "Adventitial",
                                    DS = rownames(do.call(rbind, Adventitial_LUSC_Zscores))))
LUSC_zScores.df$Subtype <- "LUSC"

#Supplementary Figure 6A####
Cox.Zscores <- rbind(reshape2::melt(LUAD_zScores.df, id.vars = c("SubPop", "DS", "Subtype"), measure.vars = "Z"),
                     reshape2::melt(LUSC_zScores.df, id.vars = c("SubPop", "DS", "Subtype"), measure.vars = "Z"))

names(Cox.Zscores)

Cox.UV.plot <- 
  Cox.Zscores %>%
  ggplot(aes(y = value, x = SubPop, fill = SubPop)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(size = 1) + 
  theme_pubr(base_size = 7) +
  scale_fill_manual(values = Fibs_col.palette) +
  facet_grid(Subtype~.) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = c(1.959964, -1.959964), linetype = "dashed") +
  geom_hline(yintercept =  c(-2.575829, 2.575829), linetype = "dashed") +
  ylab("Good Survival <- CoxPH Zscore -> Poor Survival") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0)) +
  #annotate("label", label = "*", y = c(-1.959964, 1.9596964), x = 1,
   #        size = 5, vjust = 1, label.size = NA) +
  #annotate("label", label = "**", y = c(-2.575829, 2.575829), x = 1,
  #         size = 5, vjust = 1, label.size = NA) +
  coord_flip()


p.labels <- data.frame(value = c(-1.959964, 1.9596964,-2.575829, 2.575829),
                       label = c("*", "*", "**", "**"),
                       Subtype = rep("LUSC", 4),
                       variable = rep("Myo", 4))

SupplFig_6A <- 
  Cox.UV.plot + geom_label(aes(label = label), data = p.labels,
                        x = 0.75, fill = "white")

#LUAD v LUSC Marker analysis####
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
names(Fibs.integrated_samples@meta.data)
keepCols_metaData <- c("SampleID", "nCount_RNA", "nFeature_RNA",
                       "Fibs_MajorClusters", "Dataset")
Fibs.integrated_samples@meta.data <- Fibs.integrated_samples@meta.data[, keepCols_metaData]
Fibs.integrated_samples <- AddMetaData(Fibs.integrated_samples, Fibs.Sample_MetaData)
Fibs.integrated_samples$Sample.type[Fibs.integrated_samples$Sample.type == "Tumor"] <- "Tumour"
Fibs.integrated_samples$pN <- droplevels(Fibs.integrated_samples$pN, exclude = c(NA, "NA"))
levels(Fibs.integrated_samples$pN)

Fibs.integrated_samples <- SetIdent(Fibs.integrated_samples, value = "Fibs_MajorClusters")
levels(Fibs.integrated_samples@active.ident)
Fibs.Subpop.list <- SplitObject(subset(Fibs.integrated_samples, idents = c("Alveolar", "Adventitial", "Myo")),
                                split.by = "Fibs_MajorClusters")

LUADvLUSC_markers_list <- list()
for(i in names(Fibs.Subpop.list)){
  Fibs.Subpop.list[[i]] <- SetIdent(Fibs.Subpop.list[[i]], value = "Sample.Subtype")
  LUADvLUSC_markers_list[[i]] <- FindMarkers(
    Fibs.Subpop.list[[i]], ident.1 = "LUAD", ident.2 = "LUSC", assay = "RNA", only.pos = F
  )
  LUADvLUSC_markers_list[[i]]$Gene <- rownames(LUADvLUSC_markers_list[[i]])
}
table(LUADvLUSC_markers_list$Myo$p_val_adj < 0.01)
LUADvLUSC_markers_list$Myo$label <- LUADvLUSC_markers_list$Myo$Gene
LUADvLUSC_markers_list$Myo$label[LUADvLUSC_markers_list$Myo$p_val_adj > 0.01] <- NA

#Supplementary Figure 6B####
SupplFig_6B <- 
  LUADvLUSC_markers_list$Myo %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = p_val_adj < 0.01)) + 
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  geom_point(aes(size = -log10(p_val_adj))) +
  scale_size(range = c(0.1,2)) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("Myofibroblasts") +
  xlab("LUAD v LUSC (logFC)")+
  ggrepel::geom_text_repel(aes(label = label), fontface = "italic", size = 2.5, min.segment.length = 0) +
  ylim(0,5) + #xlim(-3,3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate(geom = "label", label = "*adj.p=0.05",
           x = 3, y = -log10(0.05), size = 2.5, hjust = 0.75)
SupplFig_6B

#Kaplan Meier Analysis####
variables <- c("Myo_Fibs.pct", "Alveolar_Fibs.pct")
TCGA_LUAD <- Merged.LUADtraits %>%
  filter(Dataset == "TCGA" & Sample.Subtype == "LUAD")
LUAD.Opt_cut <- surv_cutpoint(TCGA_LUAD,
                              time = "OS_YEARS", event = "OS", 
                              variables,
                              minprop = 0.01, progressbar = TRUE)
Opt_cut.plot <- plot(LUAD.Opt_cut)

LUAD.cutpoints <- LUAD.Opt_cut$cutpoint 
Fibs_cut <- data.frame(
  Alv = factor(Hmisc::cut2(Merged.LUADtraits$Alveolar_Fibs.pct, cuts = LUAD.cutpoints[2,1]),
               labels = c("Low", "High")),
  Myo = factor(Hmisc::cut2(Merged.LUADtraits$Myo_Fibs.pct, cuts = LUAD.cutpoints[1,1]),
               labels = c("Low", "High"))
)
Merged.LUADtraits$Myo_cat <- Fibs_cut$Myo
Merged.LUADtraits$Alv_cat <- Fibs_cut$Alv

Merged.LUADtraits$Combined_strat <- paste0("Alv.", Merged.LUADtraits$Alv_cat,
                                           "_Myo.", Merged.LUADtraits$Myo_cat)


#Figure 6 A-D####
Fig_6A <- 
  data.frame(stat = LUAD.Opt_cut$Myo_Fibs.pct$stats,
           x = LUAD.Opt_cut$Myo_Fibs.pct$cuts) %>%
  ggplot(aes(x = x, y = stat, colour = x < LUAD.Opt_cut$cutpoint[1,1])) +
  geom_point(size = 0.5) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  scale_colour_manual(values =  c(pal_d3()(3)[3], "grey70")) +
  xlab("Myo (% of all fibroblasts)") +
  ylab("Standardized Log-Rank statistic") +
  geom_vline(xintercept = LUAD.Opt_cut$cutpoint[1,1], 
             linetype = "dashed") +
  annotate("label", y = 1, size = 2.5,
           label = paste("Max ranked\nCutpoint =",
                                  signif(LUAD.Opt_cut$cutpoint[1,1],3)),
           x =  LUAD.Opt_cut$cutpoint[1,1],
           label.size = NA)
Fig_6A

Fig_6B <- 
  Opt_cut.plot$Myo_Fibs.pct$distribution +
  theme_pubr(base_size = 7) +
  scale_fill_manual(values =  c(pal_d3()(3)[3], "grey70"))+
  scale_colour_manual(values =  c(pal_d3()(3)[3], "grey70")) +
  theme(legend.position = "none", plot.title = element_blank()) +
  xlab("Myo (% of all fibroblasts)")
Fig_6B  

Myo.ggsurvplot_list <- list()
for(i in levels(Merged.LUADtraits$Dataset.factor)){
  cox_res <- summary(coxph(Surv(OS_YEARS, OS) ~ Myo_cat,
                           data = Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                                                  Merged.LUADtraits$Dataset.factor == i, ]))
  Myo.ggsurvplot_list[[i]] <- ggsurvplot(survfit(Surv(OS_YEARS, OS) ~ Myo_cat,
                                                 data = Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                                                                        Merged.LUADtraits$Dataset.factor == i, ]),
                                         main = i,
                                         pval = signif(cox_res$logtest[3],3), pval.coord = c(0,0.15), pval.size = 2.5,
                                         conf.int = F, 
                                         risk.table = TRUE, risk.table.title = element_blank(),
                                         risk.table.fontsize = 2.5,tables.height = 0.4,
                                         risk.table.pos = "in",
                                         font.tickslab = 7,
                                         censor.shape = 124, censor.size = 1,
                                         legend.labs = c("Low", "High"),
                                         palette = c("grey70",pal_d3()(3)[3]),
                                         #surv.median.line = "hv",
                                         ggtheme = theme_pubr(base_size = 7),
                                         #tables.theme = theme_cleantable(),
                                         risk.table.y.col = T,
                                         legend = c(0.85,0.85), legend.title = element_blank())
  Myo.ggsurvplot_list[[i]]$plot <- Myo.ggsurvplot_list[[i]]$plot + 
    ggtitle(i) + ylab("Overall Survival Rate") +
    theme(#plot.title.position = "plot",
          axis.text.x = element_blank(), axis.title.x =  element_blank(),
          plot.margin = margin(b= 0)) +
    coord_cartesian(xlim = c(-0.3,NA))
    
  Myo.ggsurvplot_list[[i]]$table <- Myo.ggsurvplot_list[[i]]$table + xlab("Time (Years)")+
    theme(plot.margin = margin(t = 0), plot.background = element_blank()) +
    coord_cartesian(xlim = c(-0.3,NA))
   }


Fig_6C <- 
  ggarrange(Myo.ggsurvplot_list$TCGA$plot + theme(plot.title = element_blank()),
          Myo.ggsurvplot_list$TCGA$table,
          ncol = 1, nrow = 2, heights = c(0.7,0.3), align = "v",
          common.legend = T, legend = "none")
         
Fig_6C



Fig_6D <- 
  ggarrange(Myo.ggsurvplot_list$GSE72094$plot, Myo.ggsurvplot_list$GSE31210$plot,
           Myo.ggsurvplot_list$GSE68465$plot,
          Myo.ggsurvplot_list$GSE72094$table,  Myo.ggsurvplot_list$GSE31210$table,
           Myo.ggsurvplot_list$GSE68465$table,
          ncol = 3, nrow = 2, heights = c(0.7,0.3), align = "v",
          common.legend = T, legend = "none")
Fig_6D



#Figure 6 E-H####
Fig_6E <- 
  data.frame(stat = LUAD.Opt_cut$Alveolar_Fibs.pct$stats,
             x = LUAD.Opt_cut$Alveolar_Fibs.pct$cuts) %>%
  ggplot(aes(x = x, y = stat, colour = x < LUAD.Opt_cut$cutpoint[2,1])) +
  geom_point(size = 0.5) +
  theme_pubr(base_size = 7) +
  theme(legend.position = "none") +
  scale_colour_manual(values =  c(pal_d3()(3)[1], "grey70")) +
  xlab("Alveolar (% of all fibroblasts)") +
  ylab("Standardized Log-Rank statistic") +
  geom_vline(xintercept = LUAD.Opt_cut$cutpoint[2,1], 
             linetype = "dashed") +
  annotate("label", y = 1, size = 2.5,
           label = paste("Max ranked\nCutpoint =",
                         signif(LUAD.Opt_cut$cutpoint[2,1],3)),
           x =  LUAD.Opt_cut$cutpoint[2,1],
           label.size = NA)
Fig_6E

Fig_6F <- 
  Opt_cut.plot$Alveolar_Fibs.pct$distribution +
  theme_pubr(base_size = 7) +
  scale_fill_manual(values =  c(pal_d3()(3)[1], "grey70"))+
  scale_colour_manual(values =  c(pal_d3()(3)[1], "grey70")) +
  theme(legend.position = "none", plot.title = element_blank()) +
  xlab("Alveolar (% of all fibroblasts)")
Fig_6F  


Alv.ggsurvplot_list <- list()
for(i in levels(Merged.LUADtraits$Dataset.factor)){
  cox_res <- summary(coxph(Surv(OS_YEARS, OS) ~ Alv_cat,
                           data = Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                                                      Merged.LUADtraits$Dataset.factor == i, ]))
  Alv.ggsurvplot_list[[i]] <- ggsurvplot(survfit(Surv(OS_YEARS, OS) ~ Alv_cat,
                                                 data = Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                                                                            Merged.LUADtraits$Dataset.factor == i, ]),
                                         main = i,
                                         pval = T, pval.coord = c(0,0.15), pval.size = 2.5,
                                         conf.int = F, 
                                         risk.table = TRUE, risk.table.title = element_blank(),
                                         risk.table.fontsize = 2.5,tables.height = 0.4,
                                         risk.table.pos = "in",
                                         font.tickslab = 7,
                                         censor.shape = 124, censor.size = 1,
                                         legend.labs = c("Low", "High"),
                                         palette = c("grey70",pal_d3()(3)[1]),
                                         #surv.median.line = "hv",
                                         ggtheme = theme_pubr(base_size = 7),
                                         #tables.theme = theme_cleantable(),
                                         risk.table.y.col = T,
                                         legend = c(0.85,0.85), legend.title = element_blank())
  Alv.ggsurvplot_list[[i]]$plot <- Alv.ggsurvplot_list[[i]]$plot + 
    ggtitle(i) + ylab("Overall Survival Rate") +
    theme(#plot.title.position = "plot",
          axis.text.x = element_blank(), axis.title.x =  element_blank(),
          plot.margin = margin(b= 0)) +
    coord_cartesian(xlim = c(-0.3,NA))
  
  Alv.ggsurvplot_list[[i]]$table <- Alv.ggsurvplot_list[[i]]$table + xlab("Time (Years)")+
    theme(plot.margin = margin(t = 0), plot.background = element_blank()) +
    coord_cartesian(xlim = c(-0.3,NA))
}

Fig_6G <- 
  ggarrange(Alv.ggsurvplot_list$TCGA$plot + theme(plot.title = element_blank()),
            Alv.ggsurvplot_list$TCGA$table,
            ncol = 1, nrow = 2, heights = c(0.7,0.3), align = "v",
            common.legend = T, legend = "none")
Fig_6G

Fig_6H <- 
  ggarrange(Alv.ggsurvplot_list$GSE72094$plot, Alv.ggsurvplot_list$GSE31210$plot,
          Alv.ggsurvplot_list$GSE68465$plot,
          Alv.ggsurvplot_list$GSE72094$table,  Alv.ggsurvplot_list$GSE31210$table,
          Alv.ggsurvplot_list$GSE68465$table,
          ncol = 3, nrow = 2, heights = c(0.7,0.3), align = "v",
          common.legend = T, legend = "none")
Fig_6H


#Categorical analysis Multivariate####
#Figure 6I####
max_year <- 4
surv_data <- Merged.LUADtraits[Merged.LUADtraits$Dataset.factor %in% c("TCGA", "GSE68465", "GSE72094", "GSE31210"), ]
split_data <- survSplit(Surv(OS_YEARS, OS) ~ 
                          Myo_cat + Stage_4cat + Age,
                        zero = -0.1,
                        data = surv_data[!is.na(surv_data$OS_YEARS), ],
                        cut = max_year, episode = "tgroup", id = "id")
names(split_data)[1:2] <- c("Myo", "Stage")
Myo_cox.res <- coxph(Surv(OS_YEARS, OS) ~ 
                       Myo + Stage + Age,
                     data =  split_data[split_data$tgroup == 1,])
summary(Myo_cox.res)$coefficients


Fig_6I <- forest_model(Myo_cox.res,
                   format_options = forest_model_format_options(
                     banded = T, text_size = 2.5, point_size = 1)) 
Fig_6I

#Supplementary Figure 6C####
split_data <- list()
for(DS in levels(Merged.LUADtraits$Dataset.factor)[]){
  surv_data <- Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                               Merged.LUADtraits$Dataset.factor == DS, ]
  split_data[[DS]] <- survSplit(Surv(OS_YEARS, OS) ~ 
                            Myo_cat + Stage_4cat + Age,
                          zero = -0.1,
                          data = surv_data[!is.na(surv_data$OS_YEARS), ],
                          cut = max_year, episode = "tgroup", id = "id")
}
Myo_cox.res <- list()
Myo_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  Myo_cat + Stage_4cat + Age,
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Myo_cox.res[["GSE68465"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                          Myo_cat + Stage_4cat + Age,
                                        data =  split_data[["GSE68465"]][split_data[["GSE68465"]]$tgroup == 1,])
Myo_cox.res[["GSE31210"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                          Myo_cat + Stage_4cat ,
                                        data =  split_data[["GSE31210"]][split_data[["GSE31210"]]$tgroup == 1,])
Myo_cox.res[["GSE72094"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                           Myo_cat + Stage_4cat + Age,
                                         data =  split_data[["GSE72094"]][split_data[["GSE72094"]]$tgroup == 1,])

SupplFig_6C <- 
  forest_model(model_list = Myo_cox.res, panels = panels,
             format_options = forest_model_format_options(
               banded = T, text_size = 2.5, point_size = 2))
SupplFig_6C

#Figure 6J####
surv_data <- Merged.LUADtraits[Merged.LUADtraits$Dataset.factor %in% c("TCGA", "GSE68465", "GSE72094", "GSE31210"), ]
split_data <- survSplit(Surv(OS_YEARS, OS) ~ 
                          Alv_cat + Stage_4cat + Age,
                        zero = -0.1,
                        data = surv_data[!is.na(surv_data$OS_YEARS), ],
                        cut = max_year, episode = "tgroup", id = "id")
names(split_data)[1:2] <- c("Alveolar", "Stage")
Alv_cox.res <- coxph(Surv(OS_YEARS, OS) ~ 
                       Alveolar + Stage + Age,
                     data =  split_data[split_data$tgroup == 1,])

Fig_6J <- forest_model(Alv_cox.res,
                    format_options = forest_model_format_options(
                      banded = T, text_size = 2.5, point_size = 1)) 
Fig_6J

#Supplementary Figure 6D####
max_year <- 4
split_data <- list()
for(DS in levels(Merged.LUADtraits$Dataset.factor)[]){
  surv_data <- Merged.LUADtraits[Merged.LUADtraits$Sample.Subtype == "LUAD" &
                                   Merged.LUADtraits$Dataset.factor == DS, ]
  split_data[[DS]] <- survSplit(Surv(OS_YEARS, OS) ~ 
                                  Alv_cat + Stage_4cat + Age,
                                zero = -0.1,
                                data = surv_data[!is.na(surv_data$OS_YEARS), ],
                                cut = max_year, episode = "tgroup", id = "id")
}
Alv_cox.res <- list()
Alv_cox.res[["TCGA"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                  Alv_cat + Stage_4cat + Age,
                                data =  split_data[["TCGA"]][split_data[["TCGA"]]$tgroup == 1,])
Alv_cox.res[["GSE68465"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      Alv_cat + Stage_4cat + Age,
                                    data =  split_data[["GSE68465"]][split_data[["GSE68465"]]$tgroup == 1,])
Alv_cox.res[["GSE31210"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      Alv_cat + Stage_4cat ,
                                    data =  split_data[["GSE31210"]][split_data[["GSE31210"]]$tgroup == 1,])
Alv_cox.res[["GSE72094"]]  <- coxph(Surv(OS_YEARS, OS) ~ 
                                      Alv_cat + Stage_4cat + Age,
                                    data =  split_data[["GSE72094"]][split_data[["GSE72094"]]$tgroup == 1,])
SupplFig_6D <- 
  forest_model(model_list = Alv_cox.res, panels = panels,
             format_options = forest_model_format_options(
               banded = T, text_size = 2.5, point_size = 2))

SupplFig_6D

