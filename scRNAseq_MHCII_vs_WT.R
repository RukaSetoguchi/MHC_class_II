library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
set.seed(1234)

#KO1
memory_KO.data <- Read10X(data.dir = "~/data/scRNAseq/201226_NovaSeq_II_KO/filtered_feature_bc_matrix")
memory_KO <- CreateSeuratObject(counts = memory_KO.data, project = "memory_KO", min.cells = 3, min.features = 200)
memory_KO[["percent.mt"]] <- PercentageFeatureSet(memory_KO, pattern = "^mt-")
head(memory_KO@meta.data)
length(memory_KO@meta.data$nFeature_RNA)
png("VlnPlot_QC.png", width = 800, height = 500)
VlnPlot(memory_KO, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
png("FeatureScatter.png", width = 800, height = 500)
plot1 <- FeatureScatter(memory_KO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(memory_KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
memory_KO <- subset(memory_KO, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
length(memory_KO@meta.data$nFeature_RNA)
png("VlnPlot_QC_subset.png", width = 800, height = 500)
VlnPlot(memory_KO, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
memory_KO <- NormalizeData(memory_KO, normalization.method = "LogNormalize", scale.factor = 10000)
save (memory_KO, file="memory_KO_normalization.Rdata")

#WT1
memory_WT.data <- Read10X(data.dir = "~/data/scRNAseq/201226_NovaSeq_WT/filtered_feature_bc_matrix")
memory_WT <- CreateSeuratObject(counts = memory_WT.data, project = "memory_WT", min.cells = 3, min.features = 200)
memory_WT[["percent.mt"]] <- PercentageFeatureSet(memory_WT, pattern = "^mt-")
head(memory_WT@meta.data)
length(memory_WT@meta.data$nFeature_RNA)
png("VlnPlot_QC.png", width = 800, height = 500)
VlnPlot(memory_WT, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
png("FeatureScatter.png", width = 800, height = 500)
plot1 <- FeatureScatter(memory_WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(memory_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
memory_WT <- subset(memory_WT, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
length(memory_WT@meta.data$nFeature_RNA)
png("VlnPlot_QC_subset.png", width = 800, height = 500)
VlnPlot(memory_WT, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
memory_WT <- NormalizeData(memory_WT, normalization.method = "LogNormalize", scale.factor = 10000)
save (memory_WT, file="memory_WT_normalization.Rdata")

#KO2
memory_KO2.data <- Read10X(data.dir = "~/work/210426_NovaSeq_II_KO_2/filtered_feature_bc_matrix")
memory_KO2 <- CreateSeuratObject(counts = memory_KO2.data, project = "memory_KO2", min.cells = 3, min.features = 200)
memory_KO2[["percent.mt"]] <- PercentageFeatureSet(memory_KO2, pattern = "^mt-")
head(memory_KO2@meta.data)
length(memory_KO2@meta.data$nFeature_RNA)
png("VlnPlot_QC.png", width = 800, height = 500)
VlnPlot(memory_KO2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
png("FeatureScatter.png", width = 800, height = 500)
plot1 <- FeatureScatter(memory_KO2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(memory_KO2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
memory_KO2 <- subset(memory_KO2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
length(memory_KO2@meta.data$nFeature_RNA)
png("VlnPlot_QC_subset.png", width = 800, height = 500)
VlnPlot(memory_KO2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
memory_KO2 <- NormalizeData(memory_KO2, normalization.method = "LogNormalize", scale.factor = 10000)
save (memory_KO2, file="memory_KO2_normalization.Rdata")

#WT2
memory_WT2.data <- Read10X(data.dir = "~/work/210426_NovaSeq_WT_2/filtered_feature_bc_matrix")
memory_WT2 <- CreateSeuratObject(counts = memory_WT2.data, project = "memory_WT2", min.cells = 3, min.features = 200)
memory_WT2[["percent.mt"]] <- PercentageFeatureSet(memory_WT2, pattern = "^mt-")
head(memory_WT2@meta.data)
length(memory_WT2@meta.data$nFeature_RNA)
png("VlnPlot_QC.png", width = 800, height = 500)
VlnPlot(memory_WT2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
png("FeatureScatter.png", width = 800, height = 500)
plot1 <- FeatureScatter(memory_WT2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(memory_WT2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
memory_WT2 <- subset(memory_WT2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
length(memory_WT2@meta.data$nFeature_RNA)
png("VlnPlot_QC_subset.png", width = 800, height = 500)
VlnPlot(memory_WT2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
memory_WT2 <- NormalizeData(memory_WT2, normalization.method = "LogNormalize", scale.factor = 10000)
save (memory_WT2, file="memory_WT2_normalization.Rdata")


#harmony

m(list = ls())
library(tidyverse)
library(harmony)

****************
  library(effsize)

load("memory_KO_normalization.Rdata")
saveRDS(memory_KO, file = "memory_KO1.rds")

load("memory_KO2_normalization.Rdata")
saveRDS(memory_KO2, file = "memory_KO2.rds")

load("memory_WT_normalization.Rdata")
saveRDS(memory_WT, file = "memory_WT1.rds")

load("memory_WT2_normalization.Rdata")
saveRDS(memory_WT2, file = "memory_WT2.rds")

KO1 <- readRDS("memory_KO1.rds")
KO2 <- readRDS("memory_KO2.rds")
WT1 <- readRDS("memory_WT1.rds")
WT2 <- readRDS("memory_WT2.rds")

head(x = KO1[[]])

KO1$origin <- "KO1"
KO2$origin <- "KO2"
WT1$origin <- "WT1"
WT2$origin <- "WT2"

WT1$condition <- "wt"
WT2$condition <- "wt"
KO1$condition <- "ko"
KO2$condition <- "ko"

all_data <- merge(x = KO1, y = c(KO2, WT1, WT2) , add.cell.ids = (c("KO1", "KO2", "WT1", "WT2")) )
immune.combined <- all_data %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = immune.combined@var.genes, npcs = 30, verbose = FALSE)

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = immune.combined, reduction = "pca", pt.size = 0.1, group.by = "origin")
p2 <- VlnPlot(object = immune.combined, features = "PC_1", group.by = "origin", pt.size = 0.1)
plot_grid(p1,p2)

options(repr.plot.height = 2.5, repr.plot.width = 6)
immune.combined <- immune.combined %>%
  RunHarmony("origin", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(immune.combined, 'harmony')
harmony_embeddings[1:5, 1:5]


options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = immune.combined, reduction = "harmony", pt.size = 0.1, group.by = "origin")
p2 <- VlnPlot(object = immune.combined, features = "harmony_1", group.by = "origin", pt.size = 0.1)
plot_grid(p1,p2)

immune.combined <- immune.combined %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.3) %>%
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "origin", pt.size = 0.1)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.1)
plot_grid(p1, p2)

options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(immune.combined, reduction = "umap", split.by = "origin", ncol = 2, label = TRUE)


saveRDS(immune.combined, file = "221108memoryKOWT_integrate_harmony0.3.rds")

immune.combined <- readRDS("221108memoryKOWT_integrate_harmony0.3.rds")
DefaultAssay(immune.combined) <- "RNA"

#Fig.5a

p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 12)
p2

p1 <- DimPlot(immune.combined, reduction = "umap", repel = TRUE)
p1

FeaturePlot(object = immune.combined, features = c("Cd8a", "CD8b"), cols = c("grey", "red"), keep.scale = "feature") & theme(legend.position = "right")

#Fig.5b

FeaturePlot(object = immune.combined, features = c("Zeb2", "Cx3cr1", "Klrg1", "Mki67", 
                                                   "Il7r", "Cxcr3", "Ccr7", "Sell" ), cols = c("grey", "red"), keep.scale = "feature") & theme(legend.position = "right")

#Fig. 5d; cell numbers in each cluster

head(immune.combined)
immune.combined@active.ident <- immune.combined@meta.data$
head(immune.combined@active.ident)
head(immune.combined@meta.data$orig.ident)
tmp = table(immune.combined@active.ident, immune.combined@meta.data$origin)
head(tmp)
write.table(tmp, "Results/221109_res0.3.txt", row.names=T, sep="\t", quote=F)

#top30 & heatmap

all.markers <- FindAllMarkers(object = immune.combined, min.pca = 0.1, thresh.use = 0.1)
top30 <- all.markers %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC) 
write.table(top30, "221108top30.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, 
            col.names = TRUE)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) 
write.table(top10, "221108top10.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, 
            col.names = TRUE)


pdf("221108Heatmap_top10.pdf", width = 60, height = 80) 
DoHeatmap(immune.combined, features = top10$gene) +
  theme(text = element_text(size = 56))
dev.off()


#Fig. 5e; Response to type II IFN
GO_IFNg <- read.table ("GO_IFNg.txt")
GO_IFNg_df <- as.matrix(GO_IFNg)
GO_IFNg_feature <- GO_IFNg_df[, 1]

GO_IFNg_features <- list(c(GO_IFNg_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = GO_IFNg_features,
  ctrl = 5,
  name = 'GO_IFNg_features'
)

head(x = immune.combined[])
VlnPlot(immune.combined, features = 'GO_IFNg_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")


VlnPlot(WT_seu, features = 'GO_IFNg_features1', y.max = 0.8,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 
VlnPlot(KO_seu, features = 'GO_IFNg_features1', y.max = 0.8,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 

#Supple Fig. 5a; bulkRNAseq_down

down <- list(c("Dnaja1", "Il7r", "Sgms1", "Stip1", "Dtx1", "Lman2l", "Arrdc3", "Banp",
               "Ppp2r3a"))
immune.combined <- AddModuleScore(
  object = immune.combined,
  features = down,
  ctrl = 5,
  name = 'down'
)
head(x = immune.combined[])

VlnPlot(immune.combined, features = 'down1', 
        split.by = "condition")

VlnPlot(immune.combined, features = 'down1', split.by = "condition") + stat_compare_means(method = "wilcox.test")
VlnPlot(immune.combined, features = 'down1', split.by = "condition") 

stat_compare_means(method = "t.test", label = "p.signif")

VlnPlot(WT_seu, features = 'down1', y.max = 1,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 
VlnPlot(KO_seu, features = 'down1', y.max = 1,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 

##Supple Fig. 5a; bulkRNAseq_up
up_common1 <- read.table ("UP_bulk.txt")
up_common1_df <- as.matrix(up_common1)
up_common1_feature <- up_common1_df[, 1]
head(up_common1_feature)

up_common_features <- list(c(up_common1_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = up_common_features,
  ctrl = 5,
  name = 'up_features'
)

head(x = immune.combined[])
VlnPlot(immune.combined, features = 'up_features1', 
        split.by = "condition")+ stat_compare_means(method = "wilcox.test")

VlnPlot(WT_seu, features = 'up_features1', y.max = 1,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 
VlnPlot(KO_seu, features = 'up_features1', y.max = 1,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 


#Fig. 5e; cell cycle process
cell_cycle <- read.table ("GO_mitotic.txt")
cell_cycle_df <- as.matrix(cell_cycle)
cell_cycle_feature <- cell_cycle_df[, 1]

cell_cycle_features <- list(c(cell_cycle_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = cell_cycle_features,
  ctrl = 5,
  name = 'cell_cycle_features'
)

head(x = immune.combined[])
VlnPlot(immune.combined, features = 'cell_cycle_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")

VlnPlot(WT_seu, features = 'cell_cycle_features1', y.max = 1.0,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 

VlnPlot(KO_seu, features = 'cell_cycle_features1', y.max = 1.0,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 


#Fig. 5c; Tcm signature
Tcm <- read.table ("230130Tcm_common_list.txt")
Tcm_df <- as.matrix(Tcm)
Tcm_feature <- Tcm_df[, 1]
head(Tcm_feature)
Tcm_features <- list(c(Tcm_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = Tcm_features,
  ctrl = 5,
  name = 'Tcm_features'
)


head(x = immune.combined[])

FeaturePlot(immune.combined,features = 'Tcm_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

VlnPlot(immune.combined, features = 'Tcm_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")

FeaturePlot(WT_seu,features = 'Tcm_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))

FeaturePlot(KO_seu,features = 'Tcm_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))

#Fig.5c; t-TEm/LLEC signature

LLEC <- read.table ("LLECup_vs_OtherMemory.txt")
LLEC_df <- as.matrix(LLEC)
LLEC_feature <- LLEC_df[, 1]

LLEC_features <- list(c(LLEC_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = LLEC_features,
  ctrl = 5,
  name = 'LLEC_features'
)

head(x = immune.combined[])

FeaturePlot(immune.combined,features = 'LLEC_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

VlnPlot(immune.combined, features = 'LLEC_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")

FeaturePlot(WT_seu,features = 'LLEC_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

FeaturePlot(KO_seu,features = 'LLEC_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

VlnPlot(WT_seu, features = 'LLEC_features1', y.max = 0.7,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 
VlnPlot(KO_seu, features = 'LLEC_features1', y.max = 0.7,
        split.by = "condition") + stat_summary(fun = mean, geom='crossbar', width = 0.5, colour = "blue") 

#Fig. 5c; Tcm or Tem signature

TcmTem <- read.table ("230130TCMorTEM_common_list.txt")
TcmTem_df <- as.matrix(TcmTem)
TcmTem_feature <- TcmTem_df[, 1]
head(TcmTem_feature)
TcmTem_features <- list(c(TcmTem_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = TcmTem_features,
  ctrl = 5,
  name = 'TcmTem_features'
)


head(x = immune.combined[])

FeaturePlot(immune.combined,features = 'TcmTem_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

VlnPlot(immune.combined, features = 'TcmTem_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

FeaturePlot(KO_seu,features = 'TcmTem_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")), limits = c(0, 1.5))


FeaturePlot(WT_seu,features = 'TcmTem_features1', ncol = 3
)+ scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")), limits = c(0, 1.5))

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")

VlnPlot(WT_seu, features = 'TcmTem_features1',
        split.by = "condition")+ scale_y_continuous(limits = c(0, 1.5)) + stat_summary(fun = median, geom='crossbar', width = 0.5, colour = "blue") 

VlnPlot(KO_seu, features = 'TcmTem_features1',
        split.by = "condition")+ scale_y_continuous(limits = c(0, 1.5)) + stat_summary(fun = median, geom='crossbar', width = 0.5, colour = "blue") 


#Fig.5f; long-lived memory signature
memory_LL <- read.table ("230124memoryLL_common_list.txt")

memory_LL_df <- as.matrix(memory_LL)
memory_LL_feature <- memory_LL_df[, 1]
head(memory_LL_feature)

memory_LL_features <- list(c(memory_LL_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = memory_LL_features,
  ctrl = 5,
  name = 'memory_LL_features'
)

head(x = immune.combined[])


VlnPlot(immune.combined, features = 'memory_LL_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")


#Fig. 5f; Negative regulation of apoptosis in long-lived memory signature

NegRegApo <- read.table ("230317RegApo_common_list.txt")
NegRegApo_df <- as.matrix(NegRegApo)
NegRegApo_feature <- NegRegApo_df[, 1]
head(NegRegApo_feature)
NegRegApo_features <- list(c(NegRegApo_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = NegRegApo_features,
  ctrl = 5,
  name = 'NegRegApo_features'
)


head(x = immune.combined[])


VlnPlot(immune.combined, features = 'NegRegApo_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

#Fig. 5f;TE (short-lived) signature

TE <- read.table ("230130TE_common_list.txt")
TE_df <- as.matrix(TE)
TE_feature <- TE_df[, 1]
head(TE_feature)
TE_features <- list(c(TE_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = TE_features,
  ctrl = 5,
  name = 'TE_features'
)


head(x = immune.combined[])

FeaturePlot(immune.combined,features = 'TE_features1', ncol = 3
) + scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) 

VlnPlot(immune.combined, features = 'TE_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

WT_seu <- subset(immune.combined, subset = condition == "wt")
KO_seu <- subset(immune.combined, subset = condition == "ko")

VlnPlot(WT_seu, features = 'TE_features1',
        split.by = "condition")+ scale_y_continuous(limits = c(-0.3, 1.3)) + stat_summary(fun = median, geom='crossbar', width = 0.5, colour = "blue") 

VlnPlot(KO_seu, features = 'TE_features1',
        split.by = "condition")+ scale_y_continuous(limits = c(-0.3, 1.3)) + stat_summary(fun = median, geom='crossbar', width = 0.5, colour = "blue") 



#Supple Fig.5c; cell cycle process-related TE signature

TECycle <- read.table ("230317TE&cellcycle_common_list.txt")
TECycle_df <- as.matrix(TECycle)
TECycle_feature <- TECycle_df[, 1]
head(TECycle_feature)
TECycle_features <- list(c(TECycle_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = TECycle_features,
  ctrl = 5,
  name = 'TECycle_features'
)


head(x = immune.combined[])


VlnPlot(immune.combined, features = 'TECycle_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

#Supple Fig.5c: cell cycle process-unrelated TE signature

TEbutCycle <- read.table ("230317TEbut_common_list.txt")
TEbutCycle_df <- as.matrix(TEbutCycle)
TEbutCycle_feature <- TEbutCycle_df[, 1]
head(TEbutCycle_feature)
TEbutCycle_features <- list(c(TEbutCycle_df[, 1]))

immune.combined <- AddModuleScore(
  object = immune.combined,
  features = TEbutCycle_features,
  ctrl = 5,
  name = 'TEbutCycle_features'
)


head(x = immune.combined[])


VlnPlot(immune.combined, features = 'TEbutCycle_features1', 
        split.by = "condition") + stat_compare_means(method = "wilcox.test")



#Supple Fig.5b; stat1, il7r

VlnPlot(immune.combined, features = 'Stat1', y.max = 4,
        split.by = "condition") + stat_compare_means(method = "wilcox.test")

VlnPlot(immune.combined, features = 'Il7r', 
        split.by = "condition")+ stat_compare_means(method = "wilcox.test")



