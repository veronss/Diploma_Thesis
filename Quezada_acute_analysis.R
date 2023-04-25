# title: "Quezada et al. (2023) re-analysis of acute infection data"
# author: "Veronika Cimermanova"
# date: "2023-04-21"

# import packages  
library(Seurat)
library(ggplot2)
library(harmony)
library(tidyverse)
library(SingleR) # must be version 1.0.1

# For installation of the SingleR package, use:
# remotes::install_github("dviraran/SingleR", upgrade = "never")

#### import of downloaded data, QC and creation of Seurat object ####
# data import from the link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213470
# please create your own path below for downloaded files
GSE213470_UMI4Wenhao_Cl13_full_table <- read.csv("~/GSE213470_UMI4Wenhao_Cl13_full_table.csv", row.names=1)
GSE213470_Metadata4Wenhao_Cl13 <- read.csv("~/Quezada/GSE213470_Metadata4Wenhao_Cl13.csv", row.names=1)

# transpose counts
data <- t(GSE213470_UMI4Wenhao_Cl13_full_table)

# QC1: Remove mitochondrial, ribosomal, TRAV, TRBV and unannotated genes
data <- data[!(grepl('Trav', rownames(data)) |
                 grepl('Trbv', rownames(data)) |
                 grepl('mt-', rownames(data)) |
                 grepl('Rpl', rownames(data)) |
                 grepl('Rps', rownames(data)) |
                 grepl('MGP-', rownames(data))),]

# creation of Seurat object and QC2: removing genes detected in less than 3 cells and cells with less than 200 genes
seurat.obj <- CreateSeuratObject(counts = data, project = "Quezada", min.cells = 3, min.features = 200, meta.data = GSE213470_Metadata4Wenhao_Cl13)
View(seurat.obj@meta.data) # check

# check of features(number of genes) and counts (number of transcripts)
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

# QC3: filtering out low quality cells with many transcripts
seurat.obj <- subset(seurat.obj, subset = nCount_RNA < 20000)

# Extraction of only acute infection data in the seurat object
seurat_Arm <- subset(seurat.obj, virus == "LCMV_Arm")

#### Data normalisation, scaling, variable features, linear dimensionality reduction and clustering ####
seurat_Arm <- NormalizeData(seurat_Arm) # normalisation
seurat_Arm <- FindVariableFeatures(seurat_Arm) # variable features
all.genes <- rownames(seurat_Arm)
seurat_Arm <- ScaleData(seurat_Arm, features = all.genes) # scaling
seurat_Arm <- RunPCA(seurat_Arm, features = VariableFeatures(object = seurat_Arm)) # linear dimensionality reduction
ElbowPlot(seurat_Arm) # check of dimensionalit, Result: working with only principal components 1:15
seurat_Arm <- FindNeighbors(seurat_Arm, dims = 1:15)
seurat_Arm <- FindClusters(seurat_Arm, resolution = 0.1)
seurat_Arm <- RunUMAP(seurat_Arm, dims = 1:15)

# check of clustering
DimPlot(seurat_Arm, reduction = "umap")
DimPlot(seurat_Arm, reduction = "umap", label=T, group.by = "timepoint")

#### SingleR to identify cells for contamination and specification ####
counts <- GetAssayData(seurat_Arm[["RNA"]])
seu_singler <- CreateSinglerObject(counts=counts,
                                   project.name= "Quezada", # choose
                                   min.genes = 200, # ignore cells with fewer than 200 transcripts
                                   technology = "10x", # choose
                                   species = "Mouse",
                                   fine.tune = FALSE, # TRUE would take very long
                                   reduce.file.size = TRUE, # leave out less-often used fields 
                                   do.signatures = FALSE,
                                   do.main.types = TRUE,
                                   numCores = 4)

# 4 different SingleR annotations to clusters
seurat_Arm <- AddMetaData(seurat_Arm,
                          seu_singler$singler[[1]]$SingleR.single$labels,
                          col.name = "Immgen_annot_single")
seurat_Arm <- AddMetaData(seurat_Arm,
                          seu_singler$singler[[1]]$SingleR.single.main$labels,
                          col.name = "Immgen_annot_single_main")
seurat_Arm <- AddMetaData(seurat_Arm,
                          seu_singler$singler[[2]]$SingleR.single.main$labels,
                          col.name = "MouseRNAseq_single_main")
seurat_Arm <- AddMetaData(seurat_Arm,
                          seu_singler$singler[[2]]$SingleR.single$labels,
                          col.name = "MouseRNAseq_single")

# check of cell types assign to choose by which to filter
DimPlot(seurat_Arm, group.by = "MouseRNAseq_single_main") +
  DimPlot(seurat_Arm, group.by = "MouseRNAseq_single")+
  DimPlot(seurat_Arm, group.by = "Immgen_annot_single_main")

# filtering for pure T cells and re-clustering
seurat_Arm_T <- subset(seurat_Arm, Immgen_annot_single_main %in% c("T cells"))

all.genes.T <- rownames(seurat_Arm_T)
seurat_Arm_T <- ScaleData(seurat_Arm_T, features = all.genes.T)
seurat_Arm_T <- FindVariableFeatures(seurat_Arm_T)
seurat_Arm_T <- RunPCA(seurat_Arm_T, features = VariableFeatures(object = seurat_Arm_T))
seurat_Arm_T <- FindNeighbors(seurat_Arm_T, dims = 1:15)
seurat_Arm_T <- FindClusters(seurat_Arm_T, resolution = 0.1)
seurat_Arm_T <- RunUMAP(seurat_Arm_T, dims = 1:15)

# check of clustering
DimPlot(seurat_Arm_T, reduction = "umap")
DimPlot(seurat_Arm, reduction = "umap", label=T, group.by = "timepoint")

#### Integration by Harmony to reduce batch effect of timepoint collection ####
seurat_Arm_integrated <- seurat_Arm_T %>% RunHarmony(group.by.vars = "timepoint")
seurat_Arm_harmony_emb <- Embeddings(seurat_Arm_integrated, "harmony") #  harmony dimensions

# re-clustering
seurat_Arm_integrated <- seurat_Arm_integrated %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.3)

#### Final ATLAS ####
# Thesis = "Figure 12A"
DimPlot(seurat_Arm_integrated, reduction = "umap") # 6 clusters

#### Exploration of Atlas ####
## naive, effector, memory markers ##
# Thesis = "Figure 12B"
FeaturePlot(seurat_Arm_integrated, features = c("Il7r","Sell","Ccr7", "Cxcr3", "Mki67", "Klrg1", "Cx3cr1","Gzmb","Havcr2"), min.cutoff = 0)

## timepoints exploration ##
# Thesis = "Figure 12C"
my_levels <- c("d0","div1","d3","d5","d6","d7","d8","d22","d60") # order of timepoints
seurat_Arm_integrated@meta.data$timepoint <- factor(x = seurat_Arm_integrated@meta.data$timepoint, levels = my_levels)
DimPlot(seurat_Arm_integrated, reduction = "umap", label=F, group.by = "timepoint") 

# Dimplot to show cells only from the first division by creating new column to mark only first division cells
# Thesis = "Figure 12D"
seurat_Arm_integrated@meta.data$FirstD <- NA
seurat_Arm_integrated@meta.data$FirstD[seurat_Arm_integrated@meta.data$timepoint == "div1"] <- "div1"
seurat_Arm_integrated@meta.data$FirstD[seurat_Arm_integrated@meta.data$timepoint != "div1"] <- "Uns"
DimPlot(seurat_Arm_integrated, reduction = "umap", group.by = "FirstD", cols=c("darkorange1","gray80"))

# clusters quantifications based on timepoints
# Thesis = "Figure 12D"
cq <- seurat_Arm_integrated@meta.data %>% group_by(seurat_clusters, timepoint) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

plot <- ggplot(data=cq, aes(x=seurat_clusters, y=freq*100, fill=timepoint))+
  geom_bar(stat="identity")+
  theme_classic()+
  ylim(0,100)+
  ylab("cells frequency %")+
  labs(fill = "timepoints")+
  xlab("cluster")+
  theme(axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.title = element_text(size=20,),
        legend.text = element_text(size=20))
print(plot)

## Tpam signature ##
# Thesis = "Figure 12E"
# manual creation of the signature
Tpam_signature <- stringr::str_to_sentence(c("Slc7a5", "Slc3a2", "Odc1", "Srm", "Sms", "Eif5a", "Amd1", "Dhps", "Dohh", "Tnfrsf9", "Cd69","Myc"))

# calculation of Tpam signature score
Tpam_sig <- AddModuleScore(
  object = seurat_Arm_integrated,
  features = list(c(Tpam_signature)),
  search = F,
  ctrl = 50,
  nbin = 50,
  assay = "RNA",
  name = 'Tpam_signature_')

# visualisation of Tpam signature score
FeaturePlot(Tpam_sig, features = "Tpam_signature_1", min.cutoff = 0) # on Atlas
VlnPlot(Tpam_sig, features = "Tpam_signature_1", pt.size = 0) # in clusters
VlnPlot(Tpam_sig, features = "Tpam_signature_1", pt.size = 0, group.by = "timepoint") # in timepoint collections


#### session info ####
sessionInfo()


