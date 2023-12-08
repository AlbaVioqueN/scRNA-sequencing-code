#HOW TO MERGE 2 DATASETS: datadietwt1 and datadietlowiod1
library(Seurat)
library(tidyverse)
library(hexbin)
library(patchwork)
library(RSQLite)
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
install.packages('BiocManager')
BiocManager::install('limma')
install.packages("openxlsx")
library(openxlsx)


# First, establish the working directory 
getwd()
setwd("~/Desktop/1+3 PhD MRC KCL/Year 1/Rotation 1/DataDiet")

# Name the datasets - since they are .rds files, the following command is needed:
datacontrol1 <- readRDS("/Users/albavioque/Desktop/1+3 PhD MRC KCL/Year 1/Rotation 1/DataDiet/MMChow-wt1/MMchow1.rds")
datahypothyroidism1 <- readRDS("/Users/albavioque/Desktop/1+3 PhD MRC KCL/Year 1/Rotation 1/DataDiet/MMPTU-Cont2/MMPTUCont2.rds")

# Then load the raw data of each dataset
raw_datacontrol1 <- as.matrix(GetAssayData(datacontrol1, slot = "counts")[, WhichCells(datacontrol1, slot = "counts")])
raw_datahypothyroidism1<- as.matrix(GetAssayData(datahypothyroidism1, slot = "counts")[, WhichCells(datahypothyroidism1, slot = "counts")])

#Now we transform the matrix to objects
control <- CreateSeuratObject(counts = raw_datacontrol1, project = "CONTROL")
hypothyroidism <- CreateSeuratObject(counts = raw_datahypothyroidism1, project = "HYPOTHYROIDISM")
rm(raw_datahypothyroidism1)

# To confirm that both are seurat objects before proceeding merging the datasets:
control
hypothyroidism

# Before merging we need to use split.by to separate them
control$info <- "Control"
hypothyroidism$info <- "Hypothyroidism"

# Now we can marge both datasets
control1_lowiod1_merged <- merge(control, y = hypothyroidism, add.cell.ids = c("c", "m"), project = "Merged_1")

# Check again in case we want to filter again 
control1_lowiod1_merged[["percent.mt"]] <- PercentageFeatureSet(control1_lowiod1_merged, pattern = "^mt-")
view(control1_lowiod1_merged)

VlnPlot(control1_lowiod1_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

# Normalize the data
control1_lowiod1_merged <- NormalizeData(control1_lowiod1_merged)

# To find the most variable genes (the selection.method and nfeatures stay the same)
control1_lowiod1_merged <- FindVariableFeatures(control1_lowiod1_merged, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(control1_lowiod1_merged)

#Not sure which one to use of the following scaledata:
control1_lowiod1_merged <- ScaleData(control1_lowiod1_merged, features = all.genes) 
control1_lowiod1_merged <- RunPCA(control1_lowiod1_merged, features = VariableFeatures(control1_lowiod1_merged))

# Use Elbowplot to check where to cut the dimensions 
ElbowPlot(control1_lowiod1_merged, ndims = 25)
ggsave("Elbowplot.tiff", plot = last_plot(), device = NULL, width = 10, height = 7)

# Start clustering cells
control1_lowiod1_merged <- FindNeighbors(control1_lowiod1_merged, dims = 1:15)
saveRDS(object = control1_lowiod1_merged, "control1_lowiod1_merged.rds")

table(control1_lowiod1_merged$orig.ident)

# Miriam tried in her data with resolution 0.3, 0.4  0.5, but the best was 0.3
control1_lowiod1_merged_0.3 <- FindClusters(control1_lowiod1_merged, resolution = 0.3)
control1_lowiod1_merged_0.3 <- RunUMAP(control1_lowiod1_merged_0.3, dims = 1:15)
DimPlot(control1_lowiod1_merged_0.3, label = T, split.by = 'info', repel = T)

# Now we can find most upregulated genes 
control1_lowiod1_merged_0.3_markers <- FindAllMarkers(control1_lowiod1_merged_0.3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(object = control1_lowiod1_merged_0.3_markers, "control1_lowiod1_merged_0.3_markers.rds")

# Now we can check the expression of the markers of interest in each cluster, for this we'll use a ViolinPlot and FeaturePlot to see where their expression is the most enriched. 

DimPlot(control1_lowiod1_merged_0.3, label = TRUE, split.by = "info")
ggsave("merge1_labelnumbers.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

# STEM CELLS - Cluster 10
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Rbpms", "Sox9", "Mia"), split.by = "info")
ggsave("merge1_stemcells1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Sox2", "Grhl2"), split.by = "info")
ggsave("merge1_stemcells2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Rbpms", "Sox9", "Sox2", "Grhl2"), split.by = "info")
ggsave("merge1_stemcells_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#SOMATOTROPHS - Cluster 1
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Gh", "Pappa2"), split.by = "info") 
ggsave("merge1_somatotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Mmp7", "Ghrhr", "Grhl2"), split.by = "info") 
ggsave("merge1_somatotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c ("Kcnmb2", "Pou1f1"), split.by = "info") 
ggsave("merge1_somatotrophs3_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Gh", "Pappa2", "Mmp7", "Ghrhr", "Kcnmb2", "Grhl2", "Pou1f1"), split.by = "info") 
ggsave("merge1_somatotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#LACTOTROPHS - Cluster 2
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Prl", "Gpr83"), split.by = "info")
ggsave("merge1_Lactotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Edil3", "Pou1f1"), split.by = "info")
ggsave("merge1_Lactotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Prl", "Gpr83", "Edil3", "Pou1f1"), split.by = "info")
ggsave("merge1_Lactotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#THYROTROPHS - Cluster 6
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Tshb", "Trhr"), split.by = "info")
ggsave("merge1_thyrotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Dio2", "Pou1f1"), split.by = "info")
ggsave("merge1_thyrotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Tshb", "Trhr", "Dio2", "Pou1f1"), split.by = "info")
ggsave("merge1_thyrotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#MELANOTROPHS - Cluster 8
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Pomc", "Pax7"), split.by = "info")
ggsave("merge1_Melanotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Gulo", "Tbx19"), split.by = "info")
ggsave("merge1_Melanotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Pomc", "Pax7", "Gulo", "Tbx19"), split.by = "info")
ggsave("merge1_Melanotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#CORTICOTROPHS - Cluster 4
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Pomc", "Gpc5"), split.by = "info")
ggsave("merge1_Corticotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Tbx19", "Pax6"), split.by = "info")
ggsave("merge1_Corticotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Pomc", "Gpc5", "Tbx19", "Pax6"), split.by = "info")
ggsave("merge1_Corticotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#GONADOTROPHS - Cluster 7
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Lhb", "Fshb", "Foxp2"), split.by = "info")
ggsave("merge1_Gonadotrophs1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Gnrhr", "Cga"), split.by = "info")
ggsave("merge1_Gonadotrophs2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Tgfbr3l", "Nr5a1"), split.by = "info")
ggsave("merge1_Gonadotrophs3_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Lhb", "Fshb", "Foxp2", "Gnrhr", "Cga", "Tgfbr3l", "Nr5a1"), split.by = "info")
ggsave("merge1_Gonadotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#PITUICYTES - Cluster 16
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Sox2", "Nkx2-1", "Col25a1"), split.by = "info")
ggsave("merge1_Pituicytes_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Sox2", "Nkx2-1", "Col25a1"), split.by = "info")
ggsave("merge1_Pituicytes_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#ENDOTHELIAL CELLS - Cluster 15
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Pecam1", "Plvap", "Emcn"), split.by = "info")
ggsave("merge1_Endothelial1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Esam", "Adgrl4", "Stc1"), split.by = "info")
ggsave("merge1_Endothelial2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Pecam1", "Plvap", "Emcn", "Esam", "Adgrl4", "Stc1"), split.by = "info")
ggsave("merge1_Endothelial_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#MACROPHAGES - Cluster 14
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Cd163", "C1qa", "C1qc", "C1qb"), split.by = "info")
ggsave("merge1_macrophages1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Ms4a7", "Csf1r", "Cd86"), split.by = "info")
ggsave("merge1_macrophages2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Cd163", "C1qa", "C1qc", "C1qb", "Ms4a7", "Csf1r", "Cd86"), split.by = "info")
ggsave("merge1_macrophages_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#PERYCITES - Cluster 13
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Mgp", "Dcn", "Apod", "S1pr3"), split.by = "info")
ggsave("merge1_pericytes1_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(control1_lowiod1_merged_0.3, features = c("Tagln", "Pdgfra", "Bgn"), split.by = "info")
ggsave("merge1_pericytes2_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

VlnPlot(control1_lowiod1_merged_0.3, features = c("Mgp", "Dcn", "Apod", "S1pr3", "Tagln", "Pdgfra", "Bgn"), split.by = "info")
ggsave("merge1_pericytes1_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#Cluster 0: Somatotrophs 2
#Cluster 1: Somatotrophs 1
#Cluster 2: Lactotrophs 1
#Cluster 3: Lactotrophs 2
#Cluster 4: Corticotrophs
#Cluster 5: Lactotrophs 3
#Cluster 6: Thyrotrophs
#Cluster 7: Gonadotrophs
#Cluster 8: Melanotrophs 1
#Cluster 9: Melanotrophs 2
#Cluster 10: Stem cells
#Cluster 11: Lactotrophs 4
#Cluster 12: Corticotrophs + Melanotrophs
#Cluster 13: Pericytes
#Cluster 14: Macrophages 
#Cluster 15: Endothelial cells
#Cluster 16: Pituicytes

# Now that we know which cell type is each cluster, we can label them (by order):
new.cluster.ids <- c("Somatotrophs 2", "Somatotrophs 1", "Lactotrophs 1", "Lactotrophs 2", "Corticotrophs", "Lactotrophs 3", "Thyrotrophs", "Gonadotrophs", "Melanotrophs 1", "Melanotrophs 2", "Stem cells", "Lactotrophs 4", "Corticotrophs + Melanotrophs", "Pericytes", "Macrophages", "Endothelial cells", "Pituicytes")

names(new.cluster.ids) <- levels(control1_lowiod1_merged_0.3)
control1_lowiod1_merged_0.3 <- RenameIdents(control1_lowiod1_merged_0.3, new.cluster.ids)

DimPlot(control1_lowiod1_merged_0.3, label = T, pt.size = 0.6, repel = T, split.by = 'info') + NoLegend()  
ggsave("Merge1FINALUMAP_nolegend.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
DimPlot(control1_lowiod1_merged_0.3, pt.size = 0.6, split.by = 'info') 
ggsave("Merge1FINALUMAP_withlegend.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#SUBSET STEM CELLS AND THYROTROPHS CLUSTER for later, and before adding the tags to each cluster name
#Stem cells
stemcell_subset_1 <- subset(control1_lowiod1_merged_0.3, idents = "Stem cells")
stemcell_subset_1 <- NormalizeData(stemcell_subset_1)
stemcell_subset_1 <- ScaleData(stemcell_subset_1)
DimPlot(stemcell_subset_1, split.by = 'info')

#Thyrotrophs
thyrotrophs_subset_1 <- subset(control1_lowiod1_merged_0.3, idents = "Thyrotrophs")
thyrotrophs_subset_1 <- NormalizeData(thyrotrophs_subset_1)
thyrotrophs_subset_1 <- ScaleData(thyrotrophs_subset_1)
DimPlot(thyrotrophs_subset_1, split.by = 'info')


#CHECK CURRENT INDELT AND CHANGE IT TO HAVE NEW TAGS
control1_lowiod1_merged_0.3_labels <- control1_lowiod1_merged_0.3
control1_lowiod1_merged_0.3_labels@meta.data$new
control1_lowiod1_merged_0.3_labels@active.ident
#Now rename of the clusters and add the _Control or _Hypothyroidism tag
control1_lowiod1_merged_0.3_labels$different <- paste(Idents(control1_lowiod1_merged_0.3_labels), control1_lowiod1_merged_0.3_labels$info, sep = "_")
control1_lowiod1_merged_0.3_labels <- SetIdent(control1_lowiod1_merged_0.3_labels, value = control1_lowiod1_merged_0.3_labels@meta.data$different)
tail(control1_lowiod1_merged_0.3_labels@active.ident)


#CHECK STEM CELLS FOR DIFFERENTIALLY EXPRESSED MARKERS
# Stem cells vs other cells
stem.cell.markers <- FindMarkers(control1_lowiod1_merged_0.3_labels, ident.1 = c("Stem cells_Hypothyroidism"), only.pos = TRUE, min.pct = 0.3) #First the Low_iodine and second the control!!!
view(stem.cell.markers)
write.csv(stem.cell.markers, "stemcellLowIodineupregulated.csv", row.names = TRUE)

# UPREGULATED - Stem cells hypothyroidism vs stem cell control
stem.cell.lowiodine.markers <- FindMarkers(control1_lowiod1_merged_0.3_labels, ident.1 = "Stem cells_Hypothyroidism", ident.2 ="Stem cells_Control", only.pos = TRUE, min.pct = 0.30)
head(stem.cell.lowiodine.markers, 10)

#To save it, first we need to convert it into a dataframe (so we can see the column with the gene names):
stem.cell.lowiodine.markers_dataframe <- as.data.frame(stem.cell.lowiodine.markers)
stem.cell.lowiodine.markers_dataframe$Gene <- rownames(stem.cell.lowiodine.markers_dataframe)
stem.cell.lowiodine.markers_dataframe <- stem.cell.lowiodine.markers_dataframe[, c("Gene", setdiff(names(stem.cell.lowiodine.markers_dataframe), "Gene"))]
write.xlsx(stem.cell.lowiodine.markers_dataframe, "stemcell_lowiodine_dataframe_upregulated.xlsx")

view(stem.cell.lowiodine.markers_dataframe)

VlnPlot(control1_lowiod1_merged_0.3, feature = "Pomc", split.by = "info") #To make sure that the list is correct
FeaturePlot(control1_lowiod1_merged_0.3, feature = "Pomc", split.by = "info") #To make sure that the list is correct 2

# DOWNREGULATED - Stem cells hypothyroidism vs stem cell control
stem.cell.lowiodine.markers_inversed <- FindMarkers(control1_lowiod1_merged_0.3, ident.1 = "Stem cells_Control", ident.2 ="Stem cells_Low_iodine", only.pos = TRUE, min.pct = 0.30)
head(stem.cell.lowiodine.markers_inversed, 20)

#To save it, first we need to convert it into a dataframe (so we can see the column with the gene names):
stem.cell.lowiodine.markers_inversed_dataframe <- as.data.frame(stem.cell.lowiodine.markers_inversed)
stem.cell.lowiodine.markers_inversed_dataframe$Gene <- rownames(stem.cell.lowiodine.markers_inversed_dataframe)
stem.cell.lowiodine.markers_inversed_dataframe <- stem.cell.lowiodine.markers_inversed_dataframe[, c("Gene", setdiff(names(stem.cell.lowiodine.markers_inversed_dataframe), "Gene"))]
write.xlsx(stem.cell.lowiodine.markers_inversed_dataframe, "stemcell_lowiodine_inversed_dataframe_upregulated.xlsx")

view(stem.cell.lowiodine.markers_inversed_dataframe)

saveRDS(object = control1_lowiod1_merged_0.3, "control1_lowiod1_merged_0.3.rds")
saveRDS(object = control1_lowiod1_merged_0.3_labels, "control1_lowiod1_merged_0.3_labels.rds")
saveRDS(object = stem.cell.lowiodine.markers, "stem.cell.lowiodine.markers.rds")
saveRDS(object = stem.cell.lowiodine.markers_inversed, "stem.cell.lowiodine.markers_inversed.rds")
saveRDS(object = stemcell_subset_1, "stemcell_subset_1.rds")
saveRDS(object = thyrotropes_subset_1, "thyrotropes_subset_1.rds")


#Now we want to do a FeaturePlot and a VlnPlot with the most differentiated expressed genes (coding for secreted factors) in STEM CELLS population:

#UPREGULATED GENES:
#Tshb
VlnPlot(stemcell_subset_1, features =  "Tshb", split.by = "info")
ggsave("merge1_Tshb_upreg_sc_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Tshb", split.by = "info")
ggsave("merge1_Tshb_upreg_sc_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Pomc
VlnPlot(stemcell_subset_1, features =  "Pomc", split.by = "info")
ggsave("merge1_Pomc_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Pomc", split.by = "info")
ggsave("merge1_Pomc_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Cga
VlnPlot(stemcell_subset_1, features =  "Cga", split.by = "info")
ggsave("merge1_Cga_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Cga", split.by = "info")
ggsave("merge1_Cga_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Dlk1
VlnPlot(stemcell_subset_1, features =  "Dlk1", split.by = "info")
ggsave("merge1_Dlk1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Dlk1", split.by = "info")
ggsave("merge1_Dlk1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

FeaturePlot(control1_lowiod1_merged_0.3, feature =  "Dlk1", split.by = "info")
ggsave("merge1_Dlk1_upreg_merge_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Sdk1
VlnPlot(stemcell_subset_1, features =  "Sdk1", split.by = "info")
ggsave("merge1_Sdk1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Sdk1", split.by = "info")
ggsave("merge1_Sdk1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Clu
VlnPlot(stemcell_subset_1, features =  "Clu", split.by = "info")
ggsave("merge1_Clu_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Clu", split.by = "info")
ggsave("merge1_Clu_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Dscam
VlnPlot(stemcell_subset_1, features =  "Dscam", split.by = "info")
ggsave("merge1_Dscam_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Dscam", split.by = "info")
ggsave("merge1_Dscam_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Thsd4
VlnPlot(stemcell_subset_1, features =  "Thsd4", split.by = "info")
ggsave("merge1_Thsd4_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Thsd4", split.by = "info")
ggsave("merge1_Thsd4_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Mia
VlnPlot(stemcell_subset_1, features =  "Mia", split.by = "info")
ggsave("merge1_Mia_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Mia", split.by = "info")
ggsave("merge1_Mia_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Calm1
VlnPlot(stemcell_subset_1, features =  "Calm1", split.by = "info")
ggsave("merge1_Calm1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Calm1", split.by = "info")
ggsave("merge1_Calm1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Bsg
VlnPlot(stemcell_subset_1, features =  "Bsg", split.by = "info")
ggsave("merge1_Bsg_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Bsg", split.by = "info")
ggsave("merge1_Bsg_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Tenm4
VlnPlot(stemcell_subset_1, features =  "Tenm4", split.by = "info")
ggsave("merge1_Tenm4_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Tenm4", split.by = "info")
ggsave("merge1_Tenm4_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Fat3
VlnPlot(stemcell_subset_1, features =  "Fat3", split.by = "info")
ggsave("merge1_Fat3_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Fat3", split.by = "info")
ggsave("merge1_Fat3_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Cadm1
VlnPlot(stemcell_subset_1, features =  "Cadm1", split.by = "info")
ggsave("merge1_Cadm1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Cadm1", split.by = "info")
ggsave("merge1_Cadm1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Efna5
VlnPlot(stemcell_subset_1, features =  "Efna5", split.by = "info")
ggsave("merge1_Efna5_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Efna5", split.by = "info")
ggsave("merge1_Efna5_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Mdk
VlnPlot(stemcell_subset_1, features =  "Mdk", split.by = "info")
ggsave("merge1_Mdk_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Mdk", split.by = "info")
ggsave("merge1_Mdk_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Crlf1
VlnPlot(stemcell_subset_1, features =  "Crlf1", split.by = "info")
ggsave("merge1_Crlf1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Crlf1", split.by = "info")
ggsave("merge1_Crlf1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Plxnb2
VlnPlot(stemcell_subset_1, features =  "Plxnb2", split.by = "info")
ggsave("merge1_Plxnb2_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Plxnb2", split.by = "info")
ggsave("merge1_Plxnb2_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Glg1
VlnPlot(stemcell_subset_1, features =  "Glg1", split.by = "info")
ggsave("merge1_Glg1_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Glg1", split.by = "info")
ggsave("merge1_Glg1_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Pcdh9
VlnPlot(stemcell_subset_1, features =  "Pcdh9", split.by = "info")
ggsave("merge1_Pcdh9_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Pcdh9", split.by = "info")
ggsave("merge1_Pcdh9_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Psap
VlnPlot(stemcell_subset_1, features =  "Psap", split.by = "info")
ggsave("merge1_Psap_upreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Psap", split.by = "info")
ggsave("merge1_Psap_upreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
getwd()


#DOWNREGULATED GENES:
#Gh
VlnPlot(stemcell_subset_1, features =  "Gh", split.by = "info")
ggsave("merge1_Gh_downreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Gh", split.by = "info")
ggsave("merge1_Gh_downreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Prl
VlnPlot(stemcell_subset_1, features =  "Prl", split.by = "info")
ggsave("merge1_Prl_downreg_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(stemcell_subset_1, feature =  "Prl", split.by = "info")
ggsave("merge1_Prl_downreg_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#Final plot for the most upregulated and downregulated genes coding for secreted factors
DotPlot(stemcell_subset_1, features = c("Tshb", "Cga", "Dlk1", "Sdk1", "Clu", "Dscam", "Thsd4", "Bsg", "Tenm4", "Gh", "Prl"), cols = c("magenta", "blue"), split.by = "info")
ggsave("merge1_upreg_sc_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

DotPlot(stemcell_subset_1, features = c("Tshb", "Cga", "Dlk1", "Sdk1", "Clu", "Dscam", "Thsd4", "Bsg", "Tenm4", "Gh", "Prl"), cols = c("magenta", "blue"), split.by = "info") + RotatedAxis() + CoordFlip
ggsave("merge1_upreg_sc_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 4)

#CHECK WNT SIGNALLING PATHWAY
VlnPlot(stemcell_subset_1, features =  c("Wnt1", "Wnt2", "Wnt3", "Wnt4", "Wnt5", "Wnt6", "Wnt7", "Wnt8", "Wnt9", "Wnt10", "Wnt11", "Wnt12", "Wnt13", "Wnt14", "Wnt15", "Wnt16", "Wnt17", "Wnt18", "Wnt19", "Axin2", "Lef1", "Wls"), split.by = "info")
VlnPlot(stemcell_subset_1, features =  c("Axin2", "Lef1", "Wls"), split.by = "info")

#CHECK THYROTROPES FOR DIFFERENTIALLY EXPRESSED MARKERS 
#Thyrotrophs
thyrotrophs.markers <- FindMarkers(control1_lowiod1_merged_0.3, ident.1 = c("Thyrotrophs_Low_iodine"), only.pos = TRUE, min.pct = 0.3) 
view(thyrotrophs.markers)
write.csv(thyrotrophs.markers, "ThyrotrophsLowIodineupregulated.csv", row.names = TRUE)

#low iodine vs control
thyrotrophs.lowiodine.markers <- FindMarkers(control1_lowiod1_merged_0.3, ident.1 = "Thyrotrophs_Low_iodine", ident.2 ="Thyrotrophs_Control", only.pos = TRUE, min.pct = 0.30) #First the hypothyroid sample and second the control!!!
head(thyrotrophs.lowiodine.markers, 180)
write.csv(thyrotrophs.lowiodine.markers, "ThyrotrophsLowIodineupregulated.csv", row.names = TRUE)

#Now we want to check the expression of Tshb in thyrotrophs - it was the most differentially expressed in stem cells
VlnPlot(thyrotrophs_subset_1, feature = "Tshb", split.by = "info") #To make sure that the list is correct
ggsave("Tshb_thyrotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(thyrotrophs_subset_1, feature = "Tshb", split.by = "info") #To make sure that the list is correct 2
ggsave("Tshb_thyrotrophs_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Check Tshb and MiKi67 in thyrotrophs
VlnPlot(thyrotrophs_subset_1, feature = c("Tshb", "Mki67"), split.by = "info") 
ggsave("merge1_Tshb_Mki67_thyrotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(thyrotrophs_subset_1, feature = c("Tshb", "Mki67"), split.by = "info") 
ggsave("merge1_Tshb_Mki67_thyrotrophs_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)

#Check Tshr 
VlnPlot(thyrotrophs_subset_1, feature = "Tshr", split.by = "info") 
ggsave("merge1_Tshr_thyrotrophs_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)
FeaturePlot(thyrotrophs_subset_1, feature = "Tshr", split.by = "info") 
ggsave("merge1_Tshr_thyrotrophs_FeaturePlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


#Finally, we want to check the expression of Tshb, Cga, and Trhr in the all the pituitary cell types
VlnPlot(control1_lowiod1_merged_0.3, features = c("Tshb", "Cga", "Trhr"), split.by = "info", stack = TRUE, cols = "red", "blue")
ggsave("merge1_Tshb_Cga_Trhr_all_VlnPlot.tiff", plot = last_plot(), device = NULL, width = 10, height = 6)


