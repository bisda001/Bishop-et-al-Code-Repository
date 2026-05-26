
# Packages ----------------------------------------------------------------

library(tidyverse)
library(Seurat)

#renv::install("SeuratObject", rebuild = TRUE)

library(CellChat)
library(readr)


# Setup -------------------------------------------------------------------

figure_out_dir <- "~/260430-final_scRNA_analysis-output"
dir.create(figure_out_dir, showWarnings = FALSE)



# Processing --------------------------------------------------------------


Kasper.scaledata <- read_rds(file = "/media/danielbishop/T7/scRNA/Kasper.scaledata.Seurat_v4.rds")

e13.kasper.2023 <- subset(Kasper.scaledata, subset = embryonic_age == "E13.5")
rm(Kasper.scaledata)
gc()

e13.kasper.2023 <- NormalizeData(e13.kasper.2023, normalization.method = "LogNormalize", scale.factor = 10000)

e13.kasper.2023 <- FindVariableFeatures(e13.kasper.2023, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(e13.kasper.2023), 10)
top10


plot1 <- VariableFeaturePlot(e13.kasper.2023)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(e13.kasper.2023)
e13.kasper.2023 <- ScaleData(e13.kasper.2023, features = all.genes, vars.to.regress = c("sample_date", "sex", "perc_mito", "S.Score", "G2M.Score", "nCount_RNA"))


e13.kasper.2023 <- RunPCA(e13.kasper.2023, features = VariableFeatures(object = e13.kasper.2023))

DimHeatmap(e13.kasper.2023, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(e13.kasper.2023)

e13.kasper.2023 <- FindNeighbors(e13.kasper.2023, dims = 1:30)
# e13.kasper.2023 <- FindClusters(e13.kasper.2023, resolution = 0.75)

# e13.kasper.2023 |> write_rds(file = "seurat_object_end_260430.rds")
e13.kasper.2023 <- readRDS("seurat_object_end_260430.rds")
# ^ Read this object back in but check dim reductions vs orig paper - should LECs and BECs overlap?


e13.kasper.2023 <- RunUMAP(e13.kasper.2023, dims = 1:30)

DimPlot(e13.kasper.2023, reduction = "umap", group.by = "subclustering_grouped", label = TRUE)

DimPlot(e13.kasper.2023, reduction = "umap", split.by = "embryonic_age", group.by = "embryonic_age")

#Confirm that we are now correctly using data and not counts for subsequent analysis
#Adapt the below code to make a Dim Plot suitable for Figure 3


#Remaking UMAP for Figure 3

#e13.kasper.2023[[1]] %>% DimPlot()
e13.kasper.2023 %>% DimPlot()

colour_assignments <-
  c(
    "FIB Deep1-3" = "#F4D166",
    "FIB Upper1-4" = "#BF4723",
    "MUSCLE Early" = "#B3E0A6",
    "FIB Muscle1-2" = "#E69F00"  ,
    "MUSCLE Late" = "#24693D",
    "FIB Lower" = "#EC6E1C", 
    "FIB Inter1-3" = "#F8AF50",
    "VESSEL BECs" = "#EBC4E1",
    "CHOND" = "#86D0B9", 
    "IMMU Macrophages" = "#26456E", 
    "NC SchwannCells" = "#D21E1C",
    "VESSEL MuralCells" = "#E14BA8",
    "IMMU MastCells" = "#A9D2DC",
    "MUSCLE Mid" = "#60A855",
    "EPI LatePlacode" = "#BABDBC",
    "IMMU DendriticCells" = "#4993C0",
    "EPI Basal1-4" = "#B254A5" ,
    "FIB Origin1-6" = "#9E3A26",
    "NC Melanocytes" = "#F78D80",
    "EPI EarlyPlacode" = "#D5D5D5",
    "VESSEL LECs" = "#835581",
    "EPI Periderm" = "#8A9497",
    "EPI Diff" = "#222222", 
    "EPI BasalTagln" = "#59636E"
  )

e13.kasper.2023 %>% 
  DimPlot(group.by = "subclustering_grouped", label = TRUE,
          pt.size = 1.1) +
  scale_colour_manual(values = colour_assignments)
#ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_UMAP_Dim30_Rough.svg"), device = svglite::svglite, fix_text_size = FALSE)

e13.kasper.2023@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = e13.kasper.2023@reductions$umap@cell.embeddings %>% as_tibble(rownames = "barcode"),
            by = "barcode") %>%
  ggplot(aes(x = umap_1, y = umap_2, colour = subclustering_grouped)) +
  geom_point(size = 0.7) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_colour_manual(values = colour_assignments) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))
#ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_UMAP_Dim30_2.svg"), device = svglite::svglite, fix_text_size = FALSE)

# Reordering classes
e13.kasper.2023@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = e13.kasper.2023@reductions$umap@cell.embeddings %>% as_tibble(rownames = "barcode"),
            by = "barcode") %>%
  mutate(subclustering_grouped = factor(subclustering_grouped,levels = 
                                          c("CHOND", "EPI Basal1-4", "EPI EarlyPlacode", "EPI LatePlacode", "EPI Periderm", "EPI BasalTagln", "EPI Diff", 
                                            "FIB Deep1-3", "FIB Inter1-3", "FIB Muscle1-2", "FIB Lower", "FIB Upper1-4", "FIB Origin1-6",
                                            "IMMU MastCells", "IMMU DendriticCells", "IMMU Macrophages", "MUSCLE Early", "MUSCLE Mid",
                                            "MUSCLE Late", "NC Melanocytes", "NC SchwannCells", "VESSEL BECs", "VESSEL LECs", "VESSEL MuralCells"))) %>%
  ggplot(aes(x = umap_1, y = umap_2, colour = subclustering_grouped)) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_colour_manual(values = colour_assignments) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_UMAP_Dim30_2.svg"), device = svglite::svglite, fix_text_size = FALSE,
#       width = 6.5, height = 5)

##################################################################################

# Remaking Feature Plots

# renv::restore()

# library("tidyverse")
# library("Seurat")

# renv::install("SeuratObject")

library(readr)

# e13.kasper.2023 <- read_rds("seurat_object_end_260430.rds")

LayerData(e13.kasper.2023, layer = "counts")
LayerData(e13.kasper.2023, layer = "data")

DimPlot(e13.kasper.2023, reduction = "umap", group.by = "subclustering_grouped", label = TRUE)

#Igf1
FeaturePlot(e13.kasper.2023, features = "Igf1", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf1.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf1.png"), bg = "white", scale=1.5)

#Igf1r
FeaturePlot(e13.kasper.2023, features = "Igf1r", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf1r.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf1r.png"), bg = "white", scale=1.5)

#Igf2
FeaturePlot(e13.kasper.2023, features = "Igf2", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf2.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf2.png"), bg = "white", scale=1.5)

#Igf2r
FeaturePlot(e13.kasper.2023, features = "Igf2r", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf2r.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igf2r.png"), bg = "white", scale=1.5)

#Igfbp
FeaturePlot(e13.kasper.2023, features = c("Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6"), cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbps.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbps.png"), bg = "white", scale=1.5)

#Igfbp1
FeaturePlot(e13.kasper.2023, features = "Igfbp1", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp1.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp1.png"), bg = "white", scale=1.5)

#Igfbp2
FeaturePlot(e13.kasper.2023, features = "Igfbp2", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp2.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp2.png"), bg = "white", scale=1.5)

#Igfbp3
FeaturePlot(e13.kasper.2023, features = "Igfbp3", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp3.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp3.png"), bg = "white", scale=1.5)

#Igfbp4
FeaturePlot(e13.kasper.2023, features = "Igfbp4", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp4.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp4.png"), bg = "white", scale=1.5)

#Igfbp5
FeaturePlot(e13.kasper.2023, features = "Igfbp5", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp5.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp5.png"), bg = "white", scale=1.5)

#Igfbp6
FeaturePlot(e13.kasper.2023, features = "Igfbp6", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp6.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp6.png"), bg = "white", scale=1.5)

#Igfbp7
FeaturePlot(e13.kasper.2023, features = "Igfbp7", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp7.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Igfbp7.png"), bg = "white", scale=1.5)


# Flt1
FeaturePlot(e13.kasper.2023, features = "Flt1", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Flt1.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Flt1.png"), bg = "white", scale=1.5)

# VEGFR2 (Kdr)
FeaturePlot(e13.kasper.2023, features = "Kdr", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Kdr.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Kdr.png"), bg = "white", scale=1.5)

# VEGFR3 (Flt4)
FeaturePlot(e13.kasper.2023, features = "Flt4", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Flt4.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Flt4.png"), bg = "white", scale=1.5)

# VEGF-A
FeaturePlot(e13.kasper.2023, features = "Vegfa", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfa.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfa.png"), bg = "white", scale=1.5)

# VEGF-B
FeaturePlot(e13.kasper.2023, features = "Vegfb", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfb.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfb.png"), bg = "white", scale=1.5)

# VEGF-C
FeaturePlot(e13.kasper.2023, features = "Vegfc", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfc.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Vegfc.png"), bg = "white", scale=1.5)

# Nrp1
FeaturePlot(e13.kasper.2023, features = "Nrp1", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Nrp1.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Nrp1.png"), bg = "white", scale=1.5)

# Nrp2
FeaturePlot(e13.kasper.2023, features = "Nrp2", cols = c("gray88", "navy"), order = TRUE)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Nrp2.svg"), bg = "white", scale=1.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "FP1_E13_Nrp2.png"), bg = "white", scale=1.5)




##################################################################################

#Cellchat Analysis of Jacob et al 2023 Data: Creating Cellchat Object and Checking for Correct Data

# renv::restore()

# library("tidyverse")
# library("Seurat")

# renv::install("SeuratObject")

# library(readr)

# e13.kasper.2023 <- read_rds("seurat_object_end_260430.rds")

LayerData(e13.kasper.2023, layer = "counts")
LayerData(e13.kasper.2023, layer = "data")

# library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# e13.kasper.cellchat <- createCellChat(e13.kasper.2023, group.by = "subclustering_grouped")
# 
# e13.kasper.cellchat <- tryCatch(
#   createCellChat(e13.kasper.2023, group.by = "subclustering_grouped"),
#   error = function(e) e
# )
# 
# class(e13.kasper.cellchat)
# 
# "subclustering_grouped" %in% colnames(e13.kasper.2023@meta.data)
# 
# LayerData(e13.kasper.2023, layer = "counts")
# LayerData(e13.kasper.2023, layer = "data")
# 
# class(e13.kasper.cellchat)
# 
# e13.kasper.2023[["RNA"]]@data <- as.matrix(
#   LayerData(e13.kasper.2023, layer = "data")
# )
# 
# data.input <- GetAssayData(e13.kasper.2023, slot = "data")
# 
# e13.kasper.cellchat <- createCellChat(
#   object = data.input,
#   meta = e13.kasper.2023@meta.data,
#   group.by = "subclustering_grouped"
# )

#Above code caused an error, starting over to run it cleanly

# library(Seurat)

data.input <- LayerData(e13.kasper.2023, layer = "data")
class(data.input)
dim(data.input)

meta <- e13.kasper.2023@meta.data

e13.kasper.cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "subclustering_grouped"
)

#Checking that the cellchat object is correctly using the normalised data

dim(e13.kasper.cellchat@data)
e13.kasper.cellchat@data[1:10, 1:10]

seurat.norm <- LayerData(e13.kasper.2023, layer = "data")

all.equal(
  dim(e13.kasper.cellchat@data),
  dim(seurat.norm)
)

table(e13.kasper.cellchat@idents)
head(e13.kasper.cellchat@meta)
levels(e13.kasper.cellchat@idents)

summary(e13.kasper.cellchat@data)
range(e13.kasper.cellchat@data)
table(e13.kasper.cellchat@idents)

rm(data.input)
rm(seurat.norm)
rm(meta)
gc()

# Continuing Jacob 2023 Cellchat Analysis: Running Cellchat Pipeline

CellChatDB <- CellChatDB.mouse

CellChatDB.use <- CellChatDB
e13.kasper.cellchat@DB <- CellChatDB.use

e13.kasper.cellchat <- subsetData(e13.kasper.cellchat)
# future::plan("multisession", workers = 4)

# renv::status()
# renv::restore()
# 
# options(future.globals.maxSize = 2 * 1024^3)

future::plan()
future::plan("sequential")

e13.kasper.cellchat <- identifyOverExpressedGenes(e13.kasper.cellchat)
e13.kasper.cellchat <- identifyOverExpressedInteractions(e13.kasper.cellchat)

# e13.kasper.cellchat <- computeCommunProb(e13.kasper.cellchat) #Error, seemed to encounter the same error in 2023

# e13.kasper.cellchat@idents <- droplevels(e13.kasper.cellchat@idents,exclude = setdiff(levels(e13.kasper.cellchat@idents),unique(e13.kasper.cellchat@idents))) #Didnt work

unique(e13.kasper.cellchat@idents)
table(e13.kasper.cellchat@idents)

e13.kasper.cellchat@idents <- droplevels(e13.kasper.cellchat@idents)
e13.kasper.cellchat@meta$subclustering_grouped <-
  droplevels(e13.kasper.cellchat@meta$subclustering_grouped)

table(e13.kasper.cellchat@idents)
levels(e13.kasper.cellchat@idents)

e13.kasper.cellchat <- computeCommunProb(e13.kasper.cellchat)

# e13.kasper.cellchat |> write_rds(file = "e13_kasper_cellchat.rds")
# e13.kasper.cellchat <- read_rds("e13_kasper_cellchat.rds")

e13.kasper.cellchat <- filterCommunication(e13.kasper.cellchat, min.cells = 10)
 
#df.net <- subsetCommunication(e13.kasper.cellchat) #Did not work and caused an error, decided to just skip this step
#df.net

e13.kasper.cellchat <- computeCommunProbPathway(e13.kasper.cellchat)

# length(e13.kasper.cellchat@LR$LRsig)
# dim(e13.kasper.cellchat@data)
# table(e13.kasper.cellchat@idents)
# 
# is.null(e13.kasper.cellchat@net$prob)
# dim(e13.kasper.cellchat@net$prob)
# 
# length(e13.kasper.cellchat@LR$LRsig)
# dim(e13.kasper.cellchat@data)
# summary(e13.kasper.cellchat@data)


e13.kasper.cellchat <- aggregateNet(e13.kasper.cellchat)

groupSize <- as.numeric(table(e13.kasper.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(e13.kasper.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(e13.kasper.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# The following section is just trying to save the image generated by the above code

# library(svglite)
# 
# svglite::svglite(file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasapar_1_CirclePlot_Count.svg"), width = 6, height = 6)
# 
# groupSize <- as.numeric(table(e13.kasper.cellchat@idents))
# netVisual_circle(e13.kasper.cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
# 
# dev.off()
# 
# svglite::svglite(file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasapar_1_CirclePlot_Weight.svg"), width = 6, height = 6)
# 
# netVisual_circle(e13.kasper.cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
# 
# dev.off()


mat <- e13.kasper.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("IGF")
vertex.receiver = c(22, 23)
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

# library(svglite)
# 
# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasapar_3_IGF_HierarchyPlot.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# pathways.show <- c("IGF")
# vertex.receiver <- c(22, 23)
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
# 
# dev.off()

levels(e13.kasper.cellchat@idents)

 
par(mfrow=c(1,1))
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_4_CirclePlot_Top100Percent.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# par(mfrow=c(1,1))
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle")
# 
# dev.off()
       
par(mfrow=c(1,1))
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.5)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_4_CirclePlot_Top50Percent.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# par(mfrow=c(1,1))
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.5)
# 
# dev.off()

par(mfrow=c(1,1))
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.2)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_4_CirclePlot_Top20Percent.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# par(mfrow=c(1,1))
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.2)
# 
# dev.off()

par(mfrow=c(1,1))
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.1)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_4_CirclePlot_Top10Percent.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# par(mfrow=c(1,1))
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "circle", top = 0.1)
# 
# dev.off()


par(mfrow=c(1,1))
netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "chord")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_5_ChordPlot.svg"), width = 8, height = 8, fix_text_size = FALSE)
#  
# par(mfrow=c(1,1))
# netVisual_aggregate(e13.kasper.cellchat, signaling = pathways.show, layout = "chord")
# 
# dev.off()


par(mfrow=c(1,1))
netVisual_heatmap(e13.kasper.cellchat, signaling = pathways.show, color.heatmap = "Reds")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_6_Heatmap.svg"), width = 8, height = 8, fix_text_size = FALSE)
# 
# par(mfrow=c(1,1))
# netVisual_heatmap(e13.kasper.cellchat, signaling = pathways.show, color.heatmap = "Reds")
# 
# dev.off()


netAnalysis_contribution(e13.kasper.cellchat, signaling = pathways.show) 

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_7_LRContribution.svg"), fix_text_size = FALSE)
# 
# netAnalysis_contribution(e13.kasper.cellchat, signaling = pathways.show)
# 
# dev.off()


pairLR.IGF <- extractEnrichedLR(e13.kasper.cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IGF[1,]
vertex.receiver = seq(22,23) # a numeric vector
netVisual_individual(e13.kasper.cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = "hierarchy")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_8_Igf1_Igf1rHierarchyPlot.svg"), fix_text_size = FALSE)
# 
# pairLR.IGF <- extractEnrichedLR(e13.kasper.cellchat, signaling = pathways.show, geneLR.return = FALSE)
# LR.show <- pairLR.IGF[1,]
# vertex.receiver = seq(22,23) # a numeric vector
# netVisual_individual(e13.kasper.cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = "hierarchy")
# 
# dev.off()


netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_9_Igf1_Igf1rCirclePlot_Top100.svg"), fix_text_size = FALSE)
# 
# netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# 
# dev.off()

netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", top = 0.1)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_9_Igf1_Igf1rCirclePlot_Top10.svg"), fix_text_size = FALSE)
# 
# netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", top = 0.1)
# 
# dev.off()


netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_10_Igf1_Igf1rChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_individual(e13.kasper.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# 
# dev.off()


netVisual_bubble(e13.kasper.cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

netVisual_bubble(e13.kasper.cellchat, sources.use = 22, targets.use = c(1:24), signaling = "IGF", remove.isolate = FALSE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_11_BEC_Igf_BubblePlot.svg"), fix_text_size = FALSE)
# 
# netVisual_bubble(e13.kasper.cellchat, sources.use = 22, targets.use = c(1:24), signaling = "IGF", remove.isolate = FALSE)
# 
# dev.off()

netVisual_bubble(e13.kasper.cellchat, sources.use = 23, targets.use = c(1:24), signaling = "IGF", remove.isolate = FALSE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_12_LEC_Igf_BubblePlot.svg"), fix_text_size = FALSE)
# 
# netVisual_bubble(e13.kasper.cellchat, sources.use = 23, targets.use = c(1:24), signaling = "IGF", remove.isolate = FALSE)
# 
# dev.off()

netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(2), signaling = "IGF", remove.isolate = FALSE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_13_EpidermisTarget_IgfBubblePlot.svg"), fix_text_size = FALSE)
# 
# netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(2), signaling = "IGF", remove.isolate = FALSE)
# 
# dev.off()

netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(22), signaling = "IGF", remove.isolate = FALSE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_14_BECTarget_IgfBubblePlot.svg"), fix_text_size = FALSE)
# 
# netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(22), signaling = "IGF", remove.isolate = FALSE)
# 
# dev.off()

netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(23), signaling = "IGF", remove.isolate = FALSE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_15_LECTarget_IgfBubblePlot.svg"), fix_text_size = FALSE)
# 
# netVisual_bubble(e13.kasper.cellchat, sources.use = 1:24, targets.use = c(23), signaling = "IGF", remove.isolate = FALSE)
# 
# dev.off()


pairLR.use <- extractEnrichedLR(e13.kasper.cellchat, signaling = "IGF")
netVisual_bubble(e13.kasper.cellchat, sources.use = c(22, 23), targets.use = c(1:24), pairLR.use = pairLR.use, remove.isolate = TRUE)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_16_BEC+LEC_IgfBubblePlot.svg"), fix_text_size = FALSE)
# 
# pairLR.use <- extractEnrichedLR(e13.kasper.cellchat, signaling = "IGF")
# netVisual_bubble(e13.kasper.cellchat, sources.use = c(22, 23), targets.use = c(1:24), pairLR.use = pairLR.use, remove.isolate = TRUE)
# 
# dev.off()


netVisual_chord_gene(e13.kasper.cellchat, sources.use = 22,23, targets.use = c(2,22,23), lab.cex = 0.5,legend.pos.y = 30)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_17_ChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = 22,23, targets.use = c(2,22,23), lab.cex = 0.5,legend.pos.y = 30)
# 
# dev.off()

netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 2, legend.pos.x = 15)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_18_ChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 2, legend.pos.x = 15)
# 
# dev.off()

netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 22, legend.pos.x = 15)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_19_ChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 22, legend.pos.x = 15)
# 
# dev.off()

netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 23, legend.pos.x = 15)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_20_ChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 23, legend.pos.x = 15)
# 
# dev.off()


netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = c(2, 4, 8, 9, 10, 12, 13, 15, 22, 23), signaling = "IGF",legend.pos.x = 8)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_21_VESSELBEC+LEC.IgfChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = c(2, 4, 8, 9, 10, 12, 13, 15, 22, 23), signaling = "IGF",legend.pos.x = 8)
# 
# dev.off()

netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(2, 4, 8, 9, 10, 12, 13, 15, 22, 23), targets.use = 2, signaling = "IGF",legend.pos.x = 8)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_22_Epidermis.IgfChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(2, 4, 8, 9, 10, 12, 13, 15, 22, 23), targets.use = 2, signaling = "IGF",legend.pos.x = 8)
# 
# dev.off()

netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 2, slot.name = "netP", legend.pos.x = 10)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_23_AngiocrineChordPlot.svg"), fix_text_size = FALSE)
# 
# netVisual_chord_gene(e13.kasper.cellchat, sources.use = c(22,23), targets.use = 2, slot.name = "netP", legend.pos.x = 10)
# 
# dev.off()


plotGeneExpression(e13.kasper.cellchat, signaling = "IGF")
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_24_IGF_ViolinPlot.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_24_IGF_ViolinPlot.png"))
       

plotGeneExpression(e13.kasper.cellchat, signaling = "VEGF")
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_25_VEGF_ViolinPlot.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_25_VEGF_ViolinPlot.png"))


grep(pattern = "^Igfbp", x = rownames(e13.kasper.cellchat@data), value = TRUE)

plotGeneExpression(e13.kasper.cellchat, features = c("Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_26_Igfbp_ViolinPlot.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_26_Igfbp_ViolinPlot.png"))


grep(pattern = "^Vegf", x = rownames(e13.kasper.cellchat@data), value = TRUE)

plotGeneExpression(e13.kasper.cellchat, features = c("Vegfa", "Vegfb", "Vegfc"))
ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_99_Vegf_ViolinPlot.svg"))
ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_99_Vegf_ViolinPlot.png"))


plotGeneExpression(e13.kasper.cellchat, signaling = "EGF")
ggsave(filename = file.path(figure_out_dir, "E13.5_Kasper_99_EGF_ViolinPlot.svg"))
ggsave(filename = file.path(figure_out_dir, "E13.5_Kasper_99_EGF_ViolinPlot.png"))


plotGeneExpression(e13.kasper.cellchat, signaling = "EPHA")
ggsave(filename = file.path(figure_out_dir, "E13.5_Kasper_99_EPHA_ViolinPlot.svg"))
ggsave(filename = file.path(figure_out_dir, "E13.5_Kasper_99_EPHA_ViolinPlot.png"))


plotGeneExpression(e13.kasper.cellchat, signaling = "IGF", enriched.only = FALSE)



e13.kasper.cellchat.2 <- netAnalysis_computeCentrality(e13.kasper.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways


netAnalysis_signalingRole_network(e13.kasper.cellchat.2, signaling = pathways.show, width = 16, height = 5, font.size = 9)

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_27_NetworkCentrality.svg"), fix_text_size = FALSE)
# 
# netAnalysis_signalingRole_network(e13.kasper.cellchat.2, signaling = pathways.show, width = 16, height = 5, font.size = 9)
# 
# dev.off()


gg1 <- netAnalysis_signalingRole_scatter(e13.kasper.cellchat.2)
gg1 <- gg1 + labs(title = "All Signalling")
gg2 <- netAnalysis_signalingRole_scatter(e13.kasper.cellchat.2, signaling = "IGF")
gg2 <- gg2 + labs(title = "Igf Signalling")
gg1 + gg2
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_28_Visualise Senders.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_28_Visualise Senders.png"))


ht1 <- netAnalysis_signalingRole_heatmap(e13.kasper.cellchat.2, pattern = "outgoing", height = 20, font.size = 6)
ht2 <- netAnalysis_signalingRole_heatmap(e13.kasper.cellchat.2, pattern = "incoming", height = 20, font.size = 6)
ht1 + ht2

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_29_SignallingPatternHeatmap.svg"), fix_text_size = FALSE)
# ht1 + ht2
# dev.off()


ht <- netAnalysis_signalingRole_heatmap(e13.kasper.cellchat.2, pattern = "outgoing", signaling = "IGF", height = 1)
hti <- netAnalysis_signalingRole_heatmap(e13.kasper.cellchat.2, pattern = "incoming", signaling = "IGF", height = 1)
ht + hti

# svglite::svglite(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_30_SubsetPatternHeatmap.svg"), fix_text_size = FALSE)
# ht + hti
# dev.off()

 
library(NMF)
renv::install("ggalluvial")
library(ggalluvial)

packageVersion("ggplot2")
packageVersion("ggalluvial")

renv::install("assertthat")
# install.packages("pkgmaker")
renv::install("https://cran.r-project.org/src/contrib/Archive/pkgmaker/pkgmaker_0.32.10.tar.gz")
library(pkgmaker)
library(registry)
library(rngtools)
library(cluster)

selectK(e13.kasper.cellchat.2, pattern = "outgoing")
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_31_Visualise Senders.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_31_Visualise Senders.png"))


nPatterns = 3
e13.kasper.cellchat.2 <- identifyCommunicationPatterns(e13.kasper.cellchat.2, pattern = "outgoing", k = nPatterns, font.size = 4, height = 10, width = 10)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_32_Outgoing.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_32_Outgoing.png"))


netAnalysis_river(e13.kasper.cellchat.2, pattern = "outgoing",  font.size = 2.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_33_Outgoing.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_33_Outgoing.png"))


netAnalysis_dot(e13.kasper.cellchat.2, pattern = "outgoing", font.size = 7)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_34_Outgoing.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_34_Outgoing.png"))


selectK(e13.kasper.cellchat.2, pattern = "incoming")
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_35_Incoming.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_35_Incoming.png"))


nPatterns = 4
e13.kasper.cellchat.2 <- identifyCommunicationPatterns(e13.kasper.cellchat.2, pattern = "incoming", k = nPatterns)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_36_Incoming.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_36_Incoming.png"))


netAnalysis_river(e13.kasper.cellchat.2, pattern = "incoming", font.size = 2.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_37_Incoming.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_37_Incoming.png"))

netAnalysis_dot(e13.kasper.cellchat.2, pattern = "incoming", font.size = 7)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_38_Incoming.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_38_Incoming.png"))


e13.kasper.cellchat.2 <- computeNetSimilarity(e13.kasper.cellchat.2, type = "functional")

reticulate::py_install(packages = 'umap-learn')

e13.kasper.cellchat.2 <- netEmbedding(e13.kasper.cellchat.2, type = "functional")

e13.kasper.cellchat.2 <- netClustering(e13.kasper.cellchat.2, type = "functional", do.parallel = FALSE)
netVisual_embedding(e13.kasper.cellchat.2, type = "functional", label.size = 3.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_39_FunctionalSimilarity.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_39_FunctionalSimilarity.png"))


e13.kasper.cellchat.2 <- computeNetSimilarity(e13.kasper.cellchat.2, type = "structural")
e13.kasper.cellchat.2 <- netEmbedding(e13.kasper.cellchat.2, type = "structural")
e13.kasper.cellchat.2 <- netClustering(e13.kasper.cellchat.2, type = "structural", do.parallel = FALSE)
netVisual_embedding(e13.kasper.cellchat.2, type = "structural", label.size = 3.5)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_40_StructuralSimilarity.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_40_StructuralSimilarity.png"))


netVisual_embeddingZoomIn(e13.kasper.cellchat.2, type = "structural", nCol = 2)
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_41_StructuralSimilarity2.svg"))
# ggsave(filename = file.path("/media/danielbishop/T7/scRNA/IGF1_Project_Redux", "E13.5_Kasper_41_StructuralSimilarity2.png"))
