#Visualisation of CosMx data; Adapted from Satijalab script
#Tested on Nanostring CosMx NSCLC data

library(Seurat)
library(tidyverse)

nano.obj <- LoadNanostring(data.dir = "CosMx/Lung5_Rep1-Flat_files_and_images/", fov ="lung5.rep1", assay = "Nanostring")

# add in precomputed Azimuth annotations
azimuth.data <- readRDS("CosMx/nanostring_data.Rds")
nano.obj <- AddMetaData(nano.obj, metadata = azimuth.data$annotations)
nano.obj[["proj.umap"]] <- azimuth.data$umap
Idents(nano.obj) <- nano.obj$predicted.annotation.l1

# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^2)
nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10), verbose = FALSE)

# text display of annotations and prediction scores
head(slot(object = nano.obj, name = "meta.data")[2:5])

#visualize the Nanostring cells and annotations, projected onto the reference-defined UMAP
DimPlot(nano.obj)

#Visualisation of cell type and expression localization patterns
ImageDimPlot(nano.obj, fov = "lung5.rep1", axes = TRUE, cols = "glasbey")

#Visualising fewer groups
ImageDimPlot(nano.obj, fov = "lung5.rep1", cells = WhichCells(nano.obj, idents = c("Basal", "Macrophage", "Smooth Muscle", "CD4 T")), cols = c("red", "green", "blue", "orange"), size = 0.6)
#Gene expression markers visualisation
VlnPlot(nano.obj, features = "KRT17", assay = "Nanostring", layer = "counts", pt.size = 0.1, y.max = 30) +
  NoLegend()

FeaturePlot(nano.obj, features = "KRT17", max.cutoff = "q95")

p1 <- ImageFeaturePlot(nano.obj, fov = "lung5.rep1", features = "KRT17", max.cutoff = "q95")
p2 <- ImageDimPlot(nano.obj, fov = "lung5.rep1", alpha = 0.3, molecules = "KRT17", nmols = 10000) +
  NoLegend()
p1 + p2

# Plot some of the molecules which seem to display spatial correlation with each other
ImageDimPlot(nano.obj, fov = "lung5.rep1", group.by = NA, alpha = 0.3, molecules = c("KRT17", "C1QA", "IL7R", "TAGLN"), nmols = 20000)

#Zoom in on basal-rich region
basal.crop <- Crop(nano.obj[["lung5.rep1"]], x = c(159500, 164000), y = c(8700, 10500))
nano.obj[["zoom1"]] <- basal.crop
DefaultBoundary(nano.obj[["zoom1"]]) <- "segmentation"

ImageDimPlot(nano.obj, fov = "zoom1", cols = "polychrome", coord.fixed = FALSE)

# note the clouds of TPSAB1 molecules denoting mast cells
ImageDimPlot(nano.obj, fov = "zoom1", cols = "polychrome", alpha = 0.3, molecules = c("KRT17", "IL7R",
                                                                                      "TPSAB1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)