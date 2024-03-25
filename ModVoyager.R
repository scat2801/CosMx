if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(Voyager)
library(SFEData)
library(SingleCellExperiment)
library(SpatialExperiment)
BiocManager::install("DelayedArray")
library(scater) # devel version of plotExpression
library(scran)
library(bluster)
library(ggplot2)
library(patchwork)
library(stringr)
library(spdep)
library(BiocParallel)
library(BiocSingular)
theme_set(theme_bw())

BiocManager::install("SFEData")

require("SpatialFeatureExperiment")
require("SFEData")
library(Voyager)

(sfe <- HeNSCLCData())

plotGeometry(sfe, MARGIN = 2L, type = "cellSeg")
plotCellBin2D(sfe, hex = TRUE)

#QC
names(colData(sfe))

# Function to plot violin plot for distribution and spatial at once
plot_violin_spatial <- function(sfe, feature) {
  violin <- plotColData(sfe, feature, point_fun = function(...) list())
  spatial <- plotSpatialFeature(sfe, feature, colGeometryName = "centroids",
                                scattermore = TRUE)
  violin + spatial +
    plot_layout(widths = c(1, 2))
}

plot_violin_spatial(sfe, "nCounts")

summary(sfe$nCounts)

n_panel <- 960
colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel
plotColDataHistogram(sfe, c("nCounts_normed", "nGenes_normed"))

plot_violin_spatial(sfe, "nGenes")

summary(sfe$nGenes)
plotColData(sfe, x = "nCounts", y = "nGenes", bins = 100)
colData(sfe)$is_empty <- colData(sfe)$nCounts < 1
plotSpatialFeature(sfe, "is_empty", "cellSeg")

plotColData(sfe, x = "Area", y = "is_empty")
plot_violin_spatial(sfe, "Area")
plotColData(sfe, x = "nCounts", y = "Area", bins = 100) + theme_bw()

neg_inds <- str_detect(rownames(sfe), "^NegPrb")
sum(neg_inds)
colData(sfe)$prop_neg <- colSums(counts(sfe)[neg_inds,])/colData(sfe)$nCounts
plot_violin_spatial(sfe, "prop_neg")

plotColData(sfe, x = "nCounts",y = "prop_neg", bins = 100)

plotColDataHistogram(sfe, "prop_neg") +
  scale_x_log10()

# Remove low quality cells
(sfe <- sfe[,!sfe$is_empty & sfe$prop_neg < 0.1])


#Markers
plotSpatialFeature(sfe, c("AspectRatio", "Mean.DAPI", "Mean.MembraneStain", 
                          "Mean.PanCK", "Mean.CD45", "Mean.CD3"),
                   colGeometryName = "centroids", ncol = 2, scattermore = TRUE)

#Genes
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$is_neg <- neg_inds

plotRowData(sfe, x = "means", y = "vars", bins = 50) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  coord_equal()

as.data.frame(rowData(sfe)[neg_inds,]) |> 
  ggplot(aes(means, vars)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  coord_equal()

plotRowData(sfe, x = "means", y = "is_neg") +
  scale_y_log10() +
  annotation_logticks(sides = "b")

#Spatial autocorrelation in QC metrics

system.time(
  colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh",
                                                dist_type = "idw", k = 5, 
                                                style = "W")
)

features_use <- c("nCounts", "nGenes", "Area", "AspectRatio")
sfe <- colDataMoransI(sfe, features_use, colGraphName = "knn5")

colFeatureData(sfe)[features_use,]

system.time(
  sfe <- colDataUnivariate(sfe, "sp.correlogram", features = features_use,
                           colGraphName = "knn5", order = 6, zero.policy = TRUE,
                           BPPARAM = SnowParam(2))
)

plotCorrelogram(sfe, features_use)
sfe <- colDataUnivariate(sfe, "moran.plot", "nCounts", colGraphName = "knn5")
p1 <- moranPlot(sfe, "nCounts", binned = TRUE, plot_influential = FALSE)
p2 <- moranPlot(sfe, "nCounts", binned = TRUE)
p1 / p2 + plot_layout(guides = "collect")

sfe <- colDataUnivariate(sfe, "localmoran", "nCounts", colGraphName = "knn5")
plotLocalResult(sfe, "localmoran", "nCounts", colGeometryName = "cellSeg",
                divergent = TRUE, diverge_center = 0)

#Data normalisation
sfe <- logNormCounts(sfe)

#Moran’s I
sfe <- runMoransI(sfe, features = rownames(sfe), 
                  BPPARAM = SnowParam(2, progressbar = TRUE))

plotRowData(sfe, x = "moran_sample01", y = "is_neg") +
  geom_hline(yintercept = 0, linetype = 2)

top_moran <- rownames(sfe)[order(rowData(sfe)$moran_sample01, decreasing = TRUE)[1:6]]
plotSpatialFeature(sfe, top_moran, colGeometryName = "centroids", 
                   scattermore = TRUE, ncol = 2)

#Non-spatial dimension reduction and clustering
set.seed(29)
sfe <- runPCA(sfe, ncomponents = 30, scale = TRUE, BSPARAM = IrlbaParam())
ElbowPlot(sfe, ndims = 30)
plotDimLoadings(sfe, dims = 1:6)

spatialReducedDim(sfe, "PCA", 6, colGeometryName = "centroids", divergent = TRUE,
                  diverge_center = 0, ncol = 2, scattermore = TRUE)

colData(sfe)$cluster <- clusterRows(reducedDim(sfe, "PCA")[,1:15],
                                    BLUSPARAM = SNNGraphParam(
                                      cluster.fun = "leiden",
                                      cluster.args = list(
                                        resolution_parameter = 0.5,
                                        objective_function = "modularity")))


data("ditto_colors")

plotPCA(sfe, ncomponents = 4, colour_by = "cluster") +
  scale_color_manual(values = ditto_colors)

plotSpatialFeature(sfe, "cluster", colGeometryName = "cellSeg")

#Differential expression
markers <- findMarkers(sfe, groups = colData(sfe)$cluster,
                       test.type = "wilcox", pval.type = "all", direction = "up")

genes_use <- vapply(markers, function(x) rownames(x)[1], FUN.VALUE = character(1))
plotExpression(sfe, genes_use, x = "cluster", point_fun = function(...) list())

genes_use2 <- unique(unlist(lapply(markers, function(x) rownames(x)[1:5])))
plotGroupedHeatmap(sfe, genes_use2, group = "cluster", colour = scales::viridis_pal()(100))

#Local spatial statistics of marker genes
plotSpatialFeature(sfe, genes_use, colGeometryName = "centroids", ncol = 2,
                   scattermore = TRUE)

rowData(sfe)[genes_use, "moran_sample01", drop = FALSE]
sfe <- runUnivariate(sfe, "localmoran", features = genes_use, colGraphName = "knn5",
                     BPPARAM = SnowParam(2))

plotLocalResult(sfe, "localmoran", features = genes_use, 
                colGeometryName = "centroids", ncol = 2, divergent = TRUE,
                diverge_center = 0, scattermore = TRUE)

sfe <- runUnivariate(sfe, "LOSH", features = genes_use, colGraphName = "knn5",
                     BPPARAM = SnowParam(2))

plotLocalResult(sfe, "LOSH", features = genes_use, 
                colGeometryName = "centroids", ncol = 2, scattermore = TRUE)







