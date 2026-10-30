# =====================================================================
# Title: 03_prediction_mapping.R
# Project: Potential distribution models for Habitats Directive Insects (AOPK)
# Description: Load fitted model objects, restore BART memory state,
#              predict onto current or new predictor versions, and 
#              export GeoTIFFs and PDF summary plots.
# Phengaris arion (i.e., arion in the datasets) used as an example, replace
# by other species if you need
# =====================================================================

library(raster)
library(sf)
library(terra)
library(embarcadero)
library(ggplot2)
library(gridExtra)

# Load fit models
load("fitted_models.RData")

# Restore BART internal state pointers safely
if (!is.null(sdm.bart2_coarse_backup) && (is.null(sdm.bart2$fit$state) || length(sdm.bart2$fit$state) == 0)) {
  sdm.bart2 <- sdm.bart2_coarse_backup
}
invisible(sdm.bart2$fit$state)

if (!is.null(sdm.bart2.r_backup) && (is.null(sdm.bart2.r$fit$state) || length(sdm.bart2.r$fit$state) == 0)) {
  sdm.bart2.r <- sdm.bart2.r_backup
}
invisible(sdm.bart2.r$fit$state)

# Load environmental predictor layers
load("01_data_preprocessed.RData")

# =====================================================================
# 1. Coarse-Scale Predictions (~30 arc-sec)
# =====================================================================
layers <- stack(list.files(c("./rasters_kves", "./rasters_biotopes"), pattern = "tif", full.names = TRUE))
layers <- addLayer(layers, current)
layers <- addLayer(layers, dem)

# MaxEnt Prediction
arion.model <- mean(predict(data.mx, layers))
writeRaster(arion.model, "arion_potential_distribution_MaxEnt_coarse.tif", format = "GTiff", overwrite = TRUE)

# BART Prediction
sdm.bart.map <- predict(sdm.bart2, layers)
writeRaster(sdm.bart.map, "arion_potential_distribution_BART_coarse.tif", format = "GTiff", overwrite = TRUE)

# Export Coarse Distribution PDF
pdf("coarse_models_distribution.pdf")
plot(arion.model, main = "MaxEnt Coarse Model")
plot(sdm.bart.map, main = "BART Coarse Model")
dev.off()

# =====================================================================
# 2. Fine-Scale Predictions (100m Grid - Chunked Prediction)
# =====================================================================
romp <- stack(list.files("./rasters_romp", full.names = TRUE, pattern = "tif"))
climate_r <- brick("climate_100m.grd")
dem.r <- brick("dem_100m.grd")

layer.r <- stack(romp, romportl, climate_r, dem.r)
layer.r <- dropLayer(layer.r, "DEM_CZ_WGS84")

# Divide raster domain into 4 spatial tiles to manage memory during 100m BART execution
cropped1 <- c(xmin(layer.r), ((xmax(layer.r) - xmin(layer.r)) / 2) + xmin(layer.r),
              ymin(layer.r), ((ymax(layer.r) - ymin(layer.r)) / 2) + ymin(layer.r))
cropped2 <- c(xmin(layer.r), ((xmax(layer.r) - xmin(layer.r)) / 2) + xmin(layer.r),
              ((ymax(layer.r) - ymin(layer.r)) / 2) + ymin(layer.r) + 1, ymax(layer.r))
cropped3 <- c(((xmax(layer.r) - xmin(layer.r)) / 2) + xmin(layer.r) + 1, xmax(layer.r),
              ymin(layer.r), ((ymax(layer.r) - ymin(layer.r)) / 2) + ymin(layer.r))
cropped4 <- c(((xmax(layer.r) - xmin(layer.r)) / 2) + xmin(layer.r) + 1, xmax(layer.r),
              ((ymax(layer.r) - ymin(layer.r)) / 2) + ymin(layer.r) + 1, ymax(layer.r))

layer.r1 <- crop(layer.r, cropped1)
layer.r2 <- crop(layer.r, cropped2)
layer.r3 <- crop(layer.r, cropped3)
layer.r4 <- crop(layer.r, cropped4)

# Predict tile by tile
sdm.bart.map.romp1 <- predict(sdm.bart2.r, layer.r1); gc()
sdm.bart.map.romp2 <- predict(sdm.bart2.r, layer.r2); gc()
sdm.bart.map.romp3 <- predict(sdm.bart2.r, layer.r3); gc()
sdm.bart.map.romp4 <- predict(sdm.bart2.r, layer.r4); gc()

# Merge tiles back into a single surface
sdm.bart.map.romp <- merge(sdm.bart.map.romp1, sdm.bart.map.romp2, sdm.bart.map.romp3, sdm.bart.map.romp4)
crs(sdm.bart.map.romp) <- CRS("+init=epsg:3035")

writeRaster(sdm.bart.map.romp, "arion_potential_distribution_BART_100m_LAMBERT.tif", format = "GTiff", overwrite = TRUE)

# Reproject fine prediction back to WGS84
romportl.bart <- projectRaster(sdm.bart.map.romp, crs = 4326)
writeRaster(romportl.bart, "arion_potential_distribution_BART_100m_WGS84.tif", format = "GTiff", overwrite = TRUE)

print("Predictions and raster export completed successfully.")