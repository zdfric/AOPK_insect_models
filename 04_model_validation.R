# =====================================================================
# Title: 04_model_validation.R
# Project: Potential distribution models for Habitats Directive Insects (AOPK)
# Description: Cross-validation of predicted distribution models against
#              recent field data (latest survey seasons) and spatial
#              correlation analysis.
# =====================================================================

library(raster)
library(sf)
library(ENMTools)
library(terra)
library(spatialEco)

# Load saved predictions and base environmental data
load("01_data_preprocessed.RData")
romportl.bart <- raster("arion_potential_distribution_BART_100m_WGS84.tif")
arion.model <- raster("arion_potential_distribution_MaxEnt_coarse.tif")

# ---------------------------------------------------------------------
# 1. Model Suitability Thresholding & Habitat Overlay
# ---------------------------------------------------------------------
# Define threshold based on MaxEnt training sensitivity/specificity
bart_threshold <- BART.threshold

# habitats used in 95% of occassions


# Filter target key habitat categories (KVES)
kves_all <- st_read("/mnt/data/ESRI/AOPK/KVES2022/kves_all.gpkg")
kves.95<-kves_all[(kves_all$KATEGORIE %in% freq.95$Var1),]

kves.95 <- st_transform(st_as_sf(kves.95), 4326)

# Extract predicted suitability to habitat geometries
kves.95$ver.arion.romp <- raster::extract(romportl.bart, as_Spatial(kves.95), fun = mean, na.rm = TRUE)

# Retain polygons exceeding model presence threshold
kves.95_selected <- kves.95[!is.na(kves.95$ver.arion.romp) & (kves.95$ver.arion.romp > bart_threshold), ]
st_write(kves.95_selected, "arion_modelled_KVES_thresholded.gpkg", delete_dsn = TRUE)

# Load distribution data, filter by distance buffer within recent records (10 km buffer)
arion_lok.krovak <- st_as_sf(read.csv2("arion_ndop.csv"), coords = c("X", "Y"), crs = 5514)
buffer10 <- st_buffer(arion_lok.krovak, dist = 10000)
buffer10 <- st_union(buffer10)
buffer10.4326 <- st_transform(st_make_valid(buffer10), 4326)

# select habitat polygons for the records within the 10 km buffer around the known occurrence
kves.95.10km <- kves.95_selected[st_intersects(st_make_valid(kves.95_selected), buffer10.4326, sparse = FALSE), ]
st_write(kves.95.10km, "arion_modelled_KVES_10km_buffer.gpkg", delete_dsn = TRUE)

# ---------------------------------------------------------------------
# 2. Validation with Recent Field Surveys (Polygons) not included in model training
# ---------------------------------------------------------------------
arion_polygons <- st_read("arion_polygons_new2.gpkg") # New survey records
blank_raster_rompwgs84 <- romportl.bart
values(blank_raster_rompwgs84) <- NA

# Rasterize new occurrences onto model grid
arion_polygons_ras <- raster::rasterize(arion_polygons, blank_raster_rompwgs84, 
                                        field = arion_polygons$nalezen, na.rm = TRUE, fun = "mean")

writeRaster(arion_polygons_ras,"new_occurrences.gpkg",owerwrite=T)

# Local Spearman rank correlation between predictions and new observations
arion.corlocal <- corLocal(romportl.bart, arion_polygons_ras, method = "spearman", use = "complete.obs")

pdf("model_validation_correlation.pdf")
plot(arion.corlocal, main = "Local Correlation: Model Prediction vs New Seasons Field Data")
dev.off()

writeRaster(arion.corlocal, "model_validation_correlation.gpkg",owerwrite=T)

print("Model validation and spatial cross-evaluation finished.")