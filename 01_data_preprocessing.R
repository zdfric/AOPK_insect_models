# =====================================================================
# Title: 01 - Data Preprocessing and Harmonization
# Project: Potential distribution models of insect species for AOPK
# Description: Loads nationwide vector layers (KVES, biotopes, occurrences)
#              using 'sf' and 'terra' packages, cleans datasets, 
#              handles Czech diacritics, and prepares spatial predictors.
# Phengaris arion (i.e., arion in the datasets) used as an example, replace
# by other species if you need
# =====================================================================

library(terra)
library(sf)
library(vegan)
library(parallel)

ncore <- parallel::detectCores()

# Load Boundaries & Grid Templates
chko_map <- st_read("CHKO+Soutok+KrusneHory.gpkg")
cz_map <- st_read("../../../ESRI/ARCCR500/JTSK/SHAPES/kraj2000.shp")
st_crs(cz_map) <- 5514
cz_map <- st_transform(cz_map, 4326)

# Load Species Occurrence Data (e.g., Phengaris arion)
arion_lok <- read.csv2("arion_ndop.csv")
arion_lok$X <- as.numeric(arion_lok$X)
arion_lok$Y <- as.numeric(arion_lok$Y)
arion_lok <- st_as_sf(arion_lok, coords = c("X", "Y"), crs = 5514)
arion_lok <- st_transform(arion_lok, 4326)

# Filter out bad sources
if ("ZDROJ" %in% names(arion_lok)) {
  arion_lok <- arion_lok[arion_lok$ZDROJ != "Entomologický ústav Biologického centra AV ČR (2012) Databáze SOM a Entomologický ústav BC AV CR - denní motýli.", ]
}

### Coarse scale (~30 arc-sec, WGS84) 
# Nationwide KVES Habitats Processing
regions <- c("Jihocesky", "Jihomoravsky", "Kralovehradecky", "Karlovarsky", "Liberecky", 
             "Olomoucky", "Moravskoslezsky", "Plzensky", "Stredocesky", "Ustecky", 
             "Vysocina", "Zlinsky", "Praha", "Pardubicky")

if (!file.exists("kves_cely.gpkg")) {
  kves_list <- lapply(regions, function(reg) {
    st_read(paste0("/../../AOPK/KVES2022/KVES_", reg, ".shp"), options = "ENCODING=CP1250")
  })
  kves_all <- do.call(rbind, kves_list)
  
  # Fix Czech diacritics encoding
  kves_all$KATEGORIE <- gsub("<c3><a1>", "á", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c4><8d>", "č", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c3><a9>", "é", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c4><9b>", "ě", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c3><ad>", "í", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><88>", "ň", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><99>", "ř", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><a1>", "š", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><a5>", "ť", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><af>", "ů", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c3><bd>", "ý", kves_all$KATEGORIE)
  kves_all$KATEGORIE <- gsub("<c5><be>", "ž", kves_all$KATEGORIE)
  
  kves_all <- st_transform(kves_all, 4326)
  kves_all$value <- 1
  st_write(kves_all, "kves_cely.gpkg", append = FALSE)
} else {
  kves_all <- st_read("kves_cely.gpkg")
}

## Climate and Template Raster Setup (terra based)
blank_raster <- rast("../../Past_climate/current_Chelsea_0.3/bio_1.tif")
blank_raster <- crop(blank_raster, vect(kves_all))
values(blank_raster) <- 0

# Rasterize individual KVES habitat categories
kves_all_codes <- as.data.frame(unique(kves_all$KATEGORIE))
dir.create("./rasters_kves", showWarnings = FALSE)for (i in 1:nrow(kves_all_codes)) {
hab <- raster::rasterize(kves_all[kves_all$KATEGORIE == kves_all_codes[i, ], ], 
blank_raster, field = kves_all$value, background = 0, update = TRUE)
names(hab) <- paste(kves_all_codes[i, ])
writeRaster(hab, paste0("./rasters_kves/", kves_all_codes[i, ], ".tif"), overwrite = TRUE, format = "GTiff")
}
unlink("./rasters_kves/*.xml")


# Detailed Biotope Mapping and Mosaic Decomposition
if (!file.exists("biotopes_all")) {
  biotopes_list <- lapply(regions, function(reg) {
    st_read(paste0("/../../AOPK/2_mapovani", reg, ".shp"), options = "ENCODING=CP1250")
  })
  biotopes_all <- do.call(rbind, biotopes_list)
  biotopes_all <- st_transform(biotopes_all, 4326)
  biotopes_all$value <- 1
  st_write(biotopes_all, "biotopes_all.gpkg", append = FALSE)
} else {
  biotopes_all <- st_read("biotopes_all.gpkg")
}

biotopes_all$value <- 1

# Extract biotopes present at target insect locations
sf::sf_use_s2(FALSE)
arion_lok_biotopes <- over(arion_lok, sf:::as_Spatial(biotopes_all))
biotopes_all_codes <- as.data.frame(unique(arion_lok_biotopes$BIOTOP_SEZ[!is.na(arion_lok_biotopes$BIOTOP_SEZ)]))

# Process biotope mosaics and proportions
mozaiky <- as.data.frame(cbind(
unlist(strsplit(biotopes_all_codes[, 1], ",")),
unlist(regmatches(biotopes_all_codes[, 1], gregexpr("(?<=\().+?(?=\))", biotopes_all_codes[, 1], perl = TRUE)))
))
names(mozaiky) <- c("biotop", "proporce")
mozaiky$nazev.biotopu <- gsub("^(.*?)\\(.*", "\\1", unlist(strsplit(mozaiky[, 1], ",")))
mozaiky <- mozaiky[mozaiky$proporce != 0, ]
mozaiky$biotop <- trimws(mozaiky$biotop)
mozaiky$nazev.biotopu <- trimws(mozaiky$nazev.biotopu)
mozaiky$proporce <- as.numeric(mozaiky$proporce)
mozaiky <- unique(mozaiky)dir.create("./rasters_biotopes", showWarnings = FALSE)
biotopes_all$BIOTOP_SEZ <- trimws(biotopes_all$BIOTOP_SEZ)for (i in 1:length(unique(mozaiky$nazev.biotopu))) {
temp <- mozaiky[mozaiky$nazev.biotopu == unique(mozaiky$nazev.biotopu)[i], ]
biot <- raster::rasterize(biotopes_all[grep(temp$biotop[1], biotopes_all$BIOTOP_SEZ, fixed = TRUE), ],
blank_raster, field = temp$proporce[1], background = 0, update = TRUE)
if (nrow(temp) > 1) {
for (y in 2:nrow(temp)) {
biot2 <- raster::rasterize(biotopes_all[grep(temp$biotop[y], biotopes_all$BIOTOP_SEZ, fixed = TRUE), ],
blank_raster, field = temp$proporce[y], background = 0, update = TRUE)
biot <- mean(biot + biot2)
}
}
names(biot) <- paste(unique(mozaiky$nazev.biotopu)[i])
writeRaster(biot, paste0("./rasters_biotopes/", enc2utf8(unique(mozaiky$nazev.biotopu)[i]), ".tif"), overwrite = TRUE, format = "GTiff")
}
unlink("./rasters_biotopes/*.xml")

## Climate, DEM, Floods & Host Plants (Coarse Scale)
# Climate
current <- stack(list.files("../../Past_climate/current_Chelsea_0.3", full.names = TRUE, pattern = glob2rx("bio_*tif")))
current <- crop(current, kves_all)

# DEM
dem <- stack("../../DEM_CZ/DEM_CZ_WGS84.tif",
"../../DEM_CZ/DEM_CZ_WGS84_slope.tif",
"../../DEM_CZ/DEM_CZ_WGS84_roughness.tif",
"../../DEM_CZ/DEM_CZ_WGS84_northness.tif",
"../../DEM_CZ/DEM_CZ_WGS84_eastness.tif")
dem <- crop(dem, blank_raster)

# Host plant rasterization
host.plant <- read.csv2("host.plant_ndop.csv")
host.plant$X <- as.numeric(host.plant$X)
host.plant$Y <- as.numeric(host.plant$Y)
host.plant <- SpatialPointsDataFrame(host.plant[, 11:12], host.plant, proj4string = CRS("+init=epsg:5514"))
host.plant <- spTransform(host.plant, CRS("+init=epsg:4326"))
hp_ras <- raster::rasterize(host.plant, blank_raster, field = 1, background = 0, update = TRUE)
writeRaster(hp_ras, "./rasters_kves/hp.tif", overwrite = TRUE, format = "GTiff")

### Fine-Scale Dataset Alignment (DivLand / Romportl 100x100m, EPSG:3035)

romportl <- stack(list.files("../../../ESRI/Romportl", pattern = "tif", full.names = TRUE))
crs(romportl) <- CRS("+init=epsg:3035")

blank_raster.romp <- romportl[[1]]
values(blank_raster.romp) <- 0

# Vector transformations to ETRS89-LAEA
biotopes_all.romp <- st_transform(biotopes_all, 3035)
arion_lok.romp <- spTransform(arion_lok, CRS("+init=epsg:3035"))dir.create("./rasters_romp", showWarnings = FALSE)
for (i in 1:length(unique(mozaiky$nazev.biotopu))) {
temp <- mozaiky[mozaiky$nazev.biotopu == unique(mozaiky$nazev.biotopu)[i], ]
biot <- raster::rasterize(biotopes_all.romp[grep(temp$biotop[1], biotopes_all.romp$BIOTOP_SEZ, fixed = TRUE), ],
blank_raster.romp, field = temp$proporce[1], background = 0, update = TRUE)
if (nrow(temp) > 1) {
for (y in 2:nrow(temp)) {
biot2 <- raster::rasterize(biotopes_all.romp[grep(temp$biotop[y], biotopes_all.romp$BIOTOP_SEZ, fixed = TRUE), ],
blank_raster.romp, field = temp$proporce[y], background = 0, update = TRUE)
biot <- mean(biot + biot2)
}
}
names(biot) <- paste(unique(mozaiky$nazev.biotopu)[i])
writeRaster(biot, paste0("./rasters_romp/", unique(mozaiky$nazev.biotopu)[i], ".tif"), overwrite = TRUE, format = "GTiff")
}

# Reproject DEM and Climate to 100m Lambert grid

dem.r <- projectRaster(dem, crs = 3035)
dem.r <- resample(dem.r, blank_raster.romp)
writeRaster(dem.r, filename = "dem_100m.grd", overwrite = TRUE)

climate_r <- projectRaster(current, crs = 3035)
climate_r <- resample(climate_r, blank_raster.romp)
writeRaster(climate_r, filename = "climate_100m.grd", overwrite = TRUE)

save.image("01_data_preprocessed.RData")
print("Data preprocessing completed successfully.")

