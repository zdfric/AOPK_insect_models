# =====================================================================
# Title: 02_model_training.R
# Project: Potential distribution models for Habitats Directive Insects (AOPK)
# Description: Fit MaxEnt and BART models for both coarse (~30 arc-sec)
#              and fine (100m) spatial scales. Safely export model objects.
# Phengaris arion (i.e., arion in the datasets) used as an example, replace
# by other species if you need
# =====================================================================

library(terra)
library(raster)
library(sf)
library(ENMeval)
library(dismo)
library(embarcadero)

options(java.parameters = "-Xmx125000m")
library(rJava)

# Helper function for background sampling
source("../../Erebia/sample_rast.R")
load("01_data_preprocessed.RData")

# =====================================================================
# PART 1: COARSE-SCALE MODELING (~30 arc-sec: MaxEnt & BART)
# =====================================================================

# Assemble coarse predictor raster stack
layers <- stack(list.files(c("./rasters_kves", "./rasters_biotopes"), pattern = "tif", full.names = TRUE))
layers <- addLayer(layers, current)
layers <- addLayer(layers, dem)
layers <- mask(layers, cz_map)

# Sample background points
rp_ext <- sampleRast(mask(blank_raster + 1, cz_map), 4000, prob = TRUE, replace = TRUE)
colnames(rp_ext) <- c("X", "Y")

# ---------------------------------------------------------------------
# MaxEnt Tuning & Fitting via ENMeval
# ---------------------------------------------------------------------
data.enmeval <- ENMevaluate(coordinates(arion_lok), layers, bg.coords = rp_ext, 
                            RMvalues = seq(0.5, 7, 0.5), fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), 
                            method = "jackknife", algorithm = "maxent.jar")

write.table(data.enmeval@results, "enmeval_results_coarse.txt")

# Extract optimal regularization multiplier and feature classes
best_rm <- as.numeric(data.enmeval@results[data.enmeval@results$delta.AICc == 0, ]$rm[1])
best_fc <- as.character(data.enmeval@results[data.enmeval@results$delta.AICc == 0, ]$fc[1])

# Fit final MaxEnt model
nrep <- 10
data.mx <- maxent(layers, coordinates(arion_lok), removeDuplicates = TRUE, a = rp_ext,
                  args = c(paste0("replicates=", nrep), paste0("betamultiplier=", best_rm),
                         "noaddsamplestobackground", "noautofeature", "jackknife"))

write.table(data.mx@results, "arion_maxent_results_coarse.txt", sep = "\t")

# ---------------------------------------------------------------------
# BART Model Fitting (Coarse Scale)
# ---------------------------------------------------------------------
arion_bart <- as.data.frame(rbind(coordinates(arion_lok), rp_ext))
arion_bart$pres <- 0
arion_bart$pres[1:nrow(arion_lok)] <- 1
arion_bart_var <- raster::extract(layers, arion_bart[, 1:2])

# Stepwise selection and fitting
sdm.bart2 <- bart.step(y.data = arion_bart[, "pres"], x.data = as.data.frame(arion_bart_var), full = TRUE)

# Preserve BART state to avoid C++ pointer loss
sdm.bart2_coarse_backup <- sdm.bart2
invisible(sdm.bart2$fit$state)

# see the stats
pdf("arion_BART_results.pdf")
varimp(sdm.bart2,plot=T)
summary(sdm.bart2)
dev.off()

# see variables in the model
colnames(sdm.bart2$varcount)

# extract threshold value
sumBART <- capture.output(summary(sdm.bart2))
cutoff_line <- grep("Cutoff", sumBART, value = TRUE)
BART.threshold <- as.numeric(sub(".*=\\s*", "", cutoff_line))

# see variable responses
grafy.bart<-response_curve_bart(sdm.bart2,as.data.frame(arion_bart_var),colnames(sdm.bart2$varcount))
library(ggplot2)
responses<-list()
for (i in 1:length(colnames(sdm.bart2$varcount))) {
  responses[[i]]<- ggplot(as.data.frame(grafy.bart[[i]]), aes(x = value, y = mean)) + geom_point(size = 3)+
    geom_smooth(aes(ymin = q25, ymax = q975), alpha = 0.2)+
    labs(title=names(grafy[i]),y="predicted value")
}

library(gridExtra)

pdf("arion_responses_BART.pdf")
grid.arrange(grobs = responses)
dev.off()

# compare performances MaxEnt vs BART

pdf("arion_bart_vs_maxent.pdf")
plot(sdm.bart.map,arion.model,xlab="bart",ylab="maxent",maxpixels=ncell(sdm.bart.map))
abline(v=BART.threshold, col="blue") # treshold for BART
abline(h=data.mx@results[rownames(data.mx@results)=="Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",11], col="blue")
dev.off()


# =====================================================================
# PART 2: FINE-SCALE MODELING (100x100m Grid: BART Only)
# =====================================================================

# Assemble fine predictor stack (DivLand + Fine Habitat + DEM + Climate)
romp <- stack(list.files("./rasters_romp", full.names = TRUE, pattern = "tif"))
climate_r <- brick("climate_100m.grd")
dem.r <- brick("dem_100m.grd")

layer.r <- stack(romp, romportl, climate_r, dem.r)
layer.r <- dropLayer(layer.r, "DEM_CZ_WGS84")

# Sample fine-scale background points
cz_map.romp <- spTransform(cz_map, CRS("+init=epsg:3035"))
rp_ext.romp <- sampleRast(mask(blank_raster.romp + 1, cz_map.romp), 4000, prob = TRUE, replace = TRUE)
colnames(rp_ext.romp) <- c("X", "Y")

# Prepare fine-scale BART dataset
arion_bart.r <- as.data.frame(rbind(coordinates(arion_lok.romp), rp_ext.romp))
arion_bart.r$pres <- 0
arion_bart.r$pres[1:nrow(arion_lok.romp)] <- 1
arion_bart_var.r <- raster::extract(layer.r, arion_bart.r[, 1:2])

# Fit fine-scale BART model
sdm.bart2.r <- bart.step(y.data = arion_bart.r[, "pres"], x.data = as.data.frame(arion_bart_var.r), full = TRUE)

# Preserve fine-scale BART state
sdm.bart2.r_backup <- sdm.bart2.r
invisible(sdm.bart2.r$fit$state)

# see the stats
pdf("arion_BART_results_romp.pdf")
varimp(sdm.bart2.r,plot=T)
summary(sdm.bart2.r)
dev.off()

# see variables in the model
colnames(sdm.bart2.r$varcount)

# extract threshold value
sumBART <- capture.output(summary(sdm.bart2.r))
cutoff_line <- grep("Cutoff", sumBART, value = TRUE)
BART.r.threshold <- as.numeric(sub(".*=\\s*", "", cutoff_line))

# see variable responses
grafy.romportl<-response_curve_bart(sdm.bart2.r,as.data.frame(arion_bart_var.r),colnames(sdm.bart2.r$varcount))

library(ggplot2)
responses<-list()
for (i in 1:length(colnames(sdm.bart2.r$varcount))) {
  responses[[i]]<- ggplot(as.data.frame(grafy.romportl[[i]]), aes(x = value, y = mean)) + geom_point(size = 3)+
    geom_smooth(aes(ymin = q25, ymax = q975), alpha = 0.2)+
    labs(title=names(grafy.romportl[i]),y="predicted value")
}

library(gridExtra)

pdf("arion_responses_romp.pdf")
grid.arrange(grobs = responses)
dev.off()

# Save fit models to disk
save(data.mx, sdm.bart2, sdm.bart2_coarse_backup,
     sdm.bart2.r, sdm.bart2.r_backup,
     file = "fitted_models.RData")

print("All models successfully trained and saved to fitted_models.RData.")