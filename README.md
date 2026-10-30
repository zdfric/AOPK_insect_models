# AOPK_insect_models
# Species Distribution Modeling (SDM) Framework for Habitats Directive Insects

This repository contains the workflow developed for the **Nature Conservation Agency of the Czech Republic (AOPK ČR)** to model the potential distribution of 14 insect species of European importance (Habitats Directive).

The pipeline supports multi-scale modeling using two primary machine learning techniques:
- **MaxEnt** (Maximum Entropy)
- **BART** (Bayesian Additive Regression Trees)

---

## 📐 Multi-Scale Spatial Architecture

The modeling approach is split across two spatial scales:

1. **Regional / Coarse Scale (~30 arc-seconds / ~1 km)**
   - **Predictors:** CHELSA Climate (Bio1–Bio19), DEM & topographic indices, KVES (Consolidated Ecosystem Layer), and Habitat Mapping layers.
   - **Algorithms:** **MaxEnt** & **BART** (`embarcadero`).

2. **Local / Fine Scale (100 x 100 m Grid)**
   - **Predictors:** High-resolution DivLand landscape metrics ("Romportl" dataset), fine-resolution DEM/topographic variables, and targeted biotope layers (ETRS89-LAEA / EPSG:3035).
   - **Algorithms:** **BART** only (MaxEnt is omitted at this level due to excessive background point requirements and computational constraints).

---

## 📁 Repository Structure & Workflow Order

The workflow is modularized into four sequential R scripts to separate data processing, model fitting, spatial prediction, and model evaluation.

```
├── 01_data_preprocessing.R     # Loads, aligns, crops, and rasterizes environment/species data
├── 02_model_training.R         # Fits MaxEnt (coarse) and BART (coarse + fine) models
├── 03_prediction_mapping.R     # Predicts models onto spatial grids and generates map outputs
├── 04_model_validation.R       # Evaluates models against new survey seasons & field polygons
└── README.md                   # Repository documentation
```

### Script Overview

#### `01_data_preprocessing.R`
- Reads occurrence points from NDOP/internal databases.
- Cleans and standardizes Czech character encoding issues (UTF-8 / CP1250 / Latin).
- Prepares predictor raster stacks for both coarse (WGS84) and fine-resolution (EPSG:3035) scales.
- Rasterizes polygon layer mosaics (KVES & mapping biotopes) into unified single-band raster files.

#### `02_model_training.R`
- Samples background pseudo-absences across study domains.
- Optimizes MaxEnt parameters (Regularization Multiplier & Feature Classes via `ENMeval`).
- Fits stepwise BART models (`embarcadero::bart.step`).
- **BART State Preservation:** Implements robust memory/state persistence procedures (`sdm.bart2$fit$state`) to prevent memory loss during object serialization before saving to `.RData`.

#### `03_prediction_mapping.R`
- Loads saved model objects (`.RData`) independently of the training environment.
- Applies models to current layers or **new/updated data versions** (e.g., newer KVES releases).
- Handles tile-based / chunked spatial predictions (`crop`/`predict`/`merge`) to prevent RAM exhaustion during fine-scale (100m) BART predictions.
- Exports spatial predictions as GeoTIFFs and PDF visualization maps.

#### `04_model_validation.R`
- Evaluates predicted suitability maps against new field survey data (e.g., recent field seasons).
- Calculates spatial cross-correlations (`corLocal`) and threshold metrics (e.g., Maximum Training Sensitivity + Specificity Cloglog).
- Intersects suitability thresholds with habitat types (KVES) and spatial buffers.

---

## ⚙️ Installation & Requirements

Ensure you have R (>= 4.0.0) installed alongside Java (for MaxEnt). High RAM (96GB+) is recommended for the 100m BART modeling step, but at least 32GB for viewing the RData files.

### R Package Dependencies:
```R
install.packages(c("terra", "sf", "dismo", "rJava", 
                   "ENMeval", "ggplot2", "gridExtra", "vegan", 
                   "spatialEco", "gstat", "doSNOW", "parallel"))

# Install embarcadero for BART modeling
devtools::install_github("cjcarlson/embarcadero")

# Install GLOSSA for BART response curves (optional)
devtools::install_github("cjcarlson/glossa")
```

---

## 🚀 How to Run

1. **Prepare Data Inputs:** Place input shapefiles/GPKGs and rasters in their designated directory structures as referenced in `01_data_preprocessing.R`.
2. **Execute Data Pipeline:**
   ```R
   source("01_data_preprocessing.R")
   source("02_model_training.R")
   ```
3. **Re-project / Predict on New Data:**
   When updated habitat or environmental layers become available, run:
   ```R
   source("03_prediction_mapping.R")
   ```
4. **Validation:**
   To cross-evaluate predictions with new field records:
   ```R
   source("04_model_validation.R")
   ```

---

## 🔒 Notes on BART & MaxEnt Persistence in R

- **BART (`embarcadero`):** BART models stored in R objects can lose their C++ external pointers (`$fit$state`) when saved/reloaded. This workflow addresses this issue by storing backup states (`sdm.bart2.zaloha`) and calling `invisible(sdm.bart2$fit$state)` upon reloading to rebuild the object pointer tree.
- **MaxEnt (`dismo`):** MaxEnt model instances are saved directly inside `.RData` files alongside tuned hyperparameters.
