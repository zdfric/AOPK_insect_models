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
