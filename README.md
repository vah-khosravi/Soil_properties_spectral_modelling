# Soil Spectroscopy in R

An R workflow for processing Visible–Near Infrared (Vis–NIR) soil spectra and predicting soil properties using classical chemometric methods and modern machine learning algorithms.

The workflow includes spectral preprocessing, outlier detection, calibration/validation sample selection, model development, and standardized model evaluation.

The script has been written to rely primarily on base R whenever possible, using external packages only when specific functionality is required (e.g., spectral preprocessing, machine learning algorithms, and memory-based learning).

---

## Features

This repository includes:

- Data import and preparation
- Spectral preprocessing
  - Reflectance to absorbance conversion
  - Standard Normal Variate (SNV)
  - Moving average smoothing
  - Savitzky–Golay smoothing
  - First derivative transformation
- Spectral outlier detection using PCA and Mahalanobis distance
- Calibration and validation sample selection
  - Simple Random Sampling (SRS)
  - Conditioned Latin Hypercube Sampling (cLHS) (optional)   
- Predictive modelling using:
  - Principal Component Regression (PCR)
  - Partial Least Squares Regression (PLSR)
  - Cubist regression
  - Random Forest regression
  - Support Vector Regression (SVR)
  - Memory-Based Learning (MBL)
  - Artificial Neural Network (ANN)
  - Extreme Gradient Boosting (XGBoost)
- Variable importance visualization
- Standardized model evaluation
  - Coefficient of determination (R²)
  - Root Mean Square Error (RMSE)
  - Ratio of Performance to Deviation (RPD)
  - Ratio of Performance to Interquartile Distance (RPIQ)

---

## Requirements

The workflow was developed in R and relies primarily on base R, with additional packages used only where necessary.

Main packages include:
- prospectr
- signal
- Cubist
- caret
- randomForest
- kernlab
- pls
- resemble
- chls
- scatterplot3d
- nnet
- xgboost

Additional supporting packages are listed at the beginning of the script and are installed automatically if required.

---

## Input data

The script expects a CSV file containing laboratory reference values together with Vis-NIR spectral measurements.

Update the following section to match your local data location:

```r
path <- "C:/your/data/location"

soilspec <- read_csv(file.path(path, "your_data.csv"))
```

---

## Workflow

```
Raw spectra
      │
      ▼
Data preparation
      │
      ▼
Spectral preprocessing
      │
      ▼
Outlier detection
      │
      ▼
Calibration / Validation split
      │
      ▼
Model development
      ├── PCR
      ├── PLSR
      ├── Cubist
      ├── Random Forest
      ├── SVR
      ├── ANN
      ├── XGBoost
      └── MBL
      │
      ▼
Performance evaluation
      │
      ▼
Model comparison
```

---
## Output
The workflow generates:
- Preprocessed spectra
- PCA results
- Outlier identification
- Calibration and validation datasets
- Model predictions
- Variable importance plots (where applicable)
- Calibration and validation scatter plots
- Performance statistics for each model

Each model is evaluated using the same metrics:
- R²
- RMSE
- RPD
- RPIQ
making comparison among algorithms straightforward.

## Repository structure

├── README.md
├── Soil_Spectroscopy.R
├── data/
│   └── your_data.csv
├── figures/
└── LICENSE

---
