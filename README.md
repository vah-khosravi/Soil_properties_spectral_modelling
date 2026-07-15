# Soil Spectroscopy Workflow in R

An R workflow for processing Visible–Near Infrared (Vis-NIR) soil spectra and predicting Soil Organic Carbon (SOC) using several machine learning (ML) approaches.

The workflow includes spectral preprocessing, outlier detection, calibration and validation sample selection, model development, and model evaluation.

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
- Principal Component Regression (PCR)
- Partial Least Squares Regression (PLSR)
- Cubist regression
- Random Forest regression
- Support Vector Regression (SVR)
- Memory-Based Learning (MBL)
- Model evaluation using standard performance metrics

---

## Requirements

The workflow was developed in R and relies on several packages, including:

- soilspec
- prospectr
- Cubist
- caret
- randomForest
- pls
- resemble
- ggplot2

along with additional supporting packages listed in the script.

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
      └── MBL
      │
      ▼
Model evaluation
```

---

## Acknowledgements

This workflow was developed using functions and methodologies available in the **soilspec** package and other open-source R packages for soil spectroscopy and chemometric analysis. The implementation also benefited from the examples and documentation provided with the **soilspec** package.

---
