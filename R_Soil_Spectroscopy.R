###################### Soil properties spectral modeling ################


# Spectroscopic prediction of soil properties
#
# Workflow:
# 1. Data import
# 2. Spectral preprocessing
# 3. Outlier detection
# 4. Calibration and validation sampling
# 5. Principal Component Regression (PCR)
# 6. Partial Least Squares Regression (PLSR)
# 7. Cubist
# 8. Random Forest
# 9. Support Vector Regression (SVR)
# 10. Memory-Based Learning (MBL)
# 11. Artificial Neural Network (ANN)
# 12. Extreme Gradient Boosting (XGBoost)

################# Clear data, values, graphs and console ##################

rm(list=ls(all=TRUE))
graphics.off()
shell("cls")

############################ Loading required packages #############################
library(readr)
library(prospectr)
library(signal)
library(scatterplot3d)
library(clhs)
library(resemble)
library(pls)
library(caret)
library(xgboost)
library(readxl)
library(writexl)

# Install packages manually if necessary, e.g.
# install.packages("caret")
# install.packages("Cubist")

#################### Function to calculate model performance statistics #######

evaluate_model <- function(obs, pred, dataset = "") {
  rmse <- sqrt(mean((obs - pred)^2))
  r2   <- cor(obs, pred)^2
  rpd  <- sd(obs) / rmse
  rpiq <- IQR(obs) / rmse
  
  cat(
    sprintf(
      "%s\nR² = %.3f; RMSE = %.3f; RPD = %.3f; RPIQ = %.3f\n",
      dataset, r2, rmse, rpd, rpiq
    )
  )
}

############################ Plot predictions #############################

plot_predictions <- function(obs, pred, title) {
  
  plot(
    obs,
    pred,
    main = title,
    xlab = "Observed SOC",
    ylab = "Predicted SOC",
    xlim = c(0, 4),
    ylim = c(0, 4),
    pch = 16,
    cex = 0.7,
    col = rgb(0.1, 0.1, 0.8, 0.8)
  )
  
  abline(0, 1, col = "red", lwd = 2)
  
}

################### Data import and preparation ####################

# read the .csv file

path <- "C:/Enter/Your/working/directory/here"

soil_data <- read_csv(file.path(path, "Spectra.csv"))
# soil_data <- read_excel(file.path(path, "Spectra.xlsx"))

# Extract spectral variables
spec <- soil_data[, 4:ncol(soil_data)]

# Remove spectral columns from the main data frame
soil_data[, 4:ncol(soil_data)] <- NULL

# Store spectra as a matrix column
soil_data$spc <- spec

########################## Spectral Transformations ########################

# Convert reflectance to absorbance
soil_data$spcA <- log(1 / soil_data$spc)

# Standard Normal Variate (SNV) correction
soil_data$spcASnv <- standardNormalVariate(soil_data$spcA)

# Moving average smoothing
soil_data$spcAMovav <- movav(soil_data$spcASnv, w = 11)

# Savitzky–Golay smoothing
soil_data$spcAMovavSG <- savitzkyGolay(
  soil_data$spcAMovav,
  m = 0,
  p = 2,
  w = 11
)

# First derivative
soil_data$spcDerive1 <- savitzkyGolay(
  soil_data$spcA,
  m = 1,
  p = 2,
  w = 11
)

# get wavelengths directly from column names
wavelength <- as.numeric(colnames(soil_data$spc))

# assign them back (optional, but keeps consistency)
colnames(soil_data$spc) <- wavelength

# plot
matplot(x = wavelength, y = t(soil_data$spc),
        xlab = "Wavelength /nm",
        ylab = "First derivative",
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.3))

######################## Outliers detection ###############################

X <- as.matrix(soil_data$spcDerive1)

# keep only columns with at least one finite value
good_cols <- apply(X, 2, function(x) any(is.finite(x)))

sum(!good_cols)   # how many bad columns
which(!good_cols) # which columns are bad

X_clean <- X[, good_cols, drop = FALSE]

# choose the amount of variance explained
maxexplvar <- 0.99

# PCA
pca <- prcomp(
  X_clean,
  center = TRUE,
  scale. = FALSE
)

# proportion of explained variance
var_exp <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var_exp)

# number of PCs explaining 99 % of variance
nPC <- which(cum_var >= maxexplvar)[1]

# Store PCA results in a list
pcspectraA <- list(
  scores = pca$x[,1:nPC,drop=FALSE],
  loadings = pca$rotation[,1:nPC,drop=FALSE],
  n_components = nPC,
  model = pca
)

# calculate the average of the PC scores
pcspectraACentre <- colMeans(pcspectraA$scores)

# since the result of the 'colMeans' function is a vector
# reformat it to a matrix of 1 row
pcspectraACentre <- t(as.matrix(pcspectraACentre))

# compute Mahalanobis distance between scores centre and the scores of the spectra
# covariance matrix of PCA scores
S <- cov(pcspectraA$scores)

# Mahalanobis distance from the PCA center
wmahald <- mahalanobis(
  x = pcspectraA$scores,
  center = as.numeric(pcspectraACentre),
  cov = S
)

# Square root of Mahalanobis distance
wmahald <- sqrt(wmahald)

# plot the index of the spectra against the Mahalanobis distance
plot(wmahald,
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.9),
     ylab = "Mahalanobis distance")

# add a horizontal line
# visualize the spectra with Mahalanobis dissimilarity scores larger than 3
# (arbitrary threshold)
abline(h = 10, col = "red")

# obtain the indices of the outliers
indxOutM <- which(wmahald > 10)

# how many potential outliers?
length(indxOutM)

# plot the first three PCs along the identified outliers
sct3d <-scatterplot3d(pcspectraA$scores[,1],
                      pcspectraA$scores[,2],
                      pcspectraA$scores[,3],
                      xlab = "PC 1",
                      ylab = "PC 2",
                      zlab = "PC 3",
                      color=rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
                      pch = 16, grid=TRUE, angle=50)
# add location of the outliers in red
sct3d$points3d(pcspectraA$scores[indxOutM,1],
               pcspectraA$scores[indxOutM,2],
               pcspectraA$scores[indxOutM,3],
               pch = "X",
               col = "red")

soil_data <- soil_data[-c(indxOutM), ]
spec <- soil_data$spc


### compute the H distance ()
hs <- (wmahald^2 * pcspectraA$n_components)^0.5

# plot the index of the spectra against the H distance
plot(hs,
      pch = 16,
      col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
      ylab = "H distance")
# add a horizontal line
# visualize the spectra with H distance larger than 6 (arbitrary threshold)
abline(h = 80, col = "red")

# obtain the indices of the outliers
 indxOutH <- which(hs > 80)
# how many potential outliers?
length(indxOutH)

# plot the first three PCs along the identified outliers
sct3d <-scatterplot3d(pcspectraA$scores[,1],
                      pcspectraA$scores[,2],
                       pcspectraA$scores[,3],
                       xlab = "PC 1",
                       ylab = "PC 2",
                       zlab = "PC 3",
                       color=rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
                       pch = 16, grid=TRUE, angle=50)
sct3d$points3d(pcspectraA$scores[indxOutH,1],
               pcspectraA$scores[indxOutH,2],
               pcspectraA$scores[indxOutH,3],
               pch = "X",
               col = "red")

soil_data <- soil_data[-c(indxOutM), ]

############################### Sampling methods ##########################

### Simple Random Sampling
sample_size <- round(0.75 * nrow(soil_data))

# set the seed
set.seed(19101991)

# id of the rows to be used for calibration
cal_id <- sample(1:nrow(soil_data),
                size = sample_size)

# separate the dataset into calibration and validation
datC <- soil_data[cal_id,]
datV <- soil_data[-cal_id,]

# plot the value of the Total Carbon content for both calibration and validation
par(mfrow=c(1,2))

# calibration
hist(datC$`SOC`,
     main = "",
     xlab = "SOC")

# validation
hist(datV$`SOC`,
     main = "",
     xlab = "SOC")

### Conditioned Latin Hypercube Sampling

# choose the amount of variance explained
maxexplvar <- 0.99

# compute the principal components
X <- as.matrix(soil_data$spcDerive1)

good_cols <- apply(X,2,function(z) any(is.finite(z)))
X <- X[,good_cols]

pca <- prcomp(
  X,
  center = TRUE,
  scale. = FALSE
)

var_exp <- pca$sdev^2/sum(pca$sdev^2)
cum_var <- cumsum(var_exp)

nPC <- which(cum_var >= maxexplvar)[1]

pcspectraA <- list(
  scores = pca$x[,1:nPC,drop=FALSE],
  loadings = pca$rotation[,1:nPC,drop=FALSE],
  n_components = nPC,
  model = pca
)

sample_size <- round(0.75*nrow(pcspectraA$scores))

# set the seed
set.seed(19101991)

# since the clhs function accepts only 'data.frame' objects as input variables, we can
# transform our matrix of scores to 'data.frame' using the as.data.frame function
clhsS <- clhs(x = as.data.frame(pcspectraA$scores),
              size = sample_size,
              iter = 1000,
              simple = FALSE)

# verify visually that the objective function was correctly minimized
plot(clhsS)

# summary of the clhs object created with the clhs function
str(clhsS)

# plot the first three PCs
sct3d <-scatterplot3d(pcspectraA$scores[,1],
                      pcspectraA$scores[,2],
                      pcspectraA$scores[,3],
                      xlab = "PC 1",
                      ylab = "PC 2",
                      zlab = "PC 3",
                      color= rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
                      pch = 16, grid=TRUE, angle=50)

# plot the selected calibration sample
sct3d$points3d(pcspectraA$scores[,1][clhsS$index_samples],
               pcspectraA$scores[,2][clhsS$index_samples],
               pcspectraA$scores[,3][clhsS$index_samples],
               pch = 1,
               col = "red",
               cex = 1.5)

# separate the dataset into calibration and validation
cal_id <- clhsS$index_samples
datC <- soil_data[cal_id,]
datV <- soil_data[-cal_id,]

# plot the value of the Total Carbon content for both calibration and validation
par(mfrow=c(1,2))

# calibration
hist(datC$SOC,
     main = "",
     xlab = "SOC")

# validation
hist(datV$SOC,
     main = "",
     xlab = "SOC")

######################## Principal Component Regression ####################

# Perform Principal Component Analysis (PCA) on the calibration spectra
pcspectra <- prcomp(
  datC$spcDerive1,
  center = TRUE,
  scale. = TRUE
)

# Calculate the cumulative percentage of explained variance
variance <- pcspectra$sdev^2
cum_variance <- 100 * cumsum(variance) / sum(variance)

# Plot the cumulative explained variance
plot(
  cum_variance[1:30],
  type = "b",
  xlab = "Principal component",
  ylab = "Cumulative explained variance (%)"
)

# Number of principal components retained
npc <- 15

# Extract the selected principal component scores
cal_scores <- as.data.frame(pcspectra$x[, 1:npc])

# Fit the PCR model
soilCPCRModel <- lm(
  datC$SOC ~ .,
  data = cal_scores
)

# Display model summary
summary(soilCPCRModel)

# Calibration predictions
soilCPCRPred <- predict(
  soilCPCRModel,
  newdata = cal_scores
)

# Project the validation spectra into the PCA space
val_scores <- predict(
  pcspectra,
  newdata = datV$spcDerive1
)

val_scores <- as.data.frame(val_scores[, 1:npc])

# Validation predictions
soilVPCRPred <- predict(
  soilCPCRModel,
  newdata = val_scores
)

# Plot observed versus predicted values
par(mfrow = c(1,2))

plot_predictions(datC$SOC, soilCPCRPred, "Calibration")
plot_predictions(datV$SOC, soilVPCRPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCPCRPred, "Calibration")
evaluate_model(datV$SOC, soilVPCRPred, "Validation")

######################## Partial Least Squares Regression #####################

# Maximum number of latent variables
max_comp <- 30

# Fit PLSR model using repeated cross-validation
soilCPLSModel <- plsr(
  SOC ~ spcDerive1,
  data = datC,
  method = "oscorespls",
  ncomp = max_comp,
  validation = "CV",
  segments = 10,
  repeats = 5
)

# Cross-validation results
plot(
  soilCPLSModel,
  "validation",
  xlab = "Number of components",
  main = ""
)

# Select the optimal number of components
rmsep <- RMSEP(soilCPLSModel)
ncomp_opt <- which.min(rmsep$val[1, 1, -1])

cat("Optimal number of components:", ncomp_opt, "\n")

# Use the optimal number of components
nc <- ncomp_opt

# Cross-validated predictions
plot(
  soilCPLSModel,
  ncomp = nc,
  xlab = "Observed SOC",
  ylab = "Predicted SOC",
  main = ""
)

# Plot the first three loading vectors
plot(
  soilCPLSModel,
  "loadings",
  comps = 1:3,
  xlab = "Wavelength index",
  ylab = "Loading"
)

# Regression coefficients
plot(
  wavelength,
  soilCPLSModel$coefficients[, 1, nc],
  type = "l",
  xlab = "Wavelength (nm)",
  ylab = "Regression coefficient"
)
abline(h = 0, lty = 2)

# Variable Importance (VIP)

# Extract model matrices
W <- soilCPLSModel$loading.weights
Q <- soilCPLSModel$Yloadings
T_scores <- soilCPLSModel$scores

# Compute Variable Importance in Projection (VIP)
Q2 <- as.numeric(Q)^2
SSY <- Q2[1:nc] * diag(crossprod(T_scores))[1:nc]

Wnorm <- W^2 / colSums(W^2)

vip <- sqrt(
  length(wavelength) *
    rowSums(sweep(Wnorm[, 1:nc], 2, SSY, "*")) /
    sum(SSY)
)

# Plot VIP
plot(
  wavelength,
  vip,
  type = "l",
  col = "grey40",
  xlab = "Wavelength (nm)",
  ylab = "VIP"
)

abline(h = 1, col = "red", lty = 2)

# Predictions

# Calibration
soilCPLSPred <- drop(
  predict(
    soilCPLSModel,
    ncomp = nc,
    newdata = datC$spcDerive1
  )
)

# Validation
soilVPLSPred <- drop(
  predict(
    soilCPLSModel,
    ncomp = nc,
    newdata = datV$spcDerive1
  )
)

# Plot observed versus predicted values
par(mfrow = c(1, 2))

plot_predictions(datC$SOC, soilCPLSPred, "Calibration")
plot_predictions(datV$SOC, soilVPLSPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCPLSPred, "Calibration")
evaluate_model(datV$SOC, soilVPLSPred, "Validation")

############################### Cubist #######################################

# Prepare calibration and validation datasets
datC_cubist <- data.frame(
  SOC = datC$SOC,
  datC$spcDerive1
)
colnames(datC_cubist) <- c(
  "SOC",
  paste0("spec.", colnames(datC$spcDerive1))
)

datV_cubist <- data.frame(
  SOC = datV$SOC,
  datV$spcDerive1
)
colnames(datV_cubist) <- c(
  "SOC",
  paste0("spec.", colnames(datV$spcDerive1))
)

# Train Cubist model
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1
)

cubist_grid <- expand.grid(
  committees = c(1, 100),
  neighbors = c(0, 9)
)

soilCCubistModel <- train(
  SOC ~ .,
  data = datC_cubist,
  method = "cubist",
  tuneGrid = cubist_grid,
  trControl = ctrl
)

summary(soilCCubistModel)

# Predictions
soilCCubistPred <- drop(
  predict(
    soilCCubistModel,
    newdata = datC_cubist[, -1]
  )
)

soilVCubistPred <- drop(
  predict(
    soilCCubistModel,
    newdata = datV_cubist[, -1]
  )
)

# Plot observed versus predicted values
par(mfrow = c(1, 2))

plot_predictions(datC$SOC, soilCCubistPred, "Calibration")
plot_predictions(datV$SOC, soilVCubistPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCCubistPred, "Calibration")
evaluate_model(datV$SOC, soilVCubistPred, "Validation")

# Variable usage

# Variable usage statistics
cubist_usage <- soilCCubistModel$finalModel$usage

cubist_usage$total <- cubist_usage$Conditions +
  cubist_usage$Model

cubist_usage$wavelength <- as.numeric(
  sub("spec\\.", "", cubist_usage$Variable)
)

used <- subset(cubist_usage, total > 0)

# Plot spectrum and variable usage
spec_wavelength <- as.numeric(
  colnames(datC$spcDerive1)
)

plot(
  spec_wavelength,
  datC$spcDerive1[1, ],
  type = "l",
  col = "grey50",
  xlab = "Wavelength (nm)",
  ylab = "First derivative",
  xlim = range(spec_wavelength)
)

par(new = TRUE)

plot(
  used$wavelength,
  used$total,
  type = "h",
  axes = FALSE,
  xlab = "",
  ylab = "",
  col = "plum",
  lwd = 3,
  xlim = range(spec_wavelength),
  ylim = c(0, max(used$total))
)

used_conditions <- subset(
  used,
  Conditions > 0
)

lines(
  used_conditions$wavelength,
  used_conditions$Conditions,
  type = "h",
  col = "blue",
  lwd = 3
)

axis(4)

mtext(
  "Cubist usage",
  side = 4,
  line = 3
)

legend(
  "topright",
  legend = c(
    "Spectrum",
    "Model usage",
    "Condition usage"
  ),
  col = c("grey50", "plum", "blue"),
  lwd = 2,
  bty = "n"
)

########################## Random Forest #################################

# Keep only spectral variables with finite values
Xc <- as.matrix(datC$spcDerive1)
good_cols <- apply(Xc, 2, function(x) all(is.finite(x)))
Xc <- Xc[, good_cols, drop = FALSE]

# Prepare calibration data
datCSub <- data.frame(
  SOC = datC$SOC,
  Xc
)

colnames(datCSub) <- c(
  "SOC",
  paste0("spec.", colnames(Xc))
)

# Prepare validation data
Xv <- as.matrix(datV$spcDerive1)
Xv <- Xv[, good_cols, drop = FALSE]

datVSub <- data.frame(
  SOC = datV$SOC,
  Xv
)

colnames(datVSub) <- c(
  "SOC",
  paste0("spec.", colnames(Xv))
)

# Train Random Forest model
set.seed(19101991)

soilCRFModel <- train(
  SOC ~ .,
  data = datCSub,
  method = "rf",
  tuneLength = 5,
  trControl = trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 1
  )
)

# Model summary
print(soilCRFModel)

# Variable importance
rfImportance <- varImp(soilCRFModel)

print(rfImportance)

# Plot the 30 most important wavelengths
plot(rfImportance, top = 30)

# Predictions
soilCRFPred <- predict(soilCRFModel, datCSub)
soilVRFPred <- predict(soilCRFModel, datVSub)


# Plot observed versus predicted values
par(mfrow = c(1,2))

plot_predictions(datC$SOC, soilCRFPred, "Calibration")
plot_predictions(datV$SOC, soilVRFPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCRFPred, "Calibration")
evaluate_model(datV$SOC, soilVRFPred, "Validation")

######################## Support Vector Regression ########################

# Prepare calibration data
datCSub <- data.frame(
  SOC = datC$SOC,
  datC$spcDerive1
)

colnames(datCSub) <- c(
  "SOC",
  paste0("spec.", colnames(datC$spcDerive1))
)

# Prepare validation data
datVSub <- data.frame(
  SOC = datV$SOC,
  datV$spcDerive1
)

colnames(datVSub) <- c(
  "SOC",
  paste0("spec.", colnames(datV$spcDerive1))
)

# Train SVR model
set.seed(19101991)

train.control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1
)

# Hyperparameter grid
svmGrid <- expand.grid(
  sigma = 10^(-5:-1),
  C = 10^(-3:1)
)

soilCSVRModel <- train(
  SOC ~ .,
  data = datCSub,
  method = "svmRadial",
  tuneGrid = svmGrid,
  trControl = train.control
)

# Model summary
print(soilCSVRModel)

# Predictions
soilCSVRPred <- predict(soilCSVRModel, datCSub)
soilVSVRPred <- predict(soilCSVRModel, datVSub)

# Calibration and validation plots
par(mfrow = c(1, 2))

plot_predictions(datC$SOC, soilCSVRPred, "Calibration")
plot_predictions(datV$SOC, soilVSVRPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCSVRPred, "Calibration")
evaluate_model(datV$SOC, soilVSVRPred, "Validation")

########################## Memory Based Learning (MBL) ####################

# Define the neighbourhood sizes to evaluate
k2t <- seq(20, 40, by = 1)

# Nearest-neighbour optimization
soilCMBLModel <- mbl(
  Xr = datC$spcDerive1,
  Yr = datC$SOC,
  Xu = datV$spcDerive1,
  Yu = NULL,
  
  # PCA-based dissimilarity
  diss_method = diss_pca(
    center = TRUE,
    scale = FALSE
  ),
  
  # Nearest-neighbour validation
  control = mbl_control(
    validation_type = "NNv"
  ),
  
  # Test different neighbourhood sizes
  neighbors = neighbors_k(
    k = k2t
  ),
  
  # Weighted adaptive PLS
  fit_method = fit_wapls(
    min_ncomp = 4,
    max_ncomp = 17
  )
)

# Display model summary
print(soilCMBLModel)

# Select the optimal neighbourhood size
valResults <- soilCMBLModel$validation_results$nearest_neighbor_validation

rmseMBL <- valResults$rmse
neighNumber <- valResults$k

plot(
  neighNumber,
  rmseMBL,
  type = "b",
  pch = 16,
  xlab = "Number of neighbours",
  ylab = "RMSE"
)

optNn <- neighNumber[which.min(rmseMBL)]

cat(
  sprintf(
    "Optimal number of neighbours: %d\n",
    optNn
  )
)

# Final MBL model
soilVMBLModel <- mbl(
  Xr = datC$spcDerive1,
  Yr = datC$SOC,
  
  Xu = datV$spcDerive1,
  Yu = datV$SOC,
  
  diss_method = diss_pca(
    center = TRUE,
    scale = FALSE
  ),
  
  neighbors = neighbors_k(
    k = optNn
  ),
  
  fit_method = fit_wapls(
    min_ncomp = 4,
    max_ncomp = 17
  )
)

# Predictions
soilVMBLPred <-
  soilVMBLModel$results[[paste0("k_", optNn)]]$pred

# Validation plots
par(mfrow = c(1,2))

plot_predictions(datV$SOC, soilVMBLPred, "Validation")

# Model evaluation
evaluate_model(datV$SOC, soilVMBLPred, "Validation")


###################### Artificial Neural Network (ANN) ###########################
# ANN can become computationally expensive when using thousands of
# highly correlated spectral variables. To improve efficiency, Principal
# Component Analysis (PCA) is first applied to the calibration spectra, and
# the first principal components explaining 99% of the spectral variance are
# used as predictors. This substantially reduces the number of variables,
# decreases computation time, and minimizes multicollinearity while
# preserving nearly all of the spectral information.

# Perform PCA on calibration spectra
pcaANN <- prcomp(
  datC$spcDerive1,
  center = TRUE,
  scale. = TRUE
)

# Percentage of explained variance
var_exp <- pcaANN$sdev^2 / sum(pcaANN$sdev^2)
cum_var <- cumsum(var_exp)

npc <- which(cum_var >= 0.99)[1]

# Plot cumulative explained variance
plot(
  cum_var,
  type = "b",
  pch = 16,
  xlab = "Principal Component",
  ylab = "Cumulative variance explained (%)"
)
abline(v = npc, lty = 2, col = "red")

# Create calibration dataset
datCANN <- data.frame(
  SOC = datC$SOC,
  pcaANN$x[, 1:npc]
)

# Project validation spectra onto calibration PCA space
scoresV <- predict(
  pcaANN,
  newdata = datV$spcDerive1
)

datVANN <- data.frame(
  SOC = datV$SOC,
  scoresV[, 1:npc]
)

# Hyperparameter grid
annGrid <- expand.grid(
  size = c(3, 5, 7),
  decay = c(0, 0.001, 0.01)
)

# Cross-validation
set.seed(19101991)

soilCANNModel <- train(
  SOC ~ .,
  data = datCANN,
  method = "nnet",
  tuneGrid = annGrid,
  linout = TRUE,
  trace = FALSE,
  maxit = 500,
  trControl = trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 1
  )
)

# Model summary
print(soilCANNModel)

# Best tuning parameters
print(soilCANNModel$bestTune)


# Predictions
soilCANNPred <- predict(soilCANNModel, datCANN)

soilVANNPred <- predict(soilCANNModel, datVANN)

# Plot observed versus predicted values
par(mfrow = c(1,2))

plot_predictions(datC$SOC, soilCANNPred, "Calibration")
plot_predictions(datV$SOC, soilVANNPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCANNPred, "Calibration")
evaluate_model(datV$SOC, soilVANNPred, "Validation")

#################### Extreme Gradient Boosting (XGBoost) #########################
# XGBoost can become computationally expensive when using thousands of
# highly correlated spectral variables. To improve efficiency, Principal
# Component Analysis (PCA) is first applied to the calibration spectra, and
# the first principal components explaining 99% of the spectral variance are
# used as predictors. This substantially reduces the number of variables,
# decreases computation time, and minimizes multicollinearity while
# preserving nearly all of the spectral information.

# Perform PCA on calibration spectra
pca <- prcomp(
  datC$spcDerive1,
  center = TRUE,
  scale. = TRUE
)

# Percentage of explained variance
cumvar <- cumsum(pca$sdev^2) / sum(pca$sdev^2)

# Number of PCs explaining at least 99% of variance
npc <- which(cumvar >= 0.99)[1]

cat("Number of principal components used:", npc, "\n")

# Calibration scores
Xcal <- pca$x[, 1:npc]

# Validation scores
Xval <- predict(
  pca,
  newdata = datV$spcDerive1
)[, 1:npc]

# XGBoost data

datCSub <- xgb.DMatrix(
  data = Xcal,
  label = datC$SOC
)

datVSub <- xgb.DMatrix(
  data = Xval,
  label = datV$SOC
)

# Parameter selection

params <- list(
  
  booster = "gbtree",
  
  objective = "reg:squarederror",
  
  eta = 0.05,
  
  max_depth = 4,
  
  min_child_weight = 1,
  
  subsample = 0.8,
  
  colsample_bytree = 0.7
)

# Cross-validation

set.seed(19101991)

cv <- xgb.cv(
  
  params = params,
  
  data = datCSub,
  
  metrics = "rmse",
  
  nrounds = 500,
  
  nfold = 5,
  
  early_stopping_rounds = 20,
  
  verbose = 0
)

best_nrounds <- cv$early_stop$best_iteration

cat("Optimal number of boosting iterations =", best_nrounds, "\n")

# Final model

soilCXGBModel <- xgb.train(
  params = params,
  data = datCSub,
  nrounds = best_nrounds,
  verbose = 0
)

print(soilCXGBModel)

# Variable importance

importance <- xgb.importance(
  
  feature_names = colnames(Xcal),
  
  model = soilCXGBModel
  
)

print(head(importance))

xgb.plot.importance(
  importance_matrix = importance,
  top_n = 20
)

# Predictions
soilCXGBPred <- predict(soilCXGBModel, datCSub)
soilVXGBPred <- predict(soilCXGBModel, datVSub)

# Calibration and validation plots
par(mfrow = c(1,2))

plot_predictions(datC$SOC, soilCXGBPred, "Calibration")
plot_predictions(datV$SOC, soilVXGBPred, "Validation")

# Model evaluation
evaluate_model(datC$SOC, soilCXGBPred, "Calibration")
evaluate_model(datV$SOC, soilVXGBPred, "Validation")

