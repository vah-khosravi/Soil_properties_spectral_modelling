###############################################################################
# Spectroscopic prediction of soil organic carbon (SOC)
#
# This script reproduces the workflow presented in
# "A Practical Guide to Soil Visible and Near Infrared Spectroscopy".
#
# Workflow:
# 1. Data import
# 2. Spectral preprocessing
# 3. Outlier detection
# 4. Calibration/validation sampling
# 5. PCR
# 6. PLSR
# 7. Cubist
# 8. Random Forest
# 9. Support Vector Regression
# 10. Memory-Based Learning
###############################################################################

######################## Installing required packages #####################

# check if the devtools package is already installed
if (!require("devtools")) install.packages("devtools")

# install the soilspec package from GitHub
# devtools::install_github("AlexandreWadoux/soilspec")

# specify all the packages used in the book
myPackages <- c("asdreader", "RColorBrewer","wavethresh", "MASS",
                "pracma", "plyr", "signal", "SDMTools",
                "tripack", "resemble", "splancs", "sp",
                "scatterplot3d", "prospectr", "RcppArmadillo", "matrixStats",
                "clhs", "viridis", "viridisLite", "caret",
                "ggplot2", "soiltexture","randomForest","pls",
                "Cubist", "lattice", "robustbase", "mvoutlier",
                "devtools","arrow", "feather", "fst", "hexView",
                "pzfx", "readODS", "rmarkdown", "rmatio", "rio", "kernlab", "chillR",
                "beepr", "ggpubr", "Metrics", "ModelMetrics", "parallel", "caTools")

# define which packages are not installed in the current computer
notInstalled <- myPackages[!(myPackages%in%installed.packages()[ , "Package"])]

# install the missing packages
if(length(notInstalled)>0) install.packages(notInstalled)

################# Clear data, values, graphs and console ##################

rm(list=ls(all=TRUE))
graphics.off()
shell("cls")

################### Data import and preparation ####################

# Loading required packages
library(prospectr)
library(soilspec)
library(rio)
library(signal)
library(plyr)
library(resemble)
library(scatterplot3d)
library(clhs)
library(pls)
library(Cubist)
library(randomForest)
library(e1071)
library(caTools)
library(caret)
library(parallel)
library(ModelMetrics)
library(Metrics)
library(ggpubr)
library(ggplot2)
library(lattice)
library(beepr)
library(chillR)
library(kernlab)
library(readr)

# read the .csv file

path <- setwd("C:/Users/uzivatel/OneDrive - CZU v Praze/Personal/iPhD/Data/exdata")

soilspec <- read_csv(file.path(path, "All_Lab.csv"))


# put the spectra into a single dataframe
spec <- soilspec[,4:2154]

# remove the spectra from the current dataframe
soilspec [, 4:2154] <- NULL

# add the spectra to the dataframe
soilspec$spc <- spec

########################## Spectral Transformations ########################

# change from reflectance to absorbance
soilspec$spcA <- log(1/soilspec$spc)

# apply a standard normal variate transformation for baseline correction
soilspec$spcASnv <- standardNormalVariate(soilspec$spcA)

# apply a moving average window to the standard normal variate spectra
soilspec$spcAMovav <- movav(soilspec$spcASnv, w = 11)

# Savitzky-Golay smoothing
soilspec$spcAMovavSG <- savitzkyGolay(soilspec$spcAMovav, 0, 2, 11)

# Savitzky-Golay 1st derivative
soilspec$spcDerive1 <- savitzkyGolay(soilspec$spcA, 1, 2, 11)

# get wavelengths directly from column names
wavelength <- as.numeric(colnames(soilspec$spcDerive1))

# assign them back (optional, but keeps consistency)
colnames(soilspec$spcDerive1) <- wavelength

# plot
matplot(x = wavelength, y = t(soilspec$spcDerive1),
        xlab = "Wavelength /nm",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.3))

######################## Outliers detection ###############################

X <- as.matrix(soilspec$spcDerive1)

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

# object mimicking the old soilspec output
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

# print the results
print(pcspectraACentre)

# compute Mahalanobis distance between scores centre and the scores of the spectra
# covariance matrix of PCA scores
S <- cov(pcspectraA$scores)

# Mahalanobis distance from the PCA center
wmahald <- mahalanobis(
  x = pcspectraA$scores,
  center = as.numeric(pcspectraACentre),
  cov = S
)

# Match the scale of the original soilspec function
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

soilspec <- soilspec[-c(indxOutM), ]
spec <- soilspec$spc


### compute the H distance ()
# hs <- (wmahald^2 * pcspectraA$n_components)^0.5
# 
# # plot the index of the spectra against the H distance
# plot(hs,
#      pch = 16,
#      col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
#      ylab = "H distance")
# # add a horizontal line
# # visualize the spectra with H distance larger than 6 (arbitrary threshold)
# abline(h = 80, col = "red")
# 
# # obtain the indices of the outliers
# indxOutH <- which(hs > 80)
# # how many potential outliers?
# length(indxOutH)
# 
# # plot the first three PCs along the identified outliers
# sct3d <-scatterplot3d(pcspectraA$scores[,1],
#                       pcspectraA$scores[,2],
#                       pcspectraA$scores[,3],
#                       xlab = "PC 1",
#                       ylab = "PC 2",
#                       zlab = "PC 3",
#                       color=rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
#                       pch = 16, grid=TRUE, angle=50)
# sct3d$points3d(pcspectraA$scores[indxOutH,1],
#                pcspectraA$scores[indxOutH,2],
#                pcspectraA$scores[indxOutH,3],
#                pch = "X",
#                col = "red")

# soilspec <- soilspec[-c(indxOutM), ]

############################### Sampling methods ##########################

### Simple Random Sampling
SampleSize <- round(0.75*nrow(soilspec))

# set the seed
set.seed(19101991)

# id of the rows to be used for calibration
calId <- sample(1:nrow(soilspec),
                size = SampleSize)

# separate the dataset into calibration and validation
datC <- soilspec[calId,]
datV <- soilspec[-calId,]

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

# ### Conditioned Latin Hypercube Sampling
# 
# # choose the amount of variance explained
# maxexplvar <- 0.99
# 
# # compute the principal components
# X <- as.matrix(soilspec$spcDerive1)
# 
# good_cols <- apply(X,2,function(z) any(is.finite(z)))
# X <- X[,good_cols]
# 
# pca <- prcomp(
#   X,
#   center = TRUE,
#   scale. = FALSE
# )
# 
# var_exp <- pca$sdev^2/sum(pca$sdev^2)
# cum_var <- cumsum(var_exp)
# 
# nPC <- which(cum_var >= maxexplvar)[1]
# 
# pcspectraA <- list(
#   scores = pca$x[,1:nPC,drop=FALSE],
#   loadings = pca$rotation[,1:nPC,drop=FALSE],
#   n_components = nPC,
#   model = pca
# )
# 
# SampleSize <- round(0.75*nrow(pcspectraA$scores))
# 
# # set the seed
# # set.seed(19101991)
# 
# # since the clhs function accepts only 'data.frame' objects as input variables, we can
# # transform our matrix of scores to 'data.frame' using the as.data.frame function
# clhsS <- clhs(x = as.data.frame(pcspectraA$scores),
#               size = SampleSize,
#               iter = 1000,
#               simple = FALSE)
# 
# # verify visually that the objective function was correctly minimized
# plot(clhsS)
# 
# # summary of the clhs object created with the clhs function
# str(clhsS)
# 
# # plot the first three PCs
# sct3d <-scatterplot3d(pcspectraA$scores[,1],
#                       pcspectraA$scores[,2],
#                       pcspectraA$scores[,3],
#                       xlab = "PC 1",
#                       ylab = "PC 2",
#                       zlab = "PC 3",
#                       color= rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
#                       pch = 16, grid=TRUE, angle=50)
# 
# # plot the selected calibration sample
# sct3d$points3d(pcspectraA$scores[,1][clhsS$index_samples],
#                pcspectraA$scores[,2][clhsS$index_samples],
#                pcspectraA$scores[,3][clhsS$index_samples],
#                pch = 1,
#                col = "red",
#                cex = 1.5)
# 
# # separate the dataset into calibration and validation
# calId <- clhsS$index_samples
# datC <- soilspec[calId,]
# datV <- soilspec[-calId,]
# 
# # plot the value of the Total Carbon content for both calibration and validation
# par(mfrow=c(1,2))
# 
# # calibration
# hist(datC$SOC,
#      main = "",
#      xlab = "SOC")
# 
# # validation
# hist(datV$SOC,
#      main = "",
#      xlab = "SOC")

########################## Principal Component Regression ################

# first we perform a PCA on the spectra
pcspectra <- prcomp(datC$spcDerive1,
                    center = TRUE, scale = TRUE)

# calculate the percent of the variances explained by the PC
v <- pcspectra$sdev*pcspectra$sdev

# percentage of cumulative variances
cumv <- 100*cumsum(v)/sum(v)

# plot cumulative percentage of variances explained by the PCs
plot(cumv[1:30],
     type = "b",
     xlab = "PC",
     ylab = "% Cumulative variance")

# specify number of components
npc <- 15

# select PC scores
sdata <- as.data.frame(pcspectra$x[,1:npc])

# fit a linear model SOC = PC1 + PC2 + ...
soilCPcrModel <- lm(datC$SOC ~ ., data = sdata)

# obtain a summary of the fit
summary(soilCPcrModel)

# predict on the calibration dataset
soilCPcrPred <- predict(soilCPcrModel, sdata)

# predict on the validation dataset
pcspectraV <- predict(pcspectra, datV$spcDerive1)
sdataNew <- as.data.frame(pcspectraV[, 1:npc])
soilVPcrPred <- predict(soilCPcrModel, sdataNew)

par(mfrow = c(1, 2))

# plot calibration
plot(datC$SOC, soilCPcrPred,
     main = list("Calibration", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# plot validation
plot(datV$SOC, soilVPcrPred,
     main = list("Validation", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# accuracy measures for calibration
soilspec::eval(datC$SOC, soilCPcrPred, obj = "quant")

# accuracy measures for validation
soilspec::eval(datV$SOC, soilVPcrPred, obj = "quant")

 ########################## Partial Least Square Regression ################

# maximum number of components in the PLS model
maxc <- 30

# generate a PLS model based on calibration data
soilCPlsModel <- plsr(SOC ~ spcDerive1,
                      data = datC,
                      method = "oscorespls",
                      ncomp = maxc,
                      validation = "CV",
                      segments = 10,
                      repeats = 5)

# this is the plsr function, using cross validation to evaluate the RMSEP
# as a function of number of components from one until max
plot(soilCPlsModel, "val",
     main = " ",
     xlab = "Number of components")

# number of components to use
nc <- 5

# plot of cross-validated predictions
plot(soilCPlsModel,
     ncomp = nc,
     main = " ",
     xlab = "Observed",
     ylab = "Predicted")

# the three first loadings
plot(soilCPlsModel,
     "loadings",
     comps = 1:3,
     xlab = "Index of the wavelength",
     ylab = "Loading value")

# plot the coefficient
plot(wavelength, soilCPlsModel$coefficients[,1,nc],
     main = " ",
     type = "l",
     xlab = "Wavelength /nm",
     ylab = "Regression coefficient")
abline(h = 0)

# take the loadings, loading weights and scores
W <- soilCPlsModel$loading.weights
Q <- soilCPlsModel$Yloadings
TT <- soilCPlsModel$scores

# compute the variable importance, see Wold et al., (1993)
Q2 <- as.numeric(Q) * as.numeric(Q)
Q2TT <- Q2[1:nc] * diag(crossprod(TT))[1:nc]
WW <- W * W/apply(W, 2, function(x) sum(x * x))
vip <- sqrt(length(wavelength) * apply(sweep(WW[, 1:nc], 2, Q2TT, "*"),
                                 1, sum)/sum(Q2TT))

# display the variable importance
plot(wavelength, vip,
     xlab = "Wavelength /nm",
     ylab = "Importance",
     type = "l",
     lty = 1,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 1))
abline(h = 1)

# predict on the calibration dataset
soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = datC$spcDerive1)

# predict on the validation dataset
soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = datV$spcDerive1)

par(mfrow = c(1, 2))

# plot calibration
plot(datC$SOC, soilCPlsPred,
     main = list("Calibration", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# plot validation
plot(datV$SOC, soilVplsPred,
     main = list("Validation", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# accuracy measures for calibration
soilspec::eval(datC$SOC, soilCPlsPred, obj = "quant")

# accuracy measures for validation
soilspec::eval(datV$SOC, soilVplsPred, obj = "quant")

################################# Cubist ##################################

# make a Cubist model on calibration dataset

 datCSub <- data.frame(SOC = datC$SOC, datC$spcDerive1)
 colnames(datCSub) <- c("SOC", paste0("spec.", colnames(datC$spcDerive1)))

 
 datVSub <- data.frame(SOC = datV$SOC, datV$spcDerive1)
 colnames(datVSub) <- c("SOC", paste0("spec.", colnames(datV$spcDerive1)))

# Training method

 soilCCubistModel <- train(SOC~., data = datCSub, method = "cubist",
                tuneGrid = expand.grid(.committees = c(1,  100), .neighbors = c(0, 9)),
                trControl = trainControl(method = 'repeatedcv', number = 10, repeats = 1),
                tuneLength = 1) 

# summary of the model
summary(soilCCubistModel)

# predict on the calibration data
 soilCCubistPredict <- predict(soilCCubistModel, datCSub[,2:ncol(datCSub)])
 
# predict on the validation data
 soilVCubistPredict <- predict(soilCCubistModel, datVSub[,2:ncol(datVSub)])

par(mfrow = c(1, 2))

# plot calibration
plot(datCSub$SOC, soilCCubistPredict,
     main = list("Calibration", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# plot validation
plot(datVSub$SOC, soilVCubistPredict,
     main = list("Validation", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 0.5,
     pch = 16)
abline(0, 1, col = "red")

# accuracy measures for calibration
soilspec::eval(datCSub$SOC, soilCCubistPredict, obj = "quant")

# accuracy measure for validation
soilspec::eval(datVSub$SOC, soilVCubistPredict, obj = "quant")

# Plot the variables used as predictors in the Cubist model

# Extract Cubist variable usage
cubist_usage <- soilCCubistModel$finalModel$usage

# Calculate total usage (conditions + model)
cubist_usage$total <- cubist_usage$Conditions + cubist_usage$Model

# Extract wavelength numbers from variable names
cubist_usage$wavelength <- as.numeric(
  sub("spec\\.", "", cubist_usage$Variable)
)

# Keep only variables used by the model
used <- cubist_usage[cubist_usage$total > 0, ]

# Sort by importance
used <- used[order(used$total, decreasing = TRUE), ]


# Plot spectrum

spec_wavelength <- as.numeric(
  sub("spec\\.", "", colnames(datC$spcDerive1))
)

plot(
  spec_wavelength,
  datC$spcDerive1[1, ],
  type = "l",
  col = "grey50",
  xlab = "Wavelength / nm",
  ylab = "Derivative spectrum",
  xlim = c(500, 2450)
)


# Overlay Cubist usage with second axis

par(new = TRUE)

plot(
  used$wavelength,
  used$total,
  type = "h",
  col = "plum",
  lwd = 3,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = c(500, 2450),
  ylim = c(0, max(used$total))
)

used_conditions <- used[used$Conditions > 0, ]

lines(
  used_conditions$wavelength,
  used_conditions$Conditions,
  type = "h",
  col = "blue",
  lwd = 3
)

axis(
  side = 4
)

mtext(
  "Cubist usage",
  side = 4,
  line = 3
)

legend(
  "topright",
  legend = c("Spectrum", "Cubist model usage", "Cubist conditions"),
  col = c("grey50", "plum", "blue"),
  lwd = 2
)

########################## Random Forest #################################
Xc <- as.matrix(datC$spcDerive1)

# keep only complete spectral columns
good_cols <- apply(Xc, 2, function(x) all(is.finite(x)))
Xc <- Xc[, good_cols, drop = FALSE]

datCSub <- data.frame(SOC = datC$`SOC`, Xc)
colnames(datCSub) <- c("SOC", paste0("spec.", colnames(Xc)))

datCSub <- datCSub[complete.cases(datCSub), ]

# prepare the data, the column name cannot be numeric, add 'spec.' in front
datCSub <- data.frame(SOC = datC$`SOC`, Xc)
colnames(datCSub) <- c("SOC", paste0("spec.", colnames(Xc)))

soilCRFModel <- train(
  SOC ~ .,
  data = datCSub,
  method = "rf",
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 1)
)

# summary of the model
soilCRFModel

# Calculate variable importance
rfImportance <- varImp(soilCRFModel)

# Print importance values
print(rfImportance)

# Plot the 30 most important wavelengths
plot(rfImportance, top = 30)

# predict on the calibration data
soilCRFPred <- predict(soilCRFModel, datCSub)

# prepare the validation data
Xv <- as.matrix(datV$spcDerive1)
Xv <- Xv[, good_cols, drop = FALSE]

datVSub <- data.frame(SOC = datV$`SOC`, Xv)
colnames(datVSub) <- c("SOC", paste0("spec.", colnames(Xv)))

# predict on the validation data
soilVRFPred <- predict(soilCRFModel, datVSub)

par(mfrow = c(1, 2))

# plot calibration
plot(datC$`SOC`, soilCRFPred,
     main = list("Calibration", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 1,
     pch = 16)
abline(0, 1, col = "red")

# plot validation
plot(datV$`SOC`, soilVRFPred,
     main = list("Validation", cex = 1.5, col = "black", font = 1.5),
     xlab = list("Observed", cex = 1, col = "black", font = 2),
     ylab = list("Predicted", cex = 1, col = "black", font = 2),
     xlim = c(0, 4),
     ylim = c(0, 4),
     col = rgb(red = 0.1, green = 0.1, blue = 0.8, alpha = 0.8),
     cex = 1,
     pch = 16)
abline(0, 1, col = "red")

# accuracy measures for calibration
soilspec::eval(datC$`SOC`, soilCRFPred, obj = "quant")

# accuracy measures for validation
soilspec::eval(datV$`SOC`, soilVRFPred, obj = "quant")

######################## Support Vector Regression ########################

datCSub <- data.frame(SOC = datC$SOC, datC$spcDerive1)
colnames(datCSub) <- c("SOC", paste0("X.", colnames(datC$spcDerive1)))


datVSub <- data.frame(SOC = datV$SOC, datV$spcDerive1)
colnames(datVSub) <- c("SOC", paste0("X.", colnames(datV$spcDerive1)))

# Define training control
# set.seed(12) 
train.control <- trainControl(method = "repeatedcv", number = 5, repeats = 1) 

# hyperparameters
svmGrid <- expand.grid(sigma= 10^(-5:-1), C= 10^(-3:1))

# Train the model
model <- train(SOC~., data = datCSub, method = "svmRadial",
               tuneGrid = svmGrid,
               trControl = train.control)

# Summarize the results
print(model)

# plot calibration
plot(datCSub$SOC, predict(model, datCSub),
     xlab = "Observed",
     ylab = "Predicted",
     xlim = c(0, 5),
     ylim = c(0, 5),
     pch = 16)
abline(0, 1)

# plot validation
plot(datVSub$SOC, predict(model, datVSub),
     xlab = "Observed",
     ylab = "Predicted",
     xlim = c(0, 4),
     ylim = c(0, 4),
     pch = 16)
abline(0, 1)

# accuracy measures for calibration
soilspec::eval(datCSub$SOC, predict(model, datCSub), obj = "quant")

# accuracy measures for validation
soilspec::eval(datVSub$SOC, predict(model, datVSub), obj = "quant")

########################## Memory Based Learning (MBL) ####################

# Define tested neighbourhood sizes
k2t <- seq(from = 20, to = 40, by = 1)


mblResults1 <- mbl(
  Xr = datC$spcDerive1,
  Yr = datC$SOC,
  Xu = datV$spcDerive1,
  Yu = NULL,
  
  # PCA-based dissimilarity
  diss_method = diss_pca(
    center = TRUE,
    scale = FALSE
  ),
  
  # nearest neighbour validation
  control = mbl_control(
    validation_type = "NNv"
  ),
  
  # test different numbers of neighbours
  neighbors = neighbors_k(
    k = k2t
  ),
  
  # weighted adaptive PLS
  fit_method = fit_wapls(
    min_ncomp = 4,
    max_ncomp = 17
  )
)


# Display model summary
mblResults1



# Select optimal neighbours

# Extract validation results
val_results <- mblResults1$validation_results$nearest_neighbor_validation


# RMSE values
rmseMBL <- val_results$rmse


# Tested neighbour numbers
neighNumber <- val_results$k


# Plot RMSE
plot(
  neighNumber,
  rmseMBL,
  type = "b",
  pch = 16,
  xlab = "Number of neighbours",
  ylab = "RMSE"
)


# Select optimal k
optNn <- neighNumber[which.min(rmseMBL)]

cat("Optimal number of neighbours =", optNn, "\n")



# Final MBL validation model

mblResults1Val <- mbl(
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



# Extract predictions

# The prediction is stored as results$k_optNn$pred

predMBL <- mblResults1Val$results[[paste0("k_", optNn)]]$pred


# Check predictions
str(predMBL)

# Check length agreement
length(predMBL)
length(datV$SOC)

# Validation plot

plot(
  datV$SOC,
  predMBL,
  main = "MBL validation",
  xlab = "Observed SOC",
  ylab = "Predicted SOC",
  xlim = range(datV$SOC),
  ylim = range(predMBL),
  pch = 16
)

abline(
  0,
  1,
  col = "red"
)

# Accuracy statistics

soilspec::eval(
  obs = datV$SOC,
  pred = predMBL,
  obj = "quant"
)
