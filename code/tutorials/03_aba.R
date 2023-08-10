# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

# Load LiDAR data, excluding withheld flag points
las <- readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")

# Compute the mean height of points within 10x10 m pixels
hmean <- pixel_metrics(las = las, func = ~mean(Z), res = 10)
hmean
plot(hmean, col = height.colors(50))

# Compute the max height of points within 10x10 m pixels
hmax <- pixel_metrics(las = las, func = ~max(Z), res = 10)
hmax
plot(hmax, col = height.colors(50))

# Compute several metrics at once using a list
metrics <- pixel_metrics(las = las, func = ~list(hmax = max(Z), hmean = mean(Z)), res = 10)
metrics
plot(metrics, col = height.colors(50))

# Simplify computing metrics with predefined sets of metrics
metrics <- pixel_metrics(las = las, func = .stdmetrics_z, res = 10)
metrics
plot(metrics, col = height.colors(50))

# Plot a specific metric from the predefined set
plot(metrics, "zsd", col = height.colors(50))

# Generate a user-defined function to compute weighted mean
f <- function(x, weight) { sum(x*weight)/sum(weight) }

# Compute grid metrics for the user-defined function
X <- pixel_metrics(las = las, func = ~f(Z, Intensity), res = 10)

# Visualize the output
plot(X, col = height.colors(50))

##############################
##  Exercises and Questions ##
##############################

#Using:

las <- readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")

#### E1.

# Assuming that biomass is estimated using the equation `B = 0.5 * mean Z + 0.9 * 90th percentile of Z`
# applied on first returns only, map the biomass.

#### E2.

# Map the density of ground returns at a 5 m resolution with `pixel_metrics(filter = ~Classification == LASGROUND)`.

#### E3.

# Map pixels that are flat (planar) using `stdshapemetrics`. These could indicate potential roads.
