# ======================================
#         AREA-BASED APPROACH
# ======================================

# https://r-lidar.github.io/lidRbook/modeling-aba.html

# This code demonstrates an area-based approach for LiDAR data.
# Basic usage involves computing mean and max height of points within 10x10 m pixels and visualizing the results.
# The code shows how to compute multiple metrics simultaneously and use predefined metric sets.
# Advanced usage introduces user-defined metrics for more specialized calculations.
# pixel_metrics() was previously named grid_metrics() - both still work!

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries

# 1. Basic usage
# ==============

las = readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")
col = height.colors(50)

# Mean height of points within 10x10 m pixels
hmean = pixel_metrics(las = las, func = ~mean(Z), res = 10)
hmean
plot(hmean, col = col)

# Max height of points within 10x10 m pixels
hmax = pixel_metrics(las = las, func = ~max(Z), res = 10)
hmax
plot(hmax, col = col)

# We can compute several metrics at once using a list
metrics = pixel_metrics(las = las, func =  ~list(hmax = max(Z), hmean = mean(Z)), res = 10)
metrics
plot(metrics, col = col)


# For simplicity lidR proposes some sets of metrics
metrics = pixel_metrics(las = las, func = .stdmetrics_z, res = 10)
metrics
plot(metrics, col = col)

plot(metrics, "zsd", col = col)

?.stdmetrics

# 3rd party packages exist with comprehensive metric sets - https://github.com/ptompalski/lidRmetrics

# 2. Advanced usage with user defined metrics
# ===========================================

# The strength of the function is the ability to map almost everything
# https://r-lidar.github.io/lidRbook/engine2.html#engine-user-function-generic

# Generate user-defined function
f = function(x, weight) { sum(x*weight)/sum(weight) }

# Derive grid metrics for function (f)
X = pixel_metrics(las = las, ~f(Z, Intensity), 10)

# Visualize output
plot(X, col = col)

# F. Exercises and questions
# ===========================

# Using:

las = readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")

# 1. Assuming that biomass is estimated using the equation <B = 0.5 * mean Z + 0.9 * 90th percentile of Z>
#    applied on first returns only, map the biomass.

# 2. Map the density of ground returns at a 5 m resolution with pixel_metrics(filter = ~Classification == LASGROUND).

# 3. Map pixels that are flat (planar) using 'stdshapemetrics'. These could indicate potential roads.



