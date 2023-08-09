# Clear environment and specific warnings
rm(list = ls(globalenv()))
options("rgdal_show_exportToProj4_warnings" = "none")

# Load libraries
library(lidR)
library(microbenchmark)
library(terra)

# Load LiDAR data and reduce point density
las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-keep_random_fraction 0.4 -set_withheld_flag 0")

# Visualize the LiDAR point cloud
plot(las)

# Generate the CHM using a simple point-to-raster based algorithm
chm <- grid_canopy(las = las, res = 2, algorithm = p2r())

# Visualize the CHM
plot(chm, col = height.colors(50))

# The above method is strictly equivalent to using pixel_metrics to compute max height
chm <- pixel_metrics(las = las, func = ~max(Z), res = 2)

# Visualize the CHM
plot(chm, col = height.colors(50))

# However, the grid_canopy algorithm is optimized
microbenchmark::microbenchmark(canopy = grid_canopy(las = las, res = 1, algorithm = p2r()),
                               metrics = pixel_metrics(las = las, func = ~max(Z), res = 1),
                               times = 10)

# Increasing the resolution results in fewer empty pixels
chm <- grid_canopy(las = las, res = 1, algorithm = p2r())
plot(chm, col = height.colors(50))

# Using the 'subcircle' option turns each point into a disc of 8 points with a radius r
chm <- grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

# Increasing the subcircle radius, but it may not have meaningful results
chm <- grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.8))
plot(chm, col = height.colors(50))

# We can fill empty pixels using TIN interpolation
chm <- grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15, na.fill = tin()))
plot(chm, col = height.colors(50))

# Triangulation of first returns to generate the CHM
chm <- grid_canopy(las = las, res = 1, algorithm = dsmtin())
plot(chm, col = height.colors(50))

# Increasing the resolution results in a more detailed CHM
chm <- grid_canopy(las = las, res = 0.5, algorithm = dsmtin())
plot(chm, col = height.colors(50))

# Using the Khosravipour et al. pit-free algorithm with specified thresholds and maximum edge length
thresholds <- c(0, 5, 10, 20, 25, 30)
max_edge <- c(0, 1.35)
chm <- grid_canopy(las = las, res = 0.5, algorithm = pitfree(thresholds, max_edge))
plot(chm, col = height.colors(50))

# Using the 'subcircle' option with the pit-free algorithm
chm <- grid_canopy(las = las, res = 0.5, algorithm = pitfree(thresholds, max_edge, 0.1))
plot(chm, col = height.colors(50))

# Post-process the CHM using the 'terra' package and focal() function for smoothing
ker <- matrix(1, 3, 3)
schm <- terra::focal(chm, w = ker, fun = mean, na.rm = TRUE)

# Visualize the smoothed CHM
plot(schm, col = height.colors(50))
