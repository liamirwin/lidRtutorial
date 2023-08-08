# ========================
#    CANOPY HEIGHT MODEL
# ========================

# https://r-lidar.github.io/lidRbook/chm.html

# This code demonstrates the creation of a Canopy Height Model (CHM) using LiDAR data
# It shows different algorithms for generating CHMs and provides options for adjusting
# resolution, subcircle size, and filling empty pixels.

options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# libraries if not already
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  install.packages("microbenchmark")
}
library(lidR)
library(microbenchmark)

# TODO ADD Surface Model?

# The original dataset has a too high point density and provides too good output for the purpose
# of this example. I use the filter -keep_random_fraction to purposely degrade the output

las = readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-keep_random_fraction 0.4 -set_withheld_flag 0")
col = height.colors(50)

las # ~10 pts/m2 instead of 24
plot(las)

# 1. Point-to-raster based algorithm
# ==================================

# Simple method that assigns the elevation of the highest point to each pixel
chm = grid_canopy(las = las, res = 2, algorithm = p2r())
plot(chm, col = col)

# Is strictly equivalent to
chm = pixel_metrics(las = las, func = ~max(Z), res = 2)
plot(chm, col = col)

# But optimized
microbenchmark::microbenchmark(canopy = grid_canopy(las = las, res = 1, algorithm = p2r()),
                               metrics = pixel_metrics(las = las, func = ~max(Z), res = 1),
                               times = 10)

# Better resolution implies few empty pixels
chm = grid_canopy(las = las, res = 1, algorithm = p2r())
plot(chm, col = col)

chm = grid_canopy(las = las, res = 0.5, algorithm = p2r())
plot(chm, col = col)

# The option 'subcircle' turns each point into a disc of 8 points with a radius r
chm = grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = col)

# We can increase the radius but it has not necessarily any meaning
chm = grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.8))
plot(chm, col = col)

# We can fill empty pixels
chm = grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15, na.fill = tin()))
plot(chm, col = col)

# 2. Triangulation based algorithm
# ================================

# Triangulation of first returns
chm = grid_canopy(las = las, res = 1, algorithm = dsmtin())
plot(chm, col = col)

chm = grid_canopy(las = las, res = 0.5, algorithm = dsmtin())
plot(chm, col = col)

# Khosravipour et al. pitfree algorithm
thresholds = c(0,5,10,20,25, 30)
max_edge = c(0, 1.35)
chm = grid_canopy(las = las, res = 0.5, algorithm = pitfree(thresholds, max_edge))
plot(chm, col = col)

# Option 'subcircle'
chm = grid_canopy(las = las, res = 0.5, algorithm = pitfree(thresholds, max_edge, 0.1))
plot(chm, col = col)

# C. Post process
# ================

# Usually the CHM can be post-processed. Often post-processing consists in smoothing.
# lidR has no tool for that. Indeed lidR is point-cloud oriented. Once you have
# a raster, it is the user responsibility to manipulate this kind of data. The users are free
# to do whatever they want within R or within external software such as GIS tools.
#
# Here we can use the 'raster' package and the focal() function

ker <- matrix(1,3,3)
schm <- focal(chm, w = ker, fun = mean, na.rm = TRUE)
plot(schm, col = col)

# D. Questions
# ================

# No exercise. A more comprehensive exercise in next section

