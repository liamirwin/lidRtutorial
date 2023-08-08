# ======================================
#         DIGITAL TERRAIN MODEL
# ======================================

# https://r-lidar.github.io/lidRbook/dtm.html

# This code demonstrates the creation of a Digital Terrain Model (DTM) from LiDAR data.
# It covers two different algorithms for DTM generation, triangulation of ground points, and inverse-distance weighting.
# The code also showcases DTM-based normalization and point-based normalization with exercises for hands-on practice.

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)

# 1. DTM
# =======================

las = readLAS(files = "data/MixedEucaNat.laz", filter = "-set_withheld_flag 0")
col = gray(1:50/50)

plot(las)
plot(las, color = "Classification") # be patient :)

# Triangulate the ground points
dtm = grid_terrain(las = las, res = 1, algorithm = tin())

plot(dtm, col = col)
plot_dtm3d(dtm)

x = plot(las)
add_dtm3d(x, dtm)

# Inverse-distance weighting of the ground points
dtm = grid_terrain(las = las, res = 1, algorithm = knnidw())
plot(dtm, col = col)
plot_dtm3d(dtm)

# B. Normalization
# =======================

# https://r-lidar.github.io/lidRbook/norm.html

# DTM-based normalization
# ------------------------

nlas = normalize_height(las = las, dtm = dtm)
plot(nlas)

# Filter ground points
gnd = filter_ground(las = nlas)
plot(gnd)
hist(gnd$Z, breaks = seq(-1.5,1.5,0.05))

# Shortcut
nlas = las - dtm
plot(nlas)

# point-based normalization
# -------------------------

# use algorithm parameter instead of dtm
nlas = normalize_height(las = las, algorithm = tin())
plot(nlas)

gnd = filter_ground(nlas)
plot(gnd)
hist(gnd$Z, breaks = seq(-1.5,1.5,0.05))

# ==========================
# C. Exercises and questions
# ==========================

# 1. Plot and compare these two normalized point-clouds. Why do they look different? Fix that. Hint: filter.

las1 = readLAS("data/MixedEucaNat.laz", filter = "-set_withheld_flag 0")
nlas1 = normalize_height(las1, tin())
nlas2 = readLAS("data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")
plot(nlas1)
plot(nlas2)

# 2. Clip a plot somewhere in MixedEucaNat.laz (the non-normalized file).

# 3. Compute a DTM for this plot. Which method are you choosing and why?

# 4. Compute a DSM (digital surface model). Hint: Look back to how you made a CHM.

# 5. Normalize the plot.

# 6. Compute a CHM.

# 7. Estimate some metrics of interest in this plot with cloud_metrics().
