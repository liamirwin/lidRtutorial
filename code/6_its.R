# ======================================
#      INDIVIDUAL TREE SEGMENTATION
# ======================================

# https://r-lidar.github.io/lidRbook/itd-its.html

# This code demonstrates Individual Tree Segmentation (ITS) using LiDAR data.
# It covers CHM-based and point-cloud-based methods for tree detection and segmentation.
# The code also shows how to extract metrics at the tree level and visualize them.

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)
library(sf)

las = readLAS("data/MixedEucaNat_normalized.laz",  filter = "-set_withheld_flag 0")
col1 = height.colors(50)
col2 = pastel.colors(900)

# A. CHM based methods
# ====================

# 1.Build a CHM
# --------------

chm = grid_canopy(las = as, res = 0.5, algorithm = p2r(0.15))
plot(chm, col = col1)

# 2. Optionally smooth the CHM
# ----------------------------

kernel = matrix(1,3,3)
schm = raster::focal(chm, w = kernel, fun = median, na.rm = TRUE)
plot(schm, col = height.colors(30))

# 3. Tree detection
# ------------------

ttops = locate_trees(las = schm, algorithm = lmf(ws = 2.5))
ttops

plot(chm, col = col1)
plot(ttops, col = "black", add = T, cex = 0.5)

# 4. Segmentation
# -----------------

# Segment trees
las = segment_trees(las = las, algorithm = dalponte2016(chm = schm, treetops = ttops))

# Count number of trees detected and segmented
length(unique(las@data$treeID) |> na.omit())

# Visualize all trees
plot(las, color = "treeID", colorPalette = col2)

# Select trees by ID
tree25 = filter_poi(las = las, treeID == 25)
tree125 = filter_poi(las = las, treeID == 125)

plot(tree25, size = 4)
plot(tree125, size = 3)

# 5. Working with raster
# --------------------------------

# As said previously lidR deals with point-clouds. Thus, we used segment_trees with a point-cloud
# But here, from the CHM, the point-cloud is not strictly required.

trees = dalponte2016(chm = chm, treetops = ttops)() # Notice the parenthesis at the end

# Inspect
trees

plot(trees, col = col2)
plot(ttops, add = TRUE, cex = 0.5)

# B. Point-cloud based methods (no CHM)
# =====================================

# 1. Tree detection
# ------------------

# Detect
ttops = locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))

# Visualize tops
x = plot(las)
add_treetops3d(x = x, ttops = ttops, radius = .5)

# 3. Tree segmentation
# --------------------

# Segment
las = segment_trees(las = las, algorithm = li2012())

plot(las, color = "treeID", colorPalette = col2)

# This algorithm does not seem pertinent for this dataset.

# C. Extraction of metrics
# ==========================

# Redo segmentation using Dalponte's method
ttops = locate_trees(las = schm, algorithm = lmf(ws = 2.5))
las = segment_trees(las, dalponte2016(chm = schm, treetops = ttops))

plot(las, color = "treeID", colorPalette = col2)

# 1. crown_metrics() works like grid_metrics()
# ----------------------------------------

metrics = crown_metrics(las = las, func = ~list(n = length(Z)))
metrics
spplot(obj = metrics, "n", cex = 0.8)

# 2. It maps any user's function at the tree level like grid_metrics()
# -------------------------------------------------------------

# calculate convex hull of each tree and calculate  its area
f = function(x, y)
{
  ch <- chull(x,y)
  ch <- c(ch, ch[1])
  coords <- data.frame(x = x[ch], y = y[ch])
  p  = sp::Polygon(coords)
  area = p@area

  return(list(A = area))
}

# Apply user defined function
metrics = crown_metrics(las = las, func = ~f(X,Y))
metrics
spplot(obj = metrics, "A", cex = 0.8)

# 3. Some metrics are already recorded
# ------------------------------------

# Calculate tree level metrics
metrics = crown_metrics(las = las, func = .stdtreemetrics)
metrics

# Visualize individual metrics
spplot(obj = metrics, zcol = "convhull_area", cex = 0.8)
spplot(obj = metrics, zcol = "Z")

# 4. delineate_crowns: the same but with hull
# -------------------------------------------

cvx_hulls = delineate_crowns(las = las, func = .stdtreemetrics)
cvx_hulls

plot(cvx_hulls)
plot(ttops, add = TRUE, cex = 0.5)

# Visualize individual metrics based on values
spplot(obj = cvx_hulls, zcol = "convhull_area")
spplot(obj = cvx_hulls, zcol = "Z")


# E. Exercises and questions
# ==========================

# Using:
las = readLAS(files = "data/example_corrupted.laz", select = "xyz")

# 1. Run las_check() and fix the errors

# 2. Find the trees and count the trees

# 3. Compute and map the density of trees with a 10 m resolution [1]

# 4. Segment the trees

# 5. Assuming that the biomass of a tree can be estimated using the crown area and the mean Z
#   of the points with the formula <2.5 * area + 3 * mean Z>. Estimate the biomass of each tree.

# 6. Map the total biomass at a resolution of 10 m. The output is a mixed of ABA and ITS [1]

#' [1] Hint: use the 'raster' package to rasterize spatial object with the function rasterize()
