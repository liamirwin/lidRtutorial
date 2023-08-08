# ======================================
#    SELECTION OF REGIONS OF INTEREST
# ======================================

# https://r-lidar.github.io/lidRbook/engine.html#engine-clip

# This code demonstrates the selection of regions of interest from LiDAR data.
# Simple geometries like circles and rectangles are selected based on coordinates.
# Complex geometries are extracted from shapefiles to clip specific areas.

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)
library(sf)

# 1. Select simple geometries
# =============================

# Load pointcloud
las = readLAS(files = "data/MixedEucaNat_normalized.laz",  filter = "-set_withheld_flag 0")

# Inspect
las
las@header$`Number of point records`

# Check documentation
?clip_roi

# Establish coordinates
x <- 203890
y <- 7358935

# Circle based on coordinates
circle = clip_circle(las = las, xcenter = x, ycenter = y, radius = 30)

# Inspect
circle
circle@header$`Number of point records`
plot(circle)

# Rectangle based on coordinates
rect = clip_rectangle(las = las, xleft = x, ybottom = y, xright = x + 40, ytop = y + 30)
plot(rect)

# Multiple coordinates
x = runif(2, x, x)
y = runif(2, 7358900, 7359050)

plots = clip_circle(las = las, xcenter = x, ycenter = y, radius = 10)
plots

plot(plots[[1]])
plot(plots[[2]])


# 2. Extraction of complex geometries from shapefiles
# ===================================================

planting = shapefile("data/shapefiles/MixedEucaNat.shp")

plot(las@header, map = FALSE)
plot(planting, add = TRUE, col = "#08B5FF39")

eucalyptus = clip_roi(las = las, geometry = planting)


# for sf users
planting = st_read(dsn = "data/shapefiles/MixedEucaNat.shp", quiet = T)

plot(las@header, map = FALSE)
plot(planting, add = TRUE, col = "#08B5FF39")

eucalyptus = clip_roi(las, planting)

plot(eucalyptus)

# E. Exercises and questions
# ====================================

# Using:
plots = st_read("data/shapefiles/MixedEucaNatPlot.shp")
plot(las@header, map = FALSE)
plot(plots, add = TRUE)

# 1. Clip the 5 plots with a radius of 11.3 m.

# 2. Clip a transect from A c(203850, 7358950) to B c(203950, 7959000).

# 3. Clip a transect from A c(203850, 7358950) to B c(203950, 7959000) but reorient it so
# it is no longer on the XY diagonal. Hint = ?clip_transect
