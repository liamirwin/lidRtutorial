# TODO filter_poi example with buffer / classification etc.

# =========================================
#  Reading, Plotting, Querying & Validating
# =========================================

# https://r-lidar.github.io/lidRbook/index.html

# Welcome to this LiDAR data processing tutorial
# using R and the 'lidR' library! In this tutorial,
# you will learn how to read, visualize, query,
# and validate LiDAR data.

# We'll explore basic information like the header and
# tabular data and visualize the data using different
# color schemes based on attributes.

# We'll use the 'select' argument to load specific
# attributes and the 'filter' argument to load only
# points of interest or apply transformations
# on-the-fly.

# We'll validate the LiDAR data using
# the 'las_check' function on different LiDAR data
# files to ensure data integrity.

# Let's get started with processing LiDAR data
# efficiently using 'lidR' and R! Happy learning!

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)
library(sf)

# ==============
# A. Basic usage
# ==============

# 1. Load pointcloud - warning explained later
# ---------------------------------------------

las = readLAS(files = "data/MixedEucaNat_normalized.laz")

# Inspect
las
las@header
las@data

# Check file size
format(object.size(las), "Mb")

# Visualize
print(las)

# Coordinate reference system and projection
crs(x = las)
projection(x = las)
st_crs(x = las)

# Extent details
extent(x = las)
bbox(obj = las)
st_bbox(obj = las)

# 2. Visualize based on attributes and manipulate colours
# -------------------------------------------------------

plot(las)
plot(las, colorPalette = terrain.colors(50))
plot(las, color = "Intensity")
plot(las, color = "Classification") # might take a bit of time, be patient :)
plot(las, color = "ScanAngleRank", axis = TRUE, legend = T)

# ==================
# B. Optimized usage
# ==================

# 1. Use the argument 'select' to load only the attributes of interest
# ------------------------------------------------------------------
?readLAS

las = readLAS(files = "data/MixedEucaNat_normalized.laz", select = "xyz")

# Inspect
las
las@data

# Check memory size
format(object.size(las), "Mb")

# Visualize
plot(las)

# Why does this fail?
plot(las, color = "Intensity")

# 2. Use the argument 'filter' to load only points of interest
# ----------------------------------------------------------

las = readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-keep_first")
las

# Check memory size again
format(object.size(las), "Mb")

plot(las)

# 3. Use the argument 'filter' to apply a transformation to the points on-the-fly
# -----------------------------------------------------------

# Get all available options for "-filter" method
readLAS(filter = "-h")
rlas::read.las(transform = "-h")

# Load and visualize with applied filter
las = readLAS(files = "data/MixedEucaNat.laz", filter = "-keep_class 2")
plot(las)

# X. Filtering using filter_poi()
# =================================

# Filter by attributes of points
class_2 <- filter_poi(las = las, Classification == 2L)
plot(class_2)

# Combine queries - if files are buffered it will work on that attribute aswell!
first_ground <-  filter_poi(las = las, Classification == 1L & ReturnNumber == 1L)
plot(first_ground)

# 4. LAS objects validation
# =========================================

las = readLAS(files = "data/MixedEucaNat_normalized.laz")
las_check(las)

las = readLAS(files = "data/example_corrupted.laz")
plot(las)
las_check(las)

# A degenerated point in LiDAR data refers to a point
# with identical XY(Z) coordinates as another point.
# This means two or more points occupy exactly the same
# location in XY/3D space. Degenerated points can cause
# issues in tasks like creating a 2.5D digital terrain
# model, as they don't add new information and can lead
# to inconsistencies. Identifying and handling
# degenerated points appropriately is crucial for
# accurate and meaningful results.

las = readLAS(files = "data/exemple_rgb.las")
plot(las, color = "RGB")
las_check(las)

# Check your data! If you see issues from las_check() make sure to go back and make sure your processing
# steps are correct!

# =======================
# Exercises and questions
# =======================

# Using MixedEucaNat_normalized.laz

las = readLAS(files = "data/MixedEucaNat_normalized.laz")

# 1. What are withheld points? Where are they in our pointcloud?

# 2. Read the file dropping withheld points.

# 3. The withheld points seem to be legitimate points that we want to keep.
# Try to load the file including the withheld points but get rid of the warning
# (without using suppressWarnings()). Hint: Check available '-set_withheld' filters in readLAS(filter = "-h")

# 4. Load only the ground points and plot the point-cloud coloured by the return
# number of the point. Do it loading the strict minimal amount of memory (4.7
# Mb). Hint: use ?lidR::readLAS and see what 'select' options might help.








