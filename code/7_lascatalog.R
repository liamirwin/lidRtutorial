# ======================================
# LASCATALOG PROCESSING
# ======================================

# This code performs various operations on LiDAR data
# using LAScatalog functionality. We visualize and
# inspect the data, validate the files, clip the data
# based on specific coordinates, and generate a Canopy
# Height Model (CHM), compute Above Ground Biomass (ABA)
# output, detect treetops, specify processing options, and use
# parallel computing.

# https://r-lidar.github.io/lidRbook/engine.html

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)
library(sf)

# A. Basic usage
# =======================

# Read catalog from directory of files
ctg = readLAScatalog(folder = "data/Farm_A/")

# Inspect
ctg

# Visualize
plot(ctg)
plot(ctg, map = TRUE)

# Check coordinate system and extent info
crs(ctg)
projection(ctg)
st_crs(ctg)

extent(ctg)
bbox(ctg)
st_bbox(ctg)

# B. LAScatalog objects validation
# ==========================================

# Validate files in catalog
las_check(las = ctg)

# X. Indexing
# ==============================

# https://cran.r-project.org/web/packages/lidR/vignettes/lidR-computation-speed-LAScatalog.html

#####
# Delete .lax files if they already exist
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))
#####

# Indexing files
is.indexed(ctg)

# Recommended to use lasindex from LASTOOLS but an internal undocumented function does exist
# https://github.com/r-lidar/lidR/issues/250
lidR:::catalog_laxindex(ctg)

# Check indexation again
is.indexed(ctg)

# D. CHM
# ==========================================

# Generate CHM
chm = rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

# Many problems here:
# 1. Warnings? - The CHM doesnt look right
warnings()

# X. Catalog options
# ==================

# https://r-lidar.github.io/lidRbook/engine.html#engine-summary-sheet

# Options and re-tile a catalog
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"

chm = rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

# E. GENERATE ABA ON CATALOG
# ==========================================

# Set catalog options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"

# Generate ABA output and visualize
model = pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

# Repeat with different first returns only
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40 -keep_first"
model = pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

# C. Clip
# ==========================================

# Set coordinate groups
x = c(207846, 208131, 208010, 207852, 207400)
y = c(7357315, 7357537, 7357372, 7357548, 7357900)

# Visualize
plot(ctg)
points(x,y)

# Clip plots
rois = clip_circle(las = ctg, xcenter = x, ycenter = y, radius = 30)

# Visualize
plot(rois[[1]])
plot(rois[[3]])

# Validate
las_check(rois[[1]])
las_check(rois[[3]])

# X. Independent files (e.g. plots) as catalogs
# =============================================

# https://r-lidar.github.io/lidRbook/engine.html#engine-independent-files

# Read single file as catalog
ctg_non_norm <- readLAScatalog(folder = "data/MixedEucaNat.laz")

# Set options
opt_output_files(ctg_non_norm) <- paste0(tempdir(), "/test")

# Get some random plot locations and clip
x <- runif(n = 4, min = ctg_non_norm$Min.X, max = ctg_non_norm$Max.X)
y <- runif(n = 4, min = ctg_non_norm$Min.Y, max = ctg_non_norm$Max.Y)
rois <- clip_circle(las = ctg_non_norm, xcenter = x, ycenter = y, radius = 10)

# Set independent option
opt_independent_files(rois) <- TRUE

# generate plot level terrain models
dtm <- rasterize_terrain(rois, 1, tin())
plot(dtm[[1]])

#  F. ITD USING CATALOG
# ==========================================

# Set catalog options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"

# Detect treetops and visualize
ttops = locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = T, cex = 0.1, col = "black")

# Specify a few catalog options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

ttops = locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = T, cex = 0.1, col = "black")

# H. Parallel computing
# ==========================================

# Load future library
library(future)

# Specify a few catalog options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

# Visualize and
plot(ctg, chunk = TRUE)
summary(ctg)

# Single core
future::plan(sequential)

ttops = locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

# Parallel
future::plan(multisession)

ttops = locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

# NOT RUN - Ability to parallelize over a network
'future::plan(remote, workers = c("localhost", "bob@132.203.41.87", "alice@132.203.41.87"))'
'ttops = locate_trees(ctg, lmf(3, hmin = 5))'
#

# Revert to single core
future::plan(sequential)

#  H. Exercises and questions
# ==========================================

# This exercise is complex because it involves options not yet described. Be sure to use the
# lidRbook and package documentation.

# https://cran.r-project.org/web/packages/lidR/lidR.pdf
# https://r-lidar.github.io/lidRbook/index.html

#Using:
ctg = readLAScatalog(folder = "data/Farm_A/")

# 1. Generate a raster of point density for the provided catalog. Hint: Look through the documentation
# for a function that will do this!

# 2. Modify the catalog to have a point density of 10 pts/m2 using the decimate_points() function.
# If you get an error make sure to read the documentation for decimate_points() and try:
#    - using opt_output_file() to write files to a temporary directory.

# 3. Generate a raster of point density for this new decimated dataset.

# 4. Read the whole decimated catalog as a single las file.
# The catalog isnt very big - not recommended for larger datasets!

# 5. Read documentation for the catalog_retile() function and merge the decimated catalog
#    into larger tiles.

