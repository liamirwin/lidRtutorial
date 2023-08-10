# Clear environment and specific warnings
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

ctg <- readLAScatalog(folder = "data/Farm_A/")

ctg

plot(ctg)

# Interactive
plot(ctg, map = TRUE)

# coordinate system
crs(ctg)
projection(ctg)
st_crs(ctg)

# spatial extents
extent(ctg)
bbox(ctg)
st_bbox(ctg)

las_check(las = ctg)

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))

# check if files have .lax
is.indexed(ctg)

# generate index files
lidR:::catalog_laxindex(ctg)

# check if files have .lax
is.indexed(ctg)

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))

chm <- rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

# Check for warnings
warnings()

# Setting options and re-rasterizing the CHM
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"
chm <- rasterize_canopy(las = ctg, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = height.colors(50))

opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"

model <- pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40 -keep_first"
model <- pixel_metrics(las = ctg, func = ~max(Z), res = 20)
plot(model, col = height.colors(50))

# Set coordinate groups
x <- c(207846, 208131, 208010, 207852, 207400)
y <- c(7357315, 7357537, 7357372, 7357548, 7357900)

# Visualize coordinate groups
plot(ctg)
points(x, y)

# Clip plots
rois <- clip_circle(las = ctg, xcenter = x, ycenter = y, radius = 30)

# Visualize the LiDAR data with a default color palette
plot(rois[[1]], bg = "white")

# Visualize the LiDAR data with a default color palette
plot(rois[[3]], bg = "white")

las_check(rois[[1]])
las_check(rois[[3]])

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- paste0(tempdir())
file_list <- list.files(path)
delete_tif <- file_list[grep("\\.tif$", file_list)]
delete_las <- file_list[grep("\\.laz$", file_list)]
file.remove(file.path(path, delete_tif))
file.remove(file.path(path, delete_las))

# Read single file as catalog
ctg_non_norm <- readLAScatalog(folder = "data/MixedEucaNat.laz")

# Set options for output files
opt_output_files(ctg_non_norm) <- paste0(tempdir(),"/{XCENTER}_{XCENTER}")

# Write file as .laz
opt_laz_compression(ctg_non_norm) <- TRUE

# Get random plot locations and clip
x <- runif(n = 4, min = ctg_non_norm$Min.X, max = ctg_non_norm$Max.X)
y <- runif(n = 4, min = ctg_non_norm$Min.Y, max = ctg_non_norm$Max.Y)
rois <- clip_circle(las = ctg_non_norm, xcenter = x, ycenter = y, radius = 10)

# Read catalog of plots
ctg_plots <- readLAScatalog(tempdir())

# Set independent files option
opt_independent_files(ctg_plots) <- TRUE
opt_output_files(ctg_plots) <- paste0(tempdir(),"/{XCENTER}_{XCENTER}")

# Generate plot-level terrain models
rasterize_terrain(ctg_plots, 1, tin())

# Check files
path <- paste0(tempdir())
file_list <- list.files(path, full.names = TRUE)
file <- file_list[grep("\\.tif$", file_list)][[1]]

# plot dtm
plot(terra::rast(file))

# Set catalog options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"

# Detect tree tops and plot
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = TRUE, cex = 0.1, col = "black")

# Specify more options
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

# Detect treetops and plot
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))
plot(chm, col = height.colors(50))
plot(ttops, add = TRUE, cex = 0.1, col = "black")

library(future)

# Specify options
opt_filter(ctg) <- "-drop_withheld  -drop_z_below 0 -drop_z_above 40"
opt_select(ctg) <- "xyz"
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 10

# Visualize and summarize the catalog chunks
plot(ctg, chunk = TRUE)
summary(ctg)

# Process on single core
future::plan(sequential)

# Detect trees
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

# Process multi-core
future::plan(multisession)

# Detect trees
ttops <- locate_trees(las = ctg, algorithm = lmf(ws = 3, hmin = 5))

## # Example of network processing
## future::plan(remote, workers = c("localhost", "bob@132.203.41.87", "alice@132.203.41.87"))
## ttops <- locate_trees(ctg, lmf(3, hmin = 5))

# Back to single core
future::plan(sequential)

##############################
##  Exercises and Questions ##
##############################

# This exercise is complex because it involves options not yet described.
# Be sure to use the [lidRbook](https://r-lidar.github.io/lidRbook/) and 
# [package documentation](https://cran.r-project.org/web/packages/lidR/lidR.pdf).

# Using:
  
ctg <- readLAScatalog(folder = "data/Farm_A/")

#### E1.

# Generate a raster of point density for the provided catalog.
# Hint: Look through the documentation for a function that will do this!
  
#### E2.
  
# Modify the catalog to have a point density of 10 pts/m2 using the `decimate_points()` function.
# If you get an error make sure to read the documentation for `decimate_points()` and try:
# using `opt_output_file()` to write files to a temporary directory.

#### E3.

# Generate a raster of point density for this new decimated dataset.

#### E4.

# Read the whole decimated catalog as a single `las` file. The catalog isn't very big - not recommended for larger datasets!

#### E5.

# Read documentation for the `catalog_retile()` function and merge the decimated catalog into larger tiles.
