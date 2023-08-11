# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(terra)
library(future)

# Read a LAS catalog
ctg <- readLAScatalog(folder = "data/Farm_A/", filter = "-drop_withheld")

# Inspect the first LAS file in the catalog
las_file <- ctg$filename[1]
las <- readLAS(files = las_file)
las

# Visualize the LiDAR data in 3D
plot(las, bg = "white")

# Read a LAS file from the catalog and filter surface points
las_file <- ctg$filename[16]
las <- readLAS(files = las_file, filter = "-drop_z_below 0 -drop_z_above 40")
surflas <- filter_surfacepoints(las = las, res = 1)

# Visualize the LiDAR data with a default color palette
plot(las, bg = "white")

# Visualize the surface points using a default color palette
plot(surflas, bg = "white")

# Generate Area-based metrics
ri <- pixel_metrics(las = las, func = ~rumple_index(X,Y,Z), res = 10)
plot(ri)

# User-defined function for processing chunks
routine <- function(chunk){
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  
  # Perform computation
  m <- pixel_metrics(las = las, func = ~max(Z), res = 20)
  output <- terra::crop(x = m, terra::ext(chunk))
  
  return(output)
}

# Initialize parallel processing
plan(multisession)

# Specify catalog options
opt_filter(ctg) <- "-drop_withheld"

# Apply routine to catalog
out <- catalog_apply(ctg = ctg, FUN = routine)

# Inspect the output list
out[1:5]

# Use the engine-supported method for merging
options <- list(automerge = TRUE)
out <- catalog_apply(ctg = ctg, FUN = routine, .options = options)
print(out)

# User-defined function for rumple index calculation
routine_rumple <- function(chunk, res1 = 10, res2 = 1){
  las <-  readLAS(chunk)
  if (is.empty(las)) return(NULL)
  bbox <- terra::ext(x = chunk)
  
  las <- filter_surfacepoints(las, res2)
  ri  <- pixel_metrics(las = las, func = ~rumple_index(X,Y,Z), res = res1)
  
  output <- terra::crop(x = ri, y = bbox)
  return(output)
}

# Set catalog options
opt_select(ctg) <- "xyz"
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 40"
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 0

# Specify options for merging
options <- list(automerge = TRUE, alignment = 10)

# Apply the user-defined function to the catalog
ri <- catalog_apply(ctg = ctg, FUN = routine_rumple, res1 = 10, res2 = 0.5, .options = options)

# Plot the output
plot(ri, col = height.colors(50))

# User-defined routine for tree detection and metrics
routine_trees <- function(chunk) {
  # Read in the chunk and check for emptiness
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  
  # Get the chunk bounding box
  bbox <- sf::st_bbox(obj = chunk)

  # Filter surface points and create canopy height model (CHM)
  las <- filter_surfacepoints(las = las, res = 0.5)
  chm <- rasterize_canopy(las = las, res = 0.5, algorithm = p2r())

  # Detect and segment trees
  ttops <- locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))
  las_trees <- segment_trees(las = las, algorithm = dalponte2016(chm = chm, treetops = ttops))
  
  # Generate metrics for each tree
  p <- crown_metrics(las = las_trees, func = .stdtreemetrics)
  p <- sf::st_crop(x = p, y = bbox)

  # Delineate convex hulls
  m <- delineate_crowns(las_trees)
  output <- m[m$treeID %in% p$treeID,]

  return(output)
}

# Set options for the catalog
opt_chunk_buffer(ctg) <- 15
options <- list(automerge = TRUE) # Merge all outputs

# Apply the function to the catalog
m <- catalog_apply(ctg = ctg, FUN = routine_trees, .options = options)

# View and visualize the output
m
plot(m)

# End parallel processing
future::plan(sequential)

# Instructions for cleaning up any existing .lax files
# (Note: Please replace 'path' with the appropriate path)
path <- "data/Farm_A/"
file_list <- list.files(path)
delete_lax <- file_list[grep("\\.lax$", file_list)]
file.remove(file.path(path, delete_lax))

##############################
##  Exercises and Questions ##
##############################

#### E1.

#Understanding Chunk Filtering

# On line 111 (routine_trees function) `m <- m[m$treeID %in% p$treeID,]` is used. 
# Explain the purpose of this line. To understand its impact, modify the function to
# exclude this line and observe the results. You can use the `catalog_select()` 
# function to choose a subset of tiles for testing.

subctg <- catalog_select(ctg)

#### E2.

# Implement Noise Filtering
# -   Explain the purpose of the `filter_noise()` function.
# -   Create a user-defined function to apply noise filtering using the `catalog_apply()` function.
# -   Make sure to consider buffered points when using lidR's `filter_*` functions.
 
#### E3.

# Flightline Convex Hull Application

# Design an application to retrieve the convex hull of each flightline using the `concaveman::concaveman()` function
# and functions from the `sf` package.

# -   Begin by designing a test function that works on a single LAS object.
# -   Apply the function to a collection of LAS files.
# -   Visualize the results using the flightlines' shapefile.

flightlines <- st_read("data/flightlines.shp")
plot(flightlines, col = sf.colors(6, alpha = 0.5))
plot(flightlines[3,])
