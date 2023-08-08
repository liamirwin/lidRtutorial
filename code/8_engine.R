# ======================================
# LASCATALOG PROCESSING ENGINE
# ======================================

# https://r-lidar.github.io/lidRbook/engine.html
# https://r-lidar.github.io/lidRbook/engine2.html
# https://r-lidar.github.io/lidRbook/outbox.html

# This code showcases the LASCATALOG PROCESSING ENGINE, which efficiently and in
# parallel applies various functions to LiDAR catalogs. It introduces the
# catalog_apply() function for processing LiDAR data in a catalog. The code includes
# routines to detect trees and calculate metrics on the LiDAR catalog.

# Clear environment and specific warnings
options("rgdal_show_exportToProj4_warnings"="none")
rm(list = ls(globalenv()))

# Load libraries
library(lidR)
library(terra)
library(future)

# A. Example of catalog_apply()
# =====================================

# https://r-lidar.github.io/lidRbook/engine2.html#engine-catalog-apply

# Problem: How can we apply the following on a catalog?
# -----------------------------------------------------

ctg = readLAScatalog("data/Farm_A/")
las = readLAS(ctg$filename[16], filter = "-drop_withheld -drop_z_below 0 -drop_z_above 40")

surflas = filter_surfacepoints(las, 1)

plot(las)
plot(surflas)

# Calculate rumple index
?rumple_index

ri = pixel_metrics(las, ~rumple_index(X,Y,Z), 10)
plot(ri)

# Solution: LAScatalog processing engine
# --------------------------------------

# 1. Basic Usage
#===============

# Basic function
routine <- function(chunk){
  las <- readLAS(chunk)
}

# Error
catalog_apply(ctg, routine)

# Need to specify empty chunks as NULL

routine <- function(chunk){

  # read in chunk and check it is non-null
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)

  # Do some computation
  m <- pixel_metrics(las = las, func = ~max(Z),res = 20)
  output <- terra::crop(m, terra::ext(chunk))

  # output is a list
  return(output)
}

# Initialize parallel
plan(multisession)

# Apply the routine to catalog
out <- catalog_apply(ctg, routine)

# Inspect output list
out

# Use engine supported method for merging
options <- list(automerge = TRUE)
out <- catalog_apply(ctg, routine, .options = options)
print(out)


# 2. User defined functions
# =========================

# Calculating rumple index on a catalog

routine_rumple = function(chunk, res1 = 10, res2 = 1)
{
  # Read in check, check NULL status, get bbox
  las <-  readLAS(chunk)
  if (is.empty(las)) return(NULL)
  bbox <- terra::ext(chunk)

  # Calculate rumple index on surface points
  las <- filter_surfacepoints(las, res2)
  ri  <- pixel_metrics(las, ~rumple_index(X,Y,Z), res1)

  # crop to chunk extent and return output
  output <- terra::crop(ri, bbox)
  return(output)
}


opt_select(ctg) <- "xyz"
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- 0

options = list(automerge = TRUE, alignment = 10)
ri = catalog_apply(ctg, routine_rumple, res1 = 10, res2 = 0.5, .options = options)

plot(ri, col = height.colors(50))

# B. Example 2 of catalog_apply()
# =====================================

# Detecting trees and calculating metrics on a catalog

routine_trees = function(chunk) {

  # Read in check, check NULL status, get bbox
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  bbox <- st_bbox(chunk)

  # Filter surface points and generate chm
  las <- filter_surfacepoints(las, res = 0.5)
  chm <- rasterize_canopy(las = las, res = 0.5, algorithm = p2r())

  # Tree detection, segmentation, metrics
  ttops <- locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))
  las_trees <- segment_trees(las = las, algorithm = dalponte2016(chm = chm, treetops = ttops))
  p <- crown_metrics(las = las_trees, func = .stdtreemetrics)
  p <- sf::st_crop(x = p, y = bbox)

  m <- delineate_crowns(las_trees)
  output <- m[m$treeID %in% p$treeID,]

  # Remove buffer

  return(output)

}

# Set options
opt_chunk_buffer(ctg) <- 15
options = list(automerge = TRUE) # merge all outputs

# Apply function to catalog
m = catalog_apply(ctg, routine_trees, .options = options)

# Warnings - safe to ignore
warnings()

# Inspect and visualize
m
plot(m)

# End parallel
future::plan(sequential)


# C. Exercises
# ====================================

# 1. In example 2 (section B) what does last line `m <- m[m$treeID %in% p$treeID,]` do?
#    Adjust the function to not include that line to see what happens (use catalog_select() to select 4 tiles to test on)

subctg = catalog_select(ctg)

# 2. The following is a simple (and a bit naive) function to remove high noise points.
#    - Explain what this function does
#    - Create a user-defined function to apply using catalog_apply()
#    - Hint: Dont forget about buffered points... remember lidR::filter_* functions.

filter_noise = function(las, sensitivity)
{
  p95 <- grid_metrics(las, ~quantile(Z, probs = 0.95), 10)
  las <- merge_spatial(las, p95, "p95")
  las <- filter_poi(las, Z < 1+p95*sensitivity, Z > -0.5)
  las$p95 <- NULL
  return(las)
}

las = readLAS("data/Farm_A/PRJ_A_207480_7357420_g_c_d_n_u.laz")
nonoise = filter_noise(las, 1.2)

# 3. Design an application that retrieves the convex hull of each flightline (hard)
#    Use the concaveman::concaveman() function, adn functions from sf.
#    Start by designing a test function that works on a LAS object and
#    later apply on the collection. The output should look like:

flightlines = st_read("data/flightlines.shp")
plot(flightlines, col = sf.colors(6, alpha = 0.5))
plot(flightlines[3,])



