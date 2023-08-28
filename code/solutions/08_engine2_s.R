##############################
##    Exercise Solutions    ##
##############################

library(lidR)
library(sf)

################
## 08_engine2 ##
################

# 1. In example 2 (section B) what does last line `m <- m[m$treeID %in% p$treeID,]` do?
# Adjust the function to not include that line to see what happens (use `catalog_select()` to select 4 tiles to test on)

# Subset catalog
subctg = catalog_select(ctg)

# without line
routine_trees_test = function(chunk) {
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
  
  # Remove buffer
  p <- sf::st_crop(x = p, y = bbox)
  
  # Delineate crowns
  output <- delineate_crowns(las_trees)
  
  #output <- m[m$treeID %in% p$treeID,]
  
  return(output)
}

options <-  list(automerge = TRUE)
m <-  catalog_apply(subctg, routine_trees_test, .options = options)
plot(m, col = rgb(0,0,1,0.3))

ctg <-  readLAScatalog("data/Farm_A/")
opt_select(ctg) <- "xyz"
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_chunk_buffer(ctg) <- 15
opt_chunk_size(ctg) <- 0
subctg <-  catalog_select(ctg)
options <-  list(automerge = TRUE)
m <-  catalog_apply(subctg, routine_trees_test, .options = options)

plot(m, col = rgb(0,0,1,0.3))

# 2. The following is a simple (and a bit naive) function to remove high noise points.
#    - Explain what this function does
#    - Create a user-defined function to apply using catalog_apply()
#    - Hint: Dont forget about buffered points... remember lidR::filter_* functions.
filter_noise = function(las, sensitivity)
{
  p95 <- pixel_metrics(las, ~quantile(Z, probs = 0.95), 10)
  las <- merge_spatial(las, p95, "p95")
  las <- filter_poi(las, Z < 1+p95*sensitivity, Z > -0.5)
  las$p95 <- NULL
  return(las)
}

filter_noise_collection = function(cl, sensitivity)
{
  las <- readLAS(cl)
  if (is.empty(las)) return(NULL)
  las <- filter_noise(las, sensitivity)
  las <- filter_poi(las, buffer == 0L)
  return(las)
}

ctg = readLAScatalog("data/Farm_A/")
opt_select(ctg) <- "*"
opt_filter(ctg) <- "-drop_withheld -drop_"
opt_output_files(ctg) <- "{tempdir()}/*"
opt_chunk_buffer(ctg) <- 20
opt_chunk_size(ctg) <- 0

options = list(automerge = TRUE)
output = catalog_apply(ctg, filter_noise_collection, sensitivity = 1.2, .options = options)

las = readLAS(output)
plot(las)

# 3. Design an application that retrieves the convex hull of each flightline (hard)
#    Use the concaveman::concaveman() function, adn functions from sf.
#    Start by designing a test function that works on a LAS object and
#    later apply on the collection. The output should look like:

library(sf)

# Read the catalog
ctg = readLAScatalog("data/Farm_A/")

# Read a single file to perform tests
las = readLAS(ctg$filename[16], select = "xyzp", filter = "-drop_withheld -drop_z_below 0 -drop_z_above 40")

# Define a function capable of building the hull from the XY of a given PointSourceID
enveloppes = function(x,y, psi)
{
  hull = concaveman::concaveman(cbind(x,y), length_threshold = 10)
  hull = sf::st_polygon(list(hull))
  hull = sf::st_sfc(hull)
  hull = sf::st_simplify(hull, dTolerance = 1)
  hull = sf::st_sf(hull)
  hull$ID = psi[1]
  list(hull = list(hull = hull))
}

# Define a function that apply the previous function to each PointSourceID from a LAS object
flighline_polygons = function(las)
{
  u = las@data[ , enveloppes(X,Y, PointSourceID), by = PointSourceID]
  hulls = Reduce(rbind, u$hull)
  return(hulls)
}

# Test this function on a LAS
hulls = flighline_polygons(las)
plot(hulls, col = sf.colors(3, alpha = 0.5))


# It works so let make a function that works with a LAScatalog
flighline_polygons = function(las)
{
  if (is(las, "LAS"))  {
    u = las@data[ , enveloppes(X,Y, PointSourceID), by = PointSourceID]
    hulls = Reduce(rbind, u$hull)
    return(hulls)
  }
  
  if (is(las, "LAScluster")) {
    las <- readLAS(las)
    if (is.empty(las)) return(NULL)
    hulls <- flighline_polygons(las)
    return(hulls)
  }
  
  if (is(las, "LAScatalog")) {
    opt_select(las) <-  "xyzp"
    options <- list(
      need_output_file = FALSE,
      need_buffer = TRUE,
      automerge = TRUE)
    output <- catalog_apply(las, flighline_polygons, .options = options)
    hulls <- dplyr::summarise(dplyr::group_by(output, ID), ID = ID[1])
    return(hulls)
  }
  
  stop("Invalid input")
}

library(future)
future::plan(multisession)
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) = "-drop_withheld -drop_z_below 0 -drop_z_above 40"
flightlines = flighline_polygons(ctg)
plot(flightlines, col = sf.colors(6, alpha = 0.5))
