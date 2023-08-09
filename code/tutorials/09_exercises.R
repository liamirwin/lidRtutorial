##############################
##    Exercise Solutions    ##
##############################

############
## 06_its ##
############

# Using:
las = readLAS("data/example_corrupted.laz", select = "xyz")
col1 = height.colors(50)
col2 = pastel.colors(900)

# 1. Run las_check() and fix the errors

las_check(las)

las = filter_duplicates(las = las)

las_check(las)

# 2. Find the trees and count the trees

ttops = locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))
x = plot(las)
add_treetops3d(x = x, ttops = ttops)

# 3. Compute and map the density of trees with a 10 m resolution [1]

r = raster::raster(x = ttops)
res(r) = 10
r = raster::rasterize(x = ttops, y = r, "treeID", fun = 'count')
plot(ra, col = viridis::viridis(20))

# 4. Segment the trees

chm = grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = col1)
ttops = locate_trees(las = chm, algorithm = lmf(ws = 2.5))
las = segment_trees(las = las, dalponte2016(chm = chm, treetops = ttops))

plot(las, color = "treeID", pal = col2)

# 5. Assuming that a value of interest of a tree can be estimated using the crown area and the mean Z
#    of the points with the formula <2.5 * area + 3 * mean Z>. Estimate the value of interet of each tree

value_of_interest = function(x,y,z)
{
  m = stdtreemetrics(x,y,z)
  avgz = mean(z)
  v = 2.5*m$convhull_area + 3 * avgz
  return(list(V = v))
}

V = crown_metrics(las = las, func = ~value_of_interest(X,Y,Z))
spplot(obj = V, "V")

# 6. Map the total biomass at a resolution of 10 m. The output is a mixed of ABA and ITS

Vtot = rasterize(V, r, "V", fun = "sum")
plot(Vtot, col = viridis::viridis(20))


###############
## 07_engine ##
###############

#This exercise is complex because it involves options not yet described. Be sure to use the
# lidRbook and package documentation.

# https://cran.r-project.org/web/packages/lidR/lidR.pdf
# https://r-lidar.github.io/lidRbook/index.html

#Using:
ctg = readLAScatalog(folder = "data/Farm_A/")

# 1. Generate a raster of point density for the provided catalog. Hint: Look through the documentation
# for a function that will do this!

ctg = readLAScatalog("data/Farm_A/", filter = "-drop_withheld -drop_z_below 0 -drop_z_above 40")
D1 = rasterize_density(las = ctg, res = 4)
plot(D1, col = heat.colors(50))

# 2. Modify the catalog to have a point density of 10 pts/m2 using the decimate_points() function.
# If you get an error make sure to read the documentation for decimate_points() and try:
#    - using opt_output_file() to write files to a temporary directory.

# https://r-lidar.github.io/lidRbook/engine.html#engine-dtm-ondisk

newctg = decimate_points(las = ctg, algorithm = homogenize(density = 10, res = 5))
#>  Error: This function requires that the LAScatalog provides an output file template.

opt_filter(ctg) <- "-drop_withheld"
opt_output_files(ctg) <- paste0(tempdir(), "/{ORIGINALFILENAME}")
newctg = decimate_points(las = ctg, algorithm = homogenize(density = 10, res = 5))

# 3. Generate a raster of point density for this new decimated dataset.

opt_output_files(newctg) <- ""
D2 = grid_density(las = newctg, res = 4)
plot(D2, col = heat.colors(50))

# 4. Read the whole decimated catalog as a single las file.
# The catalog isn't very big - not recommended for larger data sets!

las = readLAS(newctg)
plot(las)

# 5. Read documentation for the catalog_retile() function and merge the dataset
#    into larger tiles. Use `ctg` metadata to align new chunks to the lower left corner
#    of the old ones. Hint: Visualize the chunks and use opt_chunk_* options.

opt_chunk_size(ctg) <- 280
opt_chunk_buffer(ctg) <- 0
opt_chunk_alignment(ctg) <- c(min(ctg$Min.X), min(ctg$Min.Y))
plot(ctg, chunk = T)

opt_output_files(ctg) <- "{tempdir()}/PRJ_A_{XLEFT}_{YBOTTOM}"
newctg = catalog_retile(ctg = ctg)
plot(newctg)


################
## 08_engine2 ##
################

# 1. In example 2 (section B) what does last line `m <- m[m$treeID %in% p$treeID,]` do?
#    Adjust the function to not include that line to see what happens (use catalog_select() to select 4 tiles to test on)

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

m <-  catalog_apply(subctg, routine_trees_test, .options = options)
plot(m, col = rgb(0,0,1,0.3))

ctg <-  readLAScatalog("data/Farm_A/")
opt_select(ctg) <- "xyz"
opt_filter(ctg) <- "-drop_withheld -drop_z_below 0 -drop_z_above 40"
opt_chunk_buffer(ctg) <- 15
opt_chunk_size(ctg) <- 0
subctg <-  catalog_select(ctg)
options <-  list(automerge = TRUE)
m <-  catalog_apply(subctg, my_process, .options = options)

plot(m, col = rgb(0,0,1,0.3))


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

#las = readLAS(output)
#plot(las)


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
