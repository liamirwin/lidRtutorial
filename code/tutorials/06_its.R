# Clear environment and specific warnings
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)
library(terra)

# Read in LiDAR file and set some color palettes
las <- readLAS("data/MixedEucaNat_normalized.laz",  filter = "-set_withheld_flag 0")
col1 <- height.colors(50)
col2 <- pastel.colors(900)

# Generate CHM
chm <- grid_canopy(las = las, res = 0.5, algorithm = p2r(0.15))
plot(chm, col = col1)

# Generate kernel and smooth chm
kernel <- matrix(1, 3, 3)
schm <- terra::focal(chm, w = kernel, fun = median, na.rm = TRUE)
plot(schm, col = height.colors(30))

# Detect trees
ttops <- locate_trees(las = schm, algorithm = lmf(ws = 2.5))
ttops
plot(chm, col = col1)
plot(ttops, col = "black", add = TRUE, cex = 0.5)

# Segment trees using dalponte
las <- segment_trees(las = las, algorithm = dalponte2016(chm = schm, treetops = ttops))

# Count number of trees detected and segmented
length(unique(las@data$treeID) |> na.omit())

# Visualize using intensity values as colors
plot(las, color = "treeID", bg = "white")

# Select trees by ID
tree25 <- filter_poi(las = las, treeID == 25)
tree125 <- filter_poi(las = las, treeID == 125)

# Visualize using intensity values as colors
plot(tree25, size = 4, bg = "white")

# Visualize using intensity values as colors
plot(tree125, size = 4, bg = "white")

# Generate rasterized delineation
trees <- dalponte2016(chm = chm, treetops = ttops)() # Notice the parenthesis at the end
trees

plot(trees, col = col2)
plot(ttops, add = TRUE, cex = 0.5)

# Detect trees
ttops <- locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))

x <- plot(las, bg = "white")
add_treetops3d(x = x, ttops = ttops, radius = 0.5)

# Segment using li
las <- segment_trees(las = las, algorithm = li2012())

plot(las, color = "treeID", bg = "white")

# Generate metrics for each delineated crown
metrics <- crown_metrics(las = las, func = ~list(n = length(Z)))
metrics
plot(metrics["n"], cex = 0.8)

# User defined function for area calculation
f <- function(x, y) {
  # Get xy for tree
  coords <- cbind(x, y)
  
  # Convex hull
  ch <- chull(coords)
  
  # Close coordinates
  ch <- c(ch, ch[1])
  ch_coords <- coords[ch, ]
  
  # Generate polygon
  p <- st_polygon(list(ch_coords))
  
  #calculate area
  area <- st_area(p)
  
  return(list(A = area))
}

# Apply user-defined function
metrics <- crown_metrics(las = las, func = ~f(X, Y))
metrics
plot(metrics["A"], cex = 0.8)

metrics <- crown_metrics(las = las, func = .stdtreemetrics)
metrics

# Visualize individual metrics
plot(x = metrics["convhull_area"], cex = 0.8)
plot(x = metrics["Z"], cex = 0.8)

cvx_hulls <- delineate_crowns(las = las, func = .stdtreemetrics)
cvx_hulls

plot(cvx_hulls)
plot(ttops, add = TRUE, cex = 0.5)

# Visualize individual metrics based on values
plot(x = cvx_hulls["convhull_area"])
plot(x = cvx_hulls["Z"])

##############################
##  Exercises and Questions ##
##############################

# Using:

las <- readLAS(files = "data/example_corrupted.laz", select = "xyz")

#### E1.

# Run `las_check()` and fix the errors.

#### E2.

# Find the trees and count the trees.

#### E3.

# Compute and map the density of trees with a 10 m resolution.

#### E4.

# Segment the trees.

#### E5.

# Assuming that the biomass of a tree can be estimated using the crown area
# and the mean Z of the points with the formula `2.5 * area + 3 * mean Z`, estimate the biomass of each tree.

#### E6.

# Map the total biomass at a resolution of 10 m. The output is a mix of ABA and ITS.
# Hint: use the `terra` package to rasterize spatial object with the function `rasterize()`.
