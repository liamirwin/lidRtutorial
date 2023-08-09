##############################
##    Exercise Solutions    ##
##############################

library(lidR)
library(sf)
library(terra)

############
## 06_its ##
############

# Using:
las = readLAS("data/example_corrupted.laz", select = "xyz")
col1 = height.colors(50)

# 1. Run las_check() and fix the errors

las_check(las)

las = filter_duplicates(las = las)

las_check(las)

# 2. Find the trees and count the trees

ttops = locate_trees(las = las, algorithm = lmf(ws = 3, hmin = 5))
x = plot(las)
add_treetops3d(x = x, ttops = ttops)

# 3. Compute and map the density of trees with a 10 m resolution [1]

r = terra::rast(x = ttops)
terra::res(r) <- 10
r = terra::rasterize(x = ttops, y = r, "treeID", fun = 'count')
plot(r, col = viridis::viridis(20))

# 4. Segment the trees

chm = grid_canopy(las = las, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = col1)
ttops = locate_trees(las = chm, algorithm = lmf(ws = 2.5))
las = segment_trees(las = las, dalponte2016(chm = chm, treetops = ttops))

plot(las, color = "treeID")

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
plot(x = V["V"])

# 6. Map the total biomass at a resolution of 10 m. The output is a mixed of ABA and ITS

Vtot = rasterize(V, r, "V", fun = "sum")
plot(Vtot, col = viridis::viridis(20))

