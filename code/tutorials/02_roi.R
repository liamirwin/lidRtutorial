# Clear environment and specific warnings
rm(list = ls(globalenv()))
options("rgdal_show_exportToProj4_warnings"="none")

# Load libraries
library(lidR)
library(sf)

las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")

# Inspect the header and the number of point records
las@header
las@header$`Number of point records`

# Establish coordinates
x <- 203890
y <- 7358935

# Select a circular area
circle <- clip_circle(las = las, xcenter = x, ycenter = y, radius = 30)

# Inspect the circular area and the number of point records
circle
circle@header$`Number of point records`

# Visualize the LiDAR data with a default color palette
plot(circle, bg = "white")

# Select a rectangular area
rect <- clip_rectangle(las = las, xleft = x, ybottom = y, xright = x + 40, ytop = y + 30)

# Visualize the LiDAR data with a default color palette
plot(rect, bg = "white")

# Select multiple random circular areas
x <- runif(2, x, x)
y <- runif(2, 7358900, 7359050)

plots <- clip_circle(las = las, xcenter = x, ycenter = y, radius = 10)

# Visualize the LiDAR data with a default color palette
plot(plots[[1]], bg = "white")

# Visualize the LiDAR data with a default color palette
plot(plots[[2]], bg = "white")

# Load the shapefile using sf
planting <- sf::st_read(dsn = "data/shapefiles/MixedEucaNat.shp", quiet = TRUE)

# Plot the LiDAR header information without the map
plot(las@header, map = FALSE)

# Plot the planting areas on top of the LiDAR header plot
plot(planting, add = TRUE, col = "#08B5FF39")

# Extract points within the planting areas using clip_roi()
eucalyptus <- clip_roi(las, planting)

# Plot the extracted points within the planting areas
plot(eucalyptus)

# Read the shapefile "MixedEucaNatPlot.shp" using st_read()
plots <- sf::st_read("data/shapefiles/MixedEucaNatPlot.shp")

# Plot the LiDAR header information without the map
plot(las@header, map = FALSE)

# Plot the extracted points within the planting areas
plot(plots, add = TRUE)

##############################
##  Exercises and Questions ##
##############################

# Now, let's read a shapefile called `MixedEucaNatPlot.shp` using `sf::st_read()` and plot it on top of the LiDAR header plot.

# Read the shapefile "MixedEucaNatPlot.shp" using st_read()
plots <- sf::st_read("data/shapefiles/MixedEucaNatPlot.shp")

# Plot the LiDAR header information without the map
plot(las@header, map = FALSE)

# Plot the extracted points within the planting areas
plot(plots, add = TRUE)

#### E1.

# Clip the 5 plots with a radius of 11.3 m.

#### E2.

# Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)`.

#### E3.

# Clip a transect from A `c(203850, 7358950)` to B `c(203950, 7959000)` 
# but reorient it so it is no longer on the XY diagonal. Hint = `?clip_transect`
