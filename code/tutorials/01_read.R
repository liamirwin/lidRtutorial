# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)
library(sf)

las <- readLAS(files = "data/MixedEucaNat_normalized.laz")

# Inspect header information
las@header

# Inspect attributes of the point cloud
las@data

# Check the file size of the loaded LiDAR data
format(object.size(las), "Mb")

# Visualize the LiDAR data with a default color palette
plot(las, bg = "white")

# Visualize using intensity values as colors
plot(las, color = "Intensity", bg = "white")

# Visualize using the classification attribute as colors
plot(las, color = "Classification", bg = "white")

# Visualize using the scan angle rank as colors with an axis and legend
plot(las, color = "ScanAngleRank", axis = TRUE, legend = TRUE, bg = "white")

# Load only the xyz coordinates (X, Y, Z) and ignore other attributes
las <- readLAS(files = "data/MixedEucaNat_normalized.laz", select = "xyz")

# Inspect the loaded attributes
las@data

# Check the memory size after loading only the selected attributes
format(object.size(las), "Mb")

# Load only the first return points
las <- readLAS(files = "data/MixedEucaNat_normalized.laz", filter = "-keep_first")

# Inspect the loaded points
las

# Check the memory size after loading only the filtered points
format(object.size(las), "Mb")

# Visualize the filtered points
plot(las, bg = "white")

# Load and visualize with an applied filter
las <- readLAS(files = "data/MixedEucaNat.laz", filter = "-keep_class 2")

plot(las, bg = "white")

# Filter points with Classification == 2
class_2 <- filter_poi(las = las, Classification == 2L)

# Combine queries to filter points with Classification == 1 and ReturnNumber == 1
first_ground <- filter_poi(las = las, Classification == 2L & ReturnNumber == 1L)

plot(class_2, bg = "white")

plot(first_ground, bg = "white")

# Load and validate LAS data
las <- readLAS(files = "data/MixedEucaNat_normalized.laz")
las_check(las)

# Visualize corrupted LAS data
las <- readLAS(files = "data/example_corrupted.laz")

plot(las, bg = "white")

# Validate corrupted LAS data
las_check(las)

##############################
##  Exercises and Questions ##
##############################

#Using:
  
las <- readLAS(files = "data/MixedEucaNat_normalized.laz")

#### E1.

# What are withheld points? Where are they in our point cloud?
  
#### E2.
  
# Read the file dropping withheld points.

#### E3.

# The withheld points seem to be legitimate points that we want to keep. 
# Try to load the file including the withheld points but get rid of the warning 
# (without using `suppressWarnings()`). Hint: Check available `-set_withheld` filters using `readLAS(filter = "-h")`.

#### E4.

# Load only the ground points and plot the point cloud colored by the return
# number of the point. Do it loading the strict minimal amount of memory (4.7 Mb).
# Hint: use `?lidR::readLAS` and see what `select` options might help.

## Conclusion

# This concludes our tutorial on the basic usage of the `lidR` package in R for
# processing and analyzing LiDAR data. We covered loading LiDAR data, inspecting
# and visualizing the data, selecting specific attributes, filtering points, and validating LAS objects for issues.
