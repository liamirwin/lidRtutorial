# Clear environment
rm(list = ls(globalenv()))

# Load packages
library(lidR)

# Load LiDAR data and filter out non-ground points
las <- readLAS(files = "data/MixedEucaNat.laz", filter = "-set_withheld_flag 0")

# Visualize using the classification attribute as colors
plot(las, bg = "white")

# Visualize using the classification attribute as colors
plot(las, color = "Classification", bg = "white")

# Generate a DTM using the TIN (Triangulated Irregular Network) algorithm
dtm_tin <- grid_terrain(las = las, res = 1, algorithm = tin())

# Visualize the DTM in 3D
plot_dtm3d(dtm_tin, bg = "white")

# Visualize the LiDAR data with the overlaid DTM in 3D
x <- plot(las, bg = "white")
add_dtm3d(x, dtm_tin, bg = "white")

# Generate a DTM using the IDW (Inverse-Distance Weighting) algorithm
dtm_idw <- grid_terrain(las = las, res = 1, algorithm = knnidw())

# Visualize the IDW-based DTM in 3D
plot_dtm3d(dtm_idw, bg = "white")

# Normalize the LiDAR data using DTM-based normalization
nlas_dtm <- normalize_height(las = las, algorithm = dtm_tin)

# Visualize the normalized LiDAR data
plot(nlas_dtm, bg = "white")

# Filter the normalized data to retain only ground points
gnd_dtm <- filter_ground(las = nlas_dtm)

# Visualize the filtered ground points
plot(gnd_dtm, bg = "white")

# Plot the histogram of normalized ground points' height
hist(gnd_dtm$Z, breaks = seq(-1.5, 1.5, 0.05))

# Normalize the LiDAR data using DTM-based normalization with TIN algorithm
nlas_tin <- normalize_height(las = las, algorithm = tin())

# Visualize the normalized LiDAR data using the TIN algorithm
plot(nlas_tin, bg = "white")

# Filter the normalized data (TIN-based) to retain only ground points
gnd_tin <- filter_ground(las = nlas_tin)

# Visualize the filtered ground points after TIN-based normalization
plot(gnd_tin, bg = "white")

# Plot the histogram of normalized ground points' height after TIN-based normalization
hist(gnd_tin$Z, breaks = seq(-1.5, 1.5, 0.05))

##############################
##  Exercises and Questions ##
##############################

#### E1.

# Plot and compare these two normalized point-clouds. Why do they look different? Fix that. Hint: `filter`.

# Load and visualize nlas1 and nlas2
las1 = readLAS("data/MixedEucaNat.laz", filter = "-set_withheld_flag 0")
nlas1 = normalize_height(las1, tin())
nlas2 = readLAS("data/MixedEucaNat_normalized.laz", filter = "-set_withheld_flag 0")
plot(nlas1)
plot(nlas2)

#### E2.

# Clip a plot somewhere in `MixedEucaNat.laz` (the non-normalized file).

#### E3.

# Compute a DTM for this plot. Which method are you choosing and why?
  
#### E4.
  
# Compute a DSM (digital surface model). Hint: Look back to how you made a CHM.

#### E5.

# Normalize the plot.

#### E6.

# Compute a CHM.

#### E7.

# Compute some metrics of interest in this plot with `cloud_metrics()`.
