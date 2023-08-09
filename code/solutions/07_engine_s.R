##############################
##    Exercise Solutions    ##
##############################

library(lidR)

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
