# ======================================
#    SELECTION OF REGIONS OF INTEREST
# ======================================

plots = st_read("data/shapefiles/MixedEucaNatPlot.shp")
plot(las@header, map = FALSE)
plot(plots, add = TRUE)

# 1. clip the 5 plots with a radius of 11.3 m,

inventory = clip_roi(las, plots, radius = 11.3)
plot(inventory[[2]])

# 2. Clip a transect from A c(203850, 7358950) to B c(203950, 7959000).

tr <- clip_transect(las, c(203850, 7358950), c(203950, 7359000), width = 5)
plot(tr, axis = T)

# 3. Clip a transect from A c(203850, 7358950) to B c(203950, 7959000) but reorient it so
# it is no longer on the XY diagonal. Hint = ?clip_transect

ptr <- clip_transect(las, c(203850, 7358950), c(203950, 7359000), width = 5, xz = TRUE)
plot(tr, axis = T)
plot(ptr, axis = T)
plot(ptr$X, ptr$Z, cex = 0.25, pch = 19, asp = 1)
