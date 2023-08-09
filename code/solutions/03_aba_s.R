##############################
##    Exercise Solutions    ##
##############################

library(lidR)

############
## 03_aba ##
############

las = readLAS("data/MixedEucaNat_normalized.laz", select = "*",  filter = "-set_withheld_flag 0")

# 1. Assuming that biomass is estimated using the equation <B = 0.5 * mean Z + 0.9 * 90th percentile of Z>
#    applied on first returns only, map the biomass.

B = grid_metrics(las, ~0.5*mean(Z) + 0.9*quantile(Z, probs = 0.9), 10, filter = ~ReturnNumber == 1L)
plot(B, col = height.colors(50))

B = grid_metrics(las, .stdmetrics_z, 10)
B = 0.5*B[["zmean"]] + 0.9*B[["zq90"]]
plot(B, col = height.colors(50))

grid_metrics(las, ~as.list(quantile(Z), 10))

# 2. Map the density of ground returns at a 5 m resolution with pixel_metrics(filter = ~Classification == LASGROUND).

GND = grid_metrics(las, ~length(Z)/25, res = 5, filter = ~Classification == LASGROUND)
plot(GND, col = heat.colors(50))

# 3. Map pixels that are flat (planar) using 'stdshapemetrics'. These could indicate potential roads.

m = grid_metrics(las, .stdshapemetrics, res = 3)
plot(m[["planarity"]], col = heat.colors(50))
flat = m[["planarity"]] > 0.85
plot(flat)
