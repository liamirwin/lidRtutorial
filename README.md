## Materials

This repository contains the material for a \~3 hour `lidR` tutorial workshop. You should install the material on your own machine from [this repository](https://github.com/tgoodbody/lidRtutorial). It contains the code, the shapefiles and point-clouds we will use. The workshop intends to:

-   Present an overview of what can be done with `lidR`
-   Give users an understanding of how `lidR` may fit their needs

## Requirements

### R version and Rstudio

-   You need to install a recent version of `R` i.e. `R 4.0.x` or more.
-   We will work with [Rstudio](https://www.rstudio.com/). This IDE is not mandatory to follow the workshop but is highly recommended.

### R Packages

You need to install the `lidR` package in its latest version (v \>= 3.1.3).

``` r
install.packages("lidR")
```

To run all code in the tutorial yourself, you will need to install the following libraries. You can use `lidR` without them, however.

``` r
libs <- c("geometry","viridis","future","sf","maptools","terra","mapview","mapedit","concaveman")

install.packages(libs)
```

## Estimated schedule

-   Introduction and set-up (09:00)
-   Read LAS and LAZ files (09:15)
-   Spatial queries (09:35)
-   Area-Based Approach (09:45)
-   Canopy Height Model (10:00)
-   Digital Terrain Model (10:10)

--- Break until 10:30 ---

-   Individual tree segmentation (10:30)
-   File collection processing engine (basic) (11:00)
-   File collection processing engine (advanced) (11:30)

## Resources

We strongly recommend having the following resources available to you:

-   The [`lidR` official documentation](https://cran.r-project.org/web/packages/lidR/lidR.pdf)
-   The [lidRbook](https://r-lidar.github.io/lidRbook/) of tutorials

When working on exercises:

-   [Stack Exchange with the `lidR` tag](https://gis.stackexchange.com/questions/tagged/lidr)

## `lidR`

`lidR` is an R package to work with LiDAR data developed at Laval University (Québec). It was developed & continues to be maintained by [Jean-Romain Roussel](https://github.com/Jean-Romain) and was made possible between:

-   2015 and 2018 thanks to the financial support of the AWARE project NSERC CRDPJ 462973-14; grantee Prof. Nicholas C. Coops.

-   2018 and 2021 thanks to the financial support of the Ministère des Forêts, de la Faune et des Parcs of Québec.

The current release version of `lidR` can be found on [CRAN](https://cran.r-project.org/web/packages/lidR/) and source code is hosted on [GitHub](https://github.com/r-lidar/lidR).
