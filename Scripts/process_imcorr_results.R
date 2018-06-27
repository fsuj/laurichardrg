
# Scritp to read RAW transformations from imagej-bunwarpj plugin


#Load R script with functions that you'll need...
source("/home/hex/Unikram/GEO411 - Fernerkundung und Geoinf/Photogrammetry Laurichard Rock glacier/Scripts/process_BUnwarpJ_2D_image_registrationV02.R")

# Load required rasters
setwd("/home/hex/Unikram/GEO411 - Fernerkundung und Geoinf/Photogrammetry Laurichard Rock glacier")
library(raster)
library(RSAGA)
library(rgdal)
r_source <- raster("Data/DEM-L93-16Aug2012-LIDAR.tif")
r_target <- raster("Data/DEM-L93-05Oct2017-SFMMVS.tif")

# Get the IMCORR displacement field, pre-processed by SAGA
# what did SAGA do?
# 1. calculate displacement vectors -> line shape files
# 2. convert line features to points, at the position of 2012 -> point feature shape files
# 2.a. (optional) Selection of points with a real displacement between -60° and 60° and a slope <= 0°
# 3. Inverse Distance Weighting with a 10m bandwidth -> grid

# double check....
# d_tx <- bunwarpjDisplacementField(tx_file =  "Hillshade-L93-16Aug2012-LIDAR_direct_transf_auto_raw.txt", 
#                                r_source, r_target, is_inverse=FALSE)
#d_tx <- bunwarpjDisplacementField(tx_file =  "Hillshade-L93-05Oct2017-SFMMVS_inverse_transf_auto_raw.txt",
#                                  r_source, r_target, is_inverse=TRUE)
#

rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-pts-DEM_DISP_VEC_256_64_01.sgrd", out.grids="SAGA/idw-pts-DEM_DISP_VEC_256_64_01.asc", prec=1, out.path=getwd())
d_tx <- readGDAL("SAGA/idw-pts-DEM_DISP_VEC_256_64_01.asc")

# first inspection...
plot(d_tx)
summary(d_tx$band1)
str(d_tx$band1)
hist(d_tx$band1)

raster <- raster(d_tx)

# Save d_tx for both
# saveRDS(d_tx, file = "d_tx_direct_transf_auto.RDS")
# d_tx = readRDS(file = "d_tx_direct_transf_auto.RDS")
#saveRDS(d_tx, file = "d_tx_inverse_transf_auto.RDS")
#d_tx = readRDS(file = "d_tx_inverse_transf_auto.RDS")
#r_xy <- setValues(r_source, d_tx$xy_disp)
#raster=writeRaster(r_xy, filename="xyDisp_demDirect_auto.tif", format="GTiff", overwrite=TRUE)

# only uncomment for saving the Geotiff:
#raster=writeRaster(r_xy, filename="xyDisp_demInverse_auto.tif", format="GTiff", overwrite=TRUE)
plot(raster)

