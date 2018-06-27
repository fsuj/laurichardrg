
# Scritp to read RAW transformations from imagej-bunwarpj plugin


#Load R script with functions that you'll need...
source("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Scripts/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

# Load required rasters
setwd("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data")
library(raster)
r_source <- raster("DEM-L93-16Aug2012-LIDAR.tif")
r_target <- raster("DEM-L93-05Oct2017-SFMMVS.tif")

# Get the bunwarpj transformation (i.e. displacement field)
setwd("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data/")

# double check....
# d_tx <- bunwarpjDisplacementField(tx_file =  "Hillshade-L93-16Aug2012-LIDAR_direct_transf_auto_raw.txt", 
#                                r_source, r_target, is_inverse=FALSE)

d_tx <- bunwarpjDisplacementField(tx_file =  "Hillshade-L93-05Oct2017-SFMMVS_inverse_transf_auto_raw.txt",
                                  r_source, r_target, is_inverse=TRUE)

summary(d_tx)

# Save d_tx for both
# saveRDS(d_tx, file = "d_tx_direct_transf_auto.RDS")
# d_tx = readRDS(file = "d_tx_direct_transf_auto.RDS")

saveRDS(d_tx, file = "d_tx_inverse_transf_auto.RDS")
#d_tx = readRDS(file = "d_tx_inverse_transf_auto.RDS")

r_xy <- setValues(r_source, d_tx$xy_disp)
#raster=writeRaster(r_xy, filename="xyDisp_demDirect_auto.tif", format="GTiff", overwrite=TRUE)
raster=writeRaster(r_xy, filename="xyDisp_demInverse_auto.tif", format="GTiff", overwrite=TRUE)
plot(raster)

