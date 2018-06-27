
# Scritp to read RAW transformations from imagej-bunwarpj plugin


#Load R script with functions that you'll need...
source("C:/YOUR_DIRECTORY/Ste-Flo/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

# Load required rasters
setwd("C:/YOUR_DIRECTORY/Ste-Flo/Data")
library(raster)
r_source <- raster("DEM-L93-16Aug2012-LIDAR.tif")
r_target <- raster("DEM-L93-05Oct2017-SFMMVS.tif")

# Get the bunwarpj transformation (i.e. displacement field)
setwd("D:/LOCATION_OF_TRANSFORMATIONS_RAW/")

d_tx <- bunwarpjDisplacementField(tx_file =  "bunwarpj-raw_transformation.txt", 
                               r_source, r_target, is_inverse=FALSE)

r_xyz <- setValues(r_source, d_tx$xy_disp)
writeRaster(r_xy, filename="xyzDisp_demDirect.tif", format="GTiff", overwrite=TRUE)


