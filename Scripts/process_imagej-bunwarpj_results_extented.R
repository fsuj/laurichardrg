# Scritp to read RAW transformations from imagej-bunwarpj plugin

source("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Scripts/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

library(raster)
library(Metrics)
library(caret)
library(rgdal)

#Set wdir
maindir="D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data/ImageJ"
raw = "raw"

################################################################################################################################
### 1.) Make tiffs readable in ImageJ

geoTiffForImageJ(sprintf("%s/%s/DEM-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/DEM_L93_051017_SFMMVS_IJ.tif",maindir))
geoTiffForImageJ(sprintf("%s/%s/Hillshade-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/Hillshade_L93_051017_SFMMVS_IJ.tif",maindir))
geoTiffForImageJ(sprintf("%s/%s/DEM-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/DEM_L93_160812_LIDAR_IJ.tif",maindir))
geoTiffForImageJ(sprintf("%s/%s/Hillshade-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/Hillshade_L93_160812_LIDAR_IJ.tif",maindir))

################################################################################################################################
### 2.) Process Data with Fiji/BunwrapJ
# - Convert tiffs with ArcMap to PNG
# - Opening order 2012 / 2017
# - Create Landmark and Save Landmarks (Use Hillshade) - play with different Landmarks (extreme) Test 1
# - Open resaved tiffs (2012/2017)
# - Plugins - Registration - BunwrapJ 
# - Set Source Image (2012) and Target (2017)   - play with differnt Derivates of DEM Test 1
# - Settings: 
#     - Registration Mode: Accurate
#     - Image Subsampling Factor: 0
#     - Initial Deformation: Fine
#     - Final Deformation: SuperFine
#     - Divergence Weight: 0.0
#     - Curl Weight: 0.0
#     - Landmark Weight: 1.0          - Settings to play with Test 2
#     - Image Weight: 1.0             - Settings to play with Test 2
#     - consistency Weight: 10.0
#     - Stop Threshold: 0.01
#     Check Verbose
#     Check Save Transformation
# - Open Landmarks                    Different (extreme Values), or use for "expert validation Approach"
# - Process
# - convert Transformation (landmarks) to raw
################################################################################################################################

#### 3.) Load Data

## 3.1) Load Tiffs

# List all .tifs in maindir
tiffFiles <- list.files(maindir,pattern="*.tif$", full.names = TRUE)
# Create Source and Target Stack
r_source_stack = raster::stack(tiffFiles[c(2,4)]) # DEM, Hill 2012
r_target_stack = raster::stack(tiffFiles[c(1,3)]) # DEM, Hill 2017

## 3.2) Load Validation Movement Data (Bodin et al. 2017) Calculate XY Displacement (Euclidic Distanz)

val_df = read.csv(sprintf("%s/Laurichard_2014_L93_extract_velocity.csv",maindir), header = TRUE, sep = ",")

## 3.3) Load Landmarks

landmark_60 = sprintf("%s/Landmarks_60.txt",maindir)
landmark_61 = sprintf("%s/Landmarks_61.txt",maindir)

## 3.4) Load raw Transformationfiles
rawFiles_direct = list.files(maindir,pattern="*_direct_transf_raw.txt$", full.names = TRUE)
rawFiles_inverse = list.files(maindir,pattern="*_inverse_transf_raw.txt$", full.names = TRUE)
################################################################################################################################
### 4.) Convert Data to spdf (Mapping)

# 4.1 ) convert Landmarks to Point Shapefile (Mapping)
# Extract Coords from Pixel Coords
landmark_df_60 = imagejPointsToCRS2(landmark_60, r_source_stack[[1]])
landmark_df_61 = imagejPointsToCRS2(landmark_61, r_source_stack[[1]])
# Make spdf
landmarks_60_source  = make_spdf(landmark_df_60[,c(2,3)],landmark_df_60,r_source_stack[[1]])
landmarks_60_source[,1]
landmarks_60_target  = make_spdf(landmark_df_60[,c(4,5)],landmark_df_60,r_source_stack[[1]])
landmarks_60_target[,1]
landmarks_61_source  = make_spdf(landmark_df_61[,c(2,3)],landmark_df_61,r_source_stack[[1]])
landmarks_61_source[,1]
landmarks_61_target  = make_spdf(landmark_df_61[,c(4,5)],landmark_df_61,r_source_stack[[1]])
landmarks_61_target[,1]
# Save to Shapefile

writeOGR(landmarks_60_source,dsn = maindir,layer ="Landmarks_60_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_60_target,dsn = maindir,layer ="Landmarks_60_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_61_source,dsn = maindir,layer ="Landmarks_61_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_61_target,dsn = maindir,layer ="Landmarks_61_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )

# 4.2) Convert Validation Data to Point Shapefile

val_spdf_source = make_spdf(val_df[,c(1,2)],val_df,r_source_stack[[1]])
val_spdf_source = crop(val_spdf_source, extent(r_source_stack[[1]])) # Get Points within Raster

val_spdf_target = make_spdf(val_df[,c(4,5)],val_df,r_source_stack[[1]])
val_spdf_target = crop(val_spdf_target, extent(r_source_stack[[1]])) # Get Points within Raster       ####### CHECK THIS STUPID EXTNET?!

# komisch sehen Gleich aus?! obwohl andere Koordinaten
writeOGR(val_spdf_source,dsn = maindir,layer ="Validation_points_source", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(val_spdf_target,dsn = maindir,layer ="Validation_points_target", driver="ESRI Shapefile", overwrite_layer=TRUE)

################################################################################################################################
### 5.) Process BUnwrapJ
#5.1) Direct
mybunwarpjDisplacementField = function(path,r_source,r_target, maindir,saveNAME){
  
  d_tx <- bunwarpjDisplacementField(tx_file =  path,
                                    r_source, r_target, is_inverse=FALSE)
  saveRDS(d_tx, file = paste(maindir,"/",saveNAME,".RDS", sep =""))
  
}

Sys.time()
mybunwarpjDisplacementField(rawFiles_direct[1],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[2],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_61_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[3],r_source_stack[[2]],r_target_stack[[2]], maindir ,"T1_Hill_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[4],r_source_stack[[2]],r_target_stack[[2]], maindir ,"T1_Hill_61_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[5],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_1_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[6],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_2_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[7],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_4_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[8],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_2_1_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[9],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_4_1_IJ_direct_d_tx")
Sys.time() # ca. 30 min
gc()

#5.2) Inverse
mybunwarpjDisplacementFieldinverse = function(path,r_source,r_target, maindir,saveNAME){
  
  d_tx <- bunwarpjDisplacementField(tx_file =  path,
                                    r_source, r_target, is_inverse=TRUE)
  saveRDS(d_tx, file = paste(maindir,"/",saveNAME,".RDS", sep =""))
  
}

Sys.time()
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[1],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_60_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[2],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_61_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[3],r_source_stack[[2]],r_target_stack[[2]], maindir ,"T1_Hill_60_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[4],r_source_stack[[2]],r_target_stack[[2]], maindir ,"T1_Hill_61_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[5],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_1_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[6],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_2_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[7],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_4_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[8],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_2_1_IJ_inverse_d_tx")
mybunwarpjDisplacementFieldinverse(rawFiles_inverse[9],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_4_1_IJ_inverse_d_tx")
Sys.time()
gc()

#5.3) Read in d_tx, Calculate Anual Displacement xy and xyz, Write Displacements to Raster, 
#     Plot Displacement(pdf), Calculate Stats, Write Stats csv

direct_dt_xFiles <- list.files(maindir,pattern="*direct_d_tx.RDS$", full.names = TRUE)
inverse_dt_xFiles <- list.files(maindir,pattern="*inverse_d_tx.RDS$", full.names = TRUE)

# Where to Filter?!
for (i in 1:length(direct_dt_xFiles)){
  
  # Read Transformation DataFrame
  d_tx = readRDS(file = direct_dt_xFiles[i]) 
 
  #Calculate Anual Displacement in m/year (XY,XYZ)
  d_tx$xy_disp_annual = d_tx$xy_disp / 6
  d_tx$xyz_disp_annual = d_tx$xyz_disp / 6
  head(d_tx)
  
  # Set Names for RasterStack and Outputfiles
  r_name = c(paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDisp",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDisp",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDispAnual",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDispAnual",sep ="")
             )
  #Stack Raster
  r_stack = raster::brick(list(setValues(r_source_stack[[1]], d_tx$xy_disp),setValues(r_source_stack[[1]], d_tx$xyz_disp),
                               setValues(r_source_stack[[1]], d_tx$xy_disp_annual),setValues(r_source_stack[[1]], d_tx$xyz_disp_annual)))
  names(r_stack) = r_name
  
  # Write Raster
  for (j in 1:length(names(r_stack))){
    
    writeRaster(r_stack[[j]],filename = paste(maindir,"/",names(r_stack)[j],sep=""), format="GTiff", overwrite=TRUE)
    print(paste("Write_Raster",j,sep=""))
    
  }
  
  # Plot Data RasterStack, with Points (Source,Target,Validation)
  
  my_plot_function(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xy_z_displacment.pdf",sep = ""),
                   r_stack,val_spdf_source,landmarks_60_source,landmarks_60_target,landmarks_61_source,landmarks_61_target)
  
  # Extract Values from modelled Images with Valdidation Points
  xy_disp_annual_4_val = raster::extract(r_stack[[3]],val_spdf_source, method = "simple", buffer = NULL, na.rm=FALSE, sp=TRUE)
  xyz_disp_annual_4_val = raster::extract(r_stack[[4]],val_spdf_source, method = "simple", buffer = NULL, na.rm=FALSE, sp=TRUE)
  
  
  #### Validation
  
  # Create Validation Data Frame
  full_val_spdf = cbind(val_spdf_source[,c(3,8)],xy_disp_annual_4_val[,12],xyz_disp_annual_4_val[,12])
  names(full_val_spdf)= c(colnames(val_spdf_source@data[c(3,8)]),colnames(xy_disp_annual_4_val@data[12]),colnames(xyz_disp_annual_4_val@data[12]))
  full_val_df = as.data.frame(full_val_spdf)
  
  # Calculate Displacement Error (dm-dp)
  full_val_df[paste(colnames(full_val_df)[3],"Error",sep="")] = full_val_df[colnames(full_val_df[2])] - full_val_df[colnames(full_val_df[3])] ## VORSIcHT (Bodin 3D?!)
  full_val_df[paste(colnames(full_val_df)[4],"Error",sep="")] = full_val_df[colnames(full_val_df[2])] - full_val_df[colnames(full_val_df[4])] ## VORSIcHT (Bodin 3D?!)
  
  # Plot Scatterplot
  my_scatterplot(paste(maindir,"/",colnames(full_val_df)[3],"_Scatterplot.pdf", sep =""), full_val_df,
                 colnames(full_val_df)[3],colnames(full_val_df)[2],
                 colnames(full_val_df)[4],colnames(full_val_df)[2])
  
  #### Calculate Statistics 
  # R²
  # RMSE
  # Summary
  # IQR
  # Error Displacement (inkl. Summary)
  # SD
  final_stats = my_statistic_calc(full_val_df,colnames(full_val_df)[2],colnames(full_val_df)[3],
                                  colnames(full_val_df)[4],colnames(full_val_df)[7],
                                  colnames(full_val_df)[8])  
  
  # Write Statistics to csv
  write.csv(final_stats,paste(maindir,"/",colnames(full_val_df)[3],"_Stats.csv", sep =""))
  
 }


