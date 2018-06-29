# Scritp to read RAW transformations from imagej-bunwarpj plugin

source("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Scripts/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

library(raster)
library(Metrics)
library(caret)
library(rgdal)
library(RStoolbox)

#Set wdir
maindir="D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data/ImageJ"
raw = "raw"

################################################################################################################################
### 1.) Make tiffs readable in ImageJ

# geoTiffForImageJ(sprintf("%s/%s/DEM-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/DEM_L93_051017_SFMMVS_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Hillshade-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/Hillshade_L93_051017_SFMMVS_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Slope-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/Slope_L93_051017_SFMMVS_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Aspect-L93-05Oct2017-SFMMVS.tif",maindir,raw), sprintf("%s/Aspect_L93_051017_SFMMVS_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/DEM-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/DEM_L93_160812_LIDAR_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Hillshade-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/Hillshade_L93_160812_LIDAR_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Slope-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/Slope_L93_160812_LIDAR_IJ.tif",maindir))
# geoTiffForImageJ(sprintf("%s/%s/Aspect-L93-16Aug2012-LIDAR.tif",maindir,raw), sprintf("%s/Aspect_L93_160812_LIDAR_IJ.tif",maindir))

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
tiffFiles <- list.files(maindir,pattern="*IJ.tif$", full.names = TRUE)
# Create Source and Target Stack
r_source_stack = raster::stack(tiffFiles[c(4,6,2,16,8)]) # DEM, Hill, Aspect, Slope , PC1 2012
r_target_stack = raster::stack(tiffFiles[c(3,5,1,15,7)]) # DEM, Hill, Aspect, Slope , PC1 2017

# Calculate PCA

# pca_source = rasterPCA(r_source_stack, nSamples = NULL)
# pca_target = rasterPCA(r_target_stack, nSamples = NULL)
# plot(pca_source$map)
# summary(pca_source$model)
# plot(pca_target$map)
# summary(pca_target$model)

# for (k in 1:length(names(pca_source$map))){
#   
#   writeRaster(pca_source$map[[k]],filename = paste(maindir,"/",names(pca_source$map)[k],"_",substr(names(r_source_stack)[1],5,nchar(names(r_source_stack)[1])),sep=""), format="GTiff", overwrite=TRUE)
#   print(paste("Write_Raster",k,sep=""))
#   
# }
# 
# for (l in 1:length(names(pca_target$map))){
#   
#   writeRaster(pca_target$map[[l]],filename = paste(maindir,"/",names(pca_target$map)[l],"_",substr(names(r_target_stack)[1],5,nchar(names(r_target_stack)[1])),sep=""), format="GTiff", overwrite=TRUE)
#   print(paste("Write_Raster",k,sep=""))
#   
# }


## 3.2) Load Validation Movement Data (Bodin et al. 2017) Calculate XY Displacement (Euclidic Distanz)

val_df = read.csv(sprintf("%s/Laurichard_2014_L93_extract_velocity.csv",maindir), header = TRUE, sep = ",")

## 3.3) Load Landmarks

landmark_60 = sprintf("%s/Landmarks_60.txt",maindir)
landmark_61 = sprintf("%s/Landmarks_61.txt",maindir)
landmark_val = sprintf("%s/Landmarks_val_new.txt",maindir)
#landmark_val = sprintf("%s/Landmarks_val.txt",maindir)
#landmark_val = sprintf("%s/Landmarks_val_60.txt",maindir)
## 3.4) Load raw Transformationfiles
rawFiles_direct = list.files(maindir,pattern="*_direct_transf_raw.txt$", full.names = TRUE)
rawFiles_inverse = list.files(maindir,pattern="*_inverse_transf_raw.txt$", full.names = TRUE)

##3.5) Load Contour Shapefile (Glazier)

contour_shp = readOGR(dsn = maindir, layer = "laurichard_rg_contour")


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


landmark_df_val = imagejPointsToCRS2(landmark_val, r_source_stack[[1]])

landmarks_60_val_source = make_spdf(landmark_df_val[,c(2,3)],landmark_df_val,r_source_stack[[1]])
landmarks_60_val_source[,1]
landmarks_60_val_target = make_spdf(landmark_df_val[,c(4,5)],landmark_df_val,r_source_stack[[1]])
landmarks_60_val_target[,1]

# Calculate Dispalcements XY /XYZ
landmark_df_val["XY_Displacement"] = eucDist(landmark_df_val$x_crs_Source,landmark_df_val$x_crs_Target,
                                             landmark_df_val$y_crs_Source,landmark_df_val$y_crs_Target)
landmark_df_val["XY_Displacement_year"] = landmark_df_val["XY_Displacement"] / 6
landmark_df_val["z_source"] = extract(r_source_stack[[1]],landmarks_60_val_source)
landmark_df_val["z_target"] = extract(r_target_stack[[1]],landmarks_60_val_target)
landmark_df_val["Z_Displacement"] = landmark_df_val$z_target - landmark_df_val$z_source
landmark_df_val["XYZ_Displacement"] = sqrt( (landmark_df_val$XY_Displacement)^2 + (landmark_df_val$Z_Displacement)^2 )
landmark_df_val["XYZ_Displacement_year"] = landmark_df_val["XYZ_Displacement"] / 6

landmarks_60_val_source = make_spdf(landmark_df_val[,c(2,3)],landmark_df_val,r_source_stack[[1]])
landmarks_60_val_source[,1]
landmarks_60_val_target = make_spdf(landmark_df_val[,c(4,5)],landmark_df_val,r_source_stack[[1]])
landmarks_60_val_target[,1]


# Save to Shapefile

writeOGR(landmarks_60_source,dsn = maindir,layer ="Landmarks_60_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_60_target,dsn = maindir,layer ="Landmarks_60_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_61_source,dsn = maindir,layer ="Landmarks_61_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_61_target,dsn = maindir,layer ="Landmarks_61_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_60_val_source,dsn = maindir,layer ="Landmarks_60_val_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(landmarks_60_val_target,dsn = maindir,layer ="Landmarks_60_val_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )

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

# Use different Transformations, apply on DEM
Sys.time()
mybunwarpjDisplacementField(rawFiles_direct[1],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Aspect_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[2],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[3],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_61_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[4],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Hill_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[5],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Hill_61_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[6],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_PCA1_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[7],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Slope_60_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[8],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_0_0_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[9],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_0_1_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[10],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_0_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[11],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_1_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[12],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_2_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[13],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_4_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[14],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_2_1_IJ_direct_d_tx")
mybunwarpjDisplacementField(rawFiles_direct[15],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_4_1_IJ_direct_d_tx")
Sys.time() # ca. 45 min
gc()

#5.2) Inverse
# mybunwarpjDisplacementFieldinverse = function(path,r_source,r_target, maindir,saveNAME){
#   
#   d_tx <- bunwarpjDisplacementField(tx_file =  path,
#                                     r_source, r_target, is_inverse=TRUE)
#   saveRDS(d_tx, file = paste(maindir,"/",saveNAME,".RDS", sep =""))
#   
# }
# 
# Sys.time()
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[1],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Aspect_60_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[2],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_60_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[3],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_DEM_61_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[4],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Hill_60_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[5],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Hill_61_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[6],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_PCA1_60_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[7],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T1_Slope_60_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[8],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_0_0_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[9],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_0_1_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[10],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_0_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[11],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_1_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[12],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_2_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[13],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_1_4_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[14],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_2_1_IJ_inverse_d_tx")
# mybunwarpjDisplacementFieldinverse(rawFiles_inverse[15],r_source_stack[[1]],r_target_stack[[1]], maindir ,"T2_DEM_60_4_1_IJ_inverse_d_tx")
# Sys.time()
# gc()

#5.3) Read in d_tx, Calculate Anual Displacement xy and xyz, Write Displacements to Raster, 
#     Plot Displacement(pdf), Plot (Scatterplot), Calculate Stats (filtered Raster, Validierung), Write Stats csv

direct_dt_xFiles <- list.files(maindir,pattern="*direct_d_tx.RDS$", full.names = TRUE)
inverse_dt_xFiles <- list.files(maindir,pattern="*inverse_d_tx.RDS$", full.names = TRUE)


Sys.time()
# Where to Filter?!
for (i in 1:length(direct_dt_xFiles)){
#for (i in 1:2){ 
  ########################################################### Modell Raster Stack with Transformations ##########################################
  
  # Read Transformation DataFrame
  d_tx = readRDS(file = direct_dt_xFiles[i]) 
 
  #Calculate Anual Displacement in m/year (XY,XYZ)
  d_tx$xy_disp_annual = d_tx$xy_disp / 6
  d_tx$xyz_disp_annual = d_tx$xyz_disp / 6
  #head(d_tx)
  
  # Set Names for RasterStack and Outputfiles
  r_name = c(paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDisp",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDisp",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDispAnual",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDispAnual",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Aspect",sep =""),
             paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Slope",sep ="")
             )
  #Stack Raster
  r_stack = raster::brick(list(setValues(r_source_stack[[1]], d_tx$xy_disp),setValues(r_source_stack[[1]], d_tx$xyz_disp),
                               setValues(r_source_stack[[1]], d_tx$xy_disp_annual),setValues(r_source_stack[[1]], d_tx$xyz_disp_annual),
                               setValues(r_source_stack[[1]], d_tx$aspect),setValues(r_source_stack[[1]], d_tx$slope)))  
  names(r_stack) = r_name
  ###############################################################################################################################################
  
  ######################################################### Save and Plot unfilterd Modelled Raster Data ########################################
  r_stack_unfilterd_unmasked = r_stack
  
  # Plot unfilterd Modelled Raster Stack
  plot_raster_stack(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_unfilterd_r_stack.pdf",sep = ""),r_stack_unfilterd_unmasked)
  # Plot Histogramm of Moddelled Raster Stack
  plot_raster_stack_hist(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_unfilterd_r_stack_hist.pdf",sep = ""),r_stack_unfilterd_unmasked, 300)
  ###############################################################################################################################################
  
  ######################################################### Mask and Filter Modelled Raster Data ################################################
  
  # Mask RasterStack with Contour of Glazier
  
  # Sys.time()
  # r_stack = clip_raster_Stack(r_stack,contour_shp,r_name)
  # Sys.time()
  
  library(foreach)
  library(doParallel) 
  cl<-makeCluster(4)
  registerDoParallel(cl)
  
  Sys.time()
  mask_ls = foreach(i = 1:length(names(r_stack)),.packages = "raster",.inorder = TRUE) %dopar% {
    
    mask_raster_layer(r_stack[[i]],contour_shp,r_name[i])
    
  }
  Sys.time()
  stopCluster(cl)
  gc()
  
  r_stack = raster::stack(mask_ls)
  
  ##########
  
  
  
  r_stack_unfilterd_masked = r_stack
  
  # Plot masked RasterStack
  plot_raster_stack(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_masked_r_stack.pdf",sep = ""),r_stack_unfilterd_masked)
  # Plot Histogramm of masked RasterStack
  plot_raster_stack_hist(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_masked_r_stack_hist.pdf",sep = ""),r_stack_unfilterd_masked, 300)
  
  ###############################
  
  ## Filter with Aspect 
  #r_stack = r_stack_unfilterd_masked
  r_stack = filter_raster_stack(r_stack,60,300)  #In function possibility to subset stack
  
  # Plot masked and filterd  RasterStack
  plot_raster_stack(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_masked_filterd_r_stack.pdf",sep = ""),r_stack)
  # Plot masked and filterd Histogramm of masked RasterStack
  plot_raster_stack_hist(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_masked_filtered_r_stack_hist.pdf",sep = ""),r_stack, 300)
  ###################################################################################################################
  
  ##################################### Create Summary of filterd and masked RasterStack (Write to csv) #############
  #####################################               Plot pdf with valdidation points                  #############
  r_stack = raster::brick(r_stack)
  # Summary
  summary_raster_Stack = round(summary(r_stack, na.rm = T), digits = 4)
  summary_raster_Stack = data.frame(unclass(summary_raster_Stack), check.names = FALSE, stringsAsFactors = FALSE)
  
  # IQR
  IQR_final = calc_IQR_Raster(r_stack)

  # SD
  SD_final= calc_SD_Raster(r_stack)
  
  summary_raster_Stack["IQR",] = IQR_final
  summary_raster_Stack["SD",] = SD_final
  
  write.csv(summary_raster_Stack,paste(maindir,"/",names(r_stack)[3],"_Raster_Stats.csv", sep =""))
  
  
  # Write Raster
  for (j in 1:length(names(r_stack))){
    
    writeRaster(r_stack[[j]],filename = paste(maindir,"/",names(r_stack)[j],sep=""), format="GTiff", overwrite=TRUE)
    print(paste("Write_Raster",j,sep=""))
    
  }
  
  # Plot Data RasterStack, with Points (Source,Target,Validation)
  
  my_plot_function(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xy_z_displacment.pdf",sep = ""),
                   r_stack,val_spdf_source,landmarks_60_source,landmarks_60_target,landmarks_61_source,landmarks_61_target,landmarks_60_val_source,landmarks_60_val_target)
  
  ###################################################################################################################################
  
  ################################################### Extracting Valdidation Points #################################################
  
  # Extract Values from modelled Images with Valdidation Points
  xy_disp_annual_4_val = raster::extract(r_stack[[3]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE) #Bodin XY WRONG
  xyz_disp_annual_4_val = raster::extract(r_stack[[4]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE) # Bodin XYZ WRIGHT
  
  xy_disp_annual_4_val_user = raster::extract(r_stack[[3]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE) # Manuell XY 
  xyz_disp_annual_4_val_user = raster::extract(r_stack[[4]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE) # Manuell XYZ 
  ####################################################################################################################################
  
  ########################################### Create Validatin DF for Bodin Data (3D displacement)####################################
  ########################################### Create Validiaton DF for Expert based Approach      ####################################
  
  ############ BODIN et al. 2017
  full_val_df = as.data.frame(cbind(val_spdf_source[,c(3,8)],
                                    xy_disp_annual_4_val[,12],xyz_disp_annual_4_val[,12]))
  
  # Calculate Displacement Error (dm-dp)
  full_val_df[paste(colnames(full_val_df)[3],"Error",sep="")] = full_val_df[colnames(full_val_df[2])] - full_val_df[colnames(full_val_df[3])] ## VORSIcHT (Bodin 3D?!)
  full_val_df[paste(colnames(full_val_df)[4],"Error",sep="")] = full_val_df[colnames(full_val_df[2])] - full_val_df[colnames(full_val_df[4])] ## VORSIcHT (Bodin 3D?!)
  
  full_val_df["X"] = NULL
  full_val_df["Y"] = NULL
  # Plot Scatterplot
  my_scatterplot(paste(maindir,"/",colnames(full_val_df)[3],"_Scatterplot_Bodin.pdf", sep =""), full_val_df,
                 colnames(full_val_df)[3],colnames(full_val_df)[2],
                 colnames(full_val_df)[4],colnames(full_val_df)[2],
                 "XY Displacement BunwrapJ [m/year]","XY Displacement Bodin et al. 2017 [m/year]",
                 "XYZ Displacement BunwrapJ [m/year]","XYZ Displacement Bodin et al. 2017 [m/year]")
  
  #### Calculate Statistics 
  # R²
  # RMSE
  # Summary
  # IQR
  # Error Displacement (inkl. Summary)
  # SD
  final_stats = my_statistic_calc(full_val_df,colnames(full_val_df)[2],colnames(full_val_df)[3],
                                  colnames(full_val_df)[4],colnames(full_val_df)[5],
                                  colnames(full_val_df)[6])  
  
  # Write Statistics to csv
  write.csv(final_stats,paste(maindir,"/",colnames(full_val_df)[3],"_Stats_Bodin.csv", sep =""))
  ############
  
  ############ Manuell Approach
  
  # Calclulate Validiation vor Manual Approch
  
  full_val_manuell_df = as.data.frame(cbind(xy_disp_annual_4_val_user[,c(1,7,12,13)],xyz_disp_annual_4_val_user[,13]))
  full_val_manuell_df["x_crs_Source"] = NULL
  full_val_manuell_df["y_crs_Source"] = NULL
  
  # Calculate Displacement Error 2D / 3D
  
  full_val_manuell_df[paste(colnames(full_val_manuell_df)[4],"Error",sep="")] = full_val_manuell_df[colnames(full_val_manuell_df[2])] - full_val_manuell_df[colnames(full_val_manuell_df[4])] 
  full_val_manuell_df[paste(colnames(full_val_manuell_df)[5],"Error",sep="")] = full_val_manuell_df[colnames(full_val_manuell_df[3])] - full_val_manuell_df[colnames(full_val_manuell_df[5])] 
  
  my_scatterplot(paste(maindir,"/",colnames(full_val_manuell_df)[4],"_Scatterplot_manuell.pdf", sep =""), full_val_manuell_df,
                 colnames(full_val_manuell_df)[4],colnames(full_val_manuell_df)[3],
                 colnames(full_val_manuell_df)[5],colnames(full_val_manuell_df)[3],
                 "XY Displacement BunwrapJ [m/year]","XY Displacement Expert based [m/year]",
                 "XYZ Displacement BunwrapJ [m/year]","XYZ Displacement Expert based [m/year]")
  
  final_stats_manuell = my_statistic_calc_expert(full_val_manuell_df,colnames(full_val_manuell_df)[2],colnames(full_val_manuell_df)[3],
                                                 colnames(full_val_manuell_df)[4],colnames(full_val_manuell_df)[5],
                                                 colnames(full_val_manuell_df)[6],colnames(full_val_df)[7])  
   # Write Statistics to csv
  write.csv(final_stats_manuell,paste(maindir,"/",colnames(full_val_df)[3],"_Stats_manuell.csv", sep =""))
  
  
  print(paste("----------------------- Iteration",i,"completed-----------------------", sep = ""))
 }

Sys.time()

##################################################################################################################################
# # USE of IMCORR Results
# 
# imcorr_folder = "IMCORR" 
# 
# imcorr_raster = list.files(paste(maindir,"/",imcorr_folder, sep = ""), pattern = "*.sdat$", full.names = TRUE)
# idw_pts = imcorr_raster[1:3]
# idw_corr = imcorr_raster[4:6]
# pts_imcorr = list.files(paste(maindir,"/",imcorr_folder, sep = ""), pattern = "*.shp$", full.names = TRUE)
# 
# idw_pts = raster::stack(idw_pts)
# idw_corr = raster::stack(idw_corr)
# 
# #Convert to tiff
# for (l in 1:length(names(idw_pts))){
#   
#   writeRaster(idw_pts[[l]],filename = paste(maindir,"/",imcorr_folder,"/",names(idw_pts)[l],sep=""), format="GTiff", overwrite=TRUE)
#   print(paste("Write_Raster",l,sep=""))
#   
# }
# 
# for (m in 1:length(names(idw_corr))){
#   
#   writeRaster(idw_corr[[m]],filename = paste(maindir,"/",imcorr_folder,"/",names(idw_corr)[m],sep=""), format="GTiff", overwrite=TRUE)
#   print(paste("Write_Raster",m,sep=""))
#   
# }
# # REPLACE MANUEL . with _
# 
# imcorr_raster = list.files(paste(maindir,"/",imcorr_folder, sep = ""), pattern = "*.tif$", full.names = TRUE)
# 
# geoTiffForImageJ(imcorr_raster[1], paste(substring(imcorr_raster[1],0,nchar(imcorr_raster[1])-4),"_IJ.tif", sep =""))
# geoTiffForImageJ(imcorr_raster[2], paste(substring(imcorr_raster[2],0,nchar(imcorr_raster[2])-4),"_IJ.tif", sep =""))
# geoTiffForImageJ(imcorr_raster[3], paste(substring(imcorr_raster[3],0,nchar(imcorr_raster[3])-4),"_IJ.tif", sep =""))
# geoTiffForImageJ(imcorr_raster[4], paste(substring(imcorr_raster[4],0,nchar(imcorr_raster[4])-4),"_IJ.tif", sep =""))
# geoTiffForImageJ(imcorr_raster[5], paste(substring(imcorr_raster[5],0,nchar(imcorr_raster[5])-4),"_IJ.tif", sep =""))
# geoTiffForImageJ(imcorr_raster[6], paste(substring(imcorr_raster[6],0,nchar(imcorr_raster[6])-4),"_IJ.tif", sep =""))
# 
# imcorr_raster = list.files(paste(maindir,"/",imcorr_folder, sep = ""), pattern = "*IJ.tif$", full.names = TRUE)
# 
# idw_pts = imcorr_raster[1:3]
# idw_corr = imcorr_raster[4:6]
# 
# 
# idw_pts = raster::stack(idw_pts)
# idw_corr = raster::stack(idw_corr)
# 
# # imcorrOut <- readOGR(pts_imcorr[1], substr(pts_imcorr[1],(nchar(maindir)+nchar(imcorr_folder)+3),(nchar(pts_imcorr)[1]-4)))
# # 
# # ## Estimate displacement from Aug 2012 to  Oct 2017 #######################################
# # 
# # # Create new dataframe
# # d <- data.frame(total_disp=imcorrOut$DISP)
# # d$slope <- imcorrOut$SLOPE
# # d$aspect <- imcorrOut$ASPECT
# # d$x_real <- imcorrOut$REALX
# # d$y_real <- imcorrOut$REALY
# # d$z_real <- imcorrOut$REALZ
# # d$x_target <- imcorrOut$XTARG 
# # d$y_target <- imcorrOut$YTARG
# # d$z_target <- imcorrOut$ZTARG
# # d$disp_real <- imcorrOut$DISP_REAL
# 
# 
# 
# 
# 
# plot_raster_stack(paste(maindir,"/",imcorr_folder,"/idw_pts-DEM_DISP_unfilterd.pdf", sep =""),idw_pts)
# plot_raster_stack(paste(maindir,"/",imcorr_folder,"/idw_DEM_DISP_filterd.pdf", sep =""),idw_corr)
# 
# # Create Summary of filterd IMCORR-Raster Stack
# 
# summary_raster_Stack_imcorr = summary(getValues(idw_corr), na.rm = T,digits = 20)
# summary_raster_Stack_imcorr = data.frame(unclass(summary_raster_Stack_imcorr), check.names = FALSE, stringsAsFactors = FALSE)
# 
# # IQR
# IQR_1_IMCORR = round(IQR(getValues(idw_corr[[1]]), na.rm = T, type = 7), digits = 4)
# IQR_2_IMCORR =round(IQR(getValues(idw_corr[[2]]), na.rm = T, type = 7), digits = 4)
# IQR_3_IMCORR =round(IQR(getValues(idw_corr[[2]]), na.rm = T, type = 7), digits = 4)
# 
# IQR_finale = t(as.data.frame(c(IQR_1_IMCORR,IQR_2_IMCORR,IQR_3_IMCORR)))
# 
# # SD
# SD_1_IMCORR = round(sd(getValues(idw_corr[[1]]), na.rm = T), digits = 4)
# SD_2_IMCORR = round(sd(getValues(idw_corr[[2]]), na.rm = T), digits = 4)
# SD_3_IMCORR = round(sd(getValues(idw_corr[[3]]), na.rm = T), digits = 4)
# 
# 
# SD_final = t(as.data.frame(c(SD_1_IMCORR,SD_2_IMCORR,SD_3_IMCORR)))
# summary_raster_Stack_imcorr["IQR",] = IQR_finale
# summary_raster_Stack_imcorr["SD",] = SD_final
# 
# write.csv(summary_raster_Stack,paste(maindir,"/",imcorr_folder,"/","idw_DEM_Raster_Stats_filtered.csv", sep =""))
# 
# 
# xyz_disp_annual_4_val_IMCORR_1 = raster::extract(idw_corr[[1]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# xyz_disp_annual_4_val_IMCORR_2 = raster::extract(idw_corr[[2]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# xyz_disp_annual_4_val_IMCORR_3 = raster::extract(idw_corr[[3]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# 
# xy_disp_annual_4_val_user_IMCORR_1 = raster::extract(idw_corr[[1]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# xy_disp_annual_4_val_user_IMCORR_2 = raster::extract(idw_corr[[2]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# xy_disp_annual_4_val_user_IMCORR_3 = raster::extract(idw_corr[[3]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
# 
# 
# 
# # Bodin !!!! 3D displacement
# # Create Validation Data Frame
# full_val_spdf_IMCORR = cbind(val_spdf_source[,c(3,8)],
#                       xyz_disp_annual_4_val_IMCORR_1[,12],xyz_disp_annual_4_val_IMCORR_2[,12],xyz_disp_annual_4_val_IMCORR_3[,12])
# 
# #names(full_val_spdf)= c(colnames(val_spdf_source@data[c(3,8)]),colnames(xy_disp_annual_4_val@data[12]),colnames(xyz_disp_annual_4_val@data[12]))
# full_val_df_IMCORR = as.data.frame(full_val_spdf_IMCORR)
# 
# # Calculate Displacement Error (dm-dp)
# full_val_df_IMCORR[paste(colnames(full_val_df_IMCORR)[3],"Error",sep="")] = full_val_df_IMCORR[colnames(full_val_df_IMCORR[2])] - full_val_df_IMCORR[colnames(full_val_df_IMCORR[3])] ## VORSIcHT (Bodin 3D?!)
# full_val_df_IMCORR[paste(colnames(full_val_df_IMCORR)[4],"Error",sep="")] = full_val_df_IMCORR[colnames(full_val_df_IMCORR[2])] - full_val_df_IMCORR[colnames(full_val_df_IMCORR[4])] ## VORSIcHT (Bodin 3D?!)
# full_val_df_IMCORR[paste(colnames(full_val_df_IMCORR)[5],"Error",sep="")] = full_val_df_IMCORR[colnames(full_val_df_IMCORR[2])] - full_val_df_IMCORR[colnames(full_val_df_IMCORR[5])]
# 
# 
# # Plot Scatterplot
# my_scatterplot_IMCORR(paste(maindir,"/",imcorr_folder,"/",colnames(full_val_df_IMCORR)[3],"_Scatterplot_Bodin.pdf", sep =""), full_val_df_IMCORR,
#                colnames(full_val_df_IMCORR)[3],colnames(full_val_df_IMCORR)[2],
#                colnames(full_val_df_IMCORR)[4],colnames(full_val_df_IMCORR)[2],
#                colnames(full_val_df_IMCORR)[5],colnames(full_val_df_IMCORR)[2],
#                "XYZ Displacement Bodin et al. 2017 [m/year]")
# 
# ##### AB HIER
# 
# final_stats_IMCORR = my_statistic_calc(full_val_df_IMCORR,colnames(full_val_df_IMCORR)[2],colnames(full_val_df)[3],
#                                 colnames(full_val_df)[4],colnames(full_val_df)[7],
#                                 colnames(full_val_df)[8])  
# 
# # Write Statistics to csv
# write.csv(final_stats,paste(maindir,"/",colnames(full_val_df)[3],"_Stats_Bodin.csv", sep =""))

