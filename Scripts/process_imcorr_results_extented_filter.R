# Scritp to read RAW transformations from imagej-bunwarpj plugin

setwd("/home/hex/Unikram/GEO411 - Fernerkundung und Geoinf/Photogrammetry Laurichard Rock glacier")


#Load R script with functions that you'll need...
source("Scripts/process_BUnwarpJ_2D_image_registrationV02.R")

library(raster)
library(Metrics)
library(caret)
library(rgdal)
library(RSAGA)

#Set wdir
maindir=paste(getwd(), "/Data",sep="")
raw = "raw"

################################################################################################################################
### 1.) IMCORR processing in SAGA GIS

# Geoprocessing -> Grid -> Analysis -> IMCORR or: script files from Jason Goetz
# with different search/reference window sizes and fixed gridding (0.1)
# 128/64, 256/64, 256/128
# visual inspection: 256/64 gives the best results

# resulting in line feature shapes, indicating displacement vectors

################################################################################################################################
### 2.) Post-processing for gridding
# 
################################################################################################################################

#### 3.) Load Data

## 3.1) Load Tiffs

# List all .tifs in maindir
tiffFiles <- list.files(maindir,pattern="*.tif$", full.names = TRUE)
# Create Source and Target Stack
r_source <- raster("Data/DEM-L93-16Aug2012-LIDAR.tif")
r_target <- raster("Data/DEM-L93-05Oct2017-SFMMVS.tif")
r_source_stack = raster::stack(tiffFiles[c(2,5)]) # DEM, Hill 2012
r_target_stack = raster::stack(tiffFiles[c(1,4)]) # DEM, Hill 2017

## 3.2) Load Validation Movement Data (Bodin et al. 2017) Calculate XY Displacement (Euclidic Distanz)

val_df = read.csv(sprintf("%s/Laurichard_2014_L93_extract_velocity.csv",maindir), header = TRUE, sep = ",")

## 3.3) SAGA preprocessing (convert from sgrd to asc)
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-pts-DEM_DISP_VEC_256_64_01.sgrd", out.grids="SAGA/idw_pts_DEM_DISP_VEC_256_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-pts_DEM_DISP_VEC_256_128_01.sgrd", out.grids="SAGA/idw_pts_DEM_DISP_VEC_256_128_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-pts_DEM_DISP_VEC_128_64_01.sgrd", out.grids="SAGA/idw_pts_DEM_DISP_VEC_128_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-aspectcorr-DEM_DISP_VEC_256_64_01.sgrd", out.grids="SAGA/idw_aspectcorr_DEM_DISP_VEC_256_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-aspectcorr-DEM_DISP_VEC_256_128_01.sgrd", out.grids="SAGA/idw_aspectcorr_DEM_DISP_VEC_256_128_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/idw-aspectcorr-DEM_DISP_VEC_128_64_01.sgrd", out.grids="SAGA/idw_aspectcorr_DEM_DISP_VEC_128_64_01.asc", prec=1, out.path=getwd())

rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw-DEM_DISP_VEC_256_64_01.sgrd", out.grids="SAGA/xy_idw_DEM_DISP_VEC_256_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw-DEM_DISP_VEC_256_128_01.sgrd", out.grids="SAGA/xy_idw_DEM_DISP_VEC_256_128_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw-DEM_DISP_VEC_128_64_01.sgrd", out.grids="SAGA/xy_idw_DEM_DISP_VEC_128_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw-aspectcorr-DEM_DISP_VEC_256_64_01.sgrd", out.grids="SAGA/xy_idw_aspectcorr-DEM_DISP_VEC_256_64_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw-aspectcorr_DEM_DISP_VEC_256_128_01.sgrd", out.grids="SAGA/xy_idw_aspectcorr-DEM_DISP_VEC_256_128_01.asc", prec=1, out.path=getwd())
rsaga.sgrd.to.esri(in.sgrds="SAGA/xy-idw_aspectcorr_DEM_DISP_VEC_128_64_01.sgrd", out.grids="SAGA/xy_idw_aspectcorr-DEM_DISP_VEC_128_64_01.asc", prec=1, out.path=getwd())

## 3.4) Load displacement rasters
## in the order:
## 1. XYZ displacements
## 1.1 uncorrected
## 1.2 aspect corrected
## 2. XY displacements
## 2.1 uncorrected
## 2.2 corrected
ascFiles <- list.files(paste(getwd(), "/SAGA",sep="")
  ,pattern="*.asc$"
  , full.names = TRUE)
r_displacement_stack <- raster::stack(ascFiles)

## 3.5) Load Contour Shapefile (Glazier)
contour_shp = readOGR(dsn = maindir, layer = "laurichard_rg_contour")

################################################################################################################################
### first data inspection/plot
summary(raster)
str(raster$band1@data@values)
hist(raster$band1)


################################################################################################################################
### 4.) Convert Data to spdf (Mapping)

# 4.1) Convert Validation Data to Point Shapefile

val_spdf_source = make_spdf(val_df[,c(1,2)],val_df,r_source_stack[[1]])
val_spdf_source = crop(val_spdf_source, extent(r_source_stack[[1]])) # Get Points within Raster

val_spdf_target = make_spdf(val_df[,c(4,5)],val_df,r_source_stack[[1]])
val_spdf_target = crop(val_spdf_target, extent(r_source_stack[[1]])) # Get Points within Raster       ####### TODO: CHECK THIS STUPID EXTENT !!!

# komisch sehen Gleich aus?! obwohl andere Koordinaten
writeOGR(val_spdf_source,dsn = maindir,layer ="Validation_points_source", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(val_spdf_target,dsn = maindir,layer ="Validation_points_target", driver="ESRI Shapefile", overwrite_layer=TRUE)

################################################################################################################################
### 5.) Process IMCORR

#5.3) Read in d_tx, Calculate Anual Displacement xy and xyz, Write Displacements to Raster, 


# iterate through the three window size combinations:
for (i in 1:(nlayers(r_displacement_stack)/4)){
  
  raster_xyz_aspectcorr <- r_displacement_stack[[i]]*1
  raster_xyz <- r_displacement_stack[[i+3]]*1
  raster_xy <- r_displacement_stack[[i+6]]
  raster_xy_aspectcorr <- r_displacement_stack[[i+9]]
 
  #Calculate Anual Displacement in m/year (XY,XYZ)
  raster_xyz_aspectcorr_annual <- raster_xyz_aspectcorr / 6
  raster_xyz_annual <- raster_xyz / 6
  raster_xy_aspectcorr_annual <- raster_xy_aspectcorr / 6
  raster_xy_annual <- raster_xy / 6
  
  # Set Names for RasterStack and Outputfiles
  r_name = c(paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDisp",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDisp",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyDispAnual",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xyzDispAnual",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Aspect_xy",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Aspect_xyz",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Aspect_xyDispAnnual",sep ="")
            , paste(substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_Aspect_xyzDispAnnual",sep ="")
             )
  #Stack Raster
  r_stack <- brick(list(
    raster_xy
    ,raster_xyz
    ,raster_xy_annual
    ,raster_xyz_annual
    ,raster_xy_aspectcorr
    ,raster_xyz_aspectcorr
    ,raster_xy_annual_aspectcorr
    ,raster_xyz_annual_aspectcorr
  ))
  #r_stack = raster::brick(list(setValues(r_source_stack[[1]], d_tx$xy_disp),setValues(r_source_stack[[1]], d_tx$xyz_disp),
  #                             setValues(r_source_stack[[1]], d_tx$xy_disp_annual),setValues(r_source_stack[[1]], d_tx$xyz_disp_annual),
  #                             setValues(r_source_stack[[1]], d_tx$aspect),setValues(r_source_stack[[1]], d_tx$slope)))  
  names(r_stack) = r_name
  
  # Crop with Contour of Glazier
  # contour_shp
  
  r_stack = clip_raster_Stack(r_stack,contour_shp,r_name)
  
  ############################ TODO: calculate statistics from here on and write results to file...

  
  ## Create Summary of Images write to csv
  
  # Summary
  summary_raster_Stack = round(summary(r_stack, na.rm = T), digits = 4)
  summary_raster_Stack = data.frame(unclass(summary_raster_Stack), check.names = FALSE, stringsAsFactors = FALSE)
  
  # IQR
  IQR_1 = round(IQR(r_stack[[1]]@data@values, na.rm = T, type = 7), digits = 4)
  IQR_2 =round(IQR(r_stack[[2]]@data@values, na.rm = T, type = 7), digits = 4)
  IQR_3 =round(IQR(r_stack[[3]]@data@values, na.rm = T, type = 7), digits = 4)
  IQR_4 =round(IQR(r_stack[[4]]@data@values, na.rm = T, type = 7), digits = 4)
  IQR_finale = t(as.data.frame(c(IQR_1,IQR_2,IQR_3,IQR_4)))
  
  # SD
  SD_1 = round(sd(r_stack[[1]]@data@values, na.rm = T), digits = 4)
  SD_2 = round(sd(r_stack[[2]]@data@values, na.rm = T), digits = 4)
  SD_3 = round(sd(r_stack[[3]]@data@values, na.rm = T), digits = 4)
  SD_4 = round(sd(r_stack[[4]]@data@values, na.rm = T), digits = 4)
  
  SD_final = t(as.data.frame(c(SD_1,SD_2,SD_3,SD_4)))
  summary_raster_Stack["IQR",] = IQR_finale
  summary_raster_Stack["SD",] = SD_final
  
  write.csv(summary_raster_Stack,paste(maindir,"/",names(r_stack)[3],"_Raster_Stats.csv", sep =""))
  
  
  # Write Raster
  for (j in 1:length(names(r_stack))){
    
    writeRaster(r_stack[[j]],filename = paste(maindir,"/",names(r_stack)[j],sep=""), format="GTiff", overwrite=TRUE)
    print(paste("Write_Raster",j,sep=""))
    
  }
  
  # Plot Data RasterStack, with Points (Source,Target,Validation)
  
  my_plot_function(paste(maindir,"/",substr(direct_dt_xFiles[i], nchar(maindir)+2,nchar(direct_dt_xFiles[i])-9),"_xy_z_displacment.pdf",sep = ""),
                   r_stack,val_spdf_source,landmarks_60_source,landmarks_60_target,landmarks_61_source,landmarks_61_target)
  
  # Extract Values from modelled Images with Valdidation Points
  xy_disp_annual_4_val = raster::extract(r_stack[[3]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  xyz_disp_annual_4_val = raster::extract(r_stack[[4]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  
  
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
  
  print(paste("----------------------- Iteration",i,"completed-----------------------", sep = ""))
 }

Sys.time()
