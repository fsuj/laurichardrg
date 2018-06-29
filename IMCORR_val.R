source("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Scripts/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

library(raster)
library(Metrics)
library(caret)
library(rgdal)
library(RStoolbox)

#Set wdir
maindir="D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data/ImageJ"
raw = "raw"
imcorr_folder = "IMCORR"


# Load Data
# Open idw-aspectcorr_DEM_VEC_xxx in qGIS -> speichern als tiff -> export (Remove - with _ in filename and _IJ, xy or xyz at the beginning)

imcorr_files =  list.files(sprintf("%s/%s/",maindir,imcorr_folder),pattern="*IJ.tif$", full.names = TRUE)

imcorr_r_stack = raster::stack(imcorr_files) # 1:3 = XY, 4:6 = XYZ


val_df = read.csv(sprintf("%s/Laurichard_2014_L93_extract_velocity.csv",maindir), header = TRUE, sep = ",")
landmark_val = sprintf("%s/Landmarks_val_new.txt",maindir)
contour_shp = readOGR(dsn = maindir, layer = "laurichard_rg_contour")

source_dem = raster(paste(maindir,"/","DEM_L93_160812_LIDAR_IJ.tif", sep = ""))
target_dem = raster(paste(maindir,"/","DEM_L93_051017_SFMMVS_IJ.tif", sep = ""))


imcorr_r_stack = projectRaster(from = imcorr_r_stack, to = source_dem,method = "bilinear")

#####################################################################################################################

# 4.2) Convert Validation Data to Point Shapefile / spdf

val_spdf_source = make_spdf(val_df[,c(1,2)],val_df,imcorr_r_stack[[1]])
val_spdf_source = crop(val_spdf_source, extent(imcorr_r_stack[[1]])) # Get Points within Raster

val_spdf_target = make_spdf(val_df[,c(4,5)],val_df,imcorr_r_stack[[1]])
val_spdf_target = crop(val_spdf_target, extent(imcorr_r_stack[[1]])) # Get Points within Raster       ####### CHECK THIS STUPID EXTNET?!

#########
landmark_df_val = imagejPointsToCRS2(landmark_val, imcorr_r_stack[[1]])

landmarks_60_val_source = make_spdf(landmark_df_val[,c(2,3)],landmark_df_val,imcorr_r_stack[[1]])
landmarks_60_val_source[,1]
landmarks_60_val_target = make_spdf(landmark_df_val[,c(4,5)],landmark_df_val,imcorr_r_stack[[1]])
landmarks_60_val_target[,1]


# Calculate Dispalcements XY /XYZ
landmark_df_val["XY_Displacement"] = eucDist(landmark_df_val$x_crs_Source,landmark_df_val$x_crs_Target,
                                             landmark_df_val$y_crs_Source,landmark_df_val$y_crs_Target)
landmark_df_val["XY_Displacement_year"] = landmark_df_val["XY_Displacement"] / 6
landmark_df_val["z_source"] = extract(source_dem[[1]],landmarks_60_val_source)
landmark_df_val["z_target"] = extract(target_dem[[1]],landmarks_60_val_target)
landmark_df_val["Z_Displacement"] = landmark_df_val$z_target - landmark_df_val$z_source
landmark_df_val["XYZ_Displacement"] = sqrt( (landmark_df_val$XY_Displacement)^2 + (landmark_df_val$Z_Displacement)^2 )
landmark_df_val["XYZ_Displacement_year"] = landmark_df_val["XYZ_Displacement"] / 6

landmarks_60_val_source = make_spdf(landmark_df_val[,c(2,3)],landmark_df_val,imcorr_r_stack[[1]])
landmarks_60_val_source[,1]
landmarks_60_val_target = make_spdf(landmark_df_val[,c(4,5)],landmark_df_val,imcorr_r_stack[[1]])
landmarks_60_val_target[,1]

#########
# calculate Annual Displacement from IMCORR displacements

imcorr_r_stack = imcorr_r_stack / 6

# mask Stack with contour

r_name = names(imcorr_r_stack)


library(foreach)
library(doParallel) 
cl<-makeCluster(4)
registerDoParallel(cl)

Sys.time()
mask_ls = foreach(i = 1:length(names(imcorr_r_stack)),.packages = "raster",.inorder = TRUE) %dopar% {
  
  mask_raster_layer(imcorr_r_stack[[i]],contour_shp,r_name[i])
  
}
Sys.time()
stopCluster(cl)
gc()

imcorr_r_stack = raster::stack(mask_ls)

##########
# extract Data From stack 

##################################### Create Summary of filterd and masked RasterStack (Write to csv) #############
#####################################               Plot pdf with valdidation points                  #############
imcorr_r_stack = raster::brick(imcorr_r_stack)


# Summary
summary_raster_Stack = round(summary(imcorr_r_stack, na.rm = T), digits = 4)
summary_raster_Stack = data.frame(unclass(summary_raster_Stack), check.names = FALSE, stringsAsFactors = FALSE)

# IQR
IQR_final = calc_IQR_Raster(imcorr_r_stack)

# SD
SD_final= calc_SD_Raster(imcorr_r_stack)

summary_raster_Stack["IQR",] = IQR_final
summary_raster_Stack["SD",] = SD_final

write.csv(summary_raster_Stack,paste(maindir,"/",imcorr_folder,"/",names(imcorr_r_stack)[3],"_Raster_Stats.csv", sep =""))


# Write Raster
for (j in 1:length(names(imcorr_r_stack))){
  
  writeRaster(imcorr_r_stack[[j]],filename = paste(maindir,"/",imcorr_folder,"/",names(imcorr_r_stack)[j],"_masked",sep=""), format="GTiff", overwrite=TRUE)
  print(paste("Write_Raster",j,sep=""))
  
}

#################################

################################################## Extracting Valdidation Points #################################################

# Extract Values from modelled Images with Valdidation Points
my_list_bodin_XY = list()
my_list_bodin_XYZ = list()

for ( i in 1:length(names(imcorr_r_stack))){
  
  if (i <=3) {
  
  xy_disp_annual_4_val = raster::extract(imcorr_r_stack[[i]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  name = names(imcorr_r_stack)[i]
  my_list_bodin_XY[[name]] = xy_disp_annual_4_val
  } else {
  
  xyz_disp_annual_4_val = raster::extract(imcorr_r_stack[[i]],val_spdf_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  name = names(imcorr_r_stack)[i]
  my_list_bodin_XYZ[[name]] = xyz_disp_annual_4_val
  }
  
}

#make DataFrame
bodinXY = do.call("cbind",my_list_bodin_XY)
bodinXY_df = as.data.frame(bodinXY[c("ID","X2015_2012_1","xy_idw_aspectcorr_DEM_DISP_VEC_128_64_01_IJ","xy_idw_aspectcorr_DEM_DISP_VEC_256_128_01_IJ","xy_idw_aspectcorr_DEM_DISP_VEC_256_64_01_IJ")])

bodinXYZ = do.call("cbind",my_list_bodin_XYZ)
bodinXYZ_df = as.data.frame(bodinXYZ[c("ID","X2015_2012_1","xyz_idw_aspectcorr_DEM_DISP_VEC_128_64_01_IJ","xyz_idw_aspectcorr_DEM_DISP_VEC_256_128_01_IJ","xyz_idw_aspectcorr_DEM_DISP_VEC_256_64_01_IJ")])


########
my_list_manual_XY = list()
my_list_manual_XYZ = list()

for ( i in 1:length(names(imcorr_r_stack))){
  
  if(i <= 3){
  xy_disp_annual_4_val_user = raster::extract(imcorr_r_stack[[i]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  name = names(imcorr_r_stack)[i]
  my_list_manual_XY[[name]] = xy_disp_annual_4_val_user
  } else{
  
  xyz_disp_annual_4_val_user = raster::extract(imcorr_r_stack[[i]],landmarks_60_val_source, method = "simple", buffer = NULL, na.rm=TRUE, sp=TRUE)
  name = names(imcorr_r_stack)[i]
  my_list_manual_XYZ[[name]] = xyz_disp_annual_4_val_user
  }
  
}

manualXY = do.call("cbind",my_list_manual_XY)
manualXY_df = as.data.frame(manualXY[c("Index","XY_Displacement_year","XYZ_Displacement_year","xy_idw_aspectcorr_DEM_DISP_VEC_128_64_01_IJ","xy_idw_aspectcorr_DEM_DISP_VEC_256_128_01_IJ","xy_idw_aspectcorr_DEM_DISP_VEC_256_64_01_IJ")])


manualXYZ = do.call("cbind",my_list_manual_XYZ)
manualXYZ_df = as.data.frame(manualXYZ[c("Index","XY_Displacement_year","XYZ_Displacement_year","xyz_idw_aspectcorr_DEM_DISP_VEC_128_64_01_IJ","xyz_idw_aspectcorr_DEM_DISP_VEC_256_128_01_IJ","xyz_idw_aspectcorr_DEM_DISP_VEC_256_64_01_IJ")])

################

########################################### Create Validatin DF for Bodin Data (3D displacement)####################################
########################################### Create Validiaton DF for Expert based Approach      ####################################

############ BODIN et al. 2017 XY
full_val_df_XY = bodinXY_df
full_val_df_XY["X"] = NULL
full_val_df_XY["Y"] = NULL


# Calculate Displacement Error (dm-dp)
full_val_df_XY[paste(colnames(full_val_df_XY)[3],"Error",sep="")] = full_val_df_XY[colnames(full_val_df_XY[2])] - full_val_df_XY[colnames(full_val_df_XY[3])] ## VORSIcHT (Bodin 3D?!)
full_val_df_XY[paste(colnames(full_val_df_XY)[4],"Error",sep="")] = full_val_df_XY[colnames(full_val_df_XY[2])] - full_val_df_XY[colnames(full_val_df_XY[4])] ## VORSIcHT (Bodin 3D?!)
full_val_df_XY[paste(colnames(full_val_df_XY)[5],"Error",sep="")] = full_val_df_XY[colnames(full_val_df_XY[2])] - full_val_df_XY[colnames(full_val_df_XY[5])] ## VORSIcHT (Bodin 3D?!)


# Plot Scatterplot
my_scatterplot_IMCORR(paste(maindir,"/",imcorr_folder,"/",colnames(full_val_df_XY)[3],"_Scatterplot_Bodin_XY.pdf", sep =""), full_val_df_XY,
               colnames(full_val_df_XY)[3],colnames(full_val_df_XY)[2],
               colnames(full_val_df_XY)[4],colnames(full_val_df_XY)[2],
               colnames(full_val_df_XY)[5],colnames(full_val_df_XY)[2],
               "XY Displacement IMCORR [m/year]","XY Displacement Bodin et al. 2017 [m/year]")
#################################################################################################
#### Calculate Statistics 
# R²
# RMSE
# Summary
# IQR
# Error Displacement (inkl. Summary)
# SD

final_stats = my_statistic_calc_imcorr(full_val_df_XY,colnames(full_val_df_XY)[2],
                                       colnames(full_val_df_XY)[3],colnames(full_val_df_XY)[4],
                                       colnames(full_val_df_XY)[5])  

# Write Statistics to csv
write.csv(final_stats,paste(maindir,"/",imcorr_folder,"/",colnames(full_val_df_XY)[3],"_Stats_Bodin_XY.csv", sep =""))
############
############################# Bodin et al. 2017 XYZ

full_val_df_XYZ = bodinXYZ_df
full_val_df_XYZ["X"] = NULL
full_val_df_XYZ["Y"] = NULL


# Calculate Displacement Error (dm-dp)
full_val_df_XYZ[paste(colnames(full_val_df_XYZ)[3],"Error",sep="")] = full_val_df_XYZ[colnames(full_val_df_XYZ[2])] - full_val_df_XYZ[colnames(full_val_df_XYZ[3])] ## VORSIcHT (Bodin 3D?!)
full_val_df_XYZ[paste(colnames(full_val_df_XYZ)[4],"Error",sep="")] = full_val_df_XYZ[colnames(full_val_df_XYZ[2])] - full_val_df_XYZ[colnames(full_val_df_XYZ[4])] ## VORSIcHT (Bodin 3D?!)
full_val_df_XYZ[paste(colnames(full_val_df_XYZ)[5],"Error",sep="")] = full_val_df_XYZ[colnames(full_val_df_XYZ[2])] - full_val_df_XYZ[colnames(full_val_df_XYZ[5])] ## VORSIcHT (Bodin 3D?!)


# Plot Scatterplot
my_scatterplot_IMCORR(paste(maindir,"/",imcorr_folder,"/",colnames(full_val_df_XYZ)[3],"_Scatterplot_Bodin_XYZ.pdf", sep =""), full_val_df_XYZ,
                      colnames(full_val_df_XYZ)[3],colnames(full_val_df_XYZ)[2],
                      colnames(full_val_df_XYZ)[4],colnames(full_val_df_XYZ)[2],
                      colnames(full_val_df_XYZ)[5],colnames(full_val_df_XYZ)[2],
                      "XYZ Displacement IMCORR [m/year]","XYZ Displacement Bodin et al. 2017 [m/year]")
#################################################################################################
#### Calculate Statistics 
# R²
# RMSE
# Summary
# IQR
# Error Displacement (inkl. Summary)
# SD


final_stats = my_statistic_calc_imcorr(full_val_df_XYZ,colnames(full_val_df_XYZ)[2],
                                       colnames(full_val_df_XYZ)[3],colnames(full_val_df_XYZ)[4],
                                       colnames(full_val_df_XYZ)[5])  

# Write Statistics to csv
write.csv(final_stats,paste(maindir,"/",imcorr_folder,"/",colnames(full_val_df_XYZ)[3],"_Stats_Bodin_XYZ.csv", sep =""))
#######################################################################################################

############ Manuell Approach XY

# Calclulate Validiation vor Manual Approch

full_val_manuell_df_XY = manualXY_df
full_val_manuell_df_XY["x_crs_Source"] = NULL
full_val_manuell_df_XY["y_crs_Source"] = NULL

# Calculate Displacement Error 2D / 3D

full_val_manuell_df_XY[paste(colnames(full_val_manuell_df_XY)[4],"Error",sep="")] = full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[2])] - full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[4])] 
full_val_manuell_df_XY[paste(colnames(full_val_manuell_df_XY)[5],"Error",sep="")] = full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[2])] - full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[5])] 
full_val_manuell_df_XY[paste(colnames(full_val_manuell_df_XY)[6],"Error",sep="")] = full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[2])] - full_val_manuell_df_XY[colnames(full_val_manuell_df_XY[6])] 



my_scatterplot_IMCORR(paste(maindir,"/",imcorr_folder,"/",colnames(full_val_manuell_df_XY)[4],"_Scatterplot_manual_XY.pdf", sep =""), full_val_manuell_df_XY,
                      colnames(full_val_manuell_df_XY)[4],colnames(full_val_manuell_df_XY)[2],
                      colnames(full_val_manuell_df_XY)[5],colnames(full_val_manuell_df_XY)[2],
                      colnames(full_val_manuell_df_XY)[6],colnames(full_val_manuell_df_XY)[2],
                      "XY Displacement IMCORR [m/year]","XY Displacement Expert based [m/year]")



final_stats_manuell = my_statistic_calc_expert_IMCORR(full_val_manuell_df_XY,colnames(full_val_manuell_df_XY)[2],
                                                      colnames(full_val_manuell_df_XY)[4],colnames(full_val_manuell_df_XY)[5],
                                                      colnames(full_val_manuell_df_XY)[6])  
# Write Statistics to csv
write.csv(final_stats_manuell,paste(maindir,"/",imcorr_folder,"/",colnames(full_val_manuell_df_XY)[4],"_Stats_manuell_XY.csv", sep =""))
######################

## Manuell Approach XYZ#
# Calclulate Validiation vor Manual Approch

full_val_manuell_df_XYZ = manualXYZ_df
full_val_manuell_df_XYZ["x_crs_Source"] = NULL
full_val_manuell_df_XYZ["y_crs_Source"] = NULL

# Calculate Displacement Error 2D / 3D

full_val_manuell_df_XYZ[paste(colnames(full_val_manuell_df_XYZ)[4],"Error",sep="")] = full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[3])] - full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[4])] 
full_val_manuell_df_XYZ[paste(colnames(full_val_manuell_df_XYZ)[5],"Error",sep="")] = full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[3])] - full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[5])] 
full_val_manuell_df_XYZ[paste(colnames(full_val_manuell_df_XYZ)[6],"Error",sep="")] = full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[3])] - full_val_manuell_df_XYZ[colnames(full_val_manuell_df_XYZ[6])] 



my_scatterplot_IMCORR(paste(maindir,"/",imcorr_folder,"/",colnames(full_val_manuell_df_XYZ)[4],"_Scatterplot_manual_XYZ.pdf", sep =""), full_val_manuell_df_XYZ,
                      colnames(full_val_manuell_df_XYZ)[4],colnames(full_val_manuell_df_XYZ)[3],
                      colnames(full_val_manuell_df_XYZ)[5],colnames(full_val_manuell_df_XYZ)[3],
                      colnames(full_val_manuell_df_XYZ)[6],colnames(full_val_manuell_df_XYZ)[3],
                      "XYZ Displacement IMCORR [m/year]","XYZ Displacement Expert based [m/year]")



final_stats_manuell = my_statistic_calc_expert_IMCORR(full_val_manuell_df_XYZ,colnames(full_val_manuell_df_XYZ)[3],
                                                      colnames(full_val_manuell_df_XYZ)[4],colnames(full_val_manuell_df_XYZ)[5],
                                                      colnames(full_val_manuell_df_XYZ)[6])  
# Write Statistics to csv
write.csv(final_stats_manuell,paste(maindir,"/",imcorr_folder,"/",colnames(full_val_manuell_df_XYZ)[4],"_Stats_manuell_XYZ.csv", sep =""))

