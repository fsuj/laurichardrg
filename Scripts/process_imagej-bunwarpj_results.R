# Scritp to read RAW transformations from imagej-bunwarpj plugin

source("D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Scripts/r/source/process_BUnwarpJ_2D_image_registrationV02.R")

library(raster)
library(Metrics)
library(caret)
library(rgdal)

#Set wdir
maindir="D:/Studium/Jena/M.Sc. Geoinformatik/Geo411_Landschaftsmanagment_und_Fernerkundung/Geo_411_Jason_practical/Data"

# Load Data

r_source <- raster(sprintf("%s/DEM-L93-16Aug2012-LIDAR.tif",maindir))
r_target <- raster(sprintf("%s/DEM-L93-05Oct2017-SFMMVS.tif",maindir))

points_df = imagejPointsToCRS2(sprintf("%s/Landmarks_auto.txt",maindir), r_source)
val_df = read.csv(sprintf("%s/Laurichard_2014_L93_extract_velocity.csv",maindir), header = TRUE, sep = ",")

######################################## Convert Landmarks to Shapefile

# Extract Coords from Pixel Coords
points_df = imagejPointsToCRS2(sprintf("%s/Landmarks_auto.txt",maindir), r_source)

# Make spdf
points_source  = make_spdf(points_df[,c(2,3)],points_df,r_source)
points_source[,1]
points_target  = make_spdf(points_df[,c(4,5)],points_df,r_source)
points_target[,1]

# Save to Shapefile

writeOGR(points_source,dsn = maindir,layer ="Landmarks_Points_Source", driver="ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(points_target,dsn = maindir,layer ="Landmarks_Points_Target", driver="ESRI Shapefile", overwrite_layer=TRUE )

######################################## Convert Validation Data to spdf

val_spdf = make_spdf(val_df[,c(1,2)],val_df,r_source)
val_spdf = crop(val_spdf, extent(r_source)) # Get Points within Raster

######################################## Process Data BunwrapJ

# d_tx <- bunwarpjDisplacementField(tx_file =  sprintf("%s/Hillshade-L93-16Aug2012-LIDAR_direct_transf_auto_raw.txt",maindir),
#                                r_source, r_target, is_inverse=FALSE)
d_tx <- bunwarpjDisplacementField(tx_file =  sprintf("%s/Hillshade-L93-16Aug2012-LIDAR_direct_transf_rm_raw.txt",maindir),
                                  r_source, r_target, is_inverse=FALSE)
saveRDS(d_tx, file = "d_tx_direct_transf_auto.RDS")

d_tx = readRDS(file = "d_tx_direct_transf_auto.RDS")

# d_tx <- bunwarpjDisplacementField(tx_file =  sprintf("%s/Hillshade-L93-05Oct2017-SFMMVS_inverse_transf_auto_raw.txt",maindir),
#                                   r_source, r_target, is_inverse=TRUE)
#saveRDS(d_tx, file = "d_tx_inverse_transf_auto.RDS")
#d_tx = readRDS(file = "d_tx_inverse_transf_auto.RDS")

summary(d_tx)
head(d_tx)

# Calculate Annual Displacment in m/year (XY,XYZ) ???????????

d_tx$xy_disp_annual = d_tx$xy_disp / 6
d_tx$xyz_disp_annual = d_tx$xyz_disp / 6

#Write RasterStack with XY,XYZ (2012-2017),Annual XY,XYZ (MEAN, m/Year)

r_xy <- setValues(r_source, d_tx$xy_disp)
r_xyz = setValues(r_source, d_tx$xyz_disp)
r_xy_annual = setValues(r_source, d_tx$xy_disp_annual)
r_xyz_annual = setValues(r_source, d_tx$xyz_disp_annual)

r_list = list(r_xy,r_xyz,r_xy_annual,r_xyz_annual)
r_name = c("xyDisp_demDirect_auto.tif","xyzDisp_demDirect_auto.tif","xyDisp_demDirect_auto_anual.tif","xyzDisp_demDirect_auto_anual.tif")
r_stack = raster::brick(r_list)
names(r_stack) = r_name

# Write r_stack to Files

for (i in 1:length(names(r_stack))){
  
  writeRaster(r_stack[[i]],filename = names(r_stack)[i], format="GTiff", overwrite=TRUE)
  
}

#################################### Plot Data

my_plot_function("xy_z_displacments_dem_direct.pdf",r_stack,val_spdf,points_source,points_target)

#################################### Extract Values from  Anual Raster at Valdidation Points
###########
# XY Displacement Validation or XYZ Displacement Validation?!?!?!?!?!?
############
####################################

xy_disp_annual_4_val = raster::extract(r_stack[[3]],val_spdf, method = "simple", buffer = NULL, na.rm=FALSE, sp=TRUE)
xyz_disp_annual_4_val = raster::extract(r_stack[[4]],val_spdf, method = "simple", buffer = NULL, na.rm=FALSE, sp=TRUE)

##################################### Combine all DATA
head(val_df)
head(val_spdf)
head(xy_disp_annual_4_val)
head(xyz_disp_annual_4_val)

# Create Validation Data Frame

full_val_spdf = cbind(val_spdf[,c(3,8)],xy_disp_annual_4_val[,12],xyz_disp_annual_4_val[,12])
names(full_val_spdf)= c("ID","X2015_2012_1","xy_annual","xyz_annual")

# Calculate Displacment Error (dm-dp) (xy,xyz)
full_val_df$xy_annual_error = full_val_df$X2015_2012_1 - full_val_df$xy_annual
full_val_df$xyz_annual_error = full_val_df$X2015_2012_1 - full_val_df$xyz_annual

full_val_df = as.data.frame(full_val_spdf)

##################################### Plot Full_val_df
dev.off()

par(mfrow = c(1,2))
plot(full_val_df$xy_annual,full_val_df$X2015_2012_1, main ="Scatterplot XY Displacement",
     xlab = "Displacement BUnwrap [m/year]", ylab = "Displacment Bodin et al. 2017 [m/year]",
     xlim = c(0,2.5), ylim = c(0,2.5))
abline(lm(full_val_df$X2015_2012_1 ~ full_val_df$xy_annual))
par(new = F)
plot(full_val_df$xyz_annual,full_val_df$X2015_2012_1, main ="Scatterplot XYZ Displacement",
     xlab = "Displacement BUnwrap [m/year]", ylab = "Displacment Bodin et al. 2017 [m/year]",
     xlim = c(0,2.5), ylim = c(0,2.5))
abline(lm(full_val_df$X2015_2012_1 ~ full_val_df$xyz_annual))

##################################### Calculate Statistics 
# R²
# RMSE
# Summary
# IQR
# Error Displacement (inkl. Summary)

## Calculate Statistics

# RMSE,Rsquared, MAE
xy_R_squared = caret::postResample(full_val_df$xy_annual,full_val_df$X2015_2012_1)
xyz_R_squared = caret::postResample(full_val_df$xyz_annual,full_val_df$X2015_2012_1)

# Summary
summary_df = summary(full_val_df)

# IQR (only for single Colum)
IQR_val = round(IQR(full_val_df$X2015_2012_1, na.rm =T, type = 7), digits = 4)
IQR_xy  = round(IQR(full_val_df$xy_annual, na.rm =T, type = 7), digits = 4)
IQR_xyz = round(IQR(full_val_df$xyz_annual, na.rm = T, type = 7), digits =4)
IQR_xy_error = round(IQR(full_val_df$xy_annual_error, na.rm =T, type = 7), digits = 4)
IQR_xyz_error = round(IQR(full_val_df$xyz_annual_error, na.rm =T, type = 7), digits = 4)


# Make DF with Statistics
# SUmmary and IQR
df_1 = data.frame(unclass(summary_df), check.names = FALSE, stringsAsFactors = FALSE)
df_1 = df_1[,c(2,3,7,4,8)]

# IQR
df_IQR_val = as.data.frame(IQR_val)
df_IQR_val$IQR_val_combine = c(paste0(colnames(df_IQR_val),": ",df_IQR_val[1,]))

df_IQR_xy = as.data.frame(IQR_xy)
df_IQR_xy$IQR_val_combine = c(paste0(colnames(df_IQR_xy),": ",df_IQR_xy[1,]))

df_IQR_xyz = as.data.frame(IQR_xyz)
df_IQR_xyz$IQR_val_combine = c(paste0(colnames(df_IQR_xyz),": ",df_IQR_xyz[1,]))

df_IQR_xy_error = as.data.frame(IQR_xy_error)
df_IQR_xy_error$IQR_val_combine = c(paste0(colnames(df_IQR_xy_error),": ",df_IQR_xy_error[1,]))

df_IQR_xyz_error = as.data.frame(IQR_xyz_error)
df_IQR_xyz_error$IQR_val_combine = c(paste0(colnames(df_IQR_xyz_error),": ",df_IQR_xyz_error[1,]))

df_1[7,] = c(df_IQR_val[1,2],df_IQR_xy[1,2],df_IQR_xy_error[1,2],df_IQR_xyz[1,2], df_IQR_xyz_error[1,2])


df_2 = as.data.frame(xy_R_squared)
df_2 = round(df_2, digits = 4)
df_2$xy_R_squared_combine = c(paste0(row.names(df_2),": ",df_2[,1]))
df_2[nrow(df_2)+4,] <- NA
df_2 = df_2[,2]

df_3 = as.data.frame(xyz_R_squared)
df_3 = round(df_3, digits = 4)
df_3$xyz_R_squared_combine = c(paste0(row.names(df_3),": ",df_3[,1]))
df_3[nrow(df_3)+4,] <- NA
df_3 = df_3[,2]

final_stats = cbind(df_1,df_2,df_3)
colnames(final_stats) = c("X2015_2012_1","xy_annual","xy_annual_error","xyz_annual","xyz_annual_error","xy_R_squared","xyz_R_squared") #,"5","6","7")
final_stats
