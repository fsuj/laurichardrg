# PROCESSING BUNWARPJ RAW TRANSFORMATION FOR 3D DISPLACEMENT MAPPING
# WITH DIGITAL ELEVATION MODELS

# These functions are made for mapping 3D displacements of earth-surface 
# processes (e.g. creep) from 2D image registration based on elastic 
# deformations repsented by B-splines applied using digital elevation models.
# The resulting raw transformation file from the BUnwarpJ plugin in imagej
# is converted into the local coordinate system.

# Purpose: The free form deformation using b-splines can be performed using
#          derivatives of the digital elevation model (DEM) to improve
#          the registration results (e.g. a model of hillshade, slope, 
#          surface roughness etc...). To run the bunwarpj plugin in imagej
#          the input source and target images (of equal row & col) must
#          be converted to a non-spaial image format such as a png. To
#          apply the resulting displacement field (i.e. raw trasnformation)
#          of the image registration to our geogrpahic data, the displacement
#          field must be converted to a coordinate reference system (CRS).
#          The main 'bunwarpjDisplacementField' function takes care of this.
#
#          Note the raw transformation file from the bunwarpj function
#          describes the registration by giving the registration coordinates
#          for each pixel of the source image by referencing the new row 
#          and column position. This script is designed for use in a local
#          CRS.

# https://imagej.net/BUnwarpJ

## Euclidean Distance between two points
eucDist <- function(x1, x2, y1, y2) sqrt( (x1 - x2)^2 + (y1 - y2)^2 ) 

## imageJ Coordinate transform. image to CRS ###########################

imagejPointsToCRS <- function(file_name, r){
  
  # This function can convert the points saved from imagej
  # in image coordinates (col, row) with origin (1,1)
  # to the  coordinate reference system cooresponding
  # to the image that are being used.
  
  # file_name: file name of text file with point coordinates
  #            that was saved from imageJ
  # r: Raster* object
  
  # define origin in local coord from top left corner
  x_min <- xmin(r) # cell center at this point (not corner)
  y_max <- ymax(r) # cell center (not corner)
  x_scale <- res(r)[1]
  y_scale <- res(r)[2]
  
  # Read image coordinates (col, rows) from output of SIFT in imagej
  img_crd_pnts <- read.delim(file_name, header=FALSE)
  names( img_crd_pnts) <- c("col", "row")
  
  # Shift  the origin of the image from (1,1) to (0,0)
  img_crd_pnts <-  img_crd_pnts-1
  
  # Scale 1 pixel = 0.10 m
  img_crd_pnts$x_scale <-  img_crd_pnts$col * x_scale
  img_crd_pnts$y_scale <-  img_crd_pnts$row * y_scale
  
  # Calculate position in local CRS
  #   !note: x_scale/2 is set to make the CRS origin the cell corner
  #          instead of center
  
  img_crd_pnts$x_crs <- ( x_min - x_scale/2 + img_crd_pnts$x_scale)
  img_crd_pnts$y_crs <- ( y_max + y_scale/2 - img_crd_pnts$y_scale)
  
  return(data.frame(x=img_crd_pnts$x_crs, y=img_crd_pnts$y_crs))
}



## Matrix handling ########################################

matrixToRaster <- function(m, r_ref){
  require(raster)
  # Converts a matrix to a raster object using crs and 
  # extent of a reference raster
  
  #   !matrix nrow and ncol must match ref raster
  
  # m: matrix
  # r_ref: a raster object from {raster} used as 
  #        referece for CRS and and extent
  r_m <- raster(m, xmn=r_ref@extent@xmin,
                xmx=r_ref@extent@xmax,
                ymn=r_ref@extent@ymin,
                ymx=r_ref@extent@ymax,
                crs=crs(r_ref))
  return(r_m)
}


## Angle handling #########################################

aspectCalc <- function(x,y){
  # Calculate the direction component of the vector field
  # with a geographic reference system - 0 degrees as due
  # north and 90 degrees as due east (i.e. aspect)
  
  if(x==0 & y==0){
    asp <- 0
  } else if(x==0){
    asp <- 0
  } else if(x>0 & y>=0){
    asp <- abs(atan(x/y))*180/pi
  } else if(x>0 & y<=0) {
    asp <- 180 - abs(atan(x/y))*180/pi
  } else if(x<0 & y>=0) {
    asp <- 360 - abs(atan(x/y))*180/pi
  } else if (x<0 & y<=0){
    asp <- 180 + abs(atan(x/y))*180/pi
  }
  return(asp)
}


oppDir <- function(angle){
  # Calculate direction in a reverse orientation
  
  if(angle+180 >= 360){
    opp_asp <- angle - 180
  } else{
    opp_asp <- angle + 180
  }
  return(opp_asp)
}

## Import and process BUnwarpJ raw transformation file ####
bunwarpjDisplacementField <- function(tx_file, r_source, r_target, is_inverse = FALSE){
  
  # tx_file: the raw transformation file (.txt) from BUnwarpJ
  # r_source: the DEM raster filename used as source for image registration
  # target_file: the DEM raster filename used as target for image registration
  
  # Returns a data frame
  
  ## Determine no. of columns and rows from file header
  dimensions <- as.character(read.table(tx_file, header=FALSE, nrows=2)[,1])
  dimensions <- as.numeric(gsub("^.*\\=", "",dimensions))
  n_cols <- dimensions[1]
  n_rows <- dimensions[2]
  
  ## Read the X and Y displacements from the raw transformation file
  x_raw_disp <- read.table(tx_file, header=FALSE, skip=4, nrows=n_rows)
  y_raw_disp <- read.table(tx_file, header=FALSE, skip=6+n_rows, nrows=n_rows)
  
  # Convert to matrix
  m_x_disp <- data.matrix(x_raw_disp, rownames.force=FALSE)
  
  m_y_disp <- unname(data.matrix(y_raw_disp, rownames.force=FALSE))
  
  remove(x_raw_disp, y_raw_disp, dimensions)
  
  ## Convert cell to cell description of displacement (from transformation
  #    file) to relative displacements using row and column numbers
  for (i in 1:ncol(m_x_disp)){
    m_x_disp[,i] <- m_x_disp[, i] - (i)
  }
  
  for (i in 1:nrow(m_y_disp)){
    m_y_disp[i,] <- (m_y_disp[i, ] - (i)) * -1 # Made negative to reflect UTM scale
  }
  
  ## Load source and target rasters
  require(raster)
  #r_source <- raster(source_file)
  #r_target <- raster(target_file)
  
  # Get raster resolution
  resolution <- res(r_source)[1]
  
  ## Scale matrix displacement to spatial resolution or raster cells
  m_x_disp <- m_x_disp * resolution
  m_y_disp <- m_y_disp * resolution
  
  ## Convert matrix to raster with corresponding extent and CRS
  r_x_disp <- matrixToRaster(m_x_disp, r_source)
  r_y_disp <- matrixToRaster(m_y_disp, r_source)
  
  # Reverse direction of displacements
  if(is_inverse == FALSE){
  r_x_disp <- r_x_disp * -1
  r_y_disp <- r_y_disp * -1
  }
  
  ## Convert raster to SpatialPointsDataFrame
  # Using x displacements
  spdf <- rasterToPoints(r_x_disp, spatial=TRUE)
  
  # Add y displacements
  spdf$y_disp <- values(r_y_disp)
  names(spdf@data) <- c('x_disp', 'y_disp')
  
  ## Get elevation values from source DEM
  spdf$z_source <- extract(r_source, spdf)
  
  ## Convert SpatialPointsDataFrame to data.frame
  d <- as.data.frame(spdf)
  names(d) <- c('x_disp', 'y_disp', 'z_source', 'x_source', 'y_source') 
  
  # Compute target position from x and y displacement values
  d$x_target <- d$x_disp + d$x_source
  d$y_target <- d$y_disp + d$y_source
  
  # Compute 2D (xy) displacement magnitude
  d$xy_disp <- sqrt( (d$x_disp)^2 + (d$y_disp)^2 )
  
  # Compute 2D (xy) dispalcement vector angles relative to North
  d$aspect <- mapply(aspectCalc, x=d$x_disp, y=d$y_disp)
  
  # Get elevation values from target DEM at target positions
  # i.e. elevations at the corresponding points from the displaced
  # displaced source positions
  
  xy <- data.frame(x=d$x_target, y=d$y_target)
  d$z_target <- extract(r_target, xy)

  
  # Compute z and xyz (3D) displacements
  d$z_disp <- d$z_target-d$z_source
  d$xyz_disp <- sqrt( (d$xy_disp)^2 + (d$z_disp)^2 )
  
  # Compute slope angle (degrees) of z displacement
  d$slope <- atan( d$z_disp/d$xy_disp )*180/pi
   
  
  ## Compute the values of the xyz of the registered DEM
  #https://brilliant.org/wiki/3d-coordinate-geometry-distance/
  #d$x_register <- d$x_source + sin(d$aspect*pi/180) * d$xy_disp
  #d$y_register <- d$y_source + cos(d$aspect*pi/180) * d$xy_disp
  #d$z_register <- d$z_source + tan(d$slope*pi/180) * d$xy_disp
  
  ## Estimate the performance - error (m)
  #d$z_err <- d$z_target - d$z_register
  
  ## Clean d for a nice output data.frame
  d_out <- data.frame(
    x_source=d$x_source,
    y_source=d$y_source,
    z_source=d$z_source,
    x_target=d$x_target,
    y_target=d$y_target,
    z_target=d$z_target,
    x_disp=d$x_disp,
    y_disp=d$y_disp,
    z_disp=d$z_disp,
    xy_disp = d$xy_disp, #aka total displacement for imcorr
    xyz_disp= d$xyz_disp,
    aspect=d$aspect,
    slope=d$slope
    #x_reg=d$x_register,
    #y_reg=d$y_register,
    #z_reg=d$z_register,
    #z_err=d$z_err
  )
  
  ## Description of the return - a data frame with the following columns
  # 	x_source & y_source: x & y coordinates from a CRS for the registration
  #							 source image
  #	    		   z_source: the elevation value (intensity) corresponding the source
  #				             locations
  # 	x_target & y_target: x & y coordinates from a CRS representing the
  #							 x & y transformation for every location in the
  #							 source image
  #                z_target: the elevation value (intensity) corresponding the target
  #				             locations (from the target image /elevation model)
  #                  x_disp: displacement of registered image (from source image)
  #							 in the x direction
  #                  y_disp: displacement of registered image (from source image)
  #							 in the y direction
  #                  z_disp: displacement of registered image (from source image)
  #							 in the z direction (z_target - z_source)
  #                 xy_disp: magnitude of the displacement in 2D (xy)
  #				   xyz_disp: magnitude of the displacement in 3D (xyz)
  #                  aspect: geographic orientation of vector direction with 0 degrees
  #							 being due north and 90 degrees being due east
  #                   slope: angle in degrees of displacement in z direction
  #    x_reg, y_reg & z_reg: the x,y & z coordinates of the registered image 
  #                          where the z value has been adjusted with trignometry
  #							 estimations using xy_disp, z_disp and slope angle.
  							 

  remove(d)
  return(d_out)

}


