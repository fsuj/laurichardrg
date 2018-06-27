REM @ECHO OFF

REM D:
REM D:\Software\saga-6.1.0_x64>

REM *********************************************************************
REM Calculate surface velocities from Digital elevation models (slope)
REM *********************************************************************

set PATH=D:\Software\saga-6.1.0_x64


REM SEARCH CHIP 128x128 REF CHIP 32x32 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Movement\slope_LIDAR_DTM_Aug2012.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_SFMMVS_DTM_Oct2017.tif" -DTM_1="D:\WorkingProject\Laurichard\Movement\slope_LIDAR_DTM_Aug2012.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\slope_SFMMVS_DTM_Oct2017.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts128-32-50cm.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec128-32-50cm.shp" -SEARCH_CHIPSIZE=3 -REF_CHIPSIZE=1 -GRID_SPACING=0.1



REM *********************************************************************
REM REFERFENCE CHIP 64 x 64
REM *********************************************************************

REM SEARCH CHIP 128x128 REF CHIP 64x64 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Movement\slope_LIDAR_DTM_Aug2012.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_SFMMVS_DTM_Oct2017.tif" -DTM_1="D:\WorkingProject\Laurichard\Movement\slope_LIDAR_DTM_Aug2012.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\slope_SFMMVS_DTM_Oct2017.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts128-64-50cm.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec128-64-50cm.shp" -SEARCH_CHIPSIZE=3 -REF_CHIPSIZE=2 -GRID_SPACING=0.1




