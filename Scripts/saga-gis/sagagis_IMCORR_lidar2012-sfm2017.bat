REM @ECHO OFF

REM D:
REM D:\Software\saga-6.1.0_x64>

REM *********************************************************************
REM Calculate surface velocities from Digital elevation models (slope)
REM *********************************************************************

set PATH=D:\Software\saga-6.1.0_x64

REM *********************************************************************
REM REFERFENCE CHIP 16 x 16
REM *********************************************************************

REM SEARCH CHIP 32x32 REF CHIP 16x16 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts32-16-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec32-16-10.shp" -SEARCH_CHIPSIZE=1 -REF_CHIPSIZE=0 -GRID_SPACING=5

REM SEARCH CHIP 64x64 REF CHIP 16x16 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts64-16-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec64-16-10.shp" -SEARCH_CHIPSIZE=2 -REF_CHIPSIZE=0 -GRID_SPACING=5

REM SEARCH CHIP 128x128 REF CHIP 16x16 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts128-16-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec128-16-10.shp" -SEARCH_CHIPSIZE=3 -REF_CHIPSIZE=0 -GRID_SPACING=5

REM SEARCH CHIP 128x128 REF CHIP 16x16 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts265-16-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec256-16-10.shp" -SEARCH_CHIPSIZE=4 -REF_CHIPSIZE=0 -GRID_SPACING=5

REM *********************************************************************
REM REFERFENCE CHIP 32 x 32
REM *********************************************************************

REM SEARCH CHIP 64x64 REF CHIP 32x32 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts64-32-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec64-32-10.shp" -SEARCH_CHIPSIZE=2 -REF_CHIPSIZE=1 -GRID_SPACING=5

REM SEARCH CHIP 128x128 REF CHIP 32x32 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts128-32-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec128-32-10.shp" -SEARCH_CHIPSIZE=3 -REF_CHIPSIZE=1 -GRID_SPACING=5

REM SEARCH CHIP 256x256 REF CHIP 32x32 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts256-32-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec256-32-10.shp" -SEARCH_CHIPSIZE=4 -REF_CHIPSIZE=1 -GRID_SPACING=5

REM *********************************************************************
REM REFERFENCE CHIP 64 x 64
REM *********************************************************************

REM SEARCH CHIP 128x128 REF CHIP 64x64 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts128-64-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec128-64-10.shp" -SEARCH_CHIPSIZE=3 -REF_CHIPSIZE=2 -GRID_SPACING=5

REM SEARCH CHIP 256x256 REF CHIP 64x64 GRID SPACING 5
saga_cmd grid_analysis 19 -GRID_1="D:\WorkingProject\Laurichard\Data\SLOPE_MNT50cm_LAURICHARD.tif" -GRID_2="D:\WorkingProject\Laurichard\Movement\slope_1421_20171005_50cm-match.tif" -DTM_1="D:\WorkingProject\Laurichard\Data\MNT50cm_LAURICHARD.tif" -DTM_2="D:\WorkingProject\Laurichard\Movement\dem_1421_20171005_50cm-match.tif" -CORRPOINTS="D:\WorkingProject\Laurichard\Movement\Workspace\crrPnts256-64-10.shp" -CORRLINES="D:\WorkingProject\Laurichard\Movement\Workspace\crrVec256-64-10.shp" -SEARCH_CHIPSIZE=4 -REF_CHIPSIZE=2 -GRID_SPACING=5


