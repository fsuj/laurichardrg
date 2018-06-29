## Load Data ##########################################################
setwd("D:/WorkingProject/Laurichard/Movement/Workspace/EditShape")

library(rgdal)

imcorrOut <- readOGR(".", "ed-crrPnts128-32-01")

## Estimate displacement from Aug 2012 to  Oct 2017 #######################################

# Create new dataframe
d <- data.frame(total_disp=imcorrOut$DISP)
d$slope <- imcorrOut$SLOPE
d$aspect <- imcorrOut$ASPECT
d$x_real <- imcorrOut$REALX
d$y_real <- imcorrOut$REALY
d$z_real <- imcorrOut$REALZ
d$x_target <- imcorrOut$XTARG 
d$y_target <- imcorrOut$YTARG
d$z_target <- imcorrOut$ZTARG
d$disp_real <- imcorrOut$DISP_REAL