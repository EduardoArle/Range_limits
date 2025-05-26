#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250504_GAMs/Models'

#load models
setwd(wd_models)
models <- lapply(list.files(), readRDS)
names(models) <- list.files()

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#plot A.i (minPPT relative polewardness)
plot.gam(models$minPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.04), xlim = c(0, 1)) 

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)


#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("SHAP value", side = 2, line = 6.5, cex = 2.5)

#save 800

#plot A.ii (meanPPT relative polewardness)
plot.gam(models$meanPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#dfc27d50',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.04), xlim = c(0, 1))

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#save 800

#plot A.iii (maxPPT relative polewardness)
plot.gam(models$maxPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#01665e30',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.04), xlim = c(0, 1))

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#save 800



#plot B.i (minPPT distance to range edge)
plot.gam(models$minPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250)) 

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  
mtext("SHAP value", side = 2, line = 6.5, cex = 2.5)

#save 800

#plot B.ii (meanPPT distance to range edge)
plot.gam(models$meanPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#dfc27d50',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250))

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#save 800

#plot B.iii (maxPPT distance to range edge)
plot.gam(models$maxPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#01665e30',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250))

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#save 800


