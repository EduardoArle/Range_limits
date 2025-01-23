#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250118_GAMs/Models'

#load models
setwd(wd_models)
models <- lapply(list.files(), readRDS)
names(models) <- list.files()

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#plot A.i (minT centralness)
plot.gam(models$minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.15, 0.05), xlim = c(0, 800)) 

#add axes
axis(1, pos = -0.15, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  
mtext("SHAP value", side = 2, line = 6.5, cex = 2.5)

#save 800

#plot A.ii (meanT centralness)
plot.gam(models$meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.15, 0.05), xlim = c(0, 800))

#add axes
axis(1, pos = -0.15, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#save 800

#plot A.iii (maxT centralness)
plot.gam(models$maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.15, 0.05), xlim = c(0, 800))

#add axes
axis(1, pos = -0.15, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#save 800



#plot B.i (minT relative polewardness)
plot.gam(models$minT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1)) 

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("SHAP value", side = 2, line = 6.5, cex = 2.5)

#save 800

#plot A.ii (meanT centralness)
plot.gam(models$meanT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1))

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#save 800

#plot A.iii (maxT centralness)
plot.gam(models$maxT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1))

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#save 800
