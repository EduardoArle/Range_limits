#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250504_GAMs/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Significance_GAMs'

#load models
setwd(wd_models)
models <- lapply(list.files(), readRDS)
names(models) <- list.files()

#load significance tables
setwd(wd_sig)
sigs <- lapply(list.files(), read.csv)
names(sigs) <- gsub('.csv', '', list.files())

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#plot A.i (minT relative polewardness)
plot.gam(models$minT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1),
         col = NA, lwd = 0) 

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("SHAP value", side = 2, line = 6.5, cex = 2.5)

#add significant portions

#extract one valid level
valid_species <- levels(models$minT_relPol_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  relPolewardness = sigs$minT_relPol_GS$relPolewardness,
  species = factor(valid_species,
                   levels = levels(models$minT_relPol_GS$model$species))) 


pred <- predict(models$minT_relPol_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(relPolewardness)"]

#overlay the line with color-coded significance
sig_vec <- sigs$minT_relPol_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$relPolewardness[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}

#save 800


#plot A.ii (meanT relative polewardness)
plot.gam(models$meanT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1),
         col = NA, lwd = 0)

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#add significant portions

#extract one valid level
valid_species <- levels(models$minT_relPol_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  relPolewardness = sigs$meanT_relPol_GS$relPolewardness,
  species = factor(valid_species,
                   levels = levels(models$meanT_relPol_GS$model$species))) 


pred <- predict(models$meanT_relPol_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(relPolewardness)"]

#overlay the line with color-coded significance
sig_vec <- sigs$meanT_relPol_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$relPolewardness[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}



#save 800

#plot A.iii (maxT relative polewardness)
plot.gam(models$maxT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.06, 0.06), xlim = c(0, 1),
         col = NA, lwd = 0)

#add axes
axis(1, pos = -0.06, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#add significant portions

#extract one valid level
valid_species <- levels(models$maxT_relPol_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  relPolewardness = sigs$maxT_relPol_GS$relPolewardness,
  species = factor(valid_species,
                   levels = levels(models$maxT_relPol_GS$model$species))) 


pred <- predict(models$maxT_relPol_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(relPolewardness)"]

#overlay the line with color-coded significance
sig_vec <- sigs$maxT_relPol_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$relPolewardness[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}



#save 800

#plot B.i (minT centralness 250km limit)
plot.gam(models$minT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', 
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250),
         col = NA, lwd = 0) 

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  
mtext("| SHAP value |", side = 2, line = 6.5, cex = 2.5)

#add significant portions

#extract one valid level
valid_species <- levels(models$minT_distEdge_250_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  distEdge = sigs$minT_distEdge_250_GS$distEdge,
  species = factor(valid_species,
                   levels = levels(models$minT_distEdge_250_GS$model$species))) 

pred <- predict(models$minT_distEdge_250_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(distEdge)"]

#overlay the line with color-coded significance
sig_vec <- sigs$minT_distEdge_250_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$distEdge[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}

#save 800


#plot B.ii (meanT centralness 250km limit). Either this or full dist model in final version
plot.gam(models$meanT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250),
         col = NA, lwd = 0)

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#add significant portions

#extract one valid level
valid_species <- levels(models$meanT_distEdge_250_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  distEdge = sigs$meanT_distEdge_250_GS$distEdge,
  species = factor(valid_species,
                   levels = levels(models$meanT_distEdge_250_GS$model$species))) 

pred <- predict(models$meanT_distEdge_250_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(distEdge)"]

#overlay the line with color-coded significance
sig_vec <- sigs$meanT_distEdge_250_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$distEdge[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}

#save 800


#plot B.iii (maxT centralness 250km limit). Either this or full dist model in final version
plot.gam(models$maxT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030',
         xlab = "", ylab = "",
         axes = F, xaxs = "i", yaxs = "i",
         ylim = c(-0.04, 0.02), xlim = c(0, 250),
         col = NA, lwd = 0)

#add axes
axis(1, pos = -0.04, cex.axis = 2)
axis(2, las=2, cex.axis = 2)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#add significant portions

#extract one valid level
valid_species <- levels(models$maxT_distEdge_250_GS$model$species)[1]

#use that level in your new data
newdat <- data.frame(
  distEdge = sigs$maxT_distEdge_250_GS$distEdge,
  species = factor(valid_species,
                   levels = levels(models$maxT_distEdge_250_GS$model$species))) 

pred <- predict(models$maxT_distEdge_250_GS, newdata = newdat,
                type = "terms", se.fit = F)

fit_vals <- pred[, "s(distEdge)"]

#overlay the line with color-coded significance
sig_vec <- sigs$maxT_distEdge_250_GS$sig
rle_sig <- rle(sig_vec)
idx <- cumsum(rle_sig$lengths)

start <- 1
for (i in seq_along(rle_sig$lengths)) {
  end <- idx[i]
  segment_x <- newdat$distEdge[start:end]
  segment_y <- fit_vals[start:end]
  
  col <- if (rle_sig$values[i]) "black" else "darkgray"  # or any other styling
  lty <- if (rle_sig$values[i]) 1 else 2
  lines(segment_x, segment_y, col = col, lty = lty, lwd = 2)
  
  start <- end + 1
}

#save 800



