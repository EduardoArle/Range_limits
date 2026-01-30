#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Significance_GAMs'

#install plot_gam_shade_sig
plot_gam_shade_sig <- function(model, select = 1, xlim = NULL, ylim = NULL,
                               col = NULL, alpha = 0.2, n = 200,
                               sig_x = NULL, sig = NULL,
                               sig_col = "black", nsig_col = "darkgray",
                               sig_lty = 1, nsig_lty = 2, lwd = 2,
                               add_axes = FALSE,
                               axis_pos = NULL, cex_axis = 2, las = 2,
                               xlab = NULL, ylab = NULL,
                               xlab_line = 3.8, ylab_line = 6.5,
                               cex_lab = 2.5){
  
  #extract smooth term name
  xname <- model$smooth[[select]]$term[1]
  
  #build reference newdata
  nd0 <- as.data.frame(lapply(model$model, function(z){
    if(is.numeric(z)) mean(z, na.rm = TRUE) else levels(z)[1]
  }))
  
  #x range
  if(is.null(xlim)){
    xrng <- range(model$model[[xname]], na.rm = TRUE)
  } else {
    xrng <- xlim
  }
  
  #x grid
  nd <- nd0[rep(1, n), , drop = FALSE]
  nd[[xname]] <- seq(xrng[1], xrng[2], length.out = n)
  
  #predict smooth
  pr <- predict(model, newdata = nd, type = "terms", se.fit = TRUE)
  fit <- pr$fit[, select]
  se  <- pr$se.fit[, select]
  x   <- nd[[xname]]
  
  #y range
  if(is.null(ylim)){
    yrng <- range(c(fit - 2*se, fit + 2*se), na.rm = TRUE)
  } else {
    yrng <- ylim
  }
  
  #empty plot
  plot(x, fit, type = "n",
       xlab = "", ylab = "",
       axes = FALSE, xaxs = "i", yaxs = "i",
       xlim = xrng, ylim = yrng)
  
  #shade only
  polygon(c(x, rev(x)),
          c(fit - 2*se, rev(fit + 2*se)),
          col = grDevices::adjustcolor(col, alpha.f = alpha),
          border = NA)
  
  #add significant portions (optional)
  if(!is.null(sig_x) && !is.null(sig)){
    newdat <- nd0[rep(1, length(sig_x)), , drop = FALSE]
    newdat[[xname]] <- sig_x
    
    if("species" %in% names(model$model)){
      valid_species <- levels(model$model$species)[1]
      newdat$species <- factor(valid_species,
                               levels = levels(model$model$species))
    }
    
    pred2 <- predict(model, newdata = newdat, type = "terms", se.fit = FALSE)
    fit_vals <- pred2[, select]
    
    rle_sig <- rle(sig)
    idx <- cumsum(rle_sig$lengths)
    
    start <- 1
    for(i in seq_along(rle_sig$lengths)){
      end <- idx[i]
      segment_x <- sig_x[start:end]
      segment_y <- fit_vals[start:end]
      cc <- if(rle_sig$values[i]) sig_col else nsig_col
      ll <- if(rle_sig$values[i]) sig_lty else nsig_lty
      lines(segment_x, segment_y, col = cc, lty = ll, lwd = lwd)
      start <- end + 1
    }
  }
  
  #axes
  if(add_axes){
    if(is.null(axis_pos)) axis_pos <- yrng[1]
    axis(1, pos = axis_pos, cex.axis = cex_axis)
    axis(2, las = las, cex.axis = cex_axis)
  }
  
  #labels
  if(!is.null(xlab)){
    mtext(xlab, side = 1, line = xlab_line, cex = cex_lab)
  }
  if(!is.null(ylab)){
    mtext(ylab, side = 2, line = ylab_line, cex = cex_lab)
  }
}

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#load models relPol
setwd(wd_models)
models <- lapply(list.files(pattern = 'relPol_GS'), readRDS)
names(models) <- list.files(pattern = 'relPol_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'relPol_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'relPol_GS'))


#plot A.i (minT relative polewardness)
plot_gam_shade_sig(models$minT_relPol_GS, select = 1, col = "#0000FF",
                   xlim = c(0, 1), ylim = c(-0.3, 0.3),
                   sig_x = sigs$minT_relPol_GS$relPolewardness,
                   sig = sigs$minT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.3,
                   xlab = "Relative polewardness",
                   ylab = "SHAP value")

#look at stats
summary(models$minT_relPol_GS)

#save 800


#plot A.ii (meanT relative polewardness)
plot_gam_shade_sig(models$meanT_relPol_GS, select = 1, col = "#800080",
                   xlim = c(0, 1), ylim = c(-0.3, 0.3),
                   sig_x = sigs$meanT_relPol_GS$relPolewardness,
                   sig = sigs$meanT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.3,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(models$meanT_relPol_GS)

#save 800


#plot A.iii (maxT relative polewardness)
plot_gam_shade_sig(models$maxT_relPol_GS, select = 1, col = "#FF0000",
                   xlim = c(0, 1), ylim = c(-0.3, 0.3),
                   sig_x = sigs$maxT_relPol_GS$relPolewardness,
                   sig = sigs$maxT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.3,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(models$maxT_relPol_GS)

#save 800


#load models distEdge 250
setwd(wd_models)
models <- lapply(list.files(pattern = 'distEdge_250_GS'), readRDS)
names(models) <- list.files(pattern = 'distEdge_250_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'distEdge_250_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'distEdge_250_GS'))


#plot B.i (minT centralness 250km limit)
plot_gam_shade_sig(models$minT_distEdge_250_GS, select = 1, col = "#0000FF",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$minT_distEdge_250_GS$distEdge,
                   sig = sigs$minT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "| SHAP value |")

#look at stats
summary(models$minT_distEdge_250_GS)


#save 800


#plot B.ii (meanT centralness 250km limit)
plot_gam_shade_sig(models$meanT_distEdge_250_GS, select = 1, col = "#800080",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$meanT_distEdge_250_GS$distEdge,
                   sig = sigs$meanT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$meanT_distEdge_250_GS)


#save 800

#plot B.iii (maxT centralness 250km limit)
plot_gam_shade_sig(models$maxT_distEdge_250_GS, select = 1, col = "#FF0000",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$maxT_distEdge_250_GS$distEdge,
                   sig = sigs$maxT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$maxT_distEdge_250_GS)

#save 800




##### alternative distance edge (only sps with up to 20% range touching water)

#load models distEdge 250
setwd(wd_models)
models <- lapply(list.files(pattern = 'distEdge_250_20_GS'), readRDS)
names(models) <- list.files(pattern = 'distEdge_250_20_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'distEdge_250_20_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'distEdge_250_20_GS'))


#plot B.i (minT centralness 250km limit)
plot_gam_shade_sig(models$minT_distEdge_250_20_GS, select = 1, col = "#0000FF",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$minT_distEdge_250_20_GS$distEdge,
                   sig = sigs$minT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "| SHAP value |")

#look at stats
summary(models$minT_distEdge_250_20_GS)


#save 800


#plot B.ii (meanT centralness 250km limit)
plot_gam_shade_sig(models$meanT_distEdge_250_20_GS, select = 1, col = "#800080",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$meanT_distEdge_250_20_GS$distEdge,
                   sig = sigs$meanT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$meanT_distEdge_250_20_GS)


#save 800

#plot B.iii (maxT centralness 250km limit)
plot_gam_shade_sig(models$maxT_distEdge_250_20_GS, select = 1, col = "#FF0000",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$maxT_distEdge_250_20_GS$distEdge,
                   sig = sigs$maxT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$maxT_distEdge_250_20_GS)

#save 800
















################ PPT


#load models relPol
setwd(wd_models)
models <- lapply(list.files(pattern = 'relPol_GS'), readRDS)
names(models) <- list.files(pattern = 'relPol_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'relPol_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'relPol_GS'))


#plot A.i (minPPT relative polewardness)
plot_gam_shade_sig(models$minPPT_relPol_GS, select = 1, col = "#8c510a",
                   xlim = c(0, 1), ylim = c(-0.2, 0.2),
                   sig_x = sigs$minPPT_relPol_GS$relPolewardness,
                   sig = sigs$minPPT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Relative polewardness",
                   ylab = "SHAP value")

#look at stats
summary(models$minPPT_relPol_GS)

#save 800


#plot A.ii (meanPPT relative polewardness)
plot_gam_shade_sig(models$meanPPT_relPol_GS, select = 1, col = "#dfc27d",
                   xlim = c(0, 1), ylim = c(-0.2, 0.2),
                   sig_x = sigs$meanPPT_relPol_GS$relPolewardness,
                   sig = sigs$meanPPT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(models$meanPPT_relPol_GS)

#save 800


#plot A.iii (maxPPT relative polewardness)
plot_gam_shade_sig(models$maxPPT_relPol_GS, select = 1, col = "#01665e",
                   xlim = c(0, 1), ylim = c(-0.2, 0.2),
                   sig_x = sigs$maxPPT_relPol_GS$relPolewardness,
                   sig = sigs$maxPPT_relPol_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(models$maxPPT_relPol_GS)

#save 800


#load models distEdge 250
setwd(wd_models)
models <- lapply(list.files(pattern = 'distEdge_250_GS'), readRDS)
names(models) <- list.files(pattern = 'distEdge_250_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'distEdge_250_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'distEdge_250_GS'))


#plot B.i (minPPT centralness 250km limit)
plot_gam_shade_sig(models$minPPT_distEdge_250_GS, select = 1, col = "#8c510a",
                   xlim = c(0, 250), ylim = c(-0.2, 0.2),
                   sig_x = sigs$minPPT_distEdge_250_GS$distEdge,
                   sig = sigs$minPPT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Distance to range edge",
                   ylab = "| SHAP value |")

#look at stats
summary(models$minPPT_distEdge_250_GS)


#save 800


#plot B.ii (meanPPT centralness 250km limit)
plot_gam_shade_sig(models$meanPPT_distEdge_250_GS, select = 1, col = "#dfc27d",
                   xlim = c(0, 250), ylim = c(-0.2, 0.2),
                   sig_x = sigs$meanPPT_distEdge_250_GS$distEdge,
                   sig = sigs$meanPPT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$meanPPT_distEdge_250_GS)


#save 800

#plot B.iii (maxPPT centralness 250km limit)
plot_gam_shade_sig(models$maxPPT_distEdge_250_GS, select = 1, col = "#01665e",
                   xlim = c(0, 250), ylim = c(-0.2, 0.2),
                   sig_x = sigs$maxPPT_distEdge_250_GS$distEdge,
                   sig = sigs$maxPPT_distEdge_250_GS$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$maxPPT_distEdge_250_GS)

#save 800




##### alternative distance edge (only sps with up to 20% range touching water)

#load models distEdge 250
setwd(wd_models)
models <- lapply(list.files(pattern = 'distEdge_250_20_GS'), readRDS)
names(models) <- list.files(pattern = 'distEdge_250_20_GS')

#load sigs
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'distEdge_250_20_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'distEdge_250_20_GS'))


#plot B.i (minPPT centralness 250km limit)
plot_gam_shade_sig(models$minPPT_distEdge_250_20_GS, select = 1, col = "#8c510a",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$minPPT_distEdge_250_20_GS$distEdge,
                   sig = sigs$minPPT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "| SHAP value |")

#look at stats
summary(models$minPPT_distEdge_250_20_GS)


#save 800


#plot B.ii (meanPPT centralness 250km limit)
plot_gam_shade_sig(models$meanPPT_distEdge_250_20_GS, select = 1, col = "#dfc27d",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$meanPPT_distEdge_250_20_GS$distEdge,
                   sig = sigs$meanPPT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$meanPPT_distEdge_250_20_GS)


#save 800

#plot B.iii (maxPPT centralness 250km limit)
plot_gam_shade_sig(models$maxPPT_distEdge_250_20_GS, select = 1, col = "#01665e",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sigs$maxPPT_distEdge_250_20_GS$distEdge,
                   sig = sigs$maxPPT_distEdge_250_20_GS$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(models$maxPPT_distEdge_250_20_GS)

#save 800

