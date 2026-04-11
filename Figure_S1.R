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

#load model minPPT_relPol
setwd(wd_models)
model_minPPT_relPol <- readRDS("minPPT_relPol_GS" )

#load sig
setwd(wd_sig)
sig_minPPT_relPol <- read.csv("minPPT_relPol_GS.csv")

#plot A.i (minPPT relative polewardness)
plot_gam_shade_sig(model_minPPT_relPol, select = 1, col = "#8c510a",
                   xlim = c(0, 1), ylim = c(-0.2, 0.2),
                   sig_x = sig_minPPT_relPol$relPolewardness,
                   sig = sig_minPPT_relPol$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Relative polewardness",
                   ylab = "SHAP value")

#look at stats
summary(model_minPPT_relPol)

#save 800


#load model meanPPT_relPol
setwd(wd_models)
model_meanPPT_relPol <- readRDS("meanPPT_relPol_GS")

#load sig
setwd(wd_sig)
sig_meanPPT_relPol <- read.csv("meanPPT_relPol_GS.csv")

#plot A.ii (meanPPT relative polewardness)
plot_gam_shade_sig(model_meanPPT_relPol, select = 1, col = "#dfc27d",
                   xlim = c(0, 1), ylim = c(-0.2, 0.2),
                   sig_x = sig_meanPPT_relPol$relPolewardness,
                   sig = sig_meanPPT_relPol$sig,
                   add_axes = TRUE, axis_pos = -0.2,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(model_meanPPT_relPol)

#save 800


#load model maxPPT_relPol
setwd(wd_models)
model_maxPPT_relPol <- readRDS("maxPPT_relPol_GS")

#load sig
setwd(wd_sig)
sig_maxPPT_relPol <- read.csv("maxPPT_relPol_GS.csv")

#plot A.iii (maxT relative polewardness)
plot_gam_shade_sig(model_maxPPT_relPol, select = 1, col = "#01665e",
                   xlim = c(0, 1), ylim = c(-0.3, 0.3),
                   sig_x = sig_maxPPT_relPol$relPolewardness,
                   sig = sig_maxPPT_relPol$sig,
                   add_axes = TRUE, axis_pos = -0.3,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(model_maxPPT_relPol)

#save 800


#load model distEdge 250 (only sps with up to 20% range touching water)
setwd(wd_models)
model_minPPT_distEdge <- readRDS("minPPT_distEdge_250_20_GS")

#load sig
setwd(wd_sig)
sig_minPPT_distEdge <- read.csv("minPPT_distEdge_250_20_GS.csv")

#plot B.i (minT centralness 250km limit)
plot_gam_shade_sig(model_minPPT_distEdge, select = 1, col = "#8c510a",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sig_minPPT_distEdge$distEdge,
                   sig = sig_minPPT_distEdge$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "| SHAP value |")

#look at stats
summary(model_minPPT_distEdge)

#save 800


#load model distEdge 250 (only sps with up to 20% range touching water)
setwd(wd_models)
model_meanPPT_distEdge <- readRDS("meanPPT_distEdge_250_20_GS")

#load sig
setwd(wd_sig)
sig_meanPPT_distEdge <- read.csv("meanPPT_distEdge_250_20_GS.csv")

#plot B.ii (meanT centralness 250km limit)
plot_gam_shade_sig(model_meanPPT_distEdge, select = 1, col = "#dfc27d",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sig_meanPPT_distEdge$distEdge,
                   sig = sig_meanPPT_distEdge$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(model_meanPPT_distEdge)

#save 800


#load model distEdge 250 (only sps with up to 20% range touching water)
setwd(wd_models)
model_maxPPT_distEdge <- readRDS("maxPPT_distEdge_250_20_GS")

#load sig
setwd(wd_sig)
sig_maxPPT_distEdge <- read.csv("maxPPT_distEdge_250_20_GS.csv")

#plot B.iii (maxT centralness 250km limit)
plot_gam_shade_sig(model_maxPPT_distEdge, select = 1, col = "#01665e",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sig_maxPPT_distEdge$distEdge,
                   sig = sig_maxPPT_distEdge$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(model_maxPPT_distEdge)

#save 800


### extra column with mean T

#load model meanT_relPol
setwd(wd_models)
model_meanT_relPol <- readRDS("meanT_relPol_GS")

#load sig
setwd(wd_sig)
sig_meanT_relPol <- read.csv("meanT_relPol_GS.csv")

#plot A.iv (meanPPT relative polewardness)
plot_gam_shade_sig(model_meanT_relPol, select = 1, col = "#800080",
                   xlim = c(0, 1), ylim = c(-0.3, 0.3),
                   sig_x = sig_meanT_relPol$relPolewardness,
                   sig = sig_meanT_relPol$sig,
                   add_axes = TRUE, axis_pos = -0.3,
                   xlab = "Relative polewardness",
                   ylab = "")

#look at stats
summary(model_meanT_relPol)

#save 800


#load model distEdge 250 (only sps with up to 20% range touching water)
setwd(wd_models)
model_meanT_distEdge <- readRDS("meanT_distEdge_250_20_GS")

#load sig
setwd(wd_sig)
sig_meanT_distEdge <- read.csv("meanT_distEdge_250_20_GS.csv")

#plot B.iv (meanTPP centralness 250km limit)
plot_gam_shade_sig(model_meanT_distEdge, select = 1, col = "#800080",
                   xlim = c(0, 250), ylim = c(-0.4, 0.2),
                   sig_x = sig_meanT_distEdge$distEdge,
                   sig = sig_meanT_distEdge$sig,
                   add_axes = TRUE, axis_pos = -0.4,
                   xlab = "Distance to range edge",
                   ylab = "")

#look at stats
summary(model_meanT_distEdge)

#save 800
