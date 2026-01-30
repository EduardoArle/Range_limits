#load libraries. 
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Significance_GAMs'

#load models
setwd(wd_models)
models <- lapply(list.files(pattern = 'absPol_GS'), readRDS)
names(models) <- list.files(pattern = 'absPol_GS')

#load significance tables
setwd(wd_sig)
sigs <- lapply(list.files(pattern = 'absPol_GS'), read.csv)
names(sigs) <- gsub('.csv', '', list.files(pattern = 'absPol_GS'))

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


#plot A.i (minT absolute polewardness)
plot_gam_shade_sig(models$minT_absPol_GS, select = 1, col = "#0000FF",
                   xlim = c(0, 90), ylim = c(-5, 2),
                   sig_x = sigs$minT_absPol_GS$absPolewardness,
                   sig = sigs$minT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -5,
                   xlab = "Absolute latitude",
                   ylab = "SHAP value")

#look at stats
summary(models$minT_absPol_GS)

#save 800

#plot A.ii (meanT absolute polewardness)x
plot_gam_shade_sig(models$meanT_absPol_GS, select = 1, col = "#800080",
                   xlim = c(0, 90), ylim = c(-5, 2),
                   sig_x = sigs$meanT_absPol_GS$absPolewardness,
                   sig = sigs$meanT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -5,
                   xlab = "Absolute latitude",
                   ylab = "")


#look at stats
summary(models$meanT_absPol_GS)

#save 800

#plot A.iii (maxT absolute polewardness)
plot_gam_shade_sig(models$maxT_absPol_GS, select = 1, col = "#FF0000",
                   xlim = c(0, 90), ylim = c(-5, 2),
                   sig_x = sigs$maxT_absPol_GS$absPolewardness,
                   sig = sigs$maxT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -5,
                   xlab = "Absolute latitude",
                   ylab = "")

#look at stats
summary(models$maxT_absPol_GS)

#save 800














################ PPT



#plot B.i (minPPT absolute polewardness)
plot_gam_shade_sig(models$minPPT_absPol_GS, select = 1, col = "#8c510a",
                   xlim = c(0, 90), ylim = c(-4, 2),
                   sig_x = sigs$minPPT_absPol_GS$absPolewardness,
                   sig = sigs$minPPT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -4,
                   xlab = "Absolute latitude",
                   ylab = "SHAP value")

#look at stats
summary(models$minPPT_absPol_GS)

#save 800

#plot B.ii (meanPPT absolute polewardness)
plot_gam_shade_sig(models$meanPPT_absPol_GS, select = 1, col = "#dfc27d",
                   xlim = c(0, 90), ylim = c(-4, 2),
                   sig_x = sigs$meanPPT_absPol_GS$absPolewardness,
                   sig = sigs$meanPPT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -4,
                   xlab = "Absolute latitude",
                   ylab = "")


#look at stats
summary(models$meanPPT_absPol_GS)

#save 800

#plot B.iii (maxPPT absolute polewardness)
plot_gam_shade_sig(models$maxPPT_absPol_GS, select = 1, col = "#01665e",
                   xlim = c(0, 90), ylim = c(-4, 2),
                   sig_x = sigs$maxPPT_absPol_GS$absPolewardness,
                   sig = sigs$maxPPT_absPol_GS$sig,
                   add_axes = TRUE, axis_pos = -4,
                   xlab = "Absolute latitude",
                   ylab = "")

#look at stats
summary(models$maxPPT_absPol_GS)

#save 800
