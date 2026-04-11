#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('20260221_Slopes_split.csv')

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#set colours
col_eq <- "#d73027"
col_pol <- "#4575b4"


###### MinT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_minT_relPol <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$Cor_vars_minT),]
s_minT_relPol <- s_minT_relPol[
  complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ]

#create new columns with log values (boruta does not accept it...)
s_minT_relPol$rangeSize_log10 <- log(s_minT_relPol$rangeSize, 10)
s_minT_relPol$nOcc_log <- log(s_minT_relPol$nOcc)
s_minT_relPol$nOcc_EQ_log <- log(s_minT_relPol$nOcc_EQ)
s_minT_relPol$nOcc_POL_log <- log(s_minT_relPol$nOcc_Pol)



#plot meaningful variables against slope



################
## Range size ##
################



#subset EQ and POL
s_minT_relPol_EQ <- s_minT_relPol[complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ]
s_minT_relPol_POL <- s_minT_relPol[complete.cases(s_minT_relPol$slope_POL_minT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minT_relPol_EQ$rangeSize_log10, s_minT_relPol_POL$rangeSize_log10))
y_lim <- c(-15, 5)

#plot graph with poits separated by EQ-centre and POL-centre slopes
plot(s_minT_relPol_EQ$rangeSize_log10, s_minT_relPol_EQ$slope_EQ_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minT_relPol_POL$rangeSize_log10,
       s_minT_relPol_POL$slope_POL_minT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_minT_relPol$rangeSize_log10) -
  min(s_minT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_minT_relPol$rangeSize_log10),
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_minT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)


#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minT_relPol_EQ$rangeSize_log10
x_vals_POL <- s_minT_relPol_POL$rangeSize_log10
x_range <- range(c(x_vals_EQ, x_vals_POL))

#fit linear models EQ and POL
lin_mod_minT_EQ <- lm(s_minT_relPol_EQ$slope_EQ_minT_relPol ~ x_vals_EQ,
                      weights = s_minT_relPol_EQ$nOcc_log)
lin_mod_minT_POL <- lm(s_minT_relPol_POL$slope_POL_minT_relPol ~ x_vals_POL,
                       weights = s_minT_relPol_POL$nOcc_log)

length(s_minT_relPol_EQ$slope_EQ_minT_relPol)
length(s_minT_relPol_POL$slope_POL_minT_relPol)
summary(lin_mod_minT_EQ)
summary(lin_mod_minT_POL)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred_EQ <- predict(lin_mod_minT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")
y_pred_POL <- predict(lin_mod_minT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines from the y-axis onwards
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence intervala
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################



#subset EQ only
s_minT_relPol_EQ <- s_minT_relPol[complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol_EQ$latAmplitude) #xlim
y_lim <- c(-1, 3) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_minT_relPol_EQ$latAmplitude, s_minT_relPol_EQ$slope_EQ_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) #xaxis
axis(2, pos = x_lim[1], las=2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minT_relPol_EQ$latAmplitude #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_minT_EQ <- lm(s_minT_relPol_EQ$slope_EQ_minT_relPol ~ x_vals_EQ,
                      weights = s_minT_relPol_EQ$nOcc_log) #lmEQ

length(s_minT_relPol_EQ$slope_EQ_minT_relPol) #nEQ
summary(lin_mod_minT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_minT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800



############################
### Latitudinal position ###
############################



#subset POL only
s_minT_relPol_POL <- s_minT_relPol[complete.cases(s_minT_relPol$slope_POL_minT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol_POL$rangeLoc) #xlim
x_lim[1] <- 0 #x0
y_lim <- c(-5, 5) #ylim

#plot graph (POL only)
plot(s_minT_relPol_POL$rangeLoc, s_minT_relPol_POL$slope_POL_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_pol, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add lines to fix plot ugliness
lines(c(x_lim[1], x_lim[1]), c(y_lim[1], y_lim[1] + diff(y_lim)*0.02)) #cornerY
lines(c(x_lim[1], x_lim[1] + diff(x_lim)*0.02), c(y_lim[1], y_lim[1])) #cornerX

#add axes lables 
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_POL <- s_minT_relPol_POL$rangeLoc #xvalsPOL
x_range <- range(x_vals_POL) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (POL)
lin_mod_minT_POL <- lm(s_minT_relPol_POL$slope_POL_minT_relPol ~ x_vals_POL,
                       weights = s_minT_relPol_POL$nOcc_log) #lmPOL

length(s_minT_relPol_POL$slope_POL_minT_relPol) #nPOL
summary(lin_mod_minT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_vals_POL) + range_val/100, max(x_vals_POL) - range_val/100,
             length.out = 100) #xseq
y_pred_POL <- predict(lin_mod_minT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression line
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



########################
### Median elevation ###
########################



#subset EQ only
s_minT_relPol_EQ <- s_minT_relPol[complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol_EQ$elevMedian) #xlim
y_lim <- c(-2, 2) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_minT_relPol_EQ$elevMedian, s_minT_relPol_EQ$slope_EQ_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = sort(unique(c(y_lim[1], pretty(y_lim)))),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minT_relPol_EQ$elevMedian #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_minT_EQ <- lm(s_minT_relPol_EQ$slope_EQ_minT_relPol ~ x_vals_EQ,
                      weights = s_minT_relPol_EQ$nOcc_log) #lmEQ

length(s_minT_relPol_EQ$slope_EQ_minT_relPol) #nEQ
summary(lin_mod_minT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_minT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800



#############################
### Elevational amplitude ###
#############################



#subset EQ only
s_minT_relPol_EQ <- s_minT_relPol[complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol_EQ$elevAmplitude) #xlim
y_lim <- c(-1, 2) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_minT_relPol_EQ$elevAmplitude, s_minT_relPol_EQ$slope_EQ_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minT_relPol_EQ$elevAmplitude #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_minT_EQ <- lm(s_minT_relPol_EQ$slope_EQ_minT_relPol ~ x_vals_EQ,
                      weights = s_minT_relPol_EQ$nOcc_log) #lmEQ

length(s_minT_relPol_EQ$slope_EQ_minT_relPol) #nEQ
summary(lin_mod_minT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_minT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800



#########################
### Number of records ###
#########################



#subset POL only
s_minT_relPol_POL <- s_minT_relPol[complete.cases(s_minT_relPol$slope_POL_minT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol_POL$nOcc_POL_log) #xlim
y_lim <- c(-8, 8) #ylim  (adjust if needed)

#plot graph (POL only)
plot(s_minT_relPol_POL$nOcc_POL_log,
     s_minT_relPol_POL$slope_POL_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_pol, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes labels 
#mtext('Number of records (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_POL <- s_minT_relPol_POL$nOcc_POL_log #xvalsPOL
x_range <- range(x_vals_POL) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (POL)
lin_mod_minT_POL <- lm(s_minT_relPol_POL$slope_POL_minT_relPol ~ x_vals_POL,
                       weights = s_minT_relPol_POL$nOcc_log) #lmPOL

length(s_minT_relPol_POL$slope_POL_minT_relPol) #nPOL
summary(lin_mod_minT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_vals_POL) + range_val/100,
             max(x_vals_POL) - range_val/100,
             length.out = 100) #xseq
y_pred_POL <- predict(lin_mod_minT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression line
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800





###### meanT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_meanT_relPol <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,] #corfilter
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$Cor_vars_meanT),] #corcomplete

#we do not filter by slope here (we subset EQ and POL later for plotting)

#create new columns with log values (boruta does not accept it...)
s_meanT_relPol$rangeSize_log10 <- log(s_meanT_relPol$rangeSize, 10) #lograngesize
s_meanT_relPol$nOcc_log <- log(s_meanT_relPol$nOcc) #lognocc
s_meanT_relPol$nOcc_EQ_log <- log(s_meanT_relPol$nOcc_EQ) #lognoccEQ
s_meanT_relPol$nOcc_POL_log <- log(s_meanT_relPol$nOcc_Pol) #lognoccPOL

#plot meaningful variables against slope



################
## Range size ##
################



#subset EQ and POL
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ] #subsetEQ
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanT_relPol_EQ$rangeSize_log10, s_meanT_relPol_POL$rangeSize_log10)) #xlim
y_lim <- c(-15, 5) #ylim

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanT_relPol_EQ$rangeSize_log10, s_meanT_relPol_EQ$slope_EQ_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanT_relPol_POL$rangeSize_log10,
       s_meanT_relPol_POL$slope_POL_meanT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_meanT_relPol$rangeSize_log10) -
  min(s_meanT_relPol$rangeSize_log10) #rangevals

plotting_positions <- c(min(s_meanT_relPol$rangeSize_log10),
                        min(s_meanT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_meanT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_meanT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_meanT_relPol$rangeSize_log10)) #xpos

plotting_values <- round(10 ^ plotting_positions / 1000) #xlabvals

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) #xaxis
axis(2, pos = x_lim[1], las=2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanT_relPol_EQ$rangeSize_log10 #xvalsEQ
x_vals_POL <- s_meanT_relPol_POL$rangeSize_log10 #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange

#fit linear models EQ and POL
lin_mod_meanT_EQ <- lm(s_meanT_relPol_EQ$slope_EQ_meanT_relPol ~ x_vals_EQ,
                       weights = s_meanT_relPol_EQ$nOcc_log) #lmEQ
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_EQ$slope_EQ_meanT_relPol) #nEQ
length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_EQ) #sumEQ
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression lines from the y-axis onwards
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence intervals
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



#############################
### Latitudinal amplitude ###
#############################



#subset EQ and POL
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ] #subsetEQ
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanT_relPol_EQ$latAmplitude, s_meanT_relPol_POL$latAmplitude)) #xlim
y_lim <- c(-1, 3) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanT_relPol_EQ$latAmplitude, s_meanT_relPol_EQ$slope_EQ_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanT_relPol_POL$latAmplitude,
       s_meanT_relPol_POL$slope_POL_meanT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanT_relPol_EQ$latAmplitude #xvalsEQ
x_vals_POL <- s_meanT_relPol_POL$latAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanT_EQ <- lm(s_meanT_relPol_EQ$slope_EQ_meanT_relPol ~ x_vals_EQ,
                       weights = s_meanT_relPol_EQ$nOcc_log) #lmEQ
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_EQ$slope_EQ_meanT_relPol) #nEQ
length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_EQ) #sumEQ
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence intervals
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



############################
### Latitudinal position ###
############################



#subset POL only
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol),] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol_POL$rangeLoc) #xlim
x_lim[1] <- 0 #x0
y_lim <- c(-5, 5) #ylim  (adjust if needed)

#plot graph (POL only)
plot(s_meanT_relPol_POL$rangeLoc, s_meanT_relPol_POL$slope_POL_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_pol, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_POL <- s_meanT_relPol_POL$rangeLoc #xvalsPOL
x_range <- range(x_vals_POL) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (POL)
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_vals_POL) + range_val/100, max(x_vals_POL) - range_val/100,
             length.out = 100) #xseq
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression line
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



########################
### Median elevation ###
########################



#subset EQ and POL
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ] #subsetEQ
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanT_relPol_EQ$elevMedian, s_meanT_relPol_POL$elevMedian)) #xlim
y_lim <- c(-2, 2) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanT_relPol_EQ$elevMedian, s_meanT_relPol_EQ$slope_EQ_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanT_relPol_POL$elevMedian,
       s_meanT_relPol_POL$slope_POL_meanT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanT_relPol_EQ$elevMedian #xvalsEQ
x_vals_POL <- s_meanT_relPol_POL$elevMedian #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanT_EQ <- lm(s_meanT_relPol_EQ$slope_EQ_meanT_relPol ~ x_vals_EQ,
                       weights = s_meanT_relPol_EQ$nOcc_log) #lmEQ
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_EQ$slope_EQ_meanT_relPol) #nEQ
length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_EQ) #sumEQ
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence intervals
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



#############################
### Elevational amplitude ###
#############################



#subset EQ and POL
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ] #subsetEQ
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanT_relPol_EQ$elevAmplitude,
                 s_meanT_relPol_POL$elevAmplitude)) #xlim
y_lim <- c(-1, 2) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanT_relPol_EQ$elevAmplitude, s_meanT_relPol_EQ$slope_EQ_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanT_relPol_POL$elevAmplitude,
       s_meanT_relPol_POL$slope_POL_meanT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanT_relPol_EQ$elevAmplitude #xvalsEQ
x_vals_POL <- s_meanT_relPol_POL$elevAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanT_EQ <- lm(s_meanT_relPol_EQ$slope_EQ_meanT_relPol ~ x_vals_EQ,
                       weights = s_meanT_relPol_EQ$nOcc_log) #lmEQ
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_EQ$slope_EQ_meanT_relPol) #nEQ
length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_EQ) #sumEQ
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence intervals
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



#################
### Roundness ###
#################



#subset POL only
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol_POL$roundness) #xlim
y_lim <- c(-2, 2) #ylim

#plot graph (POL only)
plot(s_meanT_relPol_POL$roundness, s_meanT_relPol_POL$slope_POL_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_pol, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_POL <- s_meanT_relPol_POL$roundness #xvalsPOL
x_range <- range(x_vals_POL) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (POL)
lin_mod_meanT_POL <- lm(s_meanT_relPol_POL$slope_POL_meanT_relPol ~ x_vals_POL,
                        weights = s_meanT_relPol_POL$nOcc_log) #lmPOL

length(s_meanT_relPol_POL$slope_POL_meanT_relPol) #nPOL
summary(lin_mod_meanT_POL) #sumPOL

#predict y-values
x_seq <- seq(min(x_vals_POL) + range_val/100, max(x_vals_POL) - range_val/100,
             length.out = 100) #xseq
y_pred_POL <- predict(lin_mod_meanT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression line
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL

#save 800



#########################
### Number of records ###
#########################



#subset EQ only
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol_EQ$nOcc_EQ_log) #xlim
y_lim <- c(-8, 8) #ylim

#plot graph (EQ only)
plot(s_meanT_relPol_EQ$nOcc_EQ_log,
     s_meanT_relPol_EQ$slope_EQ_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes labels 
#mtext('Number of records (log)', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanT_relPol_EQ$nOcc_EQ_log #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_meanT_EQ <- lm(s_meanT_relPol_EQ$slope_EQ_meanT_relPol ~ x_vals_EQ,
                       weights = s_meanT_relPol_EQ$nOcc_EQ_log) #lmEQ

length(s_meanT_relPol_EQ$slope_EQ_meanT_relPol) #nEQ
summary(lin_mod_meanT_EQ) #sumEQ

#predict y-values
x_seq <- seq(min(x_vals_EQ) + range_val/100,
             max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800





###### maxT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_maxT_relPol <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,] #corfilter
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$Cor_vars_maxT),] #corcomplete

#we do not filter by slope here (we subset EQ and POL later for plotting)

#create new columns with log values (boruta does not accept it...)
s_maxT_relPol$rangeSize_log10 <- log(s_maxT_relPol$rangeSize, 10) #lograngesize
s_maxT_relPol$nOcc_log <- log(s_maxT_relPol$nOcc) #lognocc
s_maxT_relPol$nOcc_EQ_log <- log(s_maxT_relPol$nOcc_EQ) #lognoccEQ
s_maxT_relPol$nOcc_POL_log <- log(s_maxT_relPol$nOcc_Pol) #lognoccPOL

#plot meaningful variables against slope



################
## Range size ##
################



#subset EQ and POL
s_maxT_relPol_EQ <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ] #subsetEQ
s_maxT_relPol_POL <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_POL_maxT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxT_relPol_EQ$rangeSize_log10, s_maxT_relPol_POL$rangeSize_log10)) #xlim
y_lim <- c(-15, 5) #ylim (same as minT/meanT for consistency)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxT_relPol_EQ$rangeSize_log10, s_maxT_relPol_EQ$slope_EQ_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxT_relPol_POL$rangeSize_log10,
       s_maxT_relPol_POL$slope_POL_maxT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_maxT_relPol$rangeSize_log10) -
  min(s_maxT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_maxT_relPol$rangeSize_log10),
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_maxT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) #xaxis
axis(2, pos = x_lim[1], las = 2, cex.axis = 2) #yaxis

#add axes labels
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxT_relPol_EQ$rangeSize_log10
x_vals_POL <- s_maxT_relPol_POL$rangeSize_log10
x_range <- range(c(x_vals_EQ, x_vals_POL))

#fit linear models EQ and POL
lin_mod_maxT_EQ <- lm(s_maxT_relPol_EQ$slope_EQ_maxT_relPol ~ x_vals_EQ,
                      weights = s_maxT_relPol_EQ$nOcc_log)

lin_mod_maxT_POL <- lm(s_maxT_relPol_POL$slope_POL_maxT_relPol ~ x_vals_POL,
                       weights = s_maxT_relPol_POL$nOcc_log)

length(s_maxT_relPol_EQ$slope_EQ_maxT_relPol)
length(s_maxT_relPol_POL$slope_POL_maxT_relPol)
summary(lin_mod_maxT_EQ)
summary(lin_mod_maxT_POL)

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_maxT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_maxT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence intervals
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



#############################
### Latitudinal amplitude ###
#############################



#subset EQ only
s_maxT_relPol_EQ <- s_maxT_relPol[complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol_EQ$latAmplitude) #xlim
y_lim <- c(-1, 3) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_maxT_relPol_EQ$latAmplitude, s_maxT_relPol_EQ$slope_EQ_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_maxT_relPol_EQ$latAmplitude #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_maxT_EQ <- lm(s_maxT_relPol_EQ$slope_EQ_maxT_relPol ~ x_vals_EQ,
                      weights = s_maxT_relPol_EQ$nOcc_log) #lmEQ

length(s_maxT_relPol_EQ$slope_EQ_maxT_relPol) #nEQ
summary(lin_mod_maxT_EQ) #sumEQ

#predict y-values
x_seq <- seq(min(x_vals_EQ) + range_val/100,
             max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800



########################
### Elevation median ###
########################



#subset EQ only
s_maxT_relPol_EQ <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol_EQ$elevMedian) #xlim
y_lim <- c(-2, 2) #ylim  # <-- keep or change manually for maxT

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_maxT_relPol_EQ$elevMedian, s_maxT_relPol_EQ$slope_EQ_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = sort(unique(c(y_lim[1], pretty(y_lim)))),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxT_relPol_EQ$elevMedian #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_maxT_EQ <- lm(s_maxT_relPol_EQ$slope_EQ_maxT_relPol ~ x_vals_EQ,
                      weights = s_maxT_relPol_EQ$nOcc_log) #lmEQ

length(s_maxT_relPol_EQ$slope_EQ_maxT_relPol) #nEQ
summary(lin_mod_maxT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ



############################
### Elevation amplitude  ###
############################



#subset EQ only
s_maxT_relPol_EQ <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol_EQ$elevAmplitude) #xlim
y_lim <- c(-1, 2) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_maxT_relPol_EQ$elevAmplitude, s_maxT_relPol_EQ$slope_EQ_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxT_relPol_EQ$elevAmplitude #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_maxT_EQ <- lm(s_maxT_relPol_EQ$slope_EQ_maxT_relPol ~ x_vals_EQ,
                      weights = s_maxT_relPol_EQ$nOcc_log) #lmEQ

length(s_maxT_relPol_EQ$slope_EQ_maxT_relPol) #nEQ
summary(lin_mod_maxT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800



#################
### Roundness ###
#################



#subset POL only
s_maxT_relPol_POL <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_POL_maxT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol_POL$roundness) #xlim
y_lim <- c(-2, 2) #ylim  # <-- set manually for this panel

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (POL only)
plot(s_maxT_relPol_POL$roundness, s_maxT_relPol_POL$slope_POL_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_pol, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_POL <- s_maxT_relPol_POL$roundness #xvalsPOL
x_range <- range(x_vals_POL) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (POL)
lin_mod_maxT_POL <- lm(s_maxT_relPol_POL$slope_POL_maxT_relPol ~ x_vals_POL,
                       weights = s_maxT_relPol_POL$nOcc_log) #lmPOL

length(s_maxT_relPol_POL$slope_POL_maxT_relPol) #nPOL
summary(lin_mod_maxT_POL) #sumPOL

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_POL) - range_val/100,
             length.out = 100) #xseq
y_pred_POL <- predict(lin_mod_maxT_POL, newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") #predPOL

#add regression line
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) #linePOL

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) #ciPOL


#save 800



###################
#### n Records ####
###################



#subset EQ only
s_maxT_relPol_EQ <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ] #subsetEQ

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol_EQ$nOcc_EQ_log) #xlim
y_lim <- c(-8, 8) #ylim  # <-- set manually for this panel

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph (EQ only)
plot(s_maxT_relPol_EQ$nOcc_EQ_log, s_maxT_relPol_EQ$slope_EQ_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2) #yaxis

#add axes lables 
#mtext('Number of records (EQ)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxT_relPol_EQ$nOcc_EQ_log #xvalsEQ
x_range <- range(x_vals_EQ) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear model (EQ)
lin_mod_maxT_EQ <- lm(s_maxT_relPol_EQ$slope_EQ_maxT_relPol ~ x_vals_EQ,
                      weights = s_maxT_relPol_EQ$nOcc_EQ_log) #lmEQ (weights = EQ nOcc)

length(s_maxT_relPol_EQ$slope_EQ_maxT_relPol) #nEQ
summary(lin_mod_maxT_EQ) #sumEQ

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ

#add regression line
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) #lineEQ

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) #ciEQ

#save 800




###### MinT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_minT_distEdge <- slopes[abs(slopes$Cor_vars_minT) <= 0.7, ]
s_minT_distEdge <- s_minT_distEdge[complete.cases(s_minT_distEdge$Cor_vars_minT), ]
s_minT_distEdge <- s_minT_distEdge[complete.cases(s_minT_distEdge$slope_minT_distEdge), ]

#create new columns with log values (boruta does not accept it...)
s_minT_distEdge$rangeSize_log10 <- log(s_minT_distEdge$rangeSize, 10)
s_minT_distEdge$nOcc_log <- log(s_minT_distEdge$nOcc)

#define colour
col_all <- "#303030"


#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$rangeSize_log10, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minT_distEdge$rangeSize_log10) -
  min(s_minT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_minT_distEdge$rangeSize_log10),
                        min(s_minT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_minT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_minT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_minT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)

axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_minT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$latAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$rangeLoc, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



########################
### Median elevation ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$elevMedian, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Median elevation', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Elevation amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$elevAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Elevation amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



###### MeanT vs DISTANCE TO RANGE EDGE  ######



#select only species that had lower correl between vars
s_meanT_distEdge <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_distEdge <- s_meanT_distEdge[complete.cases(s_meanT_distEdge$Cor_vars_meanT),]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$slope_meanT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_meanT_distEdge$rangeSize_log10 <- log(s_meanT_distEdge$rangeSize, 10)
s_meanT_distEdge$nOcc_log <- log(s_meanT_distEdge$nOcc)



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)  # <-- set manually for meanT if needed

#plot graph
plot(s_meanT_distEdge$rangeSize_log10, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_meanT_distEdge$rangeSize_log10) -
  min(s_meanT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_meanT_distEdge$rangeSize_log10),
                        min(s_meanT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_meanT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_meanT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_meanT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                    weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$latAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                    weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# # Define x and y limits
# x_lim <- range(s_meanT_distEdge$rangeLoc)
# y_lim <- c(-0.02, 0.02)
# 
# # Ensure x_lim starts at 0 for a clean intersection
# x_lim[1] <- 0 
# 
# #plot graph
# plot(s_meanT_distEdge$rangeLoc, s_meanT_distEdge$slope_meanT_distEdge,
#      axes = F, xaxs = "i", yaxs = "i",
#      xlab = "", ylab = "", cex = 1.5,
#      pch = 19, col = '#90909020',
#      ylim = y_lim, xlim = x_lim)
# 
# #add axes
# axis(1, pos = y_lim[1], cex.axis = 2)
# axis(2, pos = x_lim[1], las=2, cex.axis = 2)
# 
# 
# #add axes lables 
# #mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_meanT_distEdge$rangeLoc
# x_range <- range(x_vals)
# range_val <- x_range[2] - x_range[1]
# 
# #fit linear model
# lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
#                    weights = s_meanT_distEdge$nOcc_log)
# 
# length(s_meanT_distEdge$slope_meanT_distEdge)
# summary(lin_mod_meanT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
#              length.out = 100)  # Ensuring it starts at 0
# y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
#                   interval = "confidence")
# 
# #add regression line from the y-axis onwards
# lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)
# 
# #dd shaded confidence interval
# polygon(c(x_seq, rev(x_seq)),
#         c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
#         col = '#90909030',
#         border = NA)
# 
# 
# #save 800



########################
### Median elevation ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$elevMedian, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Median elevation', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                    weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
###  Elevation amplitude  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$elevAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Elevation amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                    weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#################
### Roundness ###
#################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$roundness, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                    weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



###################
#### n Records ####
###################


# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- range(s_meanT_distEdge$nOcc_log)
# y_lim <- c(-0.001, 0.001)
# 
# #plot graph
# plot(s_meanT_distEdge$nOcc_log, s_meanT_distEdge$slope_meanT_distEdge,
#      axes = F, xaxs = "i", yaxs = "i",
#      xlab = "", ylab = "", cex = 1.5,
#      pch = 19, col = '#90909020',
#      ylim = y_lim, xlim = x_lim)
# 
# 
# # Define original values and their log-transformed positions
# range_vals <- max(s_meanT_distEdge$nOcc_log) -
#   min(s_meanT_distEdge$nOcc_log)
# 
# 
# plotting_positions <- c(min(s_meanT_distEdge$nOcc_log),
#                         min(s_meanT_distEdge$nOcc_log) + range_vals/4,
#                         min(s_meanT_distEdge$nOcc_log) + range_vals/2,
#                         min(s_meanT_distEdge$nOcc_log) + range_vals/4*3,
#                         max(s_meanT_distEdge$nOcc_log))
# 
# plotting_values <- round(exp(plotting_positions))
# 
# #add axes
# axis(1, at = plotting_positions, labels = plotting_values,
#      pos = y_lim[1], cex.axis = 2)
# axis(2, pos = x_lim[1], las=2, cex.axis = 2)
# 
# 
# #add axes lables 
# #mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_meanT_distEdge$nOcc_log
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
#                    weights = s_meanT_distEdge$nOcc_log)
# 
# length(s_meanT_distEdge$slope_meanT_distEdge)
# summary(lin_mod_meanT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(min(plotting_positions) + range_vals/100,
#              max(plotting_positions) - range_vals/100,
#              length.out = 100)  # Ensuring it starts at 0
# y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
#                   interval = "confidence")
# 
# #add regression line from the y-axis onwards
# lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)
# 
# #add shaded confidence interval
# polygon(c(x_seq, rev(x_seq)),
#         c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
#         col = '#90909030',
#         border = NA)
# 
# 
# 
# #save 800




###### MaxT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_maxT_distEdge <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_distEdge <- s_maxT_distEdge[complete.cases(s_maxT_distEdge$Cor_vars_maxT),]
s_maxT_distEdge <- s_maxT_distEdge[complete.cases(s_maxT_distEdge$slope_maxT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_maxT_distEdge$rangeSize_log10 <- log(s_maxT_distEdge$rangeSize, 10)
s_maxT_distEdge$nOcc_log <- log(s_maxT_distEdge$nOcc)



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)  # <-- set manually for maxT if needed

#plot graph
plot(s_maxT_distEdge$rangeSize_log10, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_maxT_distEdge$rangeSize_log10) -
  min(s_maxT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_maxT_distEdge$rangeSize_log10),
                        min(s_maxT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_maxT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_maxT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_maxT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Range size', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$latAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$rangeLoc, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Latitudinal position', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



########################
### Median elevation ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$elevMedian, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Median elevation', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800



#################
### Roundness ###
#################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$roundness, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Range shape', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxT,
                  newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)

#save 800


#### empty plots to make labels


#########################
## Elevation amplitude ##
#########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$elevAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "00"),
     ylim = y_lim, xlim = x_lim)


#add axes lables 
mtext('Elevation amplitude', side = 1, line = 5.6, cex = 2.5)  

#save 800




#add axes lables 
mtext('Range shape', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#plot graph
plot(s_maxT_distEdge$elevAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Body mass', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#plot graph
plot(s_maxT_distEdge$elevAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Number of records', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)
