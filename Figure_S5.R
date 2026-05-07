##### important note: fix labels in x axis in nOcc and bodyMass matching 
##### what I did for rnge size

#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('20260428_Slopes_split.csv')

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

#set colours
col_eq <- "#d73027"
col_pol <- "#4575b4"

#table to store regression results
reg_results <- data.frame(response = character(),
                          gradient = character(),
                          portion = character(),
                          predictor = character(),
                          estimate = numeric(),
                          std_error = numeric(),
                          p_value = numeric(),
                          n = numeric(),
                          r_squared = numeric(),
                          adj_r_squared = numeric(),
                          weight_variable = character(),
                          stringsAsFactors = FALSE)

#function to extract model results
store_lm <- function(tab, model, response, gradient, portion,
                     predictor, n, weight_variable)
{
  sm <- summary(model)
  
  new_row <- data.frame(response = response,
                        gradient = gradient,
                        portion = portion,
                        predictor = predictor,
                        estimate = sm$coefficients[2, 1],
                        std_error = sm$coefficients[2, 2],
                        p_value = sm$coefficients[2, 4],
                        n = n,
                        r_squared = sm$r.squared,
                        adj_r_squared = sm$adj.r.squared,
                        weight_variable = weight_variable,
                        stringsAsFactors = FALSE)
  
  rbind(tab, new_row)
}


###### MinPPT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_minPPT_relPol <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_relPol <- s_minPPT_relPol[complete.cases(s_minPPT_relPol$Cor_vars_minPPT),]
s_minPPT_relPol <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]

#create new columns with log values (boruta does not accept it...)
s_minPPT_relPol$rangeSize_log10 <- log(s_minPPT_relPol$rangeSize, 10)
s_minPPT_relPol$nOcc_log <- log(s_minPPT_relPol$nOcc)
s_minPPT_relPol$nOcc_EQ_log <- log(s_minPPT_relPol$nOcc_EQ)
s_minPPT_relPol$nOcc_POL_log <- log(s_minPPT_relPol$nOcc_Pol)
s_minPPT_relPol$bodyMass_log <- log(s_minPPT_relPol$bodyMass)


#plot attributes against slope



################
## Range size ##
################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$rangeSize_log10,
                 s_minPPT_relPol_POL$rangeSize_log10))
y_lim <- c(-5, 5)

#plot graph with poits separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$rangeSize_log10, s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$rangeSize_log10,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_minPPT_relPol$rangeSize_log10) -
  min(s_minPPT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_minPPT_relPol$rangeSize_log10),
                        min(s_minPPT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_minPPT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_minPPT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_minPPT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$rangeSize_log10
x_vals_POL <- s_minPPT_relPol_POL$rangeSize_log10
x_range <- range(c(x_vals_EQ, x_vals_POL))

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeSize_log10',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeSize_log10',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
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



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ] 
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ] 

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$latAmplitude,
                 s_minPPT_relPol_POL$latAmplitude)) 
y_lim <- c(-3, 3) 

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$latAmplitude,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$latAmplitude,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1, pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) 
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$latAmplitude 
x_vals_POL <- s_minPPT_relPol_POL$latAmplitude 
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1] 

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log) 
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log) 

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'latAmplitude',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'latAmplitude',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100) 

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") 

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence") 

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6) 
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6) 

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA) 

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA) 


#save 800



############################
### Latitudinal position ###
############################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$rangeLoc,
                 s_minPPT_relPol_POL$rangeLoc))
x_lim[1] <- 0
y_lim <- c(-5, 5)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$rangeLoc,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$rangeLoc,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1, at = pretty(x_lim), pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add lines to fix plot ugliness
lines(c(x_lim[1], x_lim[1]),
      c(y_lim[1], y_lim[1] + diff(y_lim)*0.02))
lines(c(x_lim[1], x_lim[1] + diff(x_lim)*0.02),
      c(y_lim[1], y_lim[1]))

#add axes lables 
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$rangeLoc
x_vals_POL <- s_minPPT_relPol_POL$rangeLoc
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeLoc',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeLoc',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_vals_EQ) + range_val/100,
             max(x_vals_EQ) - range_val/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



########################
### Median elevation ###
########################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$elevMedian,
                 s_minPPT_relPol_POL$elevMedian))
y_lim <- c(-2, 2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$elevMedian,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$elevMedian,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2)
axis(2,
     at = sort(unique(c(y_lim[1], pretty(y_lim)))),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$elevMedian
x_vals_POL <- s_minPPT_relPol_POL$elevMedian
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevMedian',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevMedian',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
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
### Elevational amplitude ###
#############################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$elevAmplitude,
                 s_minPPT_relPol_POL$elevAmplitude))
y_lim <- c(-3, 3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$elevAmplitude,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$elevAmplitude,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2)
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$elevAmplitude
x_vals_POL <- s_minPPT_relPol_POL$elevAmplitude
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevAmplitude',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevAmplitude',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals_EQ) - range_val/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



#################
### Roundness ###
#################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$roundness,
                 s_minPPT_relPol_POL$roundness))
y_lim <- c(-2, 2)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$roundness,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$roundness,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2)
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$roundness
x_vals_POL <- s_minPPT_relPol_POL$roundness
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

#check sample sizes and inspect regression summaries (EQ and POL portions)
length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'roundness',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'roundness',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_vals_EQ) + range_val/100,
             max(x_vals_EQ) - range_val/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



#########################
### Number of records ###
#########################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                 s_minPPT_relPol_POL$nOcc_POL_log))
y_lim <- c(-3, 3)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$nOcc_EQ_log,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$nOcc_POL_log,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                    s_minPPT_relPol_POL$nOcc_POL_log)) -
  min(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
        s_minPPT_relPol_POL$nOcc_POL_log))

plotting_positions <- c(min(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                              s_minPPT_relPol_POL$nOcc_POL_log)),
                        min(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                              s_minPPT_relPol_POL$nOcc_POL_log)) + range_vals/4,
                        min(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                              s_minPPT_relPol_POL$nOcc_POL_log)) + range_vals/2,
                        min(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                              s_minPPT_relPol_POL$nOcc_POL_log)) + range_vals/4*3,
                        max(c(s_minPPT_relPol_EQ$nOcc_EQ_log,
                              s_minPPT_relPol_POL$nOcc_POL_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes labels 
#mtext('Number of records (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$nOcc_EQ_log
x_vals_POL <- s_minPPT_relPol_POL$nOcc_POL_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'nOcc_EQ_log',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'nOcc_POL_log',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)


#save 800



################
### Body mass ###
################



#subset EQ and POL
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol,
                 s_minPPT_relPol$bodyMass_log), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol,
                 s_minPPT_relPol$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_minPPT_relPol_EQ$bodyMass_log,
                 s_minPPT_relPol_POL$bodyMass_log))
y_lim <- c(-2, 2)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_minPPT_relPol_EQ$bodyMass_log,
     s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_eq, "20"),
     ylim = y_lim,
     xlim = x_lim)

#addPOLpts
points(s_minPPT_relPol_POL$bodyMass_log,
       s_minPPT_relPol_POL$slope_POL_minPPT_relPol,
       pch = 19,
       col = paste0(col_pol, "20"),
       cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_minPPT_relPol_EQ$bodyMass_log,
                    s_minPPT_relPol_POL$bodyMass_log)) -
  min(c(s_minPPT_relPol_EQ$bodyMass_log,
        s_minPPT_relPol_POL$bodyMass_log))

plotting_positions <- c(min(c(s_minPPT_relPol_EQ$bodyMass_log,
                              s_minPPT_relPol_POL$bodyMass_log)),
                        min(c(s_minPPT_relPol_EQ$bodyMass_log,
                              s_minPPT_relPol_POL$bodyMass_log)) + range_vals/4,
                        min(c(s_minPPT_relPol_EQ$bodyMass_log,
                              s_minPPT_relPol_POL$bodyMass_log)) + range_vals/2,
                        min(c(s_minPPT_relPol_EQ$bodyMass_log,
                              s_minPPT_relPol_POL$bodyMass_log)) + range_vals/4*3,
                        max(c(s_minPPT_relPol_EQ$bodyMass_log,
                              s_minPPT_relPol_POL$bodyMass_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)

axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Body mass (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_minPPT_relPol_EQ$bodyMass_log
x_vals_POL <- s_minPPT_relPol_POL$bodyMass_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_minPPT_EQ <- lm(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol ~ x_vals_EQ,
                        weights = s_minPPT_relPol_EQ$nOcc_EQ_log)

lin_mod_minPPT_POL <- lm(s_minPPT_relPol_POL$slope_POL_minPPT_relPol ~ x_vals_POL,
                         weights = s_minPPT_relPol_POL$nOcc_POL_log)

length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol)
length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol)
summary(lin_mod_minPPT_EQ)
summary(lin_mod_minPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT_EQ,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'bodyMass_log',
                        n = length(s_minPPT_relPol_EQ$slope_EQ_minPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT_POL,
                        response = 'minPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'bodyMass_log',
                        n = length(s_minPPT_relPol_POL$slope_POL_minPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_minPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_minPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



###### meanPPT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_meanPPT_relPol <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,] #corfilter
s_meanPPT_relPol <- s_meanPPT_relPol[complete.cases(s_meanPPT_relPol$Cor_vars_meanPPT),] #corcomplete

#we do not filter by slope here (we subset EQ and POL later for plotting)

#create new columns with log values (boruta does not accept it...)
s_meanPPT_relPol$rangeSize_log10 <- log(s_meanPPT_relPol$rangeSize, 10) #lograngesize
s_meanPPT_relPol$nOcc_log <- log(s_meanPPT_relPol$nOcc) #lognocc
s_meanPPT_relPol$nOcc_EQ_log <- log(s_meanPPT_relPol$nOcc_EQ) #lognoccEQ
s_meanPPT_relPol$nOcc_POL_log <- log(s_meanPPT_relPol$nOcc_Pol) #lognoccPOL
s_meanPPT_relPol$bodyMass_log <- log(s_meanPPT_relPol$bodyMass) #logbodymass

#plot variables against slope



################
## Range size ##
################



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$rangeSize_log10, s_meanPPT_relPol_POL$rangeSize_log10)) #xlim
y_lim <- c(-5, 5) #ylim

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$rangeSize_log10, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$rangeSize_log10,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_meanPPT_relPol$rangeSize_log10) -
  min(s_meanPPT_relPol$rangeSize_log10) #rangevals

plotting_positions <- c(min(s_meanPPT_relPol$rangeSize_log10),
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_meanPPT_relPol$rangeSize_log10)) #xpos

plotting_values <- round(10 ^ plotting_positions / 1000) #xlabvals

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$rangeSize_log10 #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$rangeSize_log10 #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeSize_log10',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeSize_log10',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$latAmplitude, s_meanPPT_relPol_POL$latAmplitude)) #xlim
y_lim <- c(-3, 3) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$latAmplitude, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$latAmplitude,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes (forced ticks so axes meet)
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$latAmplitude #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$latAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'latAmplitude',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'latAmplitude',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$rangeLoc, s_meanPPT_relPol_POL$rangeLoc)) #xlim
y_lim <- c(-5, 5) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$rangeLoc, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$rangeLoc,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
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
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$rangeLoc #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$rangeLoc #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeLoc',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeLoc',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



########################
### Median elevation ###
########################



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$elevMedian, s_meanPPT_relPol_POL$elevMedian)) #xlim
y_lim <- c(-2, 2) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$elevMedian, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$elevMedian,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
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
x_vals_EQ <- s_meanPPT_relPol_EQ$elevMedian #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$elevMedian #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevMedian',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevMedian',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$elevAmplitude,
                 s_meanPPT_relPol_POL$elevAmplitude)) #xlim
y_lim <- c(-3, 3) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$elevAmplitude, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$elevAmplitude,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$elevAmplitude #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$elevAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevAmplitude',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevAmplitude',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ] #subsetEQ
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$roundness, s_meanPPT_relPol_POL$roundness)) #xlim
y_lim <- c(-2, 2) #ylim

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$roundness, s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$roundness,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
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
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$roundness #xvalsEQ
x_vals_POL <- s_meanPPT_relPol_POL$roundness #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol) #nEQ
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol) #nPOL
summary(lin_mod_meanPPT_EQ) #sumEQ
summary(lin_mod_meanPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'roundness',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'roundness',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_meanPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_meanPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#########################
### Number of records ###
#########################



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ]
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                 s_meanPPT_relPol_POL$nOcc_POL_log))
y_lim <- c(-3, 3)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$nOcc_EQ_log,
     s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$nOcc_POL_log,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                    s_meanPPT_relPol_POL$nOcc_POL_log)) -
  min(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
        s_meanPPT_relPol_POL$nOcc_POL_log))

plotting_positions <- c(min(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                              s_meanPPT_relPol_POL$nOcc_POL_log)),
                        min(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                              s_meanPPT_relPol_POL$nOcc_POL_log)) + range_vals/4,
                        min(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                              s_meanPPT_relPol_POL$nOcc_POL_log)) + range_vals/2,
                        min(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                              s_meanPPT_relPol_POL$nOcc_POL_log)) + range_vals/4*3,
                        max(c(s_meanPPT_relPol_EQ$nOcc_EQ_log,
                              s_meanPPT_relPol_POL$nOcc_POL_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes labels 
#mtext('Number of records (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$nOcc_EQ_log
x_vals_POL <- s_meanPPT_relPol_POL$nOcc_POL_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log)

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol)
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol)
summary(lin_mod_meanPPT_EQ)
summary(lin_mod_meanPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'nOcc_EQ_log',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'nOcc_POL_log',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_meanPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_meanPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



################
### Body mass ###
################



#subset EQ and POL
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol,
                 s_meanPPT_relPol$bodyMass_log), ]
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol,
                 s_meanPPT_relPol$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_meanPPT_relPol_EQ$bodyMass_log,
                 s_meanPPT_relPol_POL$bodyMass_log))
y_lim <- c(-2, 2)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_meanPPT_relPol_EQ$bodyMass_log,
     s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_eq, "20"),
     ylim = y_lim,
     xlim = x_lim)

#addPOLpts
points(s_meanPPT_relPol_POL$bodyMass_log,
       s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol,
       pch = 19,
       col = paste0(col_pol, "20"),
       cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_meanPPT_relPol_EQ$bodyMass_log,
                    s_meanPPT_relPol_POL$bodyMass_log)) -
  min(c(s_meanPPT_relPol_EQ$bodyMass_log,
        s_meanPPT_relPol_POL$bodyMass_log))

plotting_positions <- c(min(c(s_meanPPT_relPol_EQ$bodyMass_log,
                              s_meanPPT_relPol_POL$bodyMass_log)),
                        min(c(s_meanPPT_relPol_EQ$bodyMass_log,
                              s_meanPPT_relPol_POL$bodyMass_log)) + range_vals/4,
                        min(c(s_meanPPT_relPol_EQ$bodyMass_log,
                              s_meanPPT_relPol_POL$bodyMass_log)) + range_vals/2,
                        min(c(s_meanPPT_relPol_EQ$bodyMass_log,
                              s_meanPPT_relPol_POL$bodyMass_log)) + range_vals/4*3,
                        max(c(s_meanPPT_relPol_EQ$bodyMass_log,
                              s_meanPPT_relPol_POL$bodyMass_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)

axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Body mass (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_meanPPT_relPol_EQ$bodyMass_log
x_vals_POL <- s_meanPPT_relPol_POL$bodyMass_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_meanPPT_EQ <- lm(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol ~ x_vals_EQ,
                         weights = s_meanPPT_relPol_EQ$nOcc_EQ_log)

lin_mod_meanPPT_POL <- lm(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol ~ x_vals_POL,
                          weights = s_meanPPT_relPol_POL$nOcc_POL_log)

length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol)
length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol)
summary(lin_mod_meanPPT_EQ)
summary(lin_mod_meanPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT_EQ,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'bodyMass_log',
                        n = length(s_meanPPT_relPol_EQ$slope_EQ_meanPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT_POL,
                        response = 'meanPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'bodyMass_log',
                        n = length(s_meanPPT_relPol_POL$slope_POL_meanPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_meanPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_meanPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800




###### maxPPT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_maxPPT_relPol <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,] #corfilter
s_maxPPT_relPol <- s_maxPPT_relPol[complete.cases(s_maxPPT_relPol$Cor_vars_maxPPT),] #corcomplete

#we do not filter by slope here (we subset EQ and POL later for plotting)

#create new columns with log values (boruta does not accept it...)
s_maxPPT_relPol$rangeSize_log10 <- log(s_maxPPT_relPol$rangeSize, 10) #lograngesize
s_maxPPT_relPol$nOcc_log <- log(s_maxPPT_relPol$nOcc) #lognocc
s_maxPPT_relPol$nOcc_EQ_log <- log(s_maxPPT_relPol$nOcc_EQ) #lognoccEQ
s_maxPPT_relPol$nOcc_POL_log <- log(s_maxPPT_relPol$nOcc_Pol) #lognoccPOL
s_maxPPT_relPol$bodyMass_log <- log(s_maxPPT_relPol$bodyMass) #logbodymass

#plot variables against slope



################
## Range size ##
################



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$rangeSize_log10, s_maxPPT_relPol_POL$rangeSize_log10)) #xlim
y_lim <- c(-5, 5) #ylim (same as minT/meanT for consistency)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$rangeSize_log10, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$rangeSize_log10,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(s_maxPPT_relPol$rangeSize_log10) -
  min(s_maxPPT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_maxPPT_relPol$rangeSize_log10),
                        min(s_maxPPT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_maxPPT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_maxPPT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_maxPPT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1] + diff(y_lim) * 0.001, cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes labels
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$rangeSize_log10
x_vals_POL <- s_maxPPT_relPol_POL$rangeSize_log10
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log)

lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log)

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol)
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol)
summary(lin_mod_maxPPT_EQ)
summary(lin_mod_maxPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeSize_log10',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeSize_log10',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_maxPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_maxPPT_POL,
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



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$latAmplitude,
                 s_maxPPT_relPol_POL$latAmplitude)) #xlim
y_lim <- c(-3, 3) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$latAmplitude, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$latAmplitude,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

#add axes
axis(1,
     at = pretty(x_lim),
     pos = y_lim[1],
     cex.axis = 2) #xaxis
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$latAmplitude #xvalsEQ
x_vals_POL <- s_maxPPT_relPol_POL$latAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol) #nEQ
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol) #nPOL
summary(lin_mod_maxPPT_EQ) #sumEQ
summary(lin_mod_maxPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'latAmplitude',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'latAmplitude',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_maxPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$rangeLoc, s_maxPPT_relPol_POL$rangeLoc)) #xlim
y_lim <- c(-5, 5) #ylim  (adjust if needed)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$rangeLoc, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$rangeLoc,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
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
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5) #ylab

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$rangeLoc #xvalsEQ
x_vals_POL <- s_maxPPT_relPol_POL$rangeLoc #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol) #nEQ
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol) #nPOL
summary(lin_mod_maxPPT_EQ) #sumEQ
summary(lin_mod_maxPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'rangeLoc',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'rangeLoc',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(x_range) + range_val/100,
             max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_maxPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



########################
### Elevation median ###
########################



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$elevMedian,
                 s_maxPPT_relPol_POL$elevMedian)) #xlim
y_lim <- c(-2, 2) #ylim  # <-- keep or change manually for maxPPT

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$elevMedian, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$elevMedian,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

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
x_vals_EQ <- s_maxPPT_relPol_EQ$elevMedian #xvalsEQ
x_vals_POL <- s_maxPPT_relPol_POL$elevMedian #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol) #nEQ
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol) #nPOL
summary(lin_mod_maxPPT_EQ) #sumEQ
summary(lin_mod_maxPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevMedian',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevMedian',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_maxPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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
### Elevation amplitude  ###
############################



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$elevAmplitude,
                 s_maxPPT_relPol_POL$elevAmplitude)) #xlim
y_lim <- c(-3, 3) #ylim

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$elevAmplitude, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$elevAmplitude,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
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
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$elevAmplitude #xvalsEQ
x_vals_POL <- s_maxPPT_relPol_POL$elevAmplitude #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol) #nEQ
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol) #nPOL
summary(lin_mod_maxPPT_EQ) #sumEQ
summary(lin_mod_maxPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'elevAmplitude',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'elevAmplitude',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_maxPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ] #subsetEQ
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ] #subsetPOL

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$roundness,
                 s_maxPPT_relPol_POL$roundness)) #xlim
y_lim <- c(-2, 2) #ylim  # <-- set manually for this panel

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 #x0

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$roundness, s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$roundness,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
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
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$roundness #xvalsEQ
x_vals_POL <- s_maxPPT_relPol_POL$roundness #xvalsPOL
x_range <- range(c(x_vals_EQ, x_vals_POL)) #xrange
range_val <- x_range[2] - x_range[1] #xrangeval

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log) #lmEQ
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log) #lmPOL

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol) #nEQ
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol) #nPOL
summary(lin_mod_maxPPT_EQ) #sumEQ
summary(lin_mod_maxPPT_POL) #sumPOL

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'roundness',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'roundness',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_range) - range_val/100,
             length.out = 100) #xseq
y_pred_EQ <- predict(lin_mod_maxPPT_EQ, newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence") #predEQ
y_pred_POL <- predict(lin_mod_maxPPT_POL, newdata = data.frame(x_vals_POL = x_seq),
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



#########################
### Number of records ###
#########################



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ]
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                 s_maxPPT_relPol_POL$nOcc_POL_log))
y_lim <- c(-3, 3)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$nOcc_EQ_log,
     s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_eq, "20"),
     ylim = y_lim, xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$nOcc_POL_log,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
       pch = 19, col = paste0(col_pol, "20"), cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                    s_maxPPT_relPol_POL$nOcc_POL_log)) -
  min(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
        s_maxPPT_relPol_POL$nOcc_POL_log))

plotting_positions <- c(min(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                              s_maxPPT_relPol_POL$nOcc_POL_log)),
                        min(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                              s_maxPPT_relPol_POL$nOcc_POL_log)) + range_vals/4,
                        min(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                              s_maxPPT_relPol_POL$nOcc_POL_log)) + range_vals/2,
                        min(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                              s_maxPPT_relPol_POL$nOcc_POL_log)) + range_vals/4*3,
                        max(c(s_maxPPT_relPol_EQ$nOcc_EQ_log,
                              s_maxPPT_relPol_POL$nOcc_POL_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2, cex.axis = 2)

#add axes labels 
#mtext('Number of records (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$nOcc_EQ_log
x_vals_POL <- s_maxPPT_relPol_POL$nOcc_POL_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log)
lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log)

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol)
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol)
summary(lin_mod_maxPPT_EQ)
summary(lin_mod_maxPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'nOcc_EQ_log',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results, lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'nOcc_POL_log',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_maxPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_maxPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800



################
### Body mass ###
################



#subset EQ and POL
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol,
                 s_maxPPT_relPol$bodyMass_log), ]
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol,
                 s_maxPPT_relPol$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(c(s_maxPPT_relPol_EQ$bodyMass_log,
                 s_maxPPT_relPol_POL$bodyMass_log))
y_lim <- c(-2, 2)

#plot graph with points separated by EQ-centre and POL-centre slopes
plot(s_maxPPT_relPol_EQ$bodyMass_log,
     s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_eq, "20"),
     ylim = y_lim,
     xlim = x_lim)

#addPOLpts
points(s_maxPPT_relPol_POL$bodyMass_log,
       s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol,
       pch = 19,
       col = paste0(col_pol, "20"),
       cex = 1.5)

# Define original values and their log-transformed positions
range_vals <- max(c(s_maxPPT_relPol_EQ$bodyMass_log,
                    s_maxPPT_relPol_POL$bodyMass_log)) -
  min(c(s_maxPPT_relPol_EQ$bodyMass_log,
        s_maxPPT_relPol_POL$bodyMass_log))

plotting_positions <- c(min(c(s_maxPPT_relPol_EQ$bodyMass_log,
                              s_maxPPT_relPol_POL$bodyMass_log)),
                        min(c(s_maxPPT_relPol_EQ$bodyMass_log,
                              s_maxPPT_relPol_POL$bodyMass_log)) + range_vals/4,
                        min(c(s_maxPPT_relPol_EQ$bodyMass_log,
                              s_maxPPT_relPol_POL$bodyMass_log)) + range_vals/2,
                        min(c(s_maxPPT_relPol_EQ$bodyMass_log,
                              s_maxPPT_relPol_POL$bodyMass_log)) + range_vals/4*3,
                        max(c(s_maxPPT_relPol_EQ$bodyMass_log,
                              s_maxPPT_relPol_POL$bodyMass_log)))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)

axis(2,
     at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1],
     las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Body mass (log)', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals_EQ <- s_maxPPT_relPol_EQ$bodyMass_log
x_vals_POL <- s_maxPPT_relPol_POL$bodyMass_log
x_range <- range(c(x_vals_EQ, x_vals_POL))
range_val <- x_range[2] - x_range[1]

#fit linear models EQ and POL
lin_mod_maxPPT_EQ <- lm(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol ~ x_vals_EQ,
                        weights = s_maxPPT_relPol_EQ$nOcc_EQ_log)

lin_mod_maxPPT_POL <- lm(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol ~ x_vals_POL,
                         weights = s_maxPPT_relPol_POL$nOcc_POL_log)

length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol)
length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol)
summary(lin_mod_maxPPT_EQ)
summary(lin_mod_maxPPT_POL)

#store regression results (EQ portion)
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT_EQ,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'EQ',
                        predictor = 'bodyMass_log',
                        n = length(s_maxPPT_relPol_EQ$slope_EQ_maxPPT_relPol),
                        weight_variable = 'nOcc_EQ_log')

#store regression results (POL portion)
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT_POL,
                        response = 'maxPPT',
                        gradient = 'relPol',
                        portion = 'POL',
                        predictor = 'bodyMass_log',
                        n = length(s_maxPPT_relPol_POL$slope_POL_maxPPT_relPol),
                        weight_variable = 'nOcc_POL_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred_EQ <- predict(lin_mod_maxPPT_EQ,
                     newdata = data.frame(x_vals_EQ = x_seq),
                     interval = "confidence")

y_pred_POL <- predict(lin_mod_maxPPT_POL,
                      newdata = data.frame(x_vals_POL = x_seq),
                      interval = "confidence")

#add regression lines
lines(x_seq, y_pred_EQ[, "fit"], col = col_eq, lwd = 6)
lines(x_seq, y_pred_POL[, "fit"], col = col_pol, lwd = 6)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred_EQ[, "lwr"], rev(y_pred_EQ[, "upr"])),
        col = paste0(col_eq, "30"),
        border = NA)

polygon(c(x_seq, rev(x_seq)),
        c(y_pred_POL[, "lwr"], rev(y_pred_POL[, "upr"])),
        col = paste0(col_pol, "30"),
        border = NA)

#save 800




###### MinPPT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_minPPT_distEdge <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7, ]
s_minPPT_distEdge <- s_minPPT_distEdge[complete.cases(s_minPPT_distEdge$Cor_vars_minPPT), ]
s_minPPT_distEdge <- s_minPPT_distEdge[complete.cases(s_minPPT_distEdge$slope_minPPT_distEdge), ]

#create new columns with log values (boruta does not accept it...)
s_minPPT_distEdge$rangeSize_log10 <- log(s_minPPT_distEdge$rangeSize, 10)
s_minPPT_distEdge$nOcc_log <- log(s_minPPT_distEdge$nOcc)
s_minPPT_distEdge$bodyMass_log <- log(s_minPPT_distEdge$bodyMass)

#define colour
col_all <- "#303030"


#plot variables against slope



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$rangeSize_log10, s_minPPT_distEdge$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minPPT_distEdge$rangeSize_log10) -
  min(s_minPPT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_minPPT_distEdge$rangeSize_log10),
                        min(s_minPPT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_minPPT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_minPPT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_minPPT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)

axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeSize_log10',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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
x_lim <- range(s_minPPT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$latAmplitude, s_minPPT_distEdge$slope_minPPT_distEdge,
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
x_vals <- s_minPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'latAmplitude',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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



############################
### Latitudinal position ###
############################



# Define x and y limits
x_lim <- range(s_minPPT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0

#plot graph
plot(s_minPPT_distEdge$rangeLoc, s_minPPT_distEdge$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeLoc',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
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
x_lim <- range(s_minPPT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$elevMedian, s_minPPT_distEdge$slope_minPPT_distEdge,
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
x_vals <- s_minPPT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevMedian',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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
x_lim <- range(s_minPPT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$elevAmplitude, s_minPPT_distEdge$slope_minPPT_distEdge,
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
x_vals <- s_minPPT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevAmplitude',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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
x_lim <- range(s_minPPT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$roundness, s_minPPT_distEdge$slope_minPPT_distEdge,
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
x_vals <- s_minPPT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'roundness',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_distEdge$nOcc_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge$nOcc_log, s_minPPT_distEdge$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_minPPT_distEdge$nOcc_log) -
  min(s_minPPT_distEdge$nOcc_log)


plotting_positions <- c(min(s_minPPT_distEdge$nOcc_log),
                        min(s_minPPT_distEdge$nOcc_log) + range_vals/4,
                        min(s_minPPT_distEdge$nOcc_log) + range_vals/2,
                        min(s_minPPT_distEdge$nOcc_log) + range_vals/4*3,
                        max(s_minPPT_distEdge$nOcc_log))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las = 2, cex.axis = 2)


#add axes lables
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge$nOcc_log)

length(s_minPPT_distEdge$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'nOcc_log',
                        n = length(s_minPPT_distEdge$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)


#save 800



#################
### Body mass ###
#################



#subset body mass
s_minPPT_distEdge_bodyMass <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$slope_minPPT_distEdge,
                 s_minPPT_distEdge$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_distEdge_bodyMass$bodyMass_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minPPT_distEdge_bodyMass$bodyMass_log,
     s_minPPT_distEdge_bodyMass$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minPPT_distEdge_bodyMass$bodyMass_log) -
  min(s_minPPT_distEdge_bodyMass$bodyMass_log)

plotting_positions <- c(min(s_minPPT_distEdge_bodyMass$bodyMass_log),
                        min(s_minPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4,
                        min(s_minPPT_distEdge_bodyMass$bodyMass_log) + range_vals/2,
                        min(s_minPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4*3,
                        max(s_minPPT_distEdge_bodyMass$bodyMass_log))

plotting_values <- round(exp(plotting_positions), 0)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge_bodyMass$bodyMass_log
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge_bodyMass$slope_minPPT_distEdge ~ x_vals,
                     weights = s_minPPT_distEdge_bodyMass$nOcc_log)

length(s_minPPT_distEdge_bodyMass$slope_minPPT_distEdge)
summary(lin_mod_minPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_minPPT,
                        response = 'minPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'bodyMass_log',
                        n = length(s_minPPT_distEdge_bodyMass$slope_minPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_minPPT,
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




###### MeanPPT vs DISTANCE TO RANGE EDGE  ######



#select only species that had lower correl between vars
s_meanPPT_distEdge <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_distEdge <- s_meanPPT_distEdge[complete.cases(s_meanPPT_distEdge$Cor_vars_meanPPT),]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$slope_meanPPT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_meanPPT_distEdge$rangeSize_log10 <- log(s_meanPPT_distEdge$rangeSize, 10)
s_meanPPT_distEdge$nOcc_log <- log(s_meanPPT_distEdge$nOcc)
s_meanPPT_distEdge$bodyMass_log <- log(s_meanPPT_distEdge$bodyMass)



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)  # <-- set manually for meanPPT if needed

#plot graph
plot(s_meanPPT_distEdge$rangeSize_log10, s_meanPPT_distEdge$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_meanPPT_distEdge$rangeSize_log10) -
  min(s_meanPPT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_meanPPT_distEdge$rangeSize_log10),
                        min(s_meanPPT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_meanPPT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_meanPPT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_meanPPT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeSize_log10',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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
x_lim <- range(s_meanPPT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge$latAmplitude, s_meanPPT_distEdge$slope_meanPPT_distEdge,
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
x_vals <- s_meanPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'latAmplitude',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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
x_lim <- range(s_meanPPT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0

#plot graph
plot(s_meanPPT_distEdge$rangeLoc, s_meanPPT_distEdge$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables
#mtext('Latitudinal position', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeLoc',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
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
x_lim <- range(s_meanPPT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge$elevMedian, s_meanPPT_distEdge$slope_meanPPT_distEdge,
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
x_vals <- s_meanPPT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevMedian',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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
x_lim <- range(s_meanPPT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge$elevAmplitude, s_meanPPT_distEdge$slope_meanPPT_distEdge,
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
x_vals <- s_meanPPT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevAmplitude',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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
x_lim <- range(s_meanPPT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge$roundness, s_meanPPT_distEdge$slope_meanPPT_distEdge,
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
x_vals <- s_meanPPT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'roundness',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_distEdge$nOcc_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge$nOcc_log, s_meanPPT_distEdge$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_meanPPT_distEdge$nOcc_log) -
  min(s_meanPPT_distEdge$nOcc_log)


plotting_positions <- c(min(s_meanPPT_distEdge$nOcc_log),
                        min(s_meanPPT_distEdge$nOcc_log) + range_vals/4,
                        min(s_meanPPT_distEdge$nOcc_log) + range_vals/2,
                        min(s_meanPPT_distEdge$nOcc_log) + range_vals/4*3,
                        max(s_meanPPT_distEdge$nOcc_log))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las = 2, cex.axis = 2)


#add axes lables
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge$nOcc_log)

length(s_meanPPT_distEdge$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'nOcc_log',
                        n = length(s_meanPPT_distEdge$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)



#save 800


#################
### Body mass ###
#################



#subset body mass
s_meanPPT_distEdge_bodyMass <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$slope_meanPPT_distEdge,
                 s_meanPPT_distEdge$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_distEdge_bodyMass$bodyMass_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanPPT_distEdge_bodyMass$bodyMass_log,
     s_meanPPT_distEdge_bodyMass$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_meanPPT_distEdge_bodyMass$bodyMass_log) -
  min(s_meanPPT_distEdge_bodyMass$bodyMass_log)

plotting_positions <- c(min(s_meanPPT_distEdge_bodyMass$bodyMass_log),
                        min(s_meanPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4,
                        min(s_meanPPT_distEdge_bodyMass$bodyMass_log) + range_vals/2,
                        min(s_meanPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4*3,
                        max(s_meanPPT_distEdge_bodyMass$bodyMass_log))

plotting_values <- round(exp(plotting_positions), 0)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge_bodyMass$bodyMass_log
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge_bodyMass$slope_meanPPT_distEdge ~ x_vals,
                      weights = s_meanPPT_distEdge_bodyMass$nOcc_log)

length(s_meanPPT_distEdge_bodyMass$slope_meanPPT_distEdge)
summary(lin_mod_meanPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_meanPPT,
                        response = 'meanPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'bodyMass_log',
                        n = length(s_meanPPT_distEdge_bodyMass$slope_meanPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_meanPPT,
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



###### MaxPPT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_maxPPT_distEdge <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_distEdge <- s_maxPPT_distEdge[complete.cases(s_maxPPT_distEdge$Cor_vars_maxPPT),]
s_maxPPT_distEdge <- s_maxPPT_distEdge[complete.cases(s_maxPPT_distEdge$slope_maxPPT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_maxPPT_distEdge$rangeSize_log10 <- log(s_maxPPT_distEdge$rangeSize, 10)
s_maxPPT_distEdge$nOcc_log <- log(s_maxPPT_distEdge$nOcc)
s_maxPPT_distEdge$bodyMass_log <- log(s_maxPPT_distEdge$bodyMass)



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.02, 0.02)  # <-- set manually for maxPPT if needed

#plot graph
plot(s_maxPPT_distEdge$rangeSize_log10, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_maxPPT_distEdge$rangeSize_log10) -
  min(s_maxPPT_distEdge$rangeSize_log10)

plotting_positions <- c(min(s_maxPPT_distEdge$rangeSize_log10),
                        min(s_maxPPT_distEdge$rangeSize_log10) + range_vals/4,
                        min(s_maxPPT_distEdge$rangeSize_log10) + range_vals/2,
                        min(s_maxPPT_distEdge$rangeSize_log10) + range_vals/4*3,
                        max(s_maxPPT_distEdge$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values, pos = y_lim[1],
     cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5), pos = x_lim[1], las = 2,
     cex.axis = 2)

#add axes lables 
mtext('Range size', side = 1, line = 5.6, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeSize_log10',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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
x_lim <- range(s_maxPPT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$latAmplitude, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
x_vals <- s_maxPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'latAmplitude',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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
x_lim <- range(s_maxPPT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$rangeLoc, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
x_vals <- s_maxPPT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'rangeLoc',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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
x_lim <- range(s_maxPPT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$elevMedian, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
x_vals <- s_maxPPT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevMedian',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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
x_lim <- range(s_maxPPT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$elevAmplitude, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
mtext('Elevation amplitude', side = 1, line = 5.6, cex = 2.5)  
# mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'elevAmplitude',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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
x_lim <- range(s_maxPPT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$roundness, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
x_vals <- s_maxPPT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxT,
                        response = 'maxT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'roundness',
                        n = length(s_maxT_distEdge$slope_maxT_distEdge),
                        weight_variable = 'nOcc_log')

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
x_lim <- range(s_maxPPT_distEdge$roundness)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$roundness, s_maxPPT_distEdge$slope_maxPPT_distEdge,
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
x_vals <- s_maxPPT_distEdge$roundness
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'roundness',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(x_vals) + range_val/100,
             max(x_vals) - range_val/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$nOcc_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge$nOcc_log, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_maxPPT_distEdge$nOcc_log) -
  min(s_maxPPT_distEdge$nOcc_log)


plotting_positions <- c(min(s_maxPPT_distEdge$nOcc_log),
                        min(s_maxPPT_distEdge$nOcc_log) + range_vals/4,
                        min(s_maxPPT_distEdge$nOcc_log) + range_vals/2,
                        min(s_maxPPT_distEdge$nOcc_log) + range_vals/4*3,
                        max(s_maxPPT_distEdge$nOcc_log))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las = 2, cex.axis = 2)


#add axes lables
mtext('Number of records', side = 1, line = 5.6, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge$nOcc_log)

length(s_maxPPT_distEdge$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'nOcc_log',
                        n = length(s_maxPPT_distEdge$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = col_all, lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = paste0(col_all, "30"),
        border = NA)



#save 800



#################
### Body mass ###
#################



#subset body mass
s_maxPPT_distEdge_bodyMass <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$slope_maxPPT_distEdge,
                 s_maxPPT_distEdge$bodyMass_log), ]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge_bodyMass$bodyMass_log)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxPPT_distEdge_bodyMass$bodyMass_log,
     s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     cex = 1.5,
     pch = 19,
     col = paste0(col_all, "20"),
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_maxPPT_distEdge_bodyMass$bodyMass_log) -
  min(s_maxPPT_distEdge_bodyMass$bodyMass_log)

plotting_positions <- c(min(s_maxPPT_distEdge_bodyMass$bodyMass_log),
                        min(s_maxPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4,
                        min(s_maxPPT_distEdge_bodyMass$bodyMass_log) + range_vals/2,
                        min(s_maxPPT_distEdge_bodyMass$bodyMass_log) + range_vals/4*3,
                        max(s_maxPPT_distEdge_bodyMass$bodyMass_log))

plotting_values <- round(exp(plotting_positions), 0)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, at = seq(y_lim[1], y_lim[2], length.out = 5),
     pos = x_lim[1], las = 2, cex.axis = 2)

#add axes lables 
mtext('Body mass', side = 1, line = 5.6, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge_bodyMass$bodyMass_log
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge ~ x_vals,
                     weights = s_maxPPT_distEdge_bodyMass$nOcc_log)

length(s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge)
summary(lin_mod_maxPPT)

#store regression results
reg_results <- store_lm(reg_results,
                        lin_mod_maxPPT,
                        response = 'maxPPT',
                        gradient = 'distEdge',
                        portion = 'all',
                        predictor = 'bodyMass_log',
                        n = length(s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge),
                        weight_variable = 'nOcc_log')

#predict y-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)

y_pred <- predict(lin_mod_maxPPT,
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


#save table with stats
setwd('/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SI/Figures/Figure_S5')

write.csv(reg_results, 'Results_lm_precipitation_slopes.csv', row.names = F)



