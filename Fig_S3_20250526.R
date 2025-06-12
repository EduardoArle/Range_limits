#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('20250531_Slopes.csv')

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

###### MinT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_minT_relPol <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$Cor_vars_minT),]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$slope_minT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_minT_relPol$rangeSize_log10 <- log(s_minT_relPol$rangeSize, 10)
s_minT_relPol$nOcc_log <- log(s_minT_relPol$nOcc)

#plot meaningful variables against slope

################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_minT_relPol$rangeSize_log10, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

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
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)



########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$latAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)



#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)



#save 800





###### MeanT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_meanT_relPol <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$Cor_vars_meanT),]
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$slope_meanT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_meanT_relPol$rangeSize_log10 <- log(s_meanT_relPol$rangeSize, 10)
s_meanT_relPol$nOcc_log <- log(s_meanT_relPol$nOcc)

#plot meaningful variables against slope


################
## Range size ##
################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10) #use minT to keep same values x
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanT_relPol$rangeSize_log10, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minnT_relPol$rangeSize_log10) -
  min(s_minnT_relPol$rangeSize_log10)

#use minT to keep same values x
plotting_positions <- c(min(s_minT_relPol$rangeSize_log10),
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_minT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)





########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$elevMedian, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$latAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$latAmplitude, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)





#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$elevAmplitude, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)


#save 800


#################
### Roundness ###
#################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanT_relPol$roundness, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00805020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#008050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00805030',
        border = NA)



#save 800


#################
### Body mass ###
#################

#select only rows with bodyMass
s_meanT_relPol_bodyMass <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanT_relPol_bodyMass$bodyMass, s_meanT_relPol_bodyMass$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00505020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol_bodyMass$slope_meanT_relPol ~ x_vals,
                    weights = s_meanT_relPol_bodyMass$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#005050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00505030',
        border = NA)


#save 800






###### MaxT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_maxT_relPol <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$Cor_vars_maxT),]
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$slope_maxT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_maxT_relPol$rangeSize_log10 <- log(s_maxT_relPol$rangeSize, 10)
s_maxT_relPol$nOcc_log <- log(s_maxT_relPol$nOcc)

#plot meaningful variables against slope


########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$elevMedian, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$elevAmplitude, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)


#save 800


#################
### Roundness ###
#################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_maxT_relPol$roundness, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00805020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#008050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00805030',
        border = NA)



#save 800



#################
### Body mass ###
#################

#select only rows with bodyMass
s_maxT_relPol_bodyMass <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_maxT_relPol_bodyMass$bodyMass, s_maxT_relPol_bodyMass$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00505020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol_bodyMass$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol_bodyMass$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#005050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00505030',
        border = NA)



#save 800



###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$nOcc_log)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_maxT_relPol$nOcc_log, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90809020',
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_maxT_relPol$nOcc_log) -
  min(s_maxT_relPol$nOcc_log)


plotting_positions <- c(min(s_maxT_relPol$nOcc_log),
                        min(s_maxT_relPol$nOcc_log) + range_vals/4,
                        min(s_maxT_relPol$nOcc_log) + range_vals/2,
                        min(s_maxT_relPol$nOcc_log) + range_vals/4*3,
                        max(s_maxT_relPol$nOcc_log))
                        

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#908090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90809030',
        border = NA)



#save 800



###### MinT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_minT_distEdge <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_distEdge <- s_minT_distEdge[
  complete.cases(s_minT_distEdge$Cor_vars_minT),]
s_minT_distEdge <- s_minT_distEdge[
  complete.cases(s_minT_distEdge$slope_minT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_minT_distEdge$rangeSize_log10 <- log(s_minT_distEdge$rangeSize, 10)
s_minT_distEdge$nOcc_log <- log(s_minT_distEdge$nOcc)

#plot meaningful variables against slope

################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_minT_distEdge$rangeSize_log10, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305020',
        border = NA)




########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$elevMedian)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$elevMedian, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$latAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)




#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, max(s_minT_distEdge$elevAmplitude))
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$elevAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)



#save 800


###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$nOcc_log)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_minT_distEdge$nOcc_log, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90809020',
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_minT_distEdge$nOcc_log) -
  min(s_minT_distEdge$nOcc_log)


plotting_positions <- c(min(s_minT_distEdge$nOcc_log),
                        min(s_minT_distEdge$nOcc_log) + range_vals/4,
                        min(s_minT_distEdge$nOcc_log) + range_vals/2,
                        min(s_minT_distEdge$nOcc_log) + range_vals/4*3,
                        max(s_minT_distEdge$nOcc_log))


plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#908090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90809030',
        border = NA)



#save 800




###### MeanT vs DISTANCE TO RANGE EDGE ######

#select only species that had lower correl between vars
s_meanT_distEdge <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$Cor_vars_meanT),]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$slope_meanT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_meanT_distEdge$rangeSize_log10 <- log(s_meanT_distEdge$rangeSize, 10)
s_meanT_distEdge$nOcc_log <- log(s_meanT_distEdge$nOcc)

#plot meaningful variables against slope

################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_meanT_distEdge$rangeSize_log10, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)


#save 800


########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$elevMedian)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$elevMedian, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$latAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)


#save 800



#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, max(s_meanT_distEdge$elevAmplitude))
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$elevAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)



#save 800


###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$nOcc_log)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_meanT_distEdge$nOcc_log, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90809020',
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_meanT_distEdge$nOcc_log) -
  min(s_meanT_distEdge$nOcc_log)


plotting_positions <- c(min(s_meanT_distEdge$nOcc_log),
                        min(s_meanT_distEdge$nOcc_log) + range_vals/4,
                        min(s_meanT_distEdge$nOcc_log) + range_vals/2,
                        min(s_meanT_distEdge$nOcc_log) + range_vals/4*3,
                        max(s_meanT_distEdge$nOcc_log))

plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#908090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90809030',
        border = NA)



#save 800




###### MaxT vs DISTANCE EDGE ######

#select only species that had lower correl between vars
s_maxT_distEdge <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_distEdge <- s_maxT_distEdge[
  complete.cases(s_maxT_distEdge$Cor_vars_maxT),]
s_maxT_distEdge <- s_maxT_distEdge[
  complete.cases(s_maxT_distEdge$slope_maxT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_maxT_distEdge$rangeSize_log10 <- log(s_maxT_distEdge$rangeSize, 10)
s_maxT_distEdge$nOcc_log <- log(s_maxT_distEdge$nOcc)

#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_maxT_distEdge$rangeSize_log10, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                    weights = s_maxT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_distEdge$latAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                    weights = s_maxT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)



#save 800


########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$elevMedian)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_distEdge$elevMedian, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



###### MinPPT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_minPPT_relPol <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_relPol <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$Cor_vars_minPPT),]
s_minPPT_relPol <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_minPPT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_minPPT_relPol$rangeSize_log10 <- log(s_minPPT_relPol$rangeSize, 10)
s_minPPT_relPol$nOcc_log <- log(s_minPPT_relPol$nOcc)

#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_relPol$rangeSize_log10)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_minPPT_relPol$rangeSize_log10, s_minPPT_relPol$slope_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

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
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_relPol$slope_minPPT_relPol ~ x_vals,
                   weights = s_minPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)



########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minPPT_relPol$elevMedian, s_minPPT_relPol$slope_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_relPol$slope_minPPT_relPol ~ x_vals,
                   weights = s_minPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_relPol$latAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minPPT_relPol$latAmplitude, s_minPPT_relPol$slope_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_relPol$slope_minPPT_relPol ~ x_vals,
                   weights = s_minPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)



#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minPPT_relPol$elevAmplitude, s_minPPT_relPol$slope_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_relPol$slope_minPPT_relPol ~ x_vals,
                   weights = s_minPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)



#save 800



###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_relPol$nOcc_log)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_minPPT_relPol$nOcc_log, s_minPPT_relPol$slope_minPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90809020',
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_minPPT_relPol$nOcc_log) -
  min(s_minPPT_relPol$nOcc_log)


plotting_positions <- c(min(s_minPPT_relPol$nOcc_log),
                        min(s_minPPT_relPol$nOcc_log) + range_vals/4,
                        min(s_minPPT_relPol$nOcc_log) + range_vals/2,
                        min(s_minPPT_relPol$nOcc_log) + range_vals/4*3,
                        max(s_minPPT_relPol$nOcc_log))


plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_relPol$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_relPol$slope_minPPT_relPol ~ x_vals,
                   weights = s_minPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#908090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90809030',
        border = NA)




###### MeanPPT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_meanPPT_relPol <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_relPol <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$Cor_vars_meanPPT),]
s_meanPPT_relPol <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_meanPPT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_meanPPT_relPol$rangeSize_log10 <- log(s_meanPPT_relPol$rangeSize, 10)
s_meanPPT_relPol$nOcc_log <- log(s_meanPPT_relPol$nOcc)

#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_relPol$rangeSize_log10)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanPPT_relPol$rangeSize_log10, s_meanPPT_relPol$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_meanPPT_relPol$rangeSize_log10) -
  min(s_meanPPT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_meanPPT_relPol$rangeSize_log10),
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_meanPPT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_meanPPT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)





########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanPPT_relPol$elevMedian, s_meanPPT_relPol$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_relPol$latAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanPPT_relPol$latAmplitude, s_meanPPT_relPol$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)





#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanPPT_relPol$elevAmplitude, s_meanPPT_relPol$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)


#save 800


#################
### Roundness ###
#################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanPPT_relPol$roundness, s_meanPPT_relPol$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00805020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#008050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00805030',
        border = NA)



#save 800


#################
### Body mass ###
#################

#select only rows with bodyMass
s_meanPPT_relPol_bodyMass <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_meanPPT_relPol_bodyMass$bodyMass, s_meanPPT_relPol_bodyMass$slope_meanPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00505020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_relPol_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_relPol_bodyMass$slope_meanPPT_relPol ~ x_vals,
                    weights = s_meanPPT_relPol_bodyMass$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#005050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00505030',
        border = NA)


#save 800





###### MaxPPT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_maxPPT_relPol <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_relPol <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$Cor_vars_maxPPT),]
s_maxPPT_relPol <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_maxPPT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_maxPPT_relPol$rangeSize_log10 <- log(s_maxPPT_relPol$rangeSize, 10)
s_maxPPT_relPol$nOcc_log <- log(s_maxPPT_relPol$nOcc)

#plot meaningful variables against slope


########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_relPol$elevMedian)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxPPT_relPol$elevMedian, s_maxPPT_relPol$slope_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_relPol$slope_maxPPT_relPol ~ x_vals,
                   weights = s_maxPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800



#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_relPol$elevAmplitude)
y_lim <- c(-0.2, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxPPT_relPol$elevAmplitude, s_maxPPT_relPol$slope_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_relPol$slope_maxPPT_relPol ~ x_vals,
                   weights = s_maxPPT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)


#save 800





#################
### Body mass ###
#################

#select only rows with bodyMass
s_maxPPT_relPol_bodyMass <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-0.2, 0.2)

#plot graph
plot(s_maxPPT_relPol_bodyMass$bodyMass, s_maxPPT_relPol_bodyMass$slope_maxPPT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00505020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_relPol_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_relPol_bodyMass$slope_maxPPT_relPol ~ x_vals,
                   weights = s_maxPPT_relPol_bodyMass$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#005050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00505030',
        border = NA)



#save 800




###### MinPPT vs DISTANCE TO RANGE EDGE  ######

#select only species that had lower correl between vars
s_minPPT_distEdge <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$Cor_vars_minPPT),]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$slope_minPPT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_minPPT_distEdge$rangeSize_log10 <- log(s_minPPT_distEdge$rangeSize, 10)
s_minPPT_distEdge$nOcc_log <- log(s_minPPT_distEdge$nOcc)



#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_minPPT_distEdge$rangeSize_log10, s_minPPT_distEdge$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                   weights = s_minPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minPPT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minPPT_distEdge$latAmplitude, s_minPPT_distEdge$slope_minPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minPPT <- lm(s_minPPT_distEdge$slope_minPPT_distEdge ~ x_vals,
                   weights = s_minPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)




###### MeanPPT vs DISTANCE TO RANGE EDGE ######

#select only species that had lower correl between vars
s_meanPPT_distEdge <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$Cor_vars_meanPPT),]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$slope_meanPPT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_meanPPT_distEdge$rangeSize_log10 <- log(s_meanPPT_distEdge$rangeSize, 10)
s_meanPPT_distEdge$nOcc_log <- log(s_meanPPT_distEdge$nOcc)



#plot meaningful variables against slope




################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_meanPPT_distEdge$rangeSize_log10, s_meanPPT_distEdge$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                    weights = s_meanPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)


#save 800




#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanPPT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanPPT_distEdge$latAmplitude, s_meanPPT_distEdge$slope_meanPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanPPT <- lm(s_meanPPT_distEdge$slope_meanPPT_distEdge ~ x_vals,
                    weights = s_meanPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)


#save 800





###### MaxPPT vs DISTANCE EDGE ######

#select only species that had lower correl between vars
s_maxPPT_distEdge <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$Cor_vars_maxPPT),]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$slope_maxPPT_distEdge),]

#create new columns with log values (boruta does not accept it...)
s_maxPPT_distEdge$rangeSize_log10 <- log(s_maxPPT_distEdge$rangeSize, 10)
s_maxPPT_distEdge$nOcc_log <- log(s_maxPPT_distEdge$nOcc)

#plot meaningful variables against slope



################
## Range size ##
################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$rangeSize_log10)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_maxPPT_distEdge$rangeSize_log10, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
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
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                   weights = s_maxPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$latAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxPPT_distEdge$latAmplitude, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#F0803020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                   weights = s_maxPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#F08030', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#F0803030',
        border = NA)



#save 800


########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$elevMedian)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxPPT_distEdge$elevMedian, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#99005020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxPPT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                   weights = s_maxPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxPPT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#990050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#99005030',
        border = NA)



#save 800


#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxPPT_distEdge$elevAmplitude)
y_lim <- c(-0.001, 0.001)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxPPT_distEdge$elevAmplitude, s_maxPPT_distEdge$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#0080FF20',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_maxPPT_distEdge$slope_maxPPT_distEdge ~ x_vals,
                   weights = s_maxPPT_distEdge$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#0080FF', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#0080FF30',
        border = NA)



#save 800



#################
### Body mass ###
#################

#select only rows with bodyMass
s_maxPPT_distEdge_bodyMass <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-0.001, 0.001)

#plot graph
plot(s_maxPPT_distEdge_bodyMass$bodyMass, s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#00505020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxPPT_distEdge_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_maxPPT_distEdge_bodyMass$slope_maxPPT_distEdge ~ x_vals,
                    weights = s_maxPPT_distEdge_bodyMass$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#005050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#00505030',
        border = NA)


#save 800


