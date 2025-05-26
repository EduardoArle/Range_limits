#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('Slopes.csv')

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

#restricted ylim

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10)
y_lim <- c(-0.3, 0.3)

#plot graph
plot(s_minT_relPol$rangeSize_log10, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#50305030',
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
mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#503050', lwd = 8)

#save 800


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10)
y_lim <- c(-0.3, 0.3)

#plot graph
plot(s_minT_relPol$rangeSize_log10, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#50305030',
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
mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

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
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#503050', lwd = 8)


#################
### Roundness ###
#################


#restricted ylim

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.3, 0.3)

#plot graph
plot(s_minT_relPol$roundness, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#00805030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#008050', lwd = 8)

#save 800


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.3, 0.3)

#plot graph
plot(s_minT_relPol$roundness, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#00805030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#008050', lwd = 8)

#save 800



########################
### Elevation median ###
########################


#restricted ylim

# Define x and y limits
x_lim <- range(s_minT_relPol$elevMedian)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#99005030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#990050', lwd = 8)

#save 800


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevMedian)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#99005030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

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
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#990050', lwd = 8)

#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim

# Define x and y limits
x_lim <- range(s_minT_relPol$latAmplitude)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#F0803030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#F08030', lwd = 8)

#save 800

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$latAmplitude)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#F0803030',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

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
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#F08030', lwd = 8)



#############################
### Elevational amplitude ###
#############################


#restricted ylim

# Define x and y limits
x_lim <- range(s_minT_relPol$elevAmplitude)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#0080FF30',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#0080FF', lwd = 8)

#save 800

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevAmplitude)
y_lim <- c(-0.3, 0.3)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     pch = 19, col = '#0080FF30',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope minT vs. relPol', side = 2, line = 6.5, cex = 2.5)

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
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq))

# Add regression line from the y-axis onwards
lines(x_seq, y_pred, col = '#0080FF', lwd = 8)

#save 800




###############################. PLOTS FROM BORUTA SCRIPT ###################


#plot meaningful variables against slope

## Range size

#lm without weights
plot(log(s_minT_relPol$rangeSize, 10), s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~
                     log(s_minT_relPol$rangeSize, 10))

abline(lin_mod_minT, col = '#503050', lwd = 3)

#restricted ylim
plot(log(s_minT_relPol$rangeSize, 10), s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~
                     log(s_minT_relPol$rangeSize, 10))

abline(lin_mod_minT, col = '#503050', lwd = 3)

#make lm with weights = nOcc
plot(log(s_minT_relPol$rangeSize, 10), s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = s_minT_relPol$nOcc)

abline(lin_mod_minT, col = '#503050', lwd = 3)

#restricted ylim
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = s_minT_relPol$nOcc)

abline(lin_mod_minT, col = '#503050', lwd = 3)


#make lm with weights = log(nOcc)
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = log(s_minT_relPol$nOcc))

abline(lin_mod_minT, col = '#503050', lwd = 3)

#restricted ylim
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = log(s_minT_relPol$nOcc))

abline(lin_mod_minT, col = '#503050', lwd = 3)




## Roundness

#lm without weights
plot(s_minT_relPol$roundness, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range roundness',
     pch = 19, col = '#10207030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~
                     log(s_minT_relPol$rangeSize, 10))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$roundness)

abline(lin_mod_minT, col = '#102070', lwd = 3)


#restricted ylim
plot(s_minT_relPol$roundness, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range roundness',
     pch = 19, col = '#10207030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$roundness)

abline(lin_mod_minT, col = '#102070', lwd = 3) 



#make lm with weights = nOcc
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = s_minT_relPol$nOcc)

abline(lin_mod_minT, col = '#503050', lwd = 3)

#restricted ylim
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = s_minT_relPol$nOcc)

abline(lin_mod_minT, col = '#503050', lwd = 3)


#make lm with weights = log(nOcc)
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = log(s_minT_relPol$nOcc))

abline(lin_mod_minT, col = '#503050', lwd = 3)

#restricted ylim
plot(s_minT_relPol$rangeSize, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030',
     ylim = c(-0.2, 0.2))

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$rangeSize,
                   weights = log(s_minT_relPol$nOcc))

abline(lin_mod_minT, col = '#503050', lwd = 3)





plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Median elevation',
     pch = 19, col = '#80800030')

plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Elevation amplitude',
     pch = 19, col = '#f0800030')

lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ s_minT_relPol$elevAmplitude,
                   weights = s_minT_relPol$nOcc)
abline(lin_mod_minT, col = '#f08000', lwd = 3)

plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     ylab = 'Slope minT_relPol', xlab = 'Latitudinal amplitude',
     pch = 19, col = '#f0806030')



#plot meaningful variables against slope

plot(s_meanT_relPol$rangeSize, s_meanT_relPol$slope_meanT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ s_meanT_relPol$rangeSize)
abline(lin_mod_meanT, col = '#503050', lwd = 3)


plot(s_meanT_relPol$bodyMass, s_meanT_relPol$slope_meanT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Body mass',
     pch = 19, col = '#80305030')

plot(s_meanT_relPol$elevMedian, s_meanT_relPol$slope_meanT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Median elevation',
     pch = 19, col = '#80800030')

plot(s_meanT_relPol$elevAmplitude, s_meanT_relPol$slope_meanT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Elevation amplitude',
     pch = 19, col = '#f0800030')

plot(s_meanT_relPol$latAmplitude, s_meanT_relPol$slope_meanT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Latitudinal amplitude',
     pch = 19, col = '#f0806030')



#plot meaningful variables against slope

plot(s_maxT_relPol$rangeSize, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope maxT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_maxT <- lm(s_meanT_relPol$slope_maxT_relPol ~ s_meanT_relPol$rangeSize)
abline(lin_mod_maxT, col = '#503050', lwd = 3)

plot(s_maxT_relPol$nOcc, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope maxT_relPol', xlab = 'n Occ',
     pch = 19, col = '#80305030')

plot(s_maxT_relPol$elevAmplitude, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Elevation amplitude',
     pch = 19, col = '#f0800030')

plot(s_maxT_relPol$latAmplitude, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Latitudinal amplitude',
     pch = 19, col = '#f0806030')


#### plot graph of slope vs range size

#### plot graph of slope vs latitudinal amplitude


### running only the species with a useful fir (lower 0.3)


#plot meaningful variables against slope

plot(s_maxT_relPol$rangeSize, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope maxT_relPol', xlab = 'Range size',
     pch = 19, col = '#50305030')

lin_mod_maxT <- lm(s_meanT_relPol$slope_maxT_relPol ~ s_meanT_relPol$rangeSize)
abline(lin_mod_maxT, col = '#503050', lwd = 3)

plot(s_maxT_relPol$nOcc, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope maxT_relPol', xlab = 'n Occ',
     pch = 19, col = '#80305030')

plot(s_maxT_relPol$elevAmplitude, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Elevation amplitude',
     pch = 19, col = '#f0800030')

plot(s_maxT_relPol$latAmplitude, s_maxT_relPol$slope_maxT_relPol,
     ylab = 'Slope meanT_relPol', xlab = 'Latitudinal amplitude',
     pch = 19, col = '#f0806030')
