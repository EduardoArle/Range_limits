#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('20260208_Slopes.csv')

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
y_lim <- c(-1, 1)

#plot graph
plot(s_minT_relPol$rangeSize_log10, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_relPol$slope_minT_relPol)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800




#############################
### Latitudinal amplitude ###
#############################




#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$latAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_relPol$slope_minT_relPol)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)

#save 800


########################
### Median elevation ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevMedian)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_relPol$slope_minT_relPol)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800





#############################
### Elevational amplitude ###
#############################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_relPol$slope_minT_relPol)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800


#################
### Body mass ###
#################



#select only rows with bodyMass
s_minT_relPol_bodyMass <- s_minT_relPol[
  complete.cases(s_minT_relPol$bodyMass),]

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1000)
y_lim <- c(-1, 1)

#plot graph
plot(s_minT_relPol_bodyMass$bodyMass, s_minT_relPol_bodyMass$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol_bodyMass$bodyMass
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol_bodyMass$slope_minT_relPol ~ x_vals,
                    weights = s_minT_relPol_bodyMass$nOcc_log)

length(s_minT_relPol_bodyMass$slope_minT_relPol)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
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
y_lim <- c(-1, 1)

#plot graph
plot(s_meanT_relPol$rangeSize_log10, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minT_relPol$rangeSize_log10) -
  min(s_minT_relPol$rangeSize_log10)

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

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800


#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$latAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$latAmplitude, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$rangeLoc)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$rangeLoc, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                    weights = s_meanT_relPol$nOcc_log)

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



########################
### Elevation median ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$elevMedian)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$elevMedian, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800




#############################
### Elevational amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$elevAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_relPol$elevAmplitude, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800


# #################
# ### Roundness ###
# #################
# 
# 
# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- c(0, 1)
# y_lim <- c(-0.2, 0.2)
# 
# #plot graph
# plot(s_meanT_relPol$roundness, s_meanT_relPol$slope_meanT_relPol,
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
# #mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_meanT_relPol$roundness
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
#                    weights = s_meanT_relPol$nOcc_log)
# 
# length(s_meanT_relPol$slope_meanT_relPol)
# summary(lin_mod_meanT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
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
# 
# 
# #################
# ### Body mass ###
# #################
# 
# #select only rows with bodyMass
# s_meanT_relPol_bodyMass <- s_meanT_relPol[
#   complete.cases(s_meanT_relPol$bodyMass),]
# 
# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- c(0, 1000)
# y_lim <- c(-0.2, 0.2)
# 
# #plot graph
# plot(s_meanT_relPol_bodyMass$bodyMass, s_meanT_relPol_bodyMass$slope_meanT_relPol,
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
# #mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_meanT_relPol_bodyMass$bodyMass
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_meanT <- lm(s_meanT_relPol_bodyMass$slope_meanT_relPol ~ x_vals,
#                     weights = s_meanT_relPol_bodyMass$nOcc_log)
# 
# length(s_meanT_relPol_bodyMass$slope_meanT_relPol)
# summary(lin_mod_meanT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
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




#save 800


###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_relPol$nOcc_log)
y_lim <- c(-1, 1)

#plot graph
plot(s_meanT_relPol$nOcc_log, s_meanT_relPol$slope_meanT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)


# Define original values and their log-transformed positions
range_vals <- max(s_meanT_relPol$nOcc_log) -
  min(s_meanT_relPol$nOcc_log)


plotting_positions <- c(min(s_meanT_relPol$nOcc_log),
                        min(s_meanT_relPol$nOcc_log) + range_vals/4,
                        min(s_meanT_relPol$nOcc_log) + range_vals/2,
                        min(s_meanT_relPol$nOcc_log) + range_vals/4*3,
                        max(s_meanT_relPol$nOcc_log))


plotting_values <- round(exp(plotting_positions))

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('n Records', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_relPol$nOcc_log
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_relPol$slope_meanT_relPol ~ x_vals,
                   weights = s_meanT_relPol$nOcc_log)

length(s_meanT_relPol$slope_meanT_relPol)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
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



################
## Range size ##
################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10) #use minT to keep same values x
y_lim <- c(-1, 1)

#plot graph
plot(s_maxT_relPol$rangeSize_log10, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_maxT_relPol$rangeSize_log10) -
  min(s_maxT_relPol$rangeSize_log10)

#use minT to keep same values x
plotting_positions <- c(min(s_maxT_relPol$rangeSize_log10),
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_maxT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_maxT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
#mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$latAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$latAmplitude, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                    weights = s_maxT_relPol$nOcc_log)

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$rangeLoc)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$rangeLoc, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_relPol$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_relPol$slope_maxT_relPol ~ x_vals,
                   weights = s_maxT_relPol$nOcc_log)

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



########################
### Elevation median ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$elevMedian)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$elevMedian, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)


#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800




#############################
### Elevational amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$elevAmplitude)
y_lim <- c(-1, 1)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_relPol$elevAmplitude, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#################
### Roundness ###
#################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-1, 1)

#plot graph
plot(s_maxT_relPol$roundness, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800



#################
### Body mass ###
#################

# #select only rows with bodyMass
# s_maxT_relPol_bodyMass <- s_maxT_relPol[
#   complete.cases(s_maxT_relPol$bodyMass),]
# 
# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- c(0, 1000)
# y_lim <- c(-0.2, 0.2)
# 
# #plot graph
# plot(s_maxT_relPol_bodyMass$bodyMass, s_maxT_relPol_bodyMass$slope_maxT_relPol,
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
# #mtext('Body mass', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_maxT_relPol_bodyMass$bodyMass
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_maxT <- lm(s_maxT_relPol_bodyMass$slope_maxT_relPol ~ x_vals,
#                     weights = s_maxT_relPol_bodyMass$nOcc_log)
# 
# length(s_maxT_relPol_bodyMass$slope_maxT_relPol)
# summary(lin_mod_maxT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(0, 100000, length.out = 100)  # Ensuring it starts at 0
# y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
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

#save 800



###################
#### n Records ####
###################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_relPol$nOcc_log)
y_lim <- c(-1, 1)

#plot graph
plot(s_maxT_relPol$nOcc_log, s_maxT_relPol$slope_maxT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_maxT_relPol$slope_maxT_relPol)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
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
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_minT_distEdge$rangeSize_log10, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)




#############################
### Latitudinal amplitude ###
#############################




#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$latAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)

#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$rangeLoc, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
                   weights = s_minT_distEdge$nOcc_log)

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800




########################
### Elevation median ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$elevMedian, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800




#############################
### Elevational amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, max(s_minT_distEdge$elevAmplitude))
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_distEdge$elevAmplitude, s_minT_distEdge$slope_minT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_minT_distEdge$slope_minT_distEdge)
summary(lin_mod_minT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800


#################
### Roundness ###
#################


# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- c(0, 1)
# y_lim <- c(-0.02, 0.02)
# 
# #plot graph
# plot(s_minT_distEdge$roundness, s_minT_distEdge$slope_minT_distEdge,
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
# #mtext('Roundness', side = 1, line = 3.8, cex = 2.5)  
# #mtext('Slope', side = 2, line = 6.5, cex = 2.5)
# 
# #define x range explicitly
# x_vals <- s_minT_distEdge$roundness
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
#                    weights = s_minT_distEdge$nOcc_log)
# 
# length(s_minT_distEdge$slope_minT_distEdge)
# summary(lin_mod_minT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
# y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
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



###################
#### n Records ####
###################


# 
# #restricted ylim weighted by nOcc
# 
# # Define x and y limits
# x_lim <- range(s_minT_distEdge$nOcc_log)
# y_lim <- c(-0.001, 0.001)
# 
# #plot graph
# plot(s_minT_distEdge$nOcc_log, s_minT_distEdge$slope_minT_distEdge,
#      axes = F, xaxs = "i", yaxs = "i",
#      xlab = "", ylab = "", cex = 1.5,
#      pch = 19, col = '#90809020',
#      ylim = y_lim, xlim = x_lim)
# 
# 
# # Define original values and their log-transformed positions
# range_vals <- max(s_minT_distEdge$nOcc_log) -
#   min(s_minT_distEdge$nOcc_log)
# 
# 
# plotting_positions <- c(min(s_minT_distEdge$nOcc_log),
#                         min(s_minT_distEdge$nOcc_log) + range_vals/4,
#                         min(s_minT_distEdge$nOcc_log) + range_vals/2,
#                         min(s_minT_distEdge$nOcc_log) + range_vals/4*3,
#                         max(s_minT_distEdge$nOcc_log))
# 
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
# x_vals <- s_minT_distEdge$nOcc_log
# x_range <- range(x_vals)
# 
# #fit linear model
# lin_mod_minT <- lm(s_minT_distEdge$slope_minT_distEdge ~ x_vals,
#                    weights = s_minT_distEdge$nOcc_log)
# 
# length(s_minT_distEdge$slope_minT_distEdge)
# summary(lin_mod_minT)
# 
# #predict y-values only for positive x-values
# x_seq <- seq(min(plotting_positions) + range_vals/100,
#              max(plotting_positions) - range_vals/100,
#              length.out = 100)  # Ensuring it starts at 0
# y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
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
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$rangeSize_log10, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$latAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
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
### Elevation median ###
########################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_meanT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$elevMedian, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800




#############################
### Elevational amplitude ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, max(s_meanT_distEdge$elevAmplitude))
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_meanT_distEdge$elevAmplitude, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800




#################
### Roundness ###
#################




#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_meanT_distEdge$roundness, s_meanT_distEdge$slope_meanT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_meanT_distEdge$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_meanT <- lm(s_meanT_distEdge$slope_meanT_distEdge ~ x_vals,
                   weights = s_meanT_distEdge$nOcc_log)

length(s_meanT_distEdge$slope_meanT_distEdge)
summary(lin_mod_meanT)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_meanT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
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
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$rangeSize_log10, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
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
mtext('Range size', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                    weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800



#############################
### Latitudinal amplitude ###
#############################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$latAmplitude)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_distEdge$latAmplitude, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                    weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800



#############################
### Latitudinal position  ###
#############################



#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$rangeLoc)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0

#plot graph
plot(s_maxT_distEdge$rangeLoc, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables
#mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$rangeLoc
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800




########################
### Elevation median ###
########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$elevMedian)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_maxT_distEdge$elevMedian, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Median elevation', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#save 800



###########################
### Elevation amplitude ###
###########################

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_maxT_distEdge$elevAmplitude)
y_lim <- c(-0.02, 0.02)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

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
mtext('Elevation amplitude', side = 1, line = 5.5, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)



#################
### Roundness ###
#################


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- c(0, 1)
y_lim <- c(-0.02, 0.02)

#plot graph
plot(s_maxT_distEdge$roundness, s_maxT_distEdge$slope_maxT_distEdge,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#90909020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables
#mtext('Roundness', side = 1, line = 3.8, cex = 2.5)
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_maxT_distEdge$roundness
x_range <- range(x_vals)

#fit linear model
lin_mod_maxT <- lm(s_maxT_distEdge$slope_maxT_distEdge ~ x_vals,
                   weights = s_maxT_distEdge$nOcc_log)

length(s_maxT_distEdge$slope_maxT_distEdge)
summary(lin_mod_maxT)

#predict y-values only for positive x-values
x_seq <- seq(0.01, 0.99, length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_maxT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#909090', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#90909030',
        border = NA)


#save 800


#### random plots to make labels

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
