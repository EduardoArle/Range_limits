#load packages
library(rnaturalearth); library(sf); library(units); library(mgcv)

#set seed for reproducibility
set.seed(69)



################################
##### Make polygon for map #####
################################



#load world map (just to have a reference system)
world <- ne_countries(returnclass = "sf")

#create polygon representing species range
poly = st_polygon(
  list(cbind(c(10,08,09,12,15,16,14,12,14,20,25,32,
               34,35,31,29,29,26,23,20,15,12,10),
             c(30,31,33,33,35,35,37,38,43,46,49,50,
               46,43,34,32,30,29,28,25,27,28,30))))

poly_sf <-  st_sfc(poly, crs = st_crs(world))
poly_sf <- st_as_sf(poly_sf)

#get min and max coord values of the range
ext <- st_bbox(poly_sf)

#create points representing species occurrences
points_sf <- st_sample(poly_sf, size = 100, type = "random")
points_sf <- st_as_sf(points_sf)

#calculate rel polewardness for all points
points_sf$relPol <- (st_coordinates(points_sf)[,2] - ext$ymin) / ext$ymin

#create box to calculate distance from edge

#calculate the size of the box to plot around the range
lon_mar <- abs(ext[3] - ext[1]) / 8
lat_mar <- abs(ext[4] - ext[2]) / 8
                          
#calculate coords for the box
b_xmin <- ext$xmin - lon_mar
b_xmax <- ext$xmax + lon_mar
b_ymin <- ext$ymin - lat_mar
b_ymax <- ext$ymax + lat_mar

#create box to calculate distance from edge
box_df <- data.frame(x = c(b_xmin,b_xmax,b_xmax,b_xmin,b_xmin),
                     y = c(b_ymax,b_ymax,b_ymin,b_ymin,b_ymax))

#make the box
box <- st_as_sfc(          
  st_bbox(st_as_sf(box_df, coords = c('x', 'y'), crs = st_crs(poly_sf))))

#cut the range out of the box
range_cut <- st_difference(box, poly_sf)

#visualise steps
plot(box)
plot(poly_sf, add = T, col = '#f23598')
plot(range_cut, add = T, col = '#89f336')
plot(st_geometry(points_sf), add = T, pch = 21, col = '#000000', bg = '#98fbcb')



#calculate the dist from each point to edge in km
points_sf$distEdge <- as.numeric(
  set_units(st_distance(points_sf, range_cut), km)[,1])

distEdge <- points_sf$distEdge


#calculate the relative polewardness of each point (0-1)

#get latitudinal range
ymax <- st_bbox(poly_sf)$ymax
ymin <- st_bbox(poly_sf)$ymin

#calculate latitudinal range
lat_range <- ymax - ymin

points_sf$relPol <- (st_coordinates(points_sf)[,2] - ymin) / lat_range

relPol <- points_sf$relPol



###################################
##### Plot latitudinal graphs #####
###################################


#prepare parametres for axes

#set number of ticks per axis
n_ticks_x <- 5
n_ticks_y <- 3

#set size of stuff
cex_axes <- 2
cex_lables <- 2.5


#generate ydata with a linear relationship plus some noise
ydata <- ((relPol + rnorm(100, mean = 0, sd = 0.1)) - 0.45)

#combine into a data frame
data <- data.frame(relPol, ydata)


###### MAP ######. FIG 1-A-i

#Legend: large ranges showing climatic extremes increase in explanatory power towards the poleward edges.

#Reference in text: "Climatic variables have higher explanatory power towards the poleward range limit (Fig. 1B)"


#set parametres for plotting
par(mar = c(2,2,2,2), pty="s", mfrow = c(1,1))

#plot small box to inform the coordinates
plot(box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#transform ydata values into something related to size
cex_pts <- (ydata + abs(min(ydata)) + 1) ^ 3

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = cex_pts,
     col = '#2c7bb680')

#save 1000



########## LINE ###########. FIG 1-A-ii


#fit the linear model
lin_mod <- lm(ydata ~ relPol, data = data)

#create a sequence of x-values within the range of the data
x_vals <- seq(min(data$relPol),
              max(data$relPol),
              length.out = 100)

#create a data frame for prediction (include 0 so that line starts from axis)
new_data <- data.frame(relPol = c(0, x_vals)) 

#predict y-values based on the model (sum min val + a bit for plotting's sake)
predicted <- predict(lin_mod, newdata = new_data) + 0.49

#make the same combina for the data$ydata
data$ydata <- data$ydata + 0.49

#set graph parametres
par(mar = c(6,9,5,5), pty="m", mfrow = c(1,1))

#determine xlim for this plot (keep ylim same as previous plot)
xlim <- c(0,1)
ylim <- c(0,1)

#plot the data points
plot(data$relPol, data$ydata, 
     pch = 19, cex = 1, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = '', xlab = '',
     ylim = ylim,
     xlim = xlim)

#set number of ticks per axis
n_ticks_x <- 3
n_ticks_y <- 2

#calculate distance between ticks
dist_tick_x <- (xlim[2] - xlim[1]) / (n_ticks_x - 1)
dist_tick_y <- (ylim[2] - ylim[1]) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes, labels = c(0, NA))


#add arrows that Shahar suggested
text(par("usr")[1] - 0.07, 0.94, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 0.0697, 0.94, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.91, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.88, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.85, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.82, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.79, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.76, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.73, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.70, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.67, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.64, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.61, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.58, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.55, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.52, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.49, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.46, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.43, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.40, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.37, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.34, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.31, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.28, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.22, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.19, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.16, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.13, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.10, expression("|"), xpd = TRUE, cex = 2)  

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("Variable contribution", side = 2, line = 6.5, cex = 2.5)

#add the restricted regression line
lines(c(0, x_vals), (predicted / 1.3) + 0.12, col = '#2c7bb6', lwd = 8)

#save 800



########## CURVE ###########. FIG 1-A-iii


# Generate 33, 34, and 33 random values for distEdgeNormal
relPol_A <- runif(33, min = 0, max = 0.3) #1st segment
relPol_B <- runif(34, min = 0.3, max = 0.7) #2nd segment
relPol_C <- runif(33, min = 0.7, max = 1) #2nd segment


# Generate ydata with a non-linear relationship plus some noise
ydata_A <-
  ((1.8 * relPol_A + rnorm(33, mean = -0.2, sd = 0.1)) - 0.25) 

ydata_B <-
  ((relPol_B + rnorm(34, mean = 0, sd = 0.3)) / 6) 

ydata_C <-
  ((1.5 *relPol_C + rnorm(33, mean = 0.2, sd = 0.1)) - 1.2) 


#combine into a data frame
data <- data.frame(relPol = c(relPol_A, relPol_B, relPol_C),
                   ydata = c(ydata_A, ydata_B, ydata_C))

#fit a gam model
gam_1 <- gam(ydata ~ s(relPol, k = 5),
             data = data,
             method="REML")

#set y and x lims for this plot 
xlim <- c(0,1)
ylim <- c(floor(min(data$ydata) * 10) / 10, - floor(min(data$ydata) * 10) / 10)

#set parametres for plotting
par(mar = c(6,9,5,5), pty="m", mfrow = c(1,1))

#plot 
plot.gam(gam_1, select = 1, residuals = F, shade = F,
         col = '#2c7bb6', se = F, lwd = 8, rug = F,
         axes = F, , xaxs = "i", yaxs = "i",
         ylab = '', xlab = '',
         ylim = ylim,
         xlim = xlim) 


# #add points
# points(data$relPol, data$ydata, # +  rnorm(100, mean = 0.1, sd = 0.1), 
#        pch = 19, cex = 1, col = '#2c7bb650')


#calculate distance between ticks
dist_tick_x <- (xlim[2] - xlim[1]) / (n_ticks_x - 1)
dist_tick_y <- (ylim[2] - ylim[1]) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes, labels = c(NA, 0, NA))


#add arrows that Shahar suggested
text(par("usr")[1] - 18.2, 0.46, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 0.0697, 0.46, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.45, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.44, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.43, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.42, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.41, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.40, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.39, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.38, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.37, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.36, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.35, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.34, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.33, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.32, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.31, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.30, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.29, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.28, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.27, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.26, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.25, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.24, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.23, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.22, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.21, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.20, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.19, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.18, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.17, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.16, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.15, expression("|"), xpd = TRUE, cex = 2)
text(par("usr")[1] - 0.0697, 0.14, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.13, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.12, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.11, expression("|"), xpd = TRUE, cex = 2)
text(par("usr")[1] - 0.0697, 0.10, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.09, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.08, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.07, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.06, expression("|"), xpd = TRUE, cex = 2)  

#add arrows that Shahar suggested
text(par("usr")[1] - 0.07, 0.5, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 0.0697, 0.5, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.45, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.4, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.35, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.3, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.2, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.15, expression("|"), xpd = TRUE, cex = 2)  


text(par("usr")[1] - 0.07, -0.5, expression("\u2193"), xpd = TRUE, cex = 3)  # Down arrow
text(par("usr")[1] - 0.0697, -0.5, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.45, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.4, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.35, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.3, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.2, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.15, expression("|"), xpd = TRUE, cex = 2)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  




# # Example data
# y <- rnorm(20)
# 
# # Create the plot without the default y-axis
# plot(y, type = "b", yaxt = "n", ylim = c(-3, 3), xlab = "Index", ylab = "")
# 
# # Custom y-axis: center '0', no labels on top and bottom
# axis(2, at = c(-3, 0, 3), labels = c("", "0", ""), las = 1)






# # Add y-axis label with arrows
# mtext(expression("Performance"~"\u2193"~"0"~"\u2191"), side = 2, line = 3)







#save 800



###################################
##### Plot dist edges graphs ######
###################################


########## LINE ###########. FIG 1-B-ii


#prepare parametres for axes

#set number of ticks per axis
n_ticks_x <- 3
n_ticks_y <- 2

#set size of stuff
cex_axes <- 2
cex_lables <- 2.5

#generate ydata with a linear relationship plus some noise
ydata <- ((relPol + rnorm(100, mean = 0, sd = 0.1)))

#combine into a data frame
data <- data.frame(relPol, ydata)

#generate ydata with a linear relationship plus some noise
ydata <- ((-0.004 * distEdge + rnorm(100, mean = 0, sd = 0.2)) + 0.6) / 1.5

#combine into a data frame
data <- data.frame(distEdge, ydata)

#fit the linear model
lin_mod <- lm(ydata ~ distEdge, data = data)

#create a sequence of x-values within the range of the data
x_vals <- seq(min(data$distEdge),
              max(data$distEdge),
              length.out = 100)

x_vals <- x_vals[x_vals < 250]

#create a data frame for prediction
new_data <- data.frame(distEdge = x_vals)

#predict y-values based on the model (sum min val + a bit for plotting's sake)
predicted <- predict(lin_mod, newdata = new_data) + 0.49

#make the same combina for the data$ydata
data$ydata <- data$ydata + 0.49

#set y and x lims
ylim <- c(0, 1)
xlim <- c(0, ceiling(max(distEdge) / 100) * 100)

#set graph parametres
par(mar = c(6,9,5,5), pty="m", mfrow = c(1,1))

#plot the data points
plot(data$distEdge, data$ydata, 
     pch = 19, cex = 1, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = '', xlab = '',
     ylim = c(0, 1),
     xlim = c(0, 250))


#calculate distance between ticks
dist_tick_x <- (250 - 0) / (n_ticks_x - 1)
dist_tick_y <- (1 - 0) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = 0, at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes, labels = c(0, NA))


#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  
mtext("Variable contribution", side = 2, line = 6.5, cex = 2.5)

#add the restricted regression line
lines(x_vals, predicted, col = '#2c7bb6', lwd = 8)

#add arrows that Shahar suggested
text(par("usr")[1] - 18.2, 0.94, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 17.9, 0.94, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.91, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.88, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.85, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.82, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.79, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.76, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.73, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.70, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.67, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.64, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.61, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.58, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.55, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.52, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.49, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.46, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.43, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.40, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.37, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.34, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.31, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.28, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.22, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.19, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.16, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.13, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.10, expression("|"), xpd = TRUE, cex = 2)  

#save 800



########## CURVE ###########. FIG 1-B-iii

#prepare parametres for axes

#set number of ticks per axis
n_ticks_x <- 3
n_ticks_y <- 2

#set size of stuff
cex_axes <- 2
cex_lables <- 2.5

set.seed(123)
x <- seq(0, 250, length.out = 200)
y <- 0.5 * exp(-x / 100) + rnorm(length(x), 0, 0.02)  # scaled for smaller range >0
y[y < 0] <- 0  # ensure no negatives

df <- data.frame(x = x, y = y)

# Fit GAM with k=3
gam_fit <- gam(y ~ s(x, k = 3), data = df, method = "REML")

x_pred <- seq(0, 250, length.out = 250)

#adjusted prediction
pred <- predict(gam_fit, newdata = data.frame(x = x_pred)) * 2

# Plot points with no axes and custom limits
plot(df$x, df$y,
     pch = 19, cex = 1, col = '#2c7bb600',
     axes = FALSE, xaxs = "i", yaxs = "i",
     xlab = '', ylab = '',
     xlim = c(0, 250), ylim = c(0, 1))

# Add GAM line on top, matching color and thicker width
lines(x_pred, pred, col = '#2c7bb6', lwd = 7)



#set y and x lims for this plot 
xlim <- c(0, 250)
ylim <- c(0, max(data$ydata))



#calculate distance between ticks
dist_tick_x <- (xlim[2] - xlim[1]) / (n_ticks_x - 1)
dist_tick_y <- (ylim[2] - ylim[1]) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 1, at = ticks_y, cex.axis = cex_axes, labels = c(0, NA))

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  

#add arrows that Shahar suggested
text(par("usr")[1] - 18.2, 0.94, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 17.9, 0.94, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.91, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.88, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.85, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.82, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.79, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.76, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.73, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.70, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.67, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.64, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.61, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.58, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.55, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.52, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.49, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.46, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.43, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.40, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.37, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.34, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.31, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.28, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 17.9, 0.22, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.19, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.16, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.13, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 17.9, 0.10, expression("|"), xpd = TRUE, cex = 2)  

#save 800




###### PLOT GRAPH SHOWING EXPECTED DIFFERENCE BETWEEN EXTREMES AND MEANS ######


#set parametres for plotting
par(mar = c(2,2,2,2), pty="s", mfrow = c(1,1))

#plot small box to inform the coordinates
plot(box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#transform ydata values into something related to size
cex_pts <- (ydata + abs(min(ydata)) + 1) ^ 3

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = cex_pts,
     col = '#2c7bb680')

#save 1000



########## LINE ###########. FIG 1-A-ii


#fit the linear model
lin_mod <- lm(ydata ~ relPol, data = data)

#create a sequence of x-values within the range of the data
x_vals <- seq(min(data$relPol),
              max(data$relPol),
              length.out = 100)

#create a data frame for prediction (include 0 so that line starts from axis)
new_data <- data.frame(relPol = c(0, x_vals)) 

#predict y-values based on the model
predicted <- predict(lin_mod, newdata = new_data)

#set graph parametres
par(mar = c(6,9,5,5), pty="m", mfrow = c(1,1))

#determine xlim for this plot (keep ylim same as previous plot)
xlim <- c(0,1)
ylim <- c(floor(min(ydata) * 10) / 10, - floor(min(ydata) * 10) / 10)

#plot the data points
plot(data$relPol, data$ydata, 
     pch = 19, cex = 1, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = '', xlab = '',
     ylim = ylim,
     xlim = xlim)


#calculate distance between ticks
dist_tick_x <- (xlim[2] - xlim[1]) / (n_ticks_x - 1)
dist_tick_y <- (ylim[2] - ylim[1]) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes, labels = c(NA, 0, NA))

#add arrows that Shahar suggested
text(par("usr")[1] - 0.07, 0.583, expression("\u2191"), xpd = TRUE, cex = 3)  # Up arrow
text(par("usr")[1] - 0.0697, 0.525, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.475, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.45, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.4, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.35, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.3, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, 0.2, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, 0.174, expression("|"), xpd = TRUE, cex = 2) 


text(par("usr")[1] - 0.07, -0.583, expression("\u2193"), xpd = TRUE, cex = 3)  # Down arrow
text(par("usr")[1] - 0.0697, -0.525, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.475, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.45, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.4, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.35, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.3, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.25, expression("|"), xpd = TRUE, cex = 2)  
text(par("usr")[1] - 0.0697, -0.2, expression("|"), xpd = TRUE, cex = 2) 
text(par("usr")[1] - 0.0697, -0.174, expression("|"), xpd = TRUE, cex = 2) 

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("Variable contribution", side = 2, line = 6.5, cex = 2.5)

#add the restricted regression line
lines(c(0, x_vals), predicted, col = '#2c7bb6', lwd = 8)
lines(c(0, x_vals), predicted / 3, col = '#ffc00c', lwd = 8)




#save 800




#generate ydata with a linear relationship plus some noise
ydata <- ((relPol + rnorm(100, mean = 0, sd = 0.1)) - 0.45)

#combine into a data frame
data <- data.frame(relPol, ydata)

#fit the linear model
lin_mod <- lm(ydata ~ relPol, data = data)

#create a sequence of x-values within the range of the data
x_vals <- seq(min(data$relPol),
              max(data$relPol),
              length.out = 100)

#create a data frame for prediction (include 0 so that line starts from axis)
new_data <- data.frame(relPol = c(0, x_vals)) 

#predict y-values based on the model
predicted <- predict(lin_mod, newdata = new_data)

#set graph parametres
par(mar = c(6,9,5,5), pty="m", mfrow = c(1,1))

#determine xlim for this plot 
xlim <- c(0,1)

#plot the data points
plot(data$relPol, data$ydata, 
     pch = 19, cex = 1, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = '', xlab = '',
     ylim = ylim,
     xlim = xlim)


#calculate distance between ticks
dist_tick_x <- (xlim[2] - xlim[1]) / (n_ticks_x - 1)
dist_tick_y <- (ylim[2] - ylim[1]) / (n_ticks_y - 1)

#calculate positions of ticks
ticks_x <- seq(xlim[1], xlim[2], dist_tick_x)
ticks_y <- seq(ylim[1], ylim[2], dist_tick_y)

#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes, labels = c(NA, 0, NA))

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  
mtext("Variable contribution", side = 2, line = 6.5, cex = 2.5)

#add the restricted regression line
lines(c(0, x_vals), predicted, col = '#2c7bb6', lwd = 8)


#save 800




###### PLOT GRAPH SHOWING EXPECTED DIFFERENCE BETWEEN EXTREMES AND MEANS ######


#plot 
plot.gam(gam_1, select = 1, residuals = F, shade = F,
         col = '#2c7bb6', se = F, lwd = 8, rug = F,
         axes = F, , xaxs = "i", yaxs = "i",
         ylab = '', xlab = '',
         ylim = ylim,
         xlim = xlim) 


#add axes
axis(1, pos = ylim[1], at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes)

#add labels
mtext("Relative polewardness", side = 1, line = 3.8, cex = 2.5)  

#save 800











###### PLOT GRAPH SHOWING THE CENTRAL GRADIENT (ALTERNATIVE MEASUREMENT) ######


#generate 33, 34, and 33 random values for distEdgeNormal (for relPol)
distEdge_A <- runif(33, min = 0, max = 175) #1st segment
distEdge_B <- runif(34, min = 175, max = 575) #2nd segment
distEdge_C <- runif(33, min = 575, max = 750) #2nd segment


#generate 33, 34, and 33 random values for distEdgeNormal (for distEdge)
distEdge_1 <- runif(33, min = 0, max = 83) #1st segment
distEdge_2 <- runif(34, min = 84, max = 166) #2nd segment
distEdge_3 <- runif(33, min = 167, max = 250) #2nd segment



#generate ydata with a linear relationship plus some noise (for distEdge)
ydata_1 <-
  ((-0.002 * distEdge_1 + rnorm(33, mean = 0.2, sd = 0.05)) + 0.13) 

ydata_B <-
  ((-0.002 * distEdge_2 + rnorm(34, mean = 0, sd = 0.3)) + 0.4) / 5

ydata_C <-
  ((-0.002 * distEdge_3 + rnorm(33, mean = -0.2, sd = 0.05)) + 1.2) 

#combine into a data frame
# data_A <- data.frame(distEdge_A, ydata_A)
# data_B <- data.frame(distEdge_B, ydata_B)
# data_C <- data.frame(distEdge_C, ydata_C)

#fit the linear model
# lin_mod_A <- lm(ydata_A ~ distEdge_A, data = data_A)
# lin_mod_B <- lm(ydata_B ~ distEdge_B, data = data_B)
# lin_mod_C <- lm(ydata_C ~ distEdge_C, data = data_C)

#combine into a data frame
data <- data.frame(distEdge = c(distEdge_A, distEdge_B, distEdge_C),
                   ydata = c(ydata_A, ydata_B, ydata_C))

data_2 <- data.frame(distEdge = c(distEdge_1, distEdge_2, distEdge_3),
                     ydata = c(ydata_1, ydata_2, ydata_3))

#fit a gam model
gam_1 <- gam(ydata ~ s(distEdge, k = 5),
                      data = data,
                      method="REML")

#plot 
plot.gam(gam_1, select = 1, residuals = F, shade = F,
         col = '#2c7bb6', se = F, lwd = 8, rug = F,
         axes = F, , xaxs = "i", yaxs = "i",
         ylab = '', xlab = '',
         ylim = ylim,
         xlim = xlim)

#add points
points(data$distEdge, data$ydata +  rnorm(100, mean = 0.1, sd = 0.1), 
       pch = 19, cex = 1, col = '#2c7bb650')

#add axes
axis(1, pos = 0, at = ticks_x, cex.axis = cex_axes)
axis(2, las = 2, at = ticks_y, cex.axis = cex_axes)

#add labels
mtext("Distance to range edge", side = 1, line = 3.8, cex = 2.5)  
#mtext("Variable contribution", side = 2, line = 6.5, cex = 2.5)

#add the restricted regression line
#lines(x_vals, predicted, col = '#2c7bb6', lwd = 7)

#save 800


###### PLOT MAP SHOWING THE CENTRAL GRADIENT ######

#Legend: "Larger ranges showing populations close to the range centre being less affected by climatic extremes than those closer to range edges."

#Reference in text: "The relationship between climatic variables and species occurrences changes between the range centre and range edges (Fig. 1A)"


#set parametres for plotting
par(mar = c(2,2,2,2), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
# plot(big_box, border = NA)

#plot small box to inform the coordinates
plot(box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#transform ydata values into something related to size
cex_pts <- (ydata + abs(min(ydata)) + 1) ^ 3

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = cex_pts,
     col = '#2c7bb680')

#save 1000











# Generate ydata with a non-linear relationship plus some noise
ydata_A <-
  ((1.2 * relPol_A + rnorm(33, mean = -0.15, sd = 0.1)) - 0.1) 

ydata_B <-
  ((relPol_B + rnorm(34, mean = 0, sd = 0.3)) / 6) 

ydata_C <-
  ((1.2 *relPol_C + rnorm(33, mean = 0.2, sd = 0.1)) - 1.1) 

#combine into a data frame
data <- data.frame(relPol = c(relPol_A, relPol_B, relPol_C),
                   ydata = c(ydata_A, ydata_B, ydata_C))

#fit a gam model
gam_2 <- gam(ydata ~ s(relPol, k = 4),
             data = data,
             method="REML")

#plot 
plot.gam(gam_2, select = 1, residuals = F, shade = F,
         col = '#ffc00c', se = F, lwd = 8, rug = F,
         axes = F, , xaxs = "i", yaxs = "i",
         ylab = '', xlab = '',
         ylim = ylim,
         xlim = xlim) 

#save 800



###### PLOT MAP SHOWING EXPECTED DIFFERENCE BETWEEN EXTREMES AND MEANS ######

#Legend: large ranges showing climatic extremes increase in explanatory power towards the edges, and means don't show a strong pattern.

#Reference in text: "Climatic variables have higher explanatory power towards the poleward range limit (Fig. 1C)"

#divide points to show extremes and means
shuffle_pts <- sample(c(1:nrow(points_sf)), size = nrow(points_sf))
extre_pts <- points_sf[shuffle_pts[c(1:50)],]
mean_pts <- points_sf[shuffle_pts[c(51:100)],]

#set parametres for plotting
par(mar = c(2,2,2,2), pty="s", mfrow = c(1,1))

#plot small box to inform the coordinates
plot(box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(extre_pts),
     add = T, pch = 19, cex = (extre_pts$relPol + 0.5) * 6,
     col = '#2c7bb680')

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(mean_pts),
     add = T, pch = 19, cex = (mean_pts$relPol + 0.3) * 3,
     col = '#ffc00c80')


#save 1000


#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)



#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(mean_pts),
     add = T, pch = 19, cex = sqrt(mean_pts$centralCex),
     col = '#f46d4380')


#save plot (width = 1000)









#######################################################################












###### PLOT GRAPH SHOWING EXPECTED DIFFERENCE BETWEEN EXTREMES AND MEANS ######

# Set seed for reproducibility
set.seed(123)

# Generate 100 random values for distEdgeNormal
distEdgeNormal <- runif(100, min = 0, max = 500)

# Generate ydata with a linear relationship plus some noise
ydata <- ((-0.002 * distEdgeNormal + rnorm(100, mean = 0, sd = 0.2)) + 0.6) / 2

# Combine into a data frame
data <- data.frame(distEdgeNormal, ydata)

# Fit the linear model
lin_mod <- lm(ydata ~ distEdgeNormal, data = data)

# Display the model summary
summary(lin_mod)

#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

# Create a sequence of x-values within the range of the data
x_vals <- seq(min(data$distEdgeNormal),
              max(data$distEdgeNormal),
              length.out = 100)

# Create a data frame for prediction
new_data <- data.frame(distEdgeNormal = x_vals)

# Predict y-values based on the model
predicted <- predict(lin_mod, newdata = new_data)

#set y and x lims
ylim <- c(-0.5, 0.5)
xlim <- c(0, 530)

# Plot the data points (with no points)
plot(data$distEdgeNormal, data$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Distance from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

# Add the restricted regression line
lines(x_vals, predicted, col = '#2c7bb6', lwd = 7)
lines(x_vals, (predicted / 1.1) - 0.07, col = '#f46d43', lwd = 7)
lines(x_vals, predicted / 6, col = '#ffc00c', lwd = 7)


# Plot the data points (with no points)
plot(data$distEdgeNormal, data$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb600',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Distance from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

# Add the restricted regression line
lines(x_vals, predicted, col = '#2c7bb6', lwd = 7)
lines(x_vals, predicted * -0.9, col = '#f46d43', lwd = 7)
lines(x_vals, predicted / 6, col = '#ffc00c', lwd = 7)








#fit linear model
lin_mod_minT <- lm(points_sf$varContPol ~ points_sf$distEdgeNormal)
lin_mod_minT_2 <- lm((points_sf$varContPol/6) + 0.2 ~ points_sf$distEdgeNormal)

abline(lin_mod_minT, col = '#9930FF', lwd = 7)
abline(lin_mod_minT_2, col = '#f46d43', lwd = 7)


#save plot (width = 1000)








