#load packages
library(rnaturalearth); library(sf); library(units)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Figures/Fig 1/Things for the figure'

#set seed
set.seed(69)



################# IMPORTANT NOTE ###################

#.           make graph fist
#.           use created values for the figure

###################################################


##### MAP #####

#load world map
world <- ne_countries(returnclass = "sf")

#create polygon representing species range
poly = st_polygon(
  list(cbind(c(10,08,09,12,15,16,14,12,14,20,25,32,
               34,35,31,29,29,26,23,20,15,12,10),
             c(30,31,33,33,35,35,37,38,43,46,49,50,
               46,43,34,32,30,29,28,25,27,28,30))))

poly_sf <-  st_sfc(poly, crs = st_crs(world))
poly_sf <- st_as_sf(poly_sf)

#create points representing species occurrences
points_sf <- st_sample(poly_sf, size = 100, type = "random")
points_sf <- st_as_sf(points_sf)

#calculate rel polewardness for all points
points_sf$relPol <- (st_coordinates(points_sf)[,2] - 25) / 25
  
#create data on var contribution for all points
points_sf$varContPol <- (points_sf$relPol +
  rnorm(nrow(points_sf), mean = 0.2, sd = 0.1)) - 0.3

summary(points_sf$varContPol)

#normalize variable contribution from 0 to 5 (for size)
points_sf$varContNormal <- 
  round(((points_sf$varContPol +
             abs(min(points_sf$varContPol))) * 2) + 0.3, 2)

summary(points_sf$varContNormal)
  
#make boxes for the map

#get min and max coord values of the range
ext <- st_bbox(poly_sf)  

#calculate the size of the small box to plot around the range
side_box <-  max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
  min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 4

#calculate central lon and central lat of the range
central_x <- (ext[3] + ext[1]) / 2
names(central_x) <- 'xmean'
central_y <- (ext[4] + ext[2]) / 2
names(central_y) <- 'ymean'

#calculate min and max lon and mean and max lat of the small box
min_x <- central_x - (side_box / 1.65)
names(min_x) <- 'xmin'
max_x <- central_x + (side_box / 1.65)
names(max_x) <- 'xmax'

min_y <- central_y - (side_box / 2)
names(min_y) <- 'ymin'
max_y <- central_y + (side_box / 2)
names(max_y) <- 'ymax'
                          
#create big box
box_df <- data.frame(x = c(min_x,max_x,max_x,min_x,min_x), 
                     y = c(max_y,max_y,min_y,min_y,max_y))

box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(box_df, coords = c('x', 'y'), crs = st_crs(poly_sf))))

#calculate the size of the box to plot in white in the background and give
#space to other stuff
side_big_box <-  max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
  min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 2

#calculate min and max lon and mean and max lat of the big box
min_x_BB <- central_x - (side_big_box / 2)
names(min_x_BB) <- 'xmin'
max_x_BB <- central_x + (side_big_box / 2)
names(max_x_BB) <- 'xmax'

min_y_BB <- central_y - (side_big_box / 1.5)
names(min_y_BB) <- 'ymin'
max_y_BB <- central_y + (side_big_box / 2)
names(max_y_BB) <- 'ymax'

#create big box
big_box_df <- data.frame(x = c(min_x_BB,max_x_BB,max_x_BB,min_x_BB,min_x_BB), 
                     y = c(max_y_BB,max_y_BB,min_y_BB,min_y_BB,max_y_BB))

big_box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(big_box_df, coords = c('x', 'y'), crs = st_crs(poly_sf))))


##### Calculate distance from the edges to show the 'centralness'

#cut the range out of the box
range_cut <- st_difference(box, poly_sf)

#calculate the dist from each point to edge in km
points_sf$distEdge <- as.numeric(
  set_units(st_distance(points_sf, range_cut), km)[,1])

#calculate index (0-1) of distance to edge
points_sf$distEdgeNormal_0 <- points_sf$distEdge / max(points_sf$distEdge)

#this option for the points uses the same calculation of the relPol but inversed
points_sf$distEdgeNormal <- 1 - points_sf$relPol

#create vector for sizes in points
points_sf$centralCex <- ((7 - log(points_sf$distEdge)) / 1.3) + 
  rnorm(nrow(points_sf), mean = 0.25, sd = 0.4)

summary(points_sf$centralCex)

#create exemplary values of var contribution showing an increase towards edges
points_sf$varContrCentral <- (points_sf$centralCex - 1.2) / 5

summary(points_sf$varContrCentral)










###### PLOT MAP SHOWING THE CENTRAL GRADIENT ######

#Legend: "Larger ranges showing populations close to the range centre being less affected by climatic extremes than those closer to range edges."

#Reference in text: "The relationship between climatic variables and species occurrences changes between the range centre and range edges (Fig. 1A)"


#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot small box to inform the coordinates
#plot(box, add = T, col = '#ffffbf20')

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = points_sf$centralCex ^ 1.4,
     col = '#2c7bb680')

#save 1000





###### PLOT GRAPH SHOWING THE CENTRAL GRADIENT ######

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

# Plot the data points again
plot(data$distEdgeNormal, data$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
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



###### PLOT GRAPH SHOWING THE CENTRAL GRADIENT (ALTERNATIVE MEASUREMENT) ######

# Set seed for reproducibility
set.seed(123)

# Generate 33, 34, and 33 random values for distEdgeNormal
distEdgeNormal_A <- runif(33, min = 0, max = 150) #1st segment
distEdgeNormal_B <- runif(34, min = 150, max = 350) #2nd segment
distEdgeNormal_C <- runif(33, min = 350, max = 500) #2nd segment


# Generate ydata with a linear relationship plus some noise
ydata_A <-
  ((-0.002 * distEdgeNormal_A + rnorm(33, mean = 0.2, sd = 0.1)) + 0.13) 

ydata_B <-
  ((-0.002 * distEdgeNormal_B + rnorm(34, mean = 0, sd = 1)) + 0.4) / 5

ydata_C <-
  ((-0.002 * distEdgeNormal_C + rnorm(33, mean = -0.2, sd = 0.1)) + 0.9) 

# Combine into a data frame
data_A <- data.frame(distEdgeNormal_A, ydata_A)
data_B <- data.frame(distEdgeNormal_B, ydata_B)
data_C <- data.frame(distEdgeNormal_C, ydata_C)

# Fit the linear model
lin_mod_A <- lm(ydata_A ~ distEdgeNormal_A, data = data_A)
lin_mod_B <- lm(ydata_B ~ distEdgeNormal_B, data = data_B)
lin_mod_C <- lm(ydata_C ~ distEdgeNormal_C, data = data_C)

# Display the model summary
summary(lin_mod_A)
summary(lin_mod_B)
summary(lin_mod_C)

#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

# Create a sequence of x-values within the range of the data
x_vals_A <- seq(min(data_A$distEdgeNormal),
                max(data_A$distEdgeNormal),
                length.out = 33)
x_vals_B <- seq(min(data_B$distEdgeNormal),
                max(data_B$distEdgeNormal),
                length.out = 34)
x_vals_C <- seq(min(data_C$distEdgeNormal),
                max(data_C$distEdgeNormal),
                length.out = 33)

# Create a data frame for prediction
new_data_A <- data.frame(distEdgeNormal_A = x_vals_A)
new_data_B <- data.frame(distEdgeNormal_B = x_vals_B)
new_data_C <- data.frame(distEdgeNormal_C = x_vals_C)

# Predict y-values based on the model
predicted_A <- predict(lin_mod_A, newdata = new_data_A)
predicted_B<- predict(lin_mod_B, newdata = new_data_B)
predicted_C <- predict(lin_mod_C, newdata = new_data_C)

#set y and x lims
ylim <- c(-0.5, 0.5)
xlim <- c(0, 530)

# Plot the data points again
plot(data_A$distEdgeNormal, data_A$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Distance from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

points(data_B$distEdgeNormal, data_B$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650')

points(data_C$distEdgeNormal, data_C$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

# Add the restricted regression line
lines(new_data_A$distEdgeNormal, predicted_A, col = '#2c7bb6', lwd = 7)
lines(new_data_B$distEdgeNormal, predicted_B, col = '#2c7bb6', lwd = 7)
lines(new_data_C$distEdgeNormal, predicted_C, col = '#2c7bb6', lwd = 7)

# For graphical purposes, I will replot joining the lines

# Plot the data points again
plot(data_A$distEdgeNormal, data_A$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Distance from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

points(data_B$distEdgeNormal, data_B$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

points(data_C$distEdgeNormal, data_C$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650'
       )

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

lines(new_data_A$distEdgeNormal, predicted_A, col = '#2c7bb6', lwd = 7)
lines(c(new_data_A$distEdgeNormal[33], new_data_C$distEdgeNormal[1]),
      c(predicted_A[33],predicted_C[1]),
      col = '#2c7bb6', lwd = 7)
lines(new_data_C$distEdgeNormal, predicted_C, col = '#2c7bb6', lwd = 7)









###### PLOT MAP SHOWING THE LATITUDINAL GRADIENT ######

#Legend: large ranges showing climatic extremes increase in explanatory power towards the poleward edges.

#Reference in text: "Climatic variables have higher explanatory power towards the poleward range limit (Fig. 1B)"

#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = points_sf$varContNormal ^ 1.6,
     col = '#2c7bb680')

#save 1000






###### PLOT GRAPH SHOWING THE LATITUDINAL GRADIENT ######

# Set seed for reproducibility
set.seed(123)

# Generate 100 random values for distEdgeNormal
relPol <- runif(100, min = 0, max = 1)

# Generate ydata with a linear relationship plus some noise
ydata <- ((relPol + rnorm(100, mean = 0, sd = 0.2))) / 2

# Combine into a data frame
data <- data.frame(relPol, ydata)

# Fit the linear model
lin_mod <- lm(ydata ~ relPol, data = data)

# Display the model summary
summary(lin_mod)

#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

# Create a sequence of x-values within the range of the data
x_vals <- seq(min(data$relPol),
              max(data$relPol),
              length.out = 100)

# Create a data frame for prediction
new_data <- data.frame(relPol = x_vals)

# Predict y-values based on the model
predicted <- predict(lin_mod, newdata = new_data)

#set y and x lims
ylim <- c(-0.5, 0.5)
xlim <- c(0, 1.05)

# Plot the data points again
plot(data$relPol, data$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Relative polewardness',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

# Add the restricted regression line
lines(x_vals, predicted, col = '#2c7bb6', lwd = 7)

#save 800






### PLOT GRAPH SHOWING THE LATITUDINAL GRADIENT (ALTERNATIVE MEASUREMENT) ###

# Set seed for reproducibility
set.seed(123)

# Generate 33, 34, and 33 random values for distEdgeNormal
relPol_A <- runif(33, min = 0, max = 0.3) #1st segment
relPol_B <- runif(34, min = 0.3, max = 0.7) #2nd segment
relPol_C <- runif(33, min = 0.7, max = 1) #2nd segment

# Generate ydata with a linear relationship plus some noise
ydata_A <-
  ((relPol_A + rnorm(33, mean = -0.2, sd = 0.1)) - 0.13) 

ydata_B <-
  ((relPol_B + rnorm(34, mean = 0, sd = 0.32)) / 4) - 0.05

ydata_C <-
  ((relPol_C + rnorm(33, mean = 0.2, sd = 0.1)) - 0.8) 


# Combine into a data frame
data_A <- data.frame(relPol_A, ydata_A)
data_B <- data.frame(relPol_B, ydata_B)
data_C <- data.frame(relPol_C, ydata_C)

# Fit the linear model
lin_mod_A <- lm(ydata_A ~ relPol_A, data = data_A)
lin_mod_B <- lm(ydata_B ~ relPol_B, data = data_B)
lin_mod_C <- lm(ydata_C ~ relPol_C, data = data_C)

# Display the model summary
summary(lin_mod_A)
summary(lin_mod_B)
summary(lin_mod_C)

#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

# Create a sequence of x-values within the range of the data
x_vals_A <- seq(min(data_A$relPo),
                max(data_A$relPo),
                length.out = 33)
x_vals_B <- seq(min(data_B$relPo),
                max(data_B$relPo),
                length.out = 34)
x_vals_C <- seq(min(data_C$relPo),
                max(data_C$relPo),
                length.out = 33)

# Create a data frame for prediction
new_data_A <- data.frame(relPol_A = x_vals_A)
new_data_B <- data.frame(relPol_B = x_vals_B)
new_data_C <- data.frame(relPol_C = x_vals_C)

# Predict y-values based on the model
predicted_A <- predict(lin_mod_A, newdata = new_data_A)
predicted_B<- predict(lin_mod_B, newdata = new_data_B)
predicted_C <- predict(lin_mod_C, newdata = new_data_C)

#set y and x lims
ylim <- c(-0.5, 0.5)
xlim <- c(0, 1.05)

# Plot the data points again
plot(data_A$relPol, data_A$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Relative polewardness',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

points(data_B$relPol, data_B$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

points(data_C$relPol, data_C$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

# Add the restricted regression line
lines(new_data_A$relPol, predicted_A, col = '#2c7bb6', lwd = 7)
lines(new_data_B$relPol, predicted_B, col = '#2c7bb6', lwd = 7)
lines(new_data_C$relPol, predicted_C, col = '#2c7bb6', lwd = 7)


# For graphical purposes, I will replot joining the lines

# Plot the data points again
plot(data_A$relPol, data_A$ydata, 
     pch = 19, cex = 0.8, col = '#2c7bb650',
     axes = F, , xaxs = "i", yaxs = "i",
     ylab = 'Variable contribution',
     xlab = 'Relative polewardness',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = ylim,
     xlim = xlim)

points(data_B$relPol, data_B$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

points(data_C$relPol, data_C$ydata, 
       pch = 19, cex = 0.8, col = '#2c7bb650')

#add axes
axis(1, pos = -0.4)
axis(2, pos = 0, las=2)

lines(new_data_A$relPol_A, predicted_A, col = '#2c7bb6', lwd = 7)
lines(c(new_data_A$relPol[33], new_data_C$relPol[1]),
      c(predicted_A[33],predicted_C[1]),
      col = '#2c7bb6', lwd = 7)
lines(new_data_C$relPol, predicted_C, col = '#2c7bb6', lwd = 7)


#save 800






###### PLOT MAP SHOWING EXPECTED DIFFERENCE BETWEEN EXTREMES AND MEANS ######

#Legend: large ranges showing climatic extremes increase in explanatory power towards the edges, and means don't show a strong pattern.

#Reference in text: "Climatic variables have higher explanatory power towards the poleward range limit (Fig. 1C)"

#divide points to show extremes and means
shuffle_pts <- sample(c(1:nrow(points_sf)), size = nrow(points_sf))
extre_pts <- points_sf[shuffle_pts[c(1:50)],]
mean_pts <- points_sf[shuffle_pts[c(51:100)],]

#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(extre_pts),
     add = T, pch = 19, cex = extre_pts$centralCex ^ 1.4,
     col = '#2c7bb680')

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(mean_pts),
     add = T, pch = 19, cex = sqrt(mean_pts$centralCex),
     col = '#f46d4380')


#save plot (width = 1000)




###### PLOT GRAPH SHOWING THE CENTRE-EDGE GRADIENT ######

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
lines(x_vals, predicted / 6, col = '#f46d43', lwd = 7)






#fit linear model
lin_mod_minT <- lm(points_sf$varContPol ~ points_sf$distEdgeNormal)
lin_mod_minT_2 <- lm((points_sf$varContPol/6) + 0.2 ~ points_sf$distEdgeNormal)

abline(lin_mod_minT, col = '#9930FF', lwd = 7)
abline(lin_mod_minT_2, col = '#f46d43', lwd = 7)


#save plot (width = 1000)








