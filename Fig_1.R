#load packages
library(rnaturalearth); library(sf); library(units)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Figures/Fig 1/Things for the figure'

#set seed
set.seed(69)

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
points_sf <- st_sample(poly_sf, size = 50, type = "random")
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
plot(box, add = T, col = '#ffffbf20')

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#ffffbf80', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = points_sf$centralCex ^ 1.4,
     col = '#2c7bb680')

#label axes
text(-5, 37.4, 'Latitude', srt = 90, cex = 1.2)
text(21.7, 16, 'Longitude', srt = 00, cex = 1.2)

#add ticks to axes
points(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)
points((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)
points(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)

points(st_bbox(box)[1] - 0.45,
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.7)
points(st_bbox(box)[1] - 0.45,
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
       pch = '—',
       cex = 0.7)
points(st_bbox(box)[1] - 0.45,
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.7)

#add values to axes
text(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)
text((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       (st_bbox(box)[1] + st_bbox(box)[3]) / 2, 1)),
     cex = 1.1)
text(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)

text(st_bbox(box)[1] - 2.75,
     st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - 2.75,
     (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
     labels = paste0(round(
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - 2.75,
     st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)

#save plot (width = 800)

###### PLOT GRAPH SHOWING THE CENTRE-EDGE GRADIENT ######

#set y and x lims
ylim <- c(-0.2, 0.8)
xlim <- c(0, 1)

#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

#minT (make points invisible)
plot(points_sf$distEdgeNormal, points_sf$varContPol, 
     pch = 19, cex = 1, col = '#FF808000',
     ylab = 'Variable contribution',
     xlab = 'Distance from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = c(ylim[1], ylim[2]))

#plot rectangle for background colour
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = '#ffffbf20') 

#fit linear model
lin_mod_minT <- lm(points_sf$varContPol ~ points_sf$distEdgeNormal)
abline(lin_mod_minT, col = '#2c7bb6', lwd = 7)

#save plot (width = 800)


#create table object to save
points_table <- cbind(st_coordinates(points_sf), st_drop_geometry(points_sf))
names(points_table)[c(1,2)] <- c('lon', 'lat')

########################

#save new coordinates for the map
setwd(wd_tables)
write.csv(points_table, 'Points_large_map.csv', row.names = F)






###### PLOT MAP SHOWING THE LATITUDINAL GRADIENT ######

#Legend: large ranges showing climatic extremes increase in explanatory power towards the poleward edges.

#Reference in text: "Climatic variables have higher explanatory power towards the poleward range limit (Fig. 1B)"

#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot small box to inform the coordinates
plot(box, add = T)

#plot the polygon
plot(poly_sf, lwd = 3, border = '#707070', col = '#F0F0F0', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_sf),
     add = T, pch = 19, cex = points_sf$varContNormal,
     col = '#0000FF80')

#label axes
text(-5, 37.4, 'Latitude', srt = 90, cex = 1.2)
text(21.7, 16, 'Longitude', srt = 00, cex = 1.2)

#add ticks to axes
points(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)
points((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)
points(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - 0.44,
       pch = '|',
       cex = 0.7)

points(st_bbox(box)[1] - 0.45,
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.7)
points(st_bbox(box)[1] - 0.45,
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
       pch = '—',
       cex = 0.7)
points(st_bbox(box)[1] - 0.45,
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.7)

#add values to axes
text(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)
text((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       (st_bbox(box)[1] + st_bbox(box)[3]) / 2, 1)),
     cex = 1.1)
text(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - 2.1,
     labels = paste0(round(
       st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)

text(st_bbox(box)[1] - 2.75,
     st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - 2.75,
     (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
     labels = paste0(round(
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - 2.75,
     st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)

#save plot (width = 800)

###### PLOT GRAPH SHOWING THE LATITUDINAL GRADIENT ######

#set y and x lims
ylim <- c(-0.2, 0.9)
xlim <- c(0, 1)


#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

#minT
plot(points_sf$relPol, points_sf$varContPol, 
     pch = 19, cex = 1, col = '#0000FF50',
     ylab = 'Variable contribution',
     xlab = 'Relative polewardness',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = c(ylim[1], ylim[2]))

#fit linear model
lin_mod_minT <- lm(points_sf$varContPol ~ points_sf$relPol)
abline(lin_mod_minT, col = '#0000FF', lwd = 7)

#save plot (width = 800)








#######################################################
######################## LOAD #########################
#######################################################

#load the coordinations when runinf from here
setwd(wd_tables)
points_table <- read.csv('Coordinates_points_large_map.csv')

#create a points_sf object with the modified longitutes
points_sf <- st_as_sf(points_table, coords = c('lon', 'lat'), crs = st_crs(world))


#create polygon representing species small range
poly_S = st_polygon(
  list(cbind(c(51, 50.7, 51, 51.2, 52.4, 53, 54.2, 55, 54.7, 53.5,
               52.2, 51),
             c(16, 16.8, 18, 18.1, 18.3, 19, 18.8, 18, 17.2, 16.2,
               15.7, 16))))

poly_S_sf <-  st_sfc(poly_S, crs = st_crs(world))
poly_S_sf <- st_as_sf(poly_S_sf)

#create points representing species occurrences
points_S_sf <- st_sample(poly_S_sf, size = 25, type = "random")
points_S_sf <- st_as_sf(points_S_sf)

#create data on var contribution for all points
points_S_sf$varCont <- rnorm(nrow(points_S_sf), mean = 0.1, sd = 0.1)

summary(points_S_sf$varCont)

#get min and max coord values of the range
ext <- st_bbox(poly_S_sf) 
minLat <- ext[2]
latAmp <- ext[4] - ext[2]

#calculate rel polewardness for all points
points_S_sf$relPol <- (st_coordinates(points_S_sf)[,2] - minLat) / latAmp

summary(points_S_sf$relPol)

#normalize variable contribution from 0 to 5 (for size)
points_S_sf$varContNormal <- 
  round(((points_S_sf$varCont +
            abs(min(points_S_sf$varCont))) * 2) + 0.7, 2)

summary(points_S_sf$varContNormal)

#make boxes for the map

#calculate the size of the small box to plot around the range
#make it big, to show that the range is small
side_box <-  (max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
  min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 4) * 3

#calculate central lon and central lat of the range
central_x <- (ext[3] + ext[1]) / 2
names(central_x) <- 'xmean'
central_y <- (ext[4] + ext[2]) / 2
names(central_y) <- 'ymean'

#calculate min and max lon and mean and max lat of the small box
min_x <- central_x - (side_box / 2.1)
names(min_x) <- 'xmin'
max_x <- central_x + (side_box / 2.1)
names(max_x) <- 'xmax'

min_y <- central_y - (side_box / 2.2)
names(min_y) <- 'ymin'
max_y <- central_y + (side_box / 2.2)
names(max_y) <- 'ymax'

#create big box
box_df <- data.frame(x = c(min_x,max_x,max_x,min_x,min_x), 
                     y = c(max_y,max_y,min_y,min_y,max_y))

box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(box_df, coords = c('x', 'y'), crs = st_crs(poly_S_sf))))

#calculate the size of the box to plot in white in the background and give
#space to other stuff
side_big_box <-  (max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
  min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 2) * 3

#calculate min and max lon and mean and max lat of the big box
min_x_BB <- central_x - (side_big_box / 1.7)
names(min_x_BB) <- 'xmin'
max_x_BB <- central_x + (side_big_box / 2)
names(max_x_BB) <- 'xmax'

min_y_BB <- central_y - (side_big_box / 1.7)
names(min_y_BB) <- 'ymin'
max_y_BB <- central_y + (side_big_box / 2)
names(max_y_BB) <- 'ymax'

#create big box
big_box_df <- data.frame(x = c(min_x_BB,max_x_BB,max_x_BB,min_x_BB,min_x_BB), 
                         y = c(max_y_BB,max_y_BB,min_y_BB,min_y_BB,max_y_BB))

big_box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(big_box_df, coords = c('x', 'y'), crs = st_crs(poly_S_sf))))


###### PLOT MAP SHOWING THE LATITUDINAL GRADIENT ######

#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot small box to inform the coordinates
plot(box, add = T)

#plot the polygon
plot(poly_S_sf, lwd = 3, border = '#707070', col = '#F0F0F0', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_S_sf),
     add = T, pch = 19, cex = points_S_sf$varContNormal,
     col = '#0000FF80')

#label axes
text(min_x_BB + (side_big_box / 19),
     (min_y + max_y) / 2,
     'Latitude', srt = 90, cex = 1.2)
text((min_x + max_x) / 2,
     min_y_BB + (side_big_box / 12),
     'Longitude', srt = 00, cex = 1.2)

#add ticks to axes
points(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)
points((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)
points(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)

points(st_bbox(box)[1] - (side_box / 100),
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.75)
points(st_bbox(box)[1] - (side_box / 100),
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
       pch = '—',
       cex = 0.75)
points(st_bbox(box)[1] - (side_box / 100),
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.75)

#add values to axes
text(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - (side_box / 17),
     labels = paste0(round(
       st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)
text((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
     st_bbox(box)[2] - (side_box / 17),
     labels = paste0(round(
       (st_bbox(box)[1] + st_bbox(box)[3]) / 2, 1)),
     cex = 1.1)
text(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - (side_box / 17),
     labels = paste0(round(
       st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)

text(st_bbox(box)[1] - (side_box / 15),
     st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - (side_box / 15),
     (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
     labels = paste0(round(
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - (side_box / 15),
     st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)

#save plot (width = 800)


###### PLOT GRAPH SHOWING THE LATITUDINAL GRADIENT ######

#set y and x lims
ylim <- c(-0.2, 0.9)
xlim <- c(0, 1)


#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

#minT
plot(points_S_sf$relPol, points_S_sf$varCont, 
     pch = 19, cex = 1, col = '#0000FF50',
     ylab = 'Variable contribution',
     xlab = 'Rel polewardness / Dist from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = c(ylim[1], ylim[2]))

#fit linear model
lin_mod_minT <- lm(points_S_sf$varCont ~ points_S_sf$relPol)
abline(lin_mod_minT, col = '#0000FF', lwd = 7)

#save plot (width = 800)


######## TROPICAL ###############

#create polygon representing species in the tropica
poly_T = st_polygon(
  list(cbind(c(20.1, 19.2, 15.3, 15, 18.4, 28, 32.2, 34, 35, 34.5,
               35.2, 40, 43, 45, 40, 28, 22, 20.1),
             c(-16, -15.8, -14, -8.2, 3.4, 19, 21.8, 21, 18, 16.2,
               14.8, 9, 2, -5.5, -12, -17, -17.3, -16))))

poly_T_sf <-  st_sfc(poly_T, crs = st_crs(world))
poly_T_sf <- st_as_sf(poly_T_sf)

plot(poly_T_sf)

#create points representing species occurrences
points_T_sf <- st_sample(poly_T_sf, size = 45, type = "random")
points_T_sf <- st_as_sf(points_T_sf)

#create data on var contribution for all points
points_T_sf$varCont <- rnorm(nrow(points_T_sf), mean = 0.2, sd = 0.2)

summary(points_T_sf$varCont)

#get min and max coord values of the range
ext <- st_bbox(poly_T_sf) 
minLat <- ext[2]
latAmp <- ext[4] - ext[2]

#calculate rel polewardness for all points
points_T_sf$relPol <- (st_coordinates(points_T_sf)[,2] - minLat) / latAmp

summary(points_T_sf$relPol)

#normalize variable contribution from 0 to 5 (for size)
points_T_sf$varContNormal <- 
  round(((points_T_sf$varCont +
            abs(min(points_T_sf$varCont))) * 2) + 0.7, 2)

summary(points_T_sf$varContNormal)

#make boxes for the map

#calculate the size of the small box to plot around the range
side_box <-  max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
                min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 4

#calculate central lon and central lat of the range
central_x <- (ext[3] + ext[1]) / 2
names(central_x) <- 'xmean'
central_y <- (ext[4] + ext[2]) / 2
names(central_y) <- 'ymean'

#calculate min and max lon and mean and max lat of the small box
min_x <- central_x - (side_box / 2)
names(min_x) <- 'xmin'
max_x <- central_x + (side_box / 2)
names(max_x) <- 'xmax'

min_y <- central_y - (side_box / 2)
names(min_y) <- 'ymin'
max_y <- central_y + (side_box / 2)
names(max_y) <- 'ymax'

#create big box
box_df <- data.frame(x = c(min_x,max_x,max_x,min_x,min_x), 
                     y = c(max_y,max_y,min_y,min_y,max_y))

box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(box_df, coords = c('x', 'y'), crs = st_crs(poly_T_sf))))

#calculate the size of the box to plot in white in the background and give
#space to other stuff
side_big_box <-  max(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) +
                    min(c(abs(ext[3] - ext[1]), abs(ext[4] - ext[2]))) / 2 

#calculate min and max lon and mean and max lat of the big box
min_x_BB <- central_x - (side_big_box / 1.675)
names(min_x_BB) <- 'xmin'
max_x_BB <- central_x + (side_big_box / 1.675)
names(max_x_BB) <- 'xmax'

min_y_BB <- central_y - (side_big_box / 1.675)
names(min_y_BB) <- 'ymin'
max_y_BB <- central_y + (side_big_box / 1.675)
names(max_y_BB) <- 'ymax'

#create big box
big_box_df <- data.frame(x = c(min_x_BB,max_x_BB,max_x_BB,min_x_BB,min_x_BB), 
                         y = c(max_y_BB,max_y_BB,min_y_BB,min_y_BB,max_y_BB))

big_box <- st_as_sfc(          #make the box    
  st_bbox(st_as_sf(big_box_df, coords = c('x', 'y'), crs = st_crs(poly_T_sf))))


###### PLOT MAP SHOWING THE LATITUDINAL GRADIENT ######

#set parametres for plotting
par(mar = c(0,0,0,0), pty="s", mfrow = c(1,1))

#plot big white box to make room for the things I need to add around
plot(big_box, border = NA)

#plot small box to inform the coordinates
plot(box, add = T)

#plot the polygon
plot(poly_T_sf, lwd = 3, border = '#707070', col = '#F0F0F0', add = T)

#plot the occurrence points (size is proportional to SHAP value)
plot(st_geometry(points_T_sf),
     add = T, pch = 19, cex = points_T_sf$varContNormal,
     col = '#0000FF80')

#label axes
text(min_x_BB + (side_big_box / 26),
     (min_y + max_y) / 2,
     'Latitude', srt = 90, cex = 1.2)
text((min_x + max_x) / 2,
     min_y_BB + (side_big_box / 22),
     'Longitude', srt = 00, cex = 1.2)

#add ticks to axes
points(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)
points((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)
points(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
       st_bbox(box)[2] - (side_box / 100),
       pch = '|',
       cex = 0.75)

points(st_bbox(box)[1] - (side_box / 100),
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.75)
points(st_bbox(box)[1] - (side_box / 100),
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
       pch = '—',
       cex = 0.75)
points(st_bbox(box)[1] - (side_box / 100),
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
       pch = '—',
       cex = 0.75)

#add values to axes
text(st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - (side_box / 16),
     labels = paste0(round(
       st_bbox(box)[1] + (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)
text((st_bbox(box)[1] + st_bbox(box)[3]) / 2,
     st_bbox(box)[2] - (side_box / 16),
     labels = paste0(round(
       (st_bbox(box)[1] + st_bbox(box)[3]) / 2, 1)),
     cex = 1.1)
text(st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10,
     st_bbox(box)[2] - (side_box / 16),
     labels = paste0(round(
       st_bbox(box)[3] - (st_bbox(box)[3] - st_bbox(box)[1]) / 10, 1)),
     cex = 1.1)

text(st_bbox(box)[1] - (side_box / 15),
     st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[2] + (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - (side_box / 15),
     (st_bbox(box)[2] + st_bbox(box)[4]) / 2,
     labels = paste0(round(
       (st_bbox(box)[2] + st_bbox(box)[4]) / 2, 1)),
     cex = 1.1, srt = 90)
text(st_bbox(box)[1] - (side_box / 15),
     st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10,
     labels = paste0(round(
       st_bbox(box)[4] - (st_bbox(box)[4] - st_bbox(box)[2]) / 10, 1)),
     cex = 1.1, srt = 90)

#add dashed line showing equator
# abline(a = 0, b = 0, lty = 'dotted', lwd = 2, col = 'red')

#save plot (width = 800)


###### PLOT GRAPH SHOWING THE LATITUDINAL GRADIENT ######

#set y and x lims
ylim <- c(-0.2, 0.9)
xlim <- c(0, 1)


#set parametres for plotting
par(mar = c(5,5,5,5), pty="s", mfrow = c(1,1))

#minT
plot(points_T_sf$relPol, points_T_sf$varCont, 
     pch = 19, cex = 1, col = '#0000FF50',
     ylab = 'Variable contribution',
     xlab = 'Rel polewardness / Dist from edge',
     cex.lab = 1.2,
     cex.axis = 1.2,
     ylim = c(ylim[1], ylim[2]))

#fit linear model
lin_mod_minT <- lm(points_T_sf$varCont ~ points_T_sf$relPol)
abline(lin_mod_minT, col = '#0000FF', lwd = 7)

#save plot (width = 800)


