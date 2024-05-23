####  This script plots Figure 3 (???) exemplifying the SHAP results per species

#load libraries
library(sf); library(rnaturalearth); library(raster)
library(rworldmap)  #change to rnaturalearth

#list wds
wd_res_shap <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Comparison'
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_elevation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/Elevation_Tozer'

#select species
species <- 'Cynomys ludovicianus'

#manually select the files for the species I need
setwd(wd_res_species)
sps_res <- lapply(list.files(pattern = species), read.csv)
names(sps_res) <- gsub(paste0(species, '_'), '',
                       gsub("_.csv", "", list.files(pattern = species)))

#get a world map from package rworldmap
world <- getMap(resolution = 'low')
world <- st_as_sf(world)
original_crs <- crs(world)
world <-  st_transform(world, crs = 3857) #to cope with the edges issue

#species range map
setwd(wd_ranges)
range <- st_read(dsn = wd_ranges, layer = species)

#select the range representing only the native range of the species
range <- range[range$legend == 'Extant (resident)',]

#unify all features
range <- st_union(range)

#get min and max coord values of the range
ext <- st_bbox(range) 

#get the values for the margin around the box
x_mar <- abs(ext[3] - ext[1]) / 3
y_mar <- abs(ext[4] - ext[2]) / 3

#make a dataframe with coordinates to create a margin around range extent
box_df <- data.frame(x = ext[c(1,3,3,1,1)] + x_mar * c(-1,1,1,-1,-1), 
                     y = ext[c(4,4,2,2,4)] + y_mar * c(1,1,-1,-1,1))

#make the box
box1 <- st_as_sfc(
  st_bbox(
    st_as_sf(box_df, coords = c('x', 'y'), crs = crs(range))))

#project the box to the same crs as the world object
box <- st_transform(box1, crs = 3857)

#crop world map to show the region where the range
region <- st_crop(world, box)

#reconvert the cropped region to WGS84
region <- st_transform(region, crs = crs(range))

#create a spatial object with the occurrences
sps_res_sp <- lapply(sps_res, function(x) {
  st_as_sf(x, crs = crs(range),coords = c('decimalLongitude', 'decimalLatitude'))})

#select only presences in the results I want to show
sps_sp_pr <- lapply(sps_res_sp, function(x) {
  x[which(x$Occurrence == 1),]})


#### PLOT
par(mfrow = c(1,1), mar = c(1,1,2,1))

#identify the limits of the contributions
lim_cont <- range(c(sps_sp_pr$minT$Min_T_SHAP,
                    sps_sp_pr$meanT$Mean_T_SHAP,
                    sps_sp_pr$maxT$Max_T_SHAP))

#make a colourRamp
colramp <- colorRampPalette(c("#FF0000", "#FFFFFF", "#0000FF"))

#populate the tables with the colours to be plotted 
sps_sp_pr$minT$colours <- colramp(100)[cut(c(-0.5,0.5,sps_sp_pr$minT$Min_T_SHAP), 
                                           breaks = 100)][-c(1,2)]

sps_sp_pr$meanT$colours <- colramp(100)[cut(c(-0.5,0.5,sps_sp_pr$meanT$Mean_T_SHAP), 
                                            breaks = 100)][-c(1,2)]

sps_sp_pr$maxT$colours <- colramp(100)[cut(c(-0.5,0.5,sps_sp_pr$maxT$Max_T_SHAP), 
                                           breaks = 100)][-c(1,2)]

#plot MIN T
plot(st_geometry(region), col = 'khaki', main = 'Min T', cex.main = 2)
plot(st_geometry(range), add = T, col = '#238b4590')
plot(st_geometry(sps_sp_pr$minT), add = T, pch = 21, 
     col = 'black', bg = sps_sp_pr$minT$colours, cex = 1.5)



# same plot but with elevation


#load elevation map
setwd(wd_elevation)
elevation <- raster('SRTM15Plus_world.tif')

#convert range from 'sfc_POLYGON' to sf object
shp_mask <- as_Spatial(range)

#crop the raster to the same extent as the mask object
ele_crop <- crop(elevation, shp_mask)

#mask elevation by species range
ele_sps <- mask(ele_crop, shp_mask)

#plot MIN T elevation
plot(st_geometry(region), col = 'khaki', main = 'Min T', cex.main = 2)
plot(ele_sps, add = T)
plot(st_geometry(sps_sp_pr$minT), add = T, pch = 21, 
     col = 'black', bg = sps_sp_pr$minT$colours, cex = 1)




# #plot MEAN T
# plot(st_geometry(region), col = 'khaki', main = 'Mean T', cex.main = 2)
# plot(st_geometry(range), add = T, col = '#238b4590')
# plot(st_geometry(sps_sp_pr$meanT), add = T, pch = 21, 
#      col = 'black', bg = sps_sp_pr$meanT$colours, cex = 1.5)
# 
# #plot MAX T
# plot(st_geometry(region), col = 'khaki', main = 'Max T', cex.main = 2)
# plot(st_geometry(range), add = T, col = '#238b4590')
# plot(st_geometry(sps_sp_pr$maxT), add = T, pch = 21, 
#      col = 'black', bg = sps_sp_pr$maxT$colours, cex = 1.5)

# par(mfrow = c(1,1), mar = c(2,2,2,2))

#plot legend
# plot(box, border = NA)
# 
# #plot map legend
# myGradientLegend(valRange = c(-0.5, 0.5), 
#                  pos=c(0.28,0,0.73,.017),
#                  color = colramp(100), 
#                  side = 1,
#                  n.seg = 1,
#                  cex = 2)
