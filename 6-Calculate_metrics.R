####  This script calculates measurements to classify the points and interpret
####  the results:
####
####  1 - Range size for each species *numeric area km^2
####
####  2 - Range shape for each species (how round is the range?) *gradient
####
####  3 - Centralness of each point (how close to the edges) *gradient
####
####  4 - Absolute polewardness of each point (dist from equator) *numeric %
####
####  5 - Relative polewardness of each point (dist from warm edge) *gradient
####
####  6 - Distance from range edges of each point *numeric km
####
####  7 - Elevation of each point *numeric m
####
####  8 - Body mass for each species *numeric g
####
####  9 - Biome for each point

#load libraries
library(sf); library(units); library(raster)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_pts_measure <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Point_and_range_measurements'
wd_elevation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_mass <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Mammal_trait_data/Smith_etal_2003_Ecology'
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'

#load elevation raster
setwd(wd_elevation)
elev <- raster('wc2.1_2.5m_elev.tif')

#load biomes shapefile
biomes <- st_read(dsn = wd_biomes, layer = 'wwf_terr_ecos')

#load mammal body mass table
setwd(wd_mass)
mass <- read.csv('Mammals_bodymass_Smmith_2003.csv')

#list species
setwd(wd_thinned_occ)
sps_list <- gsub('_thinned.csv', '', list.files())

######## Get measurements for all species ######

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #loas species occurrences
  setwd(wd_thinned_occ)
  sps_occ <- read.csv(paste0(sps, '_thinned.csv'))
  
  #select only presence occurrences
  pr_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'PRESENT'),]
  
  #select only columns we are interested in 
  
  pr_sps2 <- pr_sps[,c('species','key','decimalLongitude','decimalLatitude',
                       'datasetKey')]
  
  #check if there are absence data
  abs_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'ABSENT'),]
  if(nrow(abs_sps) != 0){
    warning(paste0('THERE IS ABSENT DATA FOR ', sps_list[i]))
  }
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #select the range representing only the native range of the species
  range <- range[range$legend == 'Extant (resident)',]
  
  #unify all features
  sps_range2 <- st_union(range)
  
  ### Calculate range area (km^2)
  pr_sps2$rangeSize <- set_units(st_area(sps_range2), km^2)
  
  ### Calculate roundness (Roundness = 4π * Area / Perimeter^2)
  area_range <- set_units(st_area(sps_range2), km^2)
  perimetre <- set_units(st_length(st_cast(sps_range2, 'MULTILINESTRING')), km)

  pr_sps2$roundness <- as.numeric((4*pi*area_range) / (perimetre^2))
  
  ### Calculate the centralness of each point
  
  #generate stratified points
  reg_pts <- st_sample(sps_range2, size = 1000, type = "regular")
  reg_pts <- st_as_sf(reg_pts, crs = crs(sps_range2))
  
  #make a box around the range
  ext <- st_bbox(sps_range2)  #get min and max coord values of the range
  
  x_mar <- abs(ext[3] - ext[1]) / 3 #get the values for the margin around the box
  y_mar <- abs(ext[4] - ext[2]) / 3
  
  box_df <- data.frame(x = ext[c(1,3,3,1,1)] + x_mar * c(-1,1,1,-1,-1), #dataframe
                       y = ext[c(4,4,2,2,4)] + y_mar * c(1,1,-1,-1,1))
  
  box <- st_as_sfc(          #make the box    
           st_bbox(
            st_as_sf(box_df, coords = c('x', 'y'), crs = crs(sps_range2))))
  
  #cut the range out of the box
  range_cut <- st_difference(box, sps_range2)
  
  #calculate the largest possible distance from the edges (in m)
  max_dist <- max(st_distance(reg_pts, range_cut))
  
  #make sf objects for the points
  pr_sps2_sf <- st_as_sf(pr_sps2, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = crs(sps_range2))
  
  #calculate the centralness (dist point to edge / max dist to edge)
  pr_sps2$centralness <- as.numeric(st_distance(pr_sps2_sf, range_cut) / max_dist)
  
  ### Calculate the dist from each point to edge in km
  pr_sps2$distEdge <- set_units(st_distance(pr_sps2_sf, range_cut), km)[,1]
  
  ### Calculate the absolute polarwardness of each point
  pr_sps2$absPolarwardness <- abs(pr_sps2$decimalLatitude) / 90
  
  ### Calculate the relative polarwardness of each point
  
  #get latitudinal range
  ymax <- st_bbox(sps_range2)$ymax
  ymin <- st_bbox(sps_range2)$ymin
  
  #flag if range is crossed by the equator (deal with this part of the code after)
  if(ymax > 0 & ymin < 0 | ymax < 0 & ymin > 0){
    pr_sps2$NOTE <- 'Crossed by Equator'
  }
  
  #identify warm and cold edges
  if(abs(ymax) > abs(ymin)){
    cold_edge <- ymax
    warm_edge <- ymin
  }else{
    cold_edge <- ymin
    warm_edge <- ymax
  }
  
  #calculate distance from cold to warm edge
  lat_range <- abs(cold_edge) - abs(warm_edge)
  
  #calculate relative distance from warm edge
  pr_sps2$relPolarwardness <- 
           1-((abs(cold_edge) - abs(pr_sps2$decimalLatitude)) / lat_range)
  
  ### Extract the elevation of each point
  pr_sps2$elevation <- extract(elev, pr_sps2_sf)
  
  ### Extract the biome info of each point
  
  #change projections to cope with edge problem
  pr_sps3_sf <-  st_transform(pr_sps2_sf, crs = 3857) 
  biomes2 <- st_transform(biomes, crs = 3857) 
  
  #extract the biome info
  biomes_sps <- st_intersection(biomes2, pr_sps3_sf)
  pr_sps2$biome <- biomes_sps$BIOME
         
  ### Get average body mass for the species
  if(sps_list[i] %in% mass$Species){
    if(length(unique(
                  mass$BodyMass[which(mass$Species == sps_list[i])])) == 1){
      pr_sps2$bodyMass <- 
        unique(mass$BodyMass[which(mass$Species == sps_list[i])])
    }else{
      pr_sps2$bodyMass <- mass$BodyMass[which(mass$Species == sps_list[i])]
    }
  }else{
    pr_sps2$bodyMass <- NA
  }
  
  t <- mass[which(mass$Species == sps_list[i]),]
  
  #save table
  setwd(wd_pts_measure)
  write.csv(pr_sps2, paste0(sps,'_point_range_metrics.csv'), row.names = F)
}


##################################
plot(sps_range2)
plot(pr_sps2_sf, add = T, pch = 19, col = 'magenta', cex = 0.4)
plot(pr_sps2_sf[which.max(pr_sps2$centralness),], add = T, pch = 19, col = 'green')
plot(pr_sps2_sf[which.min(pr_sps2$centralness),], add = T, pch = 19, col = 'red')
plot(pr_sps2_sf[which.max(pr_sps2$relPolarwardness),], add = T, pch = 19, col = 'blue')
plot(pr_sps2_sf[which.min(pr_sps2$relPolarwardness),], add = T, pch = 19, col = 'orange')

###### NOTES TO BE INCORPORATED INTO THE ANNOTATION AFTER TAKEN CARE OFF

# load th points used for shap, but don't use pa (they will not be used for the results interpretation, it would only take computational time).