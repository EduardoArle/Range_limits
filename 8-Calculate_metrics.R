####  This script calculates measurements to classify the points and interpret
####  the results:
####
####  1 - Absolute polewardness of each point (dist from equator) *numeric abs lat
####  
####  2 - Relative polewardness of each point (dist from warm edge) *gradient
####
####  3 - Distance from the edge of each point *numeric
####
####  4 - Elevation of each point *numeric m
####
####  5 - Biome for each point
####
####  6 - Range size for each species *numeric area km^2
####
####  7 - Range location for each species *numeric abs mean lat
####
####  8 - Latitudinal amplitude of each species range *numeric
####
####  9 - Range shape for each species (how round is the range?) *gradient
####
####  10 - Taxonomic order (HAS TO BE PUT BEFORE SCRIP 5)
####
####  11 - Body mass
####
####  12 - Number of presences

#load libraries
library(sf); library(units); library(raster)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Species_occ'
wd_elevation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_mass <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Mammal_trait_data/Smith_etal_2003_Ecology'
wd_tax_harm <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review'
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'
wd_pts_metrics <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables/Occurrence_metrics'

#set to avoid edge issues
sf_use_s2(FALSE)

#load elevation raster
setwd(wd_elevation)
elev <- raster('wc2.1_2.5m_elev.tif')

#load biomes shapefile
biomes <- st_read(dsn = wd_biomes, layer = 'wwf_terr_ecos')

#load taxonomic harmonisation GBIF backbone
setwd(wd_tax_harm)
tax_harm <- read.csv('Harmonised_table.csv')

#load harmonised mammal body mass table
setwd(wd_mass)
mass <- read.csv('Harmonised_Mammals_bodymass_Smmith_2003.csv')

#list species
setwd(wd_occ)
sps_list <- gsub('_occ.csv', '', list.files())

######## Get measurements for all species ######

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #load species occurrences
  setwd(wd_occ)
  sps_occ <- read.csv(paste0(sps, '_occ.csv'))
  
  ## Include number of presences
  sps_occ$nOcc <- sum(sps_occ$Occurrence)
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #work around to fix edge issue
  range <- st_make_valid(range)
  range <- st_simplify(range, dTolerance = 0.01, preserveTopology = TRUE)
  
  #unify all features
  sps_range2 <- st_union(range)
  
  #split into fragments
  frags <- st_cast(sps_range2, "POLYGON")
  
  #local equal-area CRS (metres)
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  crs_laea <- paste0("+proj=laea +lat_0=", cen[2], " +lon_0=", cen[1],
                     " +datum=WGS84 +units=m +no_defs")
  
  #areas (km2)
  frags_m <- st_transform(frags, crs_laea)
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #keep largest fragment
  frag_dom <- which.max(a_km2)
  range_dom <- frags[frag_dom, ]
  
  ## Include range area (km^2)
  sps_occ$rangeSize <- set_units(st_area(range_dom), km^2)
  
  ## Calculate roundness (Roundness = 4π * Area / Perimeter^2)
  area_range <- set_units(st_area(range_dom), km^2)
  perimetre <- set_units(st_length(st_cast(range_dom, 'MULTILINESTRING')), km)

  sps_occ$roundness <- as.numeric((4*pi*area_range) / (perimetre^2))
  
  ## Calculate the dist from each point to edge in km
  
  #make a box around the range, get min and max coord values of the range
  ext <- st_bbox(range_dom)  
  
  #get the values for the margin around the box
  x_mar <- abs(ext[3] - ext[1]) / 3 
  y_mar <- abs(ext[4] - ext[2]) / 3
  
  #dataframe
  box_df <- data.frame(x = ext[c(1,3,3,1,1)] + x_mar * c(-1,1,1,-1,-1), 
                       y = ext[c(4,4,2,2,4)] + y_mar * c(1,1,-1,-1,1))
  
  #make the box 
  box <- st_as_sfc(st_bbox(st_as_sf(box_df, coords = c('x', 'y'),
                                    crs = crs(range_dom))))
  
  #cut the range out of the box
  range_cut <- st_difference(box, range_dom)
  
  #make sf objects for the points
  sps_occ_sf <- st_as_sf(sps_occ, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = crs(sps_range2))
  
  #calculate the dist
  sps_occ$distEdge <- set_units(st_distance(sps_occ_sf, range_cut), km)[,1]
  
  ## Calculate the absolute polewardness of each point
  sps_occ$absPolewardness <- abs(sps_occ$decimalLatitude)
  
  ## Calculate the relative polewardness of each point
  
  #get latitudinal range
  ymax <- st_bbox(range_dom)$ymax
  ymin <- st_bbox(range_dom)$ymin
  
  #flag if range is crossed by the equator (deal with this part of the code after)
  if(ymax > 0 & ymin < 0 | ymax < 0 & ymin > 0){
    sps_occ$NOTE <- 'Crossed by Equator'
  }else{
    sps_occ$NOTE <- NA
  }
  
  #identify warm and cold edges
  if(abs(ymax) > abs(ymin)){
    cold_edge <- ymax
    warm_edge <- ymin
  }else{
    cold_edge <- ymin
    warm_edge <- ymax
  }
  
  #calculate dist from cold to warm edge (special calculation for cross equator)
  if(ymax > 0 & ymin < 0 | ymax < 0 & ymin > 0){
    lat_range <- max(c(abs(cold_edge), abs(warm_edge)))
  }else{
    lat_range <- abs(cold_edge - warm_edge)
  }
  
  
  #calculate relative dist from warm edge, consider the cross equator sps
  sps_occ$relPolewardness <- 
      1-((abs(cold_edge) - abs(sps_occ$decimalLatitude)) / lat_range)

  ## Extract the elevation of each point
  sps_occ$elevation <- extract(elev, sps_occ_sf)
  
  ## Extract the biome info of each point
  
  #change projections to cope with edge problem
  sps_occ2_sf <-  st_transform(sps_occ_sf, crs = 3857) 
  biomes2 <- st_transform(biomes, crs = 3857) 
  
  #extract the biome info
  biomes_sps <- st_join(sps_occ2_sf, biomes2["BIOME"], left = TRUE)
  sps_occ$biome <- biomes_sps$BIOME
         
  ##Get species order and average body mass
  
  #harmonise IUCN name to GBIF name
  gbif_name <- tax_harm$gbif_name[match(sps_list[i], tax_harm$iucn_name)]
  
  #if no GBIF name, set NA
  if(is.na(gbif_name) || !nzchar(gbif_name)){
    sps_occ$bodyMass <- NA
    sps_occ$order <- NA

  } else {
    
    #get body mass values from Smith table
    bm <- mass$BodyMass[which(mass$gbif_name == gbif_name)]
    
    #remove NA values
    bm <- bm[!is.na(bm)]
    
    #assign average body mass
    if(length(bm) == 0){
      sps_occ$bodyMass <- NA
    } else {
      sps_occ$bodyMass <- mean(bm)
    }
    
    #get order from Smith table
    ord <- tax_harm$order[match(sps_list[i], tax_harm$iucn_name)]
    ord <- ord[!is.na(ord)]
    
    if(length(ord) == 0){
      sps_occ$order <- NA
    } else {
      sps_occ$order <- ord[1]
    }
  }
  
  ## Calculate the range position (mean latitude)
  sps_occ$rangeLoc <- abs((ymax + ymin) / 2)
  
  ## Calculate the latitudinal amplitude of the range
  sps_occ$latAmpl <- abs(ymax - ymin)
  
  #save table
  setwd(wd_pts_metrics)
  write.csv(sps_occ, paste0(sps,'_point_range_metrics.csv'), row.names = F)
}
