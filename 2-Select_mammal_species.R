#load packages
library(sf); library(rworldmap)

#tells sf to use the old GEOS planar engine to avoid edge crossing error
sf_use_s2(FALSE)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_lists <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'

#install functions to be used in this script

crop_world_to_range <- function(range, world, margin_deg = 1) {
  # range: sf species range in lon/lat (WGS84)
  # world: sf world land polygons in lon/lat
  # margin_deg: degrees to expand bounding box (~1 degree ≈ 100 km)
  
  # 1. Bounding box of the species, in local metres
  bb <- st_bbox(range)
  
  # 2. Expand the longitude of bbox by the margin
  bb["xmin"] <- bb["xmin"] - margin_deg
  bb["xmax"] <- bb["xmax"] + margin_deg
  
  # 3. Expand the latitude of bbox by the margin (needs clamping)
  bb["ymin"] <- max(bb["ymin"] - margin_deg, -90)
  bb["ymax"] <- min(bb["ymax"] + margin_deg,  90)
  
  # 4. Turn expanded bbox into an sf geometry (polygon)
  bb_sfc <- st_as_sfc(bb)
  st_crs(bb_sfc) <- st_crs(world)
  
  # 5. Crop world to this expanded bbox
  world_crop <- suppressWarnings(st_intersection(world, bb_sfc))
  
  return(world_crop)
}

get_local_laea_crs <- function(range) {
# range: sf object (one species), in WGS84 (lon/lat)

# 1. Get geometry and union (in case it is MULTIPOLYGON)
range_geom   <- st_geometry(range)
range_union  <- st_union(range_geom)

# 2. Centroid of the range
range_centroid <- st_centroid(range_union)

# 3. Coordinates of centroid (lon, lat in degrees)
cent_coords <- st_coordinates(range_centroid)
lon0 <- cent_coords[1]   # longitude of range centre
lat0 <- cent_coords[2]   # latitude of range centre

# 4. Build local LAEA CRS string (in metres)
crs_laea <- paste0(
  "+proj=laea +lat_0=", lat0,
  " +lon_0=", lon0,
  " +datum=WGS84 +units=m +no_defs"
)

return(crs_laea)
}

compute_frag_dir <- function(part_geom, world_land, touches_water,
                             dfMaxLength = 2000, coast_buffer = 5000) {
  # part_geom: single POLYGON (here: range_parts[j,])
  # world_land: merged land polygon for the local area
  # touches_water: logical, T if fragment touches water (here: touch_vec[j])
  # dfMaxLength: max segment length (m) when densifying boundary
  # coast_buffer: dist (m) to define the "coastal band" around the fragment
  
  #if fragment does not touch water → 0% in all directions
  if (isFALSE(touches_water)) {
    return(c(N = 0, S = 0, E = 0, W = 0))
  }
  
  #densify boundary to get many points
  boundary <- st_boundary(part_geom)
  boundary_dense <- st_segmentize(boundary, dfMaxLength = dfMaxLength)
  boundary_pts <- st_cast(boundary_dense, "POINT")
  
  #coastal band around the fragment, then "ocean side" = band minus land
  band <- st_buffer(part_geom, dist = coast_buffer)
  
  #make sure band and world_land are valid
  band <- st_make_valid(band)
  world_land <- st_make_valid(world_land)
  
  #part of the band that is not land (i.e. towards the ocean)
  ocean_ring <- suppressWarnings(st_difference(band, world_land))
  coastal_logical <- lengths(st_intersects(boundary_pts, ocean_ring)) > 0
  coastal_pts <- boundary_pts[coastal_logical,]
  
  #coordinates of coastal points
  coords <- st_coordinates(coastal_pts)
  
  # Bounding box of the fragment
  bbox <- st_bbox(part_geom)
  xmin <- bbox["xmin"]
  xmax <- bbox["xmax"]
  ymin <- bbox["ymin"]
  ymax <- bbox["ymax"]
  
  #### visualise bbox ####
  # bbox_sfc <- sf::st_as_sfc(bbox)
  # bbox_sfc <- sf::st_set_crs(bbox_sfc, sf::st_crs(part_geom))
  # bbox_sf <- sf::st_sf(geometry = bbox_sfc)
  # plot(bbox_sf, add = T, col = NA, border = 'orange')
  #### visualise bbox ####
  
  #calc distances to each side of the bbox
  dN <- ymax - coords[, "Y"]
  dS <- coords[, "Y"] - ymin
  dE <- xmax - coords[, "X"]
  dW <- coords[, "X"] - xmin
  
  dist_mat <- cbind(dN, dS, dE, dW)
  
  #assign each coastal point to the nearest edge (1 = N, 2 = S, 3 = E, 4 = W)
  nearest_edge <- apply(dist_mat, 1, which.min)
  
  nN <- sum(nearest_edge == 1L)
  nS <- sum(nearest_edge == 2L)
  nE <- sum(nearest_edge == 3L)
  nW <- sum(nearest_edge == 4L)
  
  total <- length(nearest_edge)
  
  #get percentage
  result <- c(N = 100 * nN / total,
              S = 100 * nS / total,
              E = 100 * nE / total,
              W = 100 * nW / total)
  
  return(result)
}


#list species
setwd(wd_ranges)
sps_list <- gsub('.shp', '', list.files(pattern = '.shp'))

#load world map
world <- getMap()

#change world to sf
world <- st_as_sf(world)

#create vactors to store info about ranges
contiguous <- character() # "yes"/"no"
fragments <- numeric() # number of fragments
border_ocean <- character() # "yes"/"no" 
touch_frags <- numeric()  # number of fragments touching water
water_geoms <- character()  # list of water geometries per species (per part)
range_frag_km2 <- character() # area in km2 of each fragment
N_ocean <- character() # percentage of range limited by ocean to the North
S_ocean <- character() # percentage of range limited by ocean to the South
E_ocean <- character() # percentage of range limited by ocean to the East
W_ocean <- character() # percentage of range limited by ocean to the West

#loop through ranges
for(i in 1:length(sps_list))
{
  #load species range
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #make sure the object has only one feature
  if(nrow(range) > 1){
    stop('More than one feature')
  }
  
  #extract the geometry for the first (and only) feature
  g <- st_geometry(range)[[1]]
  
  #number of polygon parts in this MULTIPOLYGON
  n_parts <- length(g)
  
  #crop world by the extent of the range plus a 1 degree buffer
  world_crop <- crop_world_to_range(range, world)
  
  #define the best CRS for the species (in metres based on location)
  crs_sps <- get_local_laea_crs(range)
  
  #transform shp range and world to the specific crs
  range_loc <- st_transform(range, crs_sps)
  world_loc <- st_transform(world_crop, crs_sps)
  
  #clean land polygons locally
  world_loc  <- suppressWarnings(st_make_valid(world_loc))
  
  #merge all land polygons into a single geometry (ignore country borders)
  world_land <- suppressWarnings(st_union(world_loc))
  
  #nsure the merged landmass is also valid
  world_land <- suppressWarnings(st_make_valid(world_land))
  
  #split range into individual polygon parts (fragments)
  range_parts <- st_cast(range_loc, "POLYGON")
  n_parts_loc <- nrow(range_parts)  # should match n_parts, but we use this here
  
  #create temporary vectors to loop through the fragments
  touch_vec <- rep(NA, n_parts_loc)
  frag_areas_km2 <- numeric()
  
  #loop through each range portion to check whether they are limited by ocean
  for(j in 1:n_parts_loc)
  {
    #get each part of the range
    this_part <- range_parts[j,]
    
    #total area of this fragment (m^2 -> km^2)
    area_part <- as.numeric(st_area(this_part))
    frag_areas_km2[j] <- area_part / 1e6
    
    #5 km buffer around this fragment
    buf5 <- suppressWarnings(st_buffer(this_part, dist = 5000))
    
    # intersection of buffer with landmass
    land_buf <- suppressWarnings(st_intersection(buf5, world_land))
    
    # area of the buffer itself (this is NOT the fragment area)
    area_buf <- as.numeric(st_area(buf5))

    if(nrow(land_buf) == 0) {
      
      # 5 km neighbourhood has no land at all -> clearly near water
      touch_vec[j] <- TRUE
      
    } else {
      #calculate percentage of the fragment that is safely inland
      area_land_buf <- sum(as.numeric(st_area(land_buf)), na.rm = TRUE)
      
      #just a fix for tiny differences in calculations
      diff <- abs(area_buf - area_land_buf) #difference 
      tol <- 1 #tolerance: tiny area (1 m2)
      
      # That means the fragment "touches water", unless the diff is within tol
      
      if(diff > tol) {
        touch_vec[j] <- TRUE
      } else {
        touch_vec[j] <- FALSE
      }
    }
  }
  
  #minimum fragment size to keep (in km^2)
  min_frag_km2 <- 1
  
  #which fragments to keep
  keep <- frag_areas_km2 >= min_frag_km2
  
  #number of fragments to keep
  n_keep <- sum(keep)
  
  #check which fragments touch water
  idx_touch <- which(touch_vec[keep])
  
  if(length(idx_touch) == 0) {
    water_str <- "none"
  } else if(length(idx_touch) == n_parts_loc) {
    water_str <- "all"
  } else {
    water_str <- paste(idx_touch, collapse = ",")
  }
  
  #check whether any portions of the range are limited by ocean
  idx_keep <- which(keep)  
  idx_touch_keep <- idx_keep[idx_touch]
  
  #update objects eliminating what we do not wanna keep
  range_parts <- range_parts[keep,]   # keep only valid polygons
  frag_areas_km2 <- frag_areas_km2[keep]
  touch_vec <- touch_vec[keep]
  
  #loop through each range portion to check whether they are limited by ocean
  
  
  #populate vectors
  contiguous[i] <- ifelse(n_keep == 1, 'yes', 'no')
  fragments[i] <- n_keep
  border_ocean[i] <- ifelse(any(touch_vec), "yes", "no")
  touch_frags[i] <- sum(touch_vec)
  water_geoms[i] <- water_str
  range_frag_km2[i] <- paste(round(frag_areas_km2, 2), collapse = ";")
  
  
  N_ocean[i] <- 
  S_ocean[i] <- 
  E_ocean[i] <- 
  W_ocean[i] <- 
  
}

##### visualise #####

#see region on world
plot(st_geometry(world), col = 'grey90')
plot(st_geometry(world_crop), col = 'red', add = T)

#see range on region
plot(st_geometry(world_land), col = 'red')
plot(st_geometry(range_loc), col = 'green', add = T)

#see each part
plot(st_geometry(buf5), col = 'magenta', border = 'purple', add = T)

plot(st_geometry(buf50), col = 'yellow', border = 'purple', add = T)


# #borders that do not disappear in the polygon merging
# plot(st_geometry(world_loc), col = 'grey')
# plot(st_geometry(world_land), col = NA, border = 'green', add = T)





##############


#create a data frame for all species
res <- data.frame()
