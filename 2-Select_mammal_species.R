#load packages
library(sf); library(rworldmap); library(rworldxtra)

#tells sf to use the old GEOS planar engine to avoid edge crossing error
sf_use_s2(FALSE)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_plots <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Plots_ranges"
wd_tables <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables"

#install functions to be used in this script

crop_world_to_range <- function(range, world, margin_deg = 1){
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

get_local_laea_crs <- function(range){
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
                             dfMaxLength = 2000, coast_buffer = 5000){
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
  
  #% of boundary that touches water
  n_all <- nrow(boundary_pts)
  n_coastal <- sum(coastal_logical)
  prop_coast <- 100 * n_coastal / n_all  
  
  #keep only coastal points for directionality
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
  X <- coords[, 1]
  Y <- coords[, 2]
  
  dN <- ymax - Y
  dS <- Y - ymin
  dE <- xmax - X
  dW <- X - xmin
  
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
              W = 100 * nW / total,
              P = prop_coast)
  
  return(result)
}

#list species
setwd(wd_ranges)
sps_list <- gsub('.shp', '', list.files(pattern = '.shp'))

#load world map
world <- getMap(resolution = 'high') #high for analyses
world_low <- getMap() #low for faster plots

#change world to sf
world <- st_as_sf(world)
world_low <- st_as_sf(world_low)

#create vactors to store info about ranges
contiguous <- character() # "yes"/"no"
fragments <- numeric() # number of fragments
border_ocean <- character() # "yes"/"no" 
touch_frags <- numeric()  # number of fragments touching water
water_geoms <- character()  # list of water geometries per species (per part)
range_frag_km2 <- character() # area in km2 of each fragment
Perc_ocean <- character() # percentage of frag limited by ocean
N_ocean <- character() # percentage of frag limited by ocean to the North
S_ocean <- character() # percentage of frag limited by ocean to the South
E_ocean <- character() # percentage of frag limited by ocean to the East
W_ocean <- character() # percentage of frag limited by ocean to the West

#loop through ranges
for(i in 1:543)#length(sps_list))
{
  #load species range
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #make sure the object has only one feature
  range <- suppressWarnings(st_union(range))
  
  #extract the geometry for the first (and only) feature
  g <- st_geometry(range)[[1]]
  
  #number of polygon parts in this MULTIPOLYGON
  n_parts <- length(g)
  
  #define the best CRS for the species (in metres based on location)
  crs_sps <- get_local_laea_crs(range)
  
  #transform shp range to the specific crs
  range_loc <- st_transform(range, crs_sps)
  
  #split range into individual polygon parts (fragments)
  range_parts <- st_cast(st_as_sf(range_loc), "POLYGON")
  n_parts_loc <- nrow(range_parts)  # should match n_parts
  
  #------------------------------------------------------------
  # 1) fragment areas and size filter (before any world operations)
  #------------------------------------------------------------
  
  # area of each fragment (m^2 -> km^2)
  frag_areas_km2 <- as.numeric(st_area(range_parts)) / 1e6
  
  #minimum fragment size to keep (in km^2)
  min_frag_km2 <- 1
  
  #which fragments to keep
  keep <- frag_areas_km2 >= min_frag_km2
  
  #number of fragments to keep
  n_keep <- sum(keep)
  
  #if not fragments passes the size filter, skip water + direction
  if(n_keep == 0) {
    contiguous[i]     <- NA
    fragments[i]      <- 0
    border_ocean[i]   <- "no"
    touch_frags[i]    <- 0
    water_geoms[i]    <- "none"
    range_frag_km2[i] <- NA
    N_ocean[i] <- NA
    S_ocean[i] <- NA
    E_ocean[i] <- NA
    W_ocean[i] <- NA
    
    next
  }
  
  # keep only large enough fragments
  range_parts <- range_parts[keep,]
  frag_areas_km2 <- frag_areas_km2[keep] 
  
  #------------------------------------------------------------
  # 2) world cropping & transformation (only if n_keep > 0)
  #------------------------------------------------------------
  
  #crop world by the extent of the range plus a 1 degree buffer
  world_crop <- crop_world_to_range(range, world, margin_deg = 1)
  
  #transform shp of the world to the specific crs
  world_loc <- st_transform(world_crop, crs_sps)
  
  #clean land polygons locally
  world_loc  <- suppressWarnings(st_make_valid(world_loc))
  
  #merge all land polygons into a single geometry (ignore country borders)
  world_land <- suppressWarnings(st_union(world_loc))
  
  #ensure the merged landmass is also valid
  world_land <- suppressWarnings(st_make_valid(world_land))
  
  #------------------------------------------------------------
  # 3) j-loop ONLY on kept fragments to check water-touch
  #------------------------------------------------------------
  
  #create temporary vectors to loop through the fragments
  touch_vec <- rep(NA, n_keep)

  #loop through each range portion to check whether they are limited by ocean
  for(j in 1:n_keep)
  {
    #get each part of the range
    this_part <- range_parts[j,]
    
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
      area_land_buf <- sum(as.numeric(st_area(land_buf)))
      
      #just a fix for tiny differences in calculations
      diff <- abs(area_buf - area_land_buf) #difference 
      tol <- 1 #tolerance: tiny area (1 m2)
      
      #unless the diff is within tol, the fragment "touches water"
      if(diff > tol) {
        touch_vec[j] <- TRUE
      } else {
        touch_vec[j] <- FALSE
      }
    }
  }
  
  #check which kept fragments touch water
  idx_touch <- which(touch_vec)
  
  if(length(idx_touch) == 0) {
    water_str <- "none"
  } else if(length(idx_touch) == n_keep) {
    water_str <- "all"
  } else {
    water_str <- paste(idx_touch, collapse = ";")
  }
  
  #------------------------------------------------------------
  # 4) directionality – only if some kept fragment touches water
  #------------------------------------------------------------
  
  #create temporary vectors to loop through the fragments
  P_str <- NA   # % of total perimeter that is coastal
  N_str <- NA   # % of coastal part that is North
  S_str <- NA   # % of coastal part that is South
  E_str <- NA   # % of coastal part that is East
  W_str <- NA   # % of coastal part that is West
  
  if(any(touch_vec)){
    frag_P <- numeric() 
    frag_N <- numeric()
    frag_S <- numeric()
    frag_E <- numeric()
    frag_W <- numeric()
    
    for(j in 1:n_keep) {
      this_part <- range_parts[j,]
      
      #use function compute_frag_dir on each fragment touching the ocean
      dir_vec <- compute_frag_dir(part_geom = this_part,
                                  world_land = world_land,
                                  touches_water = touch_vec[j],
                                  dfMaxLength = 2000,
                                  coast_buffer = 5000)
      
      #store direction percentages for this fragment (raw numeric values)
      frag_P[j] <- dir_vec["P"]
      frag_N[j] <- dir_vec["N"]
      frag_S[j] <- dir_vec["S"]
      frag_E[j] <- dir_vec["E"]
      frag_W[j] <- dir_vec["W"]
    }
    
    #round values and append "%" to match final output format
    P_chr <- paste0(round(frag_P), "%")
    N_chr <- paste0(round(frag_N), "%")
    S_chr <- paste0(round(frag_S), "%")
    E_chr <- paste0(round(frag_E), "%")
    W_chr <- paste0(round(frag_W), "%")
    
    #collapse per-fragment values into single semicolon-separated strings
    P_str <- paste(P_chr, collapse = ";")
    N_str <- paste(N_chr, collapse = ";")
    S_str <- paste(S_chr, collapse = ";")
    E_str <- paste(E_chr, collapse = ";")
    W_str <- paste(W_chr, collapse = ";")
  }
  
  #------------------------------------------------------------
  # 5) populate per-species vectors at the very end of the i-loop
  #------------------------------------------------------------
  
  contiguous[i] <- ifelse(n_keep == 1, "yes", "no")
  fragments[i] <- n_keep
  range_frag_km2[i] <- paste(round(frag_areas_km2, 2), collapse = ";")
  border_ocean[i] <- ifelse(any(touch_vec), "yes", "no")
  touch_frags[i] <- sum(touch_vec)
  water_geoms[i] <- water_str
  Perc_ocean[i] <- P_str
  N_ocean[i] <- N_str
  S_ocean[i] <- S_str
  E_ocean[i] <- E_str
  W_ocean[i] <- W_str
  
  #------------------------------------------------------------
  # 6) plot and save maps in PDF
  #------------------------------------------------------------
  
  #plot
  setwd(wd_plots)
  
  #create PDF
  pdf(paste0(sps_list[i], ".pdf"), width = 5, height = 5)
  
  par(mfrow = c(1,2), mar = c(0,0,0,0))
  
  ## Panel 1: region on world
  
  #make sf object of bbox of the selected region
  bb <- st_bbox(world_crop)
  bb_sfc <- st_as_sfc(bb)
  st_crs(bb_sfc) <- st_crs(world)
  
  plot(st_geometry(world_low), col = 'grey90', border = 'grey70')
  plot(bb_sfc, add = T, border = "red", lwd = 1)
  title(parse(text = paste0("italic('", sps_list[i], "')")), line = -5.5)
  
  ## Panel 2: range on region
  
  #expand the bbox by 10%
  bb2 <- st_bbox(range_parts)
  pad_x <- 0.5
  pad_y <- 0.8
  
  dx <- (bb2["xmax"] - bb2["xmin"]) * pad_x
  dy <- (bb2["ymax"] - bb2["ymin"]) * pad_y
  
  xlim <- c(bb2["xmin"] - dx, bb2["xmax"] + dx)
  ylim <- c(bb2["ymin"] - dy, bb2["ymax"] + dy)
  
  plot(st_geometry(world_land), col = 'red', xlim = xlim, ylim = ylim)
  plot(st_geometry(range_loc), col = 'green', add = T)
  
  dev.off()
  
  print(i)
}


#make data.frame with results
res <- data.frame(contiguous = contiguous,
                  fragments = fragments,
                  range_frag_km2 = range_frag_km2,
                  border_ocean = border_ocean,
                  touch_frags = touch_frags,
                  water_geoms = water_geoms,
                  perc_ocean = Perc_ocean,
                  N_ocean = N_ocean,
                  S_ocean = S_ocean,
                  E_ocean = E_ocean,
                  W_ocean = W_ocean)


setwd(wd_tables)
write.csv(res, 'Partial_results_1_543.csv', row.names = F)

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
