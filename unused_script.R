##### calculate centralness from '5-Calculate metrics'

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