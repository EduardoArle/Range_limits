#load packages 
library(raster); library(sp); library(data.table); library(bRacatus)
library(CoordinateCleaner)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_variables <- '/Users/carloseduardoaribeiro/Documents/CNA/Data/Variables/wc2-5'
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_lists <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'
wd_maps <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Maps'
wd_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results'
wd_elevation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/Elevation_Tozer'

#load species list (chose with list from the folder)
setwd(wd_lists)
sps_list <- readRDS('Selected_mammal_species_12')

#load one variable to thin the records by the same resolution
setwd(wd_variables)
mean_d_range <- raster('bio2.bil')

#make ID raster to thin points by variable resolution
ID_raster <- mean_d_range
ID_raster[] <- c(1:length(ID_raster)) #has to be run twice to work for some zarbi reason

#create vectors to make a table keeping track of the species used and numbers of records
n_GBIF <- numeric()
n_recs_range <- numeric()
n_thinned_recs <- numeric()

### did not work c(19, 21)

#download and thin point records from GBIF
for(i in 1:length(sps_list))
{
  #download occ from GBIF using bRacatus
  sps_occ <- getOcc(sps_list[i])
  
  #eliminate entries without coordinates
  sps_occ_coord <- sps_occ[complete.cases(sps_occ$decimalLongitude),] 
  sps_occ_coord <- sps_occ_coord[complete.cases(sps_occ_coord$decimalLatitude),] 
  
  #populate n_GBIF column
  n_GBIF[i] <- nrow(sps_occ_coord)
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #select the range representing only the native range of the species
  range <- range[range$legend == 'Extant (resident)',]
  
  #make the range only one object, in case it is composed of more
  range <- st_union(range)
  
  #create spatial object
  sps_occ_sp <- st_as_sf(sps_occ_coord, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(range))
  
  #select only the points within the species range
  pts_range <- st_intersects(sps_occ_sp, range, sparse = F) 
  sps_occ_range <- sps_occ_coord[pts_range,] #table
  sps_occ_range_sp <- sps_occ_sp[pts_range,] #spatial object
  
  #populate n_recs_range column
  n_recs_range[i] <- nrow(sps_occ_range)
  
  #get cellID value to thin the records to max one per grid cel (variables)
  cellIDs <- extract(ID_raster, sps_occ_range_sp)
  
  #include cell value and coordinates into species data
  sps_occ_range$cellID <- cellIDs
  
  #keep only one entry per cell
  sps_occ_thin <- unique(as.data.table(sps_occ_range), by = 'cellID')
  
  #check if there are flagged occurrences with CoordinateCleaner
  sps_occ_thin_CC <- clean_coordinates(sps_occ_thin,
                              tests = c("capitals", "centroids", "equal",
                                        "gbif", "institutions", "outliers",
                                        "seas", "zeros", "duplicates"))
  
  #keep only records not flagged by CoordinateCleaner
  sps_occ_thin <- sps_occ_thin_CC[sps_occ_thin_CC$.val == T,]
  
  #populate n_thinned_recs column
  n_thinned_recs[i] <- nrow(sps_occ_thin)
  
  #save final table table
  setwd(wd_thinned_occ)
  write.csv(sps_occ_thin, paste0(sps_list[i], '_thinned.csv'), row.names = F)
  
  #save map
  setwd(wd_maps)
  
  jpeg(file=paste0(sps_list[i],"_occurrence.jpeg"))
  
  plot(world, border = NA, col = 'orange')
  plot(range, add = T, border = NA, col = 'green')
  plot(st_geometry(sps_occ_range_sp), 
       col = 'magenta', cex = 0.2, pch = 19, add = T)
  
  dev.off()
}

# make a table with the intermediate steps in data filtering

res <- data.frame(Species = sps_list,
                  n_GBIF = n_GBIF,
                  n_recs_range = n_recs_range,
                  n_thinned_recs = n_thinned_recs)


setwd(wd_res)
write.csv(res, 'Species_occurrence_4.csv', row.names = F)

######


jpeg(file=paste0(sps_list[i],"_var_contribution.jpeg"))

#plot part of the world as background
plot(world, border = NA, col = 'gray90')

#plot part of the world as background
plot(region, border = NA, col = 'gray90')

#plot points with shap result 
plot(sps_sp_pr, add = T, col = sps_sp_pr$colour, pch =19)

#make legend
points(leg_pos[1] + inter_x, 
       leg_pos[2] + inter_y, pch = 19, cex = 2, col = 'blue')
text(leg_pos[1] + (inter_x * 1.5), 
     leg_pos[2] + inter_y, 'Precipitation Seasonality', 
     cex = 1.3, pos = 4)

points(leg_pos[1] + inter_x, 
       leg_pos[2] + (2 * inter_y), pch = 19, cex = 2, col = 'orange')
text(leg_pos[1] + (inter_x * 1.5), 
     leg_pos[2] + (2 * inter_y), 'Temperature Annual Range', 
     cex = 1.3, pos = 4)

points(leg_pos[1] + inter_x, 
       leg_pos[2] + (3 * inter_y), pch = 19, cex = 2, col = 'red')
text(leg_pos[1] + (inter_x * 1.5), 
     leg_pos[2] + (3 * inter_y), 'Mean Diurnal Range', 
     cex = 1.3, pos = 4)

points(leg_pos[1] + inter_x, 
       leg_pos[2] + (4 * inter_y), pch = 19, cex = 2, col = 'yellow')
text(leg_pos[1] + (inter_x * 1.5), 
     leg_pos[2] + (4 * inter_y), 'Temperature Seasonality', 
     cex = 1.3, pos = 4)

#plot range on top of the world map
plot(sps_range2, add = T, border = 'darkgreen', lwd = 2)

#plot species name
text((leg_pos[1] + leg_pos[3]) / 2, 
     leg_pos[4] - (2 * inter_y), 
     paste0(sps_list[i]),
     font = 3,
     cex = 1.3)

dev.off()

library(rworldmap)
world <-  getMap()

plot(world, border = NA, col = 'orange')
plot(range, add = T, border = NA, col = 'green')
plot(st_geometry(sps_occ_range_sp), col = 'magenta', cex = 0.2, pch = 19, add = T)

plot(st_geometry(sps_occ_sp), col = 'magenta', cex = 0.2, pch = 19, add = T)


#####


#thin points and save 
for(i in 1:length(sps_list))
{

  

  
  #save results
  setwd(wd_thinned_occ)
  write.csv(c, paste0(sps,'.csv'), row.names = F)
  
  print(i)
}


