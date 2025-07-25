####  This script calculates the percentage of each species IUCN range map in 
####  each biome

#load libraries
library(sf); library(units)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_tab_SI <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/SI/Tables'

#load biomes shapefile
biomes <- st_read(dsn = wd_biomes, layer = 'wwf_terr_ecos')

## visually identify each biome and include name in new column ##

#create new column for biomes names
biomes$BIOME_NAME <- NA

#visualise and include name
plot(st_geometry(biomes[biomes$BIOME == 1,]), col = 'red', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 1] <-
  'Tropical and Subtropical Moist Broadleaf Forest'

plot(st_geometry(biomes[biomes$BIOME == 2,]), col = 'blue', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 2] <-
  'Tropical and Subtropical Dry Broadleaf Forest'

plot(st_geometry(biomes[biomes$BIOME == 3,]), col = 'green', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 3] <-
  'Tropical and Subtropical Dry Coniferous Forest'

plot(st_geometry(biomes[biomes$BIOME == 4,]), col = 'orange', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 4] <- 'Temperate Broadleaf and Mixed Forests'

plot(st_geometry(biomes[biomes$BIOME == 5,]), col = 'magenta', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 5] <- 'Temperate Coniferous Forest'

plot(st_geometry(biomes[biomes$BIOME == 6,]), col = 'cyan', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 6] <- 'Boreal Forest/Taiga'

plot(st_geometry(biomes[biomes$BIOME == 7,]), col = 'grey', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 7] <-
  'Tropical and Subtropical Grasslands, Savannas, and Shrublands'

plot(st_geometry(biomes[biomes$BIOME == 8,]), col = 'yellow', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 8] <-
  'Temperate Grasslands, Savannas, and Shrublands'

plot(st_geometry(biomes[biomes$BIOME == 9,]), col = 'pink', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 9] <- 'Flooded Grasslands and Savannas'

plot(st_geometry(biomes[biomes$BIOME == 10,]), col = 'purple', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 10] <- 'Montane Grasslands and Shrublands'

plot(st_geometry(biomes[biomes$BIOME == 11,]), col = 'gold', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 11] <- 'Tundra'

plot(st_geometry(biomes[biomes$BIOME == 12,]), col = 'darkgreen', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 12] <-
  'Mediterranean Forests, Woodlands, and Scrub'

plot(st_geometry(biomes[biomes$BIOME == 13,]), col = 'lightblue', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 13] <- 'Deserts and Xeric Shrublands'

plot(st_geometry(biomes[biomes$BIOME == 14,]), col = 'violet', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 14] <- 'Mangroves'

plot(st_geometry(biomes[biomes$BIOME == 98,]), col = 'darkorange', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 98] <- 'Lakes'

plot(st_geometry(biomes[biomes$BIOME == 99,]), col = 'darkcyan', border = NA)
biomes$BIOME_NAME[biomes$BIOME == 99] <- 'Rock and Ice'

#list species
setwd(wd_tab_SI)
used_sps_nOcc <- read.csv('Used_sps_nOcc.csv')
sps_list <- unique(used_species_nOcc$species)

#list biomes
biomes_list <- unique(biomes$BIOME_NAME)

#create an empty matrix with ncol == n biomes, nrow == n sps
results <- matrix(data = NA, nrow = length(sps_list), ncol = length(biomes_list))

#loop through species overlapping their range maps with biomes
for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #select the range representing only the native range of the species
  range <- range[range$legend == 'Extant (resident)',]
  
  #unify all features
  sps_range2 <- st_union(range)
  sps_range2 <- st_make_valid(sps_range2)
  
  #transform projections to deal with sf usual probles
  old_proj <- crs(biomes)
  biomes <- st_transform(biomes, crs = 3857)
  sps_range2 <- st_transform(sps_range2, crs = 3857)
  
  #calculate range area (km^2)
  sps_area <- set_units(st_area(sps_range2), km^2)
  
  #clip biomes shp by sps range
  biomes_clipped <- st_intersection(biomes, sps_range2)
  
  #calculate percentage of range in each biome
  biomes_sps <- unique(biomes_clipped$BIOME_NAME)
  
  #loop though all biomes and calculate %
  for(j in 1:length(biomes_list))
  {
    if(biomes_list[j] %in% biomes_sps){
      a <- biomes_clipped[biomes_clipped$BIOME_NAME == biomes_list[j],]
      b <- st_union(a)
      c <- set_units(st_area(b), km^2)
      results[i,j] <- as.numeric(c / sps_area) * 100
    }else{
      results[i,j] <- 0
    }
  }
  
  #reproject
  biomes <- st_transform(biomes, crs = old_proj)
  
  print(i)
}

#transform matrix into dataframe
results2 <- as.data.frame(results)

#name rows and cols
row.names(results2) <- sps_list
colnames(results2) <- biomes_list

#test
row_sums <- rowSums(results2)

#save results
setwd(wd_tab_SI)
write.csv(results2, 'Species_biome_percentage.csv')

#check how many species are represented in each biome
results3 <- results2
results3[] <- ifelse(results3[] > 0, 1, 0)
results3 <- colSums(results3)

summary(row_sums)
min(row_sums)
max(row_sums)
