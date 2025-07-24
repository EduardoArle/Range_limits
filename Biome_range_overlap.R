####  This script calculates the percentage of each species IUCN range map in 
####  each biome

#load libraries
library(sf)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'

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


