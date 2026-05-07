#load libraries
library(data.table); library(sf)

#list WDs
wd_all_sps_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'
wd_tab_SI <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SI/Tables'
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'

#load file with all results
setwd(wd_all_sps_res)
all_res <- read.csv('20260126_Results_all_sps.csv')

#select only species used in any model
all_res_models <- all_res[which(!is.na(all_res$avg_Min_T_SHAP) |
                                  !is.na(all_res$avg_Mean_T_SHAP) |
                                  !is.na(all_res$avg_Max_T_SHAP) |
                                  !is.na(all_res$avg_Min_PPT_SHAP) |
                                  !is.na(all_res$avg_Mean_PPT_SHAP) |
                                  !is.na(all_res$avg_Max_PPT_SHAP)),]


n_sps_model <- length(unique(all_res_models$sps))

#get the total number of points that were used in any models 
n_pts_models <- nrow(all_res_models)

#get the total number of species for which any model was made
n_species_all <- length(unique(all_res$sps))
n_species_all

#create a version of the file including only entries that were not too correlated
all_res_cor <- all_res_models[which(abs(all_res_models$Cor_vars_minT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_meanT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_maxT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_minPPT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_meanPPT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_maxPPT) <= 0.7),]

#n points that were used valid models (by var correl)
n_pts_val_models <- nrow(all_res_cor)
n_pts_val_models

#get the total number of species for which any model was made
n_species_model <- length(unique(all_res_cor$sps))
n_species_model

#get info on biomes represented
length(unique(all_res_cor$biome))

#make a table containing species used and number of records
used_species_nOcc <- unique(as.data.table(all_res_cor),
                            by = c('sps', 'nOcc'))

#select columns (sps, nOcc, rangeSize, order)
species_info <- used_species_nOcc[,c(1,6,7,16)]

#load biomes shapefiles
biomes <- st_read(dsn = wd_biomes, layer = 'wwf_terr_ecos')

#add biome names
biome_names <- data.frame(
  BIOME = c(1:14, 98, 99),
  BIOME_NAME = c(
    'Tropical and Subtropical Moist Broadleaf Forest',
    'Tropical and Subtropical Dry Broadleaf Forest',
    'Tropical and Subtropical Dry Coniferous Forest',
    'Temperate Broadleaf and Mixed Forests',
    'Temperate Coniferous Forest',
    'Boreal Forest/Taiga',
    'Tropical and Subtropical Grasslands, Savannas, and Shrublands',
    'Temperate Grasslands, Savannas, and Shrublands',
    'Flooded Grasslands and Savannas',
    'Montane Grasslands and Shrublands',
    'Tundra',
    'Mediterranean Forests, Woodlands, and Scrub',
    'Deserts and Xeric Shrublands',
    'Mangroves',
    'Lakes',
    'Rock and Ice'
  )
)

biomes$BIOME_NAME <- biome_names$BIOME_NAME[
  match(biomes$BIOME, biome_names$BIOME)
]

#remove non-terrestrial
biomes <- biomes[!biomes$BIOME %in% c(98,99), ]

#load species table
setwd(wd_all_sps_res)
sps_table <- read.csv('Selected_species.csv')

#number of species selected by range characteristics
n_sps_range_sel <- nrow(sps_table)

#get species list
sps_list <- species_info$sps

#create an empty table to store results
sps_biome <- data.frame()

#loop through species extracting biomes they overlap
for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #load species range
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps, quiet = TRUE)
  
  #keep native range
  range <- range[!range$legend %in% c('Introduced','Vagrant'), ]
  
  #make valid
  range <- st_make_valid(range)
  
  #dominant fragment size
  dom_frag_km2 <- as.numeric(strsplit(
    sps_table$range_frag_km2[
      match(sps, sps_table$species)], ';')[[1]][
        sps_table$dom_frag[match(sps, sps_table$species)]])
  
  #split range into polygon fragments
  frags <- st_cast(range, 'POLYGON')
  
  #local equal-area projection
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  
  crs_laea <- paste0('+proj=laea +lat_0=', cen[2],
                     ' +lon_0=', cen[1],
                     ' +datum=WGS84 +units=m +no_defs')
  
  #project fragments
  frags_m <- st_transform(frags, crs_laea)
  
  #calculate fragment areas in km2
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #select dominant fragment
  frag_pos <- which.min(abs(a_km2 - dom_frag_km2))
  
  #keep dominant fragment
  range_dom <- frags[frag_pos, ]
  
  #match CRS with biomes
  range_dom <- st_transform(range_dom, st_crs(biomes))
  
  #find overlapping biomes
  inter <- st_intersects(range_dom, biomes)
  
  #extract biome names
  biomes_i <- unique(biomes$BIOME_NAME[inter[[1]]])
  
  #store result
  tmp <- data.frame(sps = sps, biome = biomes_i)
  
  sps_biome <- rbind(sps_biome, tmp)
  
  print(i)
}



#############################
### reshape + merge table ###
#############################



#convert to presence-absence table
sps_biome$present <- 1

#wide format (species x biomes)
biome_wide <- reshape(sps_biome, idvar = 'species', timevar = 'biome',
                      direction = 'wide')

#replace NAs with zero
biome_wide[is.na(biome_wide)] <- 0

#clean column names
names(biome_wide) <- gsub('present.', '', names(biome_wide))

#merge with Table S1
species_info <- merge(species_info, biome_wide,
                      by.x = 'sps', by.y = 'species', all.x = TRUE)

#replace NA with zero
species_info[is.na(species_info)] <- 0

#save
setwd(wd_tab_SI)
write.csv(species_info, 'Table_S1.csv', row.names = F)

#calculate metrics
total_recs <- sum(species_info$nOcc)
summary(species_info$nOcc)



#############################
### species per biome ###
#############################



#select biome columns
biome_cols <- !(names(species_info) %in% c('sps','nOcc','rangeSize','order'))

#count species per biome
biome_counts <- colSums(species_info[, ..biome_cols])

#range of representation
range(biome_counts)
