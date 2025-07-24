#load packages
library(sf); library(data.table)

#list WDs
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Tables'
wd_all_sps_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_biomes <- '/Users/carloseduardoaribeiro/Documents/General data/Biomes/official'
wd_map_stuff <- '/Users/carloseduardoaribeiro/Documents/Collaborations/Rachel/Fogo'


#load file with all results
setwd(wd_all_sps_res)
all_res <- read.csv('20241210_Results_all_sps.csv')

#select only species used in any model
all_res_models <- all_res[which(!is.na(all_res$avg_Min_T_SHAP) |
                                !is.na(all_res$avg_Mean_T_SHAP) |
                                !is.na(all_res$avg_Max_T_SHAP)),]

n_sps_model <- length(unique(all_res_models$species))

#get the total number of points that were used in any models 
n_pts_models <- nrow(all_res_models)

#get the total number of species for which any model was made
n_species_all <- length(unique(all_res$species))
n_species_all

#create a version of the file including only entries that were not too correlated
all_res_cor <- all_res_models[which(abs(all_res_models$Cor_vars_minT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_meanT) <= 0.7 |
                                    abs(all_res_models$Cor_vars_maxT) <= 0.7),]

#n points that were used valid models (by var correl)
n_pts_val_models <- nrow(all_res_cor)
n_pts_val_models

#get the total number of species for which any model was made
n_species_model <- length(unique(all_res_cor$species))
n_species_model

#get info on biomes represented
length(unique(all_res_cor$biome))

#get number of species for which each model was made 


#################. FIIXXXXX ##############

# not !is.na, but by the correl


all_res_cor_minT <- all_res_cor[which(!is.na(all_res_cor$avg_Min_T_SHAP)),]
all_res_cor_meanT <- all_res_cor[which(!is.na(all_res_cor$avg_Mean_T_SHAP)),]
all_res_cor_maxT <- all_res_cor[which(!is.na(all_res_cor$avg_Max_T_SHAP)),]
all_res_cor_minPPT <- all_res_cor[which(!is.na(all_res_cor$avg_Min_PPT_SHAP)),]
all_res_cor_meanPPT <- all_res_cor[which(!is.na(all_res_cor$avg_Mean_PPT_SHAP)),]
all_res_cor_maxPPT <- all_res_cor[which(!is.na(all_res_cor$avg_Max_PPT_SHAP)),]

n_sps_minT <- length(unique(all_res_cor_minT$species))
n_sps_meanT <- length(unique(all_res_cor_meanT$species))
n_sps_maxT <- length(unique(all_res_cor_maxT$species))
n_sps_minPPT <- length(unique(all_res_cor_minPPT$species))
n_sps_meanPPT <- length(unique(all_res_cor_meanPPT$species))
n_sps_maxPPT <- length(unique(all_res_cor_maxPPT$species))

n_sps_minT
n_sps_meanT
n_sps_maxT
n_sps_minPPT
n_sps_meanPPT
n_sps_maxPPT





###########################################################################

#load mistressfile
setwd(wd_tables)
mistress_file <-  read.csv('Mistressfile.csv')

##### Count number of species with ok correl for each model #####

minT_meanPPT <- t(table(abs(mistress_file$Cor_vars_minT) <= 0.7))
meanT_meanPPT <- t(table(abs(mistress_file$Cor_vars_meanT) <= 0.7))
maxT_meanPPT <- t(table(abs(mistress_file$Cor_vars_maxT) <= 0.7))

minPPT_meanT <- t(table(abs(mistress_file$Cor_vars_minPPT) <= 0.7))
meanPPT_meanT <- t(table(abs(mistress_file$Cor_vars_meanPPT) <= 0.7))
maxPPT_meanT <- t(table(abs(mistress_file$Cor_vars_maxPPT) <= 0.7))

#count number of folds per species
n_sps_models <- sum(table(mistress_file$folds))

table(mistress_file$folds)
sum(table(mistress_file$folds))

#get minimum number of records that produced a model
  
#select only the rows in the mistress file that produced at least one model
miss_file_folds <- mistress_file[which(!is.na(mistress_file$folds)),]

#get min num occ that produced models
min(miss_file_folds$n_occ)
summary(miss_file_folds$n_occ)

#test 1

#select only the rows in the mistress file that produced at least one model
miss_file_good_models <- miss_file_folds[
  as.numeric(gsub(' %', '', mistress_file$Used_models_minT)) > 0 | 
  as.numeric(gsub(' %', '', mistress_file$Used_models_meanT)) > 0 | 
  as.numeric(gsub(' %', '', mistress_file$Used_models_maxT)) > 0 | 
  as.numeric(gsub(' %', '', mistress_file$Used_models_minPPT)) > 0 | 
  as.numeric(gsub(' %', '', mistress_file$Used_models_meanPPT)) > 0 | 
  as.numeric(gsub(' %', '', mistress_file$Used_models_maxPPT)) > 0 ,]

#eliminate NAs
miss_file_good_models <- miss_file_good_models[
  complete.cases(miss_file_good_models$Species),]

n_sps_good_model <- nrow(miss_file_good_models)

min_n_occ_good_model <- min(miss_file_good_models$n_occ)

summary(miss_file_good_models$n_occ)

#test 2

#select only the rows in the mistress file that produced at least one model
miss_file_good_models2 <- miss_file_folds[
  mistress_file$AUC_sel_minT >= 0.7 & mistress_file$TSS_sel_minT >= 0.4| 
  mistress_file$AUC_sel_meanT >= 0.7 & mistress_file$TSS_sel_meanT >= 0.4| 
  mistress_file$AUC_sel_maxT >= 0.7 & mistress_file$TSS_sel_maxT >= 0.4| 
  mistress_file$AUC_sel_minPPT >= 0.7 & mistress_file$TSS_sel_minPPT >= 0.4| 
  mistress_file$AUC_sel_meanPPT >= 0.7 & mistress_file$TSS_sel_meanPPT >= 0.4| 
  mistress_file$AUC_sel_maxPPT >= 0.7 & mistress_file$TSS_sel_maxPPT >= 0.4,]

#eliminate NAs
miss_file_good_models2 <- miss_file_good_models2[
  complete.cases(miss_file_good_models2$Species),]

n_sps_good_model2 <- nrow(miss_file_good_models2)

min_n_occ_good_model2 <- min(miss_file_good_models2$n_occ)

summary(miss_file_good_models2$n_occ)



##### Count biomes #####

#load biomes Shapefile
biomes_shp <- st_read(dsn = wd_biomes, 'wwf_terr_ecos')

#load map stuff
setwd(wd_map_stuff)
world <- readRDS("wrld.rds")
worldmapframe <- readRDS("Worldmapframe.rds")

#transform SpatialPolygons into sf object
worldmapframe <- st_as_sf(worldmapframe)

#transform the sf object and frame to the nice projection
biomes_shp2 <- st_transform(biomes_shp, crs = st_crs(world))
worldmapframe <- st_transform(worldmapframe, crs = st_crs(world))

#check unique BIOME entries in the shp
sort(unique(biomes_shp$BIOME))

#plot all features corresponding to each biome number


#make a relation table to connect biome numbers to biome names
    # I can't find a table with this info, so I made it manually based on the 
    # shapefile compared to Fig. 1 in Olson et al. 2001

#biome 1
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 1)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 1),]),
     col = '#008939', border = NA, add = F)

#biome 2 = Tropical and Subtropical Dry Broadleaf Forests
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 2)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 2),]),
     col = '#c9b600', border = NA, add = T)

#biome 3 = Tropical and Subtropical Coniferous Forests
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 3)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 3),]),
     col = '#8dd02c', border = NA, add = T)

#biome 4 = Temperate Broadleaf and Mixed Forests
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 4)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 4),]),
     col = '#007860', border = NA, add = T)

#biome 5 = Temperate Coniferous Forests
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 5)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 5),]),
     col = '#007288', border = NA, add = T)

#biome 6 = Boreal Forest / Taiga
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 6)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 6),]),
     col = '#6bc59d', border = NA, add = T)

#biome 7 = Tropical and Subtropical Grasslands, Savannas, and Shrublands
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 7)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 7),]),
     col = '#ff9f00', border = NA, add = T)

#biome 8 = Temperate Grasslands, Savannas, and Shrublands
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 8)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 8),]),
     col = '#ffd000', border = NA, add = T)

#biome 9 = Flooded Grasslands and Savannas
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 9)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 9),]),
     col = '#28d2c2', border = NA, add = T)

#biome 10 = Montane Grasslands and Shrublands
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 10)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 10),]),
     col = '#d6a36c', border = NA, add = T)

#biome 11 = Tundra
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 11)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 11),]),
     col = '#b3df8f', border = NA, add = T)

#biome 12 = Mediterranean Forests, Woodlands, and Scrub
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 12)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 12),]),
     col = '#ff0000', border = NA, add = T)

#biome 13 = Deserts and Xeric Shrublands
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 13)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 13),]),
     col = '#ff6e4e', border = NA, add = T)

#biome 14 = Mangroves
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 14)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 14),]),
     col = '#ff0089', border = NA, add = T)

#biome 98 = Lake
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 98)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 98),]),
     col = '#ffffff', border = NA, add = T)

#biome 99 = Rock and Ice
biomes_shp2$ECO_NAME[which(biomes_shp2$BIOME == 99)]
plot(st_geometry(biomes_shp2[which(biomes_shp2$BIOME == 99),]),
     col = '#dddddd', border = NA, add = T)

#plot map frame
plot(worldmapframe, add = T, lwd = 10, border = '#101010')

#save map width = 10000

#load results with all species table (including biome of each data point)
setwd(wd_all_sps_res)
all_sps_table <-  read.csv('20241210_Results_all_sps.csv')

#select only species that yielded at least one useful model (no NAs in SHAP)
sel_table_1 <- all_sps_table[-which(is.na(all_sps_table$avg_Min_T_SHAP) &
                                        is.na(all_sps_table$avg_Mean_T_SHAP) &
                                        is.na(all_sps_table$avg_Max_T_SHAP) &
                                        is.na(all_sps_table$avg_Min_PPT_SHAP) &
                                        is.na(all_sps_table$avg_Mean_PPT_SHAP) &
                                        is.na(all_sps_table$avg_Max_PPT_SHAP)),]

#select only species that yielded at least one useful model (vars not too correl)
#first eliminate rows where all correlations are NA (no models)
sel_table_2 <- sel_table_1[-which(is.na(sel_table_1$Cor_vars_minT) &
                                  is.na(sel_table_1$Cor_vars_meanT) &
                                  is.na(sel_table_1$Cor_vars_maxT) &
                                  is.na(sel_table_1$Cor_vars_minPPT) &
                                  is.na(sel_table_1$Cor_vars_meanPPT) &
                                  is.na(sel_table_1$Cor_vars_maxPPT)),]

#first eliminate rows where all correlations are higher than |0.7|
sel_table_3 <- sel_table_2[-which(abs(sel_table_2$Cor_vars_minT) >= 0.7 &
                                 abs(sel_table_2$Cor_vars_meanT) >= 0.7 &
                                 abs(sel_table_2$Cor_vars_maxT) >= 0.7 &
                                 abs(sel_table_2$Cor_vars_minPPT) >= 0.7 &
                                 abs(sel_table_2$Cor_vars_meanPPT) >= 0.7 &
                                 abs(sel_table_2$Cor_vars_maxPPT) >= 0.7),]

#count points per biome
pts_biome <- table(sel_table_3$biome)

#count species per biome
sps_biome <- unique(as.data.table(sel_table_3), by = c('species', 'biome'))
sps_biome2 <- table(sps_biome$biome)

#count species
n_sps_sel <- length(unique(sel_table_3$species))

sps_l <- unique(sel_table_3$species)

sps_l %in% miss_file_good_models2$Species


head(all_sps_table)

table(all_sps_table$biome)

