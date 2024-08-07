#load packages 
library(raster);library(data.table);library(rworldmap)
library(dplyr);library(sf);library(shapviz);library(kernelshap)
library(xgboost);library(caret)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_variables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_res_species <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Comparison'

#list species 
setwd(wd_thinned_occ)
sps_list <- gsub('_thinned.csv', '', list.files())

#load all 19 BioCLim variables (only six being used anyways)
setwd(wd_variables)
AnnualMeanTemperature <- raster('wc2.1_2.5m_bio_1.tif')
MeanDiurnalRange <- raster('wc2.1_2.5m_bio_2.tif')
Isothermality <- raster('wc2.1_2.5m_bio_3.tif')
TemperatureSeasonality <- raster('wc2.1_2.5m_bio_4.tif')
MaxTemperatureOfWarmestMonth <- raster('wc2.1_2.5m_bio_5.tif')
MinTemperatureOfColdestMonth <- raster('wc2.1_2.5m_bio_6.tif')
TemperatureAnnualRange <- raster('wc2.1_2.5m_bio_7.tif')
MeanTemperatureOfWettestQuarter <- raster('wc2.1_2.5m_bio_8.tif')
MeanTemperatureOfDriestQuarter <- raster('wc2.1_2.5m_bio_9.tif')
MeanTemperatureOfWarmestQuarter <- raster('wc2.1_2.5m_bio_10.tif')
MeanTemperatureOfColdestQuarter <- raster('wc2.1_2.5m_bio_11.tif')
AnnualPrecipitation <- raster('wc2.1_2.5m_bio_12.tif')
PrecipitationOfWettestMonth <- raster('wc2.1_2.5m_bio_13.tif')
PrecipitationOfDriestMonth <- raster('wc2.1_2.5m_bio_14.tif')
PrecipitationSeasonality  <- raster('wc2.1_2.5m_bio_15.tif')
PrecipitationOfWettestQuarter <- raster('wc2.1_2.5m_bio_16.tif')
PrecipitationOfDriestQuarter <- raster('wc2.1_2.5m_bio_17.tif')
PrecipitationOfWarmestQuarter <- raster('wc2.1_2.5m_bio_18.tif')
PrecipitationOfColdestQuarter <- raster('wc2.1_2.5m_bio_19.tif')
elevation <- raster('wc2.1_2.5m_elev.tif')


######## Run SHAP for all species ######

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #load species occurrences
  setwd(wd_thinned_occ)
  sps_occ <- read.csv(paste0(sps, '_thinned.csv'))
  
  #select only presence occurrences
  pr_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'PRESENT'),]
  
  #check if there are absence data
  abs_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'ABSENT'),]
  if(nrow(abs_sps) != 0){
    warning(paste0('THERE IS ABSENT DATA FOR ', sps_list[i]))
  }
  
  #create pseudo absences (half the number of presences)
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #select the range representing only the native range of the species
  range <- range[range$legend == 'Extant (resident)',]
  
  #unify all features
  sps_range2 <- st_union(range)
  
  #make a 50km buffer around the sps range
  sps_range_buf <- st_buffer(sps_range2, 50000)  
  
  #make spatial object
  pr_sps_sp <- st_as_sf(pr_sps, coords = c('decimalLongitude', 'decimalLatitude'),
                        crs = crs(sps_range2))
  
  #make a 50km buffer around the points to delimit the area where I don't want pseudo abs
  small_buffer <- st_buffer(pr_sps_sp, 50000)  
  
  #make a spatial polygon object with only one feature
  no_pa_area <- st_union(small_buffer)
  
  # this fixes possible 'duplicate vertex' errors
  no_pa_area <- st_make_valid(no_pa_area) 
  
  #make holes in the species range by the small buffer around points
  pa_area <- st_difference(sps_range_buf, no_pa_area)
  
  #define number of pseudo abs to be created (same as presences)
  n_pa <- nrow(pr_sps)
  
  #generate pseudo absences
  pa <- st_sample(pa_area, size = n_pa, type = "random")
  
  #get coords of pa
  pa_coords <- as.data.frame(st_coordinates(pa))
  names(pa_coords) <- c('decimalLongitude', 'decimalLatitude')
  pa_coords$Occurrence <- 0 #include column informing occ status
  pa_coords$Type <- 'Pseudo-absence'
  
  #get coords of pr
  pr_coords <- as.data.frame(st_coordinates(pr_sps_sp))
  names(pr_coords) <- c('decimalLongitude', 'decimalLatitude')
  pr_coords$Occurrence <- 1 #include column informing occ status
  pr_coords$Type <- 'Presence'
  
  #combine pseudo absences and presences
  species <- rbind(pr_coords, pa_coords)
  
  #make spatial obj again
  species_sp <- st_as_sf(species, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = crs(sps_range2))
  
  #select variables we will use
  preds <- stack(AnnualMeanTemperature, AnnualPrecipitation, #means
                 MinTemperatureOfColdestMonth, PrecipitationOfDriestMonth, #mins
                 MaxTemperatureOfWarmestMonth, PrecipitationOfWettestMonth, #maxs
                 elevation) 
  
  #extract values from each location from all variables and make a table
  vals_pts <- extract(preds, species_sp)
  tab_occ_vars <- cbind(species, vals_pts)

  #select variables we will use for each predictions
  
  #elevation and mean prec to compare temperature variables
  preds_min_temp <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_6')
  preds_mean_temp  <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
  preds_max_temp <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_5')
  
  #elevation and mean temp to compare precipitation variables
  preds_min_PPT <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_14')
  preds_mean_PPT <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_12')
  preds_max_PPT <- c('wc2.1_2.5m_elev', 'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_13')
  
  # Fit XGBoost models
  
  ## mean PPT to compare T variables importances
  
  ##############################################
  ############## mean PPT + min T ##############
  ##############################################
  
  shap_T_min <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_min_temp]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_T_min,
    nrounds = 65L
  )

  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_min_temp]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_prec_SHAP <- character()  #bio12
  MIN_temp_SHAP <- character()  #bio6
  elev_SHAP <- character() #elevation

  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
    MIN_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_6']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanPPT_minT <- cbind(Species = sps_list[i],
                        tab_occ_vars[,c(1:4,6,7,11)],
                        data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                   Min_T_SHAP = MIN_temp_SHAP,
                                   elev_SHAP = elev_SHAP))

  #save results per species
  setwd(wd_res_species)
  write.csv(meanPPT_minT, paste0(sps,'_minT_.csv'), row.names = F)
  
  
  ###############################################
  ############## mean PPT + mean T ##############
  ###############################################
  
  shap_T_mean <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_mean_temp]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_T_mean,
    nrounds = 65L
  )
  
  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_mean_temp]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_prec_SHAP <- character()  #bio12
  MEAN_temp_SHAP <- character()  #bio1
  elev_SHAP <- character() #elevation
  
  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
    MEAN_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanPPT_meanT <- cbind(Species = sps_list[i],
                         tab_occ_vars[,c(1:6,11)],
                         data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                    Mean_T_SHAP = MEAN_temp_SHAP,
                                    elev_SHAP = elev_SHAP))
  
  #save results per species
  setwd(wd_res_species)
  write.csv(meanPPT_meanT, paste0(sps,'_meanT_.csv'), row.names = F)
  
  
  ##############################################
  ############## mean PPT + max T ##############
  ##############################################
  
  shap_T_max <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_max_temp]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_T_max,
    nrounds = 65L
  )
  
  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_max_temp]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_prec_SHAP <- character()  #bio12
  MAX_temp_SHAP <- character()  #bio5
  elev_SHAP <- character() #elevation
  
  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
    MAX_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_5']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanPPT_maxT <- cbind(Species = sps_list[i],
                         tab_occ_vars[,c(1:4,6,9,11)],
                         data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                    Max_T_SHAP = MAX_temp_SHAP,
                                    elev_SHAP = elev_SHAP))
  
  #save results per species
  setwd(wd_res_species)
  write.csv(meanPPT_maxT, paste0(sps,'_maxT_.csv'), row.names = F)
  
  
  ## mean T to compare PPT variables importance
  
  ##############################################
  ############## mean T + min PPT ##############
  ##############################################
  
  shap_PPT_min <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_min_PPT]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_PPT_min,
    nrounds = 65L
  )
  
  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_min_PPT]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_temp_SHAP <- character()  #bio1
  MIN_prec_SHAP <- character()  #bio14
  elev_SHAP <- character() #elevation
  
  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
    MIN_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_14']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanT_minPPT <- cbind(Species = sps_list[i],
                        tab_occ_vars[,c(1:4,6,7,11)],
                        data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                   Min_PPT_SHAP = MIN_prec_SHAP,
                                   elev_SHAP = elev_SHAP))
  
  #save results per species
  setwd(wd_res_species)
  write.csv(meanT_minPPT, paste0(sps,'_minPPT_.csv'), row.names = F)
  
  
  ###############################################
  ############## mean T + mean PPT ##############
  ###############################################
  
  shap_PPT_mean <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_mean_PPT]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_PPT_mean,
    nrounds = 65L
  )
  
  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_mean_PPT]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_temp_SHAP <- character()  #bio1
  MEAN_prec_SHAP <- character()  #bio12
  elev_SHAP <- character() #elevation
  
  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
    MEAN_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanT_meanPPT <- cbind(Species = sps_list[i],
                         tab_occ_vars[,c(1:6,11)],
                         data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                    Mean_PPT_SHAP = MEAN_prec_SHAP,
                                    elev_SHAP = elev_SHAP))
  
  #save results per species
  setwd(wd_res_species)
  write.csv(meanT_meanPPT, paste0(sps,'_meanPPT_.csv'), row.names = F)
  
  
  ##############################################
  ############## mean T + max PPT ##############
  ##############################################
  
  shap_PPT_max <- xgb.DMatrix(data.matrix(tab_occ_vars[preds_max_PPT]), label = tab_occ_vars$Occurrence)
  
  fit <- xgb.train(
    params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
    data =   shap_PPT_max,
    nrounds = 65L
  )
  
  # We also pass feature data X with originally encoded values
  shp <- shapviz(fit, X_pred = data.matrix(tab_occ_vars[preds_max_PPT]), X = tab_occ_vars)
  
  # get the SHAP values for each variable in each prediction 
  MEAN_temp_SHAP <- character()  #bio1
  MAX_prec_SHAP <- character()  #bio13
  elev_SHAP <- character() #elevation
  
  for(j in 1:nrow(tab_occ_vars))
  {
    obj <- sv_waterfall(shp, row_id = j) +
      theme(axis.text = element_text(size = 11))
    
    MEAN_temp_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
    MAX_prec_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_13']
    elev_SHAP[j] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_elev']
    
    print(j)
  }
  
  #make a data.frame with results
  meanT_maxPPT <- cbind(Species = sps_list[i],
                        tab_occ_vars[,c(1:5,10,11)],
                        data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                   Max_PPT_SHAP = MAX_prec_SHAP,
                                   elev_SHAP))
  
  #save results per species
  setwd(wd_res_species)
  write.csv(meanT_maxPPT, paste0(sps,'_maxPPT_.csv'), row.names = F)
  
}

/Users/carloseduardoaribeiro/Documents/Files/Old files

#################


plot(sps_range2)


#make a spatial object
occ_sps_sp <- occ_sps
coordinates(occ_sps_sp) <- ~ decimalLongitude + decimalLatitude

plot(world, border = NA, col = 'gray80')
plot(occ_sps_sp, pch = 19, col = 'red', cex = 0.4, add = T)




plot(t, add = T)

class(t)
length(t)

#ger world map
world <- getMap()

class(kden.p_groups_trim)
head(kden.p_groups_trim)

class(kden.p_groups_trim[[1]])
length(kden.p_groups_trim[[1]])
head(kden.p_groups_trim[[1]])

class(kden.p_groups_trim[[2]])
length(kden.p_groups_trim[[2]])
head(kden.p_groups_trim[[2]])

class(kden.p_groups_trim[[3]])
length(kden.p_groups_trim[[3]])

class(kden.p_groups_trim[[3]][[i]])
length(kden.p_groups_trim[[3]][[i]])
head(kden.p_groups_trim[[3]][[i]])

class(kden.p_groups_trim[[3]][[i]][[1]])
length(kden.p_groups_trim[[3]][[i]][[1]])

class(kden.p_groups_trim[[3]][[i]][[2]])
length(kden.p_groups_trim[[3]][[i]][[2]])

class(kden.p_groups_trim[[3]][[i]][[1]][[1]])
nrow(kden.p_groups_trim[[3]][[i]][[1]][[1]])

class(kden.p_groups_trim[[3]][[i]][[2]][[1]])
nrow(kden.p_groups_trim[[3]][[i]][[2]][[1]])

t <- as.data.frame(kden.p_groups_trim[[3]][[i]][[1]][[1]])
t2 <- as.data.frame(kden.p_groups_trim[[3]][[i]][[2]][[1]])
t3 <- as.data.frame(kden.p_groups_trim[[3]][[i]][[3]][[1]])
t3 <- as.data.frame(kden.p_groups_trim[[3]][[i]][[93]][[1]])


plot(world)

head(t)
coordinates(t) <- ~ V1 + V2
t

plot(t, add = T)

head(t2)
coordinates(t2) <- ~ V1 + V2
t2

plot(world)
plot(t2, add = F)

head(t3)
coordinates(t3) <- ~ V1 + V2
t3

plot(world)
plot(t3, add = F)

setwd('/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP')
load('sf_world2.RData')


