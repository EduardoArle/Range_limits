#load packages
library(sf); library(raster); library(xgboost); library(PresenceAbsence)
library(dismo); library(data.table); library(plyr)

# library(data.table);library(rworldmap)
# ;library(sf);library(shapviz);library(kernelshap)
# ;library(caret);library(pROC)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_variables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_bias <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Biases'
wd_harmo <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review'
wd_res_species <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SHAP_results'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Thinned Occ'
wd_sps_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Species_occ'

#list species 
setwd(wd_tables)
sps_table <- read.csv('Selected_species.csv')
sps_list <- sps_table$species

#load lookup table for name harmonisation according to GBIF backbone
setwd(wd_harmo)
harmo <- read.csv('Harmonised_table.csv')

#load order thinned occurrences
setwd(wd_thinned_occ)
ord_occ <- lapply(list.files(), read.csv)
names(ord_occ) <- gsub('_thin.csv', '', list.files())

#load all six BioCLim variables being used
setwd(wd_variables)
AnnualMeanTemperature <- raster('wc2.1_2.5m_bio_1.tif')
MaxTemperatureOfWarmestMonth <- raster('wc2.1_2.5m_bio_5.tif')
MinTemperatureOfColdestMonth <- raster('wc2.1_2.5m_bio_6.tif')
AnnualPrecipitation <- raster('wc2.1_2.5m_bio_12.tif')
PrecipitationOfWettestMonth <- raster('wc2.1_2.5m_bio_13.tif')
PrecipitationOfDriestMonth <- raster('wc2.1_2.5m_bio_14.tif')

## List variables we will use for each predictions

#mean PPT to compare temperature variables
preds_min_temp <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_6')
preds_mean_temp  <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
preds_max_temp <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_5')

#mean T to compare precipitation variables
preds_min_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_14')
preds_mean_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_12')
preds_max_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_13')

# number of spatial folds for cross validation
k_folds <- 5           

#define parameters for xgboost
params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.1,
               max_depth = 3, subsample = 0.7, colsample_bytree = 0.7,
               min_child_weight = 0, gamma = 0)

#install modified randomPoints function
safeRandomPoints <- function(mask, n, prob = TRUE, min_n = 1,
                             verbose = FALSE){
  # mask: RasterLayer/SpatRaster
  # n: requested number of points
  # prob: pass-through to randomPoints()
  # min_n: smallest acceptable n before giving up
  # verbose: message when n is reduced
  
  # count sampleable cells (for prob=TRUE, zeros are effectively not sampleable)
  vals <- raster::getValues(mask)
  if (prob) {
    n_avail <- sum(!is.na(vals) & vals > 0)
  } else {
    n_avail <- sum(!is.na(vals))
  }
  
  if (n_avail < min_n) {
    if (verbose) message("safeRandomPoints: no available cells to sample.")
    return(NULL)
  }
  
  n_use <- min(n, n_avail)
  if (verbose && n_use < n) message("safeRandomPoints: reducing n from ", n, " to ", n_use)
  
  # attempt, and if something still goes wrong, reduce gradually until it works
  repeat {
    out <- try(
      dismo::randomPoints(mask = mask, n = n_use, prob = prob),
      silent = TRUE
    )
    if (!inherits(out, "try-error")) return(out)
    
    n_use <- floor(n_use * 0.9)  # reduce by 10% each attempt
    if (verbose) message("safeRandomPoints: retry with n = ", n_use)
    
    if (n_use < min_n) {
      if (verbose) message("safeRandomPoints: reached min_n without success.")
      return(NULL)
    }
  }
}

######## Run SHAP for all species ######

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  ## Harmonise species according to GBIF backbone
  
  #find GBIF name from crosswalk
  idx <- match(sps, harmo$iucn_name)
  gbif_name <- harmo$gbif_name[idx]
  
  #if no GBIF name, skip species
  if(is.na(gbif_name) || !nzchar(gbif_name)) next
  
  #get species order
  order <- harmo$order[idx]
  
  #select order thinned occurrences
  ord_occ2 <- ord_occ[[order]]
  
  #select only points representing the species
  sps_occ <- ord_occ2[which(ord_occ2$species == gbif_name),]
  
  ## Select only points whithin the species range
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #get area of main fragment of the range
  dom_frag_km2 <- as.numeric(strsplit(sps_table$range_frag_km2[i],
                                      ";")[[1]][sps_table$dom_frag[i]])
  
  #split into fragments (one POLYGON per fragment)
  frags <- st_cast(range, "POLYGON")
  
  #project to a local equal-area CRS so areas are meaningful
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  crs_laea <- paste0("+proj=laea +lat_0=", cen[2], " +lon_0=", cen[1],
                     " +datum=WGS84 +units=m +no_defs")
  frags_m <- st_transform(frags, crs_laea)
  
  #areas in km2
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #pick the fragment whose area is closest to dom_frag_km2
  frag_pos <- which.min(abs(a_km2 - dom_frag_km2))
  
  #keep only the dominant fragment, preserving the original attributes
  range_dom <- range
  range_dom$geometry <- st_geometry(frags[frag_pos, ])
  
  #project dominant range to metres (same CRS used for fragment areas)
  range_dom_m <- st_transform(range_dom, crs_laea)
  
  #create spatial object with sps presences
  sps_occ_sf <- st_as_sf(sps_occ, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(range_dom))
  
  #select only the points within the species range
  pts_range <- st_intersects(sps_occ_sf, range_dom, sparse = F)[, 1]
  sps_occ_range <- sps_occ[pts_range,] #table
  sps_occ_range_sf <- sps_occ_sf[pts_range,] #spatial object
  
  #check whether there are at least 10 records
  if(nrow(sps_occ_range_sf) < 11){
    next
  }
  
  ## Select pseudo-absences (same number of presences)
  
  #project presence points to metres
  sps_occ_range_sf_m <- st_transform(sps_occ_range_sf, crs_laea)
  
  #buffer in metres (50 km)
  sps_range_buf_m <- st_buffer(range_dom_m, 50000)
  
  #make 50km buffer around points to delimit area where I don't want PA
  small_buffer_m <- st_buffer(sps_occ_range_sf_m, 50000)
  
  #make a spatial polygon object with only one feature
  no_pa_area <- st_union(small_buffer_m)
  
  # this fixes possible 'duplicate vertex' errors
  no_pa_area <- st_make_valid(no_pa_area) 
  
  #make holes in the species range by the small buffer around points
  pa_area <- st_difference(sps_range_buf_m, no_pa_area)
  
  #load bias raster for this order
  setwd(wd_bias)
  bias_r <- raster(paste0(order, "_bias_0.5deg.tif"))
  
  #transfor pa_area to RGS84
  pa_area_ll <- st_transform(pa_area, st_crs(bias_r))
  
  #restrict bias raster to PA area
  bias_pa <- mask(bias_r, as(pa_area_ll, "Spatial"))
  bias_pa <- crop(bias_pa, as(pa_area_ll, "Spatial"))
  
  #crop the fine-resolution template to the same extent
  pred1_pa <- crop(AnnualMeanTemperature, extent(bias_pa))
  
  #resample bias layer to the resolution of variables
  bias_pa_fine <- resample(bias_pa, pred1_pa, method = "ngb")
  
  # number of pseudo-absences = presences inside range
  n_pa <- nrow(sps_occ_range)
  
  # sample pseudo-absences
  pa_xy <- tryCatch(safeRandomPoints(bias_pa_fine, n_pa, prob = T, min_n = 10,
                                     verbose = TRUE), error = function(e) NULL)
  
  if(is.null(pa_xy) || nrow(pa_xy) == 0){
    next
  }
  
  # convert to sf
  pa_sf <- st_as_sf(as.data.frame(pa_xy), coords = c(1, 2),
                    crs = st_crs(range))
  
  #get coords of pa
  pa_coords <- as.data.frame(st_coordinates(pa_sf))
  names(pa_coords) <- c('decimalLongitude', 'decimalLatitude')
  pa_coords$Occurrence <- 0 #include column informing occ status
  pa_coords$Type <- 'Pseudo-absence'
  
  #get coords of pr
  pr_coords <- as.data.frame(st_coordinates(sps_occ_range_sf))
  names(pr_coords) <- c('decimalLongitude', 'decimalLatitude')
  pr_coords$Occurrence <- 1 #include column informing occ status
  pr_coords$Type <- 'Presence'
  
  #combine pseudo absences and presences
  species <- rbind(pr_coords, pa_coords)
  
  #include species name
  species <- cbind(sps, species)
  
  #save table with species occ
  setwd(wd_sps_occ)
  write.csv(species, paste0(sps, '_occ.csv'), row.names = F)
  
  #create ID for points
  species$ID <- seq_len(nrow(species))
  
  #make spatial obj again
  species_sf <- st_as_sf(species, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = crs(range_dom))
  
  #select variables we will use
  preds <- stack(AnnualMeanTemperature, AnnualPrecipitation,
                 MinTemperatureOfColdestMonth, PrecipitationOfDriestMonth,
                 MaxTemperatureOfWarmestMonth, PrecipitationOfWettestMonth)
  
  #extract values from each location from all variables
  vals_pts <- extract(preds, species_sf)
  
  #make a table creating an ID for each point
  tab_occ_vars <- cbind(species, vals_pts)
  
  
  ## Build the spatial blocks for species
  
  #local projected CRS (metres) centred on this species PA area
  cen <- st_coordinates(st_centroid(st_union(sps_range_buf)))[1, ]
  crs_laea <- paste0("+proj=laea +lat_0=", cen[2], " +lon_0=", cen[1],
                     " +datum=WGS84 +units=m +no_defs")
  
  #project domain + points
  dom_m <- st_transform(sps_range_buf_m, crs_laea)
  pts_m <- st_transform(species_sf, crs_laea)
  
  #make a single polygon for the domain (avoid edge issues)
  dom_u <- st_union(dom_m)
  
  #define the block size to split the data geographically (in metres - 200 km)
  block_size <- 200000 
  
  #grid blocks (200 km squares) over the full domain extent
  grid <- st_make_grid(dom_u, cellsize = block_size, square = TRUE)
  grid <- st_sf(block_id = seq_along(grid), geometry = grid)
  
  #keep only blocks whose centroid is within the domain
  grid <- grid[st_intersects(grid, dom_u, sparse = FALSE)[, 1], ]
  
  #reduce block_size until grid has at least 5 cells
  while(nrow(grid) < 5){
    
    #divide block_size per 2
    block_size <- block_size/2
    
    #grid blocks (200 km squares) over the full domain extent
    grid <- st_make_grid(dom_u, cellsize = block_size, square = TRUE)
    grid <- st_sf(block_id = seq_along(grid), geometry = grid)
    
    #keep only blocks whose centroid is within the domain
    grid <- grid[st_within(st_centroid(grid), dom_u, sparse = FALSE)[, 1], ]
  }
  
  #assign each point (presence + PA) to a block
  pts_blk <- st_join(pts_m, grid["block_id"], left = T)
  
  #attach block_id to table
  blk_tbl <- st_drop_geometry(pts_blk)[, c("ID", "block_id")]
  pts_tbl <- merge(tab_occ_vars, blk_tbl, by = "ID", all.x = T)
  
  #remove points not assigned to any spatial block (edge cases)
  pts_tbl <- pts_tbl[!is.na(pts_tbl$block_id), ]
  
  #set seed for reproducibility
  set.seed(1)
  
  #get unique blocks
  blk_ids <- unique(pts_tbl$block_id)
  
  #assign each block to a fold
  blk_fold <- sample(rep(1:k_folds, length.out = length(blk_ids)))
  
  #attach fold to each point (by block)
  pts_tbl$fold <- blk_fold[match(pts_tbl$block_id, blk_ids)]
  
  #include block_size in the table
  pts_tbl$block_size <- block_size
  
  #make 3 replicates of the process
  for(j in 1:3)
  {
    #set seed to make block-to-fold assignment reproducible across repetitions
    set.seed(1000 + j)
    
    #get unique spatial blocks containing data
    u_blocks <- unique(pts_tbl$block_id)
    
    #randomly assign blocks to k folds
    fold_map <- data.frame(block_id = u_blocks,
                           fold = sample(rep(1:k_folds,
                                             length.out = length(u_blocks))))
    
    #attach fold ID to each point via its block ID
    pts_tbl$fold <- fold_map$fold[match(pts_tbl$block_id, fold_map$block_id)]
    
    #create directory to save each repetition for the species
    dir_species <- paste0(wd_res_species, '/', sps_list[i])
    dir.create(dir_species)
    
    #perform 5 fold cross validation for each model
    for(k in 1:k_folds)
    {
      #subset training and testing data based on spatial folds
      train_dat <- pts_tbl[pts_tbl$fold != k, ]
      test_dat <- pts_tbl[pts_tbl$fold == k, ]
      
      #skip bad folds
      if(nrow(test_dat) == 0 ||
         length(unique(test_dat$Occurrence)) < 2 ||
         length(unique(train_dat$Occurrence)) < 2) {
        next
      }
      
      ## Fit XGBoost models
      
      ## Mean PPT to compare T variables importance
      
      
      ##############################################
      ############## mean PPT + min T ##############
      ##############################################
      
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_min_temp]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest  <- xgb.DMatrix(as.matrix(test_dat[, preds_min_temp]),
                            label = test_dat$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                        observed = as.integer(train_dat$Occurrence),
                        predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                       observed = as.integer(test_dat$Occurrence),
                       predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test,  threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test,  st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test  - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_min_temp]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
                         train_dat[,preds_min_temp])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanPPT_minT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_6")],
        data.frame(
          Mean_PPT_SHAP = shap_df$wc2.1_2.5m_bio_12,
          Min_T_SHAP = shap_df$wc2.1_2.5m_bio_6,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanPPT_minT,
                paste0(sps, '_minT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
      
      
      ###############################################
      ############## mean PPT + mean T ##############
      ###############################################
      
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_mean_temp]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest <- xgb.DMatrix(as.matrix(test_dat[, preds_mean_temp]),
                            label = test_dat$Occurrence)
      
     
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                              observed = as.integer(train_dat$Occurrence),
                              predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                                  observed = as.integer(test_dat$Occurrence),
                                  predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test,  threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test,  st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test  - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_mean_temp]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
        train_dat[,preds_mean_temp])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanPPT_meanT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_1")],
        data.frame(
          Mean_PPT_SHAP = shap_df$wc2.1_2.5m_bio_12,
          Mean_T_SHAP = shap_df$wc2.1_2.5m_bio_1,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanPPT_meanT,
                paste0(sps, '_meanT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
      

      ##############################################
      ############## mean PPT + max T ##############
      ##############################################
    
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_max_temp]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest  <- xgb.DMatrix(as.matrix(test_dat[, preds_max_temp]),
                            label = test_dat$Occurrence)
      
      
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test  <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                                   observed = as.integer(train_dat$Occurrence),
                                   predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                                  observed = as.integer(test_dat$Occurrence),
                                  predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test,  threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test,  st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test  - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_max_temp]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
        train_dat[,preds_max_temp])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanPPT_maxT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_5")],
        data.frame(
          Mean_PPT_SHAP = shap_df$wc2.1_2.5m_bio_12,
          Max_T_SHAP = shap_df$wc2.1_2.5m_bio_5,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanPPT_maxT,
                paste0(sps, '_maxT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
      
      
      ## Mean PPT to compare T variables importance
      
      
      ##############################################
      ############## mean T + min PPT ##############
      ##############################################
      
      
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_min_PPT]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest  <- xgb.DMatrix(as.matrix(test_dat[, preds_min_PPT]),
                            label = test_dat$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                                   observed = as.integer(train_dat$Occurrence),
                                   predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                                  observed = as.integer(test_dat$Occurrence),
                                  predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test,  threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test,  st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test  - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_min_PPT]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
        train_dat[,preds_min_PPT])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanT_minPPT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_14")],
        data.frame(
          Mean_T_SHAP = shap_df$wc2.1_2.5m_bio_1,
          Min_PPT_SHAP = shap_df$wc2.1_2.5m_bio_14,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanT_minPPT,
                paste0(sps, '_minPPT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
      
      
      ###############################################
      ############## mean T + mean PPT ##############
      ###############################################
      
      
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_mean_PPT]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest  <- xgb.DMatrix(as.matrix(test_dat[, preds_mean_PPT]),
                            label = test_dat$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                                   observed = as.integer(train_dat$Occurrence),
                                   predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                                  observed = as.integer(test_dat$Occurrence),
                                  predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test,  st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_mean_PPT]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
        train_dat[,preds_mean_PPT])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanT_meanPPT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_12")],
        data.frame(
          Mean_T_SHAP = shap_df$wc2.1_2.5m_bio_1,
          Mean_PPT_SHAP = shap_df$wc2.1_2.5m_bio_12,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanT_meanPPT,
                paste0(sps, '_meanPPT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
      
      
      ##############################################
      ############## mean T + max PPT ##############
      ##############################################
      
      
      #prepare train data
      dtrain <- xgb.DMatrix(as.matrix(train_dat[, preds_max_PPT]),
                            label = train_dat$Occurrence)
      
      #prepare test data
      dtest  <- xgb.DMatrix(as.matrix(test_dat[, preds_max_PPT]),
                            label = test_dat$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(params = params, data = dtrain, nrounds = 100,
                       watchlist = list(train = dtrain, test = dtest),
                       verbose = 0)
      
      
      ## Get evaluation metrics
      
      #predictions
      pred_train <- predict(fit, dtrain)
      pred_test <- predict(fit, dtest)
      
      #AUC (from xgboost)
      AUC_train <- tail(fit$evaluation_log$train_auc, 1)
      AUC_test <- tail(fit$evaluation_log$test_auc, 1)
      
      #tables with observed and predicted values
      obs_pred_train <- data.frame(id = train_dat$ID,
                                   observed = as.integer(train_dat$Occurrence),
                                   predicted = pred_train)
      
      obs_pred_test <- data.frame(id = test_dat$ID,
                                  observed = as.integer(test_dat$Occurrence),
                                  predicted = pred_test)
      
      #threshold that maximises sensitivity + specificity
      th <- optimal.thresholds(obs_pred_train,
                               opt.methods = "MaxSens+Spec")$predicted
      
      #confusion matrices
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = FALSE)
      sens_test <- sensitivity(confusion_test, st.dev = FALSE)
      
      #specificity
      spec_train <- specificity(confusion_train, st.dev = FALSE)
      spec_test <- specificity(confusion_test,  st.dev = FALSE)
      
      #TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test  + spec_test - 1
      
      #get SHAP values for training data
      shap_train <- predict(fit, as.matrix(train_dat[, preds_max_PPT]),
                            predcontrib = TRUE)
      
      #drop BIAS column (keep only predictors)
      shap_train <- shap_train[,colnames(as.matrix(
        train_dat[,preds_max_PPT])), drop = FALSE]
      
      #convert SHAP matrix to data.frame
      shap_df <- as.data.frame(shap_train)
      
      #build results table (one row per training point)
      meanT_maxPPT <- cbind(
        Species = sps,
        train_dat[, c("ID",
                      "decimalLongitude", "decimalLatitude",
                      "Occurrence", "Type",
                      "wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_13")],
        data.frame(
          Mean_T_SHAP = shap_df$wc2.1_2.5m_bio_1,
          Max_PPT_SHAP = shap_df$wc2.1_2.5m_bio_13,
          iteration = j,
          fold = k,
          block_size = block_size,
          AUC_test = AUC_test,
          TSS_test = TSS_test,
          sens_test = sens_test,
          spec_test = spec_test,
          AUC_train = AUC_train,
          TSS_train = TSS_train,
          sens_train = sens_train,
          spec_train = spec_train))
      
      #save results
      setwd(dir_species)
      write.csv(meanT_maxPPT,
                paste0(sps, '_maxPPT_run_', j, '_rep_', k, '.csv'),
                row.names = F)
    }
  }
}




##### Calculate average SHAP for each point and the corresponding metrics ######

#list with all species that produced models
setwd(wd_res_species)
sps_list <- sub("^\\./", "", list.dirs(recursive = FALSE))

#create an empty objects to populate with the total number of models
tot_models <- numeric()

#create an empty objects to populate with info about the percentage of useful models
perc_used_minT <- numeric()
perc_used_meanT <- numeric()
perc_used_maxT <- numeric()

perc_used_minPPT <- numeric()
perc_used_meanPPT <- numeric()
perc_used_maxPPT <- numeric()

#evaluation metrics for the species of the used models (test metrics)

#### minT ####
AUC_all_minT <- numeric()
TSS_all_minT <- numeric()
sens_all_minT <- numeric()
spec_all_minT <- numeric()

AUC_sel_minT <- numeric()
TSS_sel_minT <- numeric()
sens_sel_minT <- numeric()
spec_sel_minT <- numeric()

SD_AUC_all_minT <- numeric()
SD_TSS_all_minT <- numeric()
SD_sens_all_minT <- numeric()
SD_spec_all_minT <- numeric()

SD_AUC_sel_minT <- numeric()
SD_TSS_sel_minT <- numeric()
SD_sens_sel_minT <- numeric()
SD_spec_sel_minT <- numeric()

#### meanT ####
AUC_all_meanT <- numeric()
TSS_all_meanT <- numeric()
sens_all_meanT <- numeric()
spec_all_meanT <- numeric()

AUC_sel_meanT <- numeric()
TSS_sel_meanT <- numeric()
sens_sel_meanT <- numeric()
spec_sel_meanT <- numeric()

SD_AUC_all_meanT <- numeric()
SD_TSS_all_meanT <- numeric()
SD_sens_all_meanT <- numeric()
SD_spec_all_meanT <- numeric()

SD_AUC_sel_meanT <- numeric()
SD_TSS_sel_meanT <- numeric()
SD_sens_sel_meanT <- numeric()
SD_spec_sel_meanT <- numeric()

#### maxT ####
AUC_all_maxT <- numeric()
TSS_all_maxT <- numeric()
sens_all_maxT <- numeric()
spec_all_maxT <- numeric()

AUC_sel_maxT <- numeric()
TSS_sel_maxT <- numeric()
sens_sel_maxT <- numeric()
spec_sel_maxT <- numeric()

SD_AUC_all_maxT <- numeric()
SD_TSS_all_maxT <- numeric()
SD_sens_all_maxT <- numeric()
SD_spec_all_maxT <- numeric()

SD_AUC_sel_maxT <- numeric()
SD_TSS_sel_maxT <- numeric()
SD_sens_sel_maxT <- numeric()
SD_spec_sel_maxT <- numeric()

#### minPPT ####
AUC_all_minPPT <- numeric()
TSS_all_minPPT <- numeric()
sens_all_minPPT <- numeric()
spec_all_minPPT <- numeric()

AUC_sel_minPPT <- numeric()
TSS_sel_minPPT <- numeric()
sens_sel_minPPT <- numeric()
spec_sel_minPPT <- numeric()

SD_AUC_all_minPPT <- numeric()
SD_TSS_all_minPPT <- numeric()
SD_sens_all_minPPT <- numeric()
SD_spec_all_minPPT <- numeric()

SD_AUC_sel_minPPT <- numeric()
SD_TSS_sel_minPPT <- numeric()
SD_sens_sel_minPPT <- numeric()
SD_spec_sel_minPPT <- numeric()

#### meanPPT ####
AUC_all_meanPPT <- numeric()
TSS_all_meanPPT <- numeric()
sens_all_meanPPT <- numeric()
spec_all_meanPPT <- numeric()

AUC_sel_meanPPT <- numeric()
TSS_sel_meanPPT <- numeric()
sens_sel_meanPPT <- numeric()
spec_sel_meanPPT <- numeric()

SD_AUC_all_meanPPT <- numeric()
SD_TSS_all_meanPPT <- numeric()
SD_sens_all_meanPPT <- numeric()
SD_spec_all_meanPPT <- numeric()

SD_AUC_sel_meanPPT <- numeric()
SD_TSS_sel_meanPPT <- numeric()
SD_sens_sel_meanPPT <- numeric()
SD_spec_sel_meanPPT <- numeric()

#### maxPPT ####
AUC_all_maxPPT <- numeric()
TSS_all_maxPPT <- numeric()
sens_all_maxPPT <- numeric()
spec_all_maxPPT <- numeric()

AUC_sel_maxPPT <- numeric()
TSS_sel_maxPPT <- numeric()
sens_sel_maxPPT <- numeric()
spec_sel_maxPPT <- numeric()

SD_AUC_all_maxPPT <- numeric()
SD_TSS_all_maxPPT <- numeric()
SD_sens_all_maxPPT <- numeric()
SD_spec_all_maxPPT <- numeric()

SD_AUC_sel_maxPPT <- numeric()
SD_TSS_sel_maxPPT <- numeric()
SD_sens_sel_maxPPT <- numeric()
SD_spec_sel_maxPPT <- numeric()


for(i in 983:length(sps_list))
{
  #select species 
  sps <- sps_list[i]
  
  #check if the species has any models
  setwd(wd_res_species)
  sps_list_models <- list.dirs(full.names = F, recursive = F)
  
  if(sps %in% sps_list_models){
    
    
    
    ##############################################
    ############## mean PPT + min T ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the min T analysis
    minT <- lapply(list.files(pattern = 'minT'), read.csv)
    
    #input the number of folds (therefore models) that could be run. This number is the same for all variable combination models because it depends only on the number of records
    tot_models[i] <- length(minT)
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_minT[i] <- sum(sapply(minT, function(x) unique(x$TSS_test)) >= 0.4 &
                            sapply(minT, function(x) unique(x$AUC_test)) >= 0.7) /
      length(minT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    minT_bind_all <- rbindlist(minT)
    
    #get metrics for all models
    AUC_all_minT[i] <- mean(minT_bind_all$AUC_test)
    TSS_all_minT[i] <- mean(minT_bind_all$TSS_test)
    sens_all_minT[i] <- mean(minT_bind_all$sens_test)
    spec_all_minT[i] <- mean(minT_bind_all$spec_test)
    
    SD_AUC_all_minT[i] <- sd(minT_bind_all$AUC_test)
    SD_TSS_all_minT[i] <- sd(minT_bind_all$TSS_test)
    SD_sens_all_minT[i] <- sd(minT_bind_all$sens_test)
    SD_spec_all_minT[i] <- sd(minT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    minT_bind <- minT_bind_all[minT_bind_all$AUC_test >= 0.7 &
                                 minT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_minT[i] <- mean(minT_bind$AUC_test)
    TSS_sel_minT[i] <- mean(minT_bind$TSS_test)
    sens_sel_minT[i] <- mean(minT_bind$sens_test)
    spec_sel_minT[i] <- mean(minT_bind$spec_test)
    
    SD_AUC_sel_minT[i] <- sd(minT_bind$AUC_test)
    SD_TSS_sel_minT[i] <- sd(minT_bind$TSS_test)
    SD_sens_sel_minT[i] <- sd(minT_bind$sens_test)
    SD_spec_sel_minT[i] <- sd(minT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(minT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    minT_average <- ddply(minT_bind, .(ID), summarise,
                          Species = unique(Species),
                          decimalLongitude = unique(decimalLongitude),
                          decimalLatitude = unique(decimalLatitude),
                          Occurrence = unique(Occurrence),
                          Type = unique(Type),
                          wc2.1_2.5m_bio_12 = unique(wc2.1_2.5m_bio_12),
                          wc2.1_2.5m_bio_6 = unique(wc2.1_2.5m_bio_6),
                          avg_Mean_PPT_SHAP = mean(Mean_PPT_SHAP),
                          SD_Mean_PPT_SHAP = sd(Mean_PPT_SHAP),
                          avg_Min_T_SHAP = mean(Min_T_SHAP),
                          SD_Min_T_SHAP = sd(Min_T_SHAP),
                          avg_AUC_test = mean(AUC_test),
                          SD_AUC_test = sd(AUC_test),
                          avg_TSS_test = mean(TSS_test),
                          SD_TSS_test = sd(TSS_test),
                          avg_sens_test = mean(sens_test),
                          SD_sens_test = sd(sens_test),
                          avg_spec_test = mean(spec_test),
                          SD_spec_test = sd(spec_test),
                          avg_AUC_train = mean(AUC_train),
                          SD_AUC_train = sd(AUC_train),
                          avg_TSS_train = mean(TSS_train),
                          SD_TSS_train = sd(TSS_train),
                          avg_sens_train = mean(sens_train),
                          SD_sens_train = sd(sens_train),
                          avg_spec_train = mean(spec_train),
                          SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    minT_average <- try(merge(minT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(minT_average, paste0(sps, '_minT', '.csv'), row.names = F)
    
    
    
    ##############################################
    ############## mean PPT + mean T ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the mean T analysis
    meanT <- lapply(list.files(pattern = 'meanT'), read.csv)
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_meanT[i] <- sum(sapply(meanT,
                                     function(x) unique(x$TSS_test)) >= 0.4 &
                              sapply(meanT,
                                     function(x) unique(x$AUC_test)) >= 0.7) /
    length(meanT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    meanT_bind_all <- rbindlist(meanT)
    
    #get metrics for all models
    AUC_all_meanT[i] <- mean(meanT_bind_all$AUC_test)
    TSS_all_meanT[i] <- mean(meanT_bind_all$TSS_test)
    sens_all_meanT[i] <- mean(meanT_bind_all$sens_test)
    spec_all_meanT[i] <- mean(meanT_bind_all$spec_test)
    
    SD_AUC_all_meanT[i] <- sd(meanT_bind_all$AUC_test)
    SD_TSS_all_meanT[i] <- sd(meanT_bind_all$TSS_test)
    SD_sens_all_meanT[i] <- sd(meanT_bind_all$sens_test)
    SD_spec_all_meanT[i] <- sd(meanT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    meanT_bind <- meanT_bind_all[meanT_bind_all$AUC_test >= 0.7 &
                                 meanT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_meanT[i] <- mean(meanT_bind$AUC_test)
    TSS_sel_meanT[i] <- mean(meanT_bind$TSS_test)
    sens_sel_meanT[i] <- mean(meanT_bind$sens_test)
    spec_sel_meanT[i] <- mean(meanT_bind$spec_test)
    
    SD_AUC_sel_meanT[i] <- sd(meanT_bind$AUC_test)
    SD_TSS_sel_meanT[i] <- sd(meanT_bind$TSS_test)
    SD_sens_sel_meanT[i] <- sd(meanT_bind$sens_test)
    SD_spec_sel_meanT[i] <- sd(meanT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(meanT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    meanT_average <- ddply(meanT_bind, .(ID), summarise,
                          Species = unique(Species),
                          decimalLongitude = unique(decimalLongitude),
                          decimalLatitude = unique(decimalLatitude),
                          Occurrence = unique(Occurrence),
                          Type = unique(Type),
                          wc2.1_2.5m_bio_12 = unique(wc2.1_2.5m_bio_12),
                          wc2.1_2.5m_bio_1 = unique(wc2.1_2.5m_bio_1),
                          avg_Mean_PPT_SHAP = mean(Mean_PPT_SHAP),
                          SD_Mean_PPT_SHAP = sd(Mean_PPT_SHAP),
                          avg_Mean_T_SHAP = mean(Mean_T_SHAP),
                          SD_Mean_T_SHAP = sd(Mean_T_SHAP),
                          avg_AUC_test = mean(AUC_test),
                          SD_AUC_test = sd(AUC_test),
                          avg_TSS_test = mean(TSS_test),
                          SD_TSS_test = sd(TSS_test),
                          avg_sens_test = mean(sens_test),
                          SD_sens_test = sd(sens_test),
                          avg_spec_test = mean(spec_test),
                          SD_spec_test = sd(spec_test),
                          avg_AUC_train = mean(AUC_train),
                          SD_AUC_train = sd(AUC_train),
                          avg_TSS_train = mean(TSS_train),
                          SD_TSS_train = sd(TSS_train),
                          avg_sens_train = mean(sens_train),
                          SD_sens_train = sd(sens_train),
                          avg_spec_train = mean(spec_train),
                          SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    meanT_average <- try(merge(meanT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(meanT_average, paste0(sps, '_meanT', '.csv'), row.names = F)
    
    
    
    ##############################################
    ############## mean PPT + max T ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the max T analysis
    maxT <- lapply(list.files(pattern = 'maxT'), read.csv)
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_maxT[i] <- sum(sapply(maxT,
                                    function(x) unique(x$TSS_test)) >= 0.4 &
                                sapply(maxT,
                                       function(x) unique(x$AUC_test)) >= 0.7) /
      length(maxT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    maxT_bind_all <- rbindlist(maxT)
    
    #get metrics for all models
    AUC_all_maxT[i] <- mean(maxT_bind_all$AUC_test)
    TSS_all_maxT[i] <- mean(maxT_bind_all$TSS_test)
    sens_all_maxT[i] <- mean(maxT_bind_all$sens_test)
    spec_all_maxT[i] <- mean(maxT_bind_all$spec_test)
    
    SD_AUC_all_maxT[i] <- sd(maxT_bind_all$AUC_test)
    SD_TSS_all_maxT[i] <- sd(maxT_bind_all$TSS_test)
    SD_sens_all_maxT[i] <- sd(maxT_bind_all$sens_test)
    SD_spec_all_maxT[i] <- sd(maxT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    maxT_bind <- maxT_bind_all[maxT_bind_all$AUC_test >= 0.7 &
                                   maxT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_maxT[i] <- mean(maxT_bind$AUC_test)
    TSS_sel_maxT[i] <- mean(maxT_bind$TSS_test)
    sens_sel_maxT[i] <- mean(maxT_bind$sens_test)
    spec_sel_maxT[i] <- mean(maxT_bind$spec_test)
    
    SD_AUC_sel_maxT[i] <- sd(maxT_bind$AUC_test)
    SD_TSS_sel_maxT[i] <- sd(maxT_bind$TSS_test)
    SD_sens_sel_maxT[i] <- sd(maxT_bind$sens_test)
    SD_spec_sel_maxT[i] <- sd(maxT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(maxT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    maxT_average <- ddply(maxT_bind, .(ID), summarise,
                           Species = unique(Species),
                           decimalLongitude = unique(decimalLongitude),
                           decimalLatitude = unique(decimalLatitude),
                           Occurrence = unique(Occurrence),
                           Type = unique(Type),
                           wc2.1_2.5m_bio_12 = unique(wc2.1_2.5m_bio_12),
                           wc2.1_2.5m_bio_5 = unique(wc2.1_2.5m_bio_5),
                           avg_Mean_PPT_SHAP = mean(Mean_PPT_SHAP),
                           SD_Mean_PPT_SHAP = sd(Mean_PPT_SHAP),
                           avg_Max_T_SHAP = mean(Max_T_SHAP),
                           SD_Max_T_SHAP = sd(Max_T_SHAP),
                           avg_AUC_test = mean(AUC_test),
                           SD_AUC_test = sd(AUC_test),
                           avg_TSS_test = mean(TSS_test),
                           SD_TSS_test = sd(TSS_test),
                           avg_sens_test = mean(sens_test),
                           SD_sens_test = sd(sens_test),
                           avg_spec_test = mean(spec_test),
                           SD_spec_test = sd(spec_test),
                           avg_AUC_train = mean(AUC_train),
                           SD_AUC_train = sd(AUC_train),
                           avg_TSS_train = mean(TSS_train),
                           SD_TSS_train = sd(TSS_train),
                           avg_sens_train = mean(sens_train),
                           SD_sens_train = sd(sens_train),
                           avg_spec_train = mean(spec_train),
                           SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    maxT_average <- try(merge(maxT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(maxT_average, paste0(sps, '_maxT', '.csv'), row.names = F)
    
    
    
    ##############################################
    ############## mean T + min PPT ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the min T analysis
    minPPT <- lapply(list.files(pattern = 'minPPT'), read.csv)
    
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_minPPT[i] <- sum(sapply(minPPT,
                                      function(x) unique(x$TSS_test)) >= 0.4 &
                               sapply(minPPT,
                                      function(x) unique(x$AUC_test)) >= 0.7) /
      length(minPPT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    minPPT_bind_all <- rbindlist(minPPT)
    
    #get metrics for all models
    AUC_all_minPPT[i] <- mean(minPPT_bind_all$AUC_test)
    TSS_all_minPPT[i] <- mean(minPPT_bind_all$TSS_test)
    sens_all_minPPT[i] <- mean(minPPT_bind_all$sens_test)
    spec_all_minPPT[i] <- mean(minPPT_bind_all$spec_test)
    
    SD_AUC_all_minPPT[i] <- sd(minPPT_bind_all$AUC_test)
    SD_TSS_all_minPPT[i] <- sd(minPPT_bind_all$TSS_test)
    SD_sens_all_minPPT[i] <- sd(minPPT_bind_all$sens_test)
    SD_spec_all_minPPT[i] <- sd(minPPT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    minPPT_bind <- minPPT_bind_all[minPPT_bind_all$AUC_test >= 0.7 &
                                 minPPT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_minPPT[i] <- mean(minPPT_bind$AUC_test)
    TSS_sel_minPPT[i] <- mean(minPPT_bind$TSS_test)
    sens_sel_minPPT[i] <- mean(minPPT_bind$sens_test)
    spec_sel_minPPT[i] <- mean(minPPT_bind$spec_test)
    
    SD_AUC_sel_minPPT[i] <- sd(minPPT_bind$AUC_test)
    SD_TSS_sel_minPPT[i] <- sd(minPPT_bind$TSS_test)
    SD_sens_sel_minPPT[i] <- sd(minPPT_bind$sens_test)
    SD_spec_sel_minPPT[i] <- sd(minPPT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(minPPT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    minPPT_average <- ddply(minPPT_bind, .(ID), summarise,
                          Species = unique(Species),
                          decimalLongitude = unique(decimalLongitude),
                          decimalLatitude = unique(decimalLatitude),
                          Occurrence = unique(Occurrence),
                          Type = unique(Type),
                          wc2.1_2.5m_bio_1 = unique(wc2.1_2.5m_bio_1),
                          wc2.1_2.5m_bio_14 = unique(wc2.1_2.5m_bio_14),
                          avg_Mean_T_SHAP = mean(Mean_T_SHAP),
                          SD_Mean_T_SHAP = sd(Mean_T_SHAP),
                          avg_Min_PPT_SHAP = mean(Min_PPT_SHAP),
                          SD_Min_PPT_SHAP = sd(Min_PPT_SHAP),
                          avg_AUC_test = mean(AUC_test),
                          SD_AUC_test = sd(AUC_test),
                          avg_TSS_test = mean(TSS_test),
                          SD_TSS_test = sd(TSS_test),
                          avg_sens_test = mean(sens_test),
                          SD_sens_test = sd(sens_test),
                          avg_spec_test = mean(spec_test),
                          SD_spec_test = sd(spec_test),
                          avg_AUC_train = mean(AUC_train),
                          SD_AUC_train = sd(AUC_train),
                          avg_TSS_train = mean(TSS_train),
                          SD_TSS_train = sd(TSS_train),
                          avg_sens_train = mean(sens_train),
                          SD_sens_train = sd(sens_train),
                          avg_spec_train = mean(spec_train),
                          SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    minPPT_average <- try(merge(minPPT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(minPPT_average, paste0(sps, '_minPPT', '.csv'), row.names = F)
    
    
    
    ##############################################
    ############## mean T + mean PPT ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the min T analysis
    meanPPT <- lapply(list.files(pattern = 'meanPPT'), read.csv)
    
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_meanPPT[i] <- sum(sapply(meanPPT,
                                       function(x) unique(x$TSS_test)) >= 0.4 &
                                 sapply(meanPPT,
                                        function(x) unique(x$AUC_test)) >= 0.7) /
      length(meanPPT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    meanPPT_bind_all <- rbindlist(meanPPT)
    
    #get metrics for all models
    AUC_all_meanPPT[i] <- mean(meanPPT_bind_all$AUC_test)
    TSS_all_meanPPT[i] <- mean(meanPPT_bind_all$TSS_test)
    sens_all_meanPPT[i] <- mean(meanPPT_bind_all$sens_test)
    spec_all_meanPPT[i] <- mean(meanPPT_bind_all$spec_test)
    
    SD_AUC_all_meanPPT[i] <- sd(meanPPT_bind_all$AUC_test)
    SD_TSS_all_meanPPT[i] <- sd(meanPPT_bind_all$TSS_test)
    SD_sens_all_meanPPT[i] <- sd(meanPPT_bind_all$sens_test)
    SD_spec_all_meanPPT[i] <- sd(meanPPT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    meanPPT_bind <- meanPPT_bind_all[meanPPT_bind_all$AUC_test >= 0.7 &
                                     meanPPT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_meanPPT[i] <- mean(meanPPT_bind$AUC_test)
    TSS_sel_meanPPT[i] <- mean(meanPPT_bind$TSS_test)
    sens_sel_meanPPT[i] <- mean(meanPPT_bind$sens_test)
    spec_sel_meanPPT[i] <- mean(meanPPT_bind$spec_test)
    
    SD_AUC_sel_meanPPT[i] <- sd(meanPPT_bind$AUC_test)
    SD_TSS_sel_meanPPT[i] <- sd(meanPPT_bind$TSS_test)
    SD_sens_sel_meanPPT[i] <- sd(meanPPT_bind$sens_test)
    SD_spec_sel_meanPPT[i] <- sd(meanPPT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(meanPPT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    meanPPT_average <- ddply(meanPPT_bind, .(ID), summarise,
                            Species = unique(Species),
                            decimalLongitude = unique(decimalLongitude),
                            decimalLatitude = unique(decimalLatitude),
                            Occurrence = unique(Occurrence),
                            Type = unique(Type),
                            wc2.1_2.5m_bio_1 = unique(wc2.1_2.5m_bio_1),
                            wc2.1_2.5m_bio_12 = unique(wc2.1_2.5m_bio_12),
                            avg_Mean_T_SHAP = mean(Mean_T_SHAP),
                            SD_Mean_T_SHAP = sd(Mean_T_SHAP),
                            avg_Mean_PPT_SHAP = mean(Mean_PPT_SHAP),
                            SD_Mean_PPT_SHAP = sd(Mean_PPT_SHAP),
                            avg_AUC_test = mean(AUC_test),
                            SD_AUC_test = sd(AUC_test),
                            avg_TSS_test = mean(TSS_test),
                            SD_TSS_test = sd(TSS_test),
                            avg_sens_test = mean(sens_test),
                            SD_sens_test = sd(sens_test),
                            avg_spec_test = mean(spec_test),
                            SD_spec_test = sd(spec_test),
                            avg_AUC_train = mean(AUC_train),
                            SD_AUC_train = sd(AUC_train),
                            avg_TSS_train = mean(TSS_train),
                            SD_TSS_train = sd(TSS_train),
                            avg_sens_train = mean(sens_train),
                            SD_sens_train = sd(sens_train),
                            avg_spec_train = mean(spec_train),
                            SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    meanPPT_average <- try(merge(meanPPT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(meanPPT_average, paste0(sps, '_meanPPT', '.csv'), row.names = F)
    
    
    
    ##############################################
    ############## mean T + max PPT ##############
    ##############################################
    
    
    
    #setwd to the folder with all results for the species
    setwd(paste0(wd_res_species, '/', sps))
    
    #select all files corresponding to the min T analysis
    maxPPT <- lapply(list.files(pattern = 'maxPPT'), read.csv)
    
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_maxPPT[i] <- sum(sapply(maxPPT,
                                      function(x) unique(x$TSS_test)) >= 0.4 &
                                  sapply(maxPPT,
                                         function(x) unique(x$AUC_test)) >= 0.7) /
      length(maxPPT) * 100
    
    #rbind all lists to calculate the average and SD of shap values per point used in training
    maxPPT_bind_all <- rbindlist(maxPPT)
    
    #get metrics for all models
    AUC_all_maxPPT[i] <- mean(maxPPT_bind_all$AUC_test)
    TSS_all_maxPPT[i] <- mean(maxPPT_bind_all$TSS_test)
    sens_all_maxPPT[i] <- mean(maxPPT_bind_all$sens_test)
    spec_all_maxPPT[i] <- mean(maxPPT_bind_all$spec_test)
    
    SD_AUC_all_maxPPT[i] <- sd(maxPPT_bind_all$AUC_test)
    SD_TSS_all_maxPPT[i] <- sd(maxPPT_bind_all$TSS_test)
    SD_sens_all_maxPPT[i] <- sd(maxPPT_bind_all$sens_test)
    SD_spec_all_maxPPT[i] <- sd(maxPPT_bind_all$spec_test)
    
    #select only models with satisfactory AUC and TSS
    maxPPT_bind <- maxPPT_bind_all[maxPPT_bind_all$AUC_test >= 0.7 &
                                       maxPPT_bind_all$TSS_test >= 0.4]
    
    #get metrics for selected models
    AUC_sel_maxPPT[i] <- mean(maxPPT_bind$AUC_test)
    TSS_sel_maxPPT[i] <- mean(maxPPT_bind$TSS_test)
    sens_sel_maxPPT[i] <- mean(maxPPT_bind$sens_test)
    spec_sel_maxPPT[i] <- mean(maxPPT_bind$spec_test)
    
    SD_AUC_sel_maxPPT[i] <- sd(maxPPT_bind$AUC_test)
    SD_TSS_sel_maxPPT[i] <- sd(maxPPT_bind$TSS_test)
    SD_sens_sel_maxPPT[i] <- sd(maxPPT_bind$sens_test)
    SD_spec_sel_maxPPT[i] <- sd(maxPPT_bind$spec_test)
    
    #count how many points (and identify which ones by ID) were used in the training of selected models
    n_models <- ddply(maxPPT_bind, .(ID), nrow)
    names(n_models)[2] <- 'n_models'
    
    #calculate average SD and number of models in which the point was included (selected models)
    maxPPT_average <- ddply(maxPPT_bind, .(ID), summarise,
                             Species = unique(Species),
                             decimalLongitude = unique(decimalLongitude),
                             decimalLatitude = unique(decimalLatitude),
                             Occurrence = unique(Occurrence),
                             Type = unique(Type),
                             wc2.1_2.5m_bio_1 = unique(wc2.1_2.5m_bio_1),
                             wc2.1_2.5m_bio_13 = unique(wc2.1_2.5m_bio_13),
                             avg_Mean_T_SHAP = mean(Mean_T_SHAP),
                             SD_Mean_T_SHAP = sd(Mean_T_SHAP),
                             avg_Max_PPT_SHAP = mean(Max_PPT_SHAP),
                             SD_Max_PPT_SHAP = sd(Max_PPT_SHAP),
                             avg_AUC_test = mean(AUC_test),
                             SD_AUC_test = sd(AUC_test),
                             avg_TSS_test = mean(TSS_test),
                             SD_TSS_test = sd(TSS_test),
                             avg_sens_test = mean(sens_test),
                             SD_sens_test = sd(sens_test),
                             avg_spec_test = mean(spec_test),
                             SD_spec_test = sd(spec_test),
                             avg_AUC_train = mean(AUC_train),
                             SD_AUC_train = sd(AUC_train),
                             avg_TSS_train = mean(TSS_train),
                             SD_TSS_train = sd(TSS_train),
                             avg_sens_train = mean(sens_train),
                             SD_sens_train = sd(sens_train),
                             avg_spec_train = mean(spec_train),
                             SD_spec_train = sd(spec_train))
    
    #merge n_models with table with averages
    maxPPT_average <- try(merge(maxPPT_average, n_models, by = 'ID'), silent = T)
    
    #save result
    setwd(wd_res_species)
    write.csv(maxPPT_average, paste0(sps, '_maxPPT', '.csv'), row.names = F)
    
    
  }else{
    
    tot_models[i] <- NA
    
    perc_used_minT[i] <- NA
    AUC_all_minT[i] <- NA
    TSS_all_minT[i] <- NA
    sens_all_minT[i] <- NA
    spec_all_minT[i] <- NA
    AUC_all_minT[i] <- NA
    TSS_all_minT[i] <- NA
    sens_all_minT[i] <- NA
    spec_all_minT[i] <- NA
    SD_AUC_all_minT[i] <- NA
    SD_TSS_all_minT[i] <- NA
    SD_sens_all_minT[i] <- NA
    SD_spec_all_minT[i] <- NA
    AUC_sel_minT[i] <- NA
    TSS_sel_minT[i] <- NA
    sens_sel_minT[i] <- NA
    spec_sel_minT[i] <- NA
    SD_AUC_sel_minT[i] <- NA
    SD_TSS_sel_minT[i] <- NA
    SD_sens_sel_minT[i] <- NA
    SD_spec_sel_minT[i] <- NA
    
    perc_used_meanT[i] <- NA
    AUC_all_meanT[i] <- NA
    TSS_all_meanT[i] <- NA
    sens_all_meanT[i] <- NA
    spec_all_meanT[i] <- NA
    AUC_all_meanT[i] <- NA
    TSS_all_meanT[i] <- NA
    sens_all_meanT[i] <- NA
    spec_all_meanT[i] <- NA
    SD_AUC_all_meanT[i] <- NA
    SD_TSS_all_meanT[i] <- NA
    SD_sens_all_meanT[i] <- NA
    SD_spec_all_meanT[i] <- NA
    AUC_sel_meanT[i] <- NA
    TSS_sel_meanT[i] <- NA
    sens_sel_meanT[i] <- NA
    spec_sel_meanT[i] <- NA
    SD_AUC_sel_meanT[i] <- NA
    SD_TSS_sel_meanT[i] <- NA
    SD_sens_sel_meanT[i] <- NA
    SD_spec_sel_meanT[i] <- NA
    
    perc_used_maxT[i] <- NA
    AUC_all_maxT[i] <- NA
    TSS_all_maxT[i] <- NA
    sens_all_maxT[i] <- NA
    spec_all_maxT[i] <- NA
    AUC_all_maxT[i] <- NA
    TSS_all_maxT[i] <- NA
    sens_all_maxT[i] <- NA
    spec_all_maxT[i] <- NA
    SD_AUC_all_maxT[i] <- NA
    SD_TSS_all_maxT[i] <- NA
    SD_sens_all_maxT[i] <- NA
    SD_spec_all_maxT[i] <- NA
    AUC_sel_maxT[i] <- NA
    TSS_sel_maxT[i] <- NA
    sens_sel_maxT[i] <- NA
    spec_sel_maxT[i] <- NA
    SD_AUC_sel_maxT[i] <- NA
    SD_TSS_sel_maxT[i] <- NA
    SD_sens_sel_maxT[i] <- NA
    SD_spec_sel_maxT[i] <- NA
    
    perc_used_minPPT[i] <- NA
    AUC_all_minPPT[i] <- NA
    TSS_all_minPPT[i] <- NA
    sens_all_minPPT[i] <- NA
    spec_all_minPPT[i] <- NA
    AUC_all_minPPT[i] <- NA
    TSS_all_minPPT[i] <- NA
    sens_all_minPPT[i] <- NA
    spec_all_minPPT[i] <- NA
    SD_AUC_all_minPPT[i] <- NA
    SD_TSS_all_minPPT[i] <- NA
    SD_sens_all_minPPT[i] <- NA
    SD_spec_all_minPPT[i] <- NA
    AUC_sel_minPPT[i] <- NA
    TSS_sel_minPPT[i] <- NA
    sens_sel_minPPT[i] <- NA
    spec_sel_minPPT[i] <- NA
    SD_AUC_sel_minPPT[i] <- NA
    SD_TSS_sel_minPPT[i] <- NA
    SD_sens_sel_minPPT[i] <- NA
    SD_spec_sel_minPPT[i] <- NA
    
    perc_used_meanPPT[i] <- NA
    AUC_all_meanPPT[i] <- NA
    TSS_all_meanPPT[i] <- NA
    sens_all_meanPPT[i] <- NA
    spec_all_meanPPT[i] <- NA
    AUC_all_meanPPT[i] <- NA
    TSS_all_meanPPT[i] <- NA
    sens_all_meanPPT[i] <- NA
    spec_all_meanPPT[i] <- NA
    SD_AUC_all_meanPPT[i] <- NA
    SD_TSS_all_meanPPT[i] <- NA
    SD_sens_all_meanPPT[i] <- NA
    SD_spec_all_meanPPT[i] <- NA
    AUC_sel_meanPPT[i] <- NA
    TSS_sel_meanPPT[i] <- NA
    sens_sel_meanPPT[i] <- NA
    spec_sel_meanPPT[i] <- NA
    SD_AUC_sel_meanPPT[i] <- NA
    SD_TSS_sel_meanPPT[i] <- NA
    SD_sens_sel_meanPPT[i] <- NA
    SD_spec_sel_meanPPT[i] <- NA
    
    perc_used_maxPPT[i] <- NA
    AUC_all_maxPPT[i] <- NA
    TSS_all_maxPPT[i] <- NA
    sens_all_maxPPT[i] <- NA
    spec_all_maxPPT[i] <- NA
    AUC_all_maxPPT[i] <- NA
    TSS_all_maxPPT[i] <- NA
    sens_all_maxPPT[i] <- NA
    spec_all_maxPPT[i] <- NA
    SD_AUC_all_maxPPT[i] <- NA
    SD_TSS_all_maxPPT[i] <- NA
    SD_sens_all_maxPPT[i] <- NA
    SD_spec_all_maxPPT[i] <- NA
    AUC_sel_maxPPT[i] <- NA
    TSS_sel_maxPPT[i] <- NA
    sens_sel_maxPPT[i] <- NA
    spec_sel_maxPPT[i] <- NA
    SD_AUC_sel_maxPPT[i] <- NA
    SD_TSS_sel_maxPPT[i] <- NA
    SD_sens_sel_maxPPT[i] <- NA
    SD_spec_sel_maxPPT[i] <- NA
  }
  print(i)
}



# make table with info for all species
info_species <- data.frame(Species = sps_list, tot_models = tot_models,
                           
    #percentage of useful models (AUC >= 0.7 and TSS >= 0.4)
    Used_models_minT =  paste0(round(perc_used_minT, 2), ' %'),
    Used_models_meanT =  paste0(round(perc_used_meanT, 2), ' %'),
    Used_models_maxT =  paste0(round(perc_used_maxT, 2), ' %'),
    Used_models_minPPT =  paste0(round(perc_used_minPPT, 2), ' %'),
    Used_models_meanPPT =  paste0(round(perc_used_meanPPT, 2), ' %'),
    Used_models_maxPPT =  paste0(round(perc_used_maxPPT, 2), ' %'),
                
    #evaluation metrics for the species of the models (test metrics)
                
    #minT
    AUC_all_minT = AUC_all_minT, AUC_sel_minT = AUC_sel_minT,
    TSS_all_minT = TSS_all_minT, TSS_sel_minT = TSS_sel_minT,
    sens_all_minT = sens_all_minT, sens_sel_minT = sens_sel_minT,
    spec_all_minT = spec_all_minT, spec_sel_minT = spec_sel_minT,
    SD_AUC_all_minT = SD_AUC_all_minT, SD_AUC_sel_minT = SD_AUC_sel_minT,
    SD_TSS_all_minT = SD_TSS_all_minT, SD_TSS_sel_minT = SD_TSS_sel_minT,
    SD_sens_all_minT = SD_sens_all_minT, SD_sens_sel_minT = SD_sens_sel_minT,
    SD_spec_all_minT = SD_spec_all_minT, SD_spec_sel_minT = SD_spec_sel_minT,
                
    #meanT
    AUC_all_meanT = AUC_all_meanT, AUC_sel_meanT = AUC_sel_meanT,
    TSS_all_meanT = TSS_all_meanT, TSS_sel_meanT = TSS_sel_meanT,
    sens_all_meanT = sens_all_meanT, sens_sel_meanT = sens_sel_meanT,
    spec_all_meanT = spec_all_meanT, spec_sel_meanT = spec_sel_meanT,
    SD_AUC_all_meanT = SD_AUC_all_meanT, SD_AUC_sel_meanT = SD_AUC_sel_meanT,
    SD_TSS_all_meanT = SD_TSS_all_meanT, SD_TSS_sel_meanT = SD_TSS_sel_meanT,
    SD_sens_all_meanT = SD_sens_all_meanT, SD_sens_sel_meanT = SD_sens_sel_meanT,
    SD_spec_all_meanT = SD_spec_all_meanT, SD_spec_sel_meanT = SD_spec_sel_meanT,
          
    #maxT
    AUC_all_maxT = AUC_all_maxT, AUC_sel_maxT = AUC_sel_maxT,
    TSS_all_maxT = TSS_all_maxT, TSS_sel_maxT = TSS_sel_maxT,
    sens_all_maxT = sens_all_maxT, sens_sel_maxT = sens_sel_maxT,
    spec_all_maxT = spec_all_maxT, spec_sel_maxT = spec_sel_maxT,
    SD_AUC_all_maxT = SD_AUC_all_maxT, SD_AUC_sel_maxT = SD_AUC_sel_maxT,
    SD_TSS_all_maxT = SD_TSS_all_maxT, SD_TSS_sel_maxT = SD_TSS_sel_maxT,
    SD_sens_all_maxT = SD_sens_all_maxT, SD_sens_sel_maxT = SD_sens_sel_maxT,
    SD_spec_all_maxT = SD_spec_all_maxT, SD_spec_sel_maxT = SD_spec_sel_maxT,
          
    #minPPT
    AUC_all_minPPT = AUC_all_minPPT, AUC_sel_minPPT = AUC_sel_minPPT,
    TSS_all_minPPT = TSS_all_minPPT, TSS_sel_minPPT = TSS_sel_minPPT,
    sens_all_minPPT = sens_all_minPPT, sens_sel_minPPT = sens_sel_minPPT,
    spec_all_minPPT = spec_all_minPPT, spec_sel_minPPT = spec_sel_minPPT,
    SD_AUC_all_minPPT = SD_AUC_all_minPPT, SD_AUC_sel_minPPT = SD_AUC_sel_minPPT,
    SD_TSS_all_minPPT = SD_TSS_all_minPPT, SD_TSS_sel_minPPT = SD_TSS_sel_minPPT,
    SD_sens_all_minPPT = SD_sens_all_minPPT, SD_sens_sel_minPPT = SD_sens_sel_minPPT,
    SD_spec_all_minPPT = SD_spec_all_minPPT, SD_spec_sel_minPPT = SD_spec_sel_minPPT,
    
    #meanPPT
    AUC_all_meanPPT = AUC_all_meanPPT, AUC_sel_meanPPT = AUC_sel_meanPPT,
    TSS_all_meanPPT = TSS_all_meanPPT, TSS_sel_meanPPT = TSS_sel_meanPPT,
    sens_all_meanPPT = sens_all_meanPPT, sens_sel_meanPPT = sens_sel_meanPPT,
    spec_all_meanPPT = spec_all_meanPPT, spec_sel_meanPPT = spec_sel_meanPPT,
    SD_AUC_all_meanPPT = SD_AUC_all_meanPPT, SD_AUC_sel_meanPPT = SD_AUC_sel_meanPPT,
    SD_TSS_all_meanPPT = SD_TSS_all_meanPPT, SD_TSS_sel_meanPPT = SD_TSS_sel_meanPPT,
    SD_sens_all_meanPPT = SD_sens_all_meanPPT, SD_sens_sel_meanPPT = SD_sens_sel_meanPPT,
    SD_spec_all_meanPPT = SD_spec_all_meanPPT, SD_spec_sel_meanPPT = SD_spec_sel_meanPPT,
    
    #maxPPT
    AUC_all_maxPPT = AUC_all_maxPPT, AUC_sel_maxPPT = AUC_sel_maxPPT,
    TSS_all_maxPPT = TSS_all_maxPPT, TSS_sel_maxPPT = TSS_sel_maxPPT,
    sens_all_maxPPT = sens_all_maxPPT, sens_sel_maxPPT = sens_sel_maxPPT,
    spec_all_maxPPT = spec_all_maxPPT, spec_sel_maxPPT = spec_sel_maxPPT,
    SD_AUC_all_maxPPT = SD_AUC_all_maxPPT, SD_AUC_sel_maxPPT = SD_AUC_sel_maxPPT,
    SD_TSS_all_maxPPT = SD_TSS_all_maxPPT, SD_TSS_sel_maxPPT = SD_TSS_sel_maxPPT,
    SD_sens_all_maxPPT = SD_sens_all_maxPPT, SD_sens_sel_maxPPT = SD_sens_sel_maxPPT,
    SD_spec_all_maxPPT = SD_spec_all_maxPPT, SD_spec_sel_maxPPT = SD_spec_sel_maxPPT)

    
#save mistressfile
setwd(wd_tables)
write.csv(info_species, 'Info_species.csv', row.names = F)
