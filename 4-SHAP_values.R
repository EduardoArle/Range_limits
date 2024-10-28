#load packages 
library(raster);library(data.table);library(rworldmap)
library(dplyr);library(sf);library(shapviz);library(kernelshap)
library(xgboost);library(caret);library(pROC);library(PresenceAbsence)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_variables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_res_species <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Comparison_20240926'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Tables'


#list species 
setwd(wd_thinned_occ)
sps_list <- gsub('_thinned.csv', '', list.files())

#load all six BioCLim variables being used
setwd(wd_variables)
AnnualMeanTemperature <- raster('wc2.1_2.5m_bio_1.tif')
MaxTemperatureOfWarmestMonth <- raster('wc2.1_2.5m_bio_5.tif')
MinTemperatureOfColdestMonth <- raster('wc2.1_2.5m_bio_6.tif')
AnnualPrecipitation <- raster('wc2.1_2.5m_bio_12.tif')
PrecipitationOfWettestMonth <- raster('wc2.1_2.5m_bio_13.tif')
PrecipitationOfDriestMonth <- raster('wc2.1_2.5m_bio_14.tif')

## List variables we will use for each predictions

#elevation and mean prec to compare temperature variables
preds_min_temp <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_6')
preds_mean_temp  <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')
preds_max_temp <- c('wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_5')

#elevation and mean temp to compare precipitation variables
preds_min_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_14')
preds_mean_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_12')
preds_max_PPT <- c('wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_13')


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
  
  ## Create pseudo absences (half the number of presences)
  
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
                 MaxTemperatureOfWarmestMonth, PrecipitationOfWettestMonth) #maxs
  
  #extract values from each location from all variables
  vals_pts <- extract(preds, species_sp)
  
  #make a table creating an ID for each point
  tab_occ_vars <- cbind(ID = c(1:nrow(species)), species, vals_pts)
  
  ## Split the data in 10 folds for training and testing
  
  #randomly shuffle the data separate presence from pseudo-absence
  shuffle_pr <- tab_occ_vars[tab_occ_vars$Type == 'Presence',]
  shuffled_pr <- shuffle_pr[sample(nrow(shuffle_pr)),]
  
  shuffle_pa <- tab_occ_vars[tab_occ_vars$Type == 'Pseudo-absence',]
  shuffled_pa <- shuffle_pa[sample(nrow(shuffle_pa)),]
  
  #create 10 equally size folds
  folds <- cut(seq(1, nrow(shuffled_pr)), breaks=10, labels=FALSE)
  
  #perform 10 fold cross validation for each model
  for(j in 1:10)
  {
    if(TRUE %in% unique(folds == j)){
      #select test data
      testIndexes <- which(folds == j)
      testData_pr <- shuffled_pr[testIndexes,]
      testData_pa <- shuffled_pa[testIndexes,]
      testData <- rbind(testData_pr, testData_pa)
      
      #select train data
      trainData_pr <- shuffled_pr[-testIndexes, ]
      trainData_pa <- shuffled_pa[-testIndexes, ]
      trainData <- rbind(trainData_pr, trainData_pa)
      
      ## Fit XGBoost models
      
      ## Mean PPT to compare T variables importance
      
      ##############################################
      ############## mean PPT + min T ##############
      ##############################################
      
      #prepare train data
      shap_T_min <- xgb.DMatrix(data.matrix(trainData[preds_min_temp]),
                                label = trainData$Occurrence)
      
      #prepare test data
      shap_T_min_test <- xgb.DMatrix(data.matrix(testData[preds_min_temp]),
                                     label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_T_min,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_T_min) 
      
      #predict values for test
      pred_test <- predict(fit, shap_T_min_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_min_temp]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_prec_SHAP <- character()  #bio12
      MIN_temp_SHAP <- character()  #bio6
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
        MIN_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_6']
        
        print(k)
      }
      
      #make a data.frame with SHAP results and model evaluations
      meanPPT_minT <- cbind(Species = sps_list[i],
                            trainData[,c('ID',
                                         'decimalLongitude', 'decimalLatitude',
                                         'Occurrence', 'Type',
                                         'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_6')],
                            data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                       Min_T_SHAP = MIN_temp_SHAP,
                                       iteration = j,
                                       AUC_test = AUC_test,
                                       TSS_test = TSS_test,
                                       sens_test = sens_test,
                                       spec_test = spec_test,
                                       AUC_train = AUC_train,
                                       TSS_train = TSS_train,
                                       sens_train = sens_train,
                                       spec_train = spec_train))
      
      #save results per species
      
      #create directory to save each repetition for the species
      dir_species <- paste0(wd_res_species, '/', sps_list[i])
      dir.create(dir_species)
      
      #save results
      setwd(dir_species)
      write.csv(meanPPT_minT, paste0(sps, '_minT_', j, '.csv'), row.names = F)
      
      
      ###############################################
      ############## mean PPT + mean T ##############
      ###############################################
      
      #prepare train data
      shap_T_mean <- xgb.DMatrix(data.matrix(trainData[preds_mean_temp]),
                                 label = trainData$Occurrence)
      
      #prepare test data
      shap_T_mean_test <- xgb.DMatrix(data.matrix(testData[preds_mean_temp]),
                                      label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_T_mean,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_T_mean) 
      
      #predict values for test
      pred_test <- predict(fit, shap_T_mean_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_mean_temp]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_prec_SHAP <- character()  #bio12
      MEAN_temp_SHAP <- character()  #bio1
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
        MEAN_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
        
        print(k)
      }
      
      #make a data.frame with SHAP results and model evaluations
      meanPPT_meanT <- cbind(Species = sps_list[i],
                             trainData[,c('ID',
                                          'decimalLongitude', 'decimalLatitude',
                                          'Occurrence', 'Type',
                                          'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1')],
                             data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                        Mean_T_SHAP = MEAN_temp_SHAP,
                                        iteration = j,
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
      write.csv(meanPPT_meanT, paste0(sps, '_meanT_', j, '.csv'), row.names = F)
      
      
      ##############################################
      ############## mean PPT + max T ##############
      ##############################################
      
      #prepare train data
      shap_T_max <- xgb.DMatrix(data.matrix(trainData[preds_max_temp]),
                                label = trainData$Occurrence)
      
      #prepare test data
      shap_T_max_test <- xgb.DMatrix(data.matrix(testData[preds_max_temp]),
                                     label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_T_max,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_T_max) 
      
      #predict values for test
      pred_test <- predict(fit, shap_T_max_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_max_temp]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_prec_SHAP <- character()  #bio12
      MAX_temp_SHAP <- character()  #bio1
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
        MAX_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_5']
        
        print(k)
      }
      
      #make a data.frame with SHAP results and model evaluations
      meanPPT_maxT <- cbind(Species = sps_list[i],
                            trainData[,c('ID',
                                         'decimalLongitude', 'decimalLatitude',
                                         'Occurrence', 'Type',
                                         'wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_5')],
                            data.frame(Mean_PPT_SHAP = MEAN_prec_SHAP,
                                       Max_T_SHAP = MAX_temp_SHAP,
                                       iteration = j,
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
      write.csv(meanPPT_maxT, paste0(sps, '_maxT_', j, '.csv'), row.names = F)
      
      
      
      
      ## mean T to compare PPT variables importance
      
      
      
      
      ##############################################
      ############## mean T + min PPT ##############
      ##############################################
      
      
      
      
      #prepare train data
      shap_PPT_min <- xgb.DMatrix(data.matrix(trainData[preds_min_PPT]),
                                  label = trainData$Occurrence)
      
      #prepare test data
      shap_PPT_min_test <- xgb.DMatrix(data.matrix(testData[preds_min_PPT]),
                                       label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_PPT_min,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_PPT_min) 
      
      #predict values for test
      pred_test <- predict(fit, shap_PPT_min_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_min_PPT]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_temp_SHAP <- character()  #bio1
      MIN_prec_SHAP <- character()  #bio14
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
        MIN_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_14']
        
        print(k)
      }
      
      
      #make a data.frame with SHAP results and model evaluations
      meanT_minPPT <- cbind(Species = sps_list[i],
                            trainData[,c('ID',
                                         'decimalLongitude', 'decimalLatitude',
                                         'Occurrence', 'Type',
                                         'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_14')],
                            data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                       Min_PPT_SHAP = MIN_prec_SHAP,
                                       iteration = j,
                                       AUC_test = AUC_test,
                                       TSS_test = TSS_test,
                                       sens_test = sens_test,
                                       spec_test = spec_test,
                                       AUC_train = AUC_train,
                                       TSS_train = TSS_train,
                                       sens_train = sens_train,
                                       spec_train = spec_train))
      
      #save results per species
      
      #save results
      setwd(dir_species)
      write.csv(meanT_minPPT, paste0(sps, '_minPPT_', j, '.csv'), row.names = F)
      
      
      
      
      ###############################################
      ############## mean T + mean PPT ##############
      ###############################################
      
      
      
      
      #prepare train data
      shap_PPT_mean <- xgb.DMatrix(data.matrix(trainData[preds_mean_PPT]),
                                   label = trainData$Occurrence)
      
      #prepare test data
      shap_PPT_mean_test <- xgb.DMatrix(data.matrix(testData[preds_mean_PPT]),
                                        label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_PPT_mean,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_PPT_mean) 
      
      #predict values for test
      pred_test <- predict(fit, shap_PPT_mean_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_mean_PPT]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_temp_SHAP <- character()  #bio1
      MEAN_prec_SHAP <- character()  #bio12
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
        MEAN_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_12']
        
        print(k)
      }
      
      #make a data.frame with SHAP results and model evaluations
      meanT_meanPPT <- cbind(Species = sps_list[i],
                             trainData[,c('ID',
                                          'decimalLongitude', 'decimalLatitude',
                                          'Occurrence', 'Type',
                                          'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_12')],
                             data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                        Mean_PPT_SHAP = MEAN_prec_SHAP,
                                        iteration = j,
                                        AUC_test = AUC_test,
                                        TSS_test = TSS_test,
                                        sens_test = sens_test,
                                        spec_test = spec_test,
                                        AUC_train = AUC_train,
                                        TSS_train = TSS_train,
                                        sens_train = sens_train,
                                        spec_train = spec_train))
      
      #save results per species
      
      #save results
      setwd(dir_species)
      write.csv(meanT_meanPPT, paste0(sps, '_meanPPT_', j, '.csv'), row.names = F)
      
      
      
      
      ##############################################
      ############## mean T + max PPT ##############
      ##############################################
      
      
      
      
      #prepare train data
      shap_PPT_max <- xgb.DMatrix(data.matrix(trainData[preds_max_PPT]),
                                  label = trainData$Occurrence)
      
      #prepare test data
      shap_PPT_max_test <- xgb.DMatrix(data.matrix(testData[preds_max_PPT]),
                                       label = testData$Occurrence)
      
      #fit the model with train data
      fit <- xgb.train(
        params = list(learning_rate = 0.1, objective = "reg:squarederror"), 
        data =   shap_PPT_max,
        nrounds = 65L
      )
      
      #predict values for train
      pred_train <- predict(fit, shap_PPT_max) 
      
      #predict values for test
      pred_test <- predict(fit, shap_PPT_max_test) 
      
      #get evaluation metrics
      
      #AUC
      AUC_train <- as.numeric(pROC::auc(trainData$Occurrence, pred_train))
      AUC_test <- as.numeric(pROC::auc(testData$Occurrence, pred_test))
      
      #sens, spec, TSS, th = max(sens+spec)
      
      #make tables with observed and predicted values
      obs_pred_train <- data.frame(plotID = c(1:length(trainData$Occurrence)),
                                   Observed = as.integer(trainData$Occurrence),
                                   Predicted = pred_train)
      
      obs_pred_test <- data.frame(plotID = c(1:length(testData$Occurrence)),
                                  Observed = as.integer(testData$Occurrence),
                                  Predicted = pred_test)
      
      #calculate th = max(sens+spec)
      th <- optimal.thresholds(obs_pred_train)[3,2]
      
      #create confusion matrix based on th = max(sens+spec)
      confusion_train <- cmx(obs_pred_train, threshold = th)
      confusion_test <- cmx(obs_pred_test, threshold = th)
      
      #calculate sensitivity
      sens_train <- sensitivity(confusion_train, st.dev = F)
      sens_test <- sensitivity(confusion_test, st.dev = F)
      
      #calculate specificity
      spec_train <- specificity(confusion_train, st.dev = F)
      spec_test <- specificity(confusion_test, st.dev = F)
      
      #calculate TSS
      TSS_train <- sens_train + spec_train - 1
      TSS_test <- sens_test + spec_test - 1
      
      # We also pass feature data X with originally encoded values
      shp <- shapviz(fit, X_pred = data.matrix(trainData[preds_max_PPT]),
                     X = trainData)
      
      # get the SHAP values for each variable in each prediction 
      MEAN_temp_SHAP <- character()  #bio1
      MAX_prec_SHAP <- character()  #bio13
      
      for(k in 1:nrow(trainData))
      {
        obj <- sv_waterfall(shp, row_id = k) +
          theme(axis.text = element_text(size = 11))
        
        MEAN_temp_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_1']
        MAX_prec_SHAP[k] <- obj$data$S[row.names(obj$data) == 'wc2.1_2.5m_bio_13']
        
        print(k)
      }
      
      
      #make a data.frame with SHAP results and model evaluations
      meanT_maxPPT <- cbind(Species = sps_list[i],
                            trainData[,c('ID',
                                         'decimalLongitude', 'decimalLatitude',
                                         'Occurrence', 'Type',
                                         'wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_13')],
                            data.frame(Mean_T_SHAP = MEAN_temp_SHAP,
                                       Max_PPT_SHAP = MAX_prec_SHAP,
                                       iteration = j,
                                       AUC_test = AUC_test,
                                       TSS_test = TSS_test,
                                       sens_test = sens_test,
                                       spec_test = spec_test,
                                       AUC_train = AUC_train,
                                       TSS_train = TSS_train,
                                       sens_train = sens_train,
                                       spec_train = spec_train))
      
      #save results per species
      
      #save results
      setwd(dir_species)
      write.csv(meanT_maxPPT, paste0(sps, '_maxPPT_', j, '.csv'), row.names = F)
    }
  }
}




######## Calculate average SHAP for each point and the corresponding metrics ######

#list with all 503 species
sps_list

#number of records for each species
setwd(wd_thinned_occ)
sps_occ <- lapply(list.files(), read.csv)
n_occ <- sapply(sps_occ, nrow)

#create an empty objects to populate with the number of folds for training and testing
folds <- numeric()

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


for(i in 63:length(sps_list))
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
    folds[i] <- length(minT)
    
    #input the percentage of models that could be used in the calculation of the SHAP. This number will vary across variable combination models because it depends on how many models were satisfactory according to AUC and TSS
    perc_used_minT[i] <- sum(sapply(minT, function(x) unique(x$TSS_test)) >= 0.4 &
                            sapply(minT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    perc_used_meanT[i] <- sum(sapply(meanT, function(x) unique(x$TSS_test)) >= 0.4 &
                               sapply(meanT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    perc_used_maxT[i] <- sum(sapply(maxT, function(x) unique(x$TSS_test)) >= 0.4 &
                                sapply(maxT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    perc_used_minPPT[i] <- sum(sapply(minPPT, function(x) unique(x$TSS_test)) >= 0.4 &
                               sapply(minPPT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    perc_used_meanPPT[i] <- sum(sapply(meanPPT, function(x) unique(x$TSS_test)) >= 0.4 &
                                 sapply(meanPPT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    perc_used_maxPPT[i] <- sum(sapply(maxPPT, function(x) unique(x$TSS_test)) >= 0.4 &
                                  sapply(maxPPT, function(x) unique(x$AUC_test)) >= 0.7) /
      folds[i] * 100
    
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
    
    folds[i] <- NA
    
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
info_species <- data.frame(Species = sps_list, n_occ = n_occ, folds = folds,
                           
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

    

#.  INCLUDE CORRELATION ANALYSIS IN FINAL TABLE.   

#load dataframe with variable correlation
setwd(wd_tables)
results_correl <- read.csv('Correlation_variables.csv')


#join correl analysis with the info species table
mistress_file <- merge(info_species, results_correl, by = c('Species', 'n_occ'))

#save mistressfile
setwd(wd_tables)
write.csv(mistress_file, 'Mistressfile.csv', row.names = F)







##################.  NOTE  #########################





###################### SCRAP #########################

#list the species that could produce any models (they have a directory inside of the main directory)
setwd(wd_res_species)
sps_list_models <- list.dirs(full.names = F, recursive = F)
