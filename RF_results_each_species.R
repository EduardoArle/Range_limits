#load packages
library(randomForest)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240603'
  
#read temperature results table

setwd(wd_tables)
res_temp <- read.csv('Temperature_Rel_Polar_all_points.csv')

#calculate species' latitudinal amplitude
res_temp$latAmplitude <- res_temp$maxLat - res_temp$minLat

#calculate species' elevation amplitude
res_temp$elevAmplitude <- res_temp$max_elev - res_temp$min_elev

#calculate species' mean latitude
res_temp$meanLat <- (res_temp$max_elev + res_temp$min_elev) / 2

# #calculate 1/SE * slope
# res_temp$one_SE_x_slope_minT_RP <- (1/res_temp$minTSE_RP) * res_temp$minTslope_RP 

#discard species with less than 10 points
res_temp <- res_temp[res_temp$n_records > 9,]

#discard species that cross the equator
res_temp_one_hemis <- res_temp[res_temp$hemisphere != 'Both',]

#run RF 10 times
minT <- randomForest(minTslope_RP ~ n_records + 
                                     order +
                                     bodymass +
                                     n_biomes +
                                     hemisphere +
                                     rangeSize + 
                                     rangeRoundness +
                                     latAmplitude +
                                     elevAmplitude +
                                     meanLat +
                                     mean_elev,
                          data = res_temp_one_hemis, ntree = 5000,
                          keep.forest = FALSE,
                          importance = TRUE,
                          na.action = na.omit)



varImpPlot(minT)

#save in the correspondent folder in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2_species'


meanT <- randomForest(meanTslope_RP ~ n_records +
                        order +
                        bodymass +
                        n_biomes +
                        hemisphere +
                        rangeSize + 
                        rangeRoundness +
                        latAmplitude +
                        elevAmplitude +
                        meanLat +
                        mean_elev,
                     data = res_temp_one_hemis, ntree = 5000,
                     keep.forest = FALSE,
                     importance = TRUE,
                     na.action = na.omit)



varImpPlot(meanT)


maxT <- randomForest(maxTslope_RP ~ n_records + 
                       order +
                       bodymass +
                       n_biomes +
                       hemisphere +
                       rangeSize + 
                       rangeRoundness +
                       latAmplitude +
                       elevAmplitude +
                       meanLat +
                       mean_elev,
                      data = res_temp_one_hemis, ntree = 5000,
                      keep.forest = FALSE,
                      importance = TRUE,
                      na.action = na.omit)



varImpPlot(maxT)


#read precipitation results table

setwd(wd_tables)
res_prec <- read.csv('Precipitation_Rel_Polar_all_points.csv')

#calculate species' latitudinal amplitude
res_prec$latAmplitude <- res_prec$maxLat - res_prec$minLat

#calculate species' elevation amplitude
res_prec$elevAmplitude <- res_prec$max_elev - res_prec$min_elev

#calculate species' mean latitude
res_prec$meanLat <- (res_prec$max_elev + res_prec$min_elev) / 2

# #calculate 1/SE * slope
# res_temp$one_SE_x_slope_minT_RP <- (1/res_temp$minTSE_RP) * res_temp$minTslope_RP 

#discard species that cross the equator
res_prec_one_hemis <- res_prec[res_prec$hemisphere != 'Both',]



minPPT <- randomForest(minPPTslope_RP ~ n_records + 
                         order +
                         bodymass +
                         n_biomes +
                         hemisphere +
                         rangeSize + 
                         rangeRoundness +
                         latAmplitude +
                         elevAmplitude +
                         meanLat +
                         mean_elev,
                     data = res_prec_one_hemis, ntree = 5000,
                     keep.forest = FALSE,
                     importance = TRUE,
                     na.action = na.omit)



varImpPlot(minPPT)


meanPPT <- randomForest(meanPPTslope_RP ~ n_records + 
                          order +
                          bodymass +
                          n_biomes +
                          hemisphere +
                          rangeSize + 
                          rangeRoundness +
                          latAmplitude +
                          elevAmplitude +
                          meanLat +
                          mean_elev,
                      data = res_prec_one_hemis, ntree = 5000,
                      keep.forest = FALSE,
                      importance = TRUE,
                      na.action = na.omit)



varImpPlot(meanPPT)


maxPPT <- randomForest(maxPPTslope_RP ~ n_records + 
                         order +
                         bodymass +
                         n_biomes +
                         hemisphere +
                         rangeSize + 
                         rangeRoundness +
                         latAmplitude +
                         elevAmplitude +
                         meanLat +
                         mean_elev,
                     data = res_prec_one_hemis, ntree = 5000,
                     keep.forest = FALSE,
                     importance = TRUE,
                     na.action = na.omit)



varImpPlot(maxPPT)



#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2_species'  1000 size


# %IncMSE is the most robust and informative measure. It is the increase in mse of predictions(estimated with out-of-bag-CV) as a result of variable j being permuted(values randomly shuffled).

# IncNodePurity relates to the loss function which by best splits are chosen. The loss function is mse for regression and gini-impurity for classification. More useful variables achieve higher increases in node purities, that is to find a split which has a high inter node 'variance' and a small intra node 'variance'. IncNodePurity is biased and should only be used if the extra computation time of calculating %IncMSE is unacceptable. 




##########################################



# weighted SHAP minTslopeRP X rangeSize

plot(res_temp_one_hemis$rangeSize, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'Range Size')

lm_weighted_minT_range_Size <- lm(minTslope_RP ~ rangeSize,
                                  weights = minTSE_RP,
                                  data = res_temp_one_hemis)

abline(lm_weighted_minT_range_Size , col = '#0000FF', lwd = 2)

summary(lm_minT_range_Size)


# SHAP minTslopeRP X log rangeSize

plot(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'log Range Size')

lm_minT_log_range_Size <- lm(minTslope_RP ~ log(rangeSize),
                             data = res_temp_one_hemis)

abline(lm_minT_log_range_Size, col = '#0000FF', lwd = 2)

summary(lm_minT_log_range_Size)


#weighted SHAP minTslopeRP X log rangeSize

plot(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP weighted',
     xlab = 'log Range Size')

lm_weigtehd_minT_log_range_Size <- lm(minTslope_RP ~ log(rangeSize),
                                   weights = minTSE_RP,
                                   data = res_temp_one_hemis)

abline(lm_weigtehd_minT_log_range_Size, col = '#0000FF', lwd = 2)

summary(lm_weigtehd_minT_log_range_Size)


##### LATITUDINAL AMPLITUDE ####

# SHAP minTslopeRP X rangeSize

plot(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'Latitudinal Amplitude')

lm_minT_lat_amplitude <- lm(minTslope_RP ~ latAmplitude,
                            data = res_temp_one_hemis)

abline(lm_minT_lat_amplitude , col = '#0000FF', lwd = 2)

summary(lm_minT_lat_amplitude)


# weighted SHAP minTslopeRP X rangeSize

plot(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP weighted',
     xlab = 'Latitudinal Amplitude')

lm_weigtehd_minT_lat_amplitude <- lm(minTslope_RP ~ latAmplitude,
                            weights = minTSE_RP,
                            data = res_temp_one_hemis)

abline(lm_weigtehd_minT_lat_amplitude, col = '#0000FF', lwd = 2)

summary(lm_weigtehd_minT_lat_amplitude)


# SHAP minTslopeRP X log lat amplitude

plot(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'log Latitudinal Amplitude')

lm_minT_log_lat_amplitude <- lm(minTslope_RP ~ log(latAmplitude),
                                data = res_temp_one_hemis)

abline(lm_minT_log_lat_amplitude , col = '#0000FF', lwd = 2)

summary(lm_minT_log_lat_amplitude)


#weighted SHAP minTslopeRP X log lat amplitude

plot(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP weighted',
     xlab = 'log Latitudinal Amplitude')

lm_weigtehd_minT_log_lat_amplitude <- lm(minTslope_RP ~ log(latAmplitude),
                                        weights = minTSE_RP,
                                        data = res_temp_one_hemis)

abline(lm_minT_log_lat_amplitude , col = '#0000FF', lwd = 2)

summary(lm_weigtehd_minT_log_lat_amplitude)



##### ELEVATION AMPLITUDE ####


# SHAP minTslopeRP X elev amplitude

plot(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'Elevation Amplitude')

lm_minT_elev_amplitude <- lm(minTslope_RP ~ elevAmplitude,
                             data = res_temp_one_hemis)

abline(lm_minT_elev_amplitude , col = '#0000FF', lwd = 2)

summary(lm_minT_elev_amplitude)


# SHAP minTslopeRP X log elev amplitude

plot(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP',
     xlab = 'log Elevation Amplitude')

lm_minT_log_elev_amplitude <- lm(minTslope_RP ~ log(elevAmplitude),
                                 data = res_temp_one_hemis)

abline(lm_minT_log_elev_amplitude, col = '#0000FF', lwd = 2)

summary(lm_minT_log_elev_amplitude)

#weighted SHAP minTslopeRP X log elev amplitude

plot(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP weighted',
     xlab = 'log Elevation Amplitude')

lm_weighted_minT_log_elev_amplitude <- lm(minTslope_RP ~ log(elevAmplitude),
                                          weights = minTSE_RP,
                                          data = res_temp_one_hemis)

abline(lm_weighted_minT_log_elev_amplitude, col = '#0000FF', lwd = 2)


summary(lm_minT_range_Size)

#weighted SHAP minTslopeRP X elev amplitude

plot(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'minT slope RP weighted',
     xlab = 'log Elevation Amplitude')

lm_weighted_minT_elev_amplitude <- lm(minTslope_RP ~ elevAmplitude,
                                          weights = minTSE_RP,
                                          data = res_temp_one_hemis)

abline(lm_weighted_minT_elev_amplitude, col = '#0000FF', lwd = 2)


summary(lm_weighted_minT_elev_amplitude)




##MAKE WEIGHTED PLOTS for each of the predictions

slopes

point levels (gam instead of linear models)

t <- lm(res_temp$minTslope_RP ~ res_temp$rangeSize)  #+ elevation + latamp + W1SD)
plot(t)
plot(res_temp$minTslope_RP ~ res_temp$rangeSize)

summary(t)
