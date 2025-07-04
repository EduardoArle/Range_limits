#load packages
library(randomForest)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'

#read results table
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#run RF minT

#select only species that had lower correl between vars
results_minT <- results[abs(results$Cor_vars_minT) <= 0.7,]
results_minT <- results_minT[complete.cases(results_minT$Cor_vars_minT),]

minT <- randomForest(avg_Min_T_SHAP ~
                        decimalLatitude +
                        species +
                        rangeSize +
                        rangeLoc +
                        roundness +
                        distEdge +
                        relPolarwardness +
                        elevation +
                        biome +
                        bodyMass +
                        nOcc +
                        order,
                      data = results_minT, ntree = 1000,
                      keep.forest = FALSE,
                      na.action = na.omit,
                      importance = TRUE)



# %IncMSE is the most robust and informative measure. It is the increase in mse of predictions(estimated with out-of-bag-CV) as a result of variable j being permuted(values randomly shuffled).

# IncNodePurity relates to the loss function which by best splits are chosen. The loss function is mse for regression and gini-impurity for classification. More useful variables achieve higher increases in node purities, that is to find a split which has a high inter node 'variance' and a small intra node 'variance'. IncNodePurity is biased and should only be used if the extra computation time of calculating %IncMSE is unacceptable. 

varImpPlot(minT)

#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2'  1000 size

minT$importanceSD
minT$importance


#run RF meanT

#select only species that had lower correl between vars
results_meanT <- results[abs(results$Cor_vars_meanT) <= 0.7,]
results_meanT <- results_meanT[complete.cases(results_meanT$Cor_vars_meanT),]

meanT <- randomForest(avg_Mean_T_SHAP ~
                       decimalLatitude +
                       species +
                       rangeSize +
                       roundness +
                       distEdge +
                       rangeLoc +
                       relPolarwardness +
                       elevation +
                       biome +
                       bodyMass +
                       nOcc +
                       order,
                     data = results_meanT, ntree = 1000,
                     keep.forest = FALSE,
                     na.action = na.omit,
                     importance = TRUE)

varImpPlot(meanT)


#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2'  1000 size

#run RF maxT

#select only species that had lower correl between vars
results_maxT <- results[abs(results$Cor_vars_maxT) <= 0.7,]
results_maxT <- results_meanT[complete.cases(results_maxT$Cor_vars_maxT),]

maxT <- randomForest(avg_Max_T_SHAP ~
                        decimalLatitude +
                        species +
                        rangeSize +
                        roundness +
                        distEdge +
                        rangeLoc +
                        relPolarwardness +
                        elevation +
                        biome +
                        bodyMass +
                        nOcc +
                        order,
                      data = results_maxT, ntree = 1000,
                      keep.forest = FALSE,
                      na.action = na.omit,
                      importance = TRUE)

varImpPlot(maxT)


#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2'  1000 size


#run RF minPPT
minPPT <- randomForest(Min_PPT_SHAP ~
                       decimalLatitude +
                       species +
                       rangeSize +
                       roundness +
                       distEdge +
                       absPolarwardness +
                       relPolarwardness +
                       elevation +
                       biome +
                       bodyMass +
                       n_occ +
                       order,
                     data = results, ntree = 1000,
                     keep.forest = FALSE,
                     na.action = na.omit,
                     importance = TRUE)

varImpPlot(minPPT)


#run RF meanPPT
meanPPT <- randomForest(Mean_PPT_SHAP ~
                         decimalLatitude +
                         species +
                         rangeSize +
                         roundness +
                         distEdge +
                         absPolarwardness +
                         relPolarwardness +
                         elevation +
                         biome +
                         bodyMass +
                         n_occ +
                         order,
                       data = results, ntree = 1000,
                       keep.forest = FALSE,
                       na.action = na.omit,
                       importance = TRUE)

varImpPlot(meanPPT)


#run RF maxPPT
maxPPT <- randomForest(Max_PPT_SHAP ~
                          decimalLatitude +
                          species +
                          rangeSize +
                          roundness +
                          distEdge +
                          absPolarwardness +
                          relPolarwardness +
                          elevation +
                          biome +
                          bodyMass +
                          n_occ +
                          order,
                        data = results, ntree = 1000,
                        keep.forest = FALSE,
                        na.action = na.omit,
                        importance = TRUE)

varImpPlot(maxPPT)

