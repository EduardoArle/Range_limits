#load packages
library(randomForest); library(glmmTMB)

#list wds

wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'

#read temperature results table

setwd(wd_tables)
results <- read.csv('Results_all_sps.csv')

#run RF minT
minT <- randomForest(Min_T_SHAP ~
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



# %IncMSE is the most robust and informative measure. It is the increase in mse of predictions(estimated with out-of-bag-CV) as a result of variable j being permuted(values randomly shuffled).

# IncNodePurity relates to the loss function which by best splits are chosen. The loss function is mse for regression and gini-impurity for classification. More useful variables achieve higher increases in node purities, that is to find a split which has a high inter node 'variance' and a small intra node 'variance'. IncNodePurity is biased and should only be used if the extra computation time of calculating %IncMSE is unacceptable. 

varImpPlot(minT)

#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2'  1000 size

minT$importanceSD
minT$importance


#run RF meanT
meanT <- randomForest(Mean_T_SHAP ~
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

varImpPlot(meanT)


#save in '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/RF_2'  1000 size

#run RF maxT
maxT <- randomForest(Max_T_SHAP ~
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

