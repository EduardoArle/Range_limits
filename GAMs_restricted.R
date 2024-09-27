#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'
wd_models_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models_restricted'

#read temperature results table

setwd(wd_tables)
results <- read.csv('Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,11)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

#select only ranges larger than 1M km2
results_more_1Mkm2 <- results[results$rangeSize >= 1000000,]

#select only ranges smaller than 1M km2
results_less_1Mkm2 <- results[results$rangeSize <= 1000000,]

#check number of occ and number os species
nrow(results_more_1Mkm2)
nrow(results_less_1Mkm2)

length(unique(results_more_1Mkm2$species))
length(unique(results_less_1Mkm2$species))

# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT (more 1,000,000)

#model G
minT_relPolewarness_G_more_1Mkm2 <- gam(Min_T_SHAP ~
                                        s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                      data = results_more_1Mkm2,
                                      method="REML",
                                      family="gaussian")


summary(minT_relPolewarness_G_more_1Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_G_more_1Mkm2, 'minT_relPolewarness_G_more_1Mkm2')


plot.gam(minT_relPolewarness_G_more_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPolewarness_GS_more_1Mkm2 <- gam(Min_T_SHAP ~
                                        s(relPolewardness, k=4, m=2) 
                                        + s(relPolewardness, species, k=4, bs="fs", m=2),
                                       data = results_more_1Mkm2,
                                       method="REML")


summary(minT_relPolewarness_GS_more_1Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_GS_more_1Mkm2, 'minT_relPolewarness_GS_more_1Mkm2')


plot.gam(minT_relPolewarness_GS_more_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### minT (less 1,000,000)

minT_relPolewarness_G_less_1Mkm2 <- gam(Min_T_SHAP ~
                                            s(relPolewardness, k=4, bs="tp") 
                                          + s(species, k=503, bs="re"),
                                          data = results_less_1Mkm2,
                                          method="REML",
                                          family="gaussian")


summary(minT_relPolewarness_G_less_1Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_G_less_1Mkm2, 'minT_relPolewarness_G_less_1Mkm2')


plot.gam(minT_relPolewarness_G_less_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPolewarness_GS_less_1Mkm2 <- gam(Min_T_SHAP ~
                                             s(relPolewardness, k=4, m=2) 
                                           + s(relPolewardness, species, k=4, bs="fs", m=2),
                                           data = results_less_1Mkm2,
                                           method="REML")


summary(minT_relPolewarness_GS_less_1Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_GS_less_1Mkm2, 'minT_relPolewarness_GS_less_1Mkm2')


plot.gam(minT_relPolewarness_GS_less_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#select only ranges larger than 100000 km2
results_more_100000km2 <- results[results$rangeSize >= 100000,]

#select only ranges smaller than 100000 km2
results_less_100000km2 <- results[results$rangeSize <= 100000,]


#check number of occ and number os species
nrow(results_more_100000km2)
nrow(results_less_100000km2)

length(unique(results_more_100000km2$species))
length(unique(results_less_100000km2$species))


### minT (more 100,000)

#model G
minT_relPolewarness_G_more_100000km2 <- gam(Min_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_more_100000km2,
                                        method="REML",
                                        family="gaussian")


summary(minT_relPolewarness_G_more_100000km2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_G_more_100000km2, 'minT_relPolewarness_G_more_100000km2')


plot.gam(minT_relPolewarness_G_more_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPolewarness_GS_more_100000km2 <- gam(Min_T_SHAP ~
                                           s(relPolewardness, k=4, m=2) 
                                         + s(relPolewardness, species, k=4, bs="fs", m=2),
                                         data = results_more_1Mkm2,
                                         method="REML")


summary(minT_relPolewarness_GS_more_1Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_GS_more_1Mkm2, 'minT_relPolewarness_GS_more_1Mkm2')


plot.gam(minT_relPolewarness_GS_more_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### minT (less 100,000)

minT_relPolewarness_G_less_100000km2 <- gam(Min_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_less_100000km2,
                                        method="REML",
                                        family="gaussian")


summary(minT_relPolewarness_G_less_100000km2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_G_less_100000km2, 'minT_relPolewarness_G_less_100000km2')


plot.gam(minT_relPolewarness_G_less_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#select only ranges larger than 2M km2
results_more_2Mkm2 <- results[results$rangeSize >= 2000000,]

#select only ranges smaller than 2M km2
results_less_2Mkm2 <- results[results$rangeSize <= 2000000,]

#check number of occ and number os species
nrow(results_more_2Mkm2)
nrow(results_less_2Mkm2)

length(unique(results_more_2Mkm2$species))
length(unique(results_less_2Mkm2$species))

### minT (more 2,000,000)

#model G
minT_relPolewarness_G_more_2Mkm2 <- gam(Min_T_SHAP ~
                                              s(relPolewardness, k=4, bs="tp") 
                                            + s(species, k=503, bs="re"),
                                            data = results_more_2Mkm2,
                                            method="REML",
                                            family="gaussian")


summary(minT_relPolewarness_G_more_2Mkm2)

setwd(wd_models_res)
saveRDS(minT_relPolewarness_G_more_2Mkm2, 'minT_relPolewarness_G_more_2Mkm2')


plot.gam(minT_relPolewarness_G_more_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




### minT (less 100,000)

minT_relPolewarness_G_less_2Mkm2 <- gam(Min_T_SHAP ~
                                              s(relPolewardness, k=4, bs="tp") 
                                            + s(species, k=503, bs="re"),
                                            data = results_less_2Mkm2,
                                            method="REML",
                                            family="gaussian")


summary(minT_relPolewarness_G_less_2Mkm2)

minsetwd(wd_models_res)
saveRDS(minT_relPolewarness_G_less_2Mkm2, 'minT_relPolewarness_G_less_2Mkm2')


plot.gam(minT_relPolewarness_G_less_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



################################################################################

#                                MEAN 

################################################################################


### meanT (more 1,000,000)

#model G
meanT_relPolewarness_G_more_1Mkm2 <- gam(Mean_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_more_1Mkm2,
                                        method="REML",
                                        family="gaussian")


summary(meanT_relPolewarness_G_more_1Mkm2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_more_1Mkm2, 'meanT_relPolewarness_G_more_1Mkm2')


plot.gam(meanT_relPolewarness_G_more_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



### meanT (less 1,000,000)

meanT_relPolewarness_G_less_1Mkm2 <- gam(Mean_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_less_1Mkm2,
                                        method="REML",
                                        family="gaussian")


summary(meanT_relPolewarness_G_less_1Mkm2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_less_1Mkm2, 'meanT_relPolewarness_G_less_1Mkm2')


plot.gam(meanT_relPolewarness_G_less_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



### meanT (more 100,000)

#model G
meanT_relPolewarness_G_more_100000km2 <- gam(Mean_T_SHAP ~
                                              s(relPolewardness, k=4, bs="tp") 
                                            + s(species, k=503, bs="re"),
                                            data = results_more_100000km2,
                                            method="REML",
                                            family="gaussian")


summary(meanT_relPolewarness_G_more_100000km2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_more_100000km2, 'meanT_relPolewarness_G_more_100000km2')


plot.gam(meanT_relPolewarness_G_more_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




### meanT (less 100,000)

meanT_relPolewarness_G_less_100000km2 <- gam(Mean_T_SHAP ~
                                              s(relPolewardness, k=4, bs="tp") 
                                            + s(species, k=503, bs="re"),
                                            data = results_less_100000km2,
                                            method="REML",
                                            family="gaussian")


summary(meanT_relPolewarness_G_less_100000km2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_less_100000km2, 'meanT_relPolewarness_G_less_100000km2')


plot.gam(meanT_relPolewarness_G_less_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


### meanT (more 2,000,000)

#model G
meanT_relPolewarness_G_more_2Mkm2 <- gam(Mean_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_more_2Mkm2,
                                        method="REML",
                                        family="gaussian")


summary(meanT_relPolewarness_G_more_2Mkm2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_more_2Mkm2, 'meanT_relPolewarness_G_more_2Mkm2')


plot.gam(meanT_relPolewarness_G_more_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




### meanT (less 100,000)

meanT_relPolewarness_G_less_2Mkm2 <- gam(Mean_T_SHAP ~
                                          s(relPolewardness, k=4, bs="tp") 
                                        + s(species, k=503, bs="re"),
                                        data = results_less_2Mkm2,
                                        method="REML",
                                        family="gaussian")


summary(meanT_relPolewarness_G_less_2Mkm2)

setwd(wd_models_res)
saveRDS(meanT_relPolewarness_G_less_2Mkm2, 'meanT_relPolewarness_G_less_2Mkm2')


plot.gam(meanT_relPolewarness_G_less_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



################################################################################

#                                MAX 

################################################################################


### maxT (more 1,000,000)

#model G
maxT_relPolewarness_G_more_1Mkm2 <- gam(Max_T_SHAP ~
                                           s(relPolewardness, k=4, bs="tp") 
                                         + s(species, k=503, bs="re"),
                                         data = results_more_1Mkm2,
                                         method="REML",
                                         family="gaussian")


summary(maxT_relPolewarness_G_more_1Mkm2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_more_1Mkm2, 'maxT_relPolewarness_G_more_1Mkm2')


plot.gam(maxT_relPolewarness_G_more_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT (less 1,000,000)

maxT_relPolewarness_G_less_1Mkm2 <- gam(Max_T_SHAP ~
                                           s(relPolewardness, k=4, bs="tp") 
                                         + s(species, k=503, bs="re"),
                                         data = results_less_1Mkm2,
                                         method="REML",
                                         family="gaussian")


summary(maxT_relPolewarness_G_less_1Mkm2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_less_1Mkm2, 'maxT_relPolewarness_G_less_1Mkm2')


plot.gam(maxT_relPolewarness_G_less_1Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



### maxT (more 100,000)

#model G
maxT_relPolewarness_G_more_100000km2 <- gam(Max_T_SHAP ~
                                               s(relPolewardness, k=4, bs="tp") 
                                             + s(species, k=503, bs="re"),
                                             data = results_more_100000km2,
                                             method="REML",
                                             family="gaussian")


summary(maxT_relPolewarness_G_more_100000km2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_more_100000km2, 'maxT_relPolewarness_G_more_100000km2')


plot.gam(maxT_relPolewarness_G_more_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




### maxT (less 100,000)

maxT_relPolewarness_G_less_100000km2 <- gam(Max_T_SHAP ~
                                               s(relPolewardness, k=4, bs="tp") 
                                             + s(species, k=503, bs="re"),
                                             data = results_less_100000km2,
                                             method="REML",
                                             family="gaussian")


summary(maxT_relPolewarness_G_less_100000km2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_less_100000km2, 'maxT_relPolewarness_G_less_100000km2')


plot.gam(maxT_relPolewarness_G_less_100000km2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


### maxT (more 2,000,000)

#model G
maxT_relPolewarness_G_more_2Mkm2 <- gam(Max_T_SHAP ~
                                           s(relPolewardness, k=4, bs="tp") 
                                         + s(species, k=503, bs="re"),
                                         data = results_more_2Mkm2,
                                         method="REML",
                                         family="gaussian")


summary(maxT_relPolewarness_G_more_2Mkm2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_more_2Mkm2, 'maxT_relPolewarness_G_more_2Mkm2')


plot.gam(maxT_relPolewarness_G_more_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800




### maxT (less 100,000)

maxT_relPolewarness_G_less_2Mkm2 <- gam(Max_T_SHAP ~
                                           s(relPolewardness, k=4, bs="tp") 
                                         + s(species, k=503, bs="re"),
                                         data = results_less_2Mkm2,
                                         method="REML",
                                         family="gaussian")


summary(maxT_relPolewarness_G_less_2Mkm2)

setwd(wd_models_res)
saveRDS(maxT_relPolewarness_G_less_2Mkm2, 'maxT_relPolewarness_G_less_2Mkm2')


plot.gam(maxT_relPolewarness_G_less_2Mkm2, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800

