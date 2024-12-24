#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20241210_GAMs/Models'

#read results table
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

###########.  GAM   ################

########################
######### PPT ##########
########################


# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPPT

#model G
minPPT_relPolewarness_G <- gam(avg_Min_PPT_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")

summary(minPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(minPPT_relPolewarness_G, 'minPPT_relPolewarness_G')

plot.gam(minPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_relPolewarness_GS <- gam(avg_Min_PPT_SHAP ~
                                s(relPolewardness, k=4, m=2) 
                              + s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(minPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(minPPT_relPolewarness_GS, 'minPPT_relPolewarness_GS')

plot.gam(minPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_relPolewarness_G <- gam(avg_Mean_PPT_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")

summary(meanPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(meanPPT_relPolewarness_G, 'meanPPT_relPolewarness_G')

plot.gam(meanPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_relPolewarness_GS <- gam(avg_Mean_PPT_SHAP ~
                                s(relPolewardness, k=4, m=2) 
                              + s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(meanPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(meanPPT_relPolewarness_GS, 'meanPPT_relPolewarness_GS')

plot.gam(meanPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxPPT

#model G
maxPPT_relPolewarness_G <- gam(avg_Max_PPT_SHAP ~
                                s(relPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")

summary(maxPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(maxPPT_relPolewarness_G, 'maxPPT_relPolewarness_G')


plot.gam(maxPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxPPT_relPolewarness_GS <- gam(avg_Max_PPT_SHAP ~
                                 s(relPolewardness, k=4, m=2) 
                               + s(relPolewardness, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")


summary(maxPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(maxPPT_relPolewarness_GS, 'maxPPT_relPolewarness_GS')


plot.gam(maxPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800



# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPPT

#model G
minPPT_absPolewarness_G <- gam(avg_Min_PPT_SHAP ~
                               s(absPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(minPPT_absPolewarness_G, 'minPPT_absPolewarness_G')


plot.gam(minPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_absPolewarness_GS <- gam(avg_Min_PPT_SHAP ~
                                s(absPolewardness, k=4, m=2) 
                              + s(absPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(minPPT_absPolewarness_GS)
 
setwd(wd_models)
saveRDS(minPPT_absPolewarness_GS, 'minPPT_absPolewarness_GS')


plot.gam(minPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800



### meanPPT

#model G
meanPPT_absPolewarness_G <- gam(avg_Mean_PPT_SHAP ~
                                s(absPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(meanPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(meanPPT_absPolewarness_G, 'meanPPT_absPolewarness_G')


plot.gam(meanPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.25, 0.12),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
meanPPT_absPolewarness_GS <- gam(avg_Mean_PPT_SHAP ~
                                 s(absPolewardness, k=4, m=2) 
                               + s(absPolewardness, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")

summary(meanPPT_absPolewarness_GS)

setwd(wd_models)
saveRDS(meanPPT_absPolewarness_GS, 'meanPPT_absPolewarness_GS')

plot.gam(meanPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.25, 0.12),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxPPT

#model G
maxPPT_absPolewarness_G <- gam(avg_Max_PPT_SHAP ~
                               s(absPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(maxPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(maxPPT_absPolewarness_G, 'maxPPT_absPolewarness_G')


plot.gam(maxPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.13, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxPPT_absPolewarness_GS <- gam(avg_Max_PPT_SHAP ~
                                s(absPolewardness, k=4, m=2) 
                              + s(absPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(maxPPT_absPolewarness_GS)

setwd(wd_models)
saveRDS(maxPPT_absPolewarness_GS, 'maxPPT_absPolewarness_GS')


plot.gam(maxPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.13, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800


# SHAP values X distance from edge (GAM)


### minPPT

#model G
minPPT_distEdge_G <- gam(avg_Min_PPT_SHAP ~
                               s(distEdge, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minPPT_distEdge_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_G, 'minPPT_distEdge_G')


plot.gam(minPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.05, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_distEdge_GS <- gam(avg_Min_PPT_SHAP ~
                                s(distEdge, k=4, m=2) 
                              + s(distEdge, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(minPPT_distEdge_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_GS, 'minPPT_distEdge_GS')


plot.gam(minPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.05, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_distEdge_G <- gam(avg_Mean_PPT_SHAP ~
                                s(distEdge, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(meanPPT_distEdge_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_G, 'meanPPT_distEdge_G')

plot.gam(meanPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
meanPPT_distEdge_GS <- gam(avg_Mean_PPT_SHAP ~
                                 s(distEdge, k=4, m=2) 
                               + s(distEdge, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")


summary(meanPPT_distEdge_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_GS, 'meanPPT_distEdge_GS')


plot.gam(meanPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxPPT

#model G
maxPPT_distEdge_G <- gam(avg_Max_PPT_SHAP ~
                               s(distEdge, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(maxPPT_distEdge_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_G, 'maxPPT_distEdge_G')

plot.gam(maxPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxPPT_distEdge_GS <- gam(avg_Max_PPT_SHAP ~
                                s(distEdge, k=4, m=2) 
                              + s(distEdge, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(maxPPT_distEdge_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_GS, 'maxPPT_distEdge_GS')


plot.gam(maxPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800










########################
######### TEMP #########
########################


# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_relPolewarness_G <- gam(avg_Min_T_SHAP ~
                                 s(relPolewardness, k=4, bs="tp") 
                               + s(species, k=503, bs="re"),
                               data = results,
                               method="REML",
                               family="gaussian")

summary(minT_relPolewarness_G)

setwd(wd_models)
saveRDS(minT_relPolewarness_G, 'minT_relPolewarness_G')

plot.gam(minT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPolewarness_GS <- gam(avg_Min_T_SHAP ~
                                  s(relPolewardness, k=4, m=2) 
                                + s(relPolewardness, species, k=4, bs="fs", m=2),
                                data = results,
                                method="REML")

summary(minT_relPolewarness_GS)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS, 'minT_relPolewarness_GS')

plot.gam(minT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT

#model G
meanT_relPolewarness_G <- gam(avg_Mean_T_SHAP ~
                                  s(relPolewardness, k=4, bs="tp") 
                                + s(species, k=503, bs="re"),
                                data = results,
                                method="REML",
                                family="gaussian")

summary(meanT_relPolewarness_G)

setwd(wd_models)
saveRDS(meanT_relPolewarness_G, 'meanT_relPolewarness_G')

plot.gam(meanT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_relPolewarness_GS <- gam(avg_Mean_T_SHAP ~
                                   s(relPolewardness, k=4, m=2) 
                                 + s(relPolewardness, species, k=4, bs="fs", m=2),
                                 data = results,
                                 method="REML")

summary(meanT_relPolewarness_GS)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS, 'meanT_relPolewarness_GS')

plot.gam(meanT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxT

#model G
maxT_relPolewarness_G <- gam(avg_Max_T_SHAP ~
                                 s(relPolewardness, k=4, bs="tp") 
                               + s(species, k=503, bs="re"),
                               data = results,
                               method="REML",
                               family="gaussian")

summary(maxT_relPolewarness_G)

setwd(wd_models)
saveRDS(maxT_relPolewarness_G, 'maxT_relPolewarness_G')


plot.gam(maxT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxT_relPolewarness_GS <- gam(avg_Max_T_SHAP ~
                                  s(relPolewardness, k=4, m=2) 
                                + s(relPolewardness, species, k=4, bs="fs", m=2),
                                data = results,
                                method="REML")


summary(maxT_relPolewarness_GS)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS, 'maxT_relPolewarness_GS')


plot.gam(maxT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800



# SHAP values X absolute polewardness (GAM)


### minT

#model G
minT_absPolewarness_G <- gam(avg_Min_T_SHAP ~
                                 s(absPolewardness, k=4, bs="tp") 
                               + s(species, k=503, bs="re"),
                               data = results,
                               method="REML",
                               family="gaussian")


summary(minT_absPolewarness_G)

setwd(wd_models)
saveRDS(minT_absPolewarness_G, 'minT_absPolewarness_G')


plot.gam(minT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.11, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_absPolewarness_GS <- gam(avg_Min_T_SHAP ~
                                  s(absPolewardness, k=4, m=2) 
                                + s(absPolewardness, species, k=4, bs="fs", m=2),
                                data = results,
                                method="REML")


summary(minT_absPolewarness_GS)

setwd(wd_models)
saveRDS(minT_absPolewarness_GS, 'minT_absPolewarness_GS')


plot.gam(minT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.11, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800



### meanT

#model G
meanT_absPolewarness_G <- gam(avg_Mean_T_SHAP ~
                                  s(absPolewardness, k=4, bs="tp") 
                                + s(species, k=503, bs="re"),
                                data = results,
                                method="REML",
                                family="gaussian")


summary(meanT_absPolewarness_G)

setwd(wd_models)
saveRDS(meanT_absPolewarness_G, 'meanT_absPolewarness_G')


plot.gam(meanT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.12, 0.12),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
meanT_absPolewarness_GS <- gam(avg_Mean_T_SHAP ~
                                   s(absPolewardness, k=4, m=2) 
                                 + s(absPolewardness, species, k=4, bs="fs", m=2),
                                 data = results,
                                 method="REML")

summary(meanT_absPolewarness_GS)

setwd(wd_models)
saveRDS(meanT_absPolewarness_GS, 'meanT_absPolewarness_GS')

plot.gam(meanT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.12, 0.12),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxT

#model G
maxT_absPolewarness_G <- gam(avg_Max_T_SHAP ~
                                 s(absPolewardness, k=4, bs="tp") 
                               + s(species, k=503, bs="re"),
                               data = results,
                               method="REML",
                               family="gaussian")


summary(maxT_absPolewarness_G)

setwd(wd_models)
saveRDS(maxT_absPolewarness_G, 'maxT_absPolewarness_G')


plot.gam(maxT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxT_absPolewarness_GS <- gam(avg_Max_T_SHAP ~
                                  s(absPolewardness, k=4, m=2) 
                                + s(absPolewardness, species, k=4, bs="fs", m=2),
                                data = results,
                                method="REML")


summary(maxT_absPolewarness_GS)

setwd(wd_models)
saveRDS(maxT_absPolewarness_GS, 'maxT_absPolewarness_GS')


plot.gam(maxT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


# SHAP values X distance from edge (GAM)


### minT

#model G
minT_distEdge_G <- gam(avg_Min_T_SHAP ~
                           s(distEdge, k=4, bs="tp") 
                         + s(species, k=503, bs="re"),
                         data = results,
                         method="REML",
                         family="gaussian")


summary(minT_distEdge_G)

setwd(wd_models)
saveRDS(minT_distEdge_G, 'minT_distEdge_G')


plot.gam(minT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_distEdge_GS <- gam(avg_Min_T_SHAP ~
                            s(distEdge, k=4, m=2) 
                          + s(distEdge, species, k=4, bs="fs", m=2),
                          data = results,
                          method="REML")


summary(minT_distEdge_GS)

setwd(wd_models)
saveRDS(minT_distEdge_GS, 'minT_distEdge_GS')


plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT

#model G
meanT_distEdge_G <- gam(avg_Mean_T_SHAP ~
                            s(distEdge, k=4, bs="tp") 
                          + s(species, k=503, bs="re"),
                          data = results,
                          method="REML",
                          family="gaussian")


summary(meanT_distEdge_G)

setwd(wd_models)
saveRDS(meanT_distEdge_G, 'meanT_distEdge_G')

plot.gam(meanT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
meanT_distEdge_GS <- gam(avg_Mean_T_SHAP ~
                             s(distEdge, k=4, m=2) 
                           + s(distEdge, species, k=4, bs="fs", m=2),
                           data = results,
                           method="REML")


summary(meanT_distEdge_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_GS, 'meanT_distEdge_GS')


plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxT

#model G
maxT_distEdge_G <- gam(avg_Max_T_SHAP ~
                           s(distEdge, k=4, bs="tp") 
                         + s(species, k=503, bs="re"),
                         data = results,
                         method="REML",
                         family="gaussian")


summary(maxT_distEdge_G)

setwd(wd_models)
saveRDS(maxT_distEdge_G, 'maxT_distEdge_G')

plot.gam(maxT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.13, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxT_distEdge_GS <- gam(avg_Max_T_SHAP ~
                            s(distEdge, k=4, m=2) 
                          + s(distEdge, species, k=4, bs="fs", m=2),
                          data = results,
                          method="REML")


summary(maxT_distEdge_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_GS, 'maxT_distEdge_GS')


plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.13, 0.02),
         cex.lab = 2, cex.axis = 1.5) #save 800



