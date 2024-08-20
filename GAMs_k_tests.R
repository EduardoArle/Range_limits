#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models_k_tests'

#read temperature results table

setwd(wd_tables)
results <- read.csv('Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,11)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

###########.  GAM   ################

# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model GS k (first term) = 5 k (second term) = 4
minT_relPolewarness_GS_k5_k4 <- gam(Min_T_SHAP ~
                                      s(relPolewardness, k=5, m=2) 
                                    + s(relPolewardness, species, k=4, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(minT_relPolewarness_GS_k5_k4)
logLik(minT_relPolewarness_GS_k5_k4)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS_k5_k4, 'minT_relPolewarness_GS_k5_k4')


plot.gam(minT_relPolewarness_GS_k5_k4, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_relPolewarness_GS_k5_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS k (first term) = 5 k (second term) = 5
minT_relPolewarness_GS_k5_k5 <- gam(Min_T_SHAP ~
                                s(relPolewardness, k=5, m=2) 
                              + s(relPolewardness, species, k=5, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(minT_relPolewarness_GS_k5_k5)
logLik(minT_relPolewarness_GS_k5_k5)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS_k5_k5, 'minT_relPolewarness_GS_k5_k5')


plot(minT_relPolewarness_GS_k5_k5, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot(minT_relPolewarness_GS_k5_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 4
minT_relPolewarness_GS_k6_k4 <- gam(Min_T_SHAP ~
                                   s(relPolewardness, k=6, m=2) 
                                 + s(relPolewardness, species, k=4, bs="fs", m=2),
                                 data = results,
                                 method="REML")


summary(minT_relPolewarness_GS_k6_k4)
logLik(minT_relPolewarness_GS_k6_k4)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS_k6_k4, 'minT_relPolewarness_GS_k6_k4')

draw(minT_relPolewarness_GS_k6_k4, page = 1)


plot.gam(minT_relPolewarness_GS_k6_k4, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_relPolewarness_GS_k6_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 5
minT_relPolewarness_GS_k6_k5 <- gam(Min_T_SHAP ~
                                      s(relPolewardness, k=6, m=2) 
                                    + s(relPolewardness, species, k=5, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(minT_relPolewarness_GS_k6_k5)
logLik(minT_relPolewarness_GS_k6_k5)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS_k6_k5, 'minT_relPolewarness_GS_k6_k5')

draw(minT_relPolewarness_GS_k6_k5, page = 1)


plot.gam(minT_relPolewarness_GS_k6_k5, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_relPolewarness_GS_k6_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 6
minT_relPolewarness_GS_k6_k6 <- gam(Min_T_SHAP ~
                                      s(relPolewardness, k=6, m=2) 
                                    + s(relPolewardness, species, k=6, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(minT_relPolewarness_GS_k6_k6)
logLik(minT_relPolewarness_GS_k6_k6)

setwd(wd_models)
saveRDS(minT_relPolewarness_GS_k6_k6, 'minT_relPolewarness_GS_k6_k6')

draw(minT_relPolewarness_GS_k6_k6, page = 1)


plot.gam(minT_relPolewarness_GS_k6_k6, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_relPolewarness_GS_k6_k6, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


### meanT

#model GS k (first term) = 5 k (second term) = 4
meanT_relPolewarness_GS_k5_k4 <- gam(Mean_T_SHAP ~
                                      s(relPolewardness, k=5, m=2) 
                                    + s(relPolewardness, species, k=4, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(meanT_relPolewarness_GS_k5_k4)
logLik(meanT_relPolewarness_GS_k5_k4)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS_k5_k4, 'meanT_relPolewarness_GS_k5_k4')


plot.gam(meanT_relPolewarness_GS_k5_k4, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k5_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800




#model GS k (first term) = 5 k (second term) = 5
meanT_relPolewarness_GS_k5_k5 <- gam(Mean_T_SHAP ~
                                      s(relPolewardness, k=5, m=2) 
                                    + s(relPolewardness, species, k=5, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(meanT_relPolewarness_GS_k5_k5)
logLik(meanT_relPolewarness_GS_k5_k5)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS_k5_k5, 'meanT_relPolewarness_GS_k5_k5')


plot.gam(meanT_relPolewarness_GS_k5_k5, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k5_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 4
meanT_relPolewarness_GS_k6_k4 <- gam(Mean_T_SHAP ~
                                      s(relPolewardness, k=6, m=2) 
                                    + s(relPolewardness, species, k=4, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(meanT_relPolewarness_GS_k6_k4)
logLik(meanT_relPolewarness_GS_k6_k4)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS_k6_k4, 'meanT_relPolewarness_GS_k6_k4')

draw(minT_relPolewarness_GS_k6_k4, page = 1)


plot.gam(meanT_relPolewarness_GS_k6_k4, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k6_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 5
meanT_relPolewarness_GS_k6_k5 <- gam(Mean_T_SHAP ~
                                      s(relPolewardness, k=6, m=2) 
                                    + s(relPolewardness, species, k=5, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(meanT_relPolewarness_GS_k6_k5)
logLik(meanT_relPolewarness_GS_k6_k5)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS_k6_k5, 'meanT_relPolewarness_GS_k6_k5')

draw(minT_relPolewarness_GS_k6_k5, page = 1)


plot.gam(meanT_relPolewarness_GS_k6_k5, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k6_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 6
meanT_relPolewarness_GS_k6_k6 <- gam(Mean_T_SHAP ~
                                      s(relPolewardness, k=6, m=2) 
                                    + s(relPolewardness, species, k=6, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(meanT_relPolewarness_GS_k6_k6)
logLik(meanT_relPolewarness_GS_k6_k6)

setwd(wd_models)
saveRDS(meanT_relPolewarness_GS_k6_k6, 'meanT_relPolewarness_GS_k6_k6')


plot.gam(meanT_relPolewarness_GS_k6_k6, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k6_k6, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800



### maxT

#model GS k (first term) = 5 k (second term) = 4
maxT_relPolewarness_GS_k5_k4 <- gam(Max_T_SHAP ~
                                       s(relPolewardness, k=5, m=2) 
                                     + s(relPolewardness, species, k=4, bs="fs", m=2),
                                     data = results,
                                     method="REML")


summary(maxT_relPolewarness_GS_k5_k4)
logLik(maxT_relPolewarness_GS_k5_k4)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k5_k4, 'maxT_relPolewarness_GS_k5_k4')


plot.gam(maxT_relPolewarness_GS_k5_k4, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k5_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800




#model GS k (first term) = 5 k (second term) = 5
maxT_relPolewarness_GS_k5_k5 <- gam(Max_T_SHAP ~
                                       s(relPolewardness, k=5, m=2) 
                                     + s(relPolewardness, species, k=5, bs="fs", m=2),
                                     data = results,
                                     method="REML")


summary(maxT_relPolewarness_GS_k5_k5)
logLik(maxT_relPolewarness_GS_k5_k5)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k5_k5, 'maxT_relPolewarness_GS_k5_k5')


plot.gam(maxT_relPolewarness_GS_k5_k5, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k5_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 4
maxT_relPolewarness_GS_k6_k4 <- gam(Max_T_SHAP ~
                                       s(relPolewardness, k=6, m=2) 
                                     + s(relPolewardness, species, k=4, bs="fs", m=2),
                                     data = results,
                                     method="REML")


summary(maxT_relPolewarness_GS_k6_k4)
logLik(maxT_relPolewarness_GS_k6_k4)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k6_k4, 'maxT_relPolewarness_GS_k6_k4')


plot.gam(maxT_relPolewarness_GS_k6_k4, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k6_k4, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 5
maxT_relPolewarness_GS_k6_k5 <- gam(Max_T_SHAP ~
                                       s(relPolewardness, k=6, m=2) 
                                     + s(relPolewardness, species, k=5, bs="fs", m=2),
                                     data = results,
                                     method="REML")


summary(maxT_relPolewarness_GS_k6_k5)
logLik(maxT_relPolewarness_GS_k6_k5)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k6_k5, 'maxT_relPolewarness_GS_k6_k5')


plot.gam(maxT_relPolewarness_GS_k6_k5, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

maxT_relPol_GS_plot.gam(maxT_relPolewarness_GS_k6_k5, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 6 k (second term) = 6
maxT_relPolewarness_GS_k6_k6 <- gam(Max_T_SHAP ~
                                       s(relPolewardness, k=6, m=2) 
                                     + s(relPolewardness, species, k=6, bs="fs", m=2),
                                     data = results,
                                     method="REML")


summary(maxT_relPolewarness_GS_k6_k6)
logLik(maxT_relPolewarness_GS_k6_k6)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k6_k6, 'maxT_relPolewarness_GS_k6_k6')


plot.gam(maxT_relPolewarness_GS_k6_k6, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k6_k6, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model GS k (first term) = 10 k (second term) = 10
maxT_relPolewarness_GS_k10_k10 <- gam(Max_T_SHAP ~
                                      s(relPolewardness, k=10, m=2) 
                                    + s(relPolewardness, species, k=10, bs="fs", m=2),
                                    data = results,
                                    method="REML")


summary(maxT_relPolewarness_GS_k10_k10)
logLik(maxT_relPolewarness_GS_k10_k10)

setwd(wd_models)
saveRDS(maxT_relPolewarness_GS_k10_k10, 'maxT_relPolewarness_GS_k10_k10')


plot.gam(maxT_relPolewarness_GS_k6_k6, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k6_k6, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#################


#load and save plots of models run with k = 4, k = 4 with the same ylim

wd_models_k4_k4 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models'

setwd(wd_models_k4_k4)

minT_relPolewarness_GS_k4_k4 <- readRDS('minT_relPolewarness_GS')
meanT_relPolewarness_GS_k4_k4 <- readRDS('meanT_relPolewarness_GS')
maxT_relPolewarness_GS_k4_k4 <- readRDS('maxT_relPolewarness_GS')

logLik(minT_relPolewarness_GS_k4_k4)
logLik(meanT_relPolewarness_GS_k4_k4)
logLik(maxT_relPolewarness_GS_k4_k4)

plot.gam(minT_relPolewarness_GS_k4_k4, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_relPolewarness_GS_k4_k4, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_relPolewarness_GS_k4_k4, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.04),
         cex.lab = 2, cex.axis = 1.5) #save 800
