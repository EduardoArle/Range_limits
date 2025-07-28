#load libraries
library(mgcv); library(itsadug); library(gratia)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20250504_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250504_GAMs/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Significance_GAMs'

#read results table
setwd(wd_tables)
results <- read.csv('20250504_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)


########################
########## T ###########
########################

#select only species that had lower correl between vars
results_minT <- results[abs(results$Cor_vars_minT) <= 0.7,]
results_minT <- results_minT[complete.cases(results_minT$Cor_vars_minT),]

results_meanT <- results[abs(results$Cor_vars_meanT) <= 0.7,]
results_meanT <- results_meanT[complete.cases(results_meanT$Cor_vars_meanT),]

results_maxT <- results[abs(results$Cor_vars_maxT) <= 0.7,]
results_maxT <- results_maxT[complete.cases(results_maxT$Cor_vars_maxT),]

#make a version of the results with absolute shap values
results_minT_abs <- results_minT
results_meanT_abs <- results_meanT
results_maxT_abs <- results_maxT

results_minT_abs$avg_Min_T_SHAP <- abs(results_minT_abs$avg_Min_T_SHAP)
results_meanT_abs$avg_Mean_T_SHAP <- abs(results_meanT_abs$avg_Mean_T_SHAP)
results_maxT_abs$avg_Max_T_SHAP <- abs(results_maxT_abs$avg_Max_T_SHAP)

#make a version of the results keeping only points up to 250km from edges
results_minT_abs_250 <- results_minT_abs[results_minT_abs$distEdge <= 250,]
results_meanT_abs_250 <- results_meanT_abs[results_meanT_abs$distEdge <= 250,]
results_maxT_abs_250 <- results_maxT_abs[results_maxT_abs$distEdge <= 250,]

#count species in each table (to be used in the latitudinal analyses)
n_sps_minT <- length(unique(results_minT$species))
n_sps_meanT <- length(unique(results_meanT$species))
n_sps_maxT <- length(unique(results_maxT$species))

#count species in each modified table (to be used in the centre vs edge analyses)
n_sps_minT_mod <- length(unique(results_minT_abs_250$species))
n_sps_meanT_mod <- length(unique(results_meanT_abs_250$species))
n_sps_maxT_mod <- length(unique(results_maxT_abs_250$species))



# SHAP values X distance from edge (GAM)


### minT

#model G
minT_distEdge_G <- gam(avg_Min_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_minT, bs="re"),
                       data = results_minT_abs,
                       method="REML",
                       family="gaussian")


summary(minT_distEdge_G)
gam.check(minT_distEdge_G)


#AIC
AIC_minT_distEdge_G <- AIC(minT_distEdge_G)

setwd(wd_models)
saveRDS(minT_distEdge_G, 'minT_distEdge_G')


plot.gam(minT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250)
minT_distEdge_250_G <- gam(avg_Min_T_SHAP ~
                           s(distEdge, k=4, bs="tp") 
                         + s(species, k = n_sps_minT_mod, bs="re"),
                           data = results_minT_abs_250,
                           method="REML",
                           family="gaussian")


summary(minT_distEdge_250_G)
gam.check(minT_distEdge_250_G)


#AIC
AIC_minT_distEdge_250_G <- AIC(minT_distEdge_250_G)

setwd(wd_models)
saveRDS(minT_distEdge_250_G, 'minT_distEdge_250_G')


plot.gam(minT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_distEdge_GS <- gam(avg_Min_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_minT_abs,
                        method="REML")


summary(minT_distEdge_GS)
gam.check(minT_distEdge_GS)

#AIC
AIC_minT_distEdge_GS <- AIC(minT_distEdge_GS)

setwd(wd_models)
saveRDS(minT_distEdge_GS, 'minT_distEdge_GS')


plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
minT_distEdge_250_GS <- gam(avg_Min_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_minT_abs_250,
                        method="REML")


summary(minT_distEdge_250_GS)
gam.check(minT_distEdge_250_GS)

#AIC
AIC_minT_distEdge_250_GS <- AIC(minT_distEdge_250_GS)

setwd(wd_models)
saveRDS(minT_distEdge_250_GS, 'minT_distEdge_250_GS')


plot.gam(minT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_distEdge_250_GS.csv', row.names = F)



### meanT

#model G
meanT_distEdge_G <- gam(avg_Mean_T_SHAP ~
                          s(distEdge, k=4, bs="tp") 
                        + s(species, k = n_sps_meanT, bs="re"),
                        data = results_meanT_abs,
                        method="REML",
                        family="gaussian")


summary(meanT_distEdge_G)

#AIC
AIC_meanT_distEdge_G <- AIC(meanT_distEdge_G)

setwd(wd_models)
saveRDS(meanT_distEdge_G, 'meanT_distEdge_G')


plot.gam(meanT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250)
meanT_distEdge_250_G <- gam(avg_Mean_T_SHAP ~
                             s(distEdge, k=4, bs="tp") 
                           + s(species, k = n_sps_meanT_mod, bs="re"),
                           data = results_meanT_abs_250,
                           method="REML",
                           family="gaussian")


summary(meanT_distEdge_250_G)
gam.check(meanT_distEdge_G)

#AIC
AIC_meanT_distEdge_250_G <- AIC(meanT_distEdge_250_G)

setwd(wd_models)
saveRDS(meanT_distEdge_250_G, 'meanT_distEdge_250_G')


plot.gam(meanT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS 
meanT_distEdge_GS <- gam(avg_Mean_T_SHAP ~
                         s(distEdge, k=4, m=2) 
                       + s(distEdge, species, k=4, bs="fs", m=2),
                         data = results_meanT_abs,
                         method="REML")


summary(meanT_distEdge_GS)
gam.check(meanT_distEdge_GS)

#AIC
AIC_meanT_distEdge_GS <- AIC(meanT_distEdge_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_GS, 'meanT_distEdge_GS')


plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
meanT_distEdge_250_GS <- gam(avg_Mean_T_SHAP ~
                           s(distEdge, k=4, m=2) 
                         + s(distEdge, species, k=4, bs="fs", m=2),
                         data = results_meanT_abs_250,
                         method="REML")


summary(meanT_distEdge_250_GS)
gam.check(meanT_distEdge_250_GS)

#AIC
AIC_meanT_distEdge_250_GS <- AIC(meanT_distEdge_250_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_250_GS, 'meanT_distEdge_250_GS')


plot.gam(meanT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(meanT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_distEdge_250_GS.csv', row.names = F)


### maxT

#model G
maxT_distEdge_G <- gam(avg_Max_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT, bs="re"),
                       data = results_maxT_abs,
                       method="REML",
                       family="gaussian")


summary(maxT_distEdge_G)
gam.check(maxT_distEdge_G)

#AIC
AIC_maxT_distEdge_G <- AIC(maxT_distEdge_G)

setwd(wd_models)
saveRDS(maxT_distEdge_G, 'maxT_distEdge_G')


plot.gam(maxT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250)
maxT_distEdge_250_G <- gam(avg_Max_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT_mod, bs="re"),
                       data = results_maxT_abs_250,
                       method="REML",
                       family="gaussian")


summary(maxT_distEdge_250_G)
gam.check(maxT_distEdge_250_G)

#AIC
AIC_maxT_distEdge_250_G <- AIC(maxT_distEdge_250_G)

setwd(wd_models)
saveRDS(maxT_distEdge_250_G, 'maxT_distEdge_250_G')


plot.gam(maxT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_distEdge_GS <- gam(avg_Max_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxT_abs,
                        method="REML")


summary(maxT_distEdge_GS)
gam.check(maxT_distEdge_GS)

#AIC
AIC_maxT_distEdge_GS <- AIC(maxT_distEdge_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_GS, 'maxT_distEdge_GS')


plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
maxT_distEdge_250_GS <- gam(avg_Max_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxT_abs_250,
                        method="REML")


summary(maxT_distEdge_250_GS)
gam.check(maxT_distEdge_250_GS)

#AIC
AIC_maxT_distEdge_250_GS <- AIC(maxT_distEdge_250_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_250_GS, 'maxT_distEdge_250_GS')


plot.gam(maxT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_distEdge_250_GS.csv', row.names = F)



# SHAP values X relative polewardness (GAM)

### minT

#model G
minT_relPol_G <- gam(avg_Min_T_SHAP ~
                         s(relPolewardness, k=4, bs="tp") 
                       + s(species, k = n_sps_minT, bs="re"),
                       data = results_minT,
                       method="REML",
                       family="gaussian")


summary(minT_relPol_G)
gam.check(minT_relPol_G)

#AIC
AIC_minT_relPol_G <- AIC(minT_relPol_G)

setwd(wd_models)
saveRDS(minT_relPol_G, 'minT_relPol_G')


plot.gam(minT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPol_GS <- gam(avg_Min_T_SHAP ~
                          s(relPolewardness, k=4, m=2) 
                        + s(relPolewardness, species, k=4, bs="fs", m=2),
                        data = results_minT,
                        method="REML")


summary(minT_relPol_GS)
gam.check(minT_relPol_GS)

#AIC
AIC_minT_relPol_GS <- AIC(minT_relPol_GS)

setwd(wd_models)
saveRDS(minT_relPol_GS, 'minT_relPol_GS')


plot.gam(minT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#test stats for significant change points from First Derivative of GAM smooth

#calculate derivatives
deriv <- derivatives(minT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_relPol_GS.csv', row.names = F)

### meanT

#model G
meanT_relPol_G <- gam(avg_Mean_T_SHAP ~
                          s(relPolewardness, k=4, bs="tp") 
                        + s(species, k = n_sps_meanT, bs="re"),
                        data = results_meanT,
                        method="REML",
                        family="gaussian")


summary(meanT_relPol_G)
gam.check(meanT_relPol_G)

#AIC
AIC_meanT_relPol_G <- AIC(meanT_relPol_G)

setwd(wd_models)
saveRDS(meanT_relPol_G, 'meanT_relPol_G')


plot.gam(meanT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_relPol_GS <- gam(avg_Mean_T_SHAP ~
                           s(relPolewardness, k=4, m=2) 
                         + s(relPolewardness, species, k=4, bs="fs", m=2),
                         data = results_meanT,
                         method="REML")


summary(meanT_relPol_GS)
gam.check(meanT_relPol_GS)

#AIC
AIC_meanT_relPol_GS <- AIC(meanT_relPol_GS)

setwd(wd_models)
saveRDS(meanT_relPol_GS, 'meanT_relPol_GS')


plot.gam(meanT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(meanT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_relPol_GS.csv', row.names = F)


### maxT

#model G
maxT_relPol_G <- gam(avg_Max_T_SHAP ~
                         s(relPolewardness, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT, bs="re"),
                       data = results_maxT,
                       method="REML",
                       family="gaussian")


summary(maxT_relPol_G)
gam.check(maxT_relPol_G)

#AIC
AIC_maxT_relPol_G <- AIC(maxT_relPol_G)

setwd(wd_models)
saveRDS(maxT_relPol_G, 'maxT_relPol_G')


plot.gam(maxT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_relPol_GS <- gam(avg_Max_T_SHAP ~
                          s(relPolewardness, k=4, m=2) 
                        + s(relPolewardness, species, k=4, bs="fs", m=2),
                        data = results_maxT,
                        method="REML")


summary(maxT_relPol_GS)
gam.check(maxT_relPol_GS)

#AIC
AIC_maxT_relPol_GS <- AIC(maxT_relPol_GS)

setwd(wd_models)
saveRDS(maxT_relPol_GS, 'maxT_relPol_GS')


plot.gam(maxT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(maxT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_relPol_GS.csv', row.names = F)


# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_absPol_G <- gam(avg_Min_T_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_minT, bs="re"),
                     data = results_minT,
                     method="REML",
                     family="gaussian")


summary(minT_absPol_G)
gam.check(minT_absPol_G)

#AIC
AIC_minT_absPol_G <- AIC(minT_absPol_G)

setwd(wd_models)
saveRDS(minT_absPol_G, 'minT_absPol_G')


plot.gam(minT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_absPol_GS <- gam(avg_Min_T_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_minT,
                      method="REML")


summary(minT_absPol_GS)
gam.check(minT_absPol_GS)

#AIC
AIC_minT_absPol_GS <- AIC(minT_absPol_GS)

setwd(wd_models)
saveRDS(minT_absPol_GS, 'minT_absPol_GS')


plot.gam(minT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_absPol_GS.csv', row.names = F)




### meanT

#model G
meanT_absPol_G <- gam(avg_Mean_T_SHAP ~
                        s(absPolewardness, k=4, bs="tp") 
                      + s(species, k = n_sps_meanT, bs="re"),
                      data = results_meanT,
                      method="REML",
                      family="gaussian")


summary(meanT_absPol_G)
gam.check(meanT_absPol_G)

#AIC
AIC_meanT_absPol_G <- AIC(meanT_absPol_G)

setwd(wd_models)
saveRDS(meanT_absPol_G, 'meanT_absPol_G')


plot.gam(meanT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_absPol_GS <- gam(avg_Mean_T_SHAP ~
                         s(absPolewardness, k=4, m=2) 
                       + s(absPolewardness, species, k=4, bs="fs", m=2),
                       data = results_meanT,
                       method="REML")


summary(meanT_absPol_GS)
gam.check(meanT_absPol_GS)

#AIC
AIC_meanT_absPol_GS <- AIC(meanT_absPol_GS)

setwd(wd_models)
saveRDS(meanT_absPol_GS, 'meanT_absPol_GS')


plot.gam(meanT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.5, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(meanT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_absPol_GS.csv', row.names = F)



### maxT

#model G
maxT_absPol_G <- gam(avg_Max_T_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_maxT, bs="re"),
                     data = results_maxT,
                     method="REML",
                     family="gaussian")


summary(maxT_absPol_G)
gam.check(maxT_absPol_G)

#AIC
AIC_maxT_absPol_G <- AIC(maxT_absPol_G)

setwd(wd_models)
saveRDS(maxT_absPol_G, 'maxT_absPol_G')


plot.gam(maxT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_absPol_GS <- gam(avg_Max_T_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_maxT,
                      method="REML")


summary(maxT_absPol_GS)
gam.check(maxT_absPol_GS)

#AIC
AIC_maxT_absPol_GS <- AIC(maxT_absPol_GS)

setwd(wd_models)
saveRDS(maxT_absPol_GS, 'maxT_absPol_GS')


plot.gam(maxT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_absPol_GS.csv', row.names = F)




