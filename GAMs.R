#load libraries
library(mgcv); library(mgcViz)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'

#read temperature results table

setwd(wd_tables)
results <- read.csv('Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,11)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

###########.  GAM   ################

# SHAP values X polewardness (GAM)

# rel polewardness


### minT
gam_minT_relPolewarness <- gam(Min_T_SHAP ~ s(relPolewardness, k = 4),
                               data = results)

summary(gam_minT_relPolewarness)

plot.gam(gam_minT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)

#including random effects
gam_minT_relPolewarness_RE <- gam(Min_T_SHAP ~ s(relPolewardness, k = 4) +
                                                  s(species, bs = "re"),
                                    data = results)

summary(gam_minT_relPolewarness_RE)

plot.gam(gam_minT_relPolewarness_RE, pages = 2, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)



### meanT
gam_meanT_relPolewarness <- gam(Mean_T_SHAP ~ s(relPolewardness, k = 4),
                                na.action = 'na.omit',
                                data = results)

summary(gam_meanT_relPolewarness)


plot.gam(gam_meanT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)

#including random effects
gam_meanT_relPolewarness_RE <- gam(Mean_T_SHAP ~ s(relPolewardness, k = 4) +
                                    s(species, bs = "re"),
                                  data = results)

summary(gam_meanT_relPolewarness_RE)

plot.gam(gam_meanT_relPolewarness_RE, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)


### maxT
gam_maxT_relPolewarness <- gam(Max_T_SHAP ~ s(relPolewardness, k = 4),
                               na.action = 'na.omit',
                               data = results)

summary(gam_maxT_relPolewarness)

plot.gam(gam_maxT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)

#including random effects
gam_maxT_relPolewarness_RE <- gam(Max_T_SHAP ~ s(relPolewardness, k = 4) +
                                     s(species, bs = "re"),
                                   data = results)

summary(gam_maxT_relPolewarness_RE)

plot.gam(gam_maxT_relPolewarness_RE, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)


# abs latitude


### minT
gam_minT_absPolewarness <- gam(Min_T_SHAP ~ s(abs(decimalLatitude), k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_minT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value', xlab = 'Latitude',
         ylim = c(-0.05, 0.05),
         cex.lab = 2, cex.axis = 2)



### meanT
gam_meanT_absPolewarness <- gam(Mean_T_SHAP ~ s(abs(decimalLatitude), k = 4),
                                na.action = 'na.omit',
                                data = results)

plot.gam(gam_meanT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value', xlab = 'Latitude',
         ylim = c(-0.05, 0.05),
         cex.lab = 2, cex.axis = 2)

### maxT
gam_maxT_absPolewarness <- gam(Max_T_SHAP ~ s(abs(decimalLatitude), k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_maxT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value', xlab = 'Latitude',
         ylim = c(-0.05, 0.05),
         cex.lab = 2, cex.axis = 2)



# dist edge 

### minT
gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge, k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_minT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value')



### meanT
gam_meanT_distEdge <- gam(Mean_T_SHAP ~ s(distEdge, k = 4),
                          na.action = 'na.omit',
                          data = results)

plot.gam(gam_meanT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value')

### maxT
gam_maxT_distEdge <- gam(Max_T_SHAP ~ s(distEdge, k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_maxT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value')


# dist edge (consider only points not too far from the edge)

results2 <- results[results$distEdge <= 100,]

### minT
gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge, k = 4),
                         na.action = 'na.omit',
                         data = results2)

plot.gam(gam_minT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value')



### meanT
gam_meanT_distEdge <- gam(Mean_T_SHAP ~ s(distEdge, k = 4),
                          na.action = 'na.omit',
                          data = results2)

plot.gam(gam_meanT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value')

### maxT
gam_maxT_distEdge <- gam(Max_T_SHAP ~ s(distEdge, k = 4),
                         na.action = 'na.omit',
                         data = results2)

plot.gam(gam_maxT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value')



# elevation

### minT
gam_minT_elev <- gam(Min_T_SHAP ~ s(elevation, k = 4),
                     na.action = 'na.omit',
                     data = results)

plot.gam(gam_minT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value')



### meanT
gam_meanT_elev <- gam(Mean_T_SHAP ~ s(elevation, k = 4),
                      na.action = 'na.omit',
                      data = results)

plot.gam(gam_meanT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value')

### maxT
gam_maxT_elev <- gam(Max_T_SHAP ~ s(elevation, k = 4),
                     na.action = 'na.omit',
                     data = results)

plot.gam(gam_maxT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value')




#####  PPT


# rel polewardness


### minPPT
gam_minPPT_relPolewarness <- gam(Min_PPT_SHAP ~ s(relPolewardness, k = 4),
                                 na.action = 'na.omit',
                                 data = results)

plot.gam(gam_minPPT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#65432130', ylab = 'SHAP value')



### meanPPT
gam_meanPPT_relPolewarness <- gam(Mean_PPT_SHAP ~ s(relPolewardness, k = 4),
                                na.action = 'na.omit',
                                data = results)

plot.gam(gam_meanPPT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#FF752D30', ylab = 'SHAP value')

### maxPPT
gam_maxPPT_relPolewarness <- gam(Max_PPT_SHAP ~ s(relPolewardness, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_maxPPT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value')


# abs polewardness


### minPPT
gam_minPPT_absPolewarness <- gam(Min_PPT_SHAP ~ s(absPolewardness, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_minPPT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#65432130', ylab = 'SHAP value')



### meanPPT
gam_meanPPT_absPolewarness <- gam(Mean_PPT_SHAP ~ s(absPolewardness, k = 4),
                                na.action = 'na.omit',
                                data = results)

plot.gam(gam_meanPPT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#FF752D30', ylab = 'SHAP value')

### maxPPT
gam_maxPPT_absPolewarness <- gam(Max_PPT_SHAP ~ s(absPolewardness, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_maxPPT_absPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value')



# dist edge

### minPPT
gam_minPPT_distEdge <- gam(Min_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_minPPT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#65432130', ylab = 'SHAP value')


### meanPPT
gam_meanPPT_distEdge <- gam(Mean_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
                          na.action = 'na.omit',
                          data = results)

plot.gam(gam_meanPPT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#FF752D30', ylab = 'SHAP value')

### maxPPT
gam_maxPPT_distEdge <- gam(Max_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_maxPPT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value')


# elevation

### minPPT
gam_minPPT_elev <- gam(Min_PPT_SHAP ~ s(elevation, k = 4),
                     na.action = 'na.omit',
                     data = results)

plot.gam(gam_minPPT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#65432130', ylab = 'SHAP value')



### meanPPT
gam_meanPPT_elev <- gam(Mean_PPT_SHAP ~ s(elevation, k = 4),
                      na.action = 'na.omit',
                      data = results)

plot.gam(gam_meanPPT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#FF752D30', ylab = 'SHAP value')

### maxPPT
gam_maxPPT_elev <- gam(Max_PPT_SHAP ~ s(elevation, k = 4),
                     na.action = 'na.omit',
                     data = results)

plot.gam(gam_maxPPT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value')



#########

install.packages("lme4")
library(lme4)
data("sleepstudy")

head(sleepstudy)

#example without random effects
gam_model <- gam(Reaction ~ s(Days, k = 5), data = sleepstudy)

summary(gam_model)

plot(gam_model, pages = 1, residuals = F, shade = T, shade.col = '#0000FF30')

#example with random effects
gam_model_re <- gam(Reaction ~ s(Days, k = 5) +
                      s(Subject, bs = "re"), data = sleepstudy)

summary(gam_model_re)

plot(gam_model_re, pages = 1, residuals = F, shade = T, shade.col = '#0000FF30')


# gam_model <- gam(Reaction ~ s(Days, k = 10), data = sleepstudy)
# summary(gam_model)
# 
# plot(gam_model, pages = 1)


###
plot.gam(gam_minT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)

