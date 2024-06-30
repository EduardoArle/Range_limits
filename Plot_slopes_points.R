#load libraries
library(mgcv)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'

#read temperature results table

setwd(wd_tables)
results <- read.csv('Results_all_sps.csv')

###### Plot the slopes per species against each variables in the RF

## Set graphical parametres

par(mar = c(5,5,4,4))

## Get slope values for y and x lims

#minT
cont_minT <- results$Min_T_SHAP

#meanT
cont_meanT <- results$Mean_T_SHAP

#maxT
cont_maxT <- results$Max_T_SHAP




##### ELEVATION ####

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(results$elevation) #no matter which

# SHAP values X elevation

#plot minT
plot(results$elevation, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Elevation',
     ylim = ylim)

lm_minT_elev <- lm(Min_T_SHAP ~ elevation,
                         data = results)

abline(lm_minT_elev, col = '#0000FF', lwd = 2)

#add meanT
points(results$elevation, results$Mean_T_SHAP,
       pch = 19, cex = 0.4, col = '#80008005')

lm_meanT_elev <- lm(Mean_T_SHAP ~ elevation,
                          data = results)

abline(lm_meanT_elev, col = '#800080', lwd = 2)

#add maxT
points(results$elevation, results$Max_T_SHAP,
       pch = 19, cex = 0.4, col = '#FF000005')

lm_maxT_elev <- lm(Max_T_SHAP ~ elevation,
                    data = results)

abline(lm_maxT_elev, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_elev)

r2_minT <- summa_minT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1], 
     expression(bold(paste('R'^2*' = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] - ((ylim[2] - ylim[1]) / 200), 
     round(r2_minT, 3), col = '#0000FF', font = 2, pos = 4)

p_value_minT <- summa_minT$coefficients[2,4]
if(p_value_minT < 0.0001){
  p_value_minT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#0000FF', pos = 4)
if(p_value_minT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
       p_value_minT, col = '#0000FF', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
       round(p_value_minT, 4), col = '#0000FF', font = 2, pos = 4)
}


summa_meanT <- summary(lm_meanT_elev)
r2_meanT <- summa_meanT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)

text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_meanT, 3), col = '#800080', font = 2, pos = 4)

p_value_meanT <- summa_meanT$coefficients[2,4]
if(p_value_meanT < 0.0001){
  p_value_meanT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#800080', pos = 4)
if(p_value_meanT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
       p_value_meanT, col = '#800080', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
       round(p_value_meanT, 4), col = '#800080', font = 2, pos = 4)
}


summa_maxT <- summary(lm_maxT_elev)
r2_maxT <- summa_maxT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_maxT, 3), col = '#FF0000', font = 2, pos = 4)


p_value_maxT <- summa_maxT$coefficients[2,4]
if(p_value_maxT < 0.0001){
  p_value_maxT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#FF0000', pos = 4)
if(p_value_maxT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - ((ylim[2] - ylim[1]) / 22), 
       p_value_maxT, col = '#FF0000', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - ((ylim[2] - ylim[1]) / 22), 
       round(p_value_maxT, 4), col = '#FF0000', font = 2, pos = 4)
}

#save in width 1000 in /Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240523/Plots_slopes/Slope_all_points

# SHAP values X elevation (GAM)

#plot minT
plot(results$elevation, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Elevation',
     ylim = ylim)

plot(results$relPolarwardness, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Elevation',
     ylim = ylim)

gam_minT_elev <- gam(Min_T_SHAP ~ elevation +
                       distEdge +
                       relPolarwardness +
                       decimalLatitude +
                       absPolarwardness +
                       biome +
                       species,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

summary(gam_minT_elev)

abline(gam_minT_elev, col = '#0000FF', lwd = 2)


#add meanT
points(results$elevation, results$Mean_T_SHAP,
       pch = 19, cex = 0.4, col = '#80008005')

gam_meanT_elev <- gam(Mean_T_SHAP ~ elevation,
                     na.action = 'na.omit',
                     data = results)

abline(gam_meanT_elev, col = '#800080', lwd = 2)


#add maxT
points(results$elevation, results$Max_T_SHAP,
       pch = 19, cex = 0.4, col = '#FF000005')

gam_maxT_elev <- gam(Max_T_SHAP ~ elevation,
                      na.action = 'na.omit',
                      data = results)

abline(gam_maxT_elev, col = '#800080', lwd = 2)


#add stats to the plot
summa_minT_gam <- summary(gam_minT_elev)

r2_minT <- summa_minT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1], 
     expression(bold(paste('R'^2*' = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] - ((ylim[2] - ylim[1]) / 200), 
     round(r2_minT, 3), col = '#0000FF', font = 2, pos = 4)

p_value_minT <- summa_minT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     round(p_value_minT, 4), col = '#0000FF', font = 2, pos = 4)

summa_meanT <- summary(lm_meanT_range_Size)
r2_meanT <- summa_meanT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_meanT, 3), col = '#800080', font = 2, pos = 4)

p_value_meanT <- summa_meanT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#800080', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     round(p_value_meanT, 4), col = '#800080', font = 2, pos = 4)

summa_maxT <- summary(lm_maxT_range_Size)
r2_maxT <- summa_maxT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_maxT, 3), col = '#FF0000', font = 2, pos = 4)


p_value_maxT <- summa_maxT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     round(p_value_maxT, 4), col = '#FF0000', font = 2, pos = 4)



##### DISTANCE FROM EDGE ####

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(results$distEdge) #no matter which

# SHAP values X dist edge

#plot minT
plot(results$distEdge, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Distance from edge',
     ylim = ylim)

lm_minT_distEdge <- lm(Min_T_SHAP ~ distEdge,
                   data = results)

abline(lm_minT_distEdge, col = '#0000FF', lwd = 2)

#add meanT
points(results$distEdge, results$Mean_T_SHAP,
       pch = 19, cex = 0.4, col = '#80008005')

lm_meanT_distEdge <- lm(Mean_T_SHAP ~ distEdge,
                    data = results)

abline(lm_meanT_distEdge, col = '#800080', lwd = 2)

#add maxT
points(results$distEdge, results$Max_T_SHAP,
       pch = 19, cex = 0.4, col = '#FF000005')

lm_maxT_distEdge <- lm(Max_T_SHAP ~ distEdge,
                   data = results)

abline(lm_maxT_distEdge, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_distEdge)

r2_minT <- summa_minT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1], 
     expression(bold(paste('R'^2*' = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] - ((ylim[2] - ylim[1]) / 200), 
     round(r2_minT, 3), col = '#0000FF', font = 2, pos = 4)

p_value_minT <- summa_minT$coefficients[2,4]
if(p_value_minT < 0.0001){
  p_value_minT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#0000FF', pos = 4)
if(p_value_minT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
       p_value_minT, col = '#0000FF', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
       round(p_value_minT, 4), col = '#0000FF', font = 2, pos = 4)
}


summa_meanT <- summary(lm_meanT_distEdge)
r2_meanT <- summa_meanT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)

text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_meanT, 3), col = '#800080', font = 2, pos = 4)

p_value_meanT <- summa_meanT$coefficients[2,4]
if(p_value_meanT < 0.0001){
  p_value_meanT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#800080', pos = 4)
if(p_value_meanT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
       p_value_meanT, col = '#800080', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
       round(p_value_meanT, 4), col = '#800080', font = 2, pos = 4)
}


summa_maxT <- summary(lm_maxT_distEdge)
r2_maxT <- summa_maxT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_maxT, 3), col = '#FF0000', font = 2, pos = 4)


p_value_maxT <- summa_maxT$coefficients[2,4]
if(p_value_maxT < 0.0001){
  p_value_maxT <- '< 0.0001'
}

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#FF0000', pos = 4)
if(p_value_maxT == '< 0.0001'){
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - ((ylim[2] - ylim[1]) / 22), 
       p_value_maxT, col = '#FF0000', font = 2, pos = 4)
}else{
  text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
       ylim[2] - ((ylim[2] - ylim[1]) / 22), 
       round(p_value_maxT, 4), col = '#FF0000', font = 2, pos = 4)
}

#save in width 1000 in /Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240523/Plots_SHAP_vars
###################

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(log(distEdge), k = 2),
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


#plot minT
plot(log(results$distEdge), results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_minT_distEdge))

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(log(distEdge), k = 2),
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


#plot meanT
plot(log(results$distEdge), results$Mean_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_meanT_distEdge))

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(log(distEdge), k = 2),
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


#plot maxT
plot(log(results$distEdge), results$Max_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_maxT_distEdge))



###########.  GAM   ################

# SHAP values X polewardness (GAM)

# rel polewardness


### minT
gam_minT_relPolewarness <- gam(Min_T_SHAP ~ s(relPolarwardness, k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_minT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30')



### meanT
gam_meanT_relPolewarness <- gam(Mean_T_SHAP ~ s(relPolarwardness, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_meanT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030')

### maxT
gam_maxT_relPolewarness <- gam(Max_T_SHAP ~ s(relPolarwardness, k = 4),
                                na.action = 'na.omit',
                                data = results)

plot.gam(gam_maxT_relPolewarness, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030')



# dist edge

### minT
gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_minT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30')



### meanT
gam_meanT_distEdge <- gam(Mean_T_SHAP ~ s(distEdge, k = 4),
                                na.action = 'na.omit',
                                data = results)

plot.gam(gam_meanT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030')

### maxT
gam_maxT_distEdge <- gam(Max_T_SHAP ~ s(distEdge, k = 4),
                               na.action = 'na.omit',
                               data = results)

plot.gam(gam_maxT_distEdge, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030')


# elevation

### minT
gam_minT_elev <- gam(Min_T_SHAP ~ s(elevation, k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_minT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30')



### meanT
gam_meanT_elev <- gam(Mean_T_SHAP ~ s(elevation, k = 4),
                          na.action = 'na.omit',
                          data = results)

plot.gam(gam_meanT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#80008030')

### maxT
gam_maxT_elev <- gam(Max_T_SHAP ~ s(elevation, k = 4),
                         na.action = 'na.omit',
                         data = results)

plot.gam(gam_maxT_elev, pages = 1, residuals = F, shade = T,
         shade.col = '#FF000030')



#plot minT
plot(results$relPolarwardness, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'relPolewarness',
     ylim = ylim)

plot()

pred_SHAP_minT <- predict(gam_minT_relPolewarness)

points(results$relPolarwardness, pred_minT, cex = 0.1)

plot(results$relPolarwardness, pred_minT, cex = 0.1)

summary(predict(gam_minT_relPolewarness))
predictions <- predict(gam_model, newdata = new_data, type = "response", 
                       se.fit = TRUE)

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(log(distEdge), k = 2),
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


#plot meanT
plot(log(results$distEdge), results$Mean_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_meanT_distEdge))

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(log(distEdge), k = 2),
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


#plot maxT
plot(log(results$distEdge), results$Max_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_maxT_distEdge))

plot(results$relPolarwardness, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Polewardness',
     ylim = ylim)

gam_minT_pole <- gam(Min_T_SHAP ~ relPolarwardness,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

gam_meanT_pole <- gam(Mean_T_SHAP ~ relPolarwardness,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

gam_maxT_pole <- gam(Max_T_SHAP ~ relPolarwardness,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)


summary(gam_minT_pole)
summary(gam_meanT_pole)
summary(gam_maxT_pole)

plot(gam_minT_pole, all.terms = TRUE)
plot(gam_meanT_pole, all.terms = TRUE)
plot(gam_maxT_pole, all.terms = TRUE)


gam.check(gam_minT_pole)

abline(gam_minT_elev, col = '#0000FF', lwd = 2)


# SHAP values X elevation (GAM)

#plot minT
plot(results$elevation, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'Elevation',
     ylim = ylim)

lines(predict(gam))


gam_minT_elev <- gam(Min_T_SHAP ~ elevation,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

gam_meanT_elev <- gam(Mean_T_SHAP ~ elevation,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

gam_maxT_elev <- gam(Max_T_SHAP ~ elevation,
                     family = gaussian(),
                     na.action = 'na.omit',
                     data = results)

summary(gam_minT_elev)
summary(gam_meanT_elev)
summary(gam_maxT_elev)

plot.gam(gam_minT_elev, all.terms = TRUE, add = T)
plot(gam_meanT_elev, all.terms = TRUE)
plot(gam_maxT_elev, all.terms = TRUE)

gam.check(gam_minT_pole)


# SHAP values X centralness (GAM)

gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge),
                         family = gaussian(),
                         k = 3, 
                         na.action = 'na.omit',
                         data = results)


#plot minT
plot(results$distEdge, results$Min_T_SHAP,
     pch = 19, cex = 0.4, col = '#0000FF05',
     ylab = 'SHAP value',
     xlab = 'distEdge',
     ylim = ylim)

lines(predict(gam_minT_distEdge))



gam_meanT_distEdge <- gam(Mean_T_SHAP ~ distEdge,
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)

gam_maxT_distEdge <- gam(Max_T_SHAP ~ distEdge,
                         family = gaussian(),
                         na.action = 'na.omit',
                         data = results)


summary(gam_minT_distEdge)
summary(gam_meanT_distEdge)
summary(gam_maxT_distEdge)


plot(gam_minT_distEdge, all.terms = TRUE)
plot(gam_meanT_distEdge, all.terms = TRUE)
plot(gam_maxT_distEdge, all.terms = TRUE)

#add meanT
points(results$elevation, results$Mean_T_SHAP,
       pch = 19, cex = 0.4, col = '#80008005')

gam_meanT_elev <- gam(Mean_T_SHAP ~ elevation,
                      na.action = 'na.omit',
                      data = results)

abline(gam_meanT_elev, col = '#800080', lwd = 2)


#add maxT
points(results$elevation, results$Max_T_SHAP,
       pch = 19, cex = 0.4, col = '#FF000005')

gam_maxT_elev <- gam(Max_T_SHAP ~ elevation,
                     na.action = 'na.omit',
                     data = results)

abline(gam_maxT_elev, col = '#800080', lwd = 2)


#add stats to the plot
summa_minT_gam <- summary(gam_minT_elev)

r2_minT <- summa_minT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1], 
     expression(bold(paste('R'^2*' = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] - ((ylim[2] - ylim[1]) / 200), 
     round(r2_minT, 3), col = '#0000FF', font = 2, pos = 4)

p_value_minT <- summa_minT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#0000FF', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - 3 * ((ylim[2] - ylim[1]) / 22), 
     round(p_value_minT, 4), col = '#0000FF', font = 2, pos = 4)

summa_meanT <- summary(lm_meanT_range_Size)
r2_meanT <- summa_meanT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_meanT, 3), col = '#800080', font = 2, pos = 4)

p_value_meanT <- summa_meanT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#800080', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - 2 * ((ylim[2] - ylim[1]) / 22), 
     round(p_value_meanT, 4), col = '#800080', font = 2, pos = 4)

summa_maxT <- summary(lm_maxT_range_Size)
r2_maxT <- summa_maxT$r.squared
text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 200), 
     round(r2_maxT, 3), col = '#FF0000', font = 2, pos = 4)


p_value_maxT <- summa_maxT$coefficients[2,4]

text(xlim[2] - ((xlim[2] - xlim[1]) / 5), 
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     expression(bold(paste('p-value = '), sep='')),
     col = '#FF0000', pos = 4)
text(xlim[2] - ((xlim[2] - xlim[1]) / 10),
     ylim[2] - ((ylim[2] - ylim[1]) / 22), 
     round(p_value_maxT, 4), col = '#FF0000', font = 2, pos = 4)