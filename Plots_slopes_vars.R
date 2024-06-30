#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240603'

#read temperature results table

setwd(wd_tables)
res_temp <- read.csv('Temperature_Rel_Polar_all_points.csv')

#fix 0s where there should be NAs in the table with results
res_temp$minTslope_RP[is.na(res_temp$minTR2_RP)] <- NA
res_temp$minTSE_RP[is.na(res_temp$minTR2_RP)] <- NA

res_temp$meanTslope_RP[is.na(res_temp$meanTR2_RP)] <- NA
res_temp$meanTSE_RP[is.na(res_temp$meanTR2_RP)] <- NA

res_temp$maxTslope_RP[is.na(res_temp$maxTR2_RP)] <- NA
res_temp$maxTSE_RP[is.na(res_temp$maxTR2_RP)] <- NA

#calculate species' latitudinal amplitude
res_temp$latAmplitude <- res_temp$maxLat - res_temp$minLat

#calculate species' elevation amplitude
res_temp$elevAmplitude <- res_temp$max_elev - res_temp$min_elev

#discard species with less than 10 points
res_temp <- res_temp[res_temp$n_records > 9,]

#discard species that cross the equator
res_temp_one_hemis <- res_temp[res_temp$hemisphere != 'Both',]


###### Plot the slopes per species against each variables in the RF

## Set graphical parametres

par(mar = c(5,5,4,4))

## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$rangeSize) #no matter which

##### RANGE SIZE ####

# SHAP slopeRP X rangeSize

#plot minT
plot(res_temp_one_hemis$rangeSize, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Range Size',
     ylim = ylim)

lm_minT_range_Size <- lm(minTslope_RP ~ rangeSize,
                         data = res_temp_one_hemis)

abline(lm_minT_range_Size, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$rangeSize, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_range_Size <- lm(meanTslope_RP ~ rangeSize,
                          data = res_temp_one_hemis)

abline(lm_meanT_range_Size, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$rangeSize, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_range_Size <- lm(maxTslope_RP ~ rangeSize,
                         data = res_temp_one_hemis)

abline(lm_maxT_range_Size, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_range_Size)

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


# SHAP slopeRP X rangeSize (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$rangeSize, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Range Size',
     ylim = ylim)

lm_minT_range_Size <- lm(minTslope_RP ~ rangeSize,
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_range_Size, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$rangeSize, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_range_Size <- lm(meanTslope_RP ~ rangeSize,
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_range_Size, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$rangeSize, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_range_Size <- lm(maxTslope_RP ~ rangeSize,
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_range_Size, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_range_Size)

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


# SHAP slopeRP X log rangeSize

#plot minT
plot(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Range Size (log)',
     ylim = ylim)

lm_minT_range_Size <- lm(minTslope_RP ~ log(rangeSize),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_range_Size, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_range_Size <- lm(meanTslope_RP ~ log(rangeSize),
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_range_Size, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_range_Size <- lm(maxTslope_RP ~ log(rangeSize),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_range_Size, col = '#FF0000', lwd = 2)


#change xlim to log(xlim)
xlim <- log(xlim)

#add stats to the plot
summa_minT <- summary(lm_minT_range_Size)

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


# SHAP slopeRP X log rangeSize (weighted by standard error)

#plot minT
plot(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Range Size (log)',
     ylim = ylim)

lm_minT_range_Size <- lm(minTslope_RP ~ log(rangeSize),
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_range_Size, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_range_Size <- lm(meanTslope_RP ~ log(rangeSize),
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_range_Size, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$rangeSize), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_range_Size <- lm(maxTslope_RP ~ log(rangeSize),
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_range_Size, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_range_Size)

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


##### LATITUDINAL AMPLITUDE ####


## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$latAmplitude) #no matter which



# SHAP slopeRP X Latitudinal Amplitude

#plot minT
plot(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Latitudinal Amplitude',
     ylim = ylim)

lm_minT_lat_ampli <- lm(minTslope_RP ~ latAmplitude,
                         data = res_temp_one_hemis)

abline(lm_minT_lat_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_lat_ampli <- lm(meanTslope_RP ~ latAmplitude,
                          data = res_temp_one_hemis)

abline(lm_meanT_lat_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_lat_ampli <- lm(maxTslope_RP ~ latAmplitude,
                         data = res_temp_one_hemis)

abline(lm_maxT_lat_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_lat_ampli)

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

summa_meanT <- summary(lm_meanT_lat_ampli)
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

summa_maxT <- summary(lm_maxT_lat_ampli)
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


# SHAP slopeRP X latitudinal amplitude (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Latitudinal Amplitude',
     ylim = ylim)

lm_minT_lat_ampli <- lm(minTslope_RP ~ latAmplitude,
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_lat_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_lat_ampli <- lm(meanTslope_RP ~ latAmplitude,
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_lat_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$latAmplitude, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_lat_ampli <- lm(maxTslope_RP ~ latAmplitude,
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_lat_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_lat_ampli)

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

summa_meanT <- summary(lm_meanT_lat_ampli)
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

summa_maxT <- summary(lm_maxT_lat_ampli)
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


# SHAP slopeRP X log latitudinal amplitude

#plot minT
plot(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Latitudinal Amplitude (log)',
     ylim = ylim)

lm_minT_lat_ampli <- lm(minTslope_RP ~ log(latAmplitude),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_lat_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_lat_ampli <- lm(meanTslope_RP ~ log(latAmplitude),
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_lat_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_lat_ampli <- lm(maxTslope_RP ~ log(latAmplitude),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_lat_ampli, col = '#FF0000', lwd = 2)


#change xlim to log(xlim)
xlim <- log(xlim)

#add stats to the plot
summa_minT <- summary(lm_minT_lat_ampli)

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

summa_meanT <- summary(lm_meanT_lat_ampli)
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

summa_maxT <- summary(lm_maxT_lat_ampli)
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


# SHAP slopeRP X latitudinal amplitude (weighted by standard error)

#plot minT
plot(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Latitudinal Amplitude (log)',
     ylim = ylim)

lm_minT_lat_ampli <- lm(minTslope_RP ~ log(latAmplitude),
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_lat_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_lat_ampli <- lm(meanTslope_RP ~ log(latAmplitude),
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_lat_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$latAmplitude), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_lat_ampli <- lm(maxTslope_RP ~ log(latAmplitude),
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_lat_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_lat_ampli)

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

summa_meanT <- summary(lm_meanT_lat_ampli)
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

summa_maxT <- summary(lm_maxT_lat_ampli)
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





##### ELEVATION AMPLITUDE ####


## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$elevAmplitude) #no matter which


# SHAP slopeRP X Elevation Amplitude

#plot minT
plot(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Elevation Amplitude',
     ylim = ylim)

lm_minT_elev_ampli <- lm(minTslope_RP ~ elevAmplitude,
                        data = res_temp_one_hemis)

abline(lm_minT_elev_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_elev_ampli <- lm(meanTslope_RP ~ elevAmplitude,
                         data = res_temp_one_hemis)

abline(lm_meanT_elev_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_elev_ampli <- lm(maxTslope_RP ~ elevAmplitude,
                        data = res_temp_one_hemis)

abline(lm_maxT_elev_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_elev_ampli)

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

summa_meanT <- summary(lm_meanT_elev_ampli)
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

summa_maxT <- summary(lm_maxT_elev_ampli)
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


# SHAP slopeRP X elevation amplitude (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Elevation Amplitude',
     ylim = ylim)

lm_minT_elev_ampli <- lm(minTslope_RP ~ elevAmplitude,
                        weights = 1 / minTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_minT_elev_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_elev_ampli <- lm(meanTslope_RP ~ elevAmplitude,
                         weights = 1 / meanTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_meanT_elev_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$elevAmplitude, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_elev_ampli <- lm(maxTslope_RP ~ elevAmplitude,
                        weights = 1 / maxTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_maxT_elev_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_elev_ampli)

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

summa_meanT <- summary(lm_meanT_elev_ampli)
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

summa_maxT <- summary(lm_maxT_elev_ampli)
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


# SHAP slopeRP X log latitudinal amplitude

#plot minT
plot(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Elevation Amplitude (log)',
     ylim = ylim)

lm_minT_elev_ampli <- lm(minTslope_RP ~ log(elevAmplitude),
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_minT_elev_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_elev_ampli <- lm(meanTslope_RP ~ log(elevAmplitude),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_meanT_elev_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_elev_ampli <- lm(maxTslope_RP ~ log(elevAmplitude),
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_maxT_elev_ampli, col = '#FF0000', lwd = 2)


#change xlim to log(xlim)
xlim <- log(xlim)

#add stats to the plot
summa_minT <- summary(lm_minT_elev_ampli)

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

summa_meanT <- summary(lm_meanT_elev_ampli)
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

summa_maxT <- summary(lm_maxT_elev_ampli)
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


# SHAP slopeRP X elevation amplitude (weighted by standard error)

#plot minT
plot(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Elevation Amplitude (log)',
     ylim = ylim)

lm_minT_elev_ampli <- lm(minTslope_RP ~ log(elevAmplitude),
                        weights = 1 / minTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_minT_elev_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_elev_ampli <- lm(meanTslope_RP ~ log(elevAmplitude),
                         weights = 1 / meanTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_meanT_elev_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$elevAmplitude), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_elev_ampli <- lm(maxTslope_RP ~ log(elevAmplitude),
                        weights = 1 / maxTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_maxT_elev_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_elev_ampli)

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

summa_meanT <- summary(lm_meanT_elev_ampli)
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

summa_maxT <- summary(lm_maxT_elev_ampli)
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


##### NUMBER OF RECORDS ####


## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$n_records) #no matter which


# SHAP slopeRP X Number of records

#plot minT
plot(res_temp_one_hemis$n_records, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Number of records',
     ylim = ylim)

lm_minT_nRec_ampli <- lm(minTslope_RP ~ n_records,
                         data = res_temp_one_hemis)

abline(lm_minT_nRec_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$n_records, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_nRec_ampli <- lm(meanTslope_RP ~ n_records,
                          data = res_temp_one_hemis)

abline(lm_meanT_nRec_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$n_records, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_nRec_ampli <- lm(maxTslope_RP ~ n_records,
                         data = res_temp_one_hemis)

abline(lm_maxT_nRec_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_nRec_ampli)

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

summa_meanT <- summary(lm_meanT_nRec_ampli)
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

summa_maxT <- summary(lm_maxT_nRec_ampli)
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


# SHAP slopeRP X Number of records (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$n_records, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Number of records',
     ylim = ylim)

lm_minT_nRec_ampli <- lm(minTslope_RP ~ n_records,
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_nRec_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$n_records, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_nRec_ampli <- lm(meanTslope_RP ~ n_records,
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_nRec_ampli, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$n_records, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_nRec_ampli <- lm(maxTslope_RP ~ n_records,
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_nRec_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_nRec_ampli)

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

summa_meanT <- summary(lm_meanT_nRec_ampli)
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

summa_maxT <- summary(lm_maxT_nRec_ampli)
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


# SHAP slopeRP X log Number of records

#plot minT
plot(log(res_temp_one_hemis$n_records), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Number of records (log)',
     ylim = ylim)

lm_minT_nRec_ampli <- lm(minTslope_RP ~ log(n_records),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_nRec_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$n_records), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_nRec_ampli <- lm(meanTslope_RP ~ log(n_records),
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_nRec_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$n_records), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_nRec_ampli <- lm(maxTslope_RP ~ log(n_records),
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_nRec_ampli, col = '#FF0000', lwd = 2)


#change xlim to log(xlim)
xlim <- log(xlim)

#add stats to the plot
summa_minT <- summary(lm_minT_nRec_ampli)

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

summa_meanT <- summary(lm_meanT_nRec_ampli)
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

summa_maxT <- summary(lm_maxT_nRec_ampli)
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


# SHAP slopeRP X Number of records (weighted by standard error)

#plot minT
plot(log(res_temp_one_hemis$n_records), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Number of records (log)',
     ylim = ylim)

lm_minT_nRec_ampli <- lm(minTslope_RP ~ log(n_records),
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_nRec_ampli, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$n_records), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_nRec_ampli <- lm(meanTslope_RP ~ log(n_records),
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_nRec_ampli, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$n_records), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_nRec_ampli <- lm(maxTslope_RP ~ log(n_records),
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_nRec_ampli, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_nRec_ampli)

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

summa_meanT <- summary(lm_meanT_nRec_ampli)
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

summa_maxT <- summary(lm_maxT_nRec_ampli)
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





##### MEAN LATITUDE ####

#make column with mean latitude
res_temp_one_hemis$meanLat <-
  (res_temp_one_hemis$maxLat + res_temp_one_hemis$minLat) / 2

## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$meanLat) #no matter which


# SHAP slopeRP X mean latitude

#plot minT
plot(res_temp_one_hemis$meanLat, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Mean Latitude',
     ylim = ylim)

lm_minT_meanLat <- lm(minTslope_RP ~ meanLat,
                         data = res_temp_one_hemis)

abline(lm_minT_meanLat, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$meanLat, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_meanLat <- lm(meanTslope_RP ~ meanLat,
                          data = res_temp_one_hemis)

abline(lm_meanT_meanLat, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$meanLat, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_meanLat <- lm(maxTslope_RP ~ meanLat,
                         data = res_temp_one_hemis)

abline(lm_maxT_meanLat, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_meanLat)

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

summa_meanT <- summary(lm_meanT_meanLat)
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

summa_maxT <- summary(lm_maxT_meanLat)
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


# SHAP slopeRP X mean Latitude (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$meanLat, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Mean Latitude',
     ylim = ylim)

lm_minT_meanLat <- lm(minTslope_RP ~ meanLat,
                         weights = 1 / minTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_minT_meanLat, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$meanLat, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_meanLat <- lm(meanTslope_RP ~ meanLat,
                          weights = 1 / meanTSE_RP,
                          na.action = 'na.omit',
                          data = res_temp_one_hemis)

abline(lm_meanT_meanLat, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$meanLat, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_meanLat <- lm(maxTslope_RP ~ meanLat,
                         weights = 1 / maxTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_maxT_meanLat, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_meanLat)

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

summa_meanT <- summary(lm_meanT_meanLat)
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

summa_maxT <- summary(lm_maxT_meanLat)
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





##### MEAN ELEVATION ####

## Get slope values for all species

#minT
cont_minT <- res_temp_one_hemis$minTslope_RP

#meanT
cont_meanT <- res_temp_one_hemis$meanTslope_RP

#maxT
cont_maxT <- res_temp_one_hemis$maxTslope_RP

## Define y and x lims

ylim <- range(c(cont_minT, cont_meanT, cont_maxT), na.rm = T)
xlim <- range(res_temp_one_hemis$mean_elev) #no matter which


# SHAP slopeRP X mean elevation

#plot minT
plot(res_temp_one_hemis$mean_elev, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Mean Elevation',
     ylim = ylim)

lm_minT_mean_elev <- lm(minTslope_RP ~ mean_elev,
                      data = res_temp_one_hemis)

abline(lm_minT_mean_elev, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$mean_elev, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_mean_elev <- lm(meanTslope_RP ~ mean_elev,
                       data = res_temp_one_hemis)

abline(lm_meanT_mean_elev, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$mean_elev, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_mean_elev <- lm(maxTslope_RP ~ mean_elev,
                      data = res_temp_one_hemis)

abline(lm_maxT_mean_elev, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_mean_elev)

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

summa_meanT <- summary(lm_meanT_mean_elev)
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

summa_maxT <- summary(lm_maxT_mean_elev)
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


# SHAP slopeRP X mean Elevation' (weighted by standard error)

#plot minT
plot(res_temp_one_hemis$mean_elev, res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Mean Elevation',
     ylim = ylim)

lm_minT_mean_elev <- lm(minTslope_RP ~ mean_elev,
                      weights = 1 / minTSE_RP,
                      na.action = 'na.omit',
                      data = res_temp_one_hemis)

abline(lm_minT_mean_elev, col = '#0000FF', lwd = 2)

#add meanT
points(res_temp_one_hemis$mean_elev, res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_mean_elev <- lm(meanTslope_RP ~ mean_elev,
                       weights = 1 / meanTSE_RP,
                       na.action = 'na.omit',
                       data = res_temp_one_hemis)

abline(lm_meanT_mean_elev, col = '#800080', lwd = 2)

#add maxT
points(res_temp_one_hemis$mean_elev, res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_mean_elev <- lm(maxTslope_RP ~ mean_elev,
                      weights = 1 / maxTSE_RP,
                      na.action = 'na.omit',
                      data = res_temp_one_hemis)

abline(lm_maxT_mean_elev, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_mean_elev)

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

summa_meanT <- summary(lm_meanT_mean_elev)
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

summa_maxT <- summary(lm_maxT_mean_elev)
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


# SHAP slopeRP X log mean elevation

#plot minT
plot(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP',
     xlab = 'Mean Elevation (log)',
     ylim = ylim)

lm_minT_mean_elev <- lm(minTslope_RP ~ log(mean_elev),
                        data = res_temp_one_hemis)

abline(lm_minT_mean_elev, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_mean_elev <- lm(meanTslope_RP ~ log(mean_elev),
                         data = res_temp_one_hemis)

abline(lm_meanT_mean_elev, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_mean_elev <- lm(maxTslope_RP ~ log(mean_elev),
                        data = res_temp_one_hemis)

abline(lm_maxT_mean_elev, col = '#FF0000', lwd = 2)

#change xlim to log(xlim)
xlim <- log(xlim)

#add stats to the plot
summa_minT <- summary(lm_minT_mean_elev)

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

summa_meanT <- summary(lm_meanT_mean_elev)
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

summa_maxT <- summary(lm_maxT_mean_elev)
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


# SHAP slopeRP X log mean Elevation' (weighted by standard error)

#plot minT
plot(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$minTslope_RP,
     pch = 19, cex = 0.4, col = '#0000FF20',
     ylab = 'slope RP (weighted)',
     xlab = 'Mean Elevation (log)',
     ylim = ylim)

lm_minT_mean_elev <- lm(minTslope_RP ~ log(mean_elev),
                        weights = 1 / minTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_minT_mean_elev, col = '#0000FF', lwd = 2)

#add meanT
points(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$meanTslope_RP,
       pch = 19, cex = 0.4, col = '#80008020')

lm_meanT_mean_elev <- lm(meanTslope_RP ~ log(mean_elev),
                         weights = 1 / meanTSE_RP,
                         na.action = 'na.omit',
                         data = res_temp_one_hemis)

abline(lm_meanT_mean_elev, col = '#800080', lwd = 2)

#add maxT
points(log(res_temp_one_hemis$mean_elev), res_temp_one_hemis$maxTslope_RP,
       pch = 19, cex = 0.4, col = '#FF000020')

lm_maxT_mean_elev <- lm(maxTslope_RP ~ log(mean_elev),
                        weights = 1 / maxTSE_RP,
                        na.action = 'na.omit',
                        data = res_temp_one_hemis)

abline(lm_maxT_mean_elev, col = '#FF0000', lwd = 2)


#add stats to the plot
summa_minT <- summary(lm_minT_mean_elev)

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

summa_meanT <- summary(lm_meanT_mean_elev)
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

summa_maxT <- summary(lm_maxT_mean_elev)
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

