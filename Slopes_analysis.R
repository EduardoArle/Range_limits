#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240523'

#read temperature results table
setwd(wd_tables)
res_temp <- read.csv('Temperature_Rel_Polar_all_points.csv')

#discard species with less than 10 records
res_temp <- res_temp[res_temp$n_records >= 10,]

#discard species that cross the equator
res_temp_one_hemis <- res_temp[res_temp$hemisphere != 'Both',]


#read precipitation results table
setwd(wd_tables)
res_prec <- read.csv('Precipitation_Rel_Polar_all_points.csv')

#discard species with less than 10 records
res_prec <- res_prec[res_prec$n_records >= 10,]

#discard species that cross the equator
res_prec_one_hemis <- res_prec[res_prec$hemisphere != 'Both',]


### check slopes temp rel polewardness
sum_minTslopeRP <- sum(res_temp_one_hemis$minTslope_RP)
sum_meanTslopeRP <- sum(res_temp_one_hemis$meanTslope_RP)
sum_maxTslopeRP <- sum(res_temp_one_hemis$maxTslope_RP)

# per p value
sum_minTslopeRP_0.1 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.1], na.rm = T)
sum_meanTslopeRP_0.1 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.1], na.rm = T)
sum_maxTslopeRP_0.1 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.1], na.rm = T)

sum_minTslopeRP_0.05 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.05], na.rm = T)
sum_meanTslopeRP_0.05 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.05], na.rm = T)
sum_maxTslopeRP_0.05 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.05], na.rm = T)

sum_minTslopeRP_0.01 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.01], na.rm = T)
sum_meanTslopeRP_0.01 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.01], na.rm = T)
sum_maxTslopeRP_0.01 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.01], na.rm = T)

sum_minTslopeRP_0.001 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.001], na.rm = T)
sum_meanTslopeRP_0.001 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.001], na.rm = T)
sum_maxTslopeRP_0.001 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.001], na.rm = T)

# per p value
sum_minTslopeRP_0.1 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.1], na.rm = T)
sum_meanTslopeRP_0.1 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.1], na.rm = T)
sum_maxTslopeRP_0.1 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.1], na.rm = T)

sum_minTslopeRP_0.05 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.05], na.rm = T)
sum_meanTslopeRP_0.05 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.05], na.rm = T)
sum_maxTslopeRP_0.05 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.05], na.rm = T)

sum_minTslopeRP_0.01 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$minTp_value_RP < 0.01], na.rm = T)
sum_meanTslopeRP_0.01 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$meanTp_value_RP < 0.01], na.rm = T)
sum_maxTslopeRP_0.01 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$maxTp_value_RP < 0.01], na.rm = T)

# per latitudinal amplitude
#classify species by latitudinal amplitude
res_temp_one_hemis$latAmplitude <-
  res_temp_one_hemis$maxLat - res_temp_one_hemis$minLat

sum_minTslopeRP_1 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$latAmplitude > 1], na.rm = T)
sum_meanTslopeRP_1 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$latAmplitude > 1], na.rm = T)
sum_maxTslopeRP_1 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$latAmplitude > 1], na.rm = T)

sum_minTslopeRP_5 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$latAmplitude > 5], na.rm = T)
sum_meanTslopeRP_5 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$latAmplitude > 5], na.rm = T)
sum_maxTslopeRP_5 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$latAmplitude > 5], na.rm = T)

sum_minTslopeRP_10 <- sum(res_temp_one_hemis$minTslope_RP[
  res_temp_one_hemis$latAmplitude > 10], na.rm = T)
sum_meanTslopeRP_10 <- sum(res_temp_one_hemis$meanTslope_RP[
  res_temp_one_hemis$latAmplitude > 10], na.rm = T)
sum_maxTslopeRP_10 <- sum(res_temp_one_hemis$maxTslope_RP[
  res_temp_one_hemis$latAmplitude > 10], na.rm = T)

# make table p-value
#set margin parametres for the plot
par(mar=c(3,3,3,3))

#make the empty plot
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
     xaxs = "i", yaxs = "i", axes = F, frame.plot = TRUE)

#make lines creating a table (cols)
lines(c(0, 0), c(0, 10), lwd = 3)
lines(c(1, 1), c(0, 10), lwd = 2)
lines(c(4, 4), c(0, 10))
lines(c(7, 7), c(0, 10))
lines(c(10, 10), c(0, 10), lwd = 3)

#make lines creating a table (rows)
lines(c(0, 10), c(0, 0), lwd = 3)
lines(c(0, 10), c(10/6, 10/6))
lines(c(0, 10), c(10/6 *2, 10/6 *2))
lines(c(0, 10), c(10/6 *3, 10/6 *3))
lines(c(0, 10), c(10/6 *4, 10/6 *4))
lines(c(0, 10), c(10/6 *5, 10/6 *5), lwd = 2)
lines(c(0, 10), c(10, 10), lwd = 3)

#add text to the table
text(0.5, 10/6 * 4.5, 'All', srt = 90, cex = 1.5, font = 2)
text(0.5, 5.8, '< 0.1', srt = 90, cex = 1.5, font = 2)
text(0.5, 4.1, '< 0.05', srt = 90, cex = 1.5, font = 2)
text(0.5, 2.5, '< 0.01', srt = 90, cex = 1.5, font = 2)
text(0.5, 0.8, '< 0.001', srt = 90, cex = 1.5, font = 2)

text(2.5, 10/6 * 5.5, 'min T', cex = 1.5, font = 2)
text(5.5, 10/6 * 5.5, 'mean T', cex = 1.5, font = 2)
text(8.5, 10/6 * 5.5, 'max T', cex = 1.5, font = 2)

text(2.5, 10/6 * 4.5, round(sum_minTslopeRP, 3), cex = 1.5)
text(5.5, 10/6 * 4.5, round(sum_meanTslopeRP, 3), cex = 1.5)
text(8.5, 10/6 * 4.5, round(sum_maxTslopeRP, 3), cex = 1.5)

text(2.5, 10/6 * 3.5, round(sum_minTslopeRP_1, 3), cex = 1.5)
text(5.5, 10/6 * 3.5, round(sum_meanTslopeRP_0.1, 3), cex = 1.5)
text(8.5, 10/6 * 3.5, round(sum_maxTslopeRP_0.1, 3), cex = 1.5)

text(2.5, 10/6 * 2.5, round(sum_minTslopeRP_0.05, 3), cex = 1.5)
text(5.5, 10/6 * 2.5, round(sum_meanTslopeRP_0.05, 3), cex = 1.5)
text(8.5, 10/6 * 2.5, round(sum_maxTslopeRP_0.05, 3), cex = 1.5)

text(2.5, 10/6 * 1.5, round(sum_minTslopeRP_0.01, 3), cex = 1.5)
text(5.5, 10/6 * 1.5, round(sum_meanTslopeRP_0.01, 3), cex = 1.5)
text(8.5, 10/6 * 1.5, round(sum_maxTslopeRP_0.01, 3), cex = 1.5)

text(2.5, 10/6 * 0.5, round(sum_minTslopeRP_0.001, 3), cex = 1.5)
text(5.5, 10/6 * 0.5, round(sum_meanTslopeRP_0.001, 3), cex = 1.5)
text(8.5, 10/6 * 0.5, round(sum_maxTslopeRP_0.001, 3), cex = 1.5)

# make table latitudinal amplitude
#set margin parametres for the plot
par(mar=c(3,3,3,3))

#make the empty plot
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
     xaxs = "i", yaxs = "i", axes = F, frame.plot = TRUE)

#make lines creating a table (cols)
lines(c(0, 0), c(0, 10), lwd = 3)
lines(c(1, 1), c(0, 10), lwd = 2)
lines(c(4, 4), c(0, 10))
lines(c(7, 7), c(0, 10))
lines(c(10, 10), c(0, 10), lwd = 3)

#make lines creating a table (rows)
lines(c(0, 10), c(0, 0), lwd = 3)
lines(c(0, 10), c(10/5, 10/5))
lines(c(0, 10), c(10/5 *2, 10/5 *2))
lines(c(0, 10), c(10/5 *3, 10/5 *3))
lines(c(0, 10), c(10/5 *4, 10/5 *4), lwd = 2)
lines(c(0, 10), c(10, 10), lwd = 3)

#add text to the table
text(0.5, 10/5 * 3.5, 'All', srt = 90, cex = 1.5, font = 2)
text(0.5, 10/5 * 2.5, '> 1', srt = 90, cex = 1.5, font = 2)
text(0.5, 10/5 * 1.5, '> 5', srt = 90, cex = 1.5, font = 2)
text(0.5, 10/5 * 0.5, '> 10', srt = 90, cex = 1.5, font = 2)

text(2.5, 10/6 * 5.5, 'min T', cex = 1.5, font = 2)
text(5.5, 10/6 * 5.5, 'mean T', cex = 1.5, font = 2)
text(8.5, 10/6 * 5.5, 'max T', cex = 1.5, font = 2)

text(2.5, 10/5 * 3.5, round(sum_minTslopeRP, 3), cex = 1.5)
text(5.5, 10/5 * 3.5, round(sum_meanTslopeRP, 3), cex = 1.5)
text(8.5, 10/5 * 3.5, round(sum_maxTslopeRP, 3), cex = 1.5)

text(2.5, 10/5 * 2.5, round(sum_minTslopeRP_1, 3), cex = 1.5)
text(5.5, 10/5 * 2.5, round(sum_meanTslopeRP_1, 3), cex = 1.5)
text(8.5, 10/5 * 2.5, round(sum_maxTslopeRP_1, 3), cex = 1.5)

text(2.5, 10/5 * 1.5, round(sum_minTslopeRP_5, 3), cex = 1.5)
text(5.5, 10/5 * 1.5, round(sum_meanTslopeRP_5, 3), cex = 1.5)
text(8.5, 10/5 * 1.5, round(sum_maxTslopeRP_5, 3), cex = 1.5)

text(2.5, 10/5 * 0.5, round(sum_minTslopeRP_10, 3), cex = 1.5)
text(5.5, 10/5 * 0.5, round(sum_meanTslopeRP_10, 3), cex = 1.5)
text(8.5, 10/5 * 0.5, round(sum_maxTslopeRP_10, 3), cex = 1.5)



### check slopes temp rel centralness
sum_minTslopeC <- sum(res_temp_one_hemis$minTslope_C)
sum_meanTslopeC <- sum(res_temp_one_hemis$meanTslope_C)
sum_maxTslopeC <- sum(res_temp_one_hemis$maxTslope_C)


### check slopes prec rel polewardness
sum_minPPTslopeRP <- sum(res_prec_one_hemis$minPPTslope_RP)
sum_meanPPTslopeRP <- sum(res_prec_one_hemis$meanPPTslope_RP)
sum_maxPPTslopeRP <- sum(res_prec_one_hemis$maxPPTslope_RP)

### check slopes prec rel centralness
sum_minPPTslopeC <- sum(res_prec_one_hemis$minPPTslope_C)
sum_meanPPTslopeC <- sum(res_prec_one_hemis$meanPPTslope_C)
sum_maxPPTslopeC <- sum(res_prec_one_hemis$maxPPTslope_C)
