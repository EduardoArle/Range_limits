#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240603'

#read temperature results table

setwd(wd_tables)
res_temp <- read.csv('Temperature_Rel_Polar_all_points.csv')

#calculate difference in slopes

res_temp$diffMeanMinTslope <- res_temp$meanTslope_RP - res_temp$minTslope_RP
res_temp$diffMeanMaxTslope <- res_temp$meanTslope_RP - res_temp$maxTslope_RP

summary(res_temp$diffMeanMinTslope)
hist(res_temp$diffMeanMinTslope)

summary(res_temp$diffMeanMaxTslope)
hist(res_temp$diffMeanMaxTslope)
