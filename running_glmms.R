#We will use this code to run GLMMs and test hypotheses regarding temperature and precipitations


#Packs
library(tidyverse)
library(glmmTMB)

#Data
shap_data=read.csv("Results_all_sps(in).csv") %>%
  mutate(abs_lat=abs(decimalLatitude),
         sqrt_n_occ=sqrt(n_occ)) %>% #add abs lat
  select(species,
         key,
         abs_lat,
         rangeSize,
         distEdge,
         relPolarwardness,
         elevation,
         biome,
         bodyMass,
         sqrt_n_occ,
         order,
         Min_T_SHAP,
         Mean_T_SHAP,
         Max_T_SHAP,
         Min_PPT_SHAP,
         Mean_PPT_SHAP,
         Max_PPT_SHAP,
         NOTE) %>% 
  filter(!NOTE%in%"Crossed by Equator")
shap_data$biome=as.factor(shap_data$biome)
#Running models

##relative polewardness X Range size
###Min_T_SHAP
Min_T_SHAP_model=glmmTMB(
                         data=shap_data,
                         formula = Min_T_SHAP~relPolarwardness*log(rangeSize)*elevation+(1|order)+(1|biome),
                         weights = sqrt_n_occ)
####model diag
sjPlot::plot_model(Min_T_SHAP_model,type="diag")
glmmTMB::diagnose(Min_T_SHAP_model)
####model summary
summary(Min_T_SHAP_model)

sjPlot::plot_model(Min_T_SHAP_model,type="pred",
                   terms=c('relPolarwardness','rangeSize[20000,1000000]','elevation[0,1000,6000]'))
##############
###Max_T_SHAP#
##############
###Min_T_SHAP
Max_T_SHAP_model=glmmTMB(
  data=shap_data,
  formula = Max_T_SHAP~relPolarwardness*log(rangeSize)*elevation+(1|order)+(1|biome),
  weights = sqrt_n_occ)
####model diag
sjPlot::plot_model(Max_T_SHAP_model,type="diag")
glmmTMB::diagnose(Max_T_SHAP_model)
####model summary
summary(Max_T_SHAP_model)

sjPlot::plot_model(Max_T_SHAP_model,type="pred",
                   terms=c('relPolarwardness','rangeSize[20000,1000000]','elevation[0,1000,6000]'))

ggplot(shap_data)+
  geom_point(aes(x=Min_T_SHAP,y=Max_T_SHAP))
