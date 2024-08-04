#load libraries
library(mgcv); library(itsadug)

data(simdat)

## Not run: 
# Model with random effect and interactions:
m1 <- gam(Y ~ s(Time, by=Group, k = 4)
          + s(Time, Subject, bs='fs', k = 4),
          data=simdat)

gam_minT_relPolewarness_RE <- gam(Min_T_SHAP ~ s(relPolewardness, k = 4) 
                                  + s(relPolewardness, species, bs = "re"),
                                  data = results)

summary(gam_minT_relPolewarness_RE)

plot.gam(m1, pages = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', all.terms = T,
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5)

inspect_random(m1, select=3)

inspect_random(gam_minT_relPolewarness_RE, select = 3)

# Simple model with smooth:
m2 <-  gam(Y ~ Group + s(Time, by=Group)
           + s(Subject, bs='re')
           + s(Subject, Time, bs='re'),
           data=simdat)

par(mfrow=c(1,2), cex=1.1)

plot(m2, select=1, shade=TRUE, rug=FALSE, ylim=c(-15,10))
abline(h=0)
plot(m2, select=2, shade=TRUE, rug=FALSE, ylim=c(-15,10))
abline(h=0)


par(mfrow=c(1,1), cex=1.1)



inspect_random(gam_minT_relPolewarness, select = 2)

children <- unique(simdat[simdat$Group=="Children", "Subject"])
adults   <- unique(simdat[simdat$Group=="Adults", "Subject"])

inspect_random(m1, select=3, main='Averages', 
               fun=mean, 
               cond=list(Subject=children))
inspect_random(m1, select=3, 
               fun=mean, cond=list(Subject=adults),
               add=TRUE, col='red', lty=5)

# add legend:
legend('bottomleft',
       legend=c('Children', 'Adults'),
       col=c('black', 'red'), lty=c(1,5),
       bty='n')

plot.gam(gam_minT_relPolewarness_RE, select=1, shade=TRUE, rug=FALSE)


class(gam_minT_relPolewarness_RE)

#### model G


CO2_modG <- gam(log(uptake) ~
                  s(log(conc), k=5, bs="tp") 
                + s(Plant_uo, k=12, bs="re"),
                data=CO2,
                method="REML",
                family="gaussian")

summary(CO2_modG)

draw(CO2_modG, page = 1)


minT_relPolewarness_G <- gam(Min_T_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minT_relPolewarness_G)

draw(minT_relPolewarness_G, page = 1)


#### model G

#test
CO2_modG <- gam(log(uptake) ~
                  s(log(conc), k=5, bs="tp") 
                + s(Plant_uo, k=12, bs="re"),
                data=CO2,
                method="REML",
                family="gaussian")

summary(CO2_modG)


draw(CO2_modG, page = 1)

#minT relPol
minT_relPolewarness_G <- gam(Min_T_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minT_relPolewarness_G)

draw(minT_relPolewarness_G, page = 1)


#meanT relPol
meanT_relPolewarness_G <- gam(Mean_T_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(meanT_relPolewarness_G)

draw(meanT_relPolewarness_G, page = 1)

#maxT relPol
maxT_relPolewarness_G <- gam(Max_T_SHAP ~
                                s(relPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(maxT_relPolewarness_G)

draw(maxT_relPolewarness_G, page = 1)


#minT absPol
minT_absPolewarness_G <- gam(Min_T_SHAP ~
                               s(absPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minT_absPolewarness_G)

draw(minT_absPolewarness_G, page = 1)


#meanT relPol
meanT_absPolewarness_G <- gam(Mean_T_SHAP ~
                                s(absPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(meanT_absPolewarness_G)

draw(meanT_absPolewarness_G, page = 1)

#maxT relPol
maxT_relPolewarness_G <- gam(Max_T_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(maxT_relPolewarness_G)

draw(maxT_relPolewarness_G, page = 1)


#### model GS

#test
CO2_modGS <- gam(log(uptake) ~
                   s(log(conc), k=5, m=2)
                 + s(log(conc), Plant_uo, k=5, bs="fs", m=2),
                 data=CO2,
                 method="REML")

summary(CO2_modGS)

draw(CO2_modGS, page = 1)

#minT relPol
minT_relPolewarness_GS <- gam(Min_T_SHAP ~
                               s(relPolewardness, k=4, m=2) 
                             + s(relPolewardness, species, k=4, bs="fs", m=2),
                             data = results,
                             method="REML")


summary(minT_relPolewarness_GS)

draw(minT_relPolewarness_GS, page = 1)

#meanT relPol
meanT_relPolewarness_GS <- gam(Mean_T_SHAP ~
                                s(relPolewardness, k=4, m=2) 
                              + s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(meanT_relPolewarness_GS)

draw(meanT_relPolewarness_GS, page = 1)

#maxT relPol
maxT_relPolewarness_GS <- gam(Max_T_SHAP ~
                                 s(relPolewardness, k=4, m=2) 
                               + s(relPolewardness, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")


summary(maxT_relPolewarness_GS)

draw(maxT_relPolewarness_GS, page = 1)


#model S

#test 
CO2_modS <- gam(log(uptake) ~
                  s(log(conc), Plant_uo, k=5, bs="fs", m=2),
                  data=CO2,
                  method="REML")

summary(CO2_modS)

draw(CO2_modS, page = 1)

#minT relPol
minT_relPolewarness_S <- gam(Min_T_SHAP ~
                              s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(minT_relPolewarness_S)

draw(minT_relPolewarness_S, page = 1)
