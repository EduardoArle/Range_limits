#load packagesx
library(Boruta)

#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Slopes'

#read slopes table
setwd(wd_slopes)
slopes <- read.csv('Slopes.csv')



###### MinT vs RELATIVE POLEWARDNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minT_relPol <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$Cor_vars_minT),]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$slope_minT_relPol),]


#transform categorical variable in factor
s_minT_relPol$order <- as.factor(s_minT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_minT_relPol$rangeSize_log10 <- log(s_minT_relPol$rangeSize, 10)
s_minT_relPol$nOcc_log <- log(s_minT_relPol$nOcc)

#run RF through boruta
boruta.minT_relPol <- Boruta(slope_minT_relPol ~
                        rangeSize_log10 +
                        rangeLoc +
                        latAmplitude +
                        roundness +
                        nOcc_log +
                        order +
                        elevMedian +
                        bodyMass +
                        elevAmplitude,
                      data = s_minT_relPol, doTrace = 2, ntree = 500,
                      maxRuns = 500)


plot(boruta.minT_relPol)
print(boruta.minT_relPol)
attStats(boruta.minT_relPol)

# save 2000




###  MEAN  ###

#select only species that had lower correl between vars
s_meanT_relPol <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$Cor_vars_meanT),]
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$slope_meanT_relPol),]

#transform categorical variable in factor
s_meanT_relPol$order <- as.factor(s_meanT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_meanT_relPol$rangeSize_log10 <- log(s_meanT_relPol$rangeSize, 10)
s_meanT_relPol$nOcc_log <- log(s_meanT_relPol$nOcc)

#run RF through boruta
boruta.meanT_relPol <- Boruta(slope_meanT_relPol ~
                               rangeSize_log10 +
                               rangeLoc +
                               latAmplitude +
                               roundness +
                               nOcc_log +
                               order +
                               elevMedian +
                               bodyMass +
                               elevAmplitude,
                             data = s_meanT_relPol, doTrace = 2, ntree = 500,
                             maxRuns = 500)


plot(boruta.meanT_relPol)
print(boruta.meanT_relPol)
attStats(boruta.meanT_relPol)


###  MAX  ###

#select only species that had lower correl between vars
s_maxT_relPol <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$Cor_vars_maxT),]
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$slope_maxT_relPol),]

#transform categorical variable in factor
s_maxT_relPol$order <- as.factor(s_maxT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_maxT_relPol$rangeSize_log10 <- log(s_maxT_relPol$rangeSize, 10)
s_maxT_relPol$nOcc_log <- log(s_maxT_relPol$nOcc)

#run RF through boruta
boruta.maxT_relPol <- Boruta(slope_maxT_relPol ~
                                rangeSize_log10 +
                                rangeLoc +
                                latAmplitude +
                                roundness +
                                nOcc_log +
                                order +
                                elevMedian +
                                bodyMass +
                                elevAmplitude,
                              data = s_maxT_relPol, doTrace = 2, ntree = 500,
                              maxRuns = 500)


plot(boruta.maxT_relPol)
print(boruta.maxT_relPol)
attStats(boruta.maxT_relPol)








###### TEMPERATURE vs CENTRALNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minT_distEdge <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_distEdge <- s_minT_distEdge[complete.cases(s_minT_distEdge$Cor_vars_minT),]
s_minT_distEdge <- s_minT_distEdge[
  complete.cases(s_minT_distEdge$slope_minT_distEdge),]


#transform categorical variable in factor
s_minT_distEdge$order <- as.factor(s_minT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_minT_distEdge$rangeSize_log10 <- log(s_minT_distEdge$rangeSize, 10)
s_minT_distEdge$nOcc_log <- log(s_minT_distEdge$nOcc)

#run RF through boruta
boruta.minT_distEdge <- Boruta(slope_minT_distEdge ~
                               rangeSize_log10 +
                               rangeLoc +
                               latAmplitude +
                               roundness +
                               nOcc_log +
                               order +
                               elevMedian +
                               bodyMass +
                               elevAmplitude,
                             data = s_minT_distEdge, doTrace = 2, ntree = 500,
                             maxRuns = 500)


plot(boruta.minT_distEdge)
print(boruta.minT_distEdge)
attStats(boruta.minT_relPol)


###  MEAN  ###

#select only species that had lower correl between vars
s_meanT_distEdge <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$Cor_vars_meanT),]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$slope_meanT_distEdge),]


#transform categorical variable in factor
s_meanT_distEdge$order <- as.factor(s_meanT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_meanT_distEdge$rangeSize_log10 <- log(s_meanT_distEdge$rangeSize, 10)
s_meanT_distEdge$nOcc_log <- log(s_meanT_distEdge$nOcc)

#run RF through boruta
boruta.meanT_distEdge <- Boruta(slope_meanT_distEdge ~
                                rangeSize_log10 +
                                rangeLoc +
                                latAmplitude +
                                roundness +
                                nOcc_log +
                                order +
                                elevMedian +
                                bodyMass +
                                elevAmplitude,
                              data = s_meanT_distEdge, doTrace = 2, ntree = 500,
                              maxRuns = 500)


plot(boruta.meanT_distEdge)
print(boruta.meanT_distEdge)
attStats(boruta.meanT_distEdge)


###  MAX  ###

#select only species that had lower correl between vars
s_maxT_distEdge <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_distEdge <- s_maxT_distEdge[complete.cases(s_maxT_distEdge$Cor_vars_maxT),]
s_maxT_distEdge <- s_maxT_distEdge[
  complete.cases(s_maxT_distEdge$slope_maxT_distEdge),]

#transform categorical variable in factor
s_maxT_distEdge$order <- as.factor(s_maxT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_maxT_distEdge$rangeSize_log10 <- log(s_maxT_distEdge$rangeSize, 10)
s_maxT_distEdge$nOcc_log <- log(s_maxT_distEdge$nOcc)

#run RF through boruta
boruta.maxT_distEdge <- Boruta(slope_maxT_distEdge ~
                               rangeSize_log10 +
                               rangeLoc +
                               latAmplitude +
                               roundness +
                               nOcc_log +
                               order +
                               elevMedian +
                               bodyMass +
                               elevAmplitude,
                             data = s_maxT_distEdge, doTrace = 2, ntree = 500,
                             maxRuns = 500)


plot(boruta.maxT_distEdge)
print(boruta.maxT_distEdge)
attStats(boruta.maxT_distEdge)




###### PRECIPITATION vs RELATIVE POLEWARDNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minPPT_relPol <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_relPol <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$Cor_vars_minPPT),]
s_minPPT_relPol <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_minPPT_relPol),]


#transform categorical variable in factor
s_minPPT_relPol$order <- as.factor(s_minPPT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_minPPT_relPol$rangeSize_log10 <- log(s_minPPT_relPol$rangeSize, 10)
s_minPPT_relPol$nOcc_log <- log(s_minPPT_relPol$nOcc)

#run RF through boruta
boruta.minPPT_relPol <- Boruta(slope_minPPT_relPol ~
                               rangeSize_log10 +
                               rangeLoc +
                               latAmplitude +
                               roundness +
                               nOcc_log +
                               order +
                               elevMedian +
                               bodyMass +
                               elevAmplitude,
                             data = s_minPPT_relPol, doTrace = 2, ntree = 500,
                             maxRuns = 500)


plot(boruta.minPPT_relPol)
print(boruta.minPPT_relPol)
attStats(boruta.minPPT_relPol)


###  MEAN  ###

#select only species that had lower correl between vars
s_meanPPT_relPol <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_relPol <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$Cor_vars_meanPPT),]
s_meanPPT_relPol <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_meanPPT_relPol),]


#transform categorical variable in factor
s_meanPPT_relPol$order <- as.factor(s_meanPPT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_meanPPT_relPol$rangeSize_log10 <- log(s_meanPPT_relPol$rangeSize, 10)
s_meanPPT_relPol$nOcc_log <- log(s_meanPPT_relPol$nOcc)

#run RF through boruta
boruta.meanPPT_relPol <- Boruta(slope_meanPPT_relPol ~
                                rangeSize_log10 +
                                rangeLoc +
                                roundness +
                                latAmplitude +
                                nOcc_log +
                                order +
                                elevMedian +
                                bodyMass +
                                elevAmplitude,
                              data = s_meanPPT_relPol, doTrace = 2, ntree = 500,
                              maxRuns = 500)


plot(boruta.meanPPT_relPol)
print(boruta.meanPPT_relPol)
attStats(boruta.meanPPT_relPol)


###  MAX  ###

#select only species that had lower correl between vars
s_maxPPT_relPol <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_relPol <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$Cor_vars_maxPPT),]
s_maxPPT_relPol <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_maxPPT_relPol),]


#transform categorical variable in factor
s_maxPPT_relPol$order <- as.factor(s_maxPPT_relPol$order)

#create new columns with log values (boruta does not accept it...)
s_maxPPT_relPol$rangeSize_log10 <- log(s_maxPPT_relPol$rangeSize, 10)
s_maxPPT_relPol$nOcc_log <- log(s_maxPPT_relPol$nOcc)

#run RF through boruta
boruta.maxPPT_relPol <- Boruta(slope_maxPPT_relPol ~
                               rangeSize_log10 +
                               rangeLoc +
                               roundness +
                               latAmplitude +
                               nOcc_log +
                               order +
                               elevMedian +
                               bodyMass +
                               elevAmplitude,
                             data = s_maxPPT_relPol, doTrace = 2, ntree = 500,
                             maxRuns = 500)


plot(boruta.maxPPT_relPol)
print(boruta.maxPPT_relPol)
attStats(boruta.maxPPT_relPol)




###### PRECIPITATION vs CENTRALNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minPPT_distEdge <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$Cor_vars_minPPT),]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$slope_minPPT_distEdge),]


#transform categorical variable in factor
s_minPPT_distEdge$order <- as.factor(s_minPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_minPPT_distEdge$rangeSize_log10 <- log(s_minPPT_distEdge$rangeSize, 10)
s_minPPT_distEdge$nOcc_log <- log(s_minPPT_distEdge$nOcc)

#run RF through boruta
boruta.minPPT_distEdge <- Boruta(slope_minPPT_distEdge ~
                                 rangeSize_log10 +
                                 rangeLoc +
                                 roundness +
                                 latAmplitude +
                                 nOcc_log +
                                 order +
                                 elevMedian +
                                 bodyMass +
                                 elevAmplitude,
                               data = s_minPPT_distEdge, doTrace = 2, ntree = 500,
                               maxRuns = 500)


plot(boruta.minPPT_distEdge)
print(boruta.minPPT_distEdge)
attStats(boruta.minPPT_distEdge)


###  MEAN  ###

#select only species that had lower correl between vars
s_meanPPT_distEdge <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$Cor_vars_meanPPT),]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$slope_meanPPT_distEdge),]


#transform categorical variable in factor
s_meanPPT_distEdge$order <- as.factor(s_meanPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_meanPPT_distEdge$rangeSize_log10 <- log(s_meanPPT_distEdge$rangeSize, 10)
s_meanPPT_distEdge$nOcc_log <- log(s_meanPPT_distEdge$nOcc)

#run RF through boruta
boruta.meanPPT_distEdge <- Boruta(slope_meanPPT_distEdge ~
                                  rangeSize_log10 +
                                  rangeLoc +
                                  roundness +
                                  latAmplitude +
                                  nOcc_log +
                                  order +
                                  elevMedian +
                                  bodyMass +
                                  elevAmplitude,
                         data = s_meanPPT_distEdge, doTrace = 2, ntree = 500,
                         maxRuns = 500)


plot(boruta.meanPPT_distEdge)
print(boruta.meanPPT_distEdge)
attStats(boruta.meanPPT_distEdge)


###  MAX  ###

#select only species that had lower correl between vars
s_maxPPT_distEdge <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$Cor_vars_maxPPT),]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$slope_maxPPT_distEdge),]


#transform categorical variable in factor
s_maxPPT_distEdge$order <- as.factor(s_maxPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_maxPPT_distEdge$rangeSize_log10 <- log(s_maxPPT_distEdge$rangeSize, 10)
s_maxPPT_distEdge$nOcc_log <- log(s_maxPPT_distEdge$nOcc)

#run RF through boruta
boruta.maxPPT_distEdge <- Boruta(slope_maxPPT_distEdge ~
                                 rangeSize_log10 +
                                 rangeLoc +
                                 roundness +
                                 latAmplitude +
                                 nOcc_log +
                                 order +
                                 elevMedian +
                                 bodyMass +
                                 elevAmplitude,
                               data = s_maxPPT_distEdge, doTrace = 2, ntree = 500,
                               maxRuns = 500)


plot(boruta.maxPPT_distEdge)
print(boruta.maxPPT_distEdge)
attStats(boruta.maxPPT_distEdge)



