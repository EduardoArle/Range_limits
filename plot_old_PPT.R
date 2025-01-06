setwd('/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models_PPT')

minPPT_relPolewarness_G <- readRDS('minPPT_relPolewarness_G')

summary(minPPT_relPolewarness_G)

plot.gam(minPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800



minPPT_relPolewarness_GS <- readRDS('minPPT_relPolewarness_GS')

summary(minPPT_relPolewarness_GS)

plot.gam(minPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800





meanPPT_relPolewarness_G <- readRDS('meanPPT_relPolewarness_G')

summary(meanPPT_relPolewarness_G)

plot.gam(meanPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



meanPPT_relPolewarness_GS <- readRDS('meanPPT_relPolewarness_GS')

summary(meanPPT_relPolewarness_GS)

plot.gam(meanPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800






maxPPT_relPolewarness_G <- readRDS('maxPPT_relPolewarness_G')

summary(maxPPT_relPolewarness_G)

plot.gam(maxPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



maxPPT_relPolewarness_GS <- readRDS('maxPPT_relPolewarness_GS')

summary(maxPPT_relPolewarness_GS)

plot.gam(maxPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800





minPPT_absPolewarness_G <- readRDS('minPPT_absPolewarness_G')

summary(minPPT_absPolewarness_G)

plot.gam(minPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800



minPPT_absPolewarness_GS <- readRDS('minPPT_absPolewarness_GS')

summary(minPPT_absPolewarness_GS)

plot.gam(minPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800





meanPPT_absPolewarness_G <- readRDS('meanPPT_absPolewarness_G')

summary(meanPPT_absPolewarness_G)

plot.gam(meanPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800



meanPPT_absPolewarness_GS <- readRDS('meanPPT_absPolewarness_GS')

summary(meanPPT_absPolewarness_GS)

plot.gam(meanPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.09, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800






maxPPT_absPolewarness_G <- readRDS('maxPPT_absPolewarness_G')

summary(maxPPT_absPolewarness_G)

plot.gam(maxPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#1a985030', ylab = 'SHAP value',
         ylim = c(-0.25, 0.12),
         cex.lab = 2, cex.axis = 1.5) #save 800



