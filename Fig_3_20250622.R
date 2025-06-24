#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Slopes'
wd_table <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Figures/Fig 3'

############################################################################# 
################################# Panel A ###################################
############################################################################# 

#read table with boruta results
setwd(wd_table)
tab <- read.csv('Covariable_importance.csv')

#add row names and delete first col
row.names(tab) <- tab[,1]
tab <- tab[,-1]

#create table with colours (option 1)
tab_col <- tab
tab_col[tab_col > -5 & tab_col < 5] <- '#ffc00c'
tab_col[tab_col == -5] <- '#990050'
tab_col[tab_col == 5] <- '#008050'

#create table with colours (option 2)
tab_col <- tab
tab_col[tab_col > -5 & tab_col < 5] <- '#00805070'
tab_col[tab_col == -5] <- '#ffffff'
tab_col[tab_col == 5] <- '#008050'

#create table with colours (option 3)
tab_col <- tab
tab_col[tab_col > -5 & tab_col < 5] <- '#ffc00c'
tab_col[tab_col == -5] <- '#ffffff'
tab_col[tab_col == 5] <- '#008050'

#restrict tab only to T results
tab_col <- tab_col[c(1:6),]

#set margins
par(mar=c(15,10,3,1))

#make the empty plot
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 5),
     xaxs = "i",yaxs = "i", axes=F, frame.plot=TRUE)

#decide how many rows and cols the table needs
rows <- nrow(tab_col)
cols <- ncol(tab_col)

#make lines creating a table (cols)
for(i in 1:(cols-1))
{
  a <- c(i*10/cols,i*10/cols,i*10/cols)
  b <- c(0,5,10)
  lines(a,b)
}

#make lines creating a table (rows)
for(i in 1:(rows-1))
{
  a <- c(0,5,10)
  b <- c(i*5/rows,i*5/rows,i*5/rows)
  lines(a,b)
}

#plot squares with colours in the results table
for(i in 1:rows)
{
  for(j in 1:cols)
  {
    points((10/cols/2)+(10/cols*(j-1)),
           (5/rows/2)+(5/rows*(i-1)),
           bg = tab_col[rows - i + 1, j],
           pch = 22, cex = 5)
  }
}

#add x axis
axis(side = 1, 
     at = seq(10/cols/2,(10/cols/2)+(10/cols*(cols-1)),by = 10/cols),
     labels = NA, cex.axis = .8, padj = 0, las =2)

#add x labels, rotate 35 degrees (srt)
text(seq(10/cols/2,(10/cols/2)+(10/cols*(cols-1)),by = 10/cols), 
     par("usr")[3]-0.35, 
     srt = 35, adj = 1, xpd = TRUE,
     labels = gsub('\\.', ' ',names(tab)), cex = 1.3)

#add y axis
# axis(side = 2, 
#      at = seq(10/rows/2,(10/rows/2)+(10/rows*(rows-1)),by = 10/rows),
#      labels = NA, cex.axis = 1, padj = 0, las =1)

#add y labels, in three parts

#horizontal text, 1 value per row 'min' , 'mean' , 'max'
text(par("usr")[3]-0.2, 
     seq(5/rows/2,(5/rows/2)+(5/rows*(rows-1)),by = 5/rows), 
     adj = 1, xpd = TRUE,
     labels = rev(gsub("^(.*?)[A-Z].*$", "\\1", row.names(tab_col))), cex = 1.3)

#vertical text, 1 value grouping every 3 rows with 'min' , 'mean' , 'max'
text(par("usr")[3]-1.4,
     seq(5/rows/2,(5/rows/2)+(5/rows*(rows-1)),by = 5/rows)[c(3,6)] +
       c(-1, -3.2)  * 1/rows,  #position adjustment term
     srt = 90, adj = 1, xpd = TRUE,
     labels = c('distEdge', 'relPol'), cex = 1.6)

#vertical text, 1 value grouping 6 rows, temperature and precipitatio 
text(par("usr")[3]-2.2,
     seq(5/rows/2,(5/rows/2)+(5/rows*(rows-1)),by = 5/rows)[5] +
       -0.8  * 1/rows,  #position adjustment term
     srt = 90, adj = 1, xpd = TRUE,
     labels = c('TEMPERATURE'), cex = 1.6)

#save 1000

setwd(wd_slopes)
slopes <- read.csv('Slopes.csv')

#set parametres for plotting
par(mar = c(7,9,5,7), pty="m", mfrow = c(1,1))

###### MinT vs RELATIVE POLEWARDNESS ######

#select only species that had lower correl between vars
s_minT_relPol <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$Cor_vars_minT),]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$slope_minT_relPol),]

#create new columns with log values (boruta does not accept it...)
s_minT_relPol$rangeSize_log10 <- log(s_minT_relPol$rangeSize, 10)
s_minT_relPol$nOcc_log <- log(s_minT_relPol$nOcc)

#plot meaningful variables against slope


############################################################################# 
################################# Panel B ###################################
############################################################################# 

## Range size ##

#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$rangeSize_log10)
y_lim <- c(-0.05, 0.2)

#plot graph
plot(s_minT_relPol$rangeSize_log10, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

# Define original values and their log-transformed positions
range_vals <- max(s_minT_relPol$rangeSize_log10) -
  min(s_minT_relPol$rangeSize_log10)

plotting_positions <- c(min(s_minT_relPol$rangeSize_log10),
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/2,
                        min(s_minT_relPol$rangeSize_log10) + range_vals/4*3,
                        max(s_minT_relPol$rangeSize_log10))

plotting_values <- round(10 ^ plotting_positions / 1000)

#add axes
axis(1, at = plotting_positions, labels = plotting_values,
     pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Range size', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$rangeSize_log10
x_range <- range(x_vals)

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(min(plotting_positions) + range_vals/100,
             max(plotting_positions) - range_vals/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)

#save 800

############################################################################# 
################################# Panel C ###################################
############################################################################# 

## Latitudinal amplitude ##


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$latAmplitude)
y_lim <- c(-0.05, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$latAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Latitudinal amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$latAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#dd shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)

#save 800



############################################################################# 
################################# Panel D ###################################
############################################################################# 

### Elevation median ###


#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevMedian)
y_lim <- c(-0.05, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevMedian, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)


#add axes lables 
mtext('Elevation median', side = 1, line = 3.8, cex = 2.5)  
mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevMedian
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)

#save 800


############################################################################# 
################################# Panel E ###################################
############################################################################# 


## Elevational amplitude ##
#restricted ylim weighted by nOcc

# Define x and y limits
x_lim <- range(s_minT_relPol$elevAmplitude)
y_lim <- c(-0.05, 0.2)

# Ensure x_lim starts at 0 for a clean intersection
x_lim[1] <- 0 

#plot graph
plot(s_minT_relPol$elevAmplitude, s_minT_relPol$slope_minT_relPol,
     axes = F, xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", cex = 1.5,
     pch = 19, col = '#50305020',
     ylim = y_lim, xlim = x_lim)

#add axes
axis(1, pos = y_lim[1], cex.axis = 2)
axis(2, pos = x_lim[1], las=2, cex.axis = 2)

#add axes lables 
mtext('Elevational amplitude', side = 1, line = 3.8, cex = 2.5)  
#mtext('Slope', side = 2, line = 6.5, cex = 2.5)

#define x range explicitly
x_vals <- s_minT_relPol$elevAmplitude
x_range <- range(x_vals)
range_val <- x_range[2] - x_range[1]

#fit linear model
lin_mod_minT <- lm(s_minT_relPol$slope_minT_relPol ~ x_vals,
                   weights = s_minT_relPol$nOcc_log)

#predict y-values only for positive x-values
x_seq <- seq(range_val/100, max(x_vals) - range_val/100,
             length.out = 100)  # Ensuring it starts at 0
y_pred <- predict(lin_mod_minT, newdata = data.frame(x_vals = x_seq),
                  interval = "confidence")

#add regression line from the y-axis onwards
lines(x_seq, y_pred[, "fit"], col = '#503050', lwd = 8)

#add shaded confidence interval
polygon(c(x_seq, rev(x_seq)),
        c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = '#50305030',
        border = NA)


