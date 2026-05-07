#list wds
wd_table <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta'


#read table with Boruta results
setwd(wd_table)
tab <- read.csv('boruta_precipitation_summary.csv')

#delete empty rows
tab <- tab[tab[,1] != "", ]

#add row names and delete first col
row.names(tab) <- tab[,1]
tab <- tab[,-1]

#create table with colours
tab_col <- tab
tab_col[tab_col < 1] <- '#ffffff'
tab_col[tab_col > 0 & tab_col < 5] <- '#00805070'
tab_col[tab_col == 5] <- '#008050'

#rowindices
rn <- row.names(tab_col)
i_rel_eq <- grep("relPol.*\\bEQ\\b", rn, ignore.case = TRUE) #indrelEQ
i_rel_pol <- grep("relPol.*\\bPOL\\b", rn, ignore.case = TRUE) #indrelPOL
i_dist <- grep("distEdge", rn, ignore.case = TRUE) #inddist

#subtables
tab_rel_eq <- tab_col[i_rel_eq,] 
tab_rel_pol <- tab_col[i_rel_pol,]
tab_dist <- tab_col[i_dist,]

#set margins
par(mar=c(15,10,3,1),
    oma=c(0,0,0,0),
    mfrow=c(1,1))

#make the empty plot
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 5),
     xaxs = "i",yaxs = "i", axes=F, frame.plot=TRUE)

#decide how many rows and cols the table needs
rows <- 6 #nrowsplot
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

#fillrectangles

#cell size in plot coordinates
cell_w <- 10/cols
cell_h <- 5/rows

#define margin and gap between split ranges
m <- 0.15 #tilemargin
g <- 0.10 #splitgap

#tileside
tile_side <- min(cell_w, cell_h) * (1 - 2*m) * 0.85 #tileshrink

#draw full cell rectangle
drawfull <- function(irow, jcol, fillcol)
{
  x0 <- (jcol-1) * cell_w
  x1 <- jcol * cell_w
  y0 <- (irow-1) * cell_h
  y1 <- irow * cell_h
  
  #tilecoords
  cx <- (x0 + x1)/2
  cy <- (y0 + y1)/2
  xL <- cx - tile_side/2
  xR <- cx + tile_side/2
  yB <- cy - tile_side/2
  yT <- cy + tile_side/2
  
  #tile
  rect(xL, yB, xR, yT, col = fillcol, border = "black")
}

#draw split cell rectangle (POL top, EQ bottom)
drawsplit <- function(irow, jcol, col_pol, col_eq)
{
  x0 <- (jcol-1) * cell_w
  x1 <- jcol * cell_w
  y0 <- (irow-1) * cell_h
  y1 <- irow * cell_h
  ym <- (y0 + y1)/2
  
  #tilecoords
  cx <- (x0 + x1)/2
  cy <- (y0 + y1)/2
  xL <- cx - tile_side/2
  xR <- cx + tile_side/2
  yB <- cy - tile_side/2
  yT <- cy + tile_side/2
  
  #splitcoords
  ym2 <- (yB + yT)/2
  gap <- g * tile_side
  
  #tophalf POL
  rect(xL, ym2 + gap/2, xR, yT, col = col_pol, border = "black")
  
  #bottomhalf EQ
  rect(xL, yB, xR, ym2 - gap/2, col = col_eq, border = "black")
}

#fill by blocks
for(j in 1:cols)
{
  #distEdge block (bottom 3 rows)
  for(k in 1:3)
  {
    irow <- k
    drawfull(irow, j, tab_dist[4-k, j])
  }
  
  #relPol block (top 3 rows)
  for(k in 1:3)
  {
    irow <- 3 + k
    drawsplit(irow, j,
              tab_rel_pol[4-k, j],  #POL top
              tab_rel_eq[4-k,  j])  #EQ bottom
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

#ylabels
ypos <- seq(5/rows/2,(5/rows/2)+(5/rows*(rows-1)),by = 5/rows)

text(par("usr")[3]-0.2,
     ypos,
     adj = 1, xpd = TRUE,
     labels = c("max","mean","min","max","mean","min"),
     cex = 1.3)

#vertical text, grouping the two row blocks
text(par("usr")[3]-1.4,
     c(1.25, 3.75),
     srt = 90, adj = 0.5, xpd = TRUE,
     labels = c('distEdge', 'relPol'), cex = 1.6)

#vertical text, full response variable label
text(par("usr")[3]-2.2,
     2.5,
     srt = 90, adj = 0.5, xpd = TRUE,
     labels = 'PRECIPITATION', cex = 1.6)

#save 1000 width, 618 height