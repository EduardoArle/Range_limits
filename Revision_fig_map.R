####  This script plots the ranges of selected species showing geographical and taxonomic distribution os our sample

#load libraries
library(sf); library(terra); library(plyr)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_elevation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_tab_SI <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/SI/Tables'
wd_map_stuff <- '/Users/carloseduardoaribeiro/Documents/Map stuff/Carstens_world_proj'

#input mammals per order
taxa <- data.frame(
  order = c("Rodentia","Chiroptera","Primates","Eulipotyphla",
            "Artiodactyla","Carnivora","Didelphimorphia","Dasyuromorphia",
            "Lagomorpha","Afrosoricida","Macroscelidea","Cingulata",
            "Pilosa","Peramelemorphia","Hyracoidea","Proboscidea",
            "Notoryctemorphia","Tubulidentata"),
 n_species = c(2550,1460,520,500,360,305,125,80,110,55,20,21,10,30,5,2,2,1)
)

#load elevation raster (just for mapping)
setwd(wd_elevation)
elev <- rast('wc2.1_2.5m_elev.tif')

#list species
setwd(wd_tab_SI)
tab <- read.csv('Table_S1.csv')
sps_list <- tab$species

#count species number in orders represented
sps_order <- ddply(tab, .(order), nrow)

#create colour palette for orders
colours <- c("#6A3D9A", "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00",
             "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
             "#CAB2D6", "#FFFF99", "#8DD3C7")

#bar plot species order
ord <- order(sps_order$V1, decreasing = F)
counts <- sps_order$V1[ord]

ord_taxa <- order(taxa$n_species, decreasing = F)
counts_taxa <- taxa$n_species[ord_taxa]

#all orders from the global table
orders_all <- taxa$order
total <- taxa$n_species

#match my counts to the order vector
my_counts <- sps_order$V1[match(orders_all, sps_order$order)]

#order we don't have become NA -> set to 0
my_counts[is.na(my_counts)] <- 0

#sort by total richness
o <- order(total, decreasing = F)
orders_sorted <- orders_all[o]
total_sorted  <- total[o]
my_sorted     <- my_counts[o]
colors_sorted <- colors[o] 

#set margins for plot
par(mar = c(5, 7, 4, 2))

#plot background bars
bp <- barplot(total_sorted,
              names.arg = orders_sorted,
              horiz = TRUE,
              las = 1,
              cex.names = 0.7,
              col = "grey90",
              border = "grey40",
              xlab = "Number of species",
              main = "Species per order")

#build a vector of my counts, aligned to orders_sorted
my_counts <- sps_order$V1[match(orders_sorted, sps_order$order)]
my_counts[is.na(my_counts)] <- 0   # orders not in your study -> 0

#colours aligned the same way
cols_match <- colours[match(orders_sorted, sps_order$order)]

#overlay plots of our data
barplot(my_counts,
        horiz  = TRUE,
        col    = cols_match,
        border = "black",
        add    = TRUE,
        axes   = FALSE)

# barplot(sps_order$V1[ord],
#               names.arg = sps_order$order[ord],
#               horiz = T,
#               las = 1,
#               col = colours,
#               xlab = "Number of species",
#               ylab = NA,
#               main = "Species per order",
#               cex.names = 0.8)

# text(counts[c(1:7)] + 4, bp[c(1:7)], labels = counts[c(1:7)], cex = 0.8)
# text(counts[c(8:12)] - 6, bp[c(8:12)], labels = counts[c(8:12)], cex = 0.8)
# text(counts[13] - 9, bp[13], labels = counts[13], cex = 0.8)

text(30, bp[c(1:5)], labels = my_counts[c(1:5)], cex = 0.7,adj = 0)
text(40, bp[c(6:8)], labels = my_counts[c(6:8)], cex = 0.7,adj = 0)
text(20, bp[c(9:15)], labels = my_counts[c(9:14)], cex = 0.7,adj = 0)
text(my_counts[16] + 20, bp[16], labels = my_counts[c(16)], cex = 0.7,adj = 0)
text(my_counts[17] + 20, bp[17], labels = my_counts[17], cex = 0.7,adj = 0)
text(my_counts[18] + 20, bp[18], labels = my_counts[18], cex = 0.7,adj = 0)




### Plot maps ###

#make table relating order to colour code (following barplot scheme)
ord_col <- data.frame(order = sps_order$order[ord], colour = colours)

#include colour code for each species based on order
tab_2 <- merge(tab, ord_col, by = 'order', sort = F)

#load map stuff
setwd(wd_map_stuff)
world <- readRDS("wrld.rds")
worldmapframe <- readRDS("Worldmapframe.rds")

#transform SpatialPolygons into sf object
worldmapframe <- st_as_sf(worldmapframe)

#change the crs of the frame and the elev raster
worldmapframe <- st_transform(worldmapframe, crs = st_crs(world))
elev <- project(elev, world)

#set margins for plot
par <- par(mar = (c(5,6,3,3)))

#plot world
plot(elev, box = NA, legend = NA, axes = NA, col = 'grey90')

#add frame
plot(worldmapframe, add = T)

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #loas species occurrences
  setwd(wd_thinned_occ)
  sps_occ <- read.csv(paste0(sps, '_thinned.csv'))
  
  #select only presence occurrences
  pr_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'PRESENT'),]
  
  #select only columns we are interested in 
  pr_sps2 <- pr_sps[,c('species','key','decimalLongitude','decimalLatitude',
                       'datasetKey')]
  
  #include number of records
  pr_sps2$nOcc <- nrow(pr_sps2)
  
  #check if there are absence data
  abs_sps <- sps_occ[which(sps_occ$occurrenceStatus == 'ABSENT'),]
  if(nrow(abs_sps) != 0){
    warning(paste0('THERE IS ABSENT DATA FOR ', sps_list[i]))
  }
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps_list[i])
  
  #select the range representing only the native range of the species
  range <- range[range$legend == 'Extant (resident)',]
  
  #unify all features
  sps_range2 <- st_union(range)
  
  #change the crs of the frame and the elev raster
  sps_range2 <- st_transform(sps_range2, crs = st_crs(world))
  
  #plot sps range map
  plot(sps_range2, add = T, col = paste0(tab_2$colour[i], '80'),
       border = tab_2$colour[i])
  
}

