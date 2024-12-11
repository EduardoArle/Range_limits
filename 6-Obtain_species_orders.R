library(taxize)

#list wds
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Thinned_occurrrences'
wd_lists <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'

#list species 
setwd(wd_thinned_occ)
sps_list <- gsub('_thinned.csv', '', list.files())

#get orders
order <- character()

for(i in 484:length(sps_list))
{
  a <- tax_name(sps_list[i], get = 'order', db = 'ncbi')
  order[i] <- a$order
  
  print(i)
}

#concatenate the sps_list and the order
sps_order <- as.data.frame(cbind(sps_list, order))
names(sps_order)[1] <-  'species'

#manually include the orders that were not founf
missing <- which(is.na(sps_order$order))

i = 1
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 2
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 3
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 4
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 5
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Artiodactyla'

i = 6
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 7
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 8
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 9
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 10
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 11
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Afrosoricida'

i = 12
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 13
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 14
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 15
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 16
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 17
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 18
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 19
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 20
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 21
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 22
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 23
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 24
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 25
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 26
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 27
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 28
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 29
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 30
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 31
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 32
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 33
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 34
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 35
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 36
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 37
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Rodentia'

i = 38
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 39
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 40
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 41
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Chiroptera'

i = 42
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 43
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 44
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 45
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 46
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 47
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 48
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 49
sps_order$species[missing[i]]
sps_order$order[missing[i]] <-  'Primates'

i = 50
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Primates'

i = 51
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Primates'

i = 52
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Primates'

i = 53
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Primates'

i = 54
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 55
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 56
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 57
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 58
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 59
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 60
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 61
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 62
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 63
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Chiroptera'

i = 64
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Eulipotyphla'

i = 65
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 66
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 67
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 68
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 69
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 70
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 71
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 72
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 73
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 74
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 75
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 76
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 77
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 78
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Didelphimorphia'

i = 79
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

i = 80
sps_order$species[missing[i]]
sps_order$order[missing[i]] <- 'Rodentia'

#write table with orders
setwd(wd_lists)
sps_order <- write.csv(sps_order, 'Species_order.csv',
                      row.names = F)

sort(unique(sps_order$order))
table((sps_order$order))
summary(sps_order$order)
