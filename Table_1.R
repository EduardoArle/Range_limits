#set margin parametres for the plot
par(mar=c(3,7,3,7))

#make the empty plot
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
     xaxs = "i", yaxs = "i", axes = F, frame.plot = TRUE)

#make lines creating a table (cols)
lines(c(0, 0), c(0, 10), lwd = 2, col = '#808080')
lines(c(3.5, 3.5), c(0, 10), lwd = 1.5, col = '#808080')
lines(c(10, 10), c(0, 10), lwd = 2,, col = '#808080')

#make lines creating a table (rows)
lines(c(0, 10), c(0, 0), lwd = 2,, col = '#808080')
lines(c(0, 10), c(10/10, 10/10), col = '#808080')
lines(c(0, 10), c(10/10 *2, 10/10 *2), col = '#808080')
lines(c(0, 10), c(10/10 *3, 10/10 *3), col = '#808080')
lines(c(0, 10), c(10/10 *4, 10/10 *4), col = '#808080')
lines(c(0, 10), c(10/10 *5, 10/10 *5), col = '#808080')
lines(c(0, 10), c(10/10 *6, 10/10 *6), col = '#808080')
lines(c(0, 10), c(10/10 *7, 10/10 *7), col = '#808080')
lines(c(0, 10), c(10/10 *8, 10/10 *8), col = '#808080')
lines(c(0, 10), c(10/10 *9, 10/10 *9), lwd = 1.5, col = '#808080')
lines(c(0, 10), c(10, 10), lwd = 2, col = '#808080')

#define sizes
cex_att <- 0.7
cex_desc <- 0.5
cex_tit <- 0.7

#add text to the table
text(0.2, 0.5, 'Taxonomic order', cex = cex_att, adj = 0)
text(0.2, 10/10 + 0.5, 'N presence records', cex = cex_att, adj = 0)
text(0.2, (10/10 * 2) + 0.5, 'Body mass', cex = cex_att, adj = 0)
text(0.2, (10/10 * 3) + 0.5, 'Range shape', cex = cex_att, adj = 0)
text(0.2, (10/10 * 4) + 0.5, 'Elevation range', cex = cex_att, adj = 0)
text(0.2, (10/10 * 5) + 0.5, 'Elevation', cex = cex_att, adj = 0)
text(0.2, (10/10 * 6) + 0.5, 'Latitudinal position', cex = cex_att, adj = 0)
text(0.2, (10/10 * 7) + 0.5, 'Latitudinal range', cex = cex_att, adj = 0)
text(0.2, (10/10 * 8) + 0.5, 'Range size', cex = cex_att, adj = 0)
text(0.2, (10/10 * 9) + 0.5, 'Species-level attributes', cex = cex_tit, adj = 0, font = 2)



'Proxy for broad phylogenetic relatedness, which may influence species’ responses to climate (Burns and Strauss 2011). Retrieved during taxonomic harmonisation using the GBIF backbone taxonomy.'

# 'Burns and Strauss 2011' 111

text(3.6, 0.75, 'Proxy for broad phylogenetic relatedness, which may influence species’',
     cex = cex_desc, adj = 0)
text(3.6, 0.5, 'responses to climate¹¹¹. Retrieved during taxonomic harmonisation using the',
     cex = cex_desc, adj = 0)
text(3.6, 0.25, 'GBIF backbone taxonomy.',
     cex = cex_desc, adj = 0)


'Number of occurrences used to fit models (log), as model performance tends to increase with sample size (Guisan et al. 2017; Moudrý et al. 2024). Fewer records generally reduce the ability to detect clear niche patterns.'

# 'Guisan et al. 2017' 15
# 'Moudrý et al. 2024' 110

text(3.6, 10/10 + 0.75, 'Number of occurrences (log) used to fit models, as model performance', cex = cex_desc, adj = 0)
text(3.6, 10/10 + 0.5, 'tends to increase with sample size¹⁵˒¹¹⁰. Fewer records generally reduce ', cex = cex_desc, adj = 0)
text(3.6, 10/10 + 0.25, 'the ability to detect clear niche patterns.',
     cex = cex_desc, adj = 0)


"Proxy for thermal sensitivity (log) (Soria et al. 2021): smaller mammals are vulnerable to cold and temperature fluctuations; larger ones conserve heat but risk overheating (Pacifici et al. 2017; Spence and Tingley 2020)."

# 'Soria et al. 2021' 107
# 'Pacifici et al. 2017' 108
# 'Spence and Tingley 2020' 109

text(3.6, (10/10 * 2) + 0.75, expression('Proxy for thermal sensitivity¹⁰⁷ (log): smaller mammals are vulnerable to'),
     cex = cex_desc, adj = 0)
text(3.6, (10/10 * 2) + 0.5, 'cold and temperature fluctuations; larger ones conserve heat but risk', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 2) + 0.25, expression('overheating¹⁰⁸˒¹⁰⁹.'), cex = cex_desc, adj = 0)


'Measures roundness as (4π × Area) / Perimeter². Higher values indicate more circular ranges, lower values suggest irregular ranges shaped by geography, potentially weakening climatic gradients^106.'

# 'Csergő et al. 2020' 106

text(3.6, (10/10 * 3) + 0.75, 'Measures roundness as (4π × Area) / Perimeter². Higher values indicate', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 3) + 0.5, 'more circular ranges, lower values suggest irregular ranges shaped by', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 3) + 0.25, 'geography, potentially weakening climatic gradients¹⁰⁶.', cex = cex_desc, adj = 0)


'Calculated as the difference between the 97.5% and 2.5% elevation quantiles across presence records. Broader ranges suggest exposure to diverse climates, potentially leading to more variable responses.'


text(3.6, (10/10 * 4) + 0.75, 'Calculated as the difference between the 97.5% and 2.5% elevation', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 4) + 0.5, 'quantiles across presence records. Broader ranges suggest exposure to', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 4) + 0.25, 'diverse climates, potentially leading to more variable responses.', cex = cex_desc, adj = 0)


'Median altitude of species’ presence records, extracted from WorldClim using the terra R package. Higher elevations are expected to impose stronger climatic constraints due to more extreme conditions.'

text(3.6, (10/10 * 5) + 0.75, 'Median altitude of species’ presence records, extracted from WorldClim', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 5) + 0.5, 'using the terra R package. Higher elevations are expected to impose', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 5) + 0.25, 'stronger climatic constraints due to more extreme conditions.', cex = cex_desc, adj = 0)


'Absolute latitude of the midpoint of the species’ range. Temperate species are predicted to be more climatically constrained than tropical ones, consistent with prior studies^12,105'

# 'Addo-Bediako 2020' 12
# 'Percino-Daniel 2021' 105

text(3.6, (10/10 * 6) + 0.75, 'Absolute latitude of the midpoint of the species’ range. Temperate species are', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 6) + 0.5, 'predicted to be more climatically constrained than tropical ones, consistent', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 6) + 0.25, 'with prior studies¹²˒¹⁰⁵.', cex = cex_desc, adj = 0)


'Difference between the northern and southernmost limits of a species range. Wider ranges suggest exposure to broader climatic conditions, potentially indicating greater ecological tolerance^104 and clearer climatic gradients.'

# 'Stevens 1989' 104

text(3.6, (10/10 * 7) + 0.75, 'Difference between the northern and southernmost limits of a range.', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 7) + 0.5, 'Wider ranges suggest exposure to broader climatic conditions, possibly', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 7) + 0.25, 'indicating higher tolerance¹⁰⁴ and clearer climatic gradients.', cex = cex_desc, adj = 0)



'Distribution area (log10, km²). Larger ranges may show stronger climatic gradients^51 but may also reflect greater ecological tolerance and weaker climatic constraints^6.'

# 'Yancovitch-Shalom 2020' 51
# 'Davies 2009' 6

text(3.6, (10/10 * 8) + 0.75, 'Distribution area (log10, km²). Larger ranges may show stronger climatic', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 8) + 0.5, 'gradients⁵¹ but may also reflect greater ecological tolerance and weaker', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 8) + 0.25, 'climatic constraints⁶.', cex = cex_desc, adj = 0)

text(3.6, (10/10 * 9) + 0.5, 'Description', cex = cex_tit, adj = 0, font = 2)

#save PDF landscape 5 x 7.2 in

