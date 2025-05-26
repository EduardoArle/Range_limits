#set margin parametres for the plot
par(mar=c(3,7,3,7))

#make the empty plot
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
     xaxs = "i", yaxs = "i", axes = F, frame.plot = TRUE)

#make lines creating a table (cols)
lines(c(0, 0), c(0, 10), lwd = 3)
lines(c(3.5, 3.5), c(0, 10), lwd = 2)
lines(c(10, 10), c(0, 10), lwd = 3)

#make lines creating a table (rows)
lines(c(0, 10), c(0, 0), lwd = 3)
lines(c(0, 10), c(10/10, 10/10))
lines(c(0, 10), c(10/10 *2, 10/10 *2))
lines(c(0, 10), c(10/10 *3, 10/10 *3))
lines(c(0, 10), c(10/10 *4, 10/10 *4))
lines(c(0, 10), c(10/10 *5, 10/10 *5))
lines(c(0, 10), c(10/10 *6, 10/10 *6))
lines(c(0, 10), c(10/10 *7, 10/10 *7))
lines(c(0, 10), c(10/10 *8, 10/10 *8))
lines(c(0, 10), c(10/10 *9, 10/10 *9), lwd = 2)
lines(c(0, 10), c(10, 10), lwd = 3)

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



'Proxy for phylogenetic relatedness, may influence species’ responses to climate (Burns & Strauss 2011). Retrieved via the taxize R package (Chamberlain & Szöcs 2013) and manually completed using IUCN Red List data (IUCN 2024).'

text(3.6, 0.75, 'Proxy for phylogenetic relatedness, may influence species’ responses',
     cex = cex_desc, adj = 0)
text(3.6, 0.5, 'to climate (Burns & Strauss 2011). Retrieved via the taxize R package',
     cex = cex_desc, adj = 0)
text(3.6, 0.25, 'and manually completed using IUCN Red List data (IUCN 2024).',
     cex = cex_desc, adj = 0)

'Number of occurrences used to fit models, as model performance tends to increase with sample size (Guisan et al. 2017; Moudrý et al. 2024). Fewer records generally reduce the ability to detect clear niche patterns.'

text(3.6, 10/10 + 0.75, 'Number of occurrence used to fit models, as model performance tends',
     cex = cex_desc, adj = 0)
text(3.6, 10/10 + 0.5, expression('to increase with sample size (Guisan' ~ italic('et al.') ~  '2017; Moudrý' ~ italic('et al.') ~ '2024).'),
     cex = cex_desc, adj = 0)
text(3.6, 10/10 + 0.25, 'Fewer records generally reduce the ability to detect clear niche patterns.',
     cex = cex_desc, adj = 0)


"Proxy for thermal sensitivity (Smith et al. 2003): smaller mammals are vulnerable to cold and temperature fluctuations; larger ones conserve heat but risk overheating (Pacifici et al. 2017; Spence and Tingley 2020)."


text(3.6, (10/10 * 2) + 0.75, expression('Proxy for thermal sensitivity (Smith' ~ italic('et al.') ~ '2003): smaller mammals are'),
     cex = cex_desc, adj = 0)
text(3.6, (10/10 * 2) + 0.5, 'vulnerable to cold and temperature fluctuations; larger ones conserve', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 2) + 0.25, expression('heat but risk overheating (Pacifici' ~ italic('et al.') ~ '2017; Spence and Tingley 2020).'), cex = cex_desc, adj = 0)


'Measures roundness as (4π × Area) / Perimeter². Higher values indicate more circular ranges, lower values suggest irregular ranges shaped by geography, potentially weakening climatic gradients (Csergő et al. 2020).'

text(3.6, (10/10 * 3) + 0.75, 'Measures roundness as (4π × Area) / Perimeter². Higher values indicate', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 3) + 0.5, 'more circular ranges, lower values suggest irregular ranges shaped by', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 3) + 0.25, expression('geography, potentially weakening climatic gradients (Csergő' ~ italic('et al.') ~ '2020).'), cex = cex_desc, adj = 0)


'Calculated as the difference between the 97.5% and 2.5% elevation quantiles across presence records. Broader ranges suggest exposure to diverse climates, potentially leading to more variable responses.'


text(3.6, (10/10 * 4) + 0.75, 'Calculated as the difference between the 97.5% and 2.5% elevation', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 4) + 0.5, 'quantiles across presence records. Broader ranges suggest exposure to diverse climates, potentially leading to more variable responses.', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 4) + 0.25, 'diverse climates, potentially leading to more variable responses.', cex = cex_desc, adj = 0)


'Median altitude of species’ presence records, extracted from WorldClim using the terra R package. Higher elevations are expected to impose stronger climatic constraints due to more extreme conditions.'

text(3.6, (10/10 * 5) + 0.75, 'Median altitude of species’ presence records, extracted from WorldClim', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 5) + 0.5, 'using the terra R package. Higher elevations are expected to impose', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 5) + 0.25, 'stronger climatic constraints due to more extreme conditions.', cex = cex_desc, adj = 0)


'Mean latitude of the range. Temperate species are predicted to be more climatically constrained than tropical ones, consistent with prior studies (Addo-Bediako et al. 2000; Stuart-Smith et al. 2017; Percino-Daniel et al. 2021).'

text(3.6, (10/10 * 6) + 0.75, 'Mean latitude of the range. Temperate species are predicted to be more', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 6) + 0.5, 'climatically constrained than tropical ones, consistent with prior studies', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 6) + 0.25, expression('(Addo-Bediako' ~ italic('et al.') ~ '2000; Percino-Daniel' ~ italic('et al.') ~ '2021).'), cex = cex_desc, adj = 0)


'Difference between the northern and southernmost limits of a species range. Wider ranges suggest exposure to broader climatic conditions, potentially indicating greater ecological tolerance (Stevens 1989) and clearer climatic gradients.'


text(3.6, (10/10 * 7) + 0.75, 'Difference between the northern and southernmost limits of a range.', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 7) + 0.5, 'Wider ranges suggest exposure to broader climatic conditions, possibly', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 7) + 0.25, 'indicating higher tolerance (Stevens 1989) and clearer climatic gradients.', cex = cex_desc, adj = 0)



'Distribution area (km²). Larger ranges may show stronger climatic gradients (Yancovitch-Shalom et al. 2020) but may also reflect greater ecological tolerance and weaker climatic constraints (Davies et al. 2009).'

text(3.6, (10/10 * 8) + 0.75, 'Distribution area (km²). Larger ranges may show stronger climatic', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 8) + 0.5, 'gradients (Yancovitch-Shalom et al. 2020) but may also reflect greater', cex = cex_desc, adj = 0)
text(3.6, (10/10 * 8) + 0.25, expression('ecological tolerance and weaker climatic constraints (Davies' ~ italic('et al.') ~ '2009).'), cex = cex_desc, adj = 0)

text(3.6, (10/10 * 9) + 0.5, 'Description', cex = cex_tit, adj = 0, font = 2)
