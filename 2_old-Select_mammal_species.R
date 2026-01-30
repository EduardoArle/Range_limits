#load packages
library(sf); library(rworldmap)

#list wds
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_lists <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'

#list species
setwd(wd_ranges)
sps_list <- gsub('.shp', '', list.files(pattern = '.shp'))

#load world map
world <- getMap()

#visualise species distributions and choose the ones that are continuous and do not border the oceans (manually)

i = 5620

range <- st_read(dsn = wd_ranges, layer = sps_list[i])

plot(st_geometry(range), add = F, col = 'red')
plot(world, border = NA, col = 'darkgreen')
plot(st_geometry(range), add = T, col = 'orange', border = NA)

sps_list[i]



#list species that may be interesting to look at

interesting_sps <- c(3, 8, 13, 14, 38, 41, 43, 44, 45, 59, 68, 69, 70, 71,
                     72, 83, 85, 87, 89, 91, 92, 93, 98, 100, 102, 107, 108, 112,
                     113, 114, 123, 129, 130, 132, 138, 144, 147, 148, 152, 155,
                     156, 159, 166, 167, 168, 180, 186, 187, 188, 195, 225, 228,
                     230, 231, 232, 233, 237, 240, 246, 248, 308, 311, 326, 330,
                     332, 337, 338, 339, 395, 401, 409, 412, 420, 460, 465, 466,
                     475, 480, 487, 499, 502, 504, 506, 507, 510, 513, 515, 516,
                     517, 518, 519, 520, 521, 531, 542, 555, 558, 565, 568, 576,
                     577, 579, 580, 589, 590, 592, 597, 601, 604, 605, 656, 658,
                     659, 660, 662, 669, 692, 699, 713, 729, 731, 732, 733, 737,
                     739, 749, 757, 761, 765, 766, 778, 811, 819, 827, 848, 850,
                     870, 881, 882, 911, 934, 939, 950, 965, 970, 974, 1018,
                     1035, 1041, 1043, 1067, 1078, 1096, 1097, 1099, 1149, 1151,
                     1152, 1155, 1156, 1157, 1161, 1162, 1163, 1164, 1169, 1175,
                     1176, 1178, 1179, 1180, 1184, 1185, 1186, 1187, 1188, 1190,
                     1192, 1194, 1197, 1202, 1204, 1218, 1219, 1220, 1221, 1232, 
                     1233, 1243, 1244, 1248, 1262, 1266, 1269, 1315, 1337, 1354, 
                     1371, 1372, 1374, 1375, 1377, 1399, 1409, 1427, 1428, 1443, 
                     1445, 1447, 1448, 1449, 1455, 1456, 1457, 1476, 1490, 1495, 
                     1496, 1498, 1499, 1501, 1502, 1505, 1531, 1534, 1553, 1569, 
                     1573, 1606, 1614, 1621, 1622, 1623, 1625, 1636, 1642, 1647, 
                     1654, 1655, 1665, 1666, 1673, 1682, 1683, 1687, 1724, 1725, 
                     1734, 1752, 1753, 1757, 1763, 1770, 1785, 1789, 1792, 1804,
                     1810, 1821, 1848, 1849, 1853, 1866, 1885, 1886, 1887, 1912, 
                     1917, 1963, 2005, 2028, 2050, 2064, 2065, 2069, 2091, 2095, 
                     2116, 2117, 2153, 2164, 2165, 2166, 2167, 2174, 2178, 2207, 
                     2214, 2220, 2221, 2238, 2256, 2265, 2271, 2272, 2274, 2275, 
                     2280, 2284, 2285, 2287, 2288, 2289, 2290, 2291, 2292, 2293, 
                     2294, 2297, 2339, 2357, 2361, 2369, 2378, 2380, 2382, 2408, 
                     2409, 2410, 2411, 2429, 2431, 2432, 2443, 2445, 2454, 2462, 
                     2476, 2490, 2520, 2543, 2549, 2553, 2556, 2562, 2566, 2567, 
                     2568, 2575, 2578, 2580, 2601, 2604, 2626, 2687, 2688, 2702, 
                     2703, 2704, 2718, 2722, 2724, 2725, 2726, 2727, 2728, 2729, 
                     2730, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2741, 2742, 
                     2816, 2818, 2827, 2833, 2838, 2861, 2875, 2935, 2950, 2952, 
                     2963, 2964, 2971, 2989, 2990, 2991, 2992, 3044, 3049, 3050, 
                     3061, 3070, 3072, 3100, 3120, 3121, 3122, 3131, 3172, 3179,
                     3184, 3186, 3240, 3262, 3286, 3313, 3314, 3317, 3321, 3327, 
                     3328, 3332, 3335, 3336, 3343, 3354, 3363, 3365, 3366, 3367, 
                     3377, 3378, 3379, 3395, 3397, 3398, 3400, 3406, 3414, 3415, 
                     3430, 3432, 3436, 3449, 3478, 3479, 3509, 3565, 3567, 3568, 
                     3569, 3574, 3575, 3576, 3577, 3578, 3581, 3582, 3584, 3587, 
                     3588, 3594, 3595, 3601, 3605, 3606, 3608, 3611, 3614, 3618, 
                     3620, 3621, 3623, 3628, 3631, 3632, 3634, 3635, 3639, 3643, 
                     3678, 3684, 3688, 3694, 3697, 3716, 3721, 3724, 3727, 3729, 
                     3753, 3805, 3844, 3850, 3851, 3876, 3888, 3893, 3921, 3922, 
                     3926, 3962, 3963, 3968, 3997, 4000, 4039, 4042, 4048, 4050,
                     4051, 4059, 4060, 4061, 4062, 4066, 4068, 4085, 4107, 4108, 
                     4109, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117, 4118, 
                     4119, 4121, 4122, 4124, 4128, 4132, 4142, 4143, 4145, 4152, 
                     4157, 4159, 4164, 4168, 4169, 4171, 4172, 4173, 4175, 4176, 
                     4177, 4178, 4179, 4180, 4181, 4182, 4184, 4187, 4188, 4189,
                     4190, 4191, 4192, 4193, 4194, 4211, 4231, 4266, 4272, 4277, 
                     4278, 4282, 4283, 4284, 4286, 4287, 4288, 4291, 4292, 4362,
                     4364, 4384, 4450, 4470, 4528, 4569, 4577, 4580, 4594, 4603, 
                     4606, 4616, 4621, 4635, 4656, 4659, 4561, 4683, 4695, 4697, 
                     4701, 4702, 4706, 4709, 4711, 4712, 4715, 4722, 4775, 4776, 
                     4778, 4779, 4780, 4782, 4787, 4788, 4792, 4795, 4797, 4803, 
                     4813, 4819, 4821, 4833, 4840, 4844, 4847, 4849, 4850, 4851, 
                     4853, 4854, 4859, 4860, 4861, 4877, 4895, 4901, 4905, 4919,
                     4933, 4951, 4969, 4970, 4974, 4977, 4980, 4988, 5001, 5002, 
                     5007, 5008, 5016, 5029, 5030, 5040, 5050, 5058, 5067, 5068, 
                     5084, 5089, 5095, 5112, 5116, 5124, 5142, 5193, 5196, 5199, 
                     5211, 5292, 5295, 5296, 5298, 5304, 5307, 5308, 5319, 5320, 
                     5321, 5322, 5323, 5324, 5325, 5326, 5327, 5328, 5329, 5330, 
                     5331, 5333, 5335, 5337, 5339, 5340, 5343, 5344, 5345, 5347, 
                     5348, 5349, 5350, 5351, 5352, 5354, 5355, 5356, 5361, 5371,
                     5373, 5375, 5376, 5379, 5380, 5382, 5383, 5385, 5386, 5403, 
                     5404, 5410, 5431, 5499, 5502, 5504, 5505, 5508, 5511, 5514, 
                     5518, 5536, 5585, 5600, 5616, 5617, 5620)


## list some chosen species (for now)
#########

sel_sps_indices <- c(13, 69, 85, 87, 92, 93, 98, 107, 108, 112, 113, 147, 167,
                     180, 195, 228, 233, 237, 248, 326, 337, 338, 339, 409, 487, 
                     504, 506, 513, 515, 516, 518, 519, 520, 590, 605, 656, 713,
                     870, 881, 950, 965, 970, 974, 1175, 1176, 1179, 1202, 1204,
                     1218, 1219)

sel_species <- c('Abrothrix andinus', 'Aethomys nyikae', 'Akodon dayi', 'Akodon fumeus', 'Akodon lutescens', 'Akodon mimus', 'Akodon orophilus', 'Akodon simulator', 'Akodon spegazzinii', 'Akodon toba', 'Akodon torques', 'Alouatta caraya', 'Alticola semicanus', 'Ammospermophilus interpres', 'Anoura aequatoris', 'Aotus lemurinus', 'Aotus vociferans', 'Apodemus alpicola', 'Apodemus pallipes', 'Ateles chamek', 'Auliscomys boliviensis', 'Auliscomys pictus', 'Auliscomys sublimis', 'Blarina hylophaga', 'Callithrix penicillata', 'Callospermophilus madrensis', 'Calomys boliviae', 'Calomys lepidus', 'Calomys sorellus', 'Calomys tener', 'Calomys venustus', 'Calomyscus bailwardi', 'Calomyscus baluchi', 'Cebuella pygmaea', 'Cebus versicolor', 'Cercopithecus wolfi', 'Chalinolobus picatus', 'Cratogeomys goldmani', 'Cricetulus longicaudatus', 'Crocidura hildegardeae', 'Crocidura latona', 'Crocidura littoralis', 'Crocidura luna', 'Ctenomys maulinus', 'Ctenomys mendocinus', 'Ctenomys opimus', 'Ctenomys yolandae', 'Cuniculus taczanowskii', 'Cynomys gunnisoni', 'Cynomys leucurus')


sel_sps_indices <- c(1220, 1221, 1232, 1233, 1243, 1244, 1248, 1262, 1266, 1269,
                     1315, 1337, 1354, 1371, 1372, 1375, 1377, 1399, 1409, 1427,
                     1428, 1443, 1445, 1447, 1448, 1449, 1455, 1456, 1457, 1476,
                     1490, 1495, 1496, 1498, 1499, 1501, 1502, 1505, 1531, 1534,
                     1553, 1569, 1573, 1606, 1614, 1621, 1622, 1623, 1625, 1636)

sel_species <- c("Cynomys ludovicianus", "Cynomys mexicanus", "Dactylomys boliviensis","Dactylomys dactylinus", "Dasycercus blythi", "Dasycercus cristicauda", "Dasymys nudipes", "Dasyprocta variegata", "Dasypus pilosus", "Dasypus yepesi", "Dermanura anderseni", "Diclidurus isabella","Dinomys branickii", "Dipodomys microps", "Dipodomys nelsoni","Dipodomys ornatus", "Dipodomys phillipsii", "Dolichotis salinicola","Dremomys gularis", "Echimys saturnus", "Echimys vieirai", "Elephantulus brachyrhynchus", "Elephantulus fuscipes","Elephantulus intufi", "Elephantulus myurus", "Elephantulus pilicaudus", "Eligmodontia moreni", "Eligmodontia morgani", "Eligmodontia puerulus", "Ellobius lutescens", "Eolagurus przewalskii", "Eospalax rothschildi", "Eospalax smithii", "Eothenomys chinensis", "Eothenomys custos", "Eothenomys miletus", "Eothenomys olitor", "Eozapus setchuanus", "Eptesicus floweri", "Eptesicus gobiensis", "Equus kiang", "Euchoreutes naso", "Eudorcas albonotata", "Euneomys mordax", "Euroscaptor grandis", "Euryoryzomys emmonsae", "Euryoryzomys lamia", "Euryoryzomys legatus", "Euryoryzomys nitidus", "Felis bieti")



sel_sps_indices <- c(1642, 1647, 1654, 1655, 1665, 1666, 1673, 1682, 1683, 1687,
                     1724, 1725, 1734, 1752, 1753, 1757, 1763, 1770, 1785, 1789,
                     1792, 1804, 1810, 1821, 1848, 1849, 1853, 1866, 1885, 1886,
                     1887, 1912, 1917, 2005, 2028, 2050, 2064, 2065, 2069, 2091, 
                     2095, 2116, 2117, 2153, 2164, 2165, 2166, 2167, 2174, 2178)

sel_species <- c("Felovia vae", "Fukomys damarensis", "Fukomys mechowii", "Fukomys ochraceocinereus", "Funisciurus bayonii", "Funisciurus carruthersi", "Funisciurus substriatus", "Galea comes", "Galea flavidens", "Galenomys garleppi", "Geomys bursarius", "Geomys knoxjonesi", "Gerbilliscus boehmi", "Gerbillus aquilus", "Gerbillus bottai", "Gerbillus cosensis", "Gerbillus gleadowi", "Gerbillus juliani", "Gerbillus pulvinatus", "Gerbillus rupicola", "Gerbillus stigmonyx", "Glauconycteris humeralis", "Glironia venusta", "Glyphonycteris behnii", "Graomys chacoensis", "Graomys domorum", "Graphiurus christyi", "Graphiurus walterverheyeni", "Handleyomys chapmani", "Handleyomys fuscatus", "Handleyomys intectus", "Heliosciurus ruwenzorii", "Hemicentetes nigriceps", "Hipposideros khasiana", "Hipposideros rotalis", "Hoolock leuconedys", "Hyladelphys kalinowskii", "Hylaeamys acritus", "Hylaeamys perenensis", "Hylomyscus denniae", "Hylomyscus kerbispeterhansi", "Hyperacrius fertilis", "Hyperacrius wynnei", "Idionycteris phyllotis", "Isothrix bistriata", "Isothrix negrensis", "Isothrix orinoci", "Isothrix pagurus", "Jaculus thaleri", "Juscelinomys huanchacae")



sel_sps_indices <- c(2207, 2214, 2220, 2221, 2238, 2256, 2265, 2271, 2272, 2274,
                     2275, 2280, 2284, 2285, 2287, 2288, 2289, 2290, 2291, 2292,
                     2293, 2294, 2297, 2339, 2357, 2361, 2369, 2378, 2380, 2382,
                     2408, 2409, 2410, 2411, 2429, 2431, 2432, 2443, 2445, 2454,
                     2462, 2476, 2490, 2520, 2543, 2549, 2553, 2556, 2562, 2566)

sel_species <- c("Kerodon acrobata", "Kunsia tomentosus", "Lagidium viscacia", "Lagidium wolffsohni", "Lasiopodomys brandtii", "Lasiurus salinae", "Lemmiscus curtatus", "Lemniscomys bellieri", "Lemniscomys griselda", "Lemniscomys linulus", "Lemniscomys macculus", "Lemniscomys zebra", "Lenoxus apicalis", "Leontocebus cruzlimai", "Leontocebus fuscus", "Leontocebus illigeri", "Leontocebus lagonotus", "Leontocebus leucogenys", "Leontocebus nigricollis", "Leontocebus nigrifrons", "Leontocebus tripartitus", "Leontocebus weddelli", "Leontopithecus chrysopygus", "Leporillus apicalis", "Lepus comus", "Lepus fagani", "Lepus oiostolus", "Lepus townsendii", "Lepus yarkandensis", "Lestoros inca", "Lonchorhina inusitata", "Lonchorhina marinkellei", "Lonchorhina orinocensis", "Lonchothrix emiliae", "Lophuromys huttereri", "Lophuromys luteogaster", "Lophuromys machangui", "Lophuromys simensis", "Lophuromys woosnami", "Lutreolina massoia", "Lyncodon patagonicus", "Macaca leucogenys", "Macaca thibetana", "Makalata macrura", "Marmosa constantiae", "Marmosa phaea", "Marmosa rubra", "Marmosops bishopi", "Marmosops impavidus", "Marmosops neblina")


sel_sps_indices <- c(2567, 2568, 2575, 2578, 2580, 2601, 2604, 2626, 2687, 2688,
                     2702, 2703, 2704, 2718, 2722, 2724, 2725, 2726, 2727, 2728,
                     2729, 2730, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2741,
                     2742, 2816, 2818, 2827, 2833, 2838, 2861, 2875, 2935, 2950,
                     2952, 2963, 2964, 2971, 2989, 2990, 2991, 2992, 3044, 3049)

sel_species <- c("Marmosops noctivagus", "Marmosops ocellatus", "Marmota broweri", "Marmota caudata", "Marmota himalayana", "Mastomys kollmannspergeri", "Mastomys shortridgei", "Mazama chunyi", "Meriones arimalius", "Meriones chengi", "Meriones vinogradovi", "Meriones zarudnyi", "Mesechinus dauuricus", "Mesomys occultus", "Micaelamys granti", "Mico acariensis", "Mico argentatus", "Mico chrysoleucos", "Mico emiliae", "Mico humeralifer", "Mico intermedius", "Mico leucippe", "Mico mauesi", "Mico melanurus", "Mico munduruku", "Mico nigriceps", "Mico rondoni", "Mico saterei", "Mico schneideri", "Microcavia niata", "Microcavia shiptoni", "Microsciurus flaviventer", "Microsciurus santanderensis", "Microtus bucharensis", "Microtus dogramacii", "Microtus guatemalensis", "Microtus qazvinensis", "Millardia kathleenae", "Molossus currentium", "Monodelphis osgoodi", "Monodelphis peruviana", "Mops congicus", "Mops demonstrator", "Mops niveiventer", "Moschus chrysogaster", "Moschus cupreus", "Moschus fuscus", "Moschus leucogaster", "Murina tubinaris", "Mus bufo")

sel_sps_indices <- c(3050, 3061, 3070, 3072, 3100, 3120, 3121, 3122, 3131, 3172,
                     3179, 3184, 3186, 3240, 3262, 3286, 3313, 3314, 3317, 3321,
                     3327, 3328, 3332, 3335, 3336, 3343, 3354, 3363, 3365, 3366,
                     3367, 3377, 3378, 3379, 3395, 3397, 3398, 3400, 3406, 3414,
                     3415, 3430, 3432, 3436, 3449, 3478, 3479, 3509, 3565, 3567)

sel_species <- c("Mus callewaerti", "Mus indutus", "Mus oubanguii", "Mus phillipsi", "Mustela nigripes", "Myomimus setzeri", "Myomyscus angolensis", "Myomyscus brockmani", "Myopterus daubentonii", "Myotis badius", "Myotis bucharensis", "Myotis ciliolabrum", "Myotis csorbai", "Myotis occultus", "Myotis sicarius", "Naemorhedus baileyi", "Neacomys minutus", "Neacomys musseri", "Neacomys spinosus", "Necromys amoenus", "Necromys punctulatus", "Necromys temchuki", "Nelsonia neotomodon", "Neodon fuscus", "Neodon irene", "Neomicroxus latebricola", "Neoromicia matroka", "Neotamias alpinus", "Neotamias bulleri", "Neotamias canipes", "Neotamias cinereicollis", "Neotamias quadrivittatus", "Neotamias ruficaudus", "Neotamias rufus", "Neotoma goldmani", "Neotoma lepida", "Neotoma leucodon", "Neotoma magister", "Neotoma stephensi", "Nephelomys keaysi", "Nephelomys levipes", "Neusticomys ferreirai", "Ningaui ridei", "Niviventer brahma", "Niviventer niviventer", "Notomys cervinus", "Notomys fuscus", "Nycteris woodi", "Ochotona curzoniae", "Ochotona erythrotis")

sel_sps_indices <- c(3568, 3569, 3574, 3575, 3576, 3577, 3578, 3581, 3582, 3584,
                     3587, 3588, 3594, 3595, 3601, 3605, 3606, 3608, 3611, 3614,
                     3618, 3620, 3621, 3623, 3628, 3631, 3632, 3634, 3635, 3639,
                     3643, 3678, 3684, 3688, 3694, 3697, 3716, 3721, 3724, 3727,
                     3729, 3753, 3805, 3844, 3850, 3851, 3876, 3888, 3893, 3921)

sel_species <- c("Ochotona forresti", "Ochotona gloveri", "Ochotona ladacensis", "Ochotona macrotis", "Ochotona mantchurica", "Ochotona nubrica", "Ochotona opaca", "Ochotona pusilla", "Ochotona roylei", "Ochotona rutila", "Ochotona thomasi", "Ochotona turuchanensis", "Octodontomys gliroides", "Octomys mimax", "Oecomys cleberi", "Oecomys paricola", "Oecomys phaeotis", "Oecomys phaeotis", "Oecomys superans", "Oenomys ornatus", "Oligoryzomys andinus", "Oligoryzomys brendae", "Oligoryzomys chacoensis", "Oligoryzomys destructor", "Oligoryzomys griseolus", "Oligoryzomys microtis", "Oligoryzomys moojeni", "Oligoryzomys rupestris", "Oligoryzomys stramineus", "Onychomys arenicola", "Oreoryzomys balneator", "Otomys anchietae", "Otomys cuanzensis", "Otomys helleri", "Otomys sloggetti", "Otomys typus", "Oxymycterus amazonicus", "Oxymycterus hiska", "Oxymycterus inca", "Oxymycterus paramensis", "Oxymycterus roberti", "Pantholops hodgsonii", "Paraxerus lucifer", "Perognathus fasciatus", "Perognathus parvus", "Peromyscus attwateri", "Peromyscus levipes", "Peromyscus nasutus", "Peromyscus polius", "Petaurista mechukaensis")


sel_sps_indices <- c(3922, 3926, 3962, 3963, 3968, 3997, 4000, 4039, 4042, 4048,
                     4050, 4051, 4059, 4060, 4061, 4062, 4066, 4068, 4085, 4107,
                     4108, 4109, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117,
                     4118, 4119, 4121, 4122, 4124, 4128, 4132, 4142, 4143, 4145,
                     4152, 4157, 4159, 4164, 4168, 4169, 4171, 4172, 4173, 4175)

sel_species <- c("Petaurista mishmiensis", "Petaurista xanthotis", "Petromyscus monticularis", "Petromyscus shortridgei", "Phaiomys leucurus", "Philander andersoni", "Philander mcilhennyi", "Phyllotis andium", "Phyllotis caprinus", "Phyllotis magister", "Phyllotis osilae", "Phyllotis wolffsohni", "Piliocolobus langi", "Piliocolobus lulindicus", "Piliocolobus oustaleti", "Piliocolobus parmentieri", "Piliocolobus semlikiensis", "Piliocolobus tholloni", "Pipistrellus inexspectatus", "Pithecia aequatorialis", "Pithecia albicans", "Pithecia cazuzai", "Pithecia chrysocephala", "Pithecia hirsuta", "Pithecia inusta", "Pithecia irrorata", "Pithecia isabela", "Pithecia milleri", "Pithecia mittermeieri", "Pithecia monachus", "Pithecia napensis", "Pithecia pissinattii", "Pithecia rylandsi", "Pithecia vanzolinii", "Planigale gilesi", "Planigale tenuirostris", "Platyrrhinus albericoi", "Platyrrhinus infuscus", "Platyrrhinus ismaeli", "Platyrrhinus masu", "Plecotus ariel", "Plecotus homochrous", "Plecotus kozlovi", "Plecotus strelkovi", "Plecotus wardi", "Plecturocebus aureipalatii", "Plecturocebus bernhardi", "Plecturocebus brunneus", "Plecturocebus caligatus", "Plecturocebus cinerascens")

sel_sps_indices <- c(4176, 4177, 4178, 4179, 4180, 4181, 4182, 4184, 4187, 4188,
                     4189, 4190, 4191, 4192, 4193, 4194, 4211, 4231, 4266, 4272,
                     4277, 4278, 4282, 4283, 4284, 4286, 4287, 4288, 4291, 4292,
                     4364, 4384, 4450, 4470, 4528, 4569, 4577, 4580, 4594, 4603,
                     4606, 4616, 4621, 4635, 4656, 4659, 4561, 4683, 4695, 4697)

sel_species <- c("Plecturocebus cupreus", "Plecturocebus discolor", "Plecturocebus donacophilus", "Plecturocebus dubius", "Plecturocebus grovesi", "Plecturocebus hoffmannsi", "Plecturocebus miltoni", "Plecturocebus moloch", "Plecturocebus ornatus", "Plecturocebus pallescens", "Plecturocebus parecis", "Plecturocebus stephennashi", "Plecturocebus toppini", "Plecturocebus urubambensis", "Plecturocebus vieirai", "Plerotes anchietae", "Poliocitellus franklinii", "Praomys misonnei", "Procapra picticaudata", "Proechimys brevicauda", "Proechimys echinothrix", "Proechimys gardneri", "Proechimys hoplomyoides", "Proechimys kulinae", "Proechimys longicaudatus", "Proechimys oconnelli", "Proechimys pattoni", "Proechimys quadruplicatus", "Proechimys simonsi", "Proechimys steerei", "Pseudomys australis", "Pseudoryzomys simplex", "Pteropus rufus", "Punomys lemminus", "Rattus pyctoris", "Reithrodontomys montanus", "Reithrodontomys zacatecae", "Rhabdomys intermedius", "Rhinolophus bocharicus", "Rhinolophus cohenae", "Rhinolophus damarensis", "Rhinolophus guineensis", "Rhinolophus huananus", "Rhinolophus mcintyrei", "Rhinolophus schnitzleri", "Rhinolophus shortridgei", "Rhinolophus silvestris", "Rhinopithecus bieti", "Rhipidomys austrinus", "Rhipidomys caucensis")

sel_sps_indices <- c(4701, 4702, 4706, 4709, 4711, 4712, 4715, 4722, 4775, 4776, 
                     4778, 4779, 4780, 4782, 4787, 4788, 4792, 4795, 4797, 4803, 
                     4813, 4819, 4821, 4833, 4840, 4844, 4847, 4849, 4850, 4851, 
                     4853, 4854, 4859, 4860, 4861, 4877, 4895, 4901, 4905, 4919,
                     4933, 4951, 4969, 4970, 4974, 4977, 4980, 4988, 5001, 5002)

sel_species <- c("Rhipidomys gardneri", "Rhipidomys ipukensis", "Rhipidomys macconnelli", "Rhipidomys modicus", "Rhipidomys ochrogaster", "Rhipidomys tribei", "Rhipidomys wetzeli", "Rhogeessa hussoni", "Saguinus imperator", "Saguinus inustus", "Saguinus leucopus", "Saguinus martinsi", "Saguinus melanoleucus", "Saguinus mystax", "Saimiri boliviensis", "Saimiri cassiquiarensis", "Saimiri ustus", "Salinomys delicatus", "Salpingotus crassicauda", "Sapajus cay", "Scapanulus oweni", "Scaptonyx fusicaudus", "Scarturus euphratica", "Sciurus alleni", "Sciurus flammifer", "Sciurus ignitus", "Sciurus nayaritensis", "Sciurus oculatus", "Sciurus pucheranii", "Sciurus pyrrhinus", "Sciurus sanborni", "Sciurus spadiceus", "Scleronycteris ega", "Scolomys melanops", "Scolomys ucayalensis", "Scotophilus ejetai", "Scutisorex somereni", "Semnopithecus hector", "Semnopithecus schistaceus", "Sicista pseudonapaea", "Sigmodon ochrognathus", "Sminthopsis ooldea", "Sorex arizonae", "Sorex asper", "Sorex buchariensis", "Sorex cansulus", "Sorex cylindricauda", "Sorex haydeni", "Sorex merriami", "Sorex milleri")

sel_sps_indices <- c(5007, 5008, 5016, 5029, 5030, 5040, 5050, 5058, 5067, 5068, 
                     5084, 5089, 5095, 5112, 5116, 5124, 5142, 5193, 5196, 5199, 
                     5211, 5292, 5295, 5296, 5298, 5304, 5307, 5308, 5319, 5320, 
                     5321, 5322, 5323, 5324, 5325, 5326, 5327, 5328, 5329, 5330, 
                     5331, 5333, 5335, 5337, 5339, 5340, 5343, 5344, 5345, 5347)

sel_species <- c("Sorex nanus", "Sorex neomexicanus", "Sorex preblei", "Sorex tenellus", "Sorex thibetanus", "Sorex yukonicus", "Spalax graecus", "Spermophilus brevicauda", "Spermophilus ralli", "Spermophilus relictus", "Steatomys bocagei", "Steatomys opimus", "Stenocephalemys ruppi", "Sturnira magna", "Sturnira oporaphilum", "Stylodipus andrewsi", "Suncus remyi", "Sylvilagus nuttallii", "Sylvilagus robustus", "Sylvilagus varynaensis", "Sylvisorex oriundus", "Taterillus congicus", "Taterillus lacustris", "Taterillus petteri", "Taterillus tranieri", "Thallomys loringi", "Thallomys shortridgei", "Thalpomys cerradensis", "Thomasomys australis", "Thomasomys baeops", "Thomasomys bombycinus", "Thomasomys caudivarius", "Thomasomys cinereiventer", "Thomasomys cinereus", "Thomasomys cinnameus", "Thomasomys contradictus", "Thomasomys daphne", "Thomasomys dispar", "Thomasomys eleusis", "Thomasomys emeritus", "Thomasomys erro", "Thomasomys gracilis", "Thomasomys hylophilus", "Thomasomys ischyrus", "Thomasomys ladewi", "Thomasomys laniger", "Thomasomys nicefori", "Thomasomys niveipes", "Thomasomys notatus", "Thomasomys oreas")

#############


sel_sps_indices <- c(5348, 5349, 5350, 5351, 5352, 5354, 5355, 5356, 5361, 5371,
                     5373, 5375, 5376, 5379, 5380, 5382, 5383, 5385, 5386, 5403, 
                     5404, 5410, 5431, 5499, 5502, 5504, 5505, 5508, 5511, 5514, 
                     5518, 5536, 5585, 5600, 5616, 5617, 5620)

sel_species <- c("Thomasomys paramorum", "Thomasomys popayanus", "Thomasomys praetor", "Thomasomys princeps", "Thomasomys pyrrhonotus", "Thomasomys silvestris", "Thomasomys taczanowskii", "Thomasomys ucucha", "Thomomys clusius", "Thrichomys inermis", "Thrichomys pachyurus", "Thylamys cinderella", "Thylamys citellus", "Thylamys karimii", "Thylamys macrurus", "Thylamys pulchellus", "Thylamys pusillus", "Thylamys velutinus", "Thylamys venustus", "Tolypeutes matacus", "Tolypeutes tricinctus", "Toromys rhipidurus", "Trachypithecus shortridgei", "Tympanoctomys barrerae", "Typhlomys cinereus", "Urocitellus armatus", "Urocitellus beldingi", "Urocitellus columbianus", "Urocitellus mollis", "Urocitellus richardsonii", "Urocricetus alticola", "Uropsilus gracilis", "Vulpes ferrilata", "Xeronycteris vieirai", "Zelotomys hildegardeae", "Zelotomys woosnami", "Zygodontomys brunneus")

#save list of selected species
setwd(wd_lists)
saveRDS(sel_species, 'Selected_mammal_species_12')
