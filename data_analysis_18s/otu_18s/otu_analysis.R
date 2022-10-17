# Analysis of eukaryotic biodiversity data using the phyloseq package in R

# OTU data analysis with the phyloseq package for statistical analysis
## According to Vaulot (2021) and McMurdie & Holmes (2013), this statistical analysis was performed.

## To install necessary packages
install.packages("pacman")
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

install.packages("reshape2")
install.packages("ape")
install.packages("gridExtra")
install.packages("plotly")
install.packages("vegan")
install.packages("dendextend")
install.packages("tidyr")
install.packages("rms")
install.packages("effectsize")
install.packages("lme4")
install.packages("picante")
install.packages("cowplot")
install.packages("here")
install.packages("lsr")
install.packages("MASS")
install.packages("rcompanion")
install.packages("ggiraphExtra")
install.packages("optimx")
install.packages("mvabund")
install.packages("gllvm")
install.packages("radiant.data ")
install.packages("tidyverse")
install.packages("MuMIn")
install.packages("corrplot")
install.packages("gclus")
install.packages("coefplot")
install.packages("jtools")
install.packages("ggstance")
install.packages("interactions")
install.packages("blmeco")
install.packages("tseries")
install.packages("plyr")
install.packages("lmerTest")
install.packages("DHARMa")
install.packages("car")

## To load the installed packages with the pacman package
pacman::p_load(pacman, dplyr, readxl, ggplot2, phyloseq, microbiome, reshape2, ape, gridExtra, plotly, vegan, dendextend, tidyr, rms, effectsize, lme4, picante, cowplot, here, lsr, MASS, rcompanion, ggiraphExtra, optimx, mvabund, gllvm, radiant.data, tidyverse, MuMIn, corrplot, gclus, coefplot, jtools, ggstance, interactions, blmeco, tseries, plyr, lmerTest, DHARMa, car)
here()

## To read the data and create phyloseq objects
otu_mat_1 <- read_excel(here("otu_18s","eucaryotic_week2_2.xlsx"), sheet = "OTU matrix")
otu_mat_2 <- read_excel(here("otu_18s","eucaryotic_week3_3.xlsx"), sheet = "OTU matrix")

tax_mat_1 <- read_excel(here("otu_18s","eucaryotic_week2_2.xlsx"), sheet = "Taxonomy table")
tax_mat_2 <- read_excel(here("otu_18s","eucaryotic_week3_3.xlsx"), sheet = "Taxonomy table")

samples_df_1 <- read_excel(here("otu_18s","eucaryotic_week2_2.xlsx"), sheet = "Samples")
samples_df_2 <- read_excel(here("otu_18s","eucaryotic_week3_3.xlsx"), sheet = "Samples")

## Phyloseq objects must have row.names
## To define the row names from the otu column
otu_mat_1 <- otu_mat_1 %>%
  tibble::column_to_rownames("otu")

otu_mat_2 <- otu_mat_2 %>%
  tibble::column_to_rownames("otu") 

## To idem for the two other matrices
## For taxa data
tax_mat_1 <- tax_mat_1 %>% 
  tibble::column_to_rownames("otu")

tax_mat_2 <- tax_mat_2 %>% 
  tibble::column_to_rownames("otu")

## For sample data
samples_df_1 <- samples_df_1 %>% 
  tibble::column_to_rownames("sample")

samples_df_2 <- samples_df_2 %>% 
  tibble::column_to_rownames("sample")

## Transformation to matrices of otu and tax tables (the sample table can be left as a data frame)
otu_mat_1 <- as.matrix(otu_mat_1)
otu_mat_2 <- as.matrix(otu_mat_2)

tax_mat_1 <- as.matrix(tax_mat_1)
tax_mat_2 <- as.matrix(tax_mat_2)

## Transforming to phyloseq objects
OTU_1 = otu_table(otu_mat_1, taxa_are_rows = TRUE)
TAX_1 = tax_table(tax_mat_1)
samples_1 = sample_data(samples_df_1)

OTU_2 = otu_table(otu_mat_2, taxa_are_rows = TRUE)
TAX_2 = tax_table(tax_mat_2)
samples_2 = sample_data(samples_df_2)

phyloseq_obj_1 <- phyloseq(OTU_1, TAX_1, samples_1)
phyloseq_obj_1

phyloseq_obj_2 <- phyloseq(OTU_2, TAX_2, samples_2)
phyloseq_obj_2

## To merge the data of two phyloseq objects from week 2 and week 3
phyloseq_object <- merge_phyloseq(phyloseq_obj_1, phyloseq_obj_2)
phyloseq_object

## To create a phy_tree slot from the original phyloseq object
random_tree = rtree(ntaxa(phyloseq_object), rooted=TRUE, tip.label=taxa_names(phyloseq_object))
plot(random_tree)

## To make a new phyloseq object with the newly created phy_tree
phyloseq_object_with_tree <- merge_phyloseq(phyloseq_object, random_tree)
phyloseq_object_with_tree



# The preprocessing of the data

## To transform the data from phyloseq_object_with_tree to relative abundance (to transform the count of samples to percent)
phyloseq_object_tsc <- transform_sample_counts(phyloseq_object_with_tree, function(x) x / sum(x))
phyloseq_object_tsc

phy_obj_tsc_otu <- as.data.frame(otu_table(phyloseq_object_tsc))

## The coverage index (The index calculates the number of groups required to have a certain proportion of the ecosystem occupied)
## This index was calculated according to Lahti & Shetty (Bioconductor, 2017) 
phy_coverage <- coverage(phyloseq_object_tsc, threshold=1.0)
phy_coverage

phy_coverage_df <- as.data.frame(phy_coverage)
phy_coverage_df
phy_coverage_df_1 <- tibble::rownames_to_column(phy_coverage_df, "samples")
phy_coverage_df_1

theme_set(theme_bw())
phy_coverage_df_1_viz <- ggplot(phy_coverage_df_1, aes(x=samples, y=phy_coverage)) + 
  geom_bar(stat="identity", width=.5, fill="tomato3") +
  ylim(0,2000)+
  labs(title="Bar chart of samples and phy_coverage", 
       subtitle="samples Vs phy_coverage", 
       caption="phy_coverage_df_1") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
phy_coverage_df_1_viz

## To transform the readings below 0.003% to 0% for removing noises on the data
phyloseq_object_tsc_minthreshold  = transform_sample_counts(phyloseq_object_with_tree, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})

phy_obj_tsc_minthres_otu <- as.data.frame(otu_table(phyloseq_object_tsc_minthreshold))

## The coverage index after removing noise below 0.003% from the data
phy_coverage_after_t <- coverage(phyloseq_object_tsc_minthreshold, threshold=1.0)
phy_coverage_after_t

phy_coverage_after_t_df <- as.data.frame(phy_coverage_after_t)
phy_coverage_after_t_df
phy_coverage_after_t_df_1 <- tibble::rownames_to_column(phy_coverage_after_t_df, "samples")
phy_coverage_after_t_df_1

theme_set(theme_bw())
phy_coverage_after_t_df_1_viz <- ggplot(phy_coverage_after_t_df_1, aes(x=samples, y=phy_coverage_after_t)) + 
  geom_bar(stat="identity", width=.5, fill="tomato3") +
  ylim(0,80)+
  labs(title="Bar chart of samples and phy_coverage_after_t", 
       subtitle="samples Vs phy_coverage_after_t", 
       caption="phy_coverage_after_t_df_1") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
phy_coverage_after_t_df_1_viz

phy_filtered <- phyloseq_object_tsc_minthreshold

## To visualize the data
sample_names(phy_filtered)
rank_names(phy_filtered)
sample_variables(phy_filtered)

## Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(phy_filtered), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered)){OTU1 <- t(OTU1)}
## Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf

## Extract taxa matrix from the phyloseq object
TAX1 = as(tax_table(phy_filtered), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered)){TAX1 <- t(TAX1)}
## Coerce to data.frame
TAXdf = as.data.frame(TAX1)
TAXdf

## Extract sample matrix from the phyloseq object
samples1 = as(sample_data(phy_filtered), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered)){samples1 <- t(samples1)}
## Coerce to data.frame
samplesdf = as.data.frame(samples1)
samplesdf

## To keep only taxa of interests (for Algae) according to phylum, genus, family and order level
sub_algae_phylum <- subset_taxa(phy_filtered, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Streptophyta", "Xanthophyceae", "Nucleariidae_and_Fonticula_group", "Phragmoplastophyta"))
sub_algae_genus <- subset_taxa(phy_filtered, genus %in% c("Messastrum", "Melosira", "Poteriospumella", "Pseudellipsoidion", "Symbiodinium", "Nasturtium", "Phalansterium", "Reynoutria", "Pleurochloris", "Mallomonas", "Makinoella", "Nuclearia", "Achnanthidium", "Acutodesmus", "Akashiwo", "Allapsa", "Ancyromonas", "Ankistrodesmus", "Aphanochaete", "Aphelidium", "Bicosoeca", "Biecheleria", "Botrydiopsis", "Botryococcus", "Bracteacoccus", "Carteria", "Chamaetrichon", "Characiopodium", "Characium", "Chlamydocapsa", "Chlamydomonas", "Chlorella", "Chlorochytrium", "Chlorococcum", "Chloromonas", "Chlorosarcinopsis", "Choricystis", "Chroomonas", "Chrysamoeba", "Chrysolepidomonas", "Closterium", "Coccomyxa", "Cocconeis", "Colponema", "Cosmarium", "Ctenocladus", "Cyclotella", "Cylindrocystis", "Cymbella", "Deasonia", "Desmodesmus", "Dictyochloropsis", "Dictyosphaerium", "Diplosphaera", "Dunaliella", "Dysmorphococcus", "Elliptochloris", "Elongatocystis", "Encyonema", " Encyonopsis", "Eucocconeis", "Eustigmatos", "Flechtneria", "Fragilaria", "Geissleria", "Geminigera", "Gloeotilopsis", "Gomphonema", "Goniomonas", "Hafniomonas", "Hemiselmis", "Heterochlorella", "Heveochlorella", "Hormotilopsis", "Kalinella", "Klebsormidium", "Kraftionema", "Kumanoa", "Lessonia", "Lobochlamys", "Lobomonas", "Luticola", "Mallomonas", "Marsupiomonas", "Microglena", "Microthamnion", "Monomastix", "Monoraphidium", "Mougeotia", "Mychonastes", "Mysteriochloris", "Navicula", "Nemalionopsis", "Neochloris", "Neochlorosarcina", "Neocystis", "Neodesmus", "Neospongiococcum", "Nephroselmis", "Nitzschia", "Ochromonas", "Oedogonium", "Olpidiopsis", "Oltmannsiellopsis", "Oogamochlamys", "Parietochloris", "Paulschulzia", "Pectinodesmus", "Pedinomonas", "Phyllosiphon", "Placoneis", "Plagioselmis", "Planktosphaeria", "Planothidium", "Pleurochrysis", "Polytoma", "Printzina", "Protosiphon", "Psammodictyon", "Pseudellipsoidion", "Pseudochlorella", "Pseudomuriella", "Pseudopediastrum", "Pyrobotrys", "Reimeria", "Rhodomonas", "Scenedesmus", "Scotinosphaera", "Sellaphora", "Sheathia", "Spirogyra", "Spongiochloris", "Spumella", "Stauridium", "Stephanodiscus", "Stichococcus", "Stigeoclonium", "Surirella", "Symbiochloris", "Symbiodinium", "Synura", "Teleaulax", "Tetracystis", "Tetranephris", "Thorea", "Trebouxia", "Trentepohlia", "Triposolenia", "Ulnaria", "Ulothrix", "Ulvella", "Uroglena", "Vitreochlamys", "Wislouchiella", "Xylochloris"))
sub_algae_family <- subset_taxa(phy_filtered, family %in% c("Achnanthidiaceae", "Actinochloridaceae", "Amphisoleniaceae", "Aphanochaetaceae", "Bacillariaceae", "Batrachospermaceae", "Biecheleriaceae", "Botryochloridaceae", "Botryococcaceae", "Bracteacoccaceae", "Catenulaceae", "Chaetopeltidaceae", "Chaetophoraceae", "Chlamydomonadaceae", "Chlorellaceae", "Chlorococcaceae", "Chromulinaceae", "Chroomonadaceae", "Chrysolepidomonadaceae", "Closteriaceae", "Coccomyxaceae", "Cocconeidaceae", "Cryptomonadaceae", "Cymbellaceae", "Desmidiaceae", "Diadesmidaceae", "Dunaliellaceae", "Eustigmataceae", "Fragilariaceae", "Geminigeraceae", "Gomphonemataceae", "Goniomonadaceae", "Gymnodiniaceae", "Hemiselmidaceae", "Heterocapsaceae", "Hydrodictyaceae", "Klebsormidiaceae", "Kraftionemaceae", "Kryptoperidiniaceae", "Lophodiniaceae", "Marsupiomonadaceae", "Melosiraceae", "Mesotaeniaceae", "Monomastigaceae", "Mrazekiidae", "Mychonastaceae", "Naviculaceae", "Neochloridaceae", "Oikomonadaceae", "Oocystaceae", "Pedinomonadaceae", "Phyllosiphonaceae", "Phytodiniaceae", "Pleurochloridaceae", "Pleurochrysidaceae", "Pseudomuriellaceae", "Pyrenomonadaceae", "Scenedesmaceae", "Schizochlamydaceae", "Scotinosphaeraceae", "Selenastraceae", "Sellaphoraceae", "Sphaerocystidae", "Spondylomoraceae", "Stephanodiscaceae", "Suessiaceae", "Surirellaceae", "Thaumatomonadida", "Thoreaceae", "Trentepohliaceae", "Ulnariaceae", "Ulotrichaceae", "Ulvellaceae", "Uronemataceae", "Volvocaceae", "Zygnemataceae"))
sub_algae_order <- subset_taxa(phy_filtered, order %in% c("Bacillariales", "Batrachospermales", "Cercomonadida", "Chaetopeltidales", "Chaetophorales", "Chlamydomonadales", "Chlorellales", "Chlorodendrales", "Chlorosarcinales", "Chromulinales", "Coccolithales", "Cocconeidales", "Cryptomonadales", "Cymbellales", "Desmidiales", "Dinophysiales", "Eustigmatales", "Fragilariales", "Gymnodiniales","Gymnodiniphycidae", "Hibberdiales", "Klebsormidiales", "Laminariales","Licmophorales", "Lophodiniales", "Marsupiomonadales", "Mastogloiales","Melosirales", "Microthamniales","Mischococcales", "Monomastigales","Naviculales", "Ochromonadales", "Oedogoniales", "Pavlovales","Pedinomonadales", "Peridiniales","Phytodiniales", "Pyrenomonadales","Scotinosphaerales", "Sphaeropleales","Suessiales", "Surirellales","Synurales", "Tetrasporales","Thalassiophysales", "Thalassiosirales", "Thaumatomonadida", "Thoreales","Trentepohliales", "Ulotrichales","Ulvales", "Zygnematales"))

sub_algae_phylum_cov <- coverage(sub_algae_phylum, threshold=1.0)
sub_algae_phylum_cov
sub_algae_genus_cov <- coverage(sub_algae_genus, threshold=1.0)
sub_algae_genus_cov
sub_algae_family_cov <- coverage(sub_algae_family, threshold=1.0)
sub_algae_family_cov
sub_algae_order_cov <- coverage(sub_algae_order, threshold=1.0)
sub_algae_order_cov

## To keep only taxa of interests (for Protists) according to phylum, genus, family and order level
sub_protist_phylum <- subset_taxa(phy_filtered, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Tubulinea", "Bicosoecida", "Discosea", "Labyrinthulomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida", "Colponemidia", "MAST-12", "Rigifilida", "Schizoplasmodiida"))
sub_protist_genus <- subset_taxa(phy_filtered, genus %in% c("Gymnophrys", "Ichthyophthirius", "Cephalodella", "Flabellula", "Bicosoeca", "Lepidodermella", "Sorodiplophrys", "Phytophthora", "Protosporangium", "Cryptosporidium", "Colpoda", "Acanthamoeba", "Acanthoeca", "Acineria", "Acorus", "Adelina", "Amoebogregarina", "Amphileptus", "Anteholosticha", "Anurofeca", "Apicystis", "Arachnula", "Arboramoeba", "Arcella", "Arcuospathidium", "Ascogregarina", "Aspidisca", "Aurigamonas", "Babesia", "Bodomorpha", "BOLA868", "Bromeliothrix", "Bryometopus", "Cavernomonas", "Ceratiomyxella", "Cercomonas", "Chaenea", "Chilodonella", "Codosiga", "Coleps", "Colpidium", "Colpoda", "Copromyxa", "Cryptosporidium", "Cyrtohymena", "Cyrtolophosis", "Desmarella", "Dictyostelium", "Difflugia", "Dileptus", "Diplophrys", "Eimeria", "Eocercomonas", "Epalxella", "Etoschophrya", "Euglypha", "Filamoeba", "Flabellula", "Flamella", "Freshwater_Choanoflagellates_2", "Frontonia", "Fuscheria", "Geneiorhynchus", "Gephyramoeba", "Gonostomum", "Gregarina", "Grellamoeba", "Hartmannella", "Hatena", "Hausmanniella", "Hemolivia", "Heteromita", "Heterophrys", "Holosticha", "Ischnamoeba", "Kahliella", "Kraken", "Leidyana", "Lembadion", "Leptomyxa", "Leptopharynx", "Leptophrys", "Limnostrombidium", "Litonotus", "Loxophyllum", "Malawimonas", "Mattesia", "Metabolomonas", "Microdiaphanosoma", "Micronuclearia", "Monas", "Mykophagophrys", "Nudifila", "Obertrumia", "Oikomonas", "Oxytricha", "Palustrimonas", "Paracercomonas", "Paraflabellula", "Parafurgasonia", "Paralagenidium", "Paramecium", "Paramoebidium", "Paraphysomonas", "Paraschneideria", "Paraurostyla", "Perisincirra", "Phalansterium", "Phascolodon", "Phialina", "Placus", "Platyophrya", "Poterioochromonas", "Poteriospumella", "Proleptomonas", "Protocyclidium", "Pseudocyrtolophosis", "Pseudodifflugia", "Pseudogastrostyla", "Pseudoplatyophrya", "Pseudouroleptus", "Pseudourostyla", "Psychodiella", "Ptolemeba", "Rhizamoeba", "Rhogostoma", "Saccamoeba", "Salpingoeca", "Spathidium", "Sphaeroeca", "Spongomonas", "Stylocephalus", "Stylonychia", "Syncystis", "Tetrahymena", "Thaumatomastix", "Thaumatomonas", "Theileria", "Tintinnidium", "Trachelius", "Trinema", "Trochilia", "Uroleptus", "Uronema", "Vampyrella", "Vermamoeba", "Viridiraptor", "Xiphocephalus"))
sub_protist_family <- subset_taxa(phy_filtered, family %in% c("Acanthamoebidae", "Acanthoecidae", "Acropisthiidae", "Actinocephalidae", "Adeleidae", "Allapsidae", "Amoebidiaceae", "Amphitraemidae", "Ancyromonadidae", "Apusomonadidae", "Arcellidae", "Babesiidae", "Bromeliophryidae", "Bryometopidae", "Cercomonadidae", "Chilodonellidae", "Codonellidae", "Colepidae", "Colpodea", "Colpodellidae", "Colpodidae", "Condylostomatidae", "Cryptosporidiidae", "Cyrtolophosididae", "Deltopylidae", "Dictyosteliaceae", "Difflugiidae", "Dileptidae", "Dysteriidae", "Echinamoebidae", "Eimeriidae", "Epalxellidae", "Euglyphida", "Euglyphidae", "Euplotia", "Flabellulidae", "Frontoniidae", "Gregarinidae", "Grossglockneriidae", "Haptoria", "Hartmannellidae", "Hausmanniellidae", "Heterophryidae", "Hibberdiaceae", "Holostichidae", "Hypotrichia", "Karyolysidae", "Kreyellidae", "Lacrymariidae", "Lecudinidae", "Leidyanidae", "Lembadionidae", "Leptomyxidae", "Lipotrophidae", "Litonotidae", "Mallomonadaceae", "Marynidae", "Microthoracidae", "Monocystidae", "Nassulidae", "Nucleariidae", "Ochromonadaceae", "Oligohymenophorea", "Olpidiopsidaceae", "Ophryocystidae", "Oxytrichidae", "Parameciidae", "Philasteridae", "Phyllopharyngea", "Placidae", "Plagiopylea", "Plasmodiophoridae", "Platyophryidae", "Prostomatea", "Pseudocharaciopsidaceae", "Pseudodifflugiidae", "Pseudourostylidae", "Rhizaspididae", "Rhogostomidae", "Salpingoecidae", "Sandonidae", "Spathidiidae", "Spirofilidae", "Spongomonadidae", "Stenophoridae", "Strombidiidae", "Stylocephalidae", "Syncystidae", "Tetrahymenidae", "Thaumatomastigidae", "Theileriidae", "Tintinnidae", "Tracheliidae", "Trachelophyllidae", "Trinematidae", "Turaniellidae", "Vampyrellidae", "Viridiraptoridae", "Woodruffiidae"))
sub_protist_order <- subset_taxa(phy_filtered, order %in% c("Arcellinida",  "Bicosoecida", "Bryometopida", "Chlamydodontida", "Choanoflagellida", "Choreotrichida", "Colpodida", "Conthreep", "Cryomonadida", "Cryptosporida", "Cyrtolophosidida", "Dysteriida", "Echinamoebida", "Eucoccidiorida", "Euglyphida", "Eugregarinorida", "Glissomonadida", "Grossglockneriida", "Haptorida","Heterotrichida", "Hymenostomatida","Leptomyxida", "Litostomatea","Longamoebia", "Microthoracida","Nassulida", "Neogregarinorida","Odontostomatida", "Peniculida","Philasterida", "Pleurostomatida", "Prorodontida", "Protosteliales","Salpingoecidae", "Silicofilosea","Spirotrichea", "Spongomonadida","Sporadotrichida", "Stichotrichida","Tintinnida", "Urostylida","Vampyrellida"))

sub_protist_phylum_cov <- coverage(sub_protist_phylum, threshold=1.0)
sub_protist_phylum_cov
sub_protist_genus_cov <- coverage(sub_protist_genus, threshold=1.0)
sub_protist_genus_cov
sub_protist_family_cov <- coverage(sub_protist_family, threshold=1.0)
sub_protist_family_cov
sub_protist_order_cov <- coverage(sub_protist_order, threshold=1.0)
sub_protist_order_cov

## Bar graphs of relative abundance

## Basic bar graphs based on phylum and order level for Algae
bar_algae_phylum <- microbiome::transform(sub_algae_phylum, "compositional")
plot_bar(bar_algae_phylum, fill="phylum", title = "The relative abundances of the algae at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_algae_order <- microbiome::transform(sub_algae_order, "compositional")
plot_bar(bar_algae_order, fill="order", title = "The relative abundances of the algae at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")

## Basic bar graphs based on phylum and order level for Protists
bar_protist_phylum <- microbiome::transform(sub_protist_phylum, "compositional")
plot_bar(bar_protist_phylum, fill="phylum", title = "The relative abundances of protists at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_protist_order <- microbiome::transform(sub_protist_order, "compositional")
plot_bar(bar_protist_order, fill="order", title = "The relative abundances of protists at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")

# According to Mariadassou et al., (2016); Gauthier & Derome, (2021), this statistical analysis was performed.
## Exploring the biodiversity : Alpha-diversity
## According to https://joey711.github.io/phyloseq/plot_richness-examples.html,
## For alpha diversity, it is advisable to prune OTUs that are not present in any of the samples.
## But that's all we can trim.
## It is tempting to trim noise, but there are many estimates of richness based on the singletons and double tons of the abundance data.
## To get a meaningful estimate, we need to leave it in the data set.

## For the algae phylum level
alph_alg_sset_phy <- subset_taxa(phyloseq_object_with_tree, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Streptophyta", "Xanthophyceae", "Nucleariidae_and_Fonticula_group", "Phragmoplastophyta"))
alph_alg_sset_phy_1 <- prune_taxa(taxa_sums(alph_alg_sset_phy) > 0, alph_alg_sset_phy)

alpha.diversity_algae_phylum <- estimate_richness(alph_alg_sset_phy_1, measures = c("Observed", "Shannon"))
alpha.diversity_algae_phylum

alpha_algae_phylum <- plot_richness(alph_alg_sset_phy_1, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
      geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_phylum)

pielou_algae_phylum <- evenness(alph_alg_sset_phy_1, 'pielou')
pielou_algae_phylum

dat_lm_algae_phylum <- data.frame(sample_data(alph_alg_sset_phy_1))

alpha_lm_algae_phylum <- cbind(alpha.diversity_algae_phylum, pielou_algae_phylum, dat_lm_algae_phylum)
alpha_lm_algae_phylum

plot_pielou_algae_phylum <- ggplot(alpha_lm_algae_phylum, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_phylum

## For the algae order level
alph_alg_sset_ord <- subset_taxa(phyloseq_object_with_tree, order %in% c("Bacillariales", "Batrachospermales", "Cercomonadida", "Chaetopeltidales", "Chaetophorales", "Chlamydomonadales", "Chlorellales", "Chlorodendrales", "Chlorosarcinales", "Chromulinales", "Coccolithales", "Cocconeidales", "Cryptomonadales", "Cymbellales", "Desmidiales", "Dinophysiales", "Eustigmatales", "Fragilariales", "Gymnodiniales","Gymnodiniphycidae", "Hibberdiales", "Klebsormidiales", "Laminariales","Licmophorales", "Lophodiniales", "Marsupiomonadales", "Mastogloiales","Melosirales", "Microthamniales","Mischococcales", "Monomastigales","Naviculales", "Ochromonadales", "Oedogoniales", "Pavlovales","Pedinomonadales", "Peridiniales","Phytodiniales", "Pyrenomonadales","Scotinosphaerales", "Sphaeropleales","Suessiales", "Surirellales","Synurales", "Tetrasporales","Thalassiophysales", "Thalassiosirales", "Thaumatomonadida", "Thoreales","Trentepohliales", "Ulotrichales","Ulvales", "Zygnematales"))
alph_alg_sset_ord_1 <- prune_taxa(taxa_sums(alph_alg_sset_ord) > 0, alph_alg_sset_ord)

alpha.diversity_algae_order <- estimate_richness(alph_alg_sset_ord_1, measures = c("Observed", "Shannon"))
alpha.diversity_algae_order

alpha_algae_order <- plot_richness(alph_alg_sset_ord_1, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_order)

pielou_algae_order <- evenness(alph_alg_sset_ord_1, 'pielou')
pielou_algae_order

dat_lm_algae_order <- data.frame(sample_data(alph_alg_sset_ord_1))

alpha_lm_algae_order <- cbind(alpha.diversity_algae_order, pielou_algae_order, dat_lm_algae_order)
alpha_lm_algae_order

plot_pielou_algae_order <- ggplot(alpha_lm_algae_order, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_order

## For the protist phylum level
alph_prot_sset_phy <- subset_taxa(phyloseq_object_with_tree, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Tubulinea", "Bicosoecida", "Discosea", "Labyrinthulomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida", "Colponemidia", "MAST-12", "Rigifilida", "Schizoplasmodiida"))
alph_prot_sset_phy_1 <- prune_taxa(taxa_sums(alph_prot_sset_phy) > 0, alph_prot_sset_phy)

alpha.diversity_protist_phylum <- estimate_richness(alph_prot_sset_phy_1, measures = c("Observed", "Shannon"))
alpha.diversity_protist_phylum

alpha_protist_phylum <- plot_richness(alph_prot_sset_phy_1, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_phylum)

pielou_protist_phylum <- evenness(alph_prot_sset_phy_1, 'pielou')
pielou_protist_phylum

dat_lm_protist_phylum <- data.frame(sample_data(alph_prot_sset_phy_1))

alpha_lm_protist_phylum <- cbind(alpha.diversity_protist_phylum, pielou_protist_phylum, dat_lm_protist_phylum)
alpha_lm_protist_phylum

plot_pielou_protist_phylum <- ggplot(alpha_lm_protist_phylum, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_phylum

## For the protist order level
alph_prot_sset_ord <- subset_taxa(phyloseq_object_with_tree, order %in% c("Arcellinida",  "Bicosoecida", "Bryometopida", "Chlamydodontida", "Choanoflagellida", "Choreotrichida", "Colpodida", "Conthreep", "Cryomonadida", "Cryptosporida", "Cyrtolophosidida", "Dysteriida", "Echinamoebida", "Eucoccidiorida", "Euglyphida", "Eugregarinorida", "Glissomonadida", "Grossglockneriida", "Haptorida","Heterotrichida", "Hymenostomatida","Leptomyxida", "Litostomatea","Longamoebia", "Microthoracida","Nassulida", "Neogregarinorida","Odontostomatida", "Peniculida","Philasterida", "Pleurostomatida", "Prorodontida", "Protosteliales","Salpingoecidae", "Silicofilosea","Spirotrichea", "Spongomonadida","Sporadotrichida", "Stichotrichida","Tintinnida", "Urostylida","Vampyrellida"))
alph_prot_sset_ord_1 <- prune_taxa(taxa_sums(alph_prot_sset_ord) > 0, alph_prot_sset_ord)

alpha.diversity_protist_order <- estimate_richness(alph_prot_sset_ord_1, measures = c("Observed", "Shannon"))
alpha.diversity_protist_order

alpha_protist_order <- plot_richness(alph_prot_sset_ord_1, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_order)

pielou_protist_order <- evenness(alph_prot_sset_ord_1, 'pielou')
pielou_protist_order

dat_lm_protist_order <- data.frame(sample_data(alph_prot_sset_ord_1))

alpha_lm_protist_order <- cbind(alpha.diversity_protist_order, pielou_protist_order, dat_lm_protist_order)
alpha_lm_protist_order

plot_pielou_protist_order <- ggplot(alpha_lm_protist_order, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_order



#################################
## Followed by https://ase.tufts.edu/bugs/guide/assets/mixed_model_guide.html
## To check over dispersion by a function
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
###################


# Mixed Models(MMs)
## According to Walker (2021) and Bolker (2018), this statistical analysis was performed.
## Also followed by https://static1.squarespace.com/static/5eb33c095018927ea433a883/t/5fbb491e773f2a6929bb97ed/1606109472637/Mixed-Models-in-R.pdf

# For the algae phylum level

## For Observed
### To check the normality of the data
hist(alpha_lm_algae_phylum$Observed)
qqnorm(alpha_lm_algae_phylum$Observed)
qqline(alpha_lm_algae_phylum$Observed)
shapiro.test(alpha_lm_algae_phylum$Observed)
### Mixed model
lmm_algae_phylum_block_alpha_Observed <- lmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), data = alpha_lm_algae_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(lmm_algae_phylum_block_alpha_Observed)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(lmm_algae_phylum_block_alpha_Observed)~as.factor(alpha_lm_algae_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(lmm_algae_phylum_block_alpha_Observed, Observed~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_phylum$Observed~fitted(lmm_algae_phylum_block_alpha_Observed))

### Normality of Residuals
qqnorm(resid(lmm_algae_phylum_block_alpha_Observed))
qqline(resid(lmm_algae_phylum_block_alpha_Observed))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(lmm_algae_phylum_block_alpha_Observed)
# Distribution of the random effects
qqnorm (ranef(lmm_algae_phylum_block_alpha_Observed)$header_tanks_block[[1]])

plot(simulateResiduals(lmm_algae_phylum_block_alpha_Observed))
Anova(lmm_algae_phylum_block_alpha_Observed)

coef(summary(lmm_algae_phylum_block_alpha_Observed))
overdisp_fun(lmm_algae_phylum_block_alpha_Observed)
print(summary(lmm_algae_phylum_block_alpha_Observed), correlation=TRUE)
eta_squared(lmm_algae_phylum_block_alpha_Observed)


## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_phylum$Shannon)
qqnorm(alpha_lm_algae_phylum$Shannon)
qqline(alpha_lm_algae_phylum$Shannon)
shapiro.test(alpha_lm_algae_phylum$Shannon)

### Mixed model
glmm_algae_phylum_block_alpha_Shannon <- glmer(Shannon ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_algae_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_algae_phylum_block_alpha_Shannon)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_algae_phylum_block_alpha_Shannon)~as.factor(alpha_lm_algae_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_algae_phylum_block_alpha_Shannon, Shannon~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_phylum$Shannon~fitted(glmm_algae_phylum_block_alpha_Shannon))

### Normality of Residuals
qqnorm(resid(glmm_algae_phylum_block_alpha_Shannon))
qqline(resid(glmm_algae_phylum_block_alpha_Shannon))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_algae_phylum_block_alpha_Shannon)
# Distribution of the random effects
qqnorm (ranef(glmm_algae_phylum_block_alpha_Shannon)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_algae_phylum_block_alpha_Shannon))
Anova(glmm_algae_phylum_block_alpha_Shannon)

coef(summary(glmm_algae_phylum_block_alpha_Shannon))
overdisp_fun(glmm_algae_phylum_block_alpha_Shannon)
print(summary(glmm_algae_phylum_block_alpha_Shannon), correlation=TRUE)
parameters::p_value(glmm_algae_phylum_block_alpha_Shannon)
eta_squared(glmm_algae_phylum_block_alpha_Shannon)

## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_phylum$pielou)
qqnorm(alpha_lm_algae_phylum$pielou)
qqline(alpha_lm_algae_phylum$pielou)
shapiro.test(alpha_lm_algae_phylum$pielou)
### Mixed model
glmm_algae_phylum_block_alpha_pielou <- glmer(pielou ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_algae_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_algae_phylum_block_alpha_pielou)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_algae_phylum_block_alpha_pielou)~as.factor(alpha_lm_algae_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_algae_phylum_block_alpha_pielou, pielou~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_phylum$pielou~fitted(glmm_algae_phylum_block_alpha_pielou))

### Normality of Residuals
qqnorm(resid(glmm_algae_phylum_block_alpha_pielou))
qqline(resid(glmm_algae_phylum_block_alpha_pielou))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_algae_phylum_block_alpha_pielou)
# Distribution of the random effects
qqnorm (ranef(glmm_algae_phylum_block_alpha_pielou)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_algae_phylum_block_alpha_pielou))
Anova(glmm_algae_phylum_block_alpha_pielou)

coef(summary(glmm_algae_phylum_block_alpha_pielou))
overdisp_fun(glmm_algae_phylum_block_alpha_pielou)
print(summary(glmm_algae_phylum_block_alpha_pielou), correlation=TRUE)
parameters::p_value(glmm_algae_phylum_block_alpha_pielou)
eta_squared(glmm_algae_phylum_block_alpha_pielou)



## For the algae order level

## For Observed
### To check the normality of the data
hist(alpha_lm_algae_order$Observed)
qqnorm(alpha_lm_algae_order$Observed)
qqline(alpha_lm_algae_order$Observed)
shapiro.test(alpha_lm_algae_order$Observed)
### Mixed model
lmm_algae_order_block_alpha_Observed <- lmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), data = alpha_lm_algae_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(lmm_algae_order_block_alpha_Observed)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(lmm_algae_order_block_alpha_Observed)~as.factor(alpha_lm_algae_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(lmm_algae_order_block_alpha_Observed, Observed~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_order$Observed~fitted(lmm_algae_order_block_alpha_Observed))

### Normality of Residuals
qqnorm(resid(lmm_algae_order_block_alpha_Observed))
qqline(resid(lmm_algae_order_block_alpha_Observed))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(lmm_algae_order_block_alpha_Observed)
# Distribution of the random effects
qqnorm (ranef(lmm_algae_order_block_alpha_Observed)$header_tanks_block[[1]])

plot(simulateResiduals(lmm_algae_order_block_alpha_Observed))
Anova(lmm_algae_order_block_alpha_Observed)

coef(summary(lmm_algae_order_block_alpha_Observed))
overdisp_fun(lmm_algae_order_block_alpha_Observed)
print(summary(lmm_algae_order_block_alpha_Observed), correlation=TRUE)
eta_squared(lmm_algae_order_block_alpha_Observed)


## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_order$Shannon)
qqnorm(alpha_lm_algae_order$Shannon)
qqline(alpha_lm_algae_order$Shannon)
shapiro.test(alpha_lm_algae_order$Shannon)

### Mixed model
glmm_algae_order_block_alpha_Shannon <- glmer(Shannon ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_algae_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_algae_order_block_alpha_Shannon)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_algae_order_block_alpha_Shannon)~as.factor(alpha_lm_algae_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_algae_order_block_alpha_Shannon, Shannon~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_order$Shannon~fitted(glmm_algae_order_block_alpha_Shannon))

### Normality of Residuals
qqnorm(resid(glmm_algae_order_block_alpha_Shannon))
qqline(resid(glmm_algae_order_block_alpha_Shannon))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_algae_order_block_alpha_Shannon)
# Distribution of the random effects
qqnorm (ranef(glmm_algae_order_block_alpha_Shannon)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_algae_order_block_alpha_Shannon))
Anova(glmm_algae_order_block_alpha_Shannon)

coef(summary(glmm_algae_order_block_alpha_Shannon))
overdisp_fun(glmm_algae_order_block_alpha_Shannon)
print(summary(glmm_algae_order_block_alpha_Shannon), correlation=TRUE)
parameters::p_value(glmm_algae_order_block_alpha_Shannon)
eta_squared(glmm_algae_order_block_alpha_Shannon)


## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_order$pielou)
qqnorm(alpha_lm_algae_order$pielou)
qqline(alpha_lm_algae_order$pielou)
shapiro.test(alpha_lm_algae_order$pielou)
### Mixed model
glmm_algae_order_block_alpha_pielou <- glmer(pielou ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_algae_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_algae_order_block_alpha_pielou)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_algae_order_block_alpha_pielou)~as.factor(alpha_lm_algae_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_algae_order_block_alpha_pielou, pielou~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_algae_order$pielou~fitted(glmm_algae_order_block_alpha_pielou))

### Normality of Residuals
qqnorm(resid(glmm_algae_order_block_alpha_pielou))
qqline(resid(glmm_algae_order_block_alpha_pielou))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_algae_order_block_alpha_pielou)
# Distribution of the random effects
qqnorm (ranef(glmm_algae_order_block_alpha_pielou)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_algae_order_block_alpha_pielou))
Anova(glmm_algae_order_block_alpha_pielou)

coef(summary(glmm_algae_order_block_alpha_pielou))
overdisp_fun(glmm_algae_order_block_alpha_pielou)
print(summary(glmm_algae_order_block_alpha_pielou), correlation=TRUE)
parameters::p_value(glmm_algae_order_block_alpha_pielou)
eta_squared(glmm_algae_order_block_alpha_pielou)


## For the protist phylum level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_phylum$Observed)
qqnorm(alpha_lm_protist_phylum$Observed)
qqline(alpha_lm_protist_phylum$Observed)
shapiro.test(alpha_lm_protist_phylum$Observed)
### Mixed model
glmm_protist_phylum_block_alpha_Observed <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_phylum_block_alpha_Observed)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_phylum_block_alpha_Observed)~as.factor(alpha_lm_protist_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_phylum_block_alpha_Observed, Observed~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_phylum$Observed~fitted(glmm_protist_phylum_block_alpha_Observed))

### Normality of Residuals
qqnorm(resid(glmm_protist_phylum_block_alpha_Observed))
qqline(resid(glmm_protist_phylum_block_alpha_Observed))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_phylum_block_alpha_Observed)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_phylum_block_alpha_Observed)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_phylum_block_alpha_Observed))
Anova(glmm_protist_phylum_block_alpha_Observed)

coef(summary(glmm_protist_phylum_block_alpha_Observed))
overdisp_fun(glmm_protist_phylum_block_alpha_Observed)
print(summary(glmm_protist_phylum_block_alpha_Observed), correlation=TRUE)
parameters::p_value(glmm_protist_phylum_block_alpha_Observed)
eta_squared(glmm_protist_phylum_block_alpha_Observed)


## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_phylum$Shannon)
qqnorm(alpha_lm_protist_phylum$Shannon)
qqline(alpha_lm_protist_phylum$Shannon)
shapiro.test(alpha_lm_protist_phylum$Shannon)
### Mixed model
glmm_protist_phylum_block_alpha_Shannon <- glmer(Shannon ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_phylum_block_alpha_Shannon)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_phylum_block_alpha_Shannon)~as.factor(alpha_lm_protist_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_phylum_block_alpha_Shannon, Shannon~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_phylum$Shannon~fitted(glmm_protist_phylum_block_alpha_Shannon))

### Normality of Residuals
qqnorm(resid(glmm_protist_phylum_block_alpha_Shannon))
qqline(resid(glmm_protist_phylum_block_alpha_Shannon))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_phylum_block_alpha_Shannon)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_phylum_block_alpha_Shannon)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_phylum_block_alpha_Shannon))
Anova(glmm_protist_phylum_block_alpha_Shannon)

coef(summary(glmm_protist_phylum_block_alpha_Shannon))
overdisp_fun(glmm_protist_phylum_block_alpha_Shannon)
print(summary(glmm_protist_phylum_block_alpha_Shannon), correlation=TRUE)
parameters::p_value(glmm_protist_phylum_block_alpha_Shannon)
eta_squared(glmm_protist_phylum_block_alpha_Shannon)


## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_phylum$pielou)
qqnorm(alpha_lm_protist_phylum$pielou)
qqline(alpha_lm_protist_phylum$pielou)
shapiro.test(alpha_lm_protist_phylum$pielou)
### Mixed model
glmm_protist_phylum_block_alpha_pielou <- glmer(pielou ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_phylum)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_phylum_block_alpha_pielou)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_phylum_block_alpha_pielou)~as.factor(alpha_lm_protist_phylum$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_phylum_block_alpha_pielou, pielou~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_phylum$pielou~fitted(glmm_protist_phylum_block_alpha_pielou))

### Normality of Residuals
qqnorm(resid(glmm_protist_phylum_block_alpha_pielou))
qqline(resid(glmm_protist_phylum_block_alpha_pielou))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_phylum_block_alpha_pielou)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_phylum_block_alpha_pielou)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_phylum_block_alpha_pielou))
Anova(glmm_protist_phylum_block_alpha_pielou)

coef(summary(glmm_protist_phylum_block_alpha_pielou))
overdisp_fun(glmm_protist_phylum_block_alpha_pielou)
print(summary(glmm_protist_phylum_block_alpha_pielou), correlation=TRUE)
parameters::p_value(glmm_protist_phylum_block_alpha_pielou)
eta_squared(glmm_protist_phylum_block_alpha_pielou)


## For the protist order level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_order$Observed)
qqnorm(alpha_lm_protist_order$Observed)
qqline(alpha_lm_protist_order$Observed)
shapiro.test(alpha_lm_protist_order$Observed)
### Mixed model
glmm_protist_order_block_alpha_Observed <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_order_block_alpha_Observed)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_order_block_alpha_Observed)~as.factor(alpha_lm_protist_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_order_block_alpha_Observed, Observed~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_order$Observed~fitted(glmm_protist_order_block_alpha_Observed))

### Normality of Residuals
qqnorm(resid(glmm_protist_order_block_alpha_Observed))
qqline(resid(glmm_protist_order_block_alpha_Observed))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_order_block_alpha_Observed)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_order_block_alpha_Observed)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_order_block_alpha_Observed))
Anova(glmm_protist_order_block_alpha_Observed)

coef(summary(glmm_protist_order_block_alpha_Observed))
overdisp_fun(glmm_protist_order_block_alpha_Observed)
print(summary(glmm_protist_order_block_alpha_Observed), correlation=TRUE)
parameters::p_value(glmm_protist_order_block_alpha_Observed)
eta_squared(glmm_protist_order_block_alpha_Observed)


## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_order$Shannon)
qqnorm(alpha_lm_protist_order$Shannon)
qqline(alpha_lm_protist_order$Shannon)
shapiro.test(alpha_lm_protist_order$Shannon)
### Mixed model
glmm_protist_order_block_alpha_Shannon <- glmer(Shannon ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_order_block_alpha_Shannon)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_order_block_alpha_Shannon)~as.factor(alpha_lm_protist_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_order_block_alpha_Shannon, Shannon~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_order$Shannon~fitted(glmm_protist_order_block_alpha_Shannon))

### Normality of Residuals
qqnorm(resid(glmm_protist_order_block_alpha_Shannon))
qqline(resid(glmm_protist_order_block_alpha_Shannon))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_order_block_alpha_Shannon)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_order_block_alpha_Shannon)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_order_block_alpha_Shannon))
Anova(glmm_protist_order_block_alpha_Shannon)

coef(summary(glmm_protist_order_block_alpha_Shannon))
overdisp_fun(glmm_protist_order_block_alpha_Shannon)
print(summary(glmm_protist_order_block_alpha_Shannon), correlation=TRUE)
parameters::p_value(glmm_protist_order_block_alpha_Shannon)
eta_squared(glmm_protist_order_block_alpha_Shannon)


## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_order$pielou)
qqnorm(alpha_lm_protist_order$pielou)
qqline(alpha_lm_protist_order$pielou)
shapiro.test(alpha_lm_protist_order$pielou)
### Mixed model
glmm_protist_order_block_alpha_pielou <- glmer(pielou ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = gaussian	(link = "identity"), data = alpha_lm_protist_order)

### To asses diagnostics with mixed-effect models

#### Test for Heteroscedasticity
plot(glmm_protist_order_block_alpha_pielou)
# Do header_tanks_block of week 2 consistently have higher variance in residuals than header_tanks_block of week 3?
plot(resid(glmm_protist_order_block_alpha_pielou)~as.factor(alpha_lm_protist_order$header_tanks_block))
abline(h=0)
# The overall fitness of the model
plot(glmm_protist_order_block_alpha_pielou, pielou~fitted(.), id=0.05, adj=-0.3)
plot(alpha_lm_protist_order$pielou~fitted(glmm_protist_order_block_alpha_pielou))

### Normality of Residuals
qqnorm(resid(glmm_protist_order_block_alpha_pielou))
qqline(resid(glmm_protist_order_block_alpha_pielou))

### Assumptions about the random effects
# Prediction of random efefcts from the model
nlme::ranef(glmm_protist_order_block_alpha_pielou)
# Distribution of the random effects
qqnorm (ranef(glmm_protist_order_block_alpha_pielou)$header_tanks_block[[1]])

plot(simulateResiduals(glmm_protist_order_block_alpha_pielou))
Anova(glmm_protist_order_block_alpha_pielou)

coef(summary(glmm_protist_order_block_alpha_pielou))
overdisp_fun(glmm_protist_order_block_alpha_pielou)
print(summary(glmm_protist_order_block_alpha_pielou), correlation=TRUE)
parameters::p_value(glmm_protist_order_block_alpha_pielou)
eta_squared(glmm_protist_order_block_alpha_pielou)


# BETA-DIVERSITY INDICES


# According to Ollberding (2019), this statistical analysis was performed.

# Permutational multivariate analysis of variance (PERMANOVA)
# Followed by https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova

## For the algae phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_algae_phylum <- phyloseq::distance(alph_alg_sset_phy_1, method = "bray")

## To make a data frame from the sample_data
perm_algae_phylum_sampledf <- data.frame(sample_data(alph_alg_sset_phy_1))

## To perform the adonis test
adonis2(bray_algae_phylum ~ nutrient*sediment*flow*time, strata = perm_algae_phylum_sampledf$week, data = perm_algae_phylum_sampledf)

## For the algae order level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_algae_order <- phyloseq::distance(alph_alg_sset_ord_1, method = "bray")

## To make a data frame from the sample_data
perm_algae_order_sampledf <- data.frame(sample_data(alph_alg_sset_ord_1))

## To perform the adonis test
adonis2(bray_algae_order ~ nutrient*sediment*flow*time, strata = perm_algae_order_sampledf$week, data = perm_algae_order_sampledf)

## For the protist phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_phylum <- phyloseq::distance(alph_prot_sset_phy_1, method = "bray")

## To make a data frame from the sample_data
perm_protist_phylum_sampledf <- data.frame(sample_data(alph_prot_sset_phy_1))

## To perform the adonis test
adonis2(bray_protist_phylum ~ nutrient*sediment*flow*time, strata = perm_protist_phylum_sampledf$week, data = perm_protist_phylum_sampledf)

## For the protist order level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_order <- phyloseq::distance(alph_prot_sset_ord_1, method = "bray")

## To make a data frame from the sample_data
perm_protist_order_sampledf <- data.frame(sample_data(alph_prot_sset_ord_1))

## To perform the adonis test
adonis2(bray_protist_order ~ nutrient*sediment*flow*time, strata = perm_protist_order_sampledf$week, data = perm_protist_order_sampledf)


# Permutational analysis of multivariate dispersion (PERMDISP)

## For the algae phylum level
dispersion_algae_phylum <- vegan::betadisper(bray_algae_phylum, phyloseq::sample_data(alph_alg_sset_phy_1)$treatment)
dispersion_algae_phylum

plot(dispersion_algae_phylum, main = "For the algae phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_phylum, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_phylum)

## For the algae order level
dispersion_algae_order <- vegan::betadisper(bray_algae_order, phyloseq::sample_data(alph_alg_sset_ord_1)$treatment)
dispersion_algae_order

plot(dispersion_algae_order, main = "For the algae order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_order, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_order)

## For the protist phylum level
dispersion_protist_phylum <- vegan::betadisper(bray_protist_phylum, phyloseq::sample_data(alph_prot_sset_phy_1)$treatment)
dispersion_protist_phylum

plot(dispersion_protist_phylum, main = "For the protist phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_phylum, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_phylum)

## For the protist order level
dispersion_protist_order <- vegan::betadisper(bray_protist_order, phyloseq::sample_data(alph_prot_sset_ord_1)$treatment)
dispersion_protist_order

plot(dispersion_protist_order, main = "For the protist order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_order, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_order)



# The PCoA analysis (Principal Coordinate Analysis)
# Following by https://joey711.github.io/phyloseq/plot_ordination-examples.html, this statistical analysis was performed.

pcoa_dist = "bray"
pcoa_ord_meths = c("PCoA")

## For the algae phylum level
ap_plist = llply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, alph_alg_sset_phy_1, pcoa_dist)

names(ap_plist) <- pcoa_ord_meths

ap_pdataframe = ldply(ap_plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("PCoA1", "PCoA2")
  return(cbind(df, x$data))
})
names(ap_pdataframe)[1] = "method"

ap_p = ggplot(ap_pdataframe, aes(PCoA1, PCoA2, color=treatment, shape=as.factor(week), fill=treatment))+
  geom_point(size=4) + geom_polygon()+ 
  facet_wrap(~method, scales="free")+ 
  scale_fill_brewer(type="qual", palette="Set1")+ 
  scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA1(27.34%) and PCoA2(20.96%)")
ap_p

## For the algae order level
ao_plist = llply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, alph_alg_sset_ord_1, pcoa_dist)

names(ao_plist) <- pcoa_ord_meths

ao_pdataframe = ldply(ao_plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("PCoA1", "PCoA2")
  return(cbind(df, x$data))
})
names(ao_pdataframe)[1] = "method"

ao_p = ggplot(ao_pdataframe, aes(PCoA1, PCoA2, color=treatment, shape=as.factor(week), fill=treatment))+
  geom_point(size=4) + geom_polygon()+ 
  facet_wrap(~method, scales="free")+ 
  scale_fill_brewer(type="qual", palette="Set1")+ 
  scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA1(26.51%) and PCoA2(22.78%)")
ao_p

## For the protist phylum level
pp_plist = llply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, alph_prot_sset_phy_1, pcoa_dist)

names(pp_plist) <- pcoa_ord_meths

pp_pdataframe = ldply(pp_plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("PCoA1", "PCoA2")
  return(cbind(df, x$data))
})
names(pp_pdataframe)[1] = "method"

pp_p = ggplot(pp_pdataframe, aes(PCoA1, PCoA2, color=treatment, shape=as.factor(week), fill=treatment))+
  geom_point(size=4) + geom_polygon()+ 
  facet_wrap(~method, scales="free")+ 
  scale_fill_brewer(type="qual", palette="Set1")+ 
  scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA1(24.01%) and PCoA2(19.02%)")
pp_p

## For the protist order level
po_plist = llply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, alph_prot_sset_ord_1, pcoa_dist)

names(po_plist) <- pcoa_ord_meths

po_pdataframe = ldply(po_plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("PCoA1", "PCoA2")
  return(cbind(df, x$data))
})
names(po_pdataframe)[1] = "method"

po_p = ggplot(po_pdataframe, aes(PCoA1, PCoA2, color=treatment, shape=as.factor(week), fill=treatment))+
  geom_point(size=4) + geom_polygon()+ 
  facet_wrap(~method, scales="free")+ 
  scale_fill_brewer(type="qual", palette="Set1")+ 
  scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA1(24.97%) and PCoA2(18.92%)")
po_p


# According to joey711 (2022), the most abundant taxa was calculated. 
## To find the most abundant taxa
## Function
most_abundant_taxa <- function(x,taxa){
  require(phyloseq)
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa))
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- data.frame(tax[k,])
  m <- data.frame(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
 
## Keeping >1% of total reads based on abundant taxa at phylum and order level of algae and protists
subset_algae_phylum <- subset_taxa(phyloseq_object_tsc, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Streptophyta", "Xanthophyceae", "Nucleariidae_and_Fonticula_group", "Phragmoplastophyta"))
prune_algae_phylum <- prune_taxa(taxa_sums(subset_algae_phylum) > 0.01, subset_algae_phylum)
ap_clr <- transform(prune_algae_phylum, 'clr')
abundant_algae_phylum <- most_abundant_taxa(ap_clr,"phylum")
abundant_algae_phylum

subset_algae_order <- subset_taxa(phyloseq_object_tsc, order %in% c("Bacillariales", "Batrachospermales", "Cercomonadida", "Chaetopeltidales", "Chaetophorales", "Chlamydomonadales", "Chlorellales", "Chlorodendrales", "Chlorosarcinales", "Chromulinales", "Coccolithales", "Cocconeidales", "Cryptomonadales", "Cymbellales", "Desmidiales", "Dinophysiales", "Eustigmatales", "Fragilariales", "Gymnodiniales","Gymnodiniphycidae", "Hibberdiales", "Klebsormidiales", "Laminariales","Licmophorales", "Lophodiniales", "Marsupiomonadales", "Mastogloiales","Melosirales", "Microthamniales","Mischococcales", "Monomastigales","Naviculales", "Ochromonadales", "Oedogoniales", "Pavlovales","Pedinomonadales", "Peridiniales","Phytodiniales", "Pyrenomonadales","Scotinosphaerales", "Sphaeropleales","Suessiales", "Surirellales","Synurales", "Tetrasporales","Thalassiophysales", "Thalassiosirales", "Thaumatomonadida", "Thoreales","Trentepohliales", "Ulotrichales","Ulvales", "Zygnematales"))
prune_algae_order <- prune_taxa(taxa_sums(subset_algae_order) > 0.01, subset_algae_order)
ao_clr <- transform(prune_algae_order, 'clr')
abundant_algae_order <- most_abundant_taxa(ao_clr,"order")
abundant_algae_order

subset_protist_phylum <- subset_taxa(phyloseq_object_tsc, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Tubulinea", "Bicosoecida", "Discosea", "Labyrinthulomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida", "Colponemidia", "MAST-12", "Rigifilida", "Schizoplasmodiida"))
prune_protist_phylum <- prune_taxa(taxa_sums(subset_protist_phylum) > 0.01, subset_protist_phylum)
pp_clr <- transform(prune_protist_phylum, 'clr')
abundant_protist_phylum <- most_abundant_taxa(pp_clr,"phylum")
abundant_protist_phylum

subset_protist_order <- subset_taxa(phyloseq_object_tsc, order %in% c("Arcellinida",  "Bicosoecida", "Bryometopida", "Chlamydodontida", "Choanoflagellida", "Choreotrichida", "Colpodida", "Conthreep", "Cryomonadida", "Cryptosporida", "Cyrtolophosidida", "Dysteriida", "Echinamoebida", "Eucoccidiorida", "Euglyphida", "Eugregarinorida", "Glissomonadida", "Grossglockneriida", "Haptorida","Heterotrichida", "Hymenostomatida","Leptomyxida", "Litostomatea","Longamoebia", "Microthoracida","Nassulida", "Neogregarinorida","Odontostomatida", "Peniculida","Philasterida", "Pleurostomatida", "Prorodontida", "Protosteliales","Salpingoecidae", "Silicofilosea","Spirotrichea", "Spongomonadida","Sporadotrichida", "Stichotrichida","Tintinnida", "Urostylida","Vampyrellida"))
prune_protist_order <- prune_taxa(taxa_sums(subset_protist_order) > 0.01, subset_protist_order)
po_clr <- transform(prune_protist_order, 'clr')
abundant_protist_order <- most_abundant_taxa(po_clr,"order")
abundant_protist_order

# Multivariate Analysis Of Variance (MANOVA)
## According to Kassambara (2017); Ben-Shachar et al., (2020), this statistical analysis was performed.

## For the most abundant taxa of the algal phylum level
manova_df_algae_phylum = as(sample_data(ap_clr), "data.frame")
manova_df_algae_phylum
manova_df_algae_phylum_1 = cbind(manova_df_algae_phylum, abundant_algae_phylum)
manova_df_algae_phylum_1
### For Bacillariophyta, Eustigmatophyceae, Phragmoplastophyta, Streptophyta and Chlorophyta 
shapiro.test(manova_df_algae_phylum_1$Bacillariophyta)

shapiro.test(manova_df_algae_phylum_1$Eustigmatophyceae)
Eustigmatophyceae_1 <- sqrt(manova_df_algae_phylum_1$Eustigmatophyceae)
shapiro.test(Eustigmatophyceae_1)

shapiro.test(manova_df_algae_phylum_1$Phragmoplastophyta)
shapiro.test(manova_df_algae_phylum_1$Streptophyta)
shapiro.test(manova_df_algae_phylum_1$Chlorophyta)

manova_df_algae_phylum_11 <- cbind(manova_df_algae_phylum_1, Eustigmatophyceae_1)
manova_df_algae_phylum_12 <- manova_df_algae_phylum_11 %>% 
  mutate_at(vars(Eustigmatophyceae_1), ~replace_na(.,mean(., na.rm = TRUE)))

manova_abundant_algae_phylum <- manova(cbind(Bacillariophyta, Eustigmatophyceae_1, Phragmoplastophyta, Streptophyta, Chlorophyta) ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_12)
manova_abundant_algae_phylum
coef(manova_abundant_algae_phylum)
summary.aov(manova_abundant_algae_phylum)
eta_squared(aov(manova_abundant_algae_phylum))

## For the most abundant taxa of the algal order level
manova_df_algae_order = as(sample_data(ao_clr), "data.frame")
manova_df_algae_order
manova_df_algae_order_1 = cbind(manova_df_algae_order, abundant_algae_order)
manova_df_algae_order_1

### For Cymbellales, Sphaeropleales, Cocconeidales, Chlamydomonadales, Cercomonadida, Ulotrichales and Melosirales
shapiro.test(manova_df_algae_order_1$Cymbellales)

shapiro.test(manova_df_algae_order_1$Sphaeropleales)
Sphaeropleales_1 <- sqrt(manova_df_algae_order_1$Sphaeropleales)
shapiro.test(Sphaeropleales_1)

shapiro.test(manova_df_algae_order_1$Cocconeidales)
Cocconeidales_1 <- sqrt(manova_df_algae_order_1$Cocconeidales)
shapiro.test(Cocconeidales_1)

shapiro.test(manova_df_algae_order_1$Chlamydomonadales)

shapiro.test(manova_df_algae_order_1$Cercomonadida)
Cercomonadida_1 <- sqrt(manova_df_algae_order_1$Cercomonadida)
shapiro.test(Cercomonadida_1)

shapiro.test(manova_df_algae_order_1$Ulotrichales)
Ulotrichales_1 <- log(manova_df_algae_order_1$Ulotrichales)
shapiro.test(Ulotrichales_1)

shapiro.test(manova_df_algae_order_1$Melosirales)
Melosirales_1 <- sqrt(manova_df_algae_order_1$Melosirales)
shapiro.test(Melosirales_1)

manova_df_algae_order_11 <- cbind(manova_df_algae_order_1, Sphaeropleales_1, Cocconeidales_1, Cercomonadida_1, Ulotrichales_1, Melosirales_1)
manova_df_algae_order_12 <- manova_df_algae_order_11 %>% 
  mutate_at(vars(Sphaeropleales_1, Cocconeidales_1, Cercomonadida_1, Ulotrichales_1, Melosirales_1), ~replace_na(.,mean(., na.rm = TRUE)))

manova_abundant_algae_order <- manova(cbind(Cymbellales, Sphaeropleales_1, Cocconeidales_1, Chlamydomonadales, Cercomonadida_1, Ulotrichales_1, Melosirales_1) ~ nutrient*sediment*flow*time, data = manova_df_algae_order_12)
manova_abundant_algae_order
coef(manova_abundant_algae_order)
summary.aov(manova_abundant_algae_order)
eta_squared(aov(manova_abundant_algae_order))

## For the most abundant taxa of the protist phylum level
manova_df_protist_phylum = as(sample_data(pp_clr), "data.frame")
manova_df_protist_phylum
manova_df_protist_phylum_1 = cbind(manova_df_protist_phylum, abundant_protist_phylum)
manova_df_protist_phylum_1

### For Apicomplexa, Cercozoa,Tubulinea, Ciliophora and Labyrinthulomycetes
shapiro.test(manova_df_protist_phylum_1$Apicomplexa)
shapiro.test(manova_df_protist_phylum_1$Cercozoa)
shapiro.test(manova_df_protist_phylum_1$Tubulinea)

shapiro.test(manova_df_protist_phylum_1$Ciliophora)
Ciliophora_1 <- sqrt(manova_df_protist_phylum_1$Ciliophora)
shapiro.test(Ciliophora_1)

shapiro.test(manova_df_protist_phylum_1$Labyrinthulomycetes)
Labyrinthulomycetes_1 <- sqrt(manova_df_protist_phylum_1$Labyrinthulomycetes)
shapiro.test(Labyrinthulomycetes_1)

manova_df_protist_phylum_11 <- cbind(manova_df_protist_phylum_1, Ciliophora_1, Labyrinthulomycetes_1)
manova_df_protist_phylum_12 <- manova_df_protist_phylum_11 %>% 
  mutate_at(vars(Ciliophora_1, Labyrinthulomycetes_1), ~replace_na(.,mean(., na.rm = TRUE)))

manova_abundant_protist_phylum <- manova(cbind(Apicomplexa, Cercozoa, Tubulinea, Ciliophora_1, Labyrinthulomycetes_1) ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_12)
manova_abundant_protist_phylum
coef(manova_abundant_protist_phylum)
summary.aov(manova_abundant_protist_phylum)
eta_squared(aov(manova_abundant_protist_phylum))

## For the most abundant taxa of the protist order level
manova_df_protist_order = as(sample_data(po_clr), "data.frame")
manova_df_protist_order
manova_df_protist_order_1 = cbind(manova_df_protist_order, abundant_protist_order)
manova_df_protist_order_1

### For Pleurostomatida, Eucoccidiorida, Cryomonadida, Conthreep, Urostylida, Echinamoebida, Euglyphida, Glissomonadida, Vampyrellida, Sporadotrichida, Stichotrichida, Cyrtolophosidida and Silicofilosea
shapiro.test(manova_df_protist_order_1$Pleurostomatida)
Pleurostomatida_1 <- sqrt(manova_df_protist_order_1$Pleurostomatida)
shapiro.test(Pleurostomatida_1)

shapiro.test(manova_df_protist_order_1$Eucoccidiorida)
Eucoccidiorida_1 <- sqrt(manova_df_protist_order_1$Eucoccidiorida)
shapiro.test(Eucoccidiorida_1)

shapiro.test(manova_df_protist_order_1$Cryomonadida)

shapiro.test(manova_df_protist_order_1$Conthreep)
Conthreep_1 <- sqrt(manova_df_protist_order_1$Conthreep)
shapiro.test(Conthreep_1)

shapiro.test(manova_df_protist_order_1$Urostylida)

shapiro.test(manova_df_protist_order_1$Echinamoebida)
Echinamoebida_1 <- log(manova_df_protist_order_1$Echinamoebida)
shapiro.test(Echinamoebida_1)

shapiro.test(manova_df_protist_order_1$Euglyphida)

shapiro.test(manova_df_protist_order_1$Glissomonadida)
Glissomonadida_1 <- sqrt(manova_df_protist_order_1$Glissomonadida)
shapiro.test(Glissomonadida_1)

shapiro.test(manova_df_protist_order_1$Vampyrellida)

shapiro.test(manova_df_protist_order_1$Sporadotrichida)

shapiro.test(manova_df_protist_order_1$Stichotrichida)
Stichotrichida_1 <- sqrt(manova_df_protist_order_1$Stichotrichida)
shapiro.test(Stichotrichida_1)

shapiro.test(manova_df_protist_order_1$Cyrtolophosidida)
shapiro.test(manova_df_protist_order_1$Silicofilosea)

manova_df_protist_order_11 <- cbind(manova_df_protist_order_1, Pleurostomatida_1, Eucoccidiorida_1, Conthreep_1, Echinamoebida_1, Glissomonadida_1, Stichotrichida_1)
manova_df_protist_order_12 <- manova_df_protist_order_11 %>% 
  mutate_at(vars(Pleurostomatida_1, Eucoccidiorida_1, Conthreep_1, Echinamoebida_1, Glissomonadida_1, Stichotrichida_1), ~replace_na(.,mean(., na.rm = TRUE)))


manova_abundant_protist_order <- manova(cbind(Pleurostomatida_1, Eucoccidiorida_1, Cryomonadida, Conthreep_1, Urostylida, Echinamoebida_1, Euglyphida, Glissomonadida_1, Vampyrellida, Sporadotrichida, Stichotrichida_1, Cyrtolophosidida, Silicofilosea) ~ nutrient*sediment*flow*time, data = manova_df_protist_order_12)
manova_abundant_protist_order
coef(manova_abundant_protist_order)
summary.aov(manova_abundant_protist_order)
eta_squared(aov(manova_abundant_protist_order))






List of references

Callahan, Benjamin. (2018). Silva taxonomic training data formatted for DADA2 (Silva version 132) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.1172783

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." _Nature Methods_, *13*, 581-583. doi:10.1038/nmeth.3869 <https://doi.org/10.1038/nmeth.3869>.

Callahan, Benjamin. (2018). DADA2 Pipeline Tutorial (1.8). dada2. https://benjjneb.github.io/dada2/tutorial_1_8.html

Vaulot, D. (2021, Feb 15). Phyloseq tutorial. https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." _PLoS ONE_, *8*(4), e61217. <http://dx.plos.org/10.1371/journal.pone.0061217>.

Martins, C., & Moreau, C. S. (2020). Influence of host phylogeny, geographical location and seed harvesting diet on the bacterial community of globally distributed Pheidole ants. PeerJ, 8, e8492.

Mariadassou, M., Bernard, M., Pascal, G., Cauquil, L., & Chaillou, S. (2016). Analysis of community composition data using phyloseq. Montpellier Décembre 2016. https://genoweb.toulouse.inra.fr/~formation/15_FROGS/8-February2017/FROGS_phyloseq_02_2017.pdf

Gauthier, J., & Derome, N. (2021). Evenness-Richness Scatter Plots: a Visual and Insightful Representation of Shannon Entropy Measurements for Ecological Community Analysis. Msphere, 6(2), e01019-20.

ZACH. (2021, Sep 29). How to Test for Normality in R (4 Methods). Statology. https://www.statology.org/test-for-normality-in-r/
  
Humboldt-Universität zu Berlin | Geography Department. (2021). Generalized linear models in R. Quantitative Methods for Geographers. https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab08_GLM1.html

Phillips, N.D. (2018, Jan 22). 15.4 Regression on non-Normal data with glm(). YaRrr! The Pirate's Guide to R. https://bookdown.org/ndphillips/YaRrr/regression-on-non-normal-data-with-glm.html

Bolker, B. (2018, Sep 25). GLMM worked examples. https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html

Walker, J.A. (2021, Sep 25). Applied Statistics for Experimental Biology. Elements of Applied Biostatistics. https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html

Ollberding, N.J. (2019, Jul 28). Introduction to the Statistical Analysis of Microbiome Data in R. Wowchemy. https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

joey711. (2022). Find the Most abundant Taxa in individual samples #847. https://github.com/joey711/phyloseq/issues/847

Kassambara, A. (2017). MANOVA Test in R: Multivariate Analysis of Variance. Statistical tools for high-throughput data analysis. STHDA. http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance#infos

Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815

Niku, J., Hui, F. K., Taskinen, S., & Warton, D. I. (2019). gllvm: Fast analysis of multivariate abundance data with generalized linear latent variable models in r. Methods in Ecology and Evolution, 10(12), 2173-2182.

Niku J, Hui FKC, Taskinen S, Warton DI (2019). "gllvm - Fast analysis of multivariate abundance data with generalized linear latent variable models in R." _Methods in Ecology and Evolution_, *10*, 2173-2182.

Lahti, L., & Shetty, S. (Bioconductor, 2017). Tools for microbiome analysis in R. Microbiome package version 1.17.43. URL: https://github.com/microbiome/microbiome/blob/master/R/coverage.R

