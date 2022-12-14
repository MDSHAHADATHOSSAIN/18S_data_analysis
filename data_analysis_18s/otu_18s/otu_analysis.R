# Analysis of eukaryotic biodiversity data using the phyloseq package in R

# OTU data analysis with the phyloseq package for statistical analysis

## According to Vaulot (2021) and McMurdie & Holmes (2013), this statistical analysis was performed.
## To install necessary packages
install.packages("pacman")
install.packages("dplyr") 
install.packages("readxl")    
install.packages("ggplot2")

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
install.packages("lmerTest")
install.packages("DHARMa")
install.packages("car")
install.packages("nlme")
install.packages("stats")

## To load the installed packages with the pacman package
pacman::p_load(pacman, dplyr, readxl, ggplot2, phyloseq, microbiome, reshape2, ape, gridExtra, plotly, vegan, dendextend, tidyr, rms, effectsize, lme4, picante, cowplot, here, lsr, MASS, rcompanion, ggiraphExtra, optimx, mvabund, gllvm, radiant.data, tidyverse, MuMIn, corrplot, gclus, coefplot, jtools, ggstance, interactions, blmeco, tseries, lmerTest, DHARMa, car, nlme, stats)
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

## To keep only taxa of interests (for Algae) according to phylum and order level
sub_algae_phylum <- subset_taxa(phyloseq_object, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Streptophyta", "Xanthophyceae", "Phragmoplastophyta", "unclassified", "unclassified_unidentified"))
sub_algae_order <- subset_taxa(phyloseq_object, order %in% c("Bacillariales", "Batrachospermales", "Cercomonadida", "Chaetopeltidales", "Chaetophorales", "Chlamydomonadales", "Chlorellales", "Chlorodendrales", "Chlorosarcinales", "Chromulinales", "Coccolithales", "Cocconeidales", "Cryptomonadales", "Cymbellales", "Desmidiales", "Dinophysiales", "Eustigmatales", "Fragilariales", "Gymnodiniales","Gymnodiniphycidae", "Hibberdiales", "Klebsormidiales", "Laminariales","Licmophorales", "Lophodiniales", "Marsupiomonadales", "Mastogloiales","Melosirales", "Microthamniales","Mischococcales", "Monomastigales","Naviculales", "Ochromonadales", "Oedogoniales", "Pavlovales","Pedinomonadales", "Peridiniales","Phytodiniales", "Pyrenomonadales","Scotinosphaerales", "Sphaeropleales","Suessiales", "Surirellales","Synurales", "Tetrasporales","Thalassiophysales", "Thalassiosirales", "Thaumatomonadida", "Thoreales","Trentepohliales", "Ulotrichales","Ulvales", "Zygnematales", "unclassified", "unclassified_Bacillariophyta", "unclassified_Bicosoecida", "unclassified_Cercomonadidae", "unclassified_Chlorophyceae", "unclassified_Chlorophyta", "unclassified_Chrysophyceae", "unclassified_Clade_L", "unclassified_Colponemidia", "unclassified_Cryptophyta", "unclassified_Dinophyceae", "unclassified_Eustigmatophyceae", "unclassified_Katablepharidophyta", "unclassified_Nephroselmidophyceae", "unclassified_Trebouxiophyceae", "unclassified_Ulvophyceae", "unclassified_unidentified"))

## To keep only taxa of interests (for Protists) according to phylum and order level
sub_protist_phylum <- subset_taxa(phyloseq_object, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Tubulinea", "Bicosoecida", "Discosea", "Labyrinthulomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida", "Colponemidia", "MAST-12", "Rigifilida", "Schizoplasmodiida", "unclassified", "unclassified_unidentified", "Nucleariidae_and_Fonticula_group", "unclassified_Eukaryota"))
sub_protist_order <- subset_taxa(phyloseq_object, order %in% c("Arcellinida",  "Bicosoecida", "Bryometopida", "Chlamydodontida", "Choanoflagellida", "Choreotrichida", "Colpodida", "Conthreep", "Cryomonadida", "Cryptosporida", "Cyrtolophosidida", "Dysteriida", "Echinamoebida", "Eucoccidiorida", "Euglyphida", "Eugregarinorida", "Glissomonadida", "Grossglockneriida", "Haptorida","Heterotrichida", "Hymenostomatida","Leptomyxida", "Litostomatea","Longamoebia", "Microthoracida","Nassulida", "Neogregarinorida","Odontostomatida", "Peniculida","Philasterida", "Pleurostomatida", "Prorodontida", "Protosteliales","Salpingoecidae", "Silicofilosea","Spirotrichea", "Spongomonadida","Sporadotrichida", "Stichotrichida","Tintinnida", "Urostylida","Vampyrellida", "unclassified_Aconoidasida", "unclassified_Apicomplexa", "unclassified_Apusomonadidae", "unclassified_Arcellinida", "unclassified_Cercozoa", "unclassified_Conoidasida", "unclassified_CV1-B1-93", "unclassified_Discosea", "unclassified_Euamoebida", "unclassified_Glissomonadida", "unclassified_Gracilipodida", "unclassified_Imbricatea", "unclassified_Labyrinthulomycetes", "unclassified_LEMD267", "unclassified_Leptomyxida", "unclassified_MAST-12C", "unclassified_Nassophorea", "unclassified_Novel_Clade_Gran-3", "unclassified_Oligohymenophorea", "unclassified_Protosporangiidae", "unclassified_Schizoplasmodiida", "unclassified_Spirotrichea", "unclassified_Thecofilosea", "unclassified", "unclassified_unidentified"))


# The pre-processing of the data

## For the algae phylum level
## To transform the data from sub_algae_phylum to relative abundance (to transform the count of samples to percent)
sub_algae_phylum_tsc <- transform_sample_counts(sub_algae_phylum, function(x) x / sum(x))
sub_algae_phylum_tsc
## Sum of the values in each column of the otu_table before removing noises on the data
sub_algae_phylum_tsc_a <- as.data.frame(otu_table(sub_algae_phylum_tsc))
sub_algae_phylum_tsc_a
colSums(sub_algae_phylum_tsc_a)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_algae_phylum_tsc_minthreshold  = transform_sample_counts(sub_algae_phylum, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the otu_table after removing noises on the data
sub_algae_phylum_tsc_minthreshold_b <- as.data.frame(otu_table(sub_algae_phylum_tsc_minthreshold))
sub_algae_phylum_tsc_minthreshold_b
colSums(sub_algae_phylum_tsc_minthreshold_b)

## For the algae order level
## To transform the data from sub_algae_order to relative abundance (to transform the count of samples to percent)
sub_algae_order_tsc <- transform_sample_counts(sub_algae_order, function(x) x / sum(x))
sub_algae_order_tsc
## Sum of the values in each column of the otu_table before removing noises on the data
sub_algae_order_tsc_a <- as.data.frame(otu_table(sub_algae_order_tsc))
sub_algae_order_tsc_a
colSums(sub_algae_phylum_tsc_a)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_algae_order_tsc_minthreshold  = transform_sample_counts(sub_algae_order, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the otu_table after removing noises on the data
sub_algae_order_tsc_minthreshold_b <- as.data.frame(otu_table(sub_algae_order_tsc_minthreshold))
sub_algae_order_tsc_minthreshold_b
colSums(sub_algae_order_tsc_minthreshold_b)

## For the protist phylum level
## To transform the data from sub_protist_phylum to relative abundance (to transform the count of samples to percent)
sub_protist_phylum_tsc <- transform_sample_counts(sub_protist_phylum, function(x) x / sum(x))
sub_protist_phylum_tsc

as.data.frame(otu_table(sub_protist_phylum_tsc))

## Sum of the values in each column of the otu_table before removing noises on the data
sub_protist_phylum_tsc_a <- as.data.frame(otu_table(sub_protist_phylum_tsc))
sub_protist_phylum_tsc_a
colSums(sub_protist_phylum_tsc_a)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_protist_phylum_tsc_minthreshold  = transform_sample_counts(sub_protist_phylum, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})

as.data.frame(otu_table(sub_protist_phylum_tsc_minthreshold))

## Sum of the values in each column of the otu_table after removing noises on the data
sub_protist_phylum_tsc_minthreshold_b <- as.data.frame(otu_table(sub_protist_phylum_tsc_minthreshold))
sub_protist_phylum_tsc_minthreshold_b
colSums(sub_protist_phylum_tsc_minthreshold_b)

## For the protist order level
## To transform the data from sub_protist_order to relative abundance (to transform the count of samples to percent)
sub_protist_order_tsc <- transform_sample_counts(sub_protist_order, function(x) x / sum(x))
sub_protist_order_tsc
## Sum of the values in each column of the otu_table before removing noises on the data
sub_protist_order_tsc_a <- as.data.frame(otu_table(sub_protist_order_tsc))
sub_protist_order_tsc_a
colSums(sub_protist_order_tsc_a)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_protist_order_tsc_minthreshold  = transform_sample_counts(sub_protist_order, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the otu_table after removing noises on the data
sub_protist_order_tsc_minthreshold_b <- as.data.frame(otu_table(sub_protist_order_tsc_minthreshold))
sub_protist_order_tsc_minthreshold_b
colSums(sub_protist_order_tsc_minthreshold_b)


# Bar graphs of relative abundance

## Basic bar graphs based on phylum and order level for Algae
bar_algae_phylum <- microbiome::transform(sub_algae_phylum_tsc_minthreshold, "compositional")
plot_bar(bar_algae_phylum, fill="phylum", title = "The relative abundances of the algae at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_algae_order <- microbiome::transform(sub_algae_order_tsc_minthreshold, "compositional")
plot_bar(bar_algae_order, fill="order", title = "The relative abundances of the algae at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")

## Basic bar graphs based on phylum and order level for Protists
bar_protist_phylum <- microbiome::transform(sub_protist_phylum_tsc_minthreshold, "compositional")
plot_bar(bar_protist_phylum, fill="phylum", title = "The relative abundances of protists at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_protist_order <- microbiome::transform(sub_protist_order_tsc_minthreshold, "compositional")
plot_bar(bar_protist_order, fill="order", title = "The relative abundances of protists at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")


# Exploring the biodiversity : Alpha-diversity

## According to Mariadassou et al., (2016), this statistical analysis was performed.
## Followed by https://joey711.github.io/phyloseq/plot_richness-examples.html,
## For alpha diversity, it is advisable to prune OTUs that are not present in any of the samples.
## But that's all we can trim. It is tempting to trim noise, but there are many estimates of richness based on the singletons and double tons of the abundance data.
## To get a meaningful estimate, we need to leave it in the data set.

## For the algae phylum level
## To prune taxa
alph_alg_sset_phy_1 <- prune_taxa(taxa_sums(sub_algae_phylum) > 0, sub_algae_phylum)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_algae_phylum <- estimate_richness(alph_alg_sset_phy_1, measures = c("Observed", "Shannon"))
alpha.diversity_algae_phylum
## To visualize the estimated Observed Richness and Shannon Diversity
p_ap = plot_richness(alph_alg_sset_phy_1, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the algal phylum level")
p_ap + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness             
pielou_algae_phylum <- evenness(alph_alg_sset_phy_1, 'pielou')
pielou_algae_phylum
## To arrange the sample data in a data frame
dat_lm_algae_phylum <- data.frame(sample_data(alph_alg_sset_phy_1))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_algae_phylum <- cbind(alpha.diversity_algae_phylum, pielou_algae_phylum, dat_lm_algae_phylum)
alpha_lm_algae_phylum
alpha_lm_algae_phylum_1 <- rownames_to_column(alpha_lm_algae_phylum, var = "samples")
alpha_lm_algae_phylum_1
## To visualize the Pielous Evenness
ggplot(alpha_lm_algae_phylum_1, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the algal phylum level")
## To mutate the variables as factor
str(alpha_lm_algae_phylum_1)
alpha_lm_algae_phylum_11 <- alpha_lm_algae_phylum_1 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(alpha_lm_algae_phylum_11)

## For the algae order level
## To prune taxa
alph_alg_sset_ord_1 <- prune_taxa(taxa_sums(sub_algae_order) > 0, sub_algae_order)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_algae_order <- estimate_richness(alph_alg_sset_ord_1, measures = c("Observed", "Shannon"))
alpha.diversity_algae_order
## To visualize the estimated Observed Richness and Shannon Diversity
p_ao = plot_richness(alph_alg_sset_ord_1, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the algal order level")
p_ao + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness 
pielou_algae_order <- evenness(alph_alg_sset_ord_1, 'pielou')
pielou_algae_order
## To arrange the sample data in a data frame
dat_lm_algae_order <- data.frame(sample_data(alph_alg_sset_ord_1))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_algae_order <- cbind(alpha.diversity_algae_order, pielou_algae_order, dat_lm_algae_order)
alpha_lm_algae_order
alpha_lm_algae_order_1 <- rownames_to_column(alpha_lm_algae_order, var = "samples")
alpha_lm_algae_order_1
## To visualize the Pielous Evenness
ggplot(alpha_lm_algae_order_1, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the algal order level")
## To mutate the variables as factor
str(alpha_lm_algae_order_1)
alpha_lm_algae_order_11 <- alpha_lm_algae_order_1 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(alpha_lm_algae_order_11)

## For the protist phylum level
## To prune taxa
alph_prot_sset_phy_1 <- prune_taxa(taxa_sums(sub_protist_phylum) > 0, sub_protist_phylum)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_protist_phylum <- estimate_richness(alph_prot_sset_phy_1, measures = c("Observed", "Shannon"))
alpha.diversity_protist_phylum
## To visualize the estimated Observed Richness and Shannon Diversity
p_pp = plot_richness(alph_prot_sset_phy_1, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the protist phylum level")
p_pp + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness
pielou_protist_phylum <- evenness(alph_prot_sset_phy_1, 'pielou')
pielou_protist_phylum
## To arrange the sample data in a data frame
dat_lm_protist_phylum <- data.frame(sample_data(alph_prot_sset_phy_1))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_protist_phylum <- cbind(alpha.diversity_protist_phylum, pielou_protist_phylum, dat_lm_protist_phylum)
alpha_lm_protist_phylum
alpha_lm_protist_phylum_1 <- rownames_to_column(alpha_lm_protist_phylum, var = "samples")
alpha_lm_protist_phylum_1
## To visualize the Pielous Evenness
ggplot(alpha_lm_protist_phylum_1, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the protist phylum level")
## To mutate the variables as factor
str(alpha_lm_protist_phylum_1)
alpha_lm_protist_phylum_11 <- alpha_lm_protist_phylum_1 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(alpha_lm_protist_phylum_11)

## For the protist order level
## To prune taxa
alph_prot_sset_ord_1 <- prune_taxa(taxa_sums(sub_protist_order) > 0, sub_protist_order)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_protist_order <- estimate_richness(alph_prot_sset_ord_1, measures = c("Observed", "Shannon"))
alpha.diversity_protist_order
## To visualize the estimated Observed Richness and Shannon Diversity
p_po = plot_richness(alph_prot_sset_ord_1, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the protist order level")
p_po + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness
pielou_protist_order <- evenness(alph_prot_sset_ord_1, 'pielou')
pielou_protist_order
## To arrange the sample data in a data frame
dat_lm_protist_order <- data.frame(sample_data(alph_prot_sset_ord_1))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_protist_order <- cbind(alpha.diversity_protist_order, pielou_protist_order, dat_lm_protist_order)
alpha_lm_protist_order
alpha_lm_protist_order_1 <- rownames_to_column(alpha_lm_protist_order, var = "samples")
alpha_lm_protist_order_1
## To visualize the Pielous Evenness
ggplot(alpha_lm_protist_order_1, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the protist order level")
## To mutate the variables as factor
str(alpha_lm_protist_order_1)
alpha_lm_protist_order_11 <- alpha_lm_protist_order_1 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(alpha_lm_protist_order_11)


# Modeling
## According to Zuur et al., (2009), this statistical analysis was performed.

# For the algae phylum level

## For Observed
### To check the normality of the data
hist((alpha_lm_algae_phylum_11$Observed))
boxplot(Observed ~ treatment, data = alpha_lm_algae_phylum_11, main= "Observed richness for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11$Observed)
qqline(alpha_lm_algae_phylum_11$Observed)
shapiro.test(alpha_lm_algae_phylum_11$Observed)
### To fit a linear model with the gls function
ap_ob.gls <- gls(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(ap_ob.gls)
summary(ap_ob.gls)
### To fit a linear mixed model
lmm_algae_phylum_Observed <- lme(log(10+Observed) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11)
plot(lmm_algae_phylum_Observed)
summary(lmm_algae_phylum_Observed)
### To compare the models
anova(ap_ob.gls, lmm_algae_phylum_Observed)
### To fit the final model
lm_algae_phylum_Observed <- lm(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(lm_algae_phylum_Observed)
leveneTest(lm_algae_phylum_Observed)
summary(lm_algae_phylum_Observed)
anova(lm_algae_phylum_Observed)
eta_squared(lm_algae_phylum_Observed)
### Bar plots of ANOVA results for algal phylum level for Observed Richness
plot_ao_ap_ob <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+Observed), fill = treatment))+
    stat_summary(geom = "bar", fun = mean)+
    stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
    facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
    theme_classic()+
    theme(legend.position = "none")+
    labs( size = 20)+
    ylab("Observed Richness")+
    xlab("Treatments")+
    ggtitle("Bar plots of ANOVA results for algal phylum level for Observed Richness")+
    theme_classic()+
    theme(strip.text.x = element_text(size = 25),
          legend.position = "none",
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x =element_text(size = 20),
          axis.title.y =element_text(size = 20),
          strip.background = element_blank())
plot_ao_ap_ob
### Significant interaction plots of algal phylum level for Observed Richness
### Followed by https://psyteachr.github.io/introdataviz/multi-part-plots.html
interaction_ap_ob <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+Observed), 
                     shape = week,
                     group = week,
                     color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal phylum level for Observed Richness")+
  theme_classic()
interaction_ap_ob
   
## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_phylum_11$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_algae_phylum_11, main= "Shannon diversity for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11$Shannon)
qqline(alpha_lm_algae_phylum_11$Shannon)
shapiro.test(alpha_lm_algae_phylum_11$Shannon)
### To fit a linear model with the gls function
ap_sh.gls <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(ap_sh.gls)
summary(ap_sh.gls)
### To fit a linear mixed model
lmm_algae_phylum_Shannon <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11)
plot(lmm_algae_phylum_Shannon)
summary(lmm_algae_phylum_Shannon)
### To compare the models
anova(ap_sh.gls, lmm_algae_phylum_Shannon)
### To fit the final model
lm_algae_phylum_Shannon <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(lm_algae_phylum_Shannon)
leveneTest(lm_algae_phylum_Shannon)
summary(lm_algae_phylum_Shannon)
anova(lm_algae_phylum_Shannon)
eta_squared(lm_algae_phylum_Shannon)
### Bar plots of ANOVA results for algal phylum level for Shannon Diversity
plot_ao_ap_sh <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for algal phylum level for Shannon Diversity")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_ap_sh
### Significant interaction plots of algal phylum level for Shannon Diversity
interaction_ap_sh <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+Shannon), 
                                                             shape = week,
                                                             group = week,
                                                             color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal phylum level for Shannon Diversity")+
  theme_classic()
interaction_ap_sh

## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_phylum_11$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_algae_phylum_11, main= "Pielous evenness for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11$pielou)
qqline(alpha_lm_algae_phylum_11$pielou)
shapiro.test(alpha_lm_algae_phylum_11$pielou)
### To fit a linear model with the gls function
ap_pi.gls <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(ap_pi.gls)
summary(ap_pi.gls)
### To fit a linear mixed model
lmm_algae_phylum_pielou <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11)
plot(lmm_algae_phylum_pielou)
summary(lmm_algae_phylum_pielou)
### To compare the models
anova(ap_pi.gls, lmm_algae_phylum_pielou)
### To fit the final model
lm_algae_phylum_pielou <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11)
plot(lm_algae_phylum_pielou)
leveneTest(lm_algae_phylum_pielou)
summary(lm_algae_phylum_pielou)
anova(lm_algae_phylum_pielou)
eta_squared(lm_algae_phylum_pielou)
### Bar plots of ANOVA results for algal phylum level for Pielou's Evenness
plot_ao_ap_pi <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+pielou), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for algal phylum level for Pielou's Evenness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_ap_pi
### Significant interaction plots of algal phylum level for Pielou's Evenness
interaction_ap_pi <- ggplot(alpha_lm_algae_phylum_11, aes(x = treatment, y = log(10+pielou), 
                                                             shape = week,
                                                             group = week,
                                                             color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal phylum level for Pielou's Evenness")+
  theme_classic()
interaction_ap_pi


# For the algae order level

## For Observed
### To check the normality of the data
hist(alpha_lm_algae_order_11$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_algae_order_11, main= "Observed richness for the algae order level")
qqnorm(alpha_lm_algae_order_11$Observed)
qqline(alpha_lm_algae_order_11$Observed)
shapiro.test(alpha_lm_algae_order_11$Observed)
### To fit a linear model with the gls function
ao_ob.gls <- gls(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(ao_ob.gls)
summary(ao_ob.gls)
### To fit a linear mixed model
lmm_algae_order_Observed <- lme(log(10+Observed) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11)
plot(lmm_algae_order_Observed)
summary(lmm_algae_order_Observed)
### To compare the models
anova(ao_ob.gls, lmm_algae_order_Observed)
### To fit the final model
lm_algae_order_Observed <- lm(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(lm_algae_order_Observed)
leveneTest(lm_algae_order_Observed)
summary(lm_algae_order_Observed)
anova(lm_algae_order_Observed)
eta_squared(lm_algae_order_Observed)
### Bar plots of ANOVA results for algal order level for Observed Richness
plot_ao_ao_ob <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+Observed), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for algal order level for Observed Richness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_ao_ob
### Significant interaction plots of algal order level for Observed Richness
interaction_ao_ob <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+Observed), 
                                                          shape = week,
                                                          group = week,
                                                          color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal order level for Observed Richness")+
  theme_classic()
interaction_ao_ob

## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_order_11$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_algae_order_11, main= "Shannon diversity for the algae order level")
qqnorm(alpha_lm_algae_order_11$Shannon)
qqline(alpha_lm_algae_order_11$Shannon)
shapiro.test(alpha_lm_algae_order_11$Shannon)
### To fit a linear model with the gls function
ao_sh.gls <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(ao_sh.gls)
summary(ao_sh.gls)
### To fit a linear mixed model
lmm_algae_order_Shannon <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11)
plot(lmm_algae_order_Shannon)
summary(lmm_algae_order_Shannon)
### To compare the models
anova(ao_sh.gls, lmm_algae_order_Shannon)
### To fit the final model
lm_algae_order_Shannon <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(lm_algae_order_Shannon)
leveneTest(lm_algae_order_Shannon)
summary(lm_algae_order_Shannon)
anova(lm_algae_order_Shannon)
eta_squared(lm_algae_order_Shannon)
### Bar plots of ANOVA results for algal order level for Shannon Diversity
plot_ao_ao_sh <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for algal order level for Shannon Diversity")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_ao_sh
### Significant interaction plots of algal order level for Shannon Diversity
interaction_ao_sh <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+Shannon), 
                                                         shape = week,
                                                         group = week,
                                                         color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal order level for Shannon Diversity")+
  theme_classic()
interaction_ao_sh

## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_order_11$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_algae_order_11, main= "Pielous evenness for the algae order level")
qqnorm(alpha_lm_algae_order_11$pielou)
qqline(alpha_lm_algae_order_11$pielou)
shapiro.test(alpha_lm_algae_order_11$pielou)
### To fit a linear model with the gls function
ao_pi.gls <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(ao_pi.gls)
summary(ao_pi.gls)
### To fit a linear mixed model
lmm_algae_order_pielou <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11)
plot(lmm_algae_order_pielou)
summary(lmm_algae_order_pielou)
### To compare the models
anova(ao_pi.gls, lmm_algae_order_pielou)
### To fit the final model
lm_algae_order_pielou <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11)
plot(lm_algae_order_pielou)
leveneTest(lm_algae_order_pielou)
summary(lm_algae_order_pielou)
anova(lm_algae_order_pielou)
eta_squared(lm_algae_order_pielou)
### Bar plots of ANOVA results for algal order level for Pielou's Evenness
plot_ao_ao_pi <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+pielou), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for algal order level for Pielou's Evenness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_ao_pi
### Significant interaction plots of algal order level for Pielou's Evenness
interaction_ao_pi <- ggplot(alpha_lm_algae_order_11, aes(x = treatment, y = log(10+pielou), 
                                                         shape = week,
                                                         group = week,
                                                         color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of algal order level for Pielou's Evenness")+
  theme_classic()
interaction_ao_pi

# For the protist phylum level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_phylum_11$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_protist_phylum_11, main= "Observed richness for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11$Observed)
qqline(alpha_lm_protist_phylum_11$Observed)
shapiro.test(alpha_lm_protist_phylum_11$Observed)
### To fit a linear nmodel with the gls function
pp_ob.gls <- gls(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(pp_ob.gls)
summary(pp_ob.gls)
### To fit a linear mixed model
lmm_protist_phylum_Observed <- lme(log(10+Observed) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11)
plot(lmm_protist_phylum_Observed)
summary(lmm_protist_phylum_Observed)
### To compare the models
anova(pp_ob.gls, lmm_protist_phylum_Observed)
### To fit the final model
lm_protist_phylum_Observed <- lm(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(lm_protist_phylum_Observed)
leveneTest(lm_protist_phylum_Observed)
summary(lm_protist_phylum_Observed)
anova(lm_protist_phylum_Observed)
eta_squared(lm_protist_phylum_Observed)
### Bar plots of ANOVA results for protist phylum level for Observed Richness
plot_ao_pp_ob <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+Observed), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist phylum level for Observed Richness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_pp_ob
### Significant interaction plots of protist phylum level for Observed Richness
interaction_pp_ob <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+Observed), 
                                                          shape = week,
                                                          group = week,
                                                          color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist phylum level for Observed Richness")+
  theme_classic()
interaction_pp_ob

## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_phylum_11$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_protist_phylum_11, main= "Shannon diversity for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11$Shannon)
qqline(alpha_lm_protist_phylum_11$Shannon)
shapiro.test(alpha_lm_protist_phylum_11$Shannon)
### To fit a linear model with the gls function
pp_sh.gls <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(pp_sh.gls)
summary(pp_sh.gls)
### To fit a linear mixed model
lmm_protist_phylum_Shannon <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11)
plot(lmm_protist_phylum_Shannon)
summary(lmm_protist_phylum_Shannon)
### To compare the models
anova(pp_sh.gls, lmm_protist_phylum_Shannon)
### To fit the final model
lm_protist_phylum_Shannon <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(lm_protist_phylum_Shannon)
leveneTest(lm_protist_phylum_Shannon)
summary(lm_protist_phylum_Shannon)
anova(lm_protist_phylum_Shannon)
eta_squared(lm_protist_phylum_Shannon)
### Bar plots of ANOVA results for protist phylum level for Shannon Diversity
plot_ao_pp_sh <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist phylum level for Shannon Diversity")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_pp_sh
### Significant interaction plots of protist phylum level for Shannon Diversity
interaction_pp_sh <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+Shannon), 
                                                          shape = week,
                                                          group = week,
                                                          color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist phylum level for Shannon Diversity")+
  theme_classic()
interaction_pp_sh

## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_phylum_11$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_protist_phylum_11, main= "Pielous evenness for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11$pielou)
qqline(alpha_lm_protist_phylum_11$pielou)
shapiro.test(alpha_lm_protist_phylum_11$pielou)
### To fit a linear model with the gls function
pp_pi.gls <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(pp_pi.gls)
summary(pp_pi.gls)
### To fit a linear mixed model
lmm_protist_phylum_pielou <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11)
plot(lmm_protist_phylum_pielou)
summary(lmm_protist_phylum_pielou)
### To compare the models
anova(pp_pi.gls, lmm_protist_phylum_pielou)
### To fit the final model
lm_protist_phylum_pielou <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11)
plot(lm_protist_phylum_pielou)
leveneTest(lm_protist_phylum_pielou)
summary(lm_protist_phylum_pielou)
anova(lm_protist_phylum_pielou)
eta_squared(lm_protist_phylum_pielou)
### Bar plots of ANOVA results for protist phylum level for Pielou's Evenness
plot_ao_pp_pi <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+pielou), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist phylum level for Pielou's Evenness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_pp_pi
### Significant interaction plots of protist phylum level for Pielou's Evenness
interaction_pp_pi <- ggplot(alpha_lm_protist_phylum_11, aes(x = treatment, y = log(10+pielou), 
                                                          shape = week,
                                                          group = week,
                                                          color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist phylum level for Pielou's Evenness")+
  theme_classic()
interaction_pp_pi

# For the protist order level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_order_11$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_protist_order_11, main= "Observed richness for the protist order level")
qqnorm(alpha_lm_protist_order_11$Observed)
qqline(alpha_lm_protist_order_11$Observed)
shapiro.test(alpha_lm_protist_order_11$Observed)
### To fit a linear model with the gls function
po_ob.gls <- gls(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(po_ob.gls)
summary(po_ob.gls)
### To fit a linear mixed model
lmm_protist_order_Observed <- lme(log(10+Observed) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11)
plot(lmm_protist_order_Observed)
summary(lmm_protist_order_Observed)
### To compare the models
anova(po_ob.gls, lmm_protist_order_Observed)
### To fit the final model
lm_protist_order_Observed <- lm(log(10+Observed) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(lm_protist_order_Observed)
leveneTest(lm_protist_order_Observed)
summary(lm_protist_order_Observed)
anova(lm_protist_order_Observed)
eta_squared(lm_protist_order_Observed)
### Bar plots of ANOVA results for protist order level for Observed Richness
plot_ao_po_ob <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+Observed), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist order level for Observed Richness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_po_ob
### Significant interaction plots of protist order level for Observed Richness
interaction_po_ob <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+Observed), 
                                                            shape = week,
                                                            group = week,
                                                            color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Observed Richness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist order level for Observed Richness")+
  theme_classic()
interaction_po_ob

## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_order_11$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_protist_order_11, main= "Shannon diversity for the protist order level")
qqnorm(alpha_lm_protist_order_11$Shannon)
qqline(alpha_lm_protist_order_11$Shannon)
shapiro.test(alpha_lm_protist_order_11$Shannon)
### To fit a linear model with the gls function
po_sh.gls <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(po_sh.gls)
summary(po_sh.gls)
### To fit a linear mixed model
lmm_protist_order_Shannon <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11)
plot(lmm_protist_order_Shannon)
summary(lmm_protist_order_Shannon)
### To compare the models
anova(po_sh.gls, lmm_protist_order_Shannon)
### To fit the final model
lm_protist_order_Shannon <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(lm_protist_order_Shannon)
leveneTest(lm_protist_order_Shannon)
summary(lm_protist_order_Shannon)
anova(lm_protist_order_Shannon)
eta_squared(lm_protist_order_Shannon)
### Bar plots of ANOVA results for protist order level for Shannon Diversity
plot_ao_po_sh <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist order level for Shannon Diversity")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_po_sh
### Significant interaction plots of protist order level for Shannon Diversity
interaction_po_sh <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+Shannon), 
                                                           shape = week,
                                                           group = week,
                                                           color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Shannon Diversity")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist order level for Shannon Diversity")+
  theme_classic()
interaction_po_sh

## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_order_11$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_protist_order_11, main= "Pielous evenness for the protist order level")
qqnorm(alpha_lm_protist_order_11$pielou)
qqline(alpha_lm_protist_order_11$pielou)
shapiro.test(alpha_lm_protist_order_11$pielou)
### To fit a linear model with the gls function
po_pi.gls <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(po_pi.gls)
summary(po_pi.gls)
### To fit a linear mixed model
lmm_protist_order_pielou <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11)
plot(lmm_protist_order_pielou)
summary(lmm_protist_order_pielou)
### To compare the models
anova(po_pi.gls, lmm_protist_order_pielou)
### To fit the final model
lm_protist_order_pielou <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11)
plot(lm_protist_order_pielou)
leveneTest(lm_protist_order_pielou)
summary(lm_protist_order_pielou)
anova(lm_protist_order_pielou)
eta_squared(lm_protist_order_pielou)
### Bar plots of ANOVA results for protist order level for Pielou's Evenness
plot_ao_po_pi <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+pielou), fill = treatment))+
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=.1,position=position_dodge(.9))+
  facet_wrap(~week, labeller = labeller(week = c("2" = "After two weeks", "3" = "After three weeks")))+
  theme_classic()+
  theme(legend.position = "none")+
  labs( size = 20)+
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Bar plots of ANOVA results for protist order level for Pielou's Evenness")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x =element_text(size = 20),
        axis.title.y =element_text(size = 20),
        strip.background = element_blank())
plot_ao_po_pi
### Significant interaction plots of protist order level for Pielou's Evenness
interaction_po_pi <- ggplot(alpha_lm_protist_order_11, aes(x = treatment, y = log(10+pielou), 
                                                           shape = week,
                                                           group = week,
                                                           color = week)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  scale_color_manual(values = c("blue", "darkorange")) +
  ylab("Pielou's Evenness")+
  xlab("Treatments")+
  ggtitle("Significant interaction plots of protist order level for Pielou's Evenness")+
  theme_classic()
interaction_po_pi



# BETA-DIVERSITY INDICES

# According to Ollberding (2019), the PERMANOVA and PERMDISP analysis were performed.

# Permutational multivariate analysis of variance (PERMANOVA)
# Followed by https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova

## For the algae phylum level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_algae_phylum <- phyloseq::distance(sub_algae_phylum_tsc_minthreshold, method = "bray")
## To make a data frame from the sample_data
perm_algae_phylum_sampledf <- data.frame(sample_data(sub_algae_phylum_tsc_minthreshold))
## To make the varisables as factors
str(perm_algae_phylum_sampledf)
perm_algae_phylum_sampledf_1 <- perm_algae_phylum_sampledf %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(perm_algae_phylum_sampledf_1)
## To perform the adonis test
adonis2(bray_algae_phylum ~ nutrient*sediment*flow*time, strata = perm_algae_phylum_sampledf_1$week, data = perm_algae_phylum_sampledf_1)

## For the algae order level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_algae_order <- phyloseq::distance(sub_algae_order_tsc_minthreshold, method = "bray")
## To make a data frame from the sample_data
perm_algae_order_sampledf <- data.frame(sample_data(sub_algae_order_tsc_minthreshold))
## To make the variables as factors
str(perm_algae_order_sampledf)
perm_algae_order_sampledf_1 <- perm_algae_order_sampledf %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(perm_algae_order_sampledf_1)
## To perform the adonis test
adonis2(bray_algae_order ~ nutrient*sediment*flow*time, strata = perm_algae_order_sampledf_1$week, data = perm_algae_order_sampledf_1)

## For the protist phylum level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_protist_phylum <- phyloseq::distance(sub_protist_phylum_tsc_minthreshold, method = "bray")
## To make a data frame from the sample_data
perm_protist_phylum_sampledf <- data.frame(sample_data(sub_protist_phylum_tsc_minthreshold))
## To make the variables as factors
str(perm_protist_phylum_sampledf)
perm_protist_phylum_sampledf_1 <- perm_protist_phylum_sampledf %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(perm_protist_phylum_sampledf_1)
## To perform the adonis test
adonis2(bray_protist_phylum ~ nutrient*sediment*flow*time, strata = perm_protist_phylum_sampledf_1$week, data = perm_protist_phylum_sampledf_1)

## For the protist order level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_protist_order <- phyloseq::distance(sub_protist_order_tsc_minthreshold, method = "bray")
## To make a data frame from the sample_data
perm_protist_order_sampledf <- data.frame(sample_data(sub_protist_order_tsc_minthreshold))
## To make the variables as factors
str(perm_protist_order_sampledf)
perm_protist_order_sampledf_1 <- perm_protist_order_sampledf %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(perm_protist_order_sampledf_1)
## To perform the adonis test
adonis2(bray_protist_order ~ nutrient*sediment*flow*time, strata = perm_protist_order_sampledf_1$week, data = perm_protist_order_sampledf_1)


# Permutational analysis of multivariate dispersion (PERMDISP)

## For the algae phylum level
dispersion_algae_phylum <- vegan::betadisper(bray_algae_phylum, phyloseq::sample_data(sub_algae_phylum_tsc_minthreshold)$treatment)
dispersion_algae_phylum
plot(dispersion_algae_phylum, main = "For the algae phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_algae_phylum, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_phylum)

## For the algae order level
dispersion_algae_order <- vegan::betadisper(bray_algae_order, phyloseq::sample_data(sub_algae_order_tsc_minthreshold)$treatment)
dispersion_algae_order
plot(dispersion_algae_order, main = "For the algae order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_algae_order, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_order)

## For the protist phylum level
dispersion_protist_phylum <- vegan::betadisper(bray_protist_phylum, phyloseq::sample_data(sub_protist_phylum_tsc_minthreshold)$treatment)
dispersion_protist_phylum
plot(dispersion_protist_phylum, main = "For the protist phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_protist_phylum, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_phylum)

## For the protist order level
dispersion_protist_order <- vegan::betadisper(bray_protist_order, phyloseq::sample_data(sub_protist_order_tsc_minthreshold)$treatment)
dispersion_protist_order
plot(dispersion_protist_order, main = "For the protist order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_protist_order, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_order)


# The PCoA analysis (Principal Coordinate Analysis)
# Followed by, https://joey711.github.io/phyloseq/plot_ordination-examples.html, this statistical analysis was performed.

pcoa_dist = "bray"
pcoa_ord_meths = c("PCoA")

## For the algae phylum level
ap_plist = lapply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, sub_algae_phylum_tsc_minthreshold, pcoa_dist)

names(ap_plist) <- pcoa_ord_meths

ap_pdataframe = plyr::ldply(ap_plist, function(x){
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
  ggtitle("PCoA1(22.78%) and PCoA2(18.34%)")
ap_p

## For the algae order level
ao_plist = lapply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, sub_algae_order_tsc_minthreshold, pcoa_dist)

names(ao_plist) <- pcoa_ord_meths

ao_pdataframe = plyr::ldply(ao_plist, function(x){
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
  ggtitle("PCoA1(20.88%) and PCoA2(17.99%)")
ao_p

## For the protist phylum level
pp_plist = lapply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, sub_protist_phylum_tsc_minthreshold, pcoa_dist)

names(pp_plist) <- pcoa_ord_meths

pp_pdataframe = plyr::ldply(pp_plist, function(x){
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
  ggtitle("PCoA1(28.36%) and PCoA2(19.51%)")
pp_p

## For the protist order level
po_plist = lapply(as.list(pcoa_ord_meths), function(i, physeq, pcoa_dist){
  ordi = ordinate(physeq, method=i, distance=pcoa_dist)
  plot_ordination(physeq, ordi, "samples", color="treatment")
}, sub_protist_order_tsc_minthreshold, pcoa_dist)

names(po_plist) <- pcoa_ord_meths

po_pdataframe = plyr::ldply(po_plist, function(x){
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
  ggtitle("PCoA1(25.76%) and PCoA2(18.13%)")
po_p


# Followed by https://github.com/joey711/phyloseq/issues/847, the most abundant taxa was calculated. 
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
prune_algae_phylum <- prune_taxa(taxa_sums(sub_algae_phylum_tsc_minthreshold) > 0.01, sub_algae_phylum_tsc_minthreshold)
abundant_algae_phylum <- most_abundant_taxa(prune_algae_phylum,"phylum")
abundant_algae_phylum

prune_algae_order <- prune_taxa(taxa_sums(sub_algae_order_tsc_minthreshold) > 0.01, sub_algae_order_tsc_minthreshold)
abundant_algae_order <- most_abundant_taxa(prune_algae_order,"order")
abundant_algae_order

prune_protist_phylum <- prune_taxa(taxa_sums(sub_protist_phylum_tsc_minthreshold) > 0.01, sub_protist_phylum_tsc_minthreshold)
abundant_protist_phylum <- most_abundant_taxa(prune_protist_phylum,"phylum")
abundant_protist_phylum

prune_protist_order <- prune_taxa(taxa_sums(sub_protist_order_tsc_minthreshold) > 0.01, sub_protist_order_tsc_minthreshold)
abundant_protist_order <- most_abundant_taxa(prune_protist_order,"order")
abundant_protist_order

# Multivariate Analysis Of Variance (MANOVA)
## According to Kassambara (2017); Ben-Shachar et al., (2020), this statistical analysis was performed.

## For the most abundant taxa of the algal phylum level
## To make a sample data frame
manova_df_algae_phylum = as(sample_data(prune_algae_phylum), "data.frame")
manova_df_algae_phylum
## To combine the sample data frame and abundant taxa data frame
manova_df_algae_phylum_1 = cbind(manova_df_algae_phylum, abundant_algae_phylum)
manova_df_algae_phylum_1

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Bacillariophyta, Chlorophyta, Streptophyta, Phragmoplastophyta and unclassified 
shapiro.test(manova_df_algae_phylum_1$Bacillariophyta)
Bacillariophyta_1 <- log10(manova_df_algae_phylum_1$Bacillariophyta)
shapiro.test(Bacillariophyta_1)

shapiro.test(manova_df_algae_phylum_1$Chlorophyta)

shapiro.test(manova_df_algae_phylum_1$Streptophyta)
Streptophyta_1 <- sqrt(manova_df_algae_phylum_1$Streptophyta)
shapiro.test(Streptophyta_1)

shapiro.test(manova_df_algae_phylum_1$Phragmoplastophyta)
Phragmoplastophyta_1 <- sqrt(manova_df_algae_phylum_1$Phragmoplastophyta)
shapiro.test(Phragmoplastophyta_1)

shapiro.test(manova_df_algae_phylum_1$unclassified)
unclassified_1 <- sqrt(manova_df_algae_phylum_1$unclassified)
shapiro.test(unclassified_1)
## To combine the transformed data
manova_df_algae_phylum_11 <- cbind(manova_df_algae_phylum_1, Bacillariophyta_1, Streptophyta_1, Phragmoplastophyta_1, unclassified_1)
## To make the variables as factors
str(manova_df_algae_phylum_11)
manova_df_algae_phylum_111 <- manova_df_algae_phylum_11 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(manova_df_algae_phylum_111)
## To do the leveneTest
leveneTest(Bacillariophyta_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
leveneTest(Chlorophyta ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
leveneTest(Streptophyta_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
leveneTest(Phragmoplastophyta_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
leveneTest(unclassified_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
## To do the MANOVA test
manova_abundant_algae_phylum <- manova(cbind(Bacillariophyta_1, Chlorophyta, Streptophyta_1, Phragmoplastophyta_1, unclassified_1) ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111)
manova_abundant_algae_phylum
summary.aov(manova_abundant_algae_phylum)
eta_squared(aov(manova_abundant_algae_phylum))

## For the most abundant taxa of the algal order level
## To make a sample data frame
manova_df_algae_order = as(sample_data(prune_algae_order), "data.frame")
manova_df_algae_order
## To combine the sample data frame and abundant taxa data frame
manova_df_algae_order_1 = cbind(manova_df_algae_order, abundant_algae_order)
manova_df_algae_order_1

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Licmophorales, unclassified_Chlorophyta, Zygnematales, Cocconeidales, unclassified, unclassified_Chlorophyceae, Chlamydomonadales, Melosirales, Sphaeropleales, Cercomonadida, unclassified_Trebouxiophyceae and Chlorellales 
shapiro.test(manova_df_algae_order_1$Licmophorales)
Licmophorales_1 <- sqrt(manova_df_algae_order_1$Licmophorales)
shapiro.test(Licmophorales_1)

shapiro.test(manova_df_algae_order_1$unclassified_Chlorophyta)

shapiro.test(manova_df_algae_order_1$Zygnematales)
Zygnematales_1 <- sqrt(manova_df_algae_order_1$Zygnematales)
shapiro.test(Zygnematales_1)

shapiro.test(manova_df_algae_order_1$Cocconeidales)
Cocconeidales_1 <- sqrt(manova_df_algae_order_1$Cocconeidales)
shapiro.test(Cocconeidales_1)

shapiro.test(manova_df_algae_order_1$unclassified)
unclassified_2 <- sqrt(manova_df_algae_order_1$unclassified)
shapiro.test(unclassified_2)

shapiro.test(manova_df_algae_order_1$unclassified_Chlorophyceae)
unclassified_Chlorophyceae_1 <- sqrt(manova_df_algae_order_1$unclassified_Chlorophyceae)
shapiro.test(unclassified_Chlorophyceae_1)

shapiro.test(manova_df_algae_order_1$Chlamydomonadales)
Chlamydomonadales_1 <- sqrt(manova_df_algae_order_1$Chlamydomonadales)
shapiro.test(Chlamydomonadales_1)

shapiro.test(manova_df_algae_order_1$Melosirales)
Melosirales_1 <- sqrt(manova_df_algae_order_1$Melosirales)
shapiro.test(Melosirales_1)

shapiro.test(manova_df_algae_order_1$Sphaeropleales)

shapiro.test(manova_df_algae_order_1$Cercomonadida)
Cercomonadida_1 <- sqrt(manova_df_algae_order_1$Cercomonadida)
shapiro.test(Cercomonadida_1)

shapiro.test(manova_df_algae_order_1$unclassified_Trebouxiophyceae)

shapiro.test(manova_df_algae_order_1$Chlorellales)
Chlorellales_1 <- sqrt(manova_df_algae_order_1$Chlorellales)
shapiro.test(Chlorellales_1)
## To combine the transformed data
manova_df_algae_order_11 <- cbind(manova_df_algae_order_1, Licmophorales_1, Zygnematales_1, Cocconeidales_1, unclassified_2, unclassified_Chlorophyceae_1, Chlamydomonadales_1, Melosirales_1, Cercomonadida_1, Chlorellales_1)
## To make the variables as factors
str(manova_df_algae_order_11)
manova_df_algae_order_111 <- manova_df_algae_order_11 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(manova_df_algae_order_111)
## To do the leveneTest
leveneTest(Licmophorales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(unclassified_Chlorophyta ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Zygnematales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Cocconeidales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(unclassified_2 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(unclassified_Chlorophyceae_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Chlamydomonadales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Melosirales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Sphaeropleales ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Cercomonadida_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(unclassified_Trebouxiophyceae ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
leveneTest(Chlorellales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
## To do the MANOVA test
manova_abundant_algae_order <- manova(cbind(Licmophorales_1, unclassified_Chlorophyta, Zygnematales_1, Cocconeidales_1, unclassified_2, unclassified_Chlorophyceae_1, Chlamydomonadales_1, Melosirales_1, Sphaeropleales, Cercomonadida_1, unclassified_Trebouxiophyceae, Chlorellales_1) ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111)
manova_abundant_algae_order
summary.aov(manova_abundant_algae_order)
eta_squared(aov(manova_abundant_algae_order))

## For the most abundant taxa of the protist phylum level
## To make a sample data frame
manova_df_protist_phylum = as(sample_data(prune_protist_phylum), "data.frame")
manova_df_protist_phylum
## To combine the sample data frame and abundant taxa data frame
manova_df_protist_phylum_1 = cbind(manova_df_protist_phylum, abundant_protist_phylum)
manova_df_protist_phylum_1

## Data normality check and transformation of data to the normality for the most abundant taxa
### For unclassified_Eukaryota, unclassified and Apicomplexa 
shapiro.test(manova_df_protist_phylum_1$unclassified_Eukaryota)

shapiro.test(manova_df_protist_phylum_1$unclassified)
unclassified_3 <- sqrt(manova_df_protist_phylum_1$unclassified)
shapiro.test(unclassified_3)

shapiro.test(manova_df_protist_phylum_1$Apicomplexa)
Apicomplexa_1 <- sqrt(manova_df_protist_phylum_1$Apicomplexa)
shapiro.test(Apicomplexa_1)
## To combine the transformed data
manova_df_protist_phylum_11 <- cbind(manova_df_protist_phylum_1, unclassified_3, Apicomplexa_1)
## To make the variables as factors
str(manova_df_protist_phylum_11)
manova_df_protist_phylum_111 <- manova_df_protist_phylum_11 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(manova_df_protist_phylum_111)
## To do the leveneTest
leveneTest(unclassified_Eukaryota ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111)
leveneTest(unclassified_3 ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111)
leveneTest(Apicomplexa_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111)
## To do the MANOVA test
manova_abundant_protist_phylum <- manova(cbind(unclassified_Eukaryota, unclassified_3, Apicomplexa_1) ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111)
manova_abundant_protist_phylum
summary.aov(manova_abundant_protist_phylum)
eta_squared(aov(manova_abundant_protist_phylum))

## For the most abundant taxa of the protist order level
## To make a sample data frame
manova_df_protist_order = as(sample_data(prune_protist_order), "data.frame")
manova_df_protist_order
## To combine the sample data frame and abundant taxa data frame
manova_df_protist_order_1 = cbind(manova_df_protist_order, abundant_protist_order)
manova_df_protist_order_1

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Euglyphida, Eucoccidiorida, Eugregarinorida, Urostylida,  unclassified, unclassified_Euamoebida, Sporadotrichida, unclassified_unidentified, Glissomonadida, unclassified_Leptomyxida, Arcellinida and Pleurostomatida
shapiro.test(manova_df_protist_order_1$Euglyphida)
Euglyphida_1 <- sqrt(manova_df_protist_order_1$Euglyphida)
shapiro.test(Euglyphida_1)

shapiro.test(manova_df_protist_order_1$Eucoccidiorida)
Eucoccidiorida_1 <- sqrt(manova_df_protist_order_1$Eucoccidiorida)
shapiro.test(Eucoccidiorida_1)

shapiro.test(manova_df_protist_order_1$Eugregarinorida)
Eugregarinorida_1 <- sqrt(manova_df_protist_order_1$Eugregarinorida)
shapiro.test(Eugregarinorida_1)

shapiro.test(manova_df_protist_order_1$Urostylida)
Urostylida_1 <- sqrt(manova_df_protist_order_1$Urostylida)
shapiro.test(Urostylida_1)

shapiro.test(manova_df_protist_order_1$unclassified)
unclassified_4 <- sqrt(manova_df_protist_order_1$unclassified)
shapiro.test(unclassified_4)

shapiro.test(manova_df_protist_order_1$unclassified_Euamoebida)
unclassified_Euamoebida_1 <- sqrt(manova_df_protist_order_1$unclassified_Euamoebida)
shapiro.test(unclassified_Euamoebida_1)

shapiro.test(manova_df_protist_order_1$Sporadotrichida)
Sporadotrichida_1 <- sqrt(manova_df_protist_order_1$Sporadotrichida)
shapiro.test(Sporadotrichida_1)

shapiro.test(manova_df_protist_order_1$unclassified_unidentified)
unclassified_unidentified_1 <- sqrt(manova_df_protist_order_1$unclassified_unidentified)
shapiro.test(unclassified_unidentified_1)

shapiro.test(manova_df_protist_order_1$Glissomonadida)
Glissomonadida_1 <- sqrt(manova_df_protist_order_1$Glissomonadida)
shapiro.test(Glissomonadida_1)

shapiro.test(manova_df_protist_order_1$unclassified_Leptomyxida)
unclassified_Leptomyxida_1 <- sqrt(manova_df_protist_order_1$unclassified_Leptomyxida)
shapiro.test(unclassified_Leptomyxida_1)

shapiro.test(manova_df_protist_order_1$Arcellinida)
Arcellinida_1 <- sqrt(manova_df_protist_order_1$Arcellinida)
shapiro.test(Arcellinida_1)

shapiro.test(manova_df_protist_order_1$Pleurostomatida)
Pleurostomatida_1 <- sqrt(manova_df_protist_order_1$Pleurostomatida)
shapiro.test(Pleurostomatida_1)
## To combine the transformed data
manova_df_protist_order_11 <- cbind(manova_df_protist_order_1, Euglyphida_1, Eucoccidiorida_1, Eugregarinorida_1, Urostylida_1,  unclassified_4, unclassified_Euamoebida_1, Sporadotrichida_1, unclassified_unidentified_1, Glissomonadida_1, unclassified_Leptomyxida_1, Arcellinida_1, Pleurostomatida_1)
## To make the variables as factors
str(manova_df_protist_order_11)
manova_df_protist_order_111 <- manova_df_protist_order_11 %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment, treatment_week), as.factor))
str(manova_df_protist_order_111)
## To do the leveneTest
leveneTest(Euglyphida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Eucoccidiorida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Eugregarinorida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Urostylida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(unclassified_4 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(unclassified_Euamoebida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Sporadotrichida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(unclassified_unidentified_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Glissomonadida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(unclassified_Leptomyxida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Arcellinida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
leveneTest(Pleurostomatida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
## To do the MANOVA test
manova_abundant_protist_order <- manova(cbind(Euglyphida_1, Eucoccidiorida_1, Eugregarinorida_1, Urostylida_1,  unclassified_4, unclassified_Euamoebida_1, Sporadotrichida_1, unclassified_unidentified_1, Glissomonadida_1, unclassified_Leptomyxida_1, Arcellinida_1, Pleurostomatida_1) ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111)
manova_abundant_protist_order
summary.aov(manova_abundant_protist_order)
eta_squared(aov(manova_abundant_protist_order))




# List of References

## Vaulot, D. (2021, Feb 15). Phyloseq tutorial. https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

## McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." _PLoS ONE_, *8*(4), e61217. <http://dx.plos.org/10.1371/journal.pone.0061217>.

## Mariadassou, M., Bernard, M., Pascal, G., Cauquil, L., & Chaillou, S. (2016). Analysis of community composition data using phyloseq. Montpellier Décembre 2016. https://genoweb.toulouse.inra.fr/~formation/15_FROGS/8-February2017/FROGS_phyloseq_02_2017.pdf

## Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M. (2009). Mixed effects models and extensions in ecology with R (Vol. 574). New York: Springer.

## Ollberding, N.J. (2019, Jul 28). Introduction to the Statistical Analysis of Microbiome Data in R. Wowchemy. https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

## Kassambara, A. (2017). MANOVA Test in R: Multivariate Analysis of Variance. Statistical tools for high-throughput data analysis. STHDA. http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance#infos

## Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815  

