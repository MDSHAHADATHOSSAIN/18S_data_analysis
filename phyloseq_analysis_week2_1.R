# Analysis of eukaryotic biodiversity data using the phyloseq package in R (week 2)

## According to Vaulot (2021), this statistical analysis was performed.

## To install necessary packages
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


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("curatedMetagenomicData")

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

## To load the installed packages with libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")
library("microbiome")
library("reshape2")
library("ape")
library("gridExtra")
library("plotly")
library("vegan")
library("dendextend")
library("tidyr")
library("rms")
library("effectsize")
library("curatedMetagenomicData")
suppressPackageStartupMessages(library(curatedMetagenomicData))
library("lme4")
library("picante")
library("cowplot")

## To read the data and create phyloseq objects
otu_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "Samples")

## Phyloseq objects must have row.names

## To define the row names from the otu column
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 

## To idem for the two other matrixes 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

## For sample data
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

## Transformation to matrices otu and tax tables (the sample table can be left as a data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

## Transforming to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

phyloseq_obj <- phyloseq(OTU, TAX, samples)
phyloseq_obj

## To normalize number of reads in each sample using median sequencing depth and filter the values that are above 1e-5
phy_norm  = transform_sample_counts(phyloseq_obj, function(x) x / sum(x) ) ## normalization
phy_filtered = filter_taxa(phy_norm, function(x) mean(x) > 1e-5, TRUE) ## filtering

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

## To visualize the data

sample_names(phy_filtered)
rank_names(phy_filtered)
sample_variables(phy_filtered)

## To keep only taxa of interests (Algae) according to phylum and genus level
sub_algae_phylum <- subset_taxa(phy_filtered, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Phragmoplastophyta", "Rigifilida", "Schizoplasmodiida", "Streptophyta", "Xanthophyceae"))
sub_algae_genus <- subset_taxa(phy_filtered, genus %in% c("Messastrum", "Melosira", "Poteriospumella", "Pseudellipsoidion", "Symbiodinium", "Nasturtium", "unclassified_Rigifilida", "Phalansterium", "Reynoutria", "Pleurochloris"))

## To keep only taxa of interests (Protists) according to phylum and genus level
sub_protist_phylum <- subset_taxa(phy_filtered, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Rotifera", "Tubulinea", "Bicosoecida", "Discosea", "Gastrotricha", "Labyrinthulomycetes", "Peronosporomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida"))
sub_protist_genus <- subset_taxa(phy_filtered, genus %in% c("Gymnophrys", "unclassified_Craspedida", "Ichthyophthirius", "Cephalodella", "Flabellula", "Bicosoeca", "unclassified_Discosea", "Lepidodermella", "Sorodiplophrys", "Phytophthora", "unclassified_Apusomonadidae", "unclassified_CV1-B1-93", "Protosporangium", "Cryptosporidium", "unclassified_LEMD267", "unclassified_Protosteliida"))

## Bar graphs of relative abundance

## Basic bar graphs based on phylum and genus level for Algae
bar_algae_phylum <- microbiome::transform(sub_algae_phylum, "compositional")
plot_bar(bar_algae_phylum, fill="phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)

bar_algae_genus <- microbiome::transform(sub_algae_genus, "compositional")
plot_bar(bar_algae_genus, fill="genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)

## Basic bar graphs based on phylum and genus level for Protists
bar_protist_phylum <- microbiome::transform(sub_protist_phylum, "compositional")
plot_bar(bar_protist_phylum, fill="phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)

bar_protist_genus <- microbiome::transform(sub_protist_genus, "compositional")
plot_bar(bar_protist_genus, fill="genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)


# According to Ollberding (2019), this statistical analysis was performed.
## Here we will create some boxplots by treatment type for algae on the selected taxa of interest and facet them by phylum based on the raw data.

## Agglomeration at phylum level and renaming
ps_sub_algae_phylum <- phyloseq::tax_glom(sub_algae_phylum, "phylum")
phyloseq::taxa_names(ps_sub_algae_phylum) <- phyloseq::tax_table(ps_sub_algae_phylum)[, "phylum"]
phyloseq::otu_table(ps_sub_algae_phylum)[1:5, 1:5]

## To melt and plot
phyloseq::psmelt(ps_sub_algae_phylum) %>%
  ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

## Here we will create some boxplots by treatment type for algae on the selected taxa of interest and facet them by genus based on the raw data.

## Agglomeration at genus level and renaming
ps_sub_algae_genus <- phyloseq::tax_glom(sub_algae_genus, "genus")
phyloseq::taxa_names(ps_sub_algae_genus) <- phyloseq::tax_table(ps_sub_algae_genus)[, "genus"]
phyloseq::otu_table(ps_sub_algae_genus)[1:5, 1:5]

## To melt and plot
phyloseq::psmelt(ps_sub_algae_genus) %>%
  ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

## Here we will create some boxplots by treatment type for protists on the selected taxa of interest and facet them by phylum based on the raw data.

## Agglomeration at phylum level and renaming
ps_sub_protist_phylum <- phyloseq::tax_glom(sub_protist_phylum, "phylum")
phyloseq::taxa_names(ps_sub_protist_phylum) <- phyloseq::tax_table(ps_sub_protist_phylum)[, "phylum"]
phyloseq::otu_table(ps_sub_protist_phylum)[1:5, 1:5]

## To melt and plot
phyloseq::psmelt(ps_sub_protist_phylum) %>%
  ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

## Here we will create some boxplots by treatment type for protists on the selected taxa of interest and facet them by genus based on the raw data.

## Agglomeration at genus level and renaming
ps_sub_protist_genus <- phyloseq::tax_glom(sub_protist_genus, "genus")
phyloseq::taxa_names(ps_sub_protist_genus) <- phyloseq::tax_table(ps_sub_protist_genus)[, "genus"]
phyloseq::otu_table(ps_sub_protist_genus)[1:5, 1:5]

## To melt and plot
phyloseq::psmelt(ps_sub_protist_genus) %>%
  ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ genus, scales = "free")


# According to Ollberding (2019), this statistical analysis was performed.
# Hierarchical clustering
## Based on taxonomic (dis)similarity, we will examine here how the samples were clustered.
## Depending on the Bray-Curtis dissimilarity, a hierarchical clustering of the samples is performed.
## The samples that have no taxa in common and have similar composition contain zero Bray-Curtis dissimilarity.
## To calculate the Bray-Curtis dissimilarity for the samples, the vegan package of community ecology is used here.
## Wards clustering is applied, where Wards clustering executes different clustered pairs based on the high level in each iteration that minimizes the increase in total variance.
## How the samples from the different treatments were clustered is assessed here by coloring the sample names.

## To extract the OTU table and calculate the BC
phy_filtered_rel_otu <- data.frame(phyloseq::otu_table(phy_filtered))
phy_filtered_rel_otu <- t(phy_filtered_rel_otu)
bc_distance <- vegan::vegdist(phy_filtered_rel_otu, method = "bray")
as.matrix(bc_distance)[1:5, 1:5]

## To save as a dendrogram
ward_dend <- as.dendrogram(hclust(bc_distance, method = "ward.D2"))

## To provide the color codes
meta_code <- data.frame(phyloseq::sample_data(phy_filtered))
color_code <- c(C = "red", N = "blue", S = "yellow", F = "green", NS = "black", NF = "purple", SF = "orange", NFS = "pink")
labels_colors(ward_dend) <- color_code[meta_code$treatment][order.dendrogram(ward_dend)]

## Plot
plot(ward_dend)
## The Bray-Curtis dissimilarity range for the sample is highly variable.
## Different samples have different compositions.
## Based on the treatment of the samples, there are some clusters in the dendogram near the tips, but there are obviously no high-level clusters.


# According to Mariadassou et al., (2016), this statistical analysis was performed.
## Exploring the biodiversity : Alpha-diversity
## The overall richness are plotted with plot_richness function
## We are interested in alpha diversity
## A good idea is to prune OTUs that are not present in any of the samples
## But we can not trim more than that
## Many richness estimates are modeled on singletons and doubletons in the abundance data
## We need to leave them in the dataset if we want a meaningful estimate

phy_alpha_prune <- prune_taxa(taxa_sums(phyloseq_obj) > 0 , phyloseq_obj)

## To try it on the trimmed phyloseq object
p1 <- plot_richness(phy_alpha_prune)
p1

## To custom it on treatment, first color on treatment, then select measures as "Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", and plot as boxplot
p2 <- plot_richness(phy_alpha_prune, color = "treatment", x = "treatment", measures = c("Observed", "Chao1", "Shannon", "Simpson", 
                                                                       "InvSimpson"))
p2

## To plot as boxplot
p3 <- p2 + geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(p3)

## The numerical values of the alpha diversities are provided by estimate_richness (applied internally by plot_richness)
alpha.diversity <- estimate_richness(phy_alpha_prune, measures = c("Observed", 
                                                            "Chao1", "Shannon", "Simpson", 
                                                            "InvSimpson"))
alpha.diversity
head(alpha.diversity)
write.table(alpha.diversity, "myfile.txt")

## A rapid ANOVA test: tests whether the observed richness is significantly different as a function of treatment
data <- cbind(sample_data(phy_alpha_prune), alpha.diversity)
phy_alpha_prune_obs.anova <- aov(Observed ~ treatment, data)
summary(phy_alpha_prune_obs.anova)
## The type of treatment has a significant impact on the richness

## A rapid ANOVA test: tests whether the Shannon richness is significantly different as a function of treatment
phy_alpha_prune_sha.anova <- aov(Shannon ~ treatment, data)
summary(phy_alpha_prune_sha.anova)
## The type of treatment has a significant impact on the richness

## A rapid ANOVA test: tests whether the Chao1 richness is significantly different as a function of treatment
phy_alpha_prune_chao1.anova <- aov(Chao1 ~ treatment, data)
summary(phy_alpha_prune_chao1.anova)
## The type of treatment has a significant impact on the richness

## A rapid ANOVA test: tests whether the Simpson richness is significantly different as a function of treatment
phy_alpha_prune_sim.anova <- aov(Simpson ~ treatment, data)
summary(phy_alpha_prune_sim.anova)
## The type of treatment has a significant impact on the richness

## A rapid ANOVA test: tests whether the InvSimpson richness is significantly different as a function of treatment
phy_alpha_prune_invsim.anova <- aov(InvSimpson ~ treatment, data)
summary(phy_alpha_prune_invsim.anova)
## The type of treatment has a significant impact on the richness


# Some ordination methods
## According to Vaulot (2021), this statistical analysis was performed.
## We will perform a multivariate analysis based on the Bray-Curtis distance and the NMDS ordination

## For the phylum level of algae
multi_ord_alg_phy <- ordinate(sub_algae_phylum, "NMDS", "bray")
## To plot the OTUs
plot_ordination(sub_algae_phylum, multi_ord_alg_phy, type="taxa", color="class", shape= "phylum", 
                title="OTUs for the algal phylum level")
## Since it is a confusing matter, we will make it simple and easy by breaking it down by taxonomic subdivisions
plot_ordination(sub_algae_phylum, multi_ord_alg_phy, type="taxa", color="class", 
                title="OTUs for the algal phylum level", label="class") + 
  facet_wrap(~phylum, 3)
## Now we will display samples and enlarge the points to make them more legible
plot_ordination(sub_algae_phylum, multi_ord_alg_phy, type="samples", color="treatment", 
                shape="phylum", title="Display of the samples") + geom_point(size=3)
## Display of samples and OTUs in 2 different panels
plot_ordination(sub_algae_phylum, multi_ord_alg_phy, type="split", color="class", 
                shape="phylum", title="Biplot of samples and taxa at phylum and class level.", label = "header_tanks_block") +  
  geom_point(size=3)


## For the genus level of algae
multi_ord_alg_gen <- ordinate(sub_algae_genus, "NMDS", "bray")
## To plot the OTUs
plot_ordination(sub_algae_genus, multi_ord_alg_gen, type="taxa", color="class", shape= "genus", 
                title="OTUs for the algal genus level")
## Since it is a confusing matter, we will make it simple and easy by breaking it down by taxonomic subdivisions
plot_ordination(sub_algae_genus, multi_ord_alg_gen, type="taxa", color="class", 
                title="OTUs for the algal genus level", label="class") + 
  facet_wrap(~genus, 3)
## Now we will display samples and enlarge the points to make them more legible
plot_ordination(sub_algae_genus, multi_ord_alg_gen, type="samples", color="treatment", 
                shape="genus", title="Display of the samples") + geom_point(size=3)
## Display of samples and OTUs in 2 different panels
plot_ordination(sub_algae_genus, multi_ord_alg_gen, type="split", color="class", 
                shape="genus", title="Biplot of samples and taxa at genus and class level.", label = "header_tanks_block") +  
  geom_point(size=3)

## For the phylum level of protists
multi_ord_prot_phy <- ordinate(sub_protist_phylum, "NMDS", "bray")
## To plot the OTUs
plot_ordination(sub_protist_phylum, multi_ord_prot_phy, type="taxa", color="class", shape= "phylum", 
                title="OTUs for the protists phylum level")
## Since it is a confusing matter, we will make it simple and easy by breaking it down by taxonomic subdivisions
plot_ordination(sub_protist_phylum, multi_ord_prot_phy, type="taxa", color="class", 
                title="OTUs for the protists phylum level", label="class") + 
  facet_wrap(~phylum, 3)
## Now we will display samples and enlarge the points to make them more legible
plot_ordination(sub_protist_phylum, multi_ord_prot_phy, type="samples", color="treatment", 
                shape="phylum", title="Display of the samples") + geom_point(size=3)
## Display of samples and OTUs in 2 different panels
plot_ordination(sub_protist_phylum, multi_ord_prot_phy, type="split", color="class", 
                shape="phylum", title="Biplot of samples and taxa at phylum and class level.", label = "header_tanks_block") +  
  geom_point(size=3)

## For the genus level of protists
multi_ord_prot_gen <- ordinate(sub_protist_genus, "NMDS", "bray")
## To plot the OTUs
plot_ordination(sub_protist_genus, multi_ord_prot_gen, type="taxa", color="class", shape= "genus", 
                title="OTUs for the protists genus level")
## Since it is a confusing matter, we will make it simple and easy by breaking it down by taxonomic subdivisions
plot_ordination(sub_protist_genus, multi_ord_prot_gen, type="taxa", color="class", 
                title="OTUs for the protists genus level", label="class") + 
  facet_wrap(~genus, 3)
## Now we will display samples and enlarge the points to make them more legible
plot_ordination(sub_protist_genus, multi_ord_prot_gen, type="samples", color="treatment", 
                shape="genus", title="Display of the samples") + geom_point(size=3)
## Display of samples and OTUs in 2 different panels
plot_ordination(sub_protist_genus, multi_ord_prot_gen, type="split", color="class", 
                shape="genus", title="Biplot of samples and taxa at genus and class level.", label = "header_tanks_block") +  
  geom_point(size=3)


## The network analysis

## The simple network analysis
plot_net(sub_algae_phylum, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="phylum")

plot_net(sub_algae_genus, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="genus")

plot_net(sub_protist_phylum, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="phylum")

plot_net(sub_protist_genus, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="genus")

## It is a bit confusing. 
## We want to make it simpler by using only the most important OTUs

total_1 = median(sample_sums(sub_algae_phylum))
total_2 = median(sample_sums(sub_algae_genus))
total_3 = median(sample_sums(sub_protist_phylum))
total_4 = median(sample_sums(sub_protist_genus))

sub_algae_phylum_1 <- filter_taxa(sub_algae_phylum, function(x) sum(x > total_1*0.20) > 0, TRUE)
sub_algae_phylum_1

sub_algae_genus_1 <- filter_taxa(sub_algae_genus, function(x) sum(x > total_2*0.20) > 0, TRUE)
sub_algae_genus_1

sub_protist_phylum_1 <- filter_taxa(sub_protist_phylum, function(x) sum(x > total_3*0.20) > 0, TRUE)
sub_protist_phylum_1

sub_protist_genus_1 <- filter_taxa(sub_protist_genus, function(x) sum(x > total_4*0.20) > 0, TRUE)
sub_protist_genus_1

plot_net(sub_algae_phylum_1, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="phylum")

plot_net(sub_algae_genus_1, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="genus")

plot_net(sub_protist_phylum_1, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="phylum")

plot_net(sub_protist_genus_1, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="genus")


# Linear Modelling
## According to Allen (2019) and Powell (2019), this statistical analysis was performed.
dat_lm <- data.frame(sample_data(phy_alpha_prune))

alpha_lm <- cbind(alpha.diversity, dat_lm)
alpha_lm

df2 = as.data.frame(sapply(alpha_lm, as.numeric))
df2
a <- cor(df2)
a
b <- as.data.frame(a)
b

## We can run the linear models with all factors manipulated 
## Here we will run linear models to determine the effects of all manipulated factors on alpha diversity and overall biodiversity of eukaryotes such as algae and protists.

## For alpha diversity
lin_mod_alpha_Observed <- lm(Observed ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, b)
summary(lin_mod_alpha_Observed)

lin_mod_alpha_Chao1 <- lm(Chao1 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, b)
summary(lin_mod_alpha_Chao1)

lin_mod_alpha_Shannon <- lm(Shannon ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, b)
summary(lin_mod_alpha_Shannon)

lin_mod_alpha_Simpson <- lm(Simpson ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, b)
summary(lin_mod_alpha_Simpson)

lin_mod_alpha_InvSimpson <- lm(InvSimpson ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, b)
summary(lin_mod_alpha_InvSimpson)

## For overall biodiversity of eukaryotes such as algae and protists

df3 = as.data.frame(sapply(samples_df, as.numeric))
df3
c <- cor(df3)
c
d <- as.data.frame(c)
d

lm_genus_alg_pro <- d
lm_genus_alg_pro

## Linear modeling for all manipulated factors of algae by genus level
lm_alg_gen_Messastrum <- lm(Messastrum ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Messastrum)

lm_alg_gen_Melosira <- lm(Melosira ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Melosira)

lm_alg_gen_Poteriospumella <- lm(Poteriospumella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Poteriospumella)

lm_alg_gen_Pseudellipsoidion <- lm(Pseudellipsoidion ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Pseudellipsoidion)

lm_alg_gen_Symbiodinium <- lm(Symbiodinium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Symbiodinium)

lm_alg_gen_Nasturtium <- lm(Nasturtium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Nasturtium)

lm_alg_gen_unclassified_Rigifilida <- lm(unclassified_Rigifilida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_unclassified_Rigifilida)

lm_alg_gen_Phalansterium <- lm(Phalansterium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Phalansterium)

lm_alg_gen_Reynoutria <- lm(Reynoutria ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Reynoutria)

lm_alg_gen_Pleurochloris <- lm(Pleurochloris ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_alg_gen_Pleurochloris)

## Linear modeling for all manipulated factors of protists by genus level
lm_prot_gen_Gymnophrys <- lm(Gymnophrys ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Gymnophrys)

lm_prot_gen_unclassified_Craspedida <- lm(unclassified_Craspedida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_Craspedida)

lm_prot_gen_Ichthyophthirius <- lm(Ichthyophthirius ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Ichthyophthirius)

lm_prot_gen_Cephalodella <- lm(Cephalodella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Cephalodella)

lm_prot_gen_Flabellula <- lm(Flabellula ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Flabellula)

lm_prot_gen_Bicosoeca <- lm(Bicosoeca ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Bicosoeca)

lm_prot_gen_unclassified_Discosea <- lm(unclassified_Discosea ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_Discosea)

lm_prot_gen_Lepidodermella <- lm(Lepidodermella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Lepidodermella)

lm_prot_gen_Sorodiplophrys <- lm(Sorodiplophrys ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Sorodiplophrys)

lm_prot_gen_Phytophthora <- lm(Phytophthora ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Phytophthora)

lm_prot_gen_unclassified_Apusomonadidae <- lm(unclassified_Apusomonadidae ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_Apusomonadidae)

lm_prot_gen_unclassified_CV1_B1_93 <- lm(unclassified_CV1_B1_93 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_CV1_B1_93)

lm_prot_gen_Protosporangium <- lm(Protosporangium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Protosporangium)

lm_prot_gen_Cryptosporidium <- lm(Cryptosporidium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_Cryptosporidium)

lm_prot_gen_unclassified_LEMD267 <- lm(unclassified_LEMD267 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_LEMD267)

lm_prot_gen_unclassified_Protosteliida <- lm(unclassified_Protosteliida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time, lm_genus_alg_pro)
summary(lm_prot_gen_unclassified_Protosteliida)


# According to Walker (2021), this statistical analysis was performed.
# Linear Mixed Models(LMM)
## Here, we will perform LMM to determine the blocking effect of header tanks on alpha diversity and total biodiversity of eukaryotes such as algae and protists.
## The environmental factors are the manipulated factors here.
## The random factor considered here is in the "header_tanks_block" column.
## For the first step, a random intercept is added for each treatment level.
## One message that may appear is "boundary (singular) fit: see help('isSingular')". 
## It does not describe that there is a problem with the fitted model.

## To fit the LMM

## For alpha diversity

## For random intercept only
lmm_block_alpha_Observed <- lmer(Observed ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = b)
summary(lmm_block_alpha_Observed)

lmm_block_alpha_Chao1 <- lmer(Chao1 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = b)
summary(lmm_block_alpha_Chao1)

lmm_block_alpha_Shannon <- lmer(Shannon ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = b)
summary(lmm_block_alpha_Shannon)

lmm_block_alpha_Simpson <- lmer(Simpson ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = b)
summary(lmm_block_alpha_Simpson)

lmm_block_alpha_InvSimpson <- lmer(InvSimpson ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = b)
summary(lmm_block_alpha_InvSimpson)

## For the total biodiversity of eukaryotes such as algae and protists

lmm_data_ap <- d

## LMM for all manipulated factors of algae by genus level for random intercept only
lmm_alg_gen_Messastrum_ri <- lmer(Messastrum ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Messastrum_ri)

lmm_alg_gen_Melosira_ri <- lmer(Melosira ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Melosira_ri)

lmm_alg_gen_Poteriospumella_ri <- lmer(Poteriospumella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Poteriospumella_ri)

lmm_alg_gen_Pseudellipsoidion_ri <- lmer(Pseudellipsoidion ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Pseudellipsoidion_ri)

lmm_alg_gen_Symbiodinium_ri <- lmer(Symbiodinium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Symbiodinium_ri)

lmm_alg_gen_Nasturtium_ri <- lmer(Nasturtium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Nasturtium_ri)

lmm_alg_gen_unclassified_Rigifilida_ri <- lmer(unclassified_Rigifilida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_unclassified_Rigifilida_ri)

lmm_alg_gen_Phalansterium_ri <- lmer(Phalansterium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Phalansterium_ri)

lmm_alg_gen_Reynoutria_ri <- lmer(Reynoutria ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Reynoutria_ri)

lmm_alg_gen_Pleurochloris_ri <- lmer(Pleurochloris ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_alg_gen_Pleurochloris_ri)

## LMM for all manipulated factors of protists by genus level for random intercept only
lmm_prot_gen_Gymnophrys_ri <- lmer(Gymnophrys ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Gymnophrys_ri)

lmm_prot_gen_unclassified_Craspedida_ri <- lmer(unclassified_Craspedida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_Craspedida_ri)

lmm_prot_gen_Ichthyophthirius_ri <- lmer(Ichthyophthirius ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Ichthyophthirius_ri)

lmm_prot_gen_Cephalodella_ri <- lmer(Cephalodella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Cephalodella_ri)

lmm_prot_gen_Flabellula_ri <- lmer(Flabellula ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Flabellula_ri)

lmm_prot_gen_Bicosoeca_ri <- lmer(Bicosoeca ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Bicosoeca_ri)

lmm_prot_gen_unclassified_Discosea_ri <- lmer(unclassified_Discosea ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_Discosea_ri)

lmm_prot_gen_Lepidodermella_ri <- lmer(Lepidodermella ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Lepidodermella_ri)

lmm_prot_gen_Sorodiplophrys_ri <- lmer(Sorodiplophrys ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Sorodiplophrys_ri)

lmm_prot_gen_Phytophthora_ri <- lmer(Phytophthora ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Phytophthora_ri)

lmm_prot_gen_unclassified_Apusomonadidae_ri <- lmer(unclassified_Apusomonadidae ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_Apusomonadidae_ri)

lmm_prot_gen_unclassified_CV1_B1_93_ri <- lmer(unclassified_CV1_B1_93 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_CV1_B1_93_ri)

lmm_prot_gen_Protosporangium_ri <- lmer(Protosporangium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Protosporangium_ri)

lmm_prot_gen_Cryptosporidium_ri <- lmer(Cryptosporidium ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_Cryptosporidium_ri)

lmm_prot_gen_unclassified_LEMD267_ri <- lmer(unclassified_LEMD267 ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_LEMD267_ri)

lmm_prot_gen_unclassified_Protosteliida_ri <- lmer(unclassified_Protosteliida ~ nutrient + sediment + flow + time + nutrient*sediment + nutrient*flow + sediment*flow + nutrient*time + sediment*time + flow*time + nutrient*sediment*flow + nutrient*sediment*time + nutrient*flow*time + sediment*flow*time + nutrient*sediment*flow*time + (1|header_tanks_block), data = lmm_data_ap)
summary(lmm_prot_gen_unclassified_Protosteliida_ri)


## According to Ollberding (2019), this statistical analysis was performed.
## Alpha-diversity

## To create a phy_tree slot from the original phyloseq object
phy_tree = rtree(ntaxa(phyloseq_obj), rooted=TRUE, tip.label=taxa_names(phyloseq_obj))
plot(phy_tree)

## To make a new phyloseq object with the newly created phy_tree
object_phy <- phyloseq(OTU, TAX, samples, phy_tree)
object_phy

## To prune OTUs that are not present in any of the samples
object_phy_new <- prune_taxa(taxa_sums(object_phy) > 0 , object_phy)

## Based on the plug-in estimates of observed richness, Shannon diversity, and phylogenetic diversity of the subsamples, we will attempt to calculate and examine differences between treatments.

ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(object_phy_new),
                         "observed" = phyloseq::estimate_richness(object_phy_new, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")
## As expected, the number of OTUs observed correlates with the total count of reads.

## We will subsample, plot and test for various group differences now.

## To Subsample the reads
phyloseq_rare <- phyloseq::rarefy_even_depth(object_phy_new, rngseed = 123, replace = FALSE) 

head(phyloseq::sample_sums(phyloseq_rare))

## To create a data.frame with adiv measures
adiv_data <- data.frame(
  "Observed" = phyloseq::estimate_richness(phyloseq_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phyloseq_rare, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(phyloseq_rare)))), tree = phyloseq::phy_tree(phyloseq_rare))[, 1],
  "Status" = phyloseq::sample_data(phyloseq_rare)$treatment)
head(adiv_data)

## To plot adiv_data measures
adiv_data %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

## To summarize
adiv_data %>%
  group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))


# According to Mariadassou et al., (2016), this statistical analysis was performed.
## BETA-DIVERSITY INDICES

## To create a phy_tree slot from the original phyloseq object
random_tree = rtree(ntaxa(phyloseq_obj), rooted=TRUE, tip.label=taxa_names(phyloseq_obj))
plot(random_tree)

## To make a new phyloseq object with the newly created phy_tree
phyloseq_obj_1 <- phyloseq(OTU, TAX, samples, random_tree)
phyloseq_obj_1

## To prune OTUs that are not present in any of the samples
phy_beta_prune <- prune_taxa(taxa_sums(phyloseq_obj_1) > 0 , phyloseq_obj_1)

## The generic method for distance
dist.bc <- phyloseq::distance(phy_beta_prune, method = "bray")
dist.bcc <- as.matrix(dist.bc)
dist.bcc
head(dist.bcc)[,1:6]

## To know all available distances
distanceMethodList
## or distance("list") in the previous version of phyloseq

## To compute other matrix distances
dist.jac <- phyloseq::distance(phy_beta_prune, method = "cc")
dist.jacc <- as.matrix(dist.jac)
dist.jacc
head(dist.jacc)[,1:6]

## The phy_tree slot is needed for unifrac and weighted unifrac matrix distances
dist.uf <- phyloseq::distance(phy_beta_prune, method = "unifrac") ## Unifrac
dist.uff <- as.matrix(dist.uf)
dist.uff
head(dist.uff)[,1:6]

dist.wuf <- phyloseq::distance(phy_beta_prune, method = "wunifrac") ## Weighted Unifrac
dist.wuff <- as.matrix(dist.wuf)
dist.wuff
head(dist.wuff)[,1:6]

## Just some codes to visualize the distances
sampleOrder <- levels(reorder(sample_names(phy_beta_prune), as.numeric(get_variable(phy_beta_prune, "treatment")))) 
plot_dist_as_heatmap <- function(dist, order = sampleOrder, title = NULL) {
  data <- melt(as(dist, "matrix"))
  colnames(data) <- c("x", "y", "distance")
  if (!is.null(order)) {
    data$x <- factor(data$x, levels = order)
    data$y <- factor(data$y, levels = order)
  }
  p4 <- ggplot(data, aes(x = x, y = y, fill = distance)) + geom_tile() 
  p5 <- p4 + theme(axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.text.y = element_blank()
  )
  p6 <- p5 + scale_fill_continuous(limits = c(0, 1))
  if (!is.null(title)) {
    p7 <- p6 + ggtitle(title)
  }
  return(p7)
}

## To compare Bray-Curtis (bray) and Jaccard (cc)
p8 <- plot_dist_as_heatmap(dist.bcc, title = "Bray-Curtis")
plot(p8)
p9 <- plot_dist_as_heatmap(dist.jacc, title = "Jaccard")
plot(p9)
grid.arrange(p8, p9, ncol = 2, nrow = 1)

## To compare unifrac and wunifrac
p10 <- plot_dist_as_heatmap(dist.uff, title = "Unifrac")
plot(p10)
p11 <- plot_dist_as_heatmap(dist.wuff, title = "wUnifrac")
plot(p11)
grid.arrange(p10, p11, ncol = 2, nrow = 1)

## To compare Bray-Curtis, Jaccard, Unifrac and wUnifrac
grid.arrange(p8, p9, p10, p11, ncol = 2, nrow = 2)



# According to Ollberding (2019), this statistical analysis was performed.
## From one microbial composition to another, beta diversity gives an estimate of similarity or dissimilarity.
## The differential pairwise Euclidean distance between samples is the general estimator of microbial beta diversity.

## We create the beta diversity ordination by applying the Aitchison distance.
## PCA (Principal Component Analysis) should only be applied to centered logarithmic ratio (CLR) transformed counts.
## The microbiome package would be applied and also used to assign a pseudo count of 1 to support the transformation.
## Because the logarithm of zero is undefined.

## CLR transformation
phy_beta_clr <- microbiome::transform(phyloseq_obj_1, "clr")
phy_beta_clr
phyloseq::otu_table(phyloseq_obj_1)[1:5, 1:5]
phyloseq::otu_table(phy_beta_clr)[1:5, 1:5]
## The values here are no longer defined as count values.
## Here we see that the dominance of each taxa relative to the geometric mean of the total taxa on a logarithmic scale.

## PCA (Princial Component Analysis)
## We perform PCA, test the relative importance of each principal component, and create the ordination.
## PCA is defined as an unsupervised learning technique that determines the similarities between different samples with large features.
## PCA helps to map the data in a lower dimensional space. 
## PCA identifies the variables called principal components (PCs) that capture the greatest amount of information possible.
## The information is considered as the amount of variation in the data.
## We can only represent the samples in a two- or three-dimensional space, where most microbiome analyses record data from the first few PCs.

## PCA via phyloseq

pca_ord_clr <- phyloseq::ordinate(phy_beta_clr, "RDA")

## To Plot the scree plot

phyloseq::plot_scree(pca_ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

## To Examine eigenvalues and % prop. variance explained

head(pca_ord_clr$CA$eig)   

sapply(pca_ord_clr$CA$eig[1:5], function(x) x / sum(pca_ord_clr$CA$eig))     
## PCA is the RDA without constraints.
## The PCs created by applying the phyloseq::ordinate function.
## The Scree Plot is used to test the exact proportion of the total variation described by each PC.
## The first PC stands out, and after that the other components gradually decrease.

## Here we will present the first two components.
## Also, the diagram is scaled so that the necessary information is displayed on each axis.

## To scale axes and plot ordination

c_clr1 <- pca_ord_clr$CA$eig[1] / sum(pca_ord_clr$CA$eig)
c_clr2 <- pca_ord_clr$CA$eig[2] / sum(pca_ord_clr$CA$eig)
phyloseq::plot_ordination(phyloseq_obj_1, pca_ord_clr, type="sample", color="treatment") + 
  geom_point(size = 2) +
  coord_fixed(c_clr2 / c_clr1) +
  stat_ellipse(aes(group = treatment), linetype = 2)

phyloseq::plot_ordination(phyloseq_obj_1, pca_ord_clr, type="sample", color="treatment") + 
  geom_point(size = 2) +
  coord_fixed(c_clr2 / c_clr1) +
  stat_ellipse(aes(group = week), linetype = 2)
## Based on the different sample types, there are some differences and overlaps in treatments that are reflected in the community.


# Permutational multivariate analysis of variance (PERMANOVA) 
## Because PCA is a tool for visualizing explanatory data, PERMANOVA can be used to test whether different samples are clustering beyond that expected from sampling variability.
## This can be done by partitioning the sums of squares of the components in and between clusters based on centroid theories.
## Various permutations of the data, such as random shuffling, are used to create the null distribution.
## ADONIS tests could be affected by differences in dispersion or spread.

## To generate the distance matrix

phy_clr_dist_matrix <- phyloseq::distance(phy_beta_clr, method = "euclidean")

# To do the ADONIS test

vegan::adonis2(phy_clr_dist_matrix ~ phyloseq::sample_data(phy_beta_clr)$treatment)

adonis2(dist.wuff ~ sample_data(phy_beta_prune)$treatment)

set.seed(1)

# Calculate bray curtis distance matrix
erie_bray <- phyloseq::distance(phy_beta_prune, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phy_beta_prune))

# Adonis test
adonis2(erie_bray ~ treatment, data = sampledf)

set.seed(1)

# Calculate bray curtis distance matrix
erie_bray_1 <- phyloseq::distance(phy_beta_clr, method = "bray")

# make a data frame from the sample_data
sampledf_1 <- data.frame(sample_data(phy_beta_clr))

# Adonis test
adonis2(erie_bray_1 ~ treatment, data = sampledf_1)




# Permutational analysis of multivariate dispersions (PERMDISP)
## Homogeneity of multivariate dispersions
## To do the dispersion test (PERMDISP) and plot

dispersion <- vegan::betadisper(phy_clr_dist_matrix, phyloseq::sample_data(phy_beta_clr)$treatment)
dispersion

plot(dispersion, main = "Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion)

## Homogeneity of dispersion test
beta_dispersion <- betadisper(erie_bray, sampledf$treatment)
permutest(beta_dispersion)

beta_dispersion_1 <- betadisper(erie_bray_1, sampledf_1$treatment)
permutest(beta_dispersion_1)


# Multivariate Analysis Of Variance (MANOVA)
## According to Flores et al., (2011), this statistical analysis was performed.
## adonis (a permutation MANOVA, also in the vegan package)
## We will use here several dependent variables for the MANOVA
manova_df = as(sample_data(phy_filtered), "data.frame")
manova_d = dist.bcc
phy_filtered_adonis = adonis2(manova_d ~ treatment + flow + sediment + nutrient, manova_df)
phy_filtered_adonis



## To provide information in the form of a phylogenetic tree on beta diversity analysis.
## The UniFrac metric includes phylogenetic information by measuring the total branch lengths that are not shared between samples divided by the total branch lengths.
## Between samples and sample types, this technique expresses differences in phylogenetic relatedness.
## Calculation of weighted and unweighted UniFrac metrics with PCoA.
## PCoA can be considered as PCA for the non-Euclidean measures.

## To subsample the reads

sub_ps_rare <- phyloseq::rarefy_even_depth(phy_beta_prune, rngseed = 123, replace = FALSE)

## To generate distances

g_ord_unifrac <- ordinate(sub_ps_rare, method = "PCoA", distance = "wunifrac") 
g_ord_unifrac_un <- ordinate(sub_ps_rare, method = "PCoA", distance = "unifrac")   

## To plot ordinations

m <- plot_ordination(sub_ps_rare, g_ord_unifrac, color = "treatment") + geom_point(size = 2)
n <- plot_ordination(sub_ps_rare, g_ord_unifrac_un, color = "treatment") + geom_point(size = 2)
cowplot::plot_grid(m, n, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))
## There is some overlap in treatment types in these two Unifrac distance plots.
## It measures the relative abundance of each taxa in the community.
## On the axis of unweighted Unifrac distance, there is a clustering that cannot be explained by treatment.



List of references

Ollberding, N.J. (2019, Jul 28). Introduction to the Statistical Analysis of Microbiome Data in R. Wowchemy. https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

Flores, G. E., Bates, S. T., Knights, D., Lauber, C. L., Stombaugh, J., Knight, R., & Fierer, N. (2011). Microbial biogeography of public restroom surfaces. PloS one, 6(11), e28132.

Allen, D.R. (2019). Ecological Diversity. RPubs by RStudio. M.Res. Marine Biology 2019. https://www.rpubs.com/roalle/mres_2019

Powell, J. (2019, Sep 23). Analyses of alpha diversity. http://www.hiercourse.com/docs/microbial/02_alphaDiversity.html

Walker, J.A. (2021, Sep 25). Applied Statistics for Experimental Biology. Elements of Applied Biostatistics. https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html

Vaulot, D. (2021, Feb 15). Phyloseq tutorial. https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#content

Mariadassou, M., Bernard, M., Pascal, G., Cauquil, L., & Chaillou, S. (2016). Analysis of community composition data using phyloseq. Montpellier Décembre 2016. https://genoweb.toulouse.inra.fr/~formation/15_FROGS/8-February2017/FROGS_phyloseq_02_2017.pdf

