# Analysis of eukaryotic biodiversity data using the phyloseq package in R (Week 2 and week 3)

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
install.packages("broom")
install.packages("interactions")
install.packages("MuMIn")
install.packages("blmeco")
install.packages("tseries")

## To load the installed packages with the pacman package
pacman::p_load(pacman, dplyr, readxl, ggplot2, phyloseq, microbiome, reshape2, ape, gridExtra, plotly, vegan, dendextend, tidyr, rms, effectsize, lme4, picante, cowplot, here, lsr, MASS, rcompanion, ggiraphExtra, optimx, mvabund, gllvm, radiant.data, tidyverse, MuMIn, corrplot, gclus, coefplot, jtools, ggstance, interactions, blmeco, tseries)
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

## To idem for the two other matrixes 
tax_mat_1 <- tax_mat_1 %>% 
  tibble::column_to_rownames("otu")

tax_mat_2 <- tax_mat_2 %>% 
  tibble::column_to_rownames("otu")

## For sample data
samples_df_1 <- samples_df_1 %>% 
  tibble::column_to_rownames("sample")

samples_df_2 <- samples_df_2 %>% 
  tibble::column_to_rownames("sample")

## Transformation to matrices otu and tax tables (the sample table can be left as a data frame)
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

## The preprocessing of the data

## log10p transformation
phyloseq_object_final <- microbiome::transform(phyloseq_object_with_tree, "log10p")
phyloseq_object_final
phyloseq::otu_table(phyloseq_object_final)[1:5, 1:5]
phyloseq::sample_data(phyloseq_object_final)[1:5, 1:5]

## To normalize or standardize the number of reads or abundances in each sample based on the median sequencing depth
mydata_total = median(sample_sums(phyloseq_object_final))
mydata_standf = function(x, t=mydata_total) round(t * (x / sum(x)))

## Then, the phyloseq_object_final data is first transformed into relative abundance data and stored under phy_transformed.
## After that, the phy_transformed object is filtered to phy_filtered to keep the OTUs with a mean >= 0.003.
phy_transformed  = transform_sample_counts(phyloseq_object_final, mydata_standf)
phy_transformed
phy_filtered = filter_taxa(phy_transformed, function(x) mean(x) >= 0.003, TRUE) 
phy_filtered

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

## To keep only taxa of interests (Algae) according to phylum and genus level
sub_algae_phylum <- subset_taxa(phy_filtered, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Phragmoplastophyta", "Rigifilida", "Schizoplasmodiida", "Streptophyta", "Xanthophyceae", "Nucleariidae_and_Fonticula_group"))
sub_algae_genus <- subset_taxa(phy_filtered, genus %in% c("Messastrum", "Melosira", "Poteriospumella", "Pseudellipsoidion", "Symbiodinium", "Nasturtium", "unclassified_Rigifilida", "Phalansterium", "Reynoutria", "Pleurochloris", "Mallomonas", "Makinoella", "Nuclearia"))

## To keep only taxa of interests (Protists) according to phylum and genus level
sub_protist_phylum <- subset_taxa(phy_filtered, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Rotifera", "Tubulinea", "Bicosoecida", "Discosea", "Gastrotricha", "Labyrinthulomycetes", "Peronosporomycetes", "Apusomonadidae", "CV1-B1-93", "Protosporangiida", "Apicomplexa", "Gracilipodida", "Protosteliida", "Colponemidia", "MAST-12"))
sub_protist_genus <- subset_taxa(phy_filtered, genus %in% c("Gymnophrys", "unclassified_Craspedida", "Ichthyophthirius", "Cephalodella", "Flabellula", "Bicosoeca", "unclassified_Discosea", "Lepidodermella", "Sorodiplophrys", "Phytophthora", "unclassified_Apusomonadidae", "unclassified_CV1-B1-93", "Protosporangium", "Cryptosporidium", "unclassified_LEMD267", "unclassified_Protosteliida", "Colpoda"))

## Bar graphs of relative abundance

## Basic bar graphs based on phylum and genus level for Algae
bar_algae_phylum <- microbiome::transform(sub_algae_phylum, "compositional")
plot_bar(bar_algae_phylum, fill="phylum", title = "The relative abundances of the algae at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_algae_genus <- microbiome::transform(sub_algae_genus, "compositional")
plot_bar(bar_algae_genus, fill="genus", title = "The relative abundances of the algae at the level of the genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

## Basic bar graphs based on phylum and genus level for Protists
bar_protist_phylum <- microbiome::transform(sub_protist_phylum, "compositional")
plot_bar(bar_protist_phylum, fill="phylum", title = "The relative abundances of protists at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

bar_protist_genus <- microbiome::transform(sub_protist_genus, "compositional")
plot_bar(bar_protist_genus, fill="genus", title = "The relative abundances of protists at the level of the genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

# According to Mariadassou et al., (2016), this statistical analysis was performed.
## Exploring the biodiversity : Alpha-diversity
## The overall richness are plotted with plot_richness function
## We are interested in alpha diversity

## For the algae phylum level

alpha.diversity_algae_phylum <- estimate_richness(sub_algae_phylum, measures = c("Observed", "Shannon"))
alpha.diversity_algae_phylum

alpha_algae_phylum <- plot_richness(sub_algae_phylum, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
      geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_phylum)

pielou_algae_phylum <- evenness(sub_algae_phylum, 'pielou')
pielou_algae_phylum

dat_lm_algae_phylum <- data.frame(sample_data(sub_algae_phylum))

alpha_lm_algae_phylum <- cbind(alpha.diversity_algae_phylum, pielou_algae_phylum, dat_lm_algae_phylum)
alpha_lm_algae_phylum

plot_pielou_algae_phylum <- ggplot(alpha_lm_algae_phylum, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_phylum                                        

## For the algae genus level

alpha.diversity_algae_genus <- estimate_richness(sub_algae_genus, measures = c("Observed", "Shannon"))
alpha.diversity_algae_genus

alpha_algae_genus <- plot_richness(sub_algae_genus, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_genus)

pielou_algae_genus <- evenness(sub_algae_genus, 'pielou')
pielou_algae_genus

dat_lm_algae_genus <- data.frame(sample_data(sub_algae_genus))

alpha_lm_algae_genus <- cbind(alpha.diversity_algae_genus, pielou_algae_genus, dat_lm_algae_genus)
alpha_lm_algae_genus

plot_pielou_algae_genus <- ggplot(alpha_lm_algae_genus, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_genus      

## For the protist phylum level

alpha.diversity_protist_phylum <- estimate_richness(sub_protist_phylum, measures = c("Observed", "Shannon"))
alpha.diversity_protist_phylum

alpha_protist_phylum <- plot_richness(sub_protist_phylum, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_phylum)

pielou_protist_phylum <- evenness(sub_protist_phylum, 'pielou')
pielou_protist_phylum

dat_lm_protist_phylum <- data.frame(sample_data(sub_protist_phylum))

alpha_lm_protist_phylum <- cbind(alpha.diversity_protist_phylum, pielou_protist_phylum, dat_lm_protist_phylum)
alpha_lm_protist_phylum

plot_pielou_protist_phylum <- ggplot(alpha_lm_protist_phylum, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_phylum 

## For the protist genus level

alpha.diversity_protist_genus <- estimate_richness(sub_protist_genus, measures = c("Observed", "Shannon"))
alpha.diversity_protist_genus

alpha_protist_genus <- plot_richness(sub_protist_genus, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_genus)

pielou_protist_genus <- evenness(sub_protist_genus, 'pielou')
pielou_protist_genus

dat_lm_protist_genus <- data.frame(sample_data(sub_protist_genus))

alpha_lm_protist_genus <- cbind(alpha.diversity_protist_genus, pielou_protist_genus, dat_lm_protist_genus)
alpha_lm_protist_genus

plot_pielou_protist_genus <- ggplot(alpha_lm_protist_genus, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_genus


# Generalized Linear Modelling (GLM)
# According to ZACH (2021), Humboldt-Universität zu Berlin | Geography Department (2021), Phillips (2018), this statistical analysis was performed.

## For the algae phylum level
Shannon_1 <- ceiling(alpha_lm_algae_phylum$Shannon)
pielou_1 <- ceiling(ifelse(test = alpha_lm_algae_phylum$pielou > 0.8966068, yes = alpha_lm_algae_phylum$pielou + 0.01, no = alpha_lm_algae_phylum$pielou))
alpha_lm_algae_phylum_1 <- cbind(alpha_lm_algae_phylum, Shannon_1, pielou_1)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_phylum_1$Observed)

glin_mod_algae_phylum_alpha_Observed <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1)
summary(glin_mod_algae_phylum_alpha_Observed)
eta_squared(aov(glin_mod_algae_phylum_alpha_Observed))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_algae_phylum_1$Shannon_1)

glin_mod_algae_phylum_alpha_Shannon <- glm(Shannon_1 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1)
summary(glin_mod_algae_phylum_alpha_Shannon)
eta_squared(aov(glin_mod_algae_phylum_alpha_Shannon))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_phylum_1$pielou_1)

glin_mod_algae_phylum_alpha_pielou <- glm(pielou_1 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1)
summary(glin_mod_algae_phylum_alpha_pielou)
eta_squared(aov(glin_mod_algae_phylum_alpha_pielou))

## For the algae genus level
alpha_lm_algae_genus
Shannon_2 <- ceiling(alpha_lm_algae_genus$Shannon)

alpha_lm_algae_genus_mu <- alpha_lm_algae_genus %>% 
  mutate(pielou = replace_na(pielou,mean(pielou, na.rm = TRUE)))

pielou_2 <- ceiling(ifelse(test = alpha_lm_algae_genus_mu$pielou > 0.8774688, yes = alpha_lm_algae_genus_mu$pielou + 0.01, no = alpha_lm_algae_genus_mu$pielou))
alpha_lm_algae_genus_1 <- cbind(alpha_lm_algae_genus, Shannon_2, pielou_2)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_algae_genus_1$Observed)

glin_mod_algae_genus_alpha_Observed <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
summary(glin_mod_algae_genus_alpha_Observed)
eta_squared(aov(glin_mod_algae_genus_alpha_Observed))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_genus_1$Shannon_2)

glin_mod_algae_genus_alpha_Shannon <- glm(Shannon_2 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
summary(glin_mod_algae_genus_alpha_Shannon)
eta_squared(aov(glin_mod_algae_genus_alpha_Shannon))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_genus_1$pielou_2)

glin_mod_algae_genus_alpha_pielou <- glm(pielou_2 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
summary(glin_mod_algae_genus_alpha_pielou)
eta_squared(aov(glin_mod_algae_genus_alpha_pielou))

## For the protist phylum level
alpha_lm_protist_phylum
Shannon_3 <- ceiling(alpha_lm_protist_phylum$Shannon)
pielou_3 <- ceiling(ifelse(test = alpha_lm_protist_phylum$pielou > 0.8817381, yes = alpha_lm_protist_phylum$pielou + 0.01, no = alpha_lm_protist_phylum$pielou))
alpha_lm_protist_phylum_1 <- cbind(alpha_lm_protist_phylum, Shannon_3, pielou_3)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_phylum_1$Observed)

glin_mod_protist_phylum_alpha_Observed <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1)
summary(glin_mod_protist_phylum_alpha_Observed)
eta_squared(aov(glin_mod_protist_phylum_alpha_Observed))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_phylum_1$Shannon_3)

glin_mod_protist_phylum_alpha_Shannon <- glm(Shannon_3 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1)
summary(glin_mod_protist_phylum_alpha_Shannon)
eta_squared(aov(glin_mod_protist_phylum_alpha_Shannon))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_protist_phylum_1$pielou_3)

glin_mod_protist_phylum_alpha_pielou <- glm(pielou_3 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1)
summary(glin_mod_protist_phylum_alpha_pielou)
eta_squared(aov(glin_mod_protist_phylum_alpha_pielou))

## For the protist genus level
alpha_lm_protist_genus
Shannon_4 <- ceiling(alpha_lm_protist_genus$Shannon)

alpha_lm_protist_genus_mu <- alpha_lm_protist_genus %>% 
  mutate(pielou = replace_na(pielou,mean(pielou, na.rm = TRUE)))

pielou_4 <- ceiling(ifelse(test = alpha_lm_protist_genus_mu$pielou > 0.8974600, yes = alpha_lm_protist_genus_mu$pielou + 0.01, no = alpha_lm_protist_genus_mu$pielou))
alpha_lm_protist_genus_1 <- cbind(alpha_lm_protist_genus, Shannon_4, pielou_4)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_genus_1$Observed)

glin_mod_protist_genus_alpha_Observed <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1)
summary(glin_mod_protist_genus_alpha_Observed)
eta_squared(aov(glin_mod_protist_genus_alpha_Observed))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_genus_1$Shannon_4)

glin_mod_protist_genus_alpha_Shannon <- glm(Shannon_4 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1)
summary(glin_mod_protist_genus_alpha_Shannon)
eta_squared(aov(glin_mod_protist_genus_alpha_Shannon))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_protist_genus_1$pielou_4)

glin_mod_protist_genus_alpha_pielou <- glm(pielou_4 ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1)
summary(glin_mod_protist_genus_alpha_pielou)
eta_squared(aov(glin_mod_protist_genus_alpha_pielou))


# Generalized Linear Mixed Models(GLMM)
## According to Walker (2021), Bolker (2018), this statistical analysis was performed.
## The GLMM not only includes the random intercept of each level of nutrient, sediment, flow and time for block effects but also the random slope for each level of nutrient, sediment, flow and time here.

## For the algae phylum level

## For Observed
glmm_algae_phylum_block_alpha_Observed_ris <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_algae_phylum_1)
print(summary(glmm_algae_phylum_block_alpha_Observed_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_Observed_ris)

plot(glmm_algae_phylum_block_alpha_Observed_ris)

qqnorm(resid(glmm_algae_phylum_block_alpha_Observed_ris))
qqline(resid(glmm_algae_phylum_block_alpha_Observed_ris))

## For Shannon
glmm_algae_phylum_block_alpha_Shannon_ris <- glmer(Shannon_1 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_phylum_1)
print(summary(glmm_algae_phylum_block_alpha_Shannon_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_Shannon_ris)

plot(glmm_algae_phylum_block_alpha_Shannon_ris)

qqnorm(resid(glmm_algae_phylum_block_alpha_Shannon_ris))
qqline(resid(glmm_algae_phylum_block_alpha_Shannon_ris))

## For Pielou
glmm_algae_phylum_block_alpha_pielou_ris <- glmer(pielou_1 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_phylum_1)
print(summary(glmm_algae_phylum_block_alpha_pielou_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_pielou_ris)

plot(glmm_algae_phylum_block_alpha_pielou_ris)

qqnorm(resid(glmm_algae_phylum_block_alpha_pielou_ris))
qqline(resid(glmm_algae_phylum_block_alpha_pielou_ris))


## For the algae genus level

## For Observed
glmm_algae_genus_block_alpha_Observed_ris <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
print(summary(glmm_algae_genus_block_alpha_Observed_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_Observed_ris)

plot(glmm_algae_genus_block_alpha_Observed_ris)

qqnorm(resid(glmm_algae_genus_block_alpha_Observed_ris))
qqline(resid(glmm_algae_genus_block_alpha_Observed_ris))

## For Shannon
glmm_algae_genus_block_alpha_Shannon_ris <- glmer(Shannon_2 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
print(summary(glmm_algae_genus_block_alpha_Shannon_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_Shannon_ris)

plot(glmm_algae_genus_block_alpha_Shannon_ris)

qqnorm(resid(glmm_algae_genus_block_alpha_Shannon_ris))
qqline(resid(glmm_algae_genus_block_alpha_Shannon_ris))

## For Pielou
glmm_algae_genus_block_alpha_pielou_ris <- glmer(pielou_2 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1)
print(summary(glmm_algae_genus_block_alpha_pielou_ris), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_pielou_ris)

plot(glmm_algae_genus_block_alpha_pielou_ris)

qqnorm(resid(glmm_algae_genus_block_alpha_pielou_ris))
qqline(resid(glmm_algae_genus_block_alpha_pielou_ris))

## For the protist phylum level

## For Observed
glmm_protist_phylum_block_alpha_Observed_ris <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_protist_phylum_1)
print(summary(glmm_protist_phylum_block_alpha_Observed_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_Observed_ris)

plot(glmm_protist_phylum_block_alpha_Observed_ris)

qqnorm(resid(glmm_protist_phylum_block_alpha_Observed_ris))
qqline(resid(glmm_protist_phylum_block_alpha_Observed_ris))

## For Shannon
glmm_protist_phylum_block_alpha_Shannon_ris <- glmer(Shannon_3 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_phylum_1)
print(summary(glmm_protist_phylum_block_alpha_Shannon_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_Shannon_ris)

plot(glmm_protist_phylum_block_alpha_Shannon_ris)

qqnorm(resid(glmm_protist_phylum_block_alpha_Shannon_ris))
qqline(resid(glmm_protist_phylum_block_alpha_Shannon_ris))

## For Pielou
glmm_protist_phylum_block_alpha_pielou_ris <- glmer(pielou_3 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_phylum_1)
print(summary(glmm_protist_phylum_block_alpha_pielou_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_pielou_ris)

plot(glmm_protist_phylum_block_alpha_pielou_ris)

qqnorm(resid(glmm_protist_phylum_block_alpha_pielou_ris))
qqline(resid(glmm_protist_phylum_block_alpha_pielou_ris))


## For the protist genus level

## For Observed
glmm_protist_genus_block_alpha_Observed_ris <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_protist_genus_1)
print(summary(glmm_protist_genus_block_alpha_Observed_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_Observed_ris)

plot(glmm_protist_genus_block_alpha_Observed_ris)

qqnorm(resid(glmm_protist_genus_block_alpha_Observed_ris))
qqline(resid(glmm_protist_genus_block_alpha_Observed_ris))

## For Shannon
glmm_protist_genus_block_alpha_Shannon_ris <- glmer(Shannon_4 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_genus_1)
print(summary(glmm_protist_genus_block_alpha_Shannon_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_Shannon_ris)

plot(glmm_protist_genus_block_alpha_Shannon_ris)

qqnorm(resid(glmm_protist_genus_block_alpha_Shannon_ris))
qqline(resid(glmm_protist_genus_block_alpha_Shannon_ris))

## For Pielou
glmm_protist_genus_block_alpha_pielou_ris <- glmer(pielou_4 ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_genus_1)
print(summary(glmm_protist_genus_block_alpha_pielou_ris), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_pielou_ris)

plot(glmm_protist_genus_block_alpha_pielou_ris)

qqnorm(resid(glmm_protist_genus_block_alpha_pielou_ris))
qqline(resid(glmm_protist_genus_block_alpha_pielou_ris))


# BETA-DIVERSITY INDICES
# According to Ollberding (2019), this statistical analysis was performed.
# Permutational multivariate analysis of variance (PERMANOVA)
## At the algal phylum level, protist phylum level, and protist genus level, the Bray-Curtis distance matrix was used. 
## For the algal genus level, the Euclidean distance matrix was used because the Bray-Curtis distance matrix could not create all distance matrix values appropriately. 

## For the algae phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_algae_phylum <- phyloseq::distance(sub_algae_phylum, method = "bray")

## To make a data frame from the sample_data
perm_algae_phylum_sampledf <- data.frame(sample_data(sub_algae_phylum))

## To perform the adonis test
adonis2(bray_algae_phylum ~ nutrient*sediment*flow*time, data = perm_algae_phylum_sampledf)

## For the algae genus level
set.seed(1)

## To calculate the euclidean distance matrix
euclidean_algae_genus <- phyloseq::distance(sub_algae_genus, method = "euclidean")

## To make a data frame from the sample_data
perm_algae_genus_sampledf <- data.frame(sample_data(sub_algae_genus))

## To perform the adonis test
adonis2(euclidean_algae_genus ~ nutrient*sediment*flow*time, data = perm_algae_genus_sampledf)

## For the protist phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_phylum <- phyloseq::distance(sub_protist_phylum, method = "bray")

## To make a data frame from the sample_data
perm_protist_phylum_sampledf <- data.frame(sample_data(sub_protist_phylum))

## To perform the adonis test
adonis2(bray_protist_phylum ~ nutrient*sediment*flow*time, data = perm_protist_phylum_sampledf)

## For the protist genus level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_genus <- phyloseq::distance(sub_protist_genus, method = "bray")

## To make a data frame from the sample_data
perm_protist_genus_sampledf <- data.frame(sample_data(sub_protist_genus))

## To perform the adonis test
adonis2(bray_protist_genus ~ nutrient*sediment*flow*time, data = perm_protist_genus_sampledf)


# According to Ollberding (2019), this statistical analysis was performed.
# Permutational analysis of multivariate dispersions (PERMDISP)
## Homogeneity of multivariate dispersions
## To do the dispersion test (PERMDISP) and plot

## For the algae phylum level
dispersion_algae_phylum <- vegan::betadisper(bray_algae_phylum, phyloseq::sample_data(sub_algae_phylum)$treatment)
dispersion_algae_phylum

plot(dispersion_algae_phylum, main = "For the algae phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_phylum, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_phylum)

## For the algae genus level
dispersion_algae_genus <- vegan::betadisper(euclidean_algae_genus, phyloseq::sample_data(sub_algae_genus)$treatment)
dispersion_algae_genus

plot(dispersion_algae_genus, main = "For the algae genus level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_genus, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_genus)

## For the protist phylum level
dispersion_protist_phylum <- vegan::betadisper(bray_protist_phylum, phyloseq::sample_data(sub_protist_phylum)$treatment)
dispersion_protist_phylum

plot(dispersion_protist_phylum, main = "For the protist phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_phylum, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_phylum)

## For the protist genus level
dispersion_protist_genus <- vegan::betadisper(bray_protist_genus, phyloseq::sample_data(sub_protist_genus)$treatment)
dispersion_protist_genus

plot(dispersion_protist_genus, main = "For the protist genus level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_genus, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_genus)


## According to Xu and Yu (2020), this statistical analysis was performed.
# The PCoA analysis (Principal Coordinate Analysis)
library(MicrobiotaProcess)

## For the algae phylum level
pcoares_algae_phylum <- get_pcoa(obj=sub_algae_phylum, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_algae_phylum <- ggordpoint(obj=pcoares_algae_phylum, biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_algae_phylum <- ggordpoint(obj=pcoares_algae_phylum, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_algae_phylum | pcoaplot2_algae_phylum

## For the algae genus level
pcoares_algae_genus <- get_pcoa(obj=sub_algae_genus, distmethod="euclidean", method="hellinger")
## For visualizing the result
pcoaplot1_algae_genus <- ggordpoint(obj=pcoares_algae_genus, biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_algae_genus <- ggordpoint(obj=pcoares_algae_genus, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_algae_genus | pcoaplot2_algae_genus

## For the protist phylum level
pcoares_protist_phylum <- get_pcoa(obj=sub_protist_phylum, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_protist_phylum <- ggordpoint(obj=pcoares_protist_phylum, biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_protist_phylum <- ggordpoint(obj=pcoares_protist_phylum, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_protist_phylum | pcoaplot2_protist_phylum

## For the protist genus level
pcoares_protist_genus <- get_pcoa(obj=sub_protist_genus, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_protist_genus <- ggordpoint(obj=pcoares_protist_genus, biplot=TRUE, speciesannot=TRUE,
                                       factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_protist_genus <- ggordpoint(obj=pcoares_protist_genus, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                       factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_protist_genus | pcoaplot2_protist_genus


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

## The most abundant taxa
abundant_algae_phylum <- most_abundant_taxa(sub_algae_phylum,"phylum")
abundant_algae_phylum

abundant_algae_genus <- most_abundant_taxa(sub_algae_genus,"genus")
abundant_algae_genus

abundant_protist_phylum <- most_abundant_taxa(sub_protist_phylum,"phylum")
abundant_protist_phylum

abundant_protist_genus <- most_abundant_taxa(sub_protist_genus,"genus")
abundant_protist_genus

# Multivariate Analysis Of Variance (MANOVA)
## According to Kassambara (2017), Ben-Shachar M (2020), this statistical analysis was performed.

## For the most abundant taxa of the algal phylum level
manova_df_algae_phylum = as(sample_data(sub_algae_phylum), "data.frame")
manova_df_algae_phylum
manova_df_algae_phylum_1 = cbind(manova_df_algae_phylum, abundant_algae_phylum)
manova_df_algae_phylum_1

### For Bacillariophyta, Chlorophyta, Streptophyta and Phragmoplastophyta
manova_abundant_algae_phylum <- manova(cbind(Bacillariophyta, Chlorophyta, Streptophyta, Phragmoplastophyta) ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_1)
manova_abundant_algae_phylum
summary.aov(manova_abundant_algae_phylum)
eta_squared(aov(manova_abundant_algae_phylum))

## For the most abundant taxa of the algal genus level
manova_df_algae_genus = as(sample_data(sub_algae_genus), "data.frame")
manova_df_algae_genus
manova_df_algae_genus_1 = cbind(manova_df_algae_genus, abundant_algae_genus)
manova_df_algae_genus_1

### For Melosira, Makinoella, Messastrum, Pseudellipsoidion and Symbiodinium
manova_abundant_algae_genus <- manova(cbind(Melosira, Makinoella, Messastrum, Pseudellipsoidion, Symbiodinium) ~ nutrient*sediment*flow*time, data = manova_df_algae_genus_1)
manova_abundant_algae_genus
summary.aov(manova_abundant_algae_genus)
eta_squared(aov(manova_abundant_algae_genus))

## For the most abundant taxa of the protist phylum level
manova_df_protist_phylum = as(sample_data(sub_protist_phylum), "data.frame")
manova_df_protist_phylum
manova_df_protist_phylum_1 = cbind(manova_df_protist_phylum, abundant_protist_phylum)
manova_df_protist_phylum_1

### For Apicomplexa, Rotifera and Cercozoa
manova_abundant_protist_phylum <- manova(cbind(Apicomplexa, Rotifera, Cercozoa) ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_1)
manova_abundant_protist_phylum
summary.aov(manova_abundant_protist_phylum)
eta_squared(aov(manova_abundant_protist_phylum))

## For the most abundant taxa of the protist genus level
manova_df_protist_genus = as(sample_data(sub_protist_genus), "data.frame")
manova_df_protist_genus
manova_df_protist_genus_1 = cbind(manova_df_protist_genus, abundant_protist_genus)
manova_df_protist_genus_1

### For Cephalodella, Cryptosporidium and Bicosoeca  
manova_abundant_protist_genus <- manova(cbind(Cephalodella, Cryptosporidium, Bicosoeca) ~ nutrient*sediment*flow*time, data = manova_df_protist_genus_1)
manova_abundant_protist_genus
summary.aov(manova_abundant_protist_genus)
eta_squared(aov(manova_abundant_protist_genus))



# Generalized linear latent variable models(GLLVMs)
## According to Niku et al., (2019), this statistical analysis was performed.

## For the algae phylum level
sap_otu <- data.frame(otu_table(sub_algae_phylum))
sap_otu

sap_tax <- data.frame(tax_table(sub_algae_phylum))
sap_tax

sap_sample <- data.frame(sample_data(sub_algae_phylum))
sap_sample

sap_tax_otu <- cbind(sap_tax["phylum"], sap_otu)
sap_tax_otu

rownames(sap_tax_otu) <- NULL
sap_tax_otu

sap_tax_otu_uq <- sap_tax_otu %>% distinct(phylum, .keep_all = TRUE)
sap_tax_otu_uq

sap_tax_otu_index <- sap_tax_otu_uq %>%
  remove_rownames() %>%
  column_to_rownames(var = 'phylum')
sap_tax_otu_index

sap_tax_otu_final <- as.data.frame(t(sap_tax_otu_index))
sap_tax_otu_final

## The ordination as model based
mydata_fitp_ap <- gllvm(sap_tax_otu_final, family = poisson(link = "log"))
summary(mydata_fitp_ap)

mydata_fit_ord_nb_ap <- gllvm(sap_tax_otu_final, family = "negative.binomial")
summary(mydata_fit_ord_nb_ap)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_ap, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_ap, var.colors = 1)

ordiplot(mydata_fitp_ap, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_ap, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_ap$params$theta) <- paste("spp", 1:ncol(mydata_fitp_ap$y), sep = "")
ordiplot(mydata_fitp_ap, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_ap_criterias_p <- NULL
for(i in 1:5){
  fiti_ap_p <- gllvm(sap_tax_otu_final, sap_sample, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                      formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_ap_criterias_p[i + 1] <- summary(fiti_ap_p)$AICc
  names(mydata_ap_criterias_p)[i + 1] = i
}
## To compare the AICc values
mydata_ap_criterias_p

mydata_fit_env_p_ap <- gllvm(sap_tax_otu_final, sap_sample, family = poisson(link = "log"), num.lv = 4,
                           formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_ap)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between algal phylums can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_ap_p <- getResidualCor(mydata_fit_env_p_ap)

corrplot(mydata_cr_ap_p[order.single(mydata_cr_ap_p), order.single(mydata_cr_ap_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables and 4 latent variables:
mydata_fit4lv_ap_p <- gllvm(sap_tax_otu_final, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_ap_p)

## The correlation matrix
mydata_cr0_ap_p <- getResidualCor(mydata_fit4lv_ap_p)
corrplot(mydata_cr0_ap_p[order.single(mydata_cr0_ap_p), order.single(mydata_cr0_ap_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual biplot can be used to visualize the residual correlations
mydata_fit_env2_ap_p <- gllvm(sap_tax_otu_final, sap_sample, family = poisson(link = "log"), num.lv = 2, 
                  formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_ap_p)

rownames(mydata_fit_env2_ap_p$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_ap_p$y), sep = "")
ordiplot(mydata_fitp_ap, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_ap_p, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_ap_p <- getResidualCov(mydata_fit_env_p_ap, adjust = 0)
mydata_rcov0_ap_p <- getResidualCov(mydata_fit4lv_ap_p, adjust = 0)
mydata_rcov0_ap_p$trace; mydata_rcov_ap_p$trace

1 - mydata_rcov_ap_p$trace / mydata_rcov0_ap_p$trace


## For the algae genus level
sag_otu <- data.frame(otu_table(sub_algae_genus))
sag_otu

sag_tax <- data.frame(tax_table(sub_algae_genus))
sag_tax

sag_sample <- data.frame(sample_data(sub_algae_genus))
sag_sample

sag_tax_otu <- cbind(sag_tax["genus"], sag_otu)
sag_tax_otu

rownames(sag_tax_otu) <- NULL
sag_tax_otu

sag_tax_otu_uq <- sag_tax_otu %>% distinct(genus, .keep_all = TRUE)
sag_tax_otu_uq

sag_tax_otu_index <- sag_tax_otu_uq %>%
  remove_rownames() %>%
  column_to_rownames(var = 'genus')
sag_tax_otu_index

sag_tax_otu_final <- as.data.frame(t(sag_tax_otu_index))
sag_tax_otu_final

## The ordination as model based
mydata_fitp_ag <- gllvm(sag_tax_otu_final, family = poisson(link = "log"))
summary(mydata_fitp_ag)

mydata_fit_ord_nb_ag <- gllvm(sag_tax_otu_final, family = "negative.binomial")
summary(mydata_fit_ord_nb_ag)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_ag, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_ag, var.colors = 1)

ordiplot(mydata_fitp_ag, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_ag, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_ag$params$theta) <- paste("spp", 1:ncol(mydata_fitp_ag$y), sep = "")
ordiplot(mydata_fitp_ag, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_ag_criterias_p <- NULL
for(i in 1:5){
  fiti_ag_p <- gllvm(sag_tax_otu_final, sag_sample, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_ag_criterias_p[i + 1] <- summary(fiti_ag_p)$AICc
  names(mydata_ag_criterias_p)[i + 1] = i
}
## To compare the AICc values
mydata_ag_criterias_p

mydata_fit_env_p_ag <- gllvm(sag_tax_otu_final, sag_sample, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_ag)

## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between algal genus can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_ag_p <- getResidualCor(mydata_fit_env_p_ag)

corrplot(mydata_cr_ag_p[order.single(mydata_cr_ag_p), order.single(mydata_cr_ag_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables and 4 latent variables:
mydata_fit4lv_ag_p <- gllvm(sag_tax_otu_final, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_ag_p)

## The correlation matrix
mydata_cr0_ag_p <- getResidualCor(mydata_fit4lv_ag_p)
corrplot(mydata_cr0_ag_p[order.single(mydata_cr0_ag_p), order.single(mydata_cr0_ag_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual biplot can be used to visualize the residual correlations
mydata_fit_env2_ag_p <- gllvm(sag_tax_otu_final, sag_sample, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_ag_p)

rownames(mydata_fit_env2_ag_p$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_ag_p$y), sep = "")
ordiplot(mydata_fitp_ag, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_ag_p, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_ag_p <- getResidualCov(mydata_fit_env_p_ag, adjust = 0)
mydata_rcov0_ag_p <- getResidualCov(mydata_fit4lv_ag_p, adjust = 0)
mydata_rcov0_ag_p$trace; mydata_rcov_ag_p$trace

1 - mydata_rcov_ag_p$trace / mydata_rcov0_ag_p$trace

## For the protist phylum level
spp_otu <- data.frame(otu_table(sub_protist_phylum))
spp_otu

spp_tax <- data.frame(tax_table(sub_protist_phylum))
spp_tax

spp_sample <- data.frame(sample_data(sub_protist_phylum))
spp_sample

spp_tax_otu <- cbind(spp_tax["phylum"], spp_otu)
spp_tax_otu

rownames(spp_tax_otu) <- NULL
spp_tax_otu

spp_tax_otu_uq <- spp_tax_otu %>% distinct(phylum, .keep_all = TRUE)
spp_tax_otu_uq

spp_tax_otu_index <- spp_tax_otu_uq %>%
  remove_rownames() %>%
  column_to_rownames(var = 'phylum')
spp_tax_otu_index

spp_tax_otu_final <- as.data.frame(t(spp_tax_otu_index))
spp_tax_otu_final

## The ordination as model based
mydata_fitp_pp <- gllvm(spp_tax_otu_final, family = poisson(link = "log"))
summary(mydata_fitp_pp)

mydata_fit_ord_nb_pp <- gllvm(spp_tax_otu_final, family = "negative.binomial")
summary(mydata_fit_ord_nb_pp)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_pp, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_pp, var.colors = 1)

ordiplot(mydata_fitp_pp, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_pp, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_pp$params$theta) <- paste("spp", 1:ncol(mydata_fitp_pp$y), sep = "")
ordiplot(mydata_fitp_pp, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_pp_criterias_p <- NULL
for(i in 1:5){
  fiti_pp_p <- gllvm(spp_tax_otu_final, spp_sample, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_pp_criterias_p[i + 1] <- summary(fiti_pp_p)$AICc
  names(mydata_pp_criterias_p)[i + 1] = i
}
## To compare the AICc values
mydata_pp_criterias_p

mydata_fit_env_p_pp <- gllvm(spp_tax_otu_final, spp_sample, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_pp)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between protist phylums can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_pp_p <- getResidualCor(mydata_fit_env_p_pp)

corrplot(mydata_cr_pp_p[order.single(mydata_cr_pp_p), order.single(mydata_cr_pp_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables and 4 latent variables:
mydata_fit4lv_pp_p <- gllvm(spp_tax_otu_final, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_pp_p)

## The correlation matrix
mydata_cr0_pp_p <- getResidualCor(mydata_fit4lv_pp_p)
corrplot(mydata_cr0_pp_p[order.single(mydata_cr0_pp_p), order.single(mydata_cr0_pp_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual biplot can be used to visualize the residual correlations
mydata_fit_env2_pp_p <- gllvm(spp_tax_otu_final, spp_sample, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_pp_p)

rownames(mydata_fit_env2_pp_p$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_pp_p$y), sep = "")
ordiplot(mydata_fitp_pp, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_pp_p, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_pp_p <- getResidualCov(mydata_fit_env_p_pp, adjust = 0)
mydata_rcov0_pp_p <- getResidualCov(mydata_fit4lv_pp_p, adjust = 0)
mydata_rcov0_pp_p$trace; mydata_rcov_pp_p$trace

1 - mydata_rcov_pp_p$trace / mydata_rcov0_pp_p$trace

## For the protist genus level
spg_otu <- data.frame(otu_table(sub_protist_genus))
spg_otu

spg_tax <- data.frame(tax_table(sub_protist_genus))
spg_tax

spg_sample <- data.frame(sample_data(sub_protist_genus))
spg_sample

spg_tax_otu <- cbind(spg_tax["genus"], spg_otu)
spg_tax_otu

rownames(spg_tax_otu) <- NULL
spg_tax_otu

spg_tax_otu_uq <- spg_tax_otu %>% distinct(genus, .keep_all = TRUE)
spg_tax_otu_uq

spg_tax_otu_index <- spg_tax_otu_uq %>%
  remove_rownames() %>%
  column_to_rownames(var = 'genus')
spg_tax_otu_index

spg_tax_otu_final <- as.data.frame(t(spg_tax_otu_index))
spg_tax_otu_final

## The ordination as model based
mydata_fitp_pg <- gllvm(spg_tax_otu_final, family = poisson(link = "log"))
summary(mydata_fitp_pg)

mydata_fit_ord_nb_pg <- gllvm(spg_tax_otu_final, family = "negative.binomial")
summary(mydata_fit_ord_nb_pg)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_pg, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_pg, var.colors = 1)

ordiplot(mydata_fitp_pg, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_pg, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_pg$params$theta) <- paste("spp", 1:ncol(mydata_fitp_pg$y), sep = "")
ordiplot(mydata_fitp_pg, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_pg_criterias_p <- NULL
for(i in 1:5){
  fiti_pg_p <- gllvm(spg_tax_otu_final, spg_sample, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_pg_criterias_p[i + 1] <- summary(fiti_pg_p)$AICc
  names(mydata_pg_criterias_p)[i + 1] = i
}
## To compare the AICc values
mydata_pg_criterias_p

mydata_fit_env_p_pg <- gllvm(spg_tax_otu_final, spg_sample, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_pg)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between protist genus can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_pg_p <- getResidualCor(mydata_fit_env_p_pg)

corrplot(mydata_cr_pg_p[order.single(mydata_cr_pg_p), order.single(mydata_cr_pg_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables and 4 latent variables:
mydata_fit4lv_pg_p <- gllvm(spg_tax_otu_final, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_pg_p)

## The correlation matrix
mydata_cr0_pg_p <- getResidualCor(mydata_fit4lv_pg_p)
corrplot(mydata_cr0_pg_p[order.single(mydata_cr0_pg_p), order.single(mydata_cr0_pg_p)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual biplot can be used to visualize the residual correlations
mydata_fit_env2_pg_p <- gllvm(spg_tax_otu_final, spg_sample, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_pg_p)

rownames(mydata_fit_env2_pg_p$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_pg_p$y), sep = "")
ordiplot(mydata_fitp_pg, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_pg_p, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_pg_p <- getResidualCov(mydata_fit_env_p_pg, adjust = 0)
mydata_rcov0_pg_p <- getResidualCov(mydata_fit4lv_pg_p, adjust = 0)
mydata_rcov0_pg_p$trace; mydata_rcov_pg_p$trace

1 - mydata_rcov_pg_p$trace / mydata_rcov0_pg_p$trace







List of references

Vaulot, D. (2021, Feb 15). Phyloseq tutorial. https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

McMurdie, P. J., & Holmes, S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PloS one, 8(4), e61217.

Mariadassou, M., Bernard, M., Pascal, G., Cauquil, L., & Chaillou, S. (2016). Analysis of community composition data using phyloseq. Montpellier Décembre 2016. https://genoweb.toulouse.inra.fr/~formation/15_FROGS/8-February2017/FROGS_phyloseq_02_2017.pdf

Ollberding, N.J. (2019, Jul 28). Introduction to the Statistical Analysis of Microbiome Data in R. Wowchemy. https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

Xu, S., & Yu, G. (2020, Nov 24). Workshop of microbiome dataset analysis using MicrobiotaProcess. MicrobiotaProcessWorkshop 0.0.0.92. https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

Kassambara, A. (2017). MANOVA Test in R: Multivariate Analysis of Variance. Statistical tools for high-throughput data analysis. STHDA. http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance#infos

Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815

ZACH. (2021, Sep 29). How to Test for Normality in R (4 Methods). Statology. https://www.statology.org/test-for-normality-in-r/
  
Humboldt-Universität zu Berlin | Geography Department. (2021). Generalized linear models in R. Quantitative Methods for Geographers. https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab08_GLM1.html

Phillips, N.D. (2018, Jan 22). 15.4 Regression on non-Normal data with glm(). YaRrr! The Pirate's Guide to R. https://bookdown.org/ndphillips/YaRrr/regression-on-non-normal-data-with-glm.html

Bolker, B. (2018, Sep 25). GLMM worked examples. https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html

Niku, J., Hui, F. K. C., Taskinen, S., and Warton, D. I. (2019). gllvm - Fast analysis of multivariate abundance data with generalized linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173--2182. https://github.com/JenniNiku/gllvm

Walker, J.A. (2021, Sep 25). Applied Statistics for Experimental Biology. Elements of Applied Biostatistics. https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/models-with-random-effects-blocking-and-pseudoreplication.html

joey711. (2022). Find the Most abundant Taxa in individual samples #847. https://github.com/joey711/phyloseq/issues/847

  