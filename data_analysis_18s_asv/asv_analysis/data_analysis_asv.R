# Analysis of eukaryotic biodiversity data using the phyloseq package in R

# ASV data analysis with the phyloseq package for statistical analysis

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
asv_mat_asv <- read_excel(here("asv_analysis","dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "ASV matrix")

tax_mat_asv <- read_excel(here("asv_analysis","dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "Taxonomy table")

samples_df_asv <- read_excel(here("asv_analysis","dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "Samples")

## Phyloseq objects must have row.names
## To define the row names from the asv column
asv_mat_asv <- asv_mat_asv %>%
  tibble::column_to_rownames("asv")

## To idem for the two other matrices
## For taxa data
tax_mat_asv <- tax_mat_asv %>% 
  tibble::column_to_rownames("asv")

## For sample data
samples_df_asv <- samples_df_asv %>% 
  tibble::column_to_rownames("sample")

## Transformation to matrices of asv and tax tables (the sample table can be left as a data frame)
asv_mat_asv <- as.matrix(asv_mat_asv)

tax_mat_asv <- as.matrix(tax_mat_asv)

## Transforming to phyloseq objects
ASV_1 = otu_table(asv_mat_asv, taxa_are_rows = TRUE)
TAX_1_asv = tax_table(tax_mat_asv)
samples_1_asv = sample_data(samples_df_asv)

phyloseq_obj_asv <- phyloseq(ASV_1, TAX_1_asv, samples_1_asv)
phyloseq_obj_asv


## Changing order of levels
sample_data(phyloseq_obj_asv)$treatment <- factor(sample_data(phyloseq_obj_asv)$treatment, levels = c("C", "N", "F", "S", "NF", "NS", "SF", "NFS"))
levels(sample_data(phyloseq_obj_asv)$treatment)

## To keep only taxa of interests (for Algae) according to phylum and order level
sub_algae_phylum_asv <- subset_taxa(phyloseq_obj_asv, Phylum %in% c("Chlorophyta_ph", "Cryptomonadales", "Dinoflagellata", "Florideophycidae", "Klebsormidiophyceae", "Ochrophyta", "Phragmoplastophyta", "Prymnesiophyceae"))
sub_algae_order_asv <- subset_taxa(phyloseq_obj_asv, Order %in% c("Bacillariophytina", "Chaetophorales", "Chlamydomonadales", "Chlorellales", "Chlorosarcinales", "Chromulinales", "Coccolithales_or", "Coscinodiscophytina", "Cryptomonadales_or", "Ctenocladales", "Desmidiales", "Eustigmatales", "Gymnodiniphycidae", "Klebsormidiales", "Microthamniales", "Nemaliophycidae_or", "Ochromonadales", "P34.45", "Peridiniphycidae", "Prorocentrales", "Scotinosphaerales", "Sphaeropleales", "Syndiniales_or", "Synurales", "Trebouxiales", "Trentepohliales", "Ulotrichales", "Ulvales", "Zygnematales"))

## To keep only taxa of interests (for Protists) according to phylum and order level
sub_protist_phylum_asv <- subset_taxa(phyloseq_obj_asv, Phylum %in% c("Ancyromonadida", "Apicomplexa", "Cercozoa", "Choanoflagellida", "Ciliophora", "Dictyostelia", "Euglenozoa", "Gracilipodida", "Kathablepharidae", "MAST-12", "MAST-2", "MAST-6", "Myxogastria", "Parabasalia", "Protalveolata", "Protosteliida", "Schizoplasmodiida", "Tubulinea", "LKM15"))
sub_protist_order_asv <- subset_taxa(phyloseq_obj_asv, Order %in% c("A31", "Ancyromonadida_or", "Cercomonadidae_or", "Coccidia", "Conthreep", "Cryomonadida", "Dictyostelia_or", "Echinamoebida", "Gracilipodida_or", "Granofilosea_or", "Kathablepharidae_or", "Litostomatea", "Marimonadida", "Myxogastria_or", "Phytomyxea_or", "Schizoplasmodiida_or", "Spirotrichea", "Stephanoecidae", "Hypocreales"))

# The pre-processing of the data

## For the algae phylum level
## To transform the data from sub_algae_phylum_asv to relative abundance (to transform the count of samples to percent)
sub_algae_phylum_tsc_asv <- transform_sample_counts(sub_algae_phylum_asv, function(x) x / sum(x))
sub_algae_phylum_tsc_asv
## Sum of the values in each column of the asv_table before removing noises on the data
sub_algae_phylum_tsc_a_asv <- as.data.frame(otu_table(sub_algae_phylum_tsc_asv))
sub_algae_phylum_tsc_a_asv
colSums(sub_algae_phylum_tsc_a_asv)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_algae_phylum_tsc_minthreshold_asv  = transform_sample_counts(sub_algae_phylum_asv, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the asv_table after removing noises on the data
sub_algae_phylum_tsc_minthreshold_b_asv <- as.data.frame(otu_table(sub_algae_phylum_tsc_minthreshold_asv))
sub_algae_phylum_tsc_minthreshold_b_asv
colSums(sub_algae_phylum_tsc_minthreshold_b_asv) ## Coverage


## For the algae order level
## To transform the data from sub_algae_order_asv to relative abundance (to transform the count of samples to percent)
sub_algae_order_tsc_asv <- transform_sample_counts(sub_algae_order_asv, function(x) x / sum(x))
sub_algae_order_tsc_asv
## Sum of the values in each column of the asv_table before removing noises on the data
sub_algae_order_tsc_a_asv <- as.data.frame(otu_table(sub_algae_order_tsc_asv))
sub_algae_order_tsc_a_asv
colSums(sub_algae_phylum_tsc_a_asv)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_algae_order_tsc_minthreshold_asv  = transform_sample_counts(sub_algae_order_asv, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the asv_table after removing noises on the data
sub_algae_order_tsc_minthreshold_b_asv <- as.data.frame(otu_table(sub_algae_order_tsc_minthreshold_asv))
sub_algae_order_tsc_minthreshold_b_asv
colSums(sub_algae_order_tsc_minthreshold_b_asv) ## Coverage

## For the protist phylum level
## To transform the data from sub_protist_phylum_asv to relative abundance (to transform the count of samples to percent)
sub_protist_phylum_tsc_asv <- transform_sample_counts(sub_protist_phylum_asv, function(x) x / sum(x))
sub_protist_phylum_tsc_asv

as.data.frame(otu_table(sub_protist_phylum_tsc_asv))

## Sum of the values in each column of the asv_table before removing noises on the data
sub_protist_phylum_tsc_a_asv <- as.data.frame(otu_table(sub_protist_phylum_tsc_asv))
sub_protist_phylum_tsc_a_asv
colSums(sub_protist_phylum_tsc_a_asv)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_protist_phylum_tsc_minthreshold_asv  = transform_sample_counts(sub_protist_phylum_asv, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})

as.data.frame(otu_table(sub_protist_phylum_tsc_minthreshold_asv))

## Sum of the values in each column of the asv_table after removing noises on the data
sub_protist_phylum_tsc_minthreshold_b_asv <- as.data.frame(otu_table(sub_protist_phylum_tsc_minthreshold_asv))
sub_protist_phylum_tsc_minthreshold_b_asv
colSums(sub_protist_phylum_tsc_minthreshold_b_asv) ## Coverage


## For the protist order level
## To transform the data from sub_protist_order_asv to relative abundance (to transform the count of samples to percent)
sub_protist_order_tsc_asv <- transform_sample_counts(sub_protist_order_asv, function(x) x / sum(x))
sub_protist_order_tsc_asv
## Sum of the values in each column of the asv_table before removing noises on the data
sub_protist_order_tsc_a_asv <- as.data.frame(otu_table(sub_protist_order_tsc_asv))
sub_protist_order_tsc_a_asv
colSums(sub_protist_order_tsc_a_asv)
## To transform the readings below 0.003% to 0% for removing noises on the data
sub_protist_order_tsc_minthreshold_asv  = transform_sample_counts(sub_protist_order_asv, function(x, minthreshold=0.003){
  x <- x / sum(x) 
  x[x < minthreshold] <- 0.0
  return(x)
})
## Sum of the values in each column of the asv_table after removing noises on the data
sub_protist_order_tsc_minthreshold_b_asv <- as.data.frame(otu_table(sub_protist_order_tsc_minthreshold_asv))
sub_protist_order_tsc_minthreshold_b_asv
colSums(sub_protist_order_tsc_minthreshold_b_asv) ## Coverage


# Bar graphs of relative abundance

## Basic bar graphs based on phylum and order level for Algae
bar_algae_phylum_asv <- microbiome::transform(sub_algae_phylum_tsc_minthreshold_asv, "compositional")
plot_bar(bar_algae_phylum_asv, fill="Phylum", title = "The relative abundances of the algae at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

bar_algae_order_asv <- microbiome::transform(sub_algae_order_tsc_minthreshold_asv, "compositional")
plot_bar(bar_algae_order_asv, fill="Order", title = "The relative abundances of the algae at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

## Basic bar graphs based on phylum and order level for Protists
bar_protist_phylum_asv <- microbiome::transform(sub_protist_phylum_tsc_minthreshold_asv, "compositional")
plot_bar(bar_protist_phylum_asv, fill="Phylum", title = "The relative abundances of protists at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

bar_protist_order_asv <- microbiome::transform(sub_protist_order_tsc_minthreshold_asv, "compositional")
plot_bar(bar_protist_order_asv, fill="Order", title = "The relative abundances of protists at the level of the order")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")


# Exploring the biodiversity : Alpha-diversity

## According to Mariadassou et al., (2016), this statistical analysis was performed.
## Followed by https://joey711.github.io/phyloseq/plot_richness-examples.html,
## For alpha diversity, it is advisable to prune ASVs that are not present in any of the samples.
## But that's all we can trim. It is tempting to trim noise, but there are many estimates of richness based on the singletons and double tons of the abundance data.
## To get a meaningful estimate, we need to leave it in the data set.

## For the algae phylum level
## To prune taxa
alph_alg_sset_phy_1_asv <- prune_taxa(taxa_sums(sub_algae_phylum_asv) > 0, sub_algae_phylum_asv)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_algae_phylum_asv <- estimate_richness(alph_alg_sset_phy_1_asv, measures = c("Observed", "Shannon"))
alpha.diversity_algae_phylum_asv
## To visualize the estimated Observed Richness and Shannon Diversity
p_ap_asv = plot_richness(alph_alg_sset_phy_1_asv, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the algal phylum level")
p_ap_asv + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness             
pielou_algae_phylum_asv <- evenness(alph_alg_sset_phy_1_asv, 'pielou')
pielou_algae_phylum_asv
## To arrange the sample data in a data frame
dat_lm_algae_phylum_asv <- data.frame(sample_data(alph_alg_sset_phy_1_asv))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_algae_phylum_asv <- cbind(alpha.diversity_algae_phylum_asv, pielou_algae_phylum_asv, dat_lm_algae_phylum_asv)
alpha_lm_algae_phylum_asv
alpha_lm_algae_phylum_1_asv <- rownames_to_column(alpha_lm_algae_phylum_asv, var = "samples")
alpha_lm_algae_phylum_1_asv
## To visualize the Pielous Evenness
ggplot(alpha_lm_algae_phylum_1_asv, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the algal phylum level")
## To mutate the variables as factor
str(alpha_lm_algae_phylum_1_asv)
alpha_lm_algae_phylum_11_asv <- alpha_lm_algae_phylum_1_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(alpha_lm_algae_phylum_11_asv)

## For the algae order level
## To prune taxa
alph_alg_sset_ord_1_asv <- prune_taxa(taxa_sums(sub_algae_order_asv) > 0, sub_algae_order_asv)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_algae_order_asv <- estimate_richness(alph_alg_sset_ord_1_asv, measures = c("Observed", "Shannon"))
alpha.diversity_algae_order_asv
## To visualize the estimated Observed Richness and Shannon Diversity
p_ao_asv = plot_richness(alph_alg_sset_ord_1_asv, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the algal order level")
p_ao_asv + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness 
pielou_algae_order_asv <- evenness(alph_alg_sset_ord_1_asv, 'pielou')
pielou_algae_order_asv
## To arrange the sample data in a data frame
dat_lm_algae_order_asv <- data.frame(sample_data(alph_alg_sset_ord_1_asv))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_algae_order_asv <- cbind(alpha.diversity_algae_order_asv, pielou_algae_order_asv, dat_lm_algae_order_asv)
alpha_lm_algae_order_asv
alpha_lm_algae_order_1_asv <- rownames_to_column(alpha_lm_algae_order_asv, var = "samples")
alpha_lm_algae_order_1_asv
## To visualize the Pielous Evenness
ggplot(alpha_lm_algae_order_1_asv, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the algal order level")
## To mutate the variables as factor
str(alpha_lm_algae_order_1_asv)
alpha_lm_algae_order_11_asv <- alpha_lm_algae_order_1_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(alpha_lm_algae_order_11_asv)

## For the protist phylum level
## To prune taxa
alph_prot_sset_phy_1_asv <- prune_taxa(taxa_sums(sub_protist_phylum_asv) > 0, sub_protist_phylum_asv)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_protist_phylum_asv <- estimate_richness(alph_prot_sset_phy_1_asv, measures = c("Observed", "Shannon"))
alpha.diversity_protist_phylum_asv
## To visualize the estimated Observed Richness and Shannon Diversity
p_pp_asv = plot_richness(alph_prot_sset_phy_1_asv, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the protist phylum level")
p_pp_asv + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness
pielou_protist_phylum_asv <- evenness(alph_prot_sset_phy_1_asv, 'pielou')
pielou_protist_phylum_asv
pielou_protist_phylum_asv$pielou[is.na(pielou_protist_phylum_asv$pielou)] <- mean(pielou_protist_phylum_asv$pielou, na.rm = TRUE)
pielou_protist_phylum_asv
## To arrange the sample data in a data frame
dat_lm_protist_phylum_asv <- data.frame(sample_data(alph_prot_sset_phy_1_asv))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_protist_phylum_asv <- cbind(alpha.diversity_protist_phylum_asv, pielou_protist_phylum_asv, dat_lm_protist_phylum_asv)
alpha_lm_protist_phylum_asv
alpha_lm_protist_phylum_1_asv <- rownames_to_column(alpha_lm_protist_phylum_asv, var = "samples")
alpha_lm_protist_phylum_1_asv
## To visualize the Pielous Evenness
ggplot(alpha_lm_protist_phylum_1_asv, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the protist phylum level")
## To mutate the variables as factor
str(alpha_lm_protist_phylum_1_asv)
alpha_lm_protist_phylum_11_asv <- alpha_lm_protist_phylum_1_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(alpha_lm_protist_phylum_11_asv)

## For the protist order level
## To prune taxa
alph_prot_sset_ord_1_asv <- prune_taxa(taxa_sums(sub_protist_order_asv) > 0, sub_protist_order_asv)
## To estimate the Observed Richness and Shannon Diversity
alpha.diversity_protist_order_asv <- estimate_richness(alph_prot_sset_ord_1_asv, measures = c("Observed", "Shannon"))
alpha.diversity_protist_order_asv
## To visualize the estimated Observed Richness and Shannon Diversity
p_po_asv = plot_richness(alph_prot_sset_ord_1_asv, x="samples", color="treatment", measures=c("Observed", "Shannon"), title = "Alpha diversity measures for the protist order level")
p_po_asv + geom_point(size=5, alpha=0.7)
## To estimate Pielous Evenness
pielou_protist_order_asv <- evenness(alph_prot_sset_ord_1_asv, 'pielou')
pielou_protist_order_asv
pielou_protist_order_asv$pielou[is.na(pielou_protist_order_asv$pielou)] <- mean(pielou_protist_order_asv$pielou, na.rm = TRUE)
pielou_protist_order_asv
## To arrange the sample data in a data frame
dat_lm_protist_order_asv <- data.frame(sample_data(alph_prot_sset_ord_1_asv))
## To combine the Observed Richness, Shannon Diversity, Pielous Evenness estimates and sample data in a data frame
alpha_lm_protist_order_asv <- cbind(alpha.diversity_protist_order_asv, pielou_protist_order_asv, dat_lm_protist_order_asv)
alpha_lm_protist_order_asv
alpha_lm_protist_order_1_asv <- rownames_to_column(alpha_lm_protist_order_asv, var = "samples")
alpha_lm_protist_order_1_asv
## To visualize the Pielous Evenness
ggplot(alpha_lm_protist_order_1_asv, aes(x=samples, y=pielou, color=treatment)) + 
  geom_boxplot()+
  ggtitle("Pielous evenness for the protist order level")
## To mutate the variables as factor
str(alpha_lm_protist_order_1_asv)
alpha_lm_protist_order_11_asv <- alpha_lm_protist_order_1_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(alpha_lm_protist_order_11_asv)


# Modeling
## According to Zuur et al., (2009), this statistical analysis was performed.

# For the algae phylum level

## For Observed
### To check the normality of the data
hist((alpha_lm_algae_phylum_11_asv$Observed))
boxplot(Observed ~ treatment, data = alpha_lm_algae_phylum_11_asv, main= "Observed richness for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11_asv$Observed)
qqline(alpha_lm_algae_phylum_11_asv$Observed)
shapiro.test(alpha_lm_algae_phylum_11_asv$Observed)
### To fit a linear model with the gls function
ap_ob.gls_asv <- gls(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(ap_ob.gls_asv)
summary(ap_ob.gls_asv)
### To fit a linear mixed model
lmm_algae_phylum_Observed_asv <- lme(Observed ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11_asv)
plot(lmm_algae_phylum_Observed_asv)
summary(lmm_algae_phylum_Observed_asv)
### To compare the models
anova(ap_ob.gls_asv, lmm_algae_phylum_Observed_asv)
### To fit the final model
lm_algae_phylum_Observed_asv <- lm(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(lm_algae_phylum_Observed_asv)
leveneTest(lm_algae_phylum_Observed_asv)
summary(lm_algae_phylum_Observed_asv)
anova(lm_algae_phylum_Observed_asv)
eta_squared(lm_algae_phylum_Observed_asv)
### Bar plots of ANOVA results for algal phylum level for Observed Richness
plot_ao_ap_ob_asv <- ggplot(alpha_lm_algae_phylum_11_asv, aes(x = treatment, y = Observed, fill = treatment))+
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
plot_ao_ap_ob_asv
### Significant interaction plots of algal phylum level for Observed Richness
apoi_t2_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 2)
apoi_t3_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(apoi_t2_asv$nutrient,apoi_t2_asv$flow,apoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal phylum level for nutrient and flow for observed richness")
interaction.plot(apoi_t3_asv$nutrient,apoi_t3_asv$flow,apoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal phylum level for nutrient and flow for observed richness")

par(mfrow=c(1,1))
interaction.plot(apoi_t2_asv$nutrient,apoi_t2_asv$sediment,apoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal phylum level for nutrient and sediment for observed richness")
interaction.plot(apoi_t3_asv$nutrient,apoi_t3_asv$sediment,apoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal phylum level for nutrient and sediment for observed richness")

par(mfrow=c(1,1))
interaction.plot(apoi_t2_asv$sediment,apoi_t2_asv$flow,apoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal phylum level for sediment and flow for observed richness")
interaction.plot(apoi_t3_asv$sediment,apoi_t3_asv$flow,apoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal phylum level for sediment and flow for observed richness")

## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_phylum_11_asv$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_algae_phylum_11_asv, main= "Shannon diversity for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11_asv$Shannon)
qqline(alpha_lm_algae_phylum_11_asv$Shannon)
shapiro.test(alpha_lm_algae_phylum_11_asv$Shannon)
### To fit a linear model with the gls function
ap_sh.gls_asv <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(ap_sh.gls_asv)
summary(ap_sh.gls_asv)
### To fit a linear mixed model
lmm_algae_phylum_Shannon_asv <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11_asv)
plot(lmm_algae_phylum_Shannon_asv)
summary(lmm_algae_phylum_Shannon_asv)
### To compare the models
anova(ap_sh.gls_asv, lmm_algae_phylum_Shannon_asv)
### To fit the final model
lm_algae_phylum_Shannon_asv <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(lm_algae_phylum_Shannon_asv)
leveneTest(lm_algae_phylum_Shannon_asv)
summary(lm_algae_phylum_Shannon_asv)
anova(lm_algae_phylum_Shannon_asv)
eta_squared(lm_algae_phylum_Shannon_asv)
### Bar plots of ANOVA results for algal phylum level for Shannon Diversity
plot_ao_ap_sh_asv <- ggplot(alpha_lm_algae_phylum_11_asv, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
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
plot_ao_ap_sh_asv
### Significant interaction plots of algal phylum level for Shannon Diversity
apsi_t2_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 2)
apsi_t3_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(apsi_t2_asv$nutrient,apsi_t2_asv$flow,apsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal phylum level for nutrient and flow for Shannon diversity")
interaction.plot(apsi_t3_asv$nutrient,apsi_t3_asv$flow,apsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal phylum level for nutrient and flow for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(apsi_t2_asv$nutrient,apsi_t2_asv$sediment,apsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal phylum level for nutrient and sediment for Shannon diversity")
interaction.plot(apsi_t3_asv$nutrient,apsi_t3_asv$sediment,apsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal phylum level for nutrient and sediment for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(apsi_t2_asv$sediment,apsi_t2_asv$flow,apsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal phylum level for sediment and flow for Shannon diversity")
interaction.plot(apsi_t3_asv$sediment,apsi_t3_asv$flow,apsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal phylum level for sediment and flow for Shannon diversity")

## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_phylum_11_asv$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_algae_phylum_11_asv, main= "Pielous evenness for the algae phylum level")
qqnorm(alpha_lm_algae_phylum_11_asv$pielou)
qqline(alpha_lm_algae_phylum_11_asv$pielou)
shapiro.test(alpha_lm_algae_phylum_11_asv$pielou)
### To fit a linear model with the gls function
ap_pi.gls_asv <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(ap_pi.gls_asv)
summary(ap_pi.gls_asv)
### To fit a linear mixed model
lmm_algae_phylum_pielou_asv <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_phylum_11_asv)
plot(lmm_algae_phylum_pielou_asv)
summary(lmm_algae_phylum_pielou_asv)
### To compare the models
anova(ap_pi.gls_asv, lmm_algae_phylum_pielou_asv)
### To fit the final model
lm_algae_phylum_pielou_asv <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_phylum_11_asv)
plot(lm_algae_phylum_pielou_asv)
leveneTest(lm_algae_phylum_pielou_asv)
summary(lm_algae_phylum_pielou_asv)
anova(lm_algae_phylum_pielou_asv)
eta_squared(lm_algae_phylum_pielou_asv)
### Bar plots of ANOVA results for algal phylum level for Pielou's Evenness
plot_ao_ap_pi_asv <- ggplot(alpha_lm_algae_phylum_11_asv, aes(x = treatment, y = log(10+pielou), fill = treatment))+
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
plot_ao_ap_pi_asv
### Significant interaction plots of algal phylum level for Pielou's Evenness
appi_t2_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 2)
appi_t3_asv <- alpha_lm_algae_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(appi_t2_asv$nutrient,appi_t2_asv$flow,appi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal phylum level for nutrient and flow for Pielous evenness")
interaction.plot(appi_t3_asv$nutrient,appi_t3_asv$flow,appi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal phylum level for nutrient and flow for Pielous evenness")

par(mfrow=c(1,1))
interaction.plot(appi_t2_asv$nutrient,appi_t2_asv$sediment,appi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal phylum level for nutrient and sediment for Pielous evenness")
interaction.plot(appi_t3_asv$nutrient,appi_t3_asv$sediment,appi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal phylum level for nutrient and sediment for Pielous evenness")

par(mfrow=c(1,1))
interaction.plot(appi_t2_asv$sediment,appi_t2_asv$flow,appi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal phylum level for sediment and flow for Pielous evenness")
interaction.plot(appi_t3_asv$sediment,appi_t3_asv$flow,appi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal phylum level for sediment and flow for Pielous evenness")


# For the algae order level

## For Observed
### To check the normality of the data
hist(alpha_lm_algae_order_11_asv$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_algae_order_11_asv, main= "Observed richness for the algae order level")
qqnorm(alpha_lm_algae_order_11_asv$Observed)
qqline(alpha_lm_algae_order_11_asv$Observed)
shapiro.test(alpha_lm_algae_order_11_asv$Observed)
### To fit a linear model with the gls function
ao_ob.gls_asv <- gls(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(ao_ob.gls_asv)
summary(ao_ob.gls_asv)
### To fit a linear mixed model
lmm_algae_order_Observed_asv <- lme(Observed ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11_asv)
plot(lmm_algae_order_Observed_asv)
summary(lmm_algae_order_Observed_asv)
### To compare the models
anova(ao_ob.gls_asv, lmm_algae_order_Observed_asv)
### To fit the final model
lm_algae_order_Observed_asv <- lm(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(lm_algae_order_Observed_asv)
leveneTest(lm_algae_order_Observed_asv)
summary(lm_algae_order_Observed_asv)
anova(lm_algae_order_Observed_asv)
eta_squared(lm_algae_order_Observed_asv)
### Bar plots of ANOVA results for algal order level for Observed Richness
plot_ao_ao_ob_asv <- ggplot(alpha_lm_algae_order_11_asv, aes(x = treatment, y = Observed, fill = treatment))+
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
plot_ao_ao_ob_asv
### Significant interaction plots of algal order level for Observed Richness
aooi_t2_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 2)
aooi_t3_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(aooi_t2_asv$nutrient,aooi_t2_asv$flow,aooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal order level for nutrient and flow for observed richness")
interaction.plot(aooi_t3_asv$nutrient,aooi_t3_asv$flow,aooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal order level for nutrient and flow for observed richness")

par(mfrow=c(1,1))
interaction.plot(aooi_t2_asv$nutrient,aooi_t2_asv$sediment,aooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal order level for nutrient and sediment for observed richness")
interaction.plot(aooi_t3_asv$nutrient,aooi_t3_asv$sediment,aooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal order level for nutrient and sediment for observed richness")

par(mfrow=c(1,1))
interaction.plot(aooi_t2_asv$sediment,aooi_t2_asv$flow,aooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for algal order level for sediment and flow for observed richness")
interaction.plot(aooi_t3_asv$sediment,aooi_t3_asv$flow,aooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for algal order level for sediment and flow for observed richness")

## For Shannon
### To check the normality of the data
hist(alpha_lm_algae_order_11_asv$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_algae_order_11_asv, main= "Shannon diversity for the algae order level")
qqnorm(alpha_lm_algae_order_11_asv$Shannon)
qqline(alpha_lm_algae_order_11_asv$Shannon)
shapiro.test(alpha_lm_algae_order_11_asv$Shannon)
### To fit a linear model with the gls function
ao_sh.gls_asv <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(ao_sh.gls_asv)
summary(ao_sh.gls_asv)
### To fit a linear mixed model
lmm_algae_order_Shannon_asv <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11_asv)
plot(lmm_algae_order_Shannon_asv)
summary(lmm_algae_order_Shannon_asv)
### To compare the models
anova(ao_sh.gls_asv, lmm_algae_order_Shannon_asv)
### To fit the final model
lm_algae_order_Shannon_asv <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(lm_algae_order_Shannon_asv)
leveneTest(lm_algae_order_Shannon_asv)
summary(lm_algae_order_Shannon_asv)
anova(lm_algae_order_Shannon_asv)
eta_squared(lm_algae_order_Shannon_asv)
### Bar plots of ANOVA results for algal order level for Shannon Diversity
plot_ao_ao_sh_asv <- ggplot(alpha_lm_algae_order_11_asv, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
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
plot_ao_ao_sh_asv
### Significant interaction plots of algal order level for Shannon Diversity
aosi_t2_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 2)
aosi_t3_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(aosi_t2_asv$nutrient,aosi_t2_asv$flow,aosi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal order level for nutrient and flow for Shannon diversity")
interaction.plot(aosi_t3_asv$nutrient,aosi_t3_asv$flow,aosi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal order level for nutrient and flow for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(aosi_t2_asv$nutrient,aosi_t2_asv$sediment,aosi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal order level for nutrient and sediment for Shannon diversity")
interaction.plot(aosi_t3_asv$nutrient,aosi_t3_asv$sediment,aosi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal order level for nutrient and sediment for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(aosi_t2_asv$sediment,aosi_t2_asv$flow,aosi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for algal order level for sediment and flow for Shannon Diversity")
interaction.plot(aosi_t3_asv$sediment,aosi_t3_asv$flow,aosi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for algal order level for sediment and flow for Shannon diversity")

## For Pielou
### To check the normality of the data
hist(alpha_lm_algae_order_11_asv$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_algae_order_11_asv, main= "Pielous evenness for the algae order level")
qqnorm(alpha_lm_algae_order_11_asv$pielou)
qqline(alpha_lm_algae_order_11_asv$pielou)
shapiro.test(alpha_lm_algae_order_11_asv$pielou)
### To fit a linear model with the gls function
ao_pi.gls_asv <- gls(pielou ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(ao_pi.gls_asv)
summary(ao_pi.gls_asv)
### To fit a linear mixed model
lmm_algae_order_pielou_asv <- lme(pielou ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_algae_order_11_asv)
plot(lmm_algae_order_pielou_asv)
summary(lmm_algae_order_pielou_asv)
### To compare the models
anova(ao_pi.gls_asv, lmm_algae_order_pielou_asv)
### To fit the final model
lm_algae_order_pielou_asv <- lm(pielou ~ nutrient*sediment*flow*time, data = alpha_lm_algae_order_11_asv)
plot(lm_algae_order_pielou_asv)
leveneTest(lm_algae_order_pielou_asv)
summary(lm_algae_order_pielou_asv)
anova(lm_algae_order_pielou_asv)
eta_squared(lm_algae_order_pielou_asv)
### Bar plots of ANOVA results for algal order level for Pielou's Evenness
plot_ao_ao_pi_asv <- ggplot(alpha_lm_algae_order_11_asv, aes(x = treatment, y = pielou, fill = treatment))+
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
plot_ao_ao_pi_asv
### Significant interaction plots of algal order level for Pielou's Evenness
aopi_t2_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 2)
aopi_t3_asv <- alpha_lm_algae_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(aopi_t2_asv$nutrient,aopi_t2_asv$flow,aopi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal order level for nutrient and flow for Pielous evenness")
interaction.plot(aopi_t3_asv$nutrient,aopi_t3_asv$flow,aopi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal order level for nutrient and flow for Pielous evenness")

par(mfrow=c(1,1))
interaction.plot(aopi_t2_asv$nutrient,aopi_t2_asv$sediment,aopi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal order level for nutrient and sediment for Pielous evenness")
interaction.plot(aopi_t3_asv$nutrient,aopi_t3_asv$sediment,aopi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal order level for nutrient and sediment for Pielous evenness")

par(mfrow=c(1,1))
interaction.plot(aopi_t2_asv$sediment,aopi_t2_asv$flow,aopi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for algal order level for sediment and flow for Pielous evenness")
interaction.plot(aopi_t3_asv$sediment,aopi_t3_asv$flow,aopi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for algal order level for sediment and flow for Pielous evenness")

# For the protist phylum level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_phylum_11_asv$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_protist_phylum_11_asv, main= "Observed richness for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11_asv$Observed)
qqline(alpha_lm_protist_phylum_11_asv$Observed)
shapiro.test(alpha_lm_protist_phylum_11_asv$Observed)
### To fit a linear nmodel with the gls function
pp_ob.gls_asv <- gls(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(pp_ob.gls_asv)
summary(pp_ob.gls_asv)
### To fit a linear mixed model
lmm_protist_phylum_Observed_asv <- lme(Observed ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11_asv)
plot(lmm_protist_phylum_Observed_asv)
summary(lmm_protist_phylum_Observed_asv)
### To compare the models
anova(pp_ob.gls_asv, lmm_protist_phylum_Observed_asv)
### To fit the final model
lm_protist_phylum_Observed_asv <- lm(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(lm_protist_phylum_Observed_asv)
leveneTest(lm_protist_phylum_Observed_asv)
summary(lm_protist_phylum_Observed_asv)
anova(lm_protist_phylum_Observed_asv)
eta_squared(lm_protist_phylum_Observed_asv)
### Bar plots of ANOVA results for protist phylum level for Observed Richness
plot_ao_pp_ob_asv <- ggplot(alpha_lm_protist_phylum_11_asv, aes(x = treatment, y = Observed, fill = treatment))+
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
plot_ao_pp_ob_asv
### Significant interaction plots of protist phylum level for Observed Richness
ppoi_t2_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 2)
ppoi_t3_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(ppoi_t2_asv$nutrient,ppoi_t2_asv$flow,ppoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist phylum level for nutrient and flow for Observed richness")
interaction.plot(ppoi_t3_asv$nutrient,ppoi_t3_asv$flow,ppoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist phylum level for nutrient and flow for Observed richness")

par(mfrow=c(1,1))
interaction.plot(ppoi_t2_asv$nutrient,ppoi_t2_asv$sediment,ppoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist phylum level for nutrient and sediment for Observed richness")
interaction.plot(ppoi_t3_asv$nutrient,ppoi_t3_asv$sediment,ppoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist phylum level for nutrient and sediment for Observed richness")

par(mfrow=c(1,1))
interaction.plot(ppoi_t2_asv$sediment,ppoi_t2_asv$flow,ppoi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist phylum level for sediment and flow for Observed richness")
interaction.plot(ppoi_t3_asv$sediment,ppoi_t3_asv$flow,ppoi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist phylum level for sediment and flow for Observed richness")

## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_phylum_11_asv$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_protist_phylum_11_asv, main= "Shannon diversity for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11_asv$Shannon)
qqline(alpha_lm_protist_phylum_11_asv$Shannon)
shapiro.test(alpha_lm_protist_phylum_11_asv$Shannon)
### To fit a linear model with the gls function
pp_sh.gls_asv <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(pp_sh.gls_asv)
summary(pp_sh.gls_asv)
### To fit a linear mixed model
lmm_protist_phylum_Shannon_asv <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11_asv)
plot(lmm_protist_phylum_Shannon_asv)
summary(lmm_protist_phylum_Shannon_asv)
### To compare the models
anova(pp_sh.gls_asv, lmm_protist_phylum_Shannon_asv)
### To fit the final model
lm_protist_phylum_Shannon_asv <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(lm_protist_phylum_Shannon_asv)
leveneTest(lm_protist_phylum_Shannon_asv)
summary(lm_protist_phylum_Shannon_asv)
anova(lm_protist_phylum_Shannon_asv)
eta_squared(lm_protist_phylum_Shannon_asv)
### Bar plots of ANOVA results for protist phylum level for Shannon Diversity
plot_ao_pp_sh_asv <- ggplot(alpha_lm_protist_phylum_11_asv, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
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
plot_ao_pp_sh_asv
### Significant interaction plots of protist phylum level for Shannon Diversity
ppsi_t2_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 2)
ppsi_t3_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(ppsi_t2_asv$nutrient,ppsi_t2_asv$flow,ppsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist phylum level for nutrient and flow for Shannon diversity")
interaction.plot(ppsi_t3_asv$nutrient,ppsi_t3_asv$flow,ppsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist phylum level for nutrient and flow for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(ppsi_t2_asv$nutrient,ppsi_t2_asv$sediment,ppsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist phylum level for nutrient and sediment for Shannon diversity")
interaction.plot(ppsi_t3_asv$nutrient,ppsi_t3_asv$sediment,ppsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist phylum level for nutrient and sediment for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(ppsi_t2_asv$sediment,ppsi_t2_asv$flow,ppsi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist phylum level for sediment and flow for Shannon diversity")
interaction.plot(ppsi_t3_asv$sediment,ppsi_t3_asv$flow,ppsi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist phylum level for sediment and flow for Shannon diversity")

## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_phylum_11_asv$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_protist_phylum_11_asv, main= "Pielous evenness for the protist phylum level")
qqnorm(alpha_lm_protist_phylum_11_asv$pielou)
qqline(alpha_lm_protist_phylum_11_asv$pielou)
shapiro.test(alpha_lm_protist_phylum_11_asv$pielou)
### To fit a linear model with the gls function
pp_pi.gls_asv <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(pp_pi.gls_asv)
summary(pp_pi.gls_asv)
### To fit a linear mixed model
lmm_protist_phylum_pielou_asv <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_phylum_11_asv)
plot(lmm_protist_phylum_pielou_asv)
summary(lmm_protist_phylum_pielou_asv)
### To compare the models
anova(pp_pi.gls_asv, lmm_protist_phylum_pielou_asv)
### To fit the final model
lm_protist_phylum_pielou_asv <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_phylum_11_asv)
plot(lm_protist_phylum_pielou_asv)
leveneTest(lm_protist_phylum_pielou_asv)
summary(lm_protist_phylum_pielou_asv)
anova(lm_protist_phylum_pielou_asv)
eta_squared(lm_protist_phylum_pielou_asv)
### Bar plots of ANOVA results for protist phylum level for Pielou's Evenness
plot_ao_pp_pi_asv <- ggplot(alpha_lm_protist_phylum_11_asv, aes(x = treatment, y = log(10+pielou), fill = treatment))+
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
plot_ao_pp_pi_asv
### Significant interaction plots of protist phylum level for Pielou's Evenness
pppi_t2_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 2)
pppi_t3_asv <- alpha_lm_protist_phylum_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(pppi_t2_asv$nutrient,pppi_t2_asv$flow,pppi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist phylum level for nutrient and flow for Pielous Evenness")
interaction.plot(pppi_t3_asv$nutrient,pppi_t3_asv$flow,pppi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist phylum level for nutrient and flow for Pielous Evenness")

par(mfrow=c(1,1))
interaction.plot(pppi_t2_asv$nutrient,pppi_t2_asv$sediment,pppi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist phylum level for nutrient and sediment for Pielous Evenness")
interaction.plot(pppi_t3_asv$nutrient,pppi_t3_asv$sediment,pppi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist phylum level for nutrient and sediment for Pielous Evenness")

par(mfrow=c(1,1))
interaction.plot(pppi_t2_asv$sediment,pppi_t2_asv$flow,pppi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist phylum level for sediment and flow for Pielous Evenness")
interaction.plot(pppi_t3_asv$sediment,pppi_t3_asv$flow,pppi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist phylum level for sediment and flow for Pielous Evenness")

# For the protist order level

## For Observed
### To check the normality of the data
hist(alpha_lm_protist_order_11_asv$Observed)
boxplot(Observed ~ treatment, data = alpha_lm_protist_order_11_asv, main= "Observed richness for the protist order level")
qqnorm(alpha_lm_protist_order_11_asv$Observed)
qqline(alpha_lm_protist_order_11_asv$Observed)
shapiro.test(alpha_lm_protist_order_11_asv$Observed)
### To fit a linear model with the gls function
po_ob.gls_asv <- gls(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(po_ob.gls_asv)
summary(po_ob.gls_asv)
### To fit a linear mixed model
lmm_protist_order_Observed_asv <- lme(Observed ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11_asv)
plot(lmm_protist_order_Observed_asv)
summary(lmm_protist_order_Observed_asv)
### To compare the models
anova(po_ob.gls_asv, lmm_protist_order_Observed_asv)
### To fit the final model
lm_protist_order_Observed_asv <- lm(Observed ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(lm_protist_order_Observed_asv)
leveneTest(lm_protist_order_Observed_asv)
summary(lm_protist_order_Observed_asv)
anova(lm_protist_order_Observed_asv)
eta_squared(lm_protist_order_Observed_asv)
### Bar plots of ANOVA results for protist order level for Observed Richness
plot_ao_po_ob_asv <- ggplot(alpha_lm_protist_order_11_asv, aes(x = treatment, y = Observed, fill = treatment))+
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
plot_ao_po_ob_asv
### Significant interaction plots of protist order level for Observed Richness
pooi_t2_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 2)
pooi_t3_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(pooi_t2_asv$nutrient,pooi_t2_asv$flow,pooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist order level for nutrient and flow for Observed richness")
interaction.plot(pooi_t3_asv$nutrient,pooi_t3_asv$flow,pooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist order level for nutrient and flow for Observed richness")

par(mfrow=c(1,1))
interaction.plot(pooi_t2_asv$nutrient,pooi_t2_asv$sediment,pooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist order level for nutrient and sediment for Observed richness")
interaction.plot(pooi_t3_asv$nutrient,pooi_t3_asv$sediment,pooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist order level for nutrient and sediment for Observed richness")

par(mfrow=c(1,1))
interaction.plot(pooi_t2_asv$sediment,pooi_t2_asv$flow,pooi_t2_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 2 interaction plot for protist order level for sediment and flow for Observed richness")
interaction.plot(pooi_t3_asv$sediment,pooi_t3_asv$flow,pooi_t3_asv$Observed,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Observed Richness",main="Week 3 interaction plot for protist order level for sediment and flow for Observed richness")

## For Shannon
### To check the normality of the data
hist(alpha_lm_protist_order_11_asv$Shannon)
boxplot(Shannon ~ treatment, data = alpha_lm_protist_order_11_asv, main= "Shannon diversity for the protist order level")
qqnorm(alpha_lm_protist_order_11_asv$Shannon)
qqline(alpha_lm_protist_order_11_asv$Shannon)
shapiro.test(alpha_lm_protist_order_11_asv$Shannon)
### To fit a linear model with the gls function
po_sh.gls_asv <- gls(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(po_sh.gls_asv)
summary(po_sh.gls_asv)
### To fit a linear mixed model
lmm_protist_order_Shannon_asv <- lme(log(10+Shannon) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11_asv)
plot(lmm_protist_order_Shannon_asv)
summary(lmm_protist_order_Shannon_asv)
### To compare the models
anova(po_sh.gls_asv, lmm_protist_order_Shannon_asv)
### To fit the final model
lm_protist_order_Shannon_asv <- lm(log(10+Shannon) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(lm_protist_order_Shannon_asv)
leveneTest(lm_protist_order_Shannon_asv)
summary(lm_protist_order_Shannon_asv)
anova(lm_protist_order_Shannon_asv)
eta_squared(lm_protist_order_Shannon_asv)
### Bar plots of ANOVA results for protist order level for Shannon Diversity
plot_ao_po_sh_asv <- ggplot(alpha_lm_protist_order_11_asv, aes(x = treatment, y = log(10+Shannon), fill = treatment))+
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
plot_ao_po_sh_asv
### Significant interaction plots of protist order level for Shannon Diversity
posi_t2_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 2)
posi_t3_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(posi_t2_asv$nutrient,posi_t2_asv$flow,posi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist order level for nutrient and flow for Shannon diversity")
interaction.plot(posi_t3_asv$nutrient,posi_t3_asv$flow,posi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist order level for nutrient and flow for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(posi_t2_asv$nutrient,posi_t2_asv$sediment,posi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist order level for nutrient and sediment for Shannon diversity")
interaction.plot(posi_t3_asv$nutrient,posi_t3_asv$sediment,posi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist order level for nutrient and sediment for Shannon diversity")

par(mfrow=c(1,1))
interaction.plot(posi_t2_asv$sediment,posi_t2_asv$flow,posi_t2_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 2 interaction plot for protist order level for sediment and flow for Shannon diversity")
interaction.plot(posi_t3_asv$sediment,posi_t3_asv$flow,posi_t3_asv$Shannon,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Shannon Diversity",main="Week 3 interaction plot for protist order level for sediment and flow for Shannon diversity")

## For Pielou
### To check the normality of the data
hist(alpha_lm_protist_order_11_asv$pielou)
boxplot(pielou ~ treatment, data = alpha_lm_protist_order_11_asv, main= "Pielous evenness for the protist order level")
qqnorm(alpha_lm_protist_order_11_asv$pielou)
qqline(alpha_lm_protist_order_11_asv$pielou)
shapiro.test(alpha_lm_protist_order_11_asv$pielou)
### To fit a linear model with the gls function
po_pi.gls_asv <- gls(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(po_pi.gls_asv)
summary(po_pi.gls_asv)
### To fit a linear mixed model
lmm_protist_order_pielou_asv <- lme(log(10+pielou) ~ nutrient*sediment*flow*time, random = ~ 1 | header_tanks_block, method = "REML",  data = alpha_lm_protist_order_11_asv)
plot(lmm_protist_order_pielou_asv)
summary(lmm_protist_order_pielou_asv)
### To compare the models
anova(po_pi.gls_asv, lmm_protist_order_pielou_asv)
### To fit the final model
lm_protist_order_pielou_asv <- lm(log(10+pielou) ~ nutrient*sediment*flow*time, data = alpha_lm_protist_order_11_asv)
plot(lm_protist_order_pielou_asv)
leveneTest(lm_protist_order_pielou_asv)
summary(lm_protist_order_pielou_asv)
anova(lm_protist_order_pielou_asv)
eta_squared(lm_protist_order_pielou_asv)
### Bar plots of ANOVA results for protist order level for Pielou's Evenness
plot_ao_po_pi_asv <- ggplot(alpha_lm_protist_order_11_asv, aes(x = treatment, y = log(10+pielou), fill = treatment))+
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
plot_ao_po_pi_asv
### Significant interaction plots of protist order level for Pielou's Evenness
popi_t2_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 2)
popi_t3_asv <- alpha_lm_protist_order_11_asv %>% filter(week == 3)

par(mfrow=c(1,1))
interaction.plot(popi_t2_asv$nutrient,popi_t2_asv$flow,popi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist order level for nutrient and flow for Pielous Evenness")
interaction.plot(popi_t3_asv$nutrient,popi_t3_asv$flow,popi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist order level for nutrient and flow for Pielous Evenness")

par(mfrow=c(1,1))
interaction.plot(popi_t2_asv$nutrient,popi_t2_asv$sediment,popi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist order level for nutrient and sediment for Pielous Evenness")
interaction.plot(popi_t3_asv$nutrient,popi_t3_asv$sediment,popi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist order level for nutrient and sediment for Pielous Evenness")

par(mfrow=c(1,1))
interaction.plot(popi_t2_asv$sediment,popi_t2_asv$flow,popi_t2_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 2 interaction plot for protist order level for sediment and flow for Pielous Evenness")
interaction.plot(popi_t3_asv$sediment,popi_t3_asv$flow,popi_t3_asv$pielou,type = "b",col=c(1:3),
                 leg.bty ="o",leg.bg="beige", lwd=2, pch=c(18,24,22),
                 xlab="Nutrient", ylab="Pielous Evenness",main="Week 3 interaction plot for protist order level for sediment and flow for Pielous Evenness")



# BETA-DIVERSITY INDICES

# According to Ollberding (2019), the PERMANOVA and PERMDISP analysis were performed.

# Permutational multivariate analysis of variance (PERMANOVA)
# Followed by https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova

## For the algae phylum level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_algae_phylum_asv <- phyloseq::distance(sub_algae_phylum_tsc_minthreshold_asv, method = "bray")
## To make a data frame from the sample_data
perm_algae_phylum_sampledf_asv <- data.frame(sample_data(sub_algae_phylum_tsc_minthreshold_asv))
## To make the varisables as factors
str(perm_algae_phylum_sampledf_asv)
perm_algae_phylum_sampledf_1_asv <- perm_algae_phylum_sampledf_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(perm_algae_phylum_sampledf_1_asv)
## To perform the adonis test
adonis2(bray_algae_phylum_asv ~ nutrient*sediment*flow*time, strata = perm_algae_phylum_sampledf_1_asv$week, data = perm_algae_phylum_sampledf_1_asv)

## For the algae order level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_algae_order_asv <- phyloseq::distance(sub_algae_order_tsc_minthreshold_asv, method = "bray")
## To make a data frame from the sample_data
perm_algae_order_sampledf_asv <- data.frame(sample_data(sub_algae_order_tsc_minthreshold_asv))
## To make the variables as factors
str(perm_algae_order_sampledf_asv)
perm_algae_order_sampledf_1_asv <- perm_algae_order_sampledf_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(perm_algae_order_sampledf_1_asv)
## To perform the adonis test
adonis2(bray_algae_order_asv ~ nutrient*sediment*flow*time, strata = perm_algae_order_sampledf_1_asv$week, data = perm_algae_order_sampledf_1_asv)

## For the protist phylum level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_protist_phylum_asv <- phyloseq::distance(sub_protist_phylum_tsc_minthreshold_asv, method = "bray")
## To make a data frame from the sample_data
perm_protist_phylum_sampledf_asv <- data.frame(sample_data(sub_protist_phylum_tsc_minthreshold_asv))
## To make the variables as factors
str(perm_protist_phylum_sampledf_asv)
perm_protist_phylum_sampledf_1_asv <- perm_protist_phylum_sampledf_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(perm_protist_phylum_sampledf_1_asv)
## To perform the adonis test
adonis2(bray_protist_phylum_asv ~ nutrient*sediment*flow*time, strata = perm_protist_phylum_sampledf_1_asv$week, data = perm_protist_phylum_sampledf_1_asv)

## For the protist order level
set.seed(1)
## To calculate the Bray-Curtis distance matrix
bray_protist_order_asv <- phyloseq::distance(sub_protist_order_tsc_minthreshold_asv, method = "bray")
## To make a data frame from the sample_data
perm_protist_order_sampledf_asv <- data.frame(sample_data(sub_protist_order_tsc_minthreshold_asv))
## To make the variables as factors
str(perm_protist_order_sampledf_asv)
perm_protist_order_sampledf_1_asv <- perm_protist_order_sampledf_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(perm_protist_order_sampledf_1_asv)
## To perform the adonis test
adonis2(bray_protist_order_asv ~ nutrient*sediment*flow*time, strata = perm_protist_order_sampledf_1_asv$week, data = perm_protist_order_sampledf_1_asv)


# Permutational analysis of multivariate dispersion (PERMDISP)

## For the algae phylum level
dispersion_algae_phylum_asv <- vegan::betadisper(bray_algae_phylum_asv, phyloseq::sample_data(sub_algae_phylum_tsc_minthreshold_asv)$treatment)
dispersion_algae_phylum_asv
plot(dispersion_algae_phylum_asv, main = "For the algae phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_algae_phylum_asv, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_phylum_asv)

## For the algae order level
dispersion_algae_order_asv <- vegan::betadisper(bray_algae_order_asv, phyloseq::sample_data(sub_algae_order_tsc_minthreshold_asv)$treatment)
dispersion_algae_order_asv
plot(dispersion_algae_order_asv, main = "For the algae order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_algae_order_asv, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_order_asv)

## For the protist phylum level
dispersion_protist_phylum_asv <- vegan::betadisper(bray_protist_phylum_asv, phyloseq::sample_data(sub_protist_phylum_tsc_minthreshold_asv)$treatment)
dispersion_protist_phylum_asv
plot(dispersion_protist_phylum_asv, main = "For the protist phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_protist_phylum_asv, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_phylum_asv)

## For the protist order level
dispersion_protist_order_asv <- vegan::betadisper(bray_protist_order_asv, phyloseq::sample_data(sub_protist_order_tsc_minthreshold_asv)$treatment)
dispersion_protist_order_asv
plot(dispersion_protist_order_asv, main = "For the protist order level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")
boxplot(dispersion_protist_order_asv, main = "", xlab = "")
## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_order_asv)


# The PCoA analysis (Principal Coordinate Analysis)
# Followed by, https://bioconductor.statistik.tu-dortmund.de/packages/3.4/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html, this statistical analysis was performed.

## For the algae phylum level
ap.pcoa_asv = ordinate(sub_algae_phylum_tsc_minthreshold_asv, method="PCoA", distance="bray")
ap_p_asv = plot_ordination(sub_algae_phylum_tsc_minthreshold_asv, ap.pcoa_asv, "samples", color="treatment")+
  geom_line() + geom_point(size=5)+ geom_path()+
  facet_wrap(~week, scales="free")+
  ggtitle("PCoA of the different treatments for week 2 and week 3 for the algal phylum level")
ap_p_asv

## For the algae order level
ao.pcoa_asv = ordinate(sub_algae_order_tsc_minthreshold_asv, method="PCoA", distance="bray")
ao_p_asv = plot_ordination(sub_algae_order_tsc_minthreshold_asv, ao.pcoa_asv, "samples", color="treatment")+
  geom_line() + geom_point(size=5)+ geom_path()+
  facet_wrap(~week, scales="free")+
  ggtitle("PCoA of the different treatments for week 2 and week 3 for the algal order level")
ao_p_asv

## For the protist phylum level
pp.pcoa_asv = ordinate(sub_protist_phylum_tsc_minthreshold_asv, method="PCoA", distance="bray")
pp_p_asv = plot_ordination(sub_protist_phylum_tsc_minthreshold_asv, pp.pcoa_asv, "samples", color="treatment")+
  geom_line() + geom_point(size=5)+ geom_path()+
  facet_wrap(~week, scales="free")+
  ggtitle("PCoA of the different treatments for week 2 and week 3 for the protist phylum level")
pp_p_asv

## For the protist order level
po.pcoa_asv = ordinate(sub_protist_order_tsc_minthreshold_asv, method="PCoA", distance="bray")
po_p_asv = plot_ordination(sub_protist_order_tsc_minthreshold_asv, po.pcoa_asv, "samples", color="treatment")+
  geom_line() + geom_point(size=5)+ geom_path()+
  facet_wrap(~week, scales="free")+
  ggtitle("PCoA of the different treatments for week 2 and week 3 for the protist order level")
po_p_asv


# Followed by https://github.com/joey711/phyloseq/issues/847, the most abundant taxa was calculated. 
## To find the most abundant taxa
## Function
most_abundant_taxa_asv <- function(x,taxa){
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
prune_algae_phylum_asv <- prune_taxa(taxa_sums(sub_algae_phylum_tsc_minthreshold_asv) > 0.01, sub_algae_phylum_tsc_minthreshold_asv)
abundant_algae_phylum_asv <- most_abundant_taxa_asv(prune_algae_phylum_asv,"Phylum")
abundant_algae_phylum_asv

prune_algae_order_asv <- prune_taxa(taxa_sums(sub_algae_order_tsc_minthreshold_asv) > 0.01, sub_algae_order_tsc_minthreshold_asv)
abundant_algae_order_asv <- most_abundant_taxa_asv(prune_algae_order_asv,"Order")
abundant_algae_order_asv

prune_protist_phylum_asv <- prune_taxa(taxa_sums(sub_protist_phylum_tsc_minthreshold_asv) > 0.01, sub_protist_phylum_tsc_minthreshold_asv)
abundant_protist_phylum_asv <- most_abundant_taxa_asv(prune_protist_phylum_asv,"Phylum")
abundant_protist_phylum_asv

prune_protist_order_asv <- prune_taxa(taxa_sums(sub_protist_order_tsc_minthreshold_asv) > 0.01, sub_protist_order_tsc_minthreshold_asv)
abundant_protist_order_asv <- most_abundant_taxa_asv(prune_protist_order_asv,"Order")
abundant_protist_order_asv

# Multivariate Analysis Of Variance (MANOVA)
## According to Kassambara (2017); Ben-Shachar et al., (2020), this statistical analysis was performed.

## For the most abundant taxa of the algal phylum level
## To make a sample data frame
manova_df_algae_phylum_asv = as(sample_data(prune_algae_phylum_asv), "data.frame")
manova_df_algae_phylum_asv
## To combine the sample data frame and abundant taxa data frame
manova_df_algae_phylum_1_asv = cbind(manova_df_algae_phylum_asv, abundant_algae_phylum_asv)
manova_df_algae_phylum_1_asv

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Ochrophyta, Chlorophyta_ph and Phragmoplastophyta 
shapiro.test(manova_df_algae_phylum_1_asv$Ochrophyta)
Ochrophyta_1 <- log10(manova_df_algae_phylum_1_asv$Ochrophyta)
shapiro.test(Ochrophyta_1)

shapiro.test(manova_df_algae_phylum_1_asv$Chlorophyta_ph)

shapiro.test(manova_df_algae_phylum_1_asv$Phragmoplastophyta)
Phragmoplastophyta_1 <- sqrt(manova_df_algae_phylum_1_asv$Phragmoplastophyta)
shapiro.test(Phragmoplastophyta_1)

## To combine the transformed data
manova_df_algae_phylum_11_asv <- cbind(manova_df_algae_phylum_1_asv, Ochrophyta_1, Phragmoplastophyta_1)
## To make the variables as factors
str(manova_df_algae_phylum_11_asv)
manova_df_algae_phylum_111_asv <- manova_df_algae_phylum_11_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(manova_df_algae_phylum_111_asv)
## To do the leveneTest
leveneTest(Ochrophyta_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111_asv)
leveneTest(Chlorophyta_ph ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111_asv)
leveneTest(Phragmoplastophyta_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111_asv)
## To do the MANOVA test
manova_abundant_algae_phylum_asv <- manova(cbind(Ochrophyta_1, Chlorophyta_ph, Phragmoplastophyta_1) ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_111_asv)
manova_abundant_algae_phylum_asv
summary.aov(manova_abundant_algae_phylum_asv)
eta_squared(aov(manova_abundant_algae_phylum_asv))

## For the most abundant taxa of the algal order level
## To make a sample data frame
manova_df_algae_order_asv = as(sample_data(prune_algae_order_asv), "data.frame")
manova_df_algae_order_asv
## To combine the sample data frame and abundant taxa data frame
manova_df_algae_order_1_asv = cbind(manova_df_algae_order_asv, abundant_algae_order_asv)
manova_df_algae_order_1_asv

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Sphaeropleales, Bacillariophytina, Chromulinales, Chlorellales, Chlamydomonadales and Coscinodiscophytina  
shapiro.test(manova_df_algae_order_1_asv$Sphaeropleales)

shapiro.test(manova_df_algae_order_1_asv$Bacillariophytina)
Bacillariophytina_1 <- sqrt(manova_df_algae_order_1_asv$Bacillariophytina)
shapiro.test(Bacillariophytina_1)

shapiro.test(manova_df_algae_order_1_asv$Chromulinales)

shapiro.test(manova_df_algae_order_1_asv$Chlorellales)
Chlorellales_1 <- sqrt(manova_df_algae_order_1_asv$Chlorellales)
shapiro.test(Chlorellales_1)

shapiro.test(manova_df_algae_order_1_asv$Chlamydomonadales)
Chlamydomonadales_1 <- sqrt(manova_df_algae_order_1_asv$Chlamydomonadales)
shapiro.test(Chlamydomonadales_1)

shapiro.test(manova_df_algae_order_1_asv$Coscinodiscophytina)
Coscinodiscophytina_1 <- sqrt(manova_df_algae_order_1_asv$Coscinodiscophytina)
shapiro.test(Coscinodiscophytina_1)

## To combine the transformed data
manova_df_algae_order_11_asv <- cbind(manova_df_algae_order_1_asv, Bacillariophytina_1, Chlorellales_1, Chlamydomonadales_1, Coscinodiscophytina_1)
## To make the variables as factors
str(manova_df_algae_order_11_asv)
manova_df_algae_order_111_asv <- manova_df_algae_order_11_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(manova_df_algae_order_111_asv)
## To do the leveneTest
leveneTest(Sphaeropleales ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
leveneTest(Bacillariophytina_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
leveneTest(Chromulinales ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
leveneTest(Chlorellales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
leveneTest(Chlamydomonadales_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
leveneTest(Coscinodiscophytina_1 ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
## To do the MANOVA test
manova_abundant_algae_order_asv <- manova(cbind(Sphaeropleales, Bacillariophytina_1, Chromulinales, Chlorellales_1, Chlamydomonadales_1, Coscinodiscophytina_1) ~ nutrient*sediment*flow*time, data = manova_df_algae_order_111_asv)
manova_abundant_algae_order_asv
summary.aov(manova_abundant_algae_order_asv)
eta_squared(aov(manova_abundant_algae_order_asv))

## For the most abundant taxa of the protist phylum level
## To make a sample data frame
manova_df_protist_phylum_asv = as(sample_data(prune_protist_phylum_asv), "data.frame")
manova_df_protist_phylum_asv
## To combine the sample data frame and abundant taxa data frame
manova_df_protist_phylum_1_asv = cbind(manova_df_protist_phylum_asv, abundant_protist_phylum_asv)
manova_df_protist_phylum_1_asv

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Cercozoa, LKM15, Ciliophora and Tubulinea 
shapiro.test(manova_df_protist_phylum_1_asv$Cercozoa)

shapiro.test(manova_df_protist_phylum_1_asv$LKM15)

shapiro.test(manova_df_protist_phylum_1_asv$Ciliophora)
Ciliophora_1 <- sqrt(manova_df_protist_phylum_1_asv$Ciliophora)
shapiro.test(Ciliophora_1)

shapiro.test(manova_df_protist_phylum_1_asv$Tubulinea)
Tubulinea_1 <- sqrt(manova_df_protist_phylum_1_asv$Tubulinea)
shapiro.test(Tubulinea_1)
## To combine the transformed data
manova_df_protist_phylum_11_asv <- cbind(manova_df_protist_phylum_1_asv, Ciliophora_1, Tubulinea_1)
## To make the variables as factors
str(manova_df_protist_phylum_11_asv)
manova_df_protist_phylum_111_asv <- manova_df_protist_phylum_11_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(manova_df_protist_phylum_111_asv)
## To do the leveneTest
leveneTest(Cercozoa ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111_asv)
leveneTest(LKM15 ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111_asv)
leveneTest(Ciliophora_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111_asv)
leveneTest(Tubulinea_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111_asv)
## To do the MANOVA test
manova_abundant_protist_phylum_asv <- manova(cbind(Cercozoa, LKM15, Ciliophora_1, Tubulinea_1) ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_111_asv)
manova_abundant_protist_phylum_asv
summary.aov(manova_abundant_protist_phylum_asv)
eta_squared(aov(manova_abundant_protist_phylum_asv))

## For the most abundant taxa of the protist order level
## To make a sample data frame
manova_df_protist_order_asv = as(sample_data(prune_protist_order_asv), "data.frame")
manova_df_protist_order_asv
## To combine the sample data frame and abundant taxa data frame
manova_df_protist_order_1_asv = cbind(manova_df_protist_order_asv, abundant_protist_order_asv)
manova_df_protist_order_1_asv

## Data normality check and transformation of data to the normality for the most abundant taxa
### For Litostomatea, Spirotrichea, Cryomonadida, Echinamoebida, Hypocreales, Granofilosea_or, Conthreep and Cercomonadidae_or
shapiro.test(manova_df_protist_order_1_asv$Litostomatea)
Litostomatea_1 <- sqrt(manova_df_protist_order_1_asv$Litostomatea)
shapiro.test(Litostomatea_1)

shapiro.test(manova_df_protist_order_1_asv$Spirotrichea)
Spirotrichea_1 <- sqrt(manova_df_protist_order_1_asv$Spirotrichea)
shapiro.test(Spirotrichea_1)

shapiro.test(manova_df_protist_order_1_asv$Cryomonadida)
Cryomonadida_1 <- sqrt(manova_df_protist_order_1_asv$Cryomonadida)
shapiro.test(Cryomonadida_1)

shapiro.test(manova_df_protist_order_1_asv$Echinamoebida)
Echinamoebida_1 <- sqrt(manova_df_protist_order_1_asv$Echinamoebida)
shapiro.test(Echinamoebida_1)

shapiro.test(manova_df_protist_order_1_asv$Hypocreales)
Hypocreales_1 <- sqrt(manova_df_protist_order_1_asv$Hypocreales)
shapiro.test(Hypocreales_1)

shapiro.test(manova_df_protist_order_1_asv$Granofilosea_or)
Granofilosea_or_1 <- sqrt(manova_df_protist_order_1_asv$Granofilosea_or)
shapiro.test(Granofilosea_or_1)

shapiro.test(manova_df_protist_order_1_asv$Conthreep)
Conthreep_1 <- sqrt(manova_df_protist_order_1_asv$Conthreep)
shapiro.test(Conthreep_1)

shapiro.test(manova_df_protist_order_1_asv$Cercomonadidae_or)
Cercomonadidae_or_1 <- sqrt(manova_df_protist_order_1_asv$Cercomonadidae_or)
shapiro.test(Cercomonadidae_or_1)

## To combine the transformed data
manova_df_protist_order_11_asv <- cbind(manova_df_protist_order_1_asv, Litostomatea_1, Spirotrichea_1, Cryomonadida_1, Echinamoebida_1, Hypocreales_1, Granofilosea_or_1, Conthreep_1, Cercomonadidae_or_1)
## To make the variables as factors
str(manova_df_protist_order_11_asv)
manova_df_protist_order_111_asv <- manova_df_protist_order_11_asv %>% mutate(across(c(week, flow, sediment, nutrient, time, header_tanks_block, treatment_week), as.factor))
str(manova_df_protist_order_111_asv)
## To do the leveneTest
leveneTest(Litostomatea_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Spirotrichea_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Cryomonadida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Echinamoebida_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Hypocreales_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Granofilosea_or_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Conthreep_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
leveneTest(Cercomonadidae_or_1 ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
## To do the MANOVA test
manova_abundant_protist_order_asv <- manova(cbind(Litostomatea_1, Spirotrichea_1, Cryomonadida_1, Echinamoebida_1, Hypocreales_1, Granofilosea_or_1, Conthreep_1, Cercomonadidae_or_1) ~ nutrient*sediment*flow*time, data = manova_df_protist_order_111_asv)
manova_abundant_protist_order_asv
summary.aov(manova_abundant_protist_order_asv)
eta_squared(aov(manova_abundant_protist_order_asv))




# List of References

## Vaulot, D. (2021, Feb 15). Phyloseq tutorial. https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

## McMurdie PJ, Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." _PLoS ONE_, *8*(4), e61217. <http://dx.plos.org/10.1371/journal.pone.0061217>.

## Mariadassou, M., Bernard, M., Pascal, G., Cauquil, L., & Chaillou, S. (2016). Analysis of community composition data using phyloseq. Montpellier Décembre 2016. https://genoweb.toulouse.inra.fr/~formation/15_FROGS/8-February2017/FROGS_phyloseq_02_2017.pdf

## Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M. (2009). Mixed effects models and extensions in ecology with R (Vol. 574). New York: Springer.

## Ollberding, N.J. (2019, Jul 28). Introduction to the Statistical Analysis of Microbiome Data in R. Wowchemy. https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

## Kassambara, A. (2017). MANOVA Test in R: Multivariate Analysis of Variance. Statistical tools for high-throughput data analysis. STHDA. http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance#infos

## Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815  

