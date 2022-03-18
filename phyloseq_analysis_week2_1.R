## Eucaryotic data analysis with phyloseq (week 2)

## To install necessary packages
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics

source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")  

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

## To read the data and create phyloseq objects
otu_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week2/eucaryotic_week2.xlsx", sheet = "Samples")

## Phyloseq objects need to have row.names

## To define the row names from the otu column
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 

## To idem for the two other matrixes 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

## To transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

## To transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

phyloseq_obj <- phyloseq(OTU, TAX, samples)
phyloseq_obj

## To normalize number of reads in each sample using median sequencing depth and filter the values that are above 1e-5
phy_norm  = transform_sample_counts(phyloseq_obj, function(x) x / sum(x) ) ## normalization
phy_filtered = filter_taxa(phy_norm, function(x) mean(x) > 1e-5, TRUE) ## filtering

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(phy_filtered), "matrix")
# transpose if necessary
if(taxa_are_rows(phy_filtered)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf

# Extract taxa matrix from the phyloseq object
TAX1 = as(tax_table(phy_filtered), "matrix")
# transpose if necessary
if(taxa_are_rows(phy_filtered)){TAX1 <- t(TAX1)}
# Coerce to data.frame
TAXdf = as.data.frame(TAX1)
TAXdf

# Extract sample matrix from the phyloseq object
samples1 = as(sample_data(phy_filtered), "matrix")
# transpose if necessary
if(taxa_are_rows(phy_filtered)){samples1 <- t(samples1)}
# Coerce to data.frame
samplesdf = as.data.frame(samples1)
samplesdf

## To visualize data

sample_names(phy_filtered)

rank_names(phy_filtered)

sample_variables(phy_filtered)



## To keep only taxa of interests (Algae) according to phylum and genus level
sub_algae_phylum <- subset_taxa(phy_filtered, phylum %in% c("Chlorophyta", "Bacillariophyta", "Ochrophyta", "Eustigmatophyceae", "Dinoflagellata", "Phragmoplastophyta", "Rigifilida", "Schizoplasmodiida", "Streptophyta", "Xanthophyceae"))
sub_algae_genus <- subset_taxa(phy_filtered, genus %in% c("Messastrum", "Melosira", "Poteriospumella", "Pseudellipsoidion", "Symbiodinium", "Nasturtium", "unclassified_Rigifilida", "Phalansterium", "Reynoutria", "Botryochloridaceae"))

## To keep only taxa of interests (Protists) according to phylum and genus level
sub_protist_phylum <- subset_taxa(phy_filtered, phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Rotifera", "Tubulinea", "Bicosoecida", "Discosea", "Gastrotricha", "Labyrinthulomycetes", "Peronosporomycetes"))
sub_protist_genus <- subset_taxa(phy_filtered, genus %in% c("Gymnophrys", "unclassified_Craspedida", "Gonostomum", "Cephalodella", "Flabellula", "Bicosoeca", "unclassified_Discosea", "Lepidodermella", "Sorodiplophrys", "Phytophthora"))

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




## Exploring biodiversity : £\-diversity
## Richness are plotted with plot_richness
# We are interested in alpha diversity
# A good idea is to prune OTUs that are not present in any of the samples
# But we can not trim more than that
# Many richness estimates are modeled on singletons and doubletons in the abundance data
# We need to leave them in the dataset if we want a meaningful estimate

phy_alpha_prune <- prune_taxa(taxa_sums(phyloseq_obj) > 0 , phyloseq_obj)

# To try it on trimmed phyloseq object
p1 <- plot_richness(phy_alpha_prune)
p1
# To custom it on treatment, color treatment, select measures as "Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", and plot as boxplot
p2 <- plot_richness(phy_alpha_prune, color = "treatment", x = "treatment", measures = c("Observed", "Chao1", "Shannon", "Simpson", 
                                                                       "InvSimpson"))
p2
# To plot as boxplot
p3 <- p2 + geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(p3)



# The numeric values of £\-diversities are given by estimate_richness (used internally by plot_richness)
alpha.diversity <- estimate_richness(phy_alpha_prune, measures = c("Observed", 
                                                            "Chao1", "Shannon", "Simpson", 
                                                            "InvSimpson"))
alpha.diversity
head(alpha.diversity)

write.table(alpha.diversity, "myfile.txt")

# A quick ANOVA : tests if observed richness is significantly different in function of treatment
data <- cbind(sample_data(phy_alpha_prune), alpha.diversity)
phy_alpha_prune_obs.anova <- aov(Observed ~ treatment, data)
summary(phy_alpha_prune_obs.anova)
# There is a significant effect of treatment type on richness

# A quick ANOVA : tests if Shannon indices is significantly different in function of treatment
phy_alpha_prune_sha.anova <- aov(Shannon ~ treatment, data)
summary(phy_alpha_prune_sha.anova)

# There is a significant effect of treatment type on richness

# A quick ANOVA : tests if Chao1 indices is significantly different in function of treatment
phy_alpha_prune_chao1.anova <- aov(Chao1 ~ treatment, data)
summary(phy_alpha_prune_chao1.anova)

# There is a significant effect of treatment type on richness

# A quick ANOVA : tests if Simpson indices is significantly different in function of treatment
phy_alpha_prune_sim.anova <- aov(Simpson ~ treatment, data)
summary(phy_alpha_prune_sim.anova)

# There is a significant effect of treatment type on richness

# A quick ANOVA : tests if InvSimpson indices is significantly different in function of treatment
phy_alpha_prune_invsim.anova <- aov(InvSimpson ~ treatment, data)
summary(phy_alpha_prune_invsim.anova)
# treatment effect on InvSimpson diversity is not significant




####################


### BETA-DIVERSITY INDICES

# To create a phy_tree slot from the original phyloseq object
random_tree = rtree(ntaxa(phyloseq_obj), rooted=TRUE, tip.label=taxa_names(phyloseq_obj))
plot(random_tree)

# To make a new phyloseq object with the newly created phy_tree
phyloseq_obj_1 <- phyloseq(OTU, TAX, samples, random_tree)
phyloseq_obj_1

# To prune OTUs that are not present in any of the samples
phy_beta_prune <- prune_taxa(taxa_sums(phyloseq_obj_1) > 0 , phyloseq_obj_1)

## Generic method for distance
dist.bc <- distance(phy_beta_prune, method = "bray") ## Bray-Curtis
dist.bc
## all available distances
distanceMethodList
## or distance("list") in previous version of phyloseq

## compute other matrix distances
dist.jac <- distance(phy_beta_prune, "cc")
dist.jac

# phy_tree slot is needed for unifrac and weighted unifrac matrix distances
dist.uf <- distance(phy_beta_prune, method = "unifrac") ## Unifrac
dist.uf
dist.wuf <- distance(phy_beta_prune, method = "wunifrac") ## Weighted Unifrac
dist.wuf

## some code to visualize distances
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

## Compare Bray-Curtis (bray) and Jaccard (cc)
p8 <- plot_dist_as_heatmap(dist.bc, title = "Bray-Curtis")
plot(p8)
p9 <- plot_dist_as_heatmap(dist.jac, title = "Jaccard")
plot(p9)
grid.arrange(p8, p9, ncol = 2, nrow = 1)

## Compare unifrac and wunifrac
p10 <- plot_dist_as_heatmap(dist.uf, title = "Unifrac")
plot(p10)
p11 <- plot_dist_as_heatmap(dist.wuf, title = "wUnifrac")
plot(p11)
grid.arrange(p10, p11, ncol = 2, nrow = 1)


## compare Bray-Curtis, Jaccard, Unifrac and wUnifrac
grid.arrange(p8, p9, p10, p11, ncol = 2, nrow = 2)

