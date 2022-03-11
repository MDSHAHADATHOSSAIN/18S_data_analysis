## Eucaryotic data analysis with phyloseq (week 3)

## To install necessary packages
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics

source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")  

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

## To load the installed packages with libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble") 

## To read the data and create phyloseq objects
otu_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week3/eucaryotic_week3.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week3/eucaryotic_week3.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("D:/HIWI/eukaryote_otu/shahadat_euK/eucaryotes_analysis/Week3/eucaryotic_week3.xlsx", sheet = "Samples")

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

carbom <- phyloseq(OTU, TAX, samples)
carbom

## To visualize data

sample_names(carbom)

rank_names(carbom)

sample_variables(carbom)

## To keep only samples to be analyzed
carbom <- subset_samples(carbom, treatment %in% c("C", "N", "S", "F", "NS", "NF", "SF", "NFS"))
carbom

## To keep only taxa of interests (mostly eucaryotes)
carbom <- subset_taxa(carbom, phylum %in% c("Bacillariophyta", " Apicomplexa", "Ascomycota", "Cryptophyta", 
                                            "Chlorophyta", "LKM15", "Cercozoa", "Ochrophyta", "Xanthophyceae", "Streptophyta"))
carbom <- subset_taxa(carbom, !(genus %in% c("Melosira", "Verticillium", "Messastrum", "unclassified_Glissomonadida", "Pleurochloris")))
carbom

## To normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
## To keep number of reads that are >= 0.003
carbom = filter_taxa(carbom, function(x) mean(x) >= 0.003, TRUE)

## Bar graphs

## Basic bar graph based on phylum
plot_bar(carbom, fill = "phylum")

## To make the bargraph nicer by removing OTUs boundaries. This is done by adding ggplot2 modifier
plot_bar(carbom, fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

## To regroup together according to treatment_week
carbom_fraction <- merge_samples(carbom, "treatment_week")
plot_bar(carbom_fraction, fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")

## To keep only Chlorophyta and use color according to genus. Do separate panels week and treatment_week samples
carbom_chloro <- subset_taxa(carbom, phylum %in% c("Chlorophyta"))
plot_bar(carbom_chloro, x="genus", fill = "genus", facet_grid = week~treatment_week) +
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

## Heatmaps

## A basic heatmap using the default parameters
plot_heatmap(carbom, method = "NMDS", distance = "bray")

## The heatmap is very cluttered. It is better to only consider the most abundant OTUs for heatmaps
## One can only take OTUs that represent at least 20% of reads in at least one sample. 
## But here at least 3% of reads considered because only few taxa left if 20% considered
## To remember that we normalized all the sampples to median number of reads (total)
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.030) > 0, TRUE)
carbom_abund

otu_table(carbom_abund)[1:48, 1:32]

plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

## It is possible to use different distances and different multivaraite methods.
## For example Jaccard distance and MDS and label OTUs with genus, order by genus.
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "genus", taxa.order = "genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

## Many different built-in distances can be used
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

## Another strategy is to do a heatmap for a specific taxonomy group
## For example we can taget the Chlorophyta and then label the OTUs using the genus
plot_heatmap(carbom_chloro, method = "NMDS", distance = "bray", 
             taxa.label = "genus", taxa.order = "genus", 
             low="beige", high="red", na.value="beige")

## Alpha diversity
## Plot Chao1 richness estimator and Shannon diversity estimator
plot_richness(carbom, measures=c("Chao1", "Shannon"))

## Regroup together samples from the same treatment_week
plot_richness(carbom, measures=c("Chao1", "Shannon"), x="week", color="treatment_week")

## Ordination
## To do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
carbom.ord <- ordinate(carbom, "NMDS", "bray")

## Plot OTUs
plot_ordination(carbom, carbom.ord, type="taxa", color="class", shape= "phylum", 
                title="OTUs")

## It is a bit confusing, so make it more easy to visualize by breaking according to taxonomic division
plot_ordination(carbom, carbom.ord, type="taxa", color="class", 
                title="OTUs", label="class") + 
  facet_wrap(~phylum, 3)

## To display samples and enlarge the points to make it more easy to read
plot_ordination(carbom, carbom.ord, type="samples", color="treatment_week", 
                shape="treatment", title="Samples") + geom_point(size=3)

## To display both samples and OTUs but in 2 different panels
plot_ordination(carbom, carbom.ord, type="split", color="phylum", 
                shape="treatment", title="biplot", label = "station") +  
  geom_point(size=3)

## Network analysis
## Simple network analysis
plot_net(carbom, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="phylum", point_label="genus")

## This is quite confusing. Let us make it more simple by using only major OTUs
plot_net(carbom_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="phylum", point_label="genus")

