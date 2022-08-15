# A workflow of DADA2 for Big Data (based on versions 1.4, 1.16 and 1.8)
# The workflow of DADA2 was performed according to Callahan et al., (2016); Callahan (2018). 

# Analysis with DADA2 for the 2_Merged_data
## The starting point of this workflow is paired-end Fastq files (Illumina sequencing) demultiplexed from samples.
## The barcodes or adapters had already been removed from the files.
## The final product of this analysis is an amplicon sequence variant (ASV) table, which can be expressed as the higher-resolution analog of a conventional OTU table in which the number of occurrences or times of each amplicon sequence variant for each of the samples is recorded and stored.
## The taxonomy of the output sequences is also assigned and described here. 
## The data obtained can be imported into the well-known R package phyloseq for further analysis.

# That's the starting point!
## This DADA2 pipeline expects the sequencing data to conform to certain characteristics: (1) The samples must be demultiplexed to split into separate fastq files per sample, (2) The non-biological nucleotides must be removed, such as primers, adapters, linkers, etc., (3) For paired-end sequencing data, the reads in the forward and reverse fastq files should be in a matching order.  

# Get ready to see the beauty of the DADA2 workflow now!

## To install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

install.packages("pacman")
install.packages("here")
install.packages("xlsx")
install.packages("phyloseq")
install.packages("ggplot2")
install.packages("readxl")
install.packages("dplyr")
install.packages("tibble")
install.packages("microbiome")
install.packages("reshape2")
install.packages("ape")
install.packages("gridExtra")
install.packages("plotly")
install.packages("vegan")
install.packages("dendextend")
install.packages("tidyr")
install.packages("rms")
install.packages("effectsize")
install.packages("curatedMetagenomicData")
install.packages("lme4")
install.packages("picante")
install.packages("cowplot")
install.packages("ShortRead")
install.packages("lsr")
install.packages("MASS")
insertClassMethods("rcompanion")
install.packages("ggiraphExtra")
install.packages("optimx")
install.packages("mvabund")
install.packages("gllvm")
install.packages("radiant.data")
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

## To load the installed packages with the pacman package
pacman::p_load(dada2, pacman, here, xlsx, phyloseq, ggplot2, readxl, dplyr, tibble, microbiome, reshape2, ape, gridExtra, plotly, vegan, dendextend, tidyr, rms, effectsize, curatedMetagenomicData, lme4, picante, cowplot, ShortRead, lsr, MASS, rcompanion, ggiraphExtra, optimx, mvabund, gllvm, radiant.data, tidyverse, MuMIn, corrplot, gclus, coefplot, jtools, ggstance, interactions, blmeco, tseries)
here()

## To load the directory where the demultiplexed fastq files are located
mydata_path <- here("2_Merged_data_d")

## The inclusion of the filtered files in the filtered/ subdirectory folder
mydata_filtpath <- file.path(mydata_path, "filtered_d")

## To change if there are different file extensions
## The format of the file names is expected as follows: Sample_Group_Sample_Barcode_Raw_num.fastq.
mydata_fns <- list.files(mydata_path, pattern=".fastq") 
mydata_fns

## To inspect the quality profiles read
plotQualityProfile(here("2_Merged_data_d","A1_18S191096_A1_CCAACA_55154.fastq"))
### In each base position, the gray scale is a heat map for the frequency of each quality score.
### The green line indicates the median quality score in each position.
### The orange line shows the quartiles of the distribution of the quality score.
### The red line represents the scale proportion of reads that reach at least to this position.
### At position 240, the reads are truncated (the last few nucleotides are trimmed) to avoid the less well-controlled errors that can occur there.

## Now comes the most interesting step, filtering and trimming!
filterAndTrim(file.path(mydata_path,mydata_fns), file.path(mydata_filtpath,mydata_fns), 
              truncLen=240, maxN=0, maxEE=10, truncQ=1, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=FALSE)
### Here, the standard filter parameters were used. 
### The maximum number of expected errors allowed in a read was set by the maxEE parameter, which can be considered as a good filter compared to a simple averaging of the quality score.

mydata_out <- filterAndTrim(file.path(mydata_path,mydata_fns), file.path(mydata_filtpath,mydata_fns), 
                            truncLen=240, maxN=0, maxEE=10, truncQ=1, rm.phix=TRUE,
                            compress=TRUE, verbose=TRUE, multithread=FALSE)
mydata_out

## To load the directory with the filtered fastq files
mydata_filtpath_filtered <- here("2_Merged_data_d/filtered_d")

## To change if there are different file extensions
mydata_filts_filtered <- list.files(mydata_filtpath_filtered, pattern=".fastq", full.names=TRUE)

## Assumed file name = Sample_Group_Sample_Barcode_Raw_num.fastq; for extracting the sample names
mydata_sample.names_filtered <- sapply(strsplit(basename(mydata_filts_filtered), "_"), `[`, 1)

names(mydata_filts_filtered) <- mydata_sample.names_filtered

## How to learn the error rates!
### Each amplicon dataset has different types of error rates and the DADA2 algorithm uses a parametric error model (err).
### The learnErrors process learns this error model from the data sets by alternately estimating the error rates and inferring the sample composition until they converge to a jointly and common consistent solution.
### As with the various machine learning problems, the algorithm begins with an initial guess using the highest possible error rates in the data, e.g., the error rates counted when the most frequently occurring abundant sequence is correct and all others are errors.
set.seed(100)
mydata_err_filtered <- learnErrors(mydata_filts_filtered, nbases = 1e8, multithread=TRUE, randomize=TRUE)

## The visualization of estimated error rates
plotErrors(mydata_err_filtered, nominalQ=TRUE)
### The error rates of all possible transitions (A???C, A???G, .) are very well shown here. 
### The points represent the observed error rates for each consensus quality score. 
### The black line expresses the estimated error rates after the actual convergence of the machine learning algorithms. 
### The red lines here describe the error rates expected based on the nominal definition of the Q-score. 
### The black line of estimated error rates better fits the observed rates (points), with error rates decreasing as expected as the quality score increases.
### We can proceed to the next steps, because everything seems reasonable.

## The dereplication process and inference of sequence variants
### Dereplication process combines all similar sequencing reads into unique sequences and associates them with an abundance equal to the number of reads with that unique sequence.
### By eliminating the redundant comparisons, the dereplication process reduces the computation time.
### The consensus quality profile of a unique sequence expresses the average values of the position qualities obtained from the dereplicated reads.
### The accuracy of DADA2 increases when the quality profile informs the error model for the subsequent sampling inference step.
### The core sample inference algorithm is applied here to the dereplicated data.

mydata_dds <- vector("list", length(mydata_sample.names_filtered))
names(mydata_dds) <- mydata_sample.names_filtered
for(sam in mydata_sample.names_filtered) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(mydata_filts_filtered[[sam]]) # Dereplication 
  mydata_dds[[sam]] <- dada(derep, err=mydata_err_filtered, multithread=TRUE) # Inference
}

## To inspect the object of the returned Dada class
mydata_dds[[1]]
### From the 20030 input unique sequences in the A1 sample, the DADA2 algorithm inferred 209 true sequence variants.
### We can inspect the remaining samples in a similar way.

## To construct the sequence table
## Now it is possible to start creating an amplicon sequence variant (ASV) table, which is a higher resolution version of the OTU table created using the traditional methods.
mydata_seqtab <- makeSequenceTable(mydata_dds)
dim(mydata_seqtab)

## For inspecting the overall distribution of sequence lengths
table(nchar(getSequences(mydata_seqtab)))
### The sequence table is a matrix in which the rows are named by samples and the columns are named by sequence variants.
### There are 7198 ASVs in the sequence table.

## To load the directory where the sequence table will be stored
saveRDS(mydata_seqtab, here("2_Merged_data_d/filtered_d/mydata_seqtab_asv_merged.rds")) 

## The final output of the counting matrix with 64 samples in the rows and the ASV in the columns is now loaded into the R object (mydata_foo_asv) and serialized.
mydata_foo_asv <- readRDS(here("2_Merged_data_d/filtered_d/mydata_seqtab_asv_merged.rds"))
mydata_foo_asv

## It sounds interesting to remove the chimeras!
### Compared to dealing with OTUs, the accuracy of amplicon sequence variants after the denoising process helps to identify chimeras easily.
### If the chimeric sequences are correctly reconstructed from the combination of left and right segments of two more abundant parent sequences, they can be identified. 
mydata_seqtab_removed_chimeras <- removeBimeraDenovo(mydata_foo_asv, method="consensus", multithread=TRUE)
mydata_seqtab_removed_chimeras

dim(mydata_seqtab_removed_chimeras)

sum(mydata_seqtab_removed_chimeras)/sum(mydata_foo_asv)

## To track reads through the pipeline
### We need to check the number of reads that were successful in each step of the pipeline to see the final progress.
mydata_dadaFs <- dada(derep, mydata_err_filtered, multithread=TRUE)

mydata_getN <- function(x) sum(getUniques(x))
mydata_track <- cbind(mydata_out, mydata_getN(mydata_dadaFs), rowSums(mydata_seqtab_removed_chimeras))
colnames(mydata_track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(mydata_track) <- mydata_sample.names_filtered
head(mydata_track)
mydata_track
## Looks good!

## Now comes the most important part! The assignment of taxonomies!
### Assignment of taxonomies to sequence variants specifically based on 18S amplicon sequencing is very common at this stage.
### The DADA2 package uses the native Bayesian classification method to assign taxonomies.
### First, the assignTaxonomy function loads a set of sequences as inputs to classify them. 
### Second, it also loads a trained set of reference sequences whose taxonomy is actually known. 
### Finally, the function outputs the assignment of taxonomies with the lowest minBoot bootstrap confidence.
### According to Callahan (2018), the Fasta file silva_nr_v132_train_set.fa.gz is used for taxonomy assignment and the file silva_species_assignment_v132.fa.gz is used for species assignment.
set.seed(100)
mydata_taxa <- assignTaxonomy(mydata_seqtab_removed_chimeras, tryRC=TRUE, here("silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
mydata_taxa

## Assigning species
set.seed(100)
mydata_taxa_with_species <- addSpecies(mydata_taxa, tryRC=TRUE, here("silva_species_assignment_v132.fa.gz"))
mydata_taxa_with_species

## Let's look at the taxonomy and species assignment together
mydata_taxa_with_species.print <- mydata_taxa_with_species # To remove sequence rownames for display only
rownames(mydata_taxa_with_species.print) <- NULL
head(mydata_taxa_with_species.print)
mydata_taxa_with_species.print

## To write to the disk
### To load the directory where the sequence table will be stored
saveRDS(mydata_seqtab_removed_chimeras, here("2_Merged_data_d/filtered_d/mydata_seqtab_final_asv_merged_d.rds"))
### To load the directory where the taxa table will be stored
saveRDS(mydata_taxa, here("2_Merged_data_d/filtered_d/mydata_taxa_final_asv_merged_d.rds"))
### To load the directory where the taxa with species table will be stored
saveRDS(mydata_taxa_with_species, here("2_Merged_data_d/filtered_d/mydata_taxa_with_species_final_asv_merged_d.rds")) 

### To load the directory where the sequence table will be stored in excel format
write.xlsx(mydata_seqtab_removed_chimeras, file = here("2_Merged_data_d/filtered_d/mydata_seqtab_final_asv_merged_d.xlsx"),
           append = FALSE)
### To load the directory where the taxa table will be stored in excel format
write.xlsx(mydata_taxa, file = here("2_Merged_data_d/filtered_d/mydata_taxa_final_asv_merged_d.xlsx"),
           append = FALSE)
### To load the directory where the taxa with species table will be stored in excel format
write.xlsx(mydata_taxa_with_species, here("2_Merged_data_d/filtered_d/mydata_taxa_with_species_final_asv_merged_d.xlsx"),
           append = FALSE)



# ASV data analysis with the phyloseq package for statistical analysis

## According to Vaulot (2021) and McMurdie & Holmes (2013), this statistical analysis was performed.

## To read the data and create phyloseq objects
asv_mat_1_dada <- read_excel(here("2_Merged_data_d/filtered_d", "dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "ASV matrix")

tax_mat_1_dada <- read_excel(here("2_Merged_data_d/filtered_d", "dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "Taxonomy table")

samples_df_1_dada <- read_excel(here("2_Merged_data_d/filtered_d", "dada2_asv_eucaryotic_merged_d.xlsx"), sheet = "Samples")

## Phyloseq objects must have row.names
## To define the row names from the asv column
asv_mat_1_dada <- asv_mat_1_dada %>%
  tibble::column_to_rownames("asv")

## To idem for the two other matrices
## For taxa data
tax_mat_1_dada <- tax_mat_1_dada %>% 
  tibble::column_to_rownames("asv")

## For sample data
samples_df_1_dada <- samples_df_1_dada %>% 
  tibble::column_to_rownames("sample")

## Transformation to matrices of asv and tax tables (the sample table can be left as a data frame)
asv_mat_1_dada <- as.matrix(asv_mat_1_dada)

tax_mat_1_dada <- as.matrix(tax_mat_1_dada)

## Transforming to phyloseq objects
ASV_1_dada = otu_table(asv_mat_1_dada, taxa_are_rows = TRUE)
TAX_1_dada = tax_table(tax_mat_1_dada)
samples_1_dada = sample_data(samples_df_1_dada)

phyloseq_obj_1_dada <- phyloseq(ASV_1_dada, TAX_1_dada, samples_1_dada)
phyloseq_obj_1_dada

## To create a phy_tree slot from the original phyloseq object
random_tree_dada = rtree(ntaxa(phyloseq_obj_1_dada), rooted=TRUE, tip.label=taxa_names(phyloseq_obj_1_dada))
plot(random_tree_dada)

## To make a new phyloseq object with the newly created phy_tree
phyloseq_object_with_tree_dada <- merge_phyloseq(phyloseq_obj_1_dada, random_tree_dada)
phyloseq_object_with_tree_dada

## The preprocessing of the data
## log10p transformation according to Martins & Moreau (2020).
phyloseq_object_final_dada <- microbiome::transform(phyloseq_object_with_tree_dada, "log10p")
phyloseq_object_final_dada
phyloseq::otu_table(phyloseq_object_final_dada)[1:5, 1:5]
phyloseq::sample_data(phyloseq_object_final_dada)[1:5, 1:5]

## To normalize the number of reads or abundances in each sample based on the median sequencing depth
mydata_total_dada = median(sample_sums(phyloseq_object_final_dada))
mydata_standf_dada = function(x, t=mydata_total_dada) round(t * (x / sum(x)))

## Then, the data of phyloseq_object_final_dada is first transformed into relative abundance data and stored under phy_transformed_dada.
## After that, the phy_transformed_dada object is filtered to phy_filtered_dada to keep the ASVs with a mean >= 0.003.
phy_transformed_dada  = transform_sample_counts(phyloseq_object_final_dada, mydata_standf_dada)
phy_transformed_dada
phy_filtered_dada = filter_taxa(phy_transformed_dada, function(x) mean(x) >= 0.003, TRUE) 
phy_filtered_dada

## To visualize the data
sample_names(phy_filtered_dada)
rank_names(phy_filtered_dada)
sample_variables(phy_filtered_dada)

## Extract abundance matrix from the phyloseq object
ASV1_dada = as(otu_table(phy_filtered_dada), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered_dada)){ASV1_dada <- t(ASV1_dada)}
## Coerce to data.frame
ASVdf_dada = as.data.frame(ASV1_dada)
ASVdf_dada

## Extract taxa matrix from the phyloseq object
TAX1_dada = as(tax_table(phy_filtered_dada), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered_dada)){TAX1_dada <- t(TAX1_dada)}
## Coerce to data.frame
TAXdf_dada = as.data.frame(TAX1_dada)
TAXdf_dada

## Extract sample matrix from the phyloseq object
samples1_dada = as(sample_data(phy_filtered_dada), "matrix")
## Transpose if necessary
if(taxa_are_rows(phy_filtered_dada)){samples1_dada <- t(samples1_dada)}
## Coerce to data.frame
samplesdf_dada = as.data.frame(samples1_dada)
samplesdf_dada

## To keep only taxa of interests (for Algae) according to phylum and genus level
sub_algae_phylum_dada <- subset_taxa(phy_filtered_dada, Phylum %in% c("Ochrophyta", "Dinoflagellata", "Phragmoplastophyta"))
sub_algae_genus_dada <- subset_taxa(phy_filtered_dada, Genus %in% c("Messastrum", "Poteriospumella", "Pseudellipsoidion", "Symbiodinium", "Achnanthidium", "Ancyromonas", "Ankistrodesmus", "Aphanochaete", "Aphelidium", "Carteria", "Characiopodium", "Characium", "Chlamydomonas", "Chlorella", "Chlorochytrium", "Chlorococcum", "Chlorosarcinopsis", "Chroomonas", "Chrysamoeba", "Closterium", "Coccomyxa", "Cocconeis", "Cymbella", "Deasonia", "Desmodesmus", "Dictyochloropsis", "Elliptochloris", "Eustigmatos", "Fragilaria", "Gloeotilopsis", "Gomphonema", "Hemiselmis", "Heterochlorella", "Heveochlorella", "Klebsormidium", "Lobomonas", "Microglena", "Monoraphidium", "Mougeotia", "Mychonastes", "Navicula", "Nemalionopsis", "Neochlorosarcina", "Neospongiococcum", "Nitzschia", "Ochromonas", "Parietochloris", "Plagioselmis", "Pleurochrysis", "Polytoma", "Protosiphon", "Pseudellipsoidion", "Rhodomonas", "Scenedesmus", "Scotinosphaera", "Sellaphora", "Spumella", "Stigeoclonium", "Symbiodinium", "Synura", "Teleaulax", "Tetracystis", "Tetranephris", "Thorea", "Trentepohlia", "Ulnaria", "Ulvella"))

## To keep only taxa of interests (for Protists) according to phylum and genus level
sub_protist_phylum_dada <- subset_taxa(phy_filtered_dada, Phylum %in% c("Cercozoa", "Choanoflagellida", "Ciliophora", "Tubulinea", "Apicomplexa", "Gracilipodida", "Protosteliida", "MAST-12", "Schizoplasmodiida"))
sub_protist_genus_dada <- subset_taxa(phy_filtered_dada, Genus %in% c("Gymnophrys", "Ichthyophthirius", "Phytophthora", "Amphileptus", "Anteholosticha", "Arcuospathidium", "Ceratiomyxella", "Cercomonas", "Chaenea", "Colpidium", "Cyrtohymena", "Dictyostelium", "Dileptus", "Flamella", "Gonostomum", "Holosticha", "Lembadion", "Litonotus", "Loxophyllum", "Monas", "Oikomonas", "Oxytricha", "Paracercomonas", "Paraphysomonas", "Phialina", "Poteriospumella", "Pseudourostyla", "Rhogostoma", "Spathidium", "Tetrahymena", "Uroleptus", "Uronema", "Vermamoeba"))

## Bar graphs of relative abundance

## Basic bar graphs based on phylum and genus level for Algae
bar_algae_phylum_dada <- microbiome::transform(sub_algae_phylum_dada, "compositional")
plot_bar(bar_algae_phylum_dada, fill="Phylum", title = "The relative abundances of the algae at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

bar_algae_genus_dada <- microbiome::transform(sub_algae_genus_dada, "compositional")
plot_bar(bar_algae_genus_dada, fill="Genus", title = "The relative abundances of the algae at the level of the genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

## Basic bar graphs based on phylum and genus level for Protists
bar_protist_phylum_dada <- microbiome::transform(sub_protist_phylum_dada, "compositional")
plot_bar(bar_protist_phylum_dada, fill="Phylum", title = "The relative abundances of protists at the level of the phylum")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

bar_protist_genus_dada <- microbiome::transform(sub_protist_genus_dada, "compositional")
plot_bar(bar_protist_genus_dada, fill="Genus", title = "The relative abundances of protists at the level of the genus")+ 
  facet_wrap(week~treatment, scales = "free_x", nrow = 1)+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

# According to Mariadassou et al., (2016); Gauthier & Derome, (2021), this statistical analysis was performed.
## Exploring the biodiversity : Alpha-diversity

## For the algae phylum level

alpha.diversity_algae_phylum_dada <- estimate_richness(sub_algae_phylum_dada, measures = c("Observed", "Shannon"))
alpha.diversity_algae_phylum_dada

alpha_algae_phylum_dada <- plot_richness(sub_algae_phylum_dada, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_phylum_dada)

pielou_algae_phylum_dada <- evenness(sub_algae_phylum_dada, 'pielou')
pielou_algae_phylum_dada

dat_lm_algae_phylum_dada <- data.frame(sample_data(sub_algae_phylum_dada))

alpha_lm_algae_phylum_dada <- cbind(alpha.diversity_algae_phylum_dada, pielou_algae_phylum_dada, dat_lm_algae_phylum_dada)
alpha_lm_algae_phylum_dada

plot_pielou_algae_phylum_dada <- ggplot(alpha_lm_algae_phylum_dada, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_phylum_dada                                        

## For the algae genus level

alpha.diversity_algae_genus_dada <- estimate_richness(sub_algae_genus_dada, measures = c("Observed", "Shannon"))
alpha.diversity_algae_genus_dada

alpha_algae_genus_dada <- plot_richness(sub_algae_genus_dada, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_algae_genus_dada)

pielou_algae_genus_dada <- evenness(sub_algae_genus_dada, 'pielou')
pielou_algae_genus_dada

dat_lm_algae_genus_dada <- data.frame(sample_data(sub_algae_genus_dada))

alpha_lm_algae_genus_dada <- cbind(alpha.diversity_algae_genus_dada, pielou_algae_genus_dada, dat_lm_algae_genus_dada)
alpha_lm_algae_genus_dada

plot_pielou_algae_genus_dada <- ggplot(alpha_lm_algae_genus_dada, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_algae_genus_dada      

## For the protist phylum level

alpha.diversity_protist_phylum_dada <- estimate_richness(sub_protist_phylum_dada, measures = c("Observed", "Shannon"))
alpha.diversity_protist_phylum_dada

alpha_protist_phylum_dada <- plot_richness(sub_protist_phylum_dada, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_phylum_dada)

pielou_protist_phylum_dada <- evenness(sub_protist_phylum_dada, 'pielou')
pielou_protist_phylum_dada

dat_lm_protist_phylum_dada <- data.frame(sample_data(sub_protist_phylum_dada))

alpha_lm_protist_phylum_dada <- cbind(alpha.diversity_protist_phylum_dada, pielou_protist_phylum_dada, dat_lm_protist_phylum_dada)
alpha_lm_protist_phylum_dada

plot_pielou_protist_phylum_dada <- ggplot(alpha_lm_protist_phylum_dada, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_phylum_dada 

## For the protist genus level

alpha.diversity_protist_genus_dada <- estimate_richness(sub_protist_genus_dada, measures = c("Observed", "Shannon"))
alpha.diversity_protist_genus_dada

alpha_protist_genus_dada <- plot_richness(sub_protist_genus_dada, color = "treatment", x = "treatment", measures = c("Observed", "Shannon"))+
  geom_boxplot(aes(fill = treatment), alpha=0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())
plot(alpha_protist_genus_dada)

pielou_protist_genus_dada <- evenness(sub_protist_genus_dada, 'pielou')
pielou_protist_genus_dada

dat_lm_protist_genus_dada <- data.frame(sample_data(sub_protist_genus_dada))

alpha_lm_protist_genus_dada <- cbind(alpha.diversity_protist_genus_dada, pielou_protist_genus_dada, dat_lm_protist_genus_dada)
alpha_lm_protist_genus_dada

plot_pielou_protist_genus_dada <- ggplot(alpha_lm_protist_genus_dada, aes(treatment,pielou)) +
  geom_boxplot(aes(fill = treatment)) +
  ylab("Pielou's evenness") +
  theme_classic()+
  stat_summary(fun.y=mean, geom = "point", shape = 5, size = 4) + scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
plot_pielou_protist_genus_dada


# Generalized Linear Modelling (GLM)
# According to ZACH (2021), Humboldt-Universität zu Berlin | Geography Department (2021) and Phillips (2018), this statistical analysis was performed.

## For the algae phylum level
Shannon_1_dada <- ceiling(alpha_lm_algae_phylum_dada$Shannon)

alpha_lm_algae_phylum_dada_mu <- alpha_lm_algae_phylum_dada %>% 
  mutate(pielou = replace_na(pielou,mean(pielou, na.rm = TRUE)))

pielou_1_dada <- ceiling(ifelse(test = alpha_lm_algae_phylum_dada_mu$pielou > 0.8966068, yes = alpha_lm_algae_phylum_dada_mu$pielou + 0.01, no = alpha_lm_algae_phylum_dada_mu$pielou))
alpha_lm_algae_phylum_1_dada <- cbind(alpha_lm_algae_phylum_dada, Shannon_1_dada, pielou_1_dada)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_phylum_1_dada$Observed)

glin_mod_algae_phylum_alpha_Observed_dada <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1_dada)
summary(glin_mod_algae_phylum_alpha_Observed_dada)
eta_squared(aov(glin_mod_algae_phylum_alpha_Observed_dada))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_algae_phylum_1_dada$Shannon_1_dada)

glin_mod_algae_phylum_alpha_Shannon_dada <- glm(Shannon_1_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1_dada)
summary(glin_mod_algae_phylum_alpha_Shannon_dada)
eta_squared(aov(glin_mod_algae_phylum_alpha_Shannon_dada))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_phylum_1_dada$pielou_1_dada)

glin_mod_algae_phylum_alpha_pielou_dada <- glm(pielou_1_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_phylum_1_dada)
summary(glin_mod_algae_phylum_alpha_pielou_dada)
eta_squared(aov(glin_mod_algae_phylum_alpha_pielou_dada))

## For the algae genus level
alpha_lm_algae_genus_dada
Shannon_2_dada <- ceiling(alpha_lm_algae_genus_dada$Shannon)

pielou_2_dada <- ceiling(ifelse(test = alpha_lm_algae_genus_dada$pielou > 0.8774688, yes = alpha_lm_algae_genus_dada$pielou + 0.01, no = alpha_lm_algae_genus_dada$pielou))
alpha_lm_algae_genus_1_dada <- cbind(alpha_lm_algae_genus_dada, Shannon_2_dada, pielou_2_dada)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_algae_genus_1_dada$Observed)

glin_mod_algae_genus_alpha_Observed_dada <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
summary(glin_mod_algae_genus_alpha_Observed_dada)
eta_squared(aov(glin_mod_algae_genus_alpha_Observed_dada))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_genus_1_dada$Shannon_2_dada)

glin_mod_algae_genus_alpha_Shannon_dada <- glm(Shannon_2_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
summary(glin_mod_algae_genus_alpha_Shannon_dada)
eta_squared(aov(glin_mod_algae_genus_alpha_Shannon_dada))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_algae_genus_1_dada$pielou_2_dada)

glin_mod_algae_genus_alpha_pielou_dada <- glm(pielou_2_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
summary(glin_mod_algae_genus_alpha_pielou_dada)
eta_squared(aov(glin_mod_algae_genus_alpha_pielou_dada))

## For the protist phylum level
alpha_lm_protist_phylum_dada
Shannon_3_dada <- ceiling(alpha_lm_protist_phylum_dada$Shannon)
pielou_3_dada <- ceiling(ifelse(test = alpha_lm_protist_phylum_dada$pielou > 0.8817381, yes = alpha_lm_protist_phylum_dada$pielou + 0.01, no = alpha_lm_protist_phylum_dada$pielou))
alpha_lm_protist_phylum_1_dada <- cbind(alpha_lm_protist_phylum_dada, Shannon_3_dada, pielou_3_dada)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_phylum_1_dada$Observed)

glin_mod_protist_phylum_alpha_Observed_dada <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1_dada)
summary(glin_mod_protist_phylum_alpha_Observed_dada)
eta_squared(aov(glin_mod_protist_phylum_alpha_Observed_dada))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_phylum_1_dada$Shannon_3_dada)

glin_mod_protist_phylum_alpha_Shannon_dada <- glm(Shannon_3_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1_dada)
summary(glin_mod_protist_phylum_alpha_Shannon_dada)
eta_squared(aov(glin_mod_protist_phylum_alpha_Shannon_dada))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_protist_phylum_1_dada$pielou_3_dada)

glin_mod_protist_phylum_alpha_pielou_dada <- glm(pielou_3_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_phylum_1_dada)
summary(glin_mod_protist_phylum_alpha_pielou_dada)
eta_squared(aov(glin_mod_protist_phylum_alpha_pielou_dada))

## For the protist genus level
alpha_lm_protist_genus_dada
Shannon_4_dada <- ceiling(alpha_lm_protist_genus_dada$Shannon)

pielou_4_dada <- ceiling(ifelse(test = alpha_lm_protist_genus_dada$pielou > 0.8974600, yes = alpha_lm_protist_genus_dada$pielou + 0.01, no = alpha_lm_protist_genus_dada$pielou))
alpha_lm_protist_genus_1_dada <- cbind(alpha_lm_protist_genus_dada, Shannon_4_dada, pielou_4_dada)

## For Observed
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_genus_1_dada$Observed)

glin_mod_protist_genus_alpha_Observed_dada <- glm(Observed ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1_dada)
summary(glin_mod_protist_genus_alpha_Observed_dada)
eta_squared(aov(glin_mod_protist_genus_alpha_Observed_dada))

## For Shannon
## The test of normality of the data by the Shapiro-Wilk Test 
shapiro.test(alpha_lm_protist_genus_1_dada$Shannon_4_dada)

glin_mod_protist_genus_alpha_Shannon_dada <- glm(Shannon_4_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1_dada)
summary(glin_mod_protist_genus_alpha_Shannon_dada)
eta_squared(aov(glin_mod_protist_genus_alpha_Shannon_dada))

## For Pielou
## The test of normality of the data by the Shapiro-Wilk Test
shapiro.test(alpha_lm_protist_genus_1_dada$pielou_4_dada)

glin_mod_protist_genus_alpha_pielou_dada <- glm(pielou_4_dada ~ nutrient*sediment*flow*time, family = poisson(link = "log"), data = alpha_lm_protist_genus_1_dada)
summary(glin_mod_protist_genus_alpha_pielou_dada)
eta_squared(aov(glin_mod_protist_genus_alpha_pielou_dada))


# Generalized Linear Mixed Models(GLMM)
## According to Walker (2021) and Bolker (2018), this statistical analysis was performed.

## For the algae phylum level

## For Observed
glmm_algae_phylum_block_alpha_Observed_ris_dada <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_algae_phylum_1_dada)
print(summary(glmm_algae_phylum_block_alpha_Observed_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_Observed_ris_dada)

plot(glmm_algae_phylum_block_alpha_Observed_ris_dada)

qqnorm(resid(glmm_algae_phylum_block_alpha_Observed_ris_dada))
qqline(resid(glmm_algae_phylum_block_alpha_Observed_ris_dada))

## For Shannon
glmm_algae_phylum_block_alpha_Shannon_ris_dada <- glmer(Shannon_1_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_phylum_1_dada)
print(summary(glmm_algae_phylum_block_alpha_Shannon_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_Shannon_ris_dada)

plot(glmm_algae_phylum_block_alpha_Shannon_ris_dada)

qqnorm(resid(glmm_algae_phylum_block_alpha_Shannon_ris_dada))
qqline(resid(glmm_algae_phylum_block_alpha_Shannon_ris_dada))

## For Pielou
glmm_algae_phylum_block_alpha_pielou_ris_dada <- glmer(pielou_1_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_phylum_1_dada)
print(summary(glmm_algae_phylum_block_alpha_pielou_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_phylum_block_alpha_pielou_ris_dada)

plot(glmm_algae_phylum_block_alpha_pielou_ris_dada)

qqnorm(resid(glmm_algae_phylum_block_alpha_pielou_ris_dada))
qqline(resid(glmm_algae_phylum_block_alpha_pielou_ris_dada))


## For the algae genus level

## For Observed
glmm_algae_genus_block_alpha_Observed_ris_dada <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
print(summary(glmm_algae_genus_block_alpha_Observed_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_Observed_ris_dada)

plot(glmm_algae_genus_block_alpha_Observed_ris_dada)

qqnorm(resid(glmm_algae_genus_block_alpha_Observed_ris_dada))
qqline(resid(glmm_algae_genus_block_alpha_Observed_ris_dada))

## For Shannon
glmm_algae_genus_block_alpha_Shannon_ris_dada <- glmer(Shannon_2_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
print(summary(glmm_algae_genus_block_alpha_Shannon_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_Shannon_ris_dada)

plot(glmm_algae_genus_block_alpha_Shannon_ris_dada)

qqnorm(resid(glmm_algae_genus_block_alpha_Shannon_ris_dada))
qqline(resid(glmm_algae_genus_block_alpha_Shannon_ris_dada))

## For Pielou
glmm_algae_genus_block_alpha_pielou_ris_dada <- glmer(pielou_2_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_algae_genus_1_dada)
print(summary(glmm_algae_genus_block_alpha_pielou_ris_dada), correlation=TRUE)
r.squaredGLMM(glmm_algae_genus_block_alpha_pielou_ris_dada)

plot(glmm_algae_genus_block_alpha_pielou_ris_dada)

qqnorm(resid(glmm_algae_genus_block_alpha_pielou_ris_dada))
qqline(resid(glmm_algae_genus_block_alpha_pielou_ris_dada))

## For the protist phylum level

## For Observed
glmm_protist_phylum_block_alpha_Observed_ris_dada <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_protist_phylum_1_dada)
print(summary(glmm_protist_phylum_block_alpha_Observed_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_Observed_ris_dada)

plot(glmm_protist_phylum_block_alpha_Observed_ris_dada)

qqnorm(resid(glmm_protist_phylum_block_alpha_Observed_ris_dada))
qqline(resid(glmm_protist_phylum_block_alpha_Observed_ris_dada))

## For Shannon
glmm_protist_phylum_block_alpha_Shannon_ris_dada <- glmer(Shannon_3_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_phylum_1_dada)
print(summary(glmm_protist_phylum_block_alpha_Shannon_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_Shannon_ris_dada)

plot(glmm_protist_phylum_block_alpha_Shannon_ris_dada)

qqnorm(resid(glmm_protist_phylum_block_alpha_Shannon_ris_dada))
qqline(resid(glmm_protist_phylum_block_alpha_Shannon_ris_dada))

## For Pielou
glmm_protist_phylum_block_alpha_pielou_ris_dada <- glmer(pielou_3_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_phylum_1_dada)
print(summary(glmm_protist_phylum_block_alpha_pielou_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_phylum_block_alpha_pielou_ris_dada)

plot(glmm_protist_phylum_block_alpha_pielou_ris_dada)

qqnorm(resid(glmm_protist_phylum_block_alpha_pielou_ris_dada))
qqline(resid(glmm_protist_phylum_block_alpha_pielou_ris_dada))


## For the protist genus level

## For Observed
glmm_protist_genus_block_alpha_Observed_ris_dada <- glmer(Observed ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data = alpha_lm_protist_genus_1_dada)
print(summary(glmm_protist_genus_block_alpha_Observed_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_Observed_ris_dada)

plot(glmm_protist_genus_block_alpha_Observed_ris_dada)

qqnorm(resid(glmm_protist_genus_block_alpha_Observed_ris_dada))
qqline(resid(glmm_protist_genus_block_alpha_Observed_ris_dada))

## For Shannon
glmm_protist_genus_block_alpha_Shannon_ris_dada <- glmer(Shannon_4_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_genus_1_dada)
print(summary(glmm_protist_genus_block_alpha_Shannon_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_Shannon_ris_dada)

plot(glmm_protist_genus_block_alpha_Shannon_ris_dada)

qqnorm(resid(glmm_protist_genus_block_alpha_Shannon_ris_dada))
qqline(resid(glmm_protist_genus_block_alpha_Shannon_ris_dada))

## For Pielou
glmm_protist_genus_block_alpha_pielou_ris_dada <- glmer(pielou_4_dada ~ nutrient*sediment*flow*time + (1|header_tanks_block), family = poisson(link = "log"), data = alpha_lm_protist_genus_1_dada)
print(summary(glmm_protist_genus_block_alpha_pielou_ris_dada), correlation = TRUE)
r.squaredGLMM(glmm_protist_genus_block_alpha_pielou_ris_dada)

plot(glmm_protist_genus_block_alpha_pielou_ris_dada)

qqnorm(resid(glmm_protist_genus_block_alpha_pielou_ris_dada))
qqline(resid(glmm_protist_genus_block_alpha_pielou_ris_dada))


# BETA-DIVERSITY INDICES

# According to Ollberding (2019), this statistical analysis was performed.
# Permutational multivariate analysis of variance (PERMANOVA)

## For the algae phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_algae_phylum_dada <- phyloseq::distance(sub_algae_phylum_dada, method = "bray")

## To make a data frame from the sample_data
perm_algae_phylum_sampledf_dada <- data.frame(sample_data(sub_algae_phylum_dada))

## To perform the adonis test
adonis2(bray_algae_phylum_dada ~ nutrient*sediment*flow*time, data = perm_algae_phylum_sampledf_dada)

## For the algae genus level
set.seed(1)

## To calculate the euclidean distance matrix
euclidean_algae_genus_dada <- phyloseq::distance(sub_algae_genus_dada, method = "euclidean")

## To make a data frame from the sample_data
perm_algae_genus_sampledf_dada <- data.frame(sample_data(sub_algae_genus_dada))

## To perform the adonis test
adonis2(euclidean_algae_genus_dada ~ nutrient*sediment*flow*time, data = perm_algae_genus_sampledf_dada)

## For the protist phylum level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_phylum_dada <- phyloseq::distance(sub_protist_phylum_dada, method = "bray")

## To make a data frame from the sample_data
perm_protist_phylum_sampledf_dada <- data.frame(sample_data(sub_protist_phylum_dada))

## To perform the adonis test
adonis2(bray_protist_phylum_dada ~ nutrient*sediment*flow*time, data = perm_protist_phylum_sampledf_dada)

## For the protist genus level
set.seed(1)

## To calculate the Bray-Curtis distance matrix
bray_protist_genus_dada <- phyloseq::distance(sub_protist_genus_dada, method = "bray")

## To make a data frame from the sample_data
perm_protist_genus_sampledf_dada <- data.frame(sample_data(sub_protist_genus_dada))

## To perform the adonis test
adonis2(bray_protist_genus_dada ~ nutrient*sediment*flow*time, data = perm_protist_genus_sampledf_dada)


# Permutational analysis of multivariate dispersion (PERMDISP)

## For the algae phylum level
dispersion_algae_phylum_dada <- vegan::betadisper(bray_algae_phylum_dada, phyloseq::sample_data(sub_algae_phylum_dada)$treatment)
dispersion_algae_phylum_dada

plot(dispersion_algae_phylum_dada, main = "For the algae phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_phylum_dada, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_phylum_dada)

## For the algae genus level
dispersion_algae_genus_dada <- vegan::betadisper(euclidean_algae_genus_dada, phyloseq::sample_data(sub_algae_genus_dada)$treatment)
dispersion_algae_genus_dada

plot(dispersion_algae_genus_dada, main = "For the algae genus level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_algae_genus_dada, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_algae_genus_dada)

## For the protist phylum level
dispersion_protist_phylum_dada <- vegan::betadisper(bray_protist_phylum_dada, phyloseq::sample_data(sub_protist_phylum_dada)$treatment)
dispersion_protist_phylum_dada

plot(dispersion_protist_phylum_dada, main = "For the protist phylum level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_phylum_dada, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_phylum_dada)

## For the protist genus level
dispersion_protist_genus_dada <- vegan::betadisper(bray_protist_genus_dada, phyloseq::sample_data(sub_protist_genus_dada)$treatment)
dispersion_protist_genus_dada

plot(dispersion_protist_genus_dada, main = "For the protist genus level: Ordination Centroids and Dispersion Labeled for Different Treatments: Aitchison Distance", sub = "")

boxplot(dispersion_protist_genus_dada, main = "", xlab = "")

## Permutation test for homogeneity of multivariate dispersions
permutest(dispersion_protist_genus_dada)


# According to Xu & Yu (2020), this statistical analysis was performed.
# The PCoA analysis (Principal Coordinate Analysis)
library(MicrobiotaProcess)

## For the algae phylum level
pcoares_algae_phylum_dada <- get_pcoa(obj=sub_algae_phylum_dada, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_algae_phylum_dada <- ggordpoint(obj=pcoares_algae_phylum_dada, biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_algae_phylum_dada <- ggordpoint(obj=pcoares_algae_phylum_dada, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                     factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_algae_phylum_dada | pcoaplot2_algae_phylum_dada

## For the algae genus level
pcoares_algae_genus_dada <- get_pcoa(obj=sub_algae_genus_dada, distmethod="euclidean", method="hellinger")
## For visualizing the result
pcoaplot1_algae_genus_dada <- ggordpoint(obj=pcoares_algae_genus_dada, biplot=TRUE, speciesannot=TRUE,
                                    factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_algae_genus_dada <- ggordpoint(obj=pcoares_algae_genus_dada, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                    factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_algae_genus_dada | pcoaplot2_algae_genus_dada

## For the protist phylum level
pcoares_protist_phylum_dada <- get_pcoa(obj=sub_protist_phylum_dada, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_protist_phylum_dada <- ggordpoint(obj=pcoares_protist_phylum_dada, biplot=TRUE, speciesannot=TRUE,
                                       factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_protist_phylum_dada <- ggordpoint(obj=pcoares_protist_phylum_dada, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                       factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_protist_phylum_dada | pcoaplot2_protist_phylum_dada

## For the protist genus level
pcoares_protist_genus_dada <- get_pcoa(obj=sub_protist_genus_dada, distmethod="bray", method="hellinger")
## For visualizing the result
pcoaplot1_protist_genus_dada <- ggordpoint(obj=pcoares_protist_genus_dada, biplot=TRUE, speciesannot=TRUE,
                                      factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
## The first and third principal co-ordinates
pcoaplot2_protist_genus_dada <- ggordpoint(obj=pcoares_protist_genus_dada, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                                      factorNames=c("treatment"), ellipse=TRUE) +
  scale_color_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green")) +
  scale_fill_manual(values=c("#D55E00", "#F0E442", "#E69F00", "#0072B2", "#999999", "#56B4E9", "red", "green"))
pcoaplot1_protist_genus_dada | pcoaplot2_protist_genus_dada


# According to joey711 (2022), the most abundant taxa was calculated. 
## To find the most abundant taxa
## Function
most_abundant_taxa_dada <- function(x,taxa){
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
abundant_algae_phylum_dada <- most_abundant_taxa_dada(sub_algae_phylum_dada,"Phylum")
abundant_algae_phylum_dada

abundant_algae_genus_dada <- most_abundant_taxa_dada(sub_algae_genus_dada,"Genus")
abundant_algae_genus_dada

abundant_protist_phylum_dada <- most_abundant_taxa_dada(sub_protist_phylum_dada,"Phylum")
abundant_protist_phylum_dada

abundant_protist_genus_dada <- most_abundant_taxa_dada(sub_protist_genus_dada,"Genus")
abundant_protist_genus_dada

# Multivariate Analysis Of Variance (MANOVA)
## According to Kassambara (2017); Ben-Shachar et al., (2020), this statistical analysis was performed.

## For the most abundant taxa of the algal phylum level
manova_df_algae_phylum_dada = as(sample_data(sub_algae_phylum_dada), "data.frame")
manova_df_algae_phylum_dada
manova_df_algae_phylum_1_dada = cbind(manova_df_algae_phylum_dada, abundant_algae_phylum_dada)
manova_df_algae_phylum_1_dada

### For Ochrophyta and Phragmoplastophyta
manova_abundant_algae_phylum_dada <- manova(cbind(Ochrophyta, Phragmoplastophyta) ~ nutrient*sediment*flow*time, data = manova_df_algae_phylum_1_dada)
manova_abundant_algae_phylum_dada
summary.aov(manova_abundant_algae_phylum_dada)
eta_squared(aov(manova_abundant_algae_phylum_dada))

## For the most abundant taxa of the algal genus level
manova_df_algae_genus_dada = as(sample_data(sub_algae_genus_dada), "data.frame")
manova_df_algae_genus_dada
manova_df_algae_genus_1_dada = cbind(manova_df_algae_genus_dada, abundant_algae_genus_dada)
manova_df_algae_genus_1_dada

### For Cymbella, Desmodesmus, Achnanthidium, Spumella, Chlorella, Pseudellipsoidion, Navicula, Chlamydomonas, Monoraphidium and Rhodomonas
manova_abundant_algae_genus_dada <- manova(cbind(Cymbella, Desmodesmus, Achnanthidium, Spumella, Chlorella, Pseudellipsoidion, Navicula, Chlamydomonas, Monoraphidium, Rhodomonas) ~ nutrient*sediment*flow*time, data = manova_df_algae_genus_1_dada)
manova_abundant_algae_genus_dada
summary.aov(manova_abundant_algae_genus_dada)
eta_squared(aov(manova_abundant_algae_genus_dada))

## For the most abundant taxa of the protist phylum level
manova_df_protist_phylum_dada = as(sample_data(sub_protist_phylum_dada), "data.frame")
manova_df_protist_phylum_dada
manova_df_protist_phylum_1_dada = cbind(manova_df_protist_phylum_dada, abundant_protist_phylum_dada)
manova_df_protist_phylum_1_dada

### For Cercozoa, Apicomplexa and Tubulinea
manova_abundant_protist_phylum_dada <- manova(cbind(Cercozoa, Apicomplexa, Tubulinea) ~ nutrient*sediment*flow*time, data = manova_df_protist_phylum_1_dada)
manova_abundant_protist_phylum_dada
summary.aov(manova_abundant_protist_phylum_dada)
eta_squared(aov(manova_abundant_protist_phylum_dada))

## For the most abundant taxa of the protist genus level
manova_df_protist_genus_dada = as(sample_data(sub_protist_genus_dada), "data.frame")
manova_df_protist_genus_dada
manova_df_protist_genus_1_dada = cbind(manova_df_protist_genus_dada, abundant_protist_genus_dada)
manova_df_protist_genus_1_dada

### For Rhogostoma, Dictyostelium, Vermamoeba, Gymnophrys, Paracercomonas and Anteholosticha  
manova_abundant_protist_genus_dada <- manova(cbind(Rhogostoma, Dictyostelium, Vermamoeba, Gymnophrys, Paracercomonas, Anteholosticha) ~ nutrient*sediment*flow*time, data = manova_df_protist_genus_1_dada)
manova_abundant_protist_genus_dada
summary.aov(manova_abundant_protist_genus_dada)
eta_squared(aov(manova_abundant_protist_genus_dada))


# Generalized linear latent variable models(GLLVMs)
## According to Niku et al., (2019), this statistical analysis was performed.

## For the algae phylum level
sap_asv_dada <- data.frame(otu_table(sub_algae_phylum_dada))
sap_asv_dada

sap_tax_dada <- data.frame(tax_table(sub_algae_phylum_dada))
sap_tax_dada

sap_sample_dada <- data.frame(sample_data(sub_algae_phylum_dada))
sap_sample_dada

sap_tax_asv_dada <- cbind(sap_tax_dada["Phylum"], sap_asv_dada)
sap_tax_asv_dada

rownames(sap_tax_asv_dada) <- NULL
sap_tax_asv_dada

sap_tax_asv_uq_dada <- sap_tax_asv_dada %>% distinct(Phylum, .keep_all = TRUE)
sap_tax_asv_uq_dada

sap_tax_asv_index_dada <- sap_tax_asv_uq_dada %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Phylum')
sap_tax_asv_index_dada

sap_tax_asv_final_dada <- as.data.frame(t(sap_tax_asv_index_dada))
sap_tax_asv_final_dada

## The ordination as model based
mydata_fitp_ap_dada <- gllvm(sap_tax_asv_final_dada, family = poisson(link = "log"), starting.val="random")
summary(mydata_fitp_ap_dada)

mydata_fit_ord_nb_ap_dada <- gllvm(sap_tax_asv_final_dada, family = "negative.binomial", starting.val="random")
summary(mydata_fit_ord_nb_ap_dada)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_ap_dada, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_ap_dada, var.colors = 1)

ordiplot(mydata_fitp_ap_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_ap_dada, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_ap_dada$params$theta) <- paste("spp", 1:ncol(mydata_fitp_ap_dada$y), sep = "")
ordiplot(mydata_fitp_ap_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_ap_criterias_p_dada <- NULL
for(i in 1:5){
  fiti_ap_p_dada <- gllvm(sap_tax_asv_final_dada, sap_sample_dada, family = poisson(link = "log"), num.lv = 1, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_ap_criterias_p_dada[i + 1] <- summary(fiti_ap_p_dada)$AICc
  names(mydata_ap_criterias_p_dada)[i + 1] = i
}
## To compare the AICc values
mydata_ap_criterias_p_dada

mydata_fit_env_p_ap_dada <- gllvm(sap_tax_asv_final_dada, sap_sample_dada, family = poisson(link = "log"), num.lv = 1, 
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_ap_dada)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between algal phylums can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_ap_p_dada <- getResidualCor(mydata_fit_env_p_ap_dada)

corrplot(mydata_cr_ap_p_dada[order.single(mydata_cr_ap_p_dada), order.single(mydata_cr_ap_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables
mydata_fit4lv_ap_p_dada <- gllvm(sap_tax_asv_final_dada, family = poisson(link = "log"), num.lv = 1, seed = 1234)
summary(mydata_fit4lv_ap_p_dada)

## The correlation matrix
mydata_cr0_ap_p_dada <- getResidualCor(mydata_fit4lv_ap_p_dada)
corrplot(mydata_cr0_ap_p_dada[order.single(mydata_cr0_ap_p_dada), order.single(mydata_cr0_ap_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual bi plot can be used to visualize the residual correlations
mydata_fit_env2_ap_p_dada <- gllvm(sap_tax_asv_final_dada, sap_sample_dada, family = poisson(link = "log"), num.lv = 1,  
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_ap_p_dada)

rownames(mydata_fit_env2_ap_p_dada$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_ap_p_dada$y), sep = "")
ordiplot(mydata_fitp_ap_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_ap_p_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_ap_p_dada <- getResidualCov(mydata_fit_env_p_ap_dada, adjust = 0)
mydata_rcov0_ap_p_dada <- getResidualCov(mydata_fit4lv_ap_p_dada, adjust = 0)
mydata_rcov0_ap_p_dada$trace; mydata_rcov_ap_p_dada$trace

1 - mydata_rcov_ap_p_dada$trace / mydata_rcov0_ap_p_dada$trace


## For the algae genus level
sag_asv_dada <- data.frame(otu_table(sub_algae_genus_dada))
sag_asv_dada

sag_tax_dada <- data.frame(tax_table(sub_algae_genus_dada))
sag_tax_dada

sag_sample_dada <- data.frame(sample_data(sub_algae_genus_dada))
sag_sample_dada

sag_tax_asv_dada <- cbind(sag_tax_dada["Genus"], sag_asv_dada)
sag_tax_asv_dada

rownames(sag_tax_asv_dada) <- NULL
sag_tax_asv_dada

sag_tax_asv_uq_dada <- sag_tax_asv_dada %>% distinct(Genus, .keep_all = TRUE)
sag_tax_asv_uq_dada

sag_tax_asv_index_dada <- sag_tax_asv_uq_dada %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Genus')
sag_tax_asv_index_dada

sag_tax_asv_final_dada <- as.data.frame(t(sag_tax_asv_index_dada))
sag_tax_asv_final_dada

## The ordination as model based
mydata_fitp_ag_dada <- gllvm(sag_tax_asv_final_dada, family = poisson(link = "log"))
summary(mydata_fitp_ag_dada)

mydata_fit_ord_nb_ag_dada <- gllvm(sag_tax_asv_final_dada, family = "negative.binomial")
summary(mydata_fit_ord_nb_ag_dada)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_ag_dada, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_ag_dada, var.colors = 1)

ordiplot(mydata_fitp_ag_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_ag_dada, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_ag_dada$params$theta) <- paste("spp", 1:ncol(mydata_fitp_ag_dada$y), sep = "")
ordiplot(mydata_fitp_ag_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_ag_criterias_p_dada <- NULL
for(i in 1:5){
  fiti_ag_p_dada <- gllvm(sag_tax_asv_final_dada, sag_sample_dada, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_ag_criterias_p_dada[i + 1] <- summary(fiti_ag_p_dada)$AICc
  names(mydata_ag_criterias_p_dada)[i + 1] = i
}
## To compare the AICc values
mydata_ag_criterias_p_dada

mydata_fit_env_p_ag_dada <- gllvm(sag_tax_asv_final_dada, sag_sample_dada, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_ag_dada)

## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between algal genus can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_ag_p_dada <- getResidualCor(mydata_fit_env_p_ag_dada)

corrplot(mydata_cr_ag_p_dada[order.single(mydata_cr_ag_p_dada), order.single(mydata_cr_ag_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables
mydata_fit4lv_ag_p_dada <- gllvm(sag_tax_asv_final_dada, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_ag_p_dada)

## The correlation matrix
mydata_cr0_ag_p_dada <- getResidualCor(mydata_fit4lv_ag_p_dada)
corrplot(mydata_cr0_ag_p_dada[order.single(mydata_cr0_ag_p_dada), order.single(mydata_cr0_ag_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual bi plot can be used to visualize the residual correlations
mydata_fit_env2_ag_p_dada <- gllvm(sag_tax_asv_final_dada, sag_sample_dada, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_ag_p_dada)

rownames(mydata_fit_env2_ag_p_dada$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_ag_p_dada$y), sep = "")
ordiplot(mydata_fitp_ag_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_ag_p_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_ag_p_dada <- getResidualCov(mydata_fit_env_p_ag_dada, adjust = 0)
mydata_rcov0_ag_p_dada <- getResidualCov(mydata_fit4lv_ag_p_dada, adjust = 0)
mydata_rcov0_ag_p_dada$trace; mydata_rcov_ag_p_dada$trace

1 - mydata_rcov_ag_p_dada$trace / mydata_rcov0_ag_p_dada$trace

## For the protist phylum level
spp_asv_dada <- data.frame(otu_table(sub_protist_phylum_dada))
spp_asv_dada

spp_tax_dada <- data.frame(tax_table(sub_protist_phylum_dada))
spp_tax_dada

spp_sample_dada <- data.frame(sample_data(sub_protist_phylum_dada))
spp_sample_dada

spp_tax_asv_dada <- cbind(spp_tax_dada["Phylum"], spp_asv_dada)
spp_tax_asv_dada

rownames(spp_tax_asv_dada) <- NULL
spp_tax_asv_dada

spp_tax_asv_uq_dada <- spp_tax_asv_dada %>% distinct(Phylum, .keep_all = TRUE)
spp_tax_asv_uq_dada

spp_tax_asv_index_dada <- spp_tax_asv_uq_dada %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Phylum')
spp_tax_asv_index_dada

spp_tax_asv_final_dada <- as.data.frame(t(spp_tax_asv_index_dada))
spp_tax_asv_final_dada

## The ordination as model based
mydata_fitp_pp_dada <- gllvm(spp_tax_asv_final_dada, family = poisson(link = "log"))
summary(mydata_fitp_pp_dada)

mydata_fit_ord_nb_pp_dada <- gllvm(spp_tax_asv_final_dada, family = "negative.binomial")
summary(mydata_fit_ord_nb_pp_dada)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_pp_dada, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_pp_dada, var.colors = 1)

ordiplot(mydata_fitp_pp_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_pp_dada, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_pp_dada$params$theta) <- paste("spp", 1:ncol(mydata_fitp_pp_dada$y), sep = "")
ordiplot(mydata_fitp_pp_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_pp_criterias_p_dada <- NULL
for(i in 1:5){
  fiti_pp_p_dada <- gllvm(spp_tax_asv_final_dada, spp_sample_dada, family = poisson(link = "log"), num.lv = 4, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_pp_criterias_p_dada[i + 1] <- summary(fiti_pp_p_dada)$AICc
  names(mydata_pp_criterias_p_dada)[i + 1] = i
}
## To compare the AICc values
mydata_pp_criterias_p_dada

mydata_fit_env_p_pp_dada <- gllvm(spp_tax_asv_final_dada, spp_sample_dada, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_pp_dada)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between protist phylums can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_pp_p_dada <- getResidualCor(mydata_fit_env_p_pp_dada)

corrplot(mydata_cr_pp_p_dada[order.single(mydata_cr_pp_p_dada), order.single(mydata_cr_pp_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables
mydata_fit4lv_pp_p_dada <- gllvm(spp_tax_asv_final_dada, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_pp_p_dada)

## The correlation matrix
mydata_cr0_pp_p_dada <- getResidualCor(mydata_fit4lv_pp_p_dada)
corrplot(mydata_cr0_pp_p_dada[order.single(mydata_cr0_pp_p_dada), order.single(mydata_cr0_pp_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual bi plot can be used to visualize the residual correlations
mydata_fit_env2_pp_p_dada <- gllvm(spp_tax_asv_final_dada, spp_sample_dada, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_pp_p_dada)

rownames(mydata_fit_env2_pp_p_dada$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_pp_p_dada$y), sep = "")
ordiplot(mydata_fitp_pp_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_pp_p_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_pp_p_dada <- getResidualCov(mydata_fit_env_p_pp_dada, adjust = 0)
mydata_rcov0_pp_p_dada <- getResidualCov(mydata_fit4lv_pp_p_dada, adjust = 0)
mydata_rcov0_pp_p_dada$trace; mydata_rcov_pp_p_dada$trace

1 - mydata_rcov_pp_p_dada$trace / mydata_rcov0_pp_p_dada$trace

## For the protist genus level
spg_asv_dada <- data.frame(otu_table(sub_protist_genus_dada))
spg_asv_dada

spg_tax_dada <- data.frame(tax_table(sub_protist_genus_dada))
spg_tax_dada

spg_sample_dada <- data.frame(sample_data(sub_protist_genus_dada))
spg_sample_dada

spg_tax_asv_dada <- cbind(spg_tax_dada["Genus"], spg_asv_dada)
spg_tax_asv_dada

rownames(spg_tax_asv_dada) <- NULL
spg_tax_asv_dada

spg_tax_asv_uq_dada <- spg_tax_asv_dada %>% distinct(Genus, .keep_all = TRUE)
spg_tax_asv_uq_dada

spg_tax_asv_index_dada <- spg_tax_asv_uq_dada %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Genus')
spg_tax_asv_index_dada

spg_tax_asv_final_dada <- as.data.frame(t(spg_tax_asv_index_dada))
spg_tax_asv_final_dada

## The ordination as model based
mydata_fitp_pg_dada <- gllvm(spg_tax_asv_final_dada, family = poisson(link = "log"))
summary(mydata_fitp_pg_dada)

mydata_fit_ord_nb_pg_dada <- gllvm(spg_tax_asv_final_dada, family = "negative.binomial")
summary(mydata_fit_ord_nb_pg_dada)

## To plot residuals for the Poisson model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fitp_pg_dada, var.colors = 1)

## To plot residuals for the NB model
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(mydata_fit_ord_nb_pg_dada, var.colors = 1)

ordiplot(mydata_fitp_pg_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(mydata_fitp_pg_dada, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

rownames(mydata_fitp_pg_dada$params$theta) <- paste("spp", 1:ncol(mydata_fitp_pg_dada$y), sep = "")
ordiplot(mydata_fitp_pg_dada, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-2, 1.6), 
         main = "Biplot", jitter = TRUE, cex.spp = 0.8)

## For modeling with the environmental variables
mydata_pg_criterias_p_dada <- NULL
for(i in 1:5){
  fiti_pg_p_dada <- gllvm(spg_tax_asv_final_dada, spg_sample_dada, family = poisson(link = "log"), num.lv = i, sd.errors = FALSE,
                     formula = ~ nutrient*sediment*flow*time, seed = 1234)
  mydata_pg_criterias_p_dada[i + 1] <- summary(fiti_pg_p_dada)$AICc
  names(mydata_pg_criterias_p_dada)[i + 1] = i
}
## To compare the AICc values
mydata_pg_criterias_p_dada

mydata_fit_env_p_pg_dada <- gllvm(spg_tax_asv_final_dada, spg_sample_dada, family = poisson(link = "log"), num.lv = 4,
                             formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env_p_pg_dada)


## To study the co-occurrence patterns
## The correlation of the response variables may be induced by the latent variables.
## In this way, correlation patterns between protist genus can also be estimated.
## The extents can also be described by the manipulative environmental factors.

## The residual correlation matrix:
mydata_cr_pg_p_dada <- getResidualCor(mydata_fit_env_p_pg_dada)

corrplot(mydata_cr_pg_p_dada[order.single(mydata_cr_pg_p_dada), order.single(mydata_cr_pg_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## To fit GLLVM without environmental variables
mydata_fit4lv_pg_p_dada <- gllvm(spg_tax_asv_final_dada, family = poisson(link = "log"), num.lv = 4, seed = 1234)
summary(mydata_fit4lv_pg_p_dada)

## The correlation matrix
mydata_cr0_pg_p_dada <- getResidualCor(mydata_fit4lv_pg_p_dada)
corrplot(mydata_cr0_pg_p_dada[order.single(mydata_cr0_pg_p_dada), order.single(mydata_cr0_pg_p_dada)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

## The residual bi plot can be used to visualize the residual correlations
mydata_fit_env2_pg_p_dada <- gllvm(spg_tax_asv_final_dada, spg_sample_dada, family = poisson(link = "log"), num.lv = 2, 
                              formula = ~ nutrient*sediment*flow*time, seed = 1234)
summary(mydata_fit_env2_pg_p_dada)

rownames(mydata_fit_env2_pg_p_dada$params$theta) <- paste("sp", 1:ncol(mydata_fit_env2_pg_p_dada$y), sep = "")
ordiplot(mydata_fitp_pg_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(a)")
ordiplot(mydata_fit_env2_pg_p_dada, biplot = TRUE, ind.spp = 15, jitter = TRUE, cex.spp = 1,
         xlim = c(-4, 3.5), ylim = c(-2.5, 2), main = "(b)")

## The amount of variation in the data associated to environmental manipulative factors can be calculated by the getResidualCov() function 
mydata_rcov_pg_p_dada <- getResidualCov(mydata_fit_env_p_pg_dada, adjust = 0)
mydata_rcov0_pg_p_dada <- getResidualCov(mydata_fit4lv_pg_p_dada, adjust = 0)
mydata_rcov0_pg_p_dada$trace; mydata_rcov_pg_p_dada$trace

1 - mydata_rcov_pg_p_dada$trace / mydata_rcov0_pg_p_dada$trace














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

Xu, S., & Yu, G. (2020, Nov 24). Workshop of microbiome dataset analysis using MicrobiotaProcess. MicrobiotaProcessWorkshop 0.0.0.92. https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

joey711. (2022). Find the Most abundant Taxa in individual samples #847. https://github.com/joey711/phyloseq/issues/847

Kassambara, A. (2017). MANOVA Test in R: Multivariate Analysis of Variance. Statistical tools for high-throughput data analysis. STHDA. http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance#infos

Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815

Niku, J., Hui, F. K., Taskinen, S., & Warton, D. I. (2019). gllvm: Fast analysis of multivariate abundance data with generalized linear latent variable models in r. Methods in Ecology and Evolution, 10(12), 2173-2182.

Niku J, Hui FKC, Taskinen S, Warton DI (2019). "gllvm - Fast analysis of multivariate abundance data with generalized linear latent variable models in R." _Methods in Ecology and Evolution_, *10*, 2173-2182.



