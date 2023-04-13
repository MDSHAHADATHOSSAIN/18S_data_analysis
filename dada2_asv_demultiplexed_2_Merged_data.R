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







## List of references

## Callahan, Benjamin. (2018). Silva taxonomic training data formatted for DADA2 (Silva version 132) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.1172783

## Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." _Nature Methods_, *13*, 581-583. doi:10.1038/nmeth.3869 <https://doi.org/10.1038/nmeth.3869>.

## Callahan, Benjamin. (2018). DADA2 Pipeline Tutorial (1.8). dada2. https://benjjneb.github.io/dada2/tutorial_1_8.html

