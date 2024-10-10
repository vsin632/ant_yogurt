#####################################
# author: Veronica M. Sinotte
# title: Ant yogurt 16S data analysis
# date: 2024-Oct-9
# output: R file
# R version: "R version 4.3.1 (2023-06-16)"
# manuscript: Making yogurt with the ant holobiont reveals bacteria, acids, and enzymes for food fermentation
# note: co-author VRV also contributed to the analysis, and some final figures were made by her in GraphPad Prisim, these cases are noted
######################################

######################################
#### 16S metabarcoding analysis
######################################

#this includes filter, trim, align, taxa assignment, cleaning up data including removing contaminants, mitochondira, chlorplasts
# further assessment of mock community positive contorl, pruning some reads, examining and removing positive yogurt/supernatant controls for future analyses

######################
#Notes of Sequencing:

# 250 bp pair-end sequencing, with primer and index removed about 220 bp forward and reverse
# Primers 341F, 806R
# amplicon size 465 bp, primer and index removed amplicon approx. 430 bp
# overlap 35 bp 

####################
#Load Libraries and Directory
#may require downloads through Bioconductor

#clear workspace
rm(list=ls())

#load libraries
library(dada2)
packageVersion("dada2")
#version 1.30.0
library(phyloseq)
packageVersion("phyloseq")
#version ‘1.46.0’
library(Biostrings)
packageVersion("BioStrings")
#version ‘2.70.1’
library(ggplot2)
packageVersion("ggplot2")
#version ‘3.4.4’
library(tidyverse)
packageVersion("tidyverse")
#version ‘2.0.0’
#install.packages("DECIPHER")
library(DECIPHER)
packageVersion("DECIPHER")
#version ‘2.0.0’
library("ShortRead")
packageVersion("ShortRead")
#version ‘1.60.0’
library(genefilter)
packageVersion("genefilter")
library(ComplexHeatmap)
library(microbiome)
library(microViz)
library(xlsx)
library(vegan)


#Follows the pipeline described in Callahan et al. 2016
#https://f1000research.com/articles/5-1492
#tutorial here
#https://benjjneb.github.io/dada2/tutorial.html

#set path to data files
path<- "/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes"

#check files in path
list.files(path)
#files look good
#total sample is 80

##################
#Identify Forward and Reverse

#perform string manipulation to get forward and reverse fastq.gz files
#files have the format
#forward primer: SampleX_1.fastq.gz
#reverse primer: SampleX_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

#extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####################
#Inspect Read Profiles

#visualize forward reads 
#note numbers in brackets refer to samples order in the data frame
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnFs[11:20])
plotQualityProfile(fnFs[21:30])
plotQualityProfile(fnFs[31:40])
plotQualityProfile(fnFs[41:50])
plotQualityProfile(fnFs[51:60])
plotQualityProfile(fnFs[61:70])
plotQualityProfile(fnFs[71:80])
#Notes: very high quality overall for Illumina, seems every now and there is a base with very low quality after 150 bp
#Trim: not even really needed
#Lower sequencing quality: A10, A11, A26, A24, A4, A7, B64
#OBS: read count low for: 
# A24 (19,121 reads)
# A18 (40,757 reads)
# A6 (6,671 reads)
# A75 (9,704 reads)
# A77 (4,944 reads)
# A80 (27,797 reads)
# A86 (29,914 reads)

#visualize reverse reads reads 
#note numbers in brackets refer to samples order in the data frame
plotQualityProfile(fnRs[1:10])
plotQualityProfile(fnRs[11:20])
plotQualityProfile(fnRs[21:30])
plotQualityProfile(fnRs[31:40])
plotQualityProfile(fnRs[41:50])
plotQualityProfile(fnRs[51:60])
plotQualityProfile(fnRs[61:70])
plotQualityProfile(fnRs[71:80])
#Notes: again super high quality for Illumina, may good until 250 bp, some have quality drop to around QS30 at 180 bp
#Trim:Not necessarily needed
#Lower sequencing quality after 180: A10, A11, A110, A112, A12, A2, A26, A4, A53, A57, A52, A67, A71, A83, B63. B64
# the lower sequencing quality is only that detailed in the notes above, not bad
#OBS: read count low for: 
# A24 (19,121 reads)
# A18 (40,757 reads)
# A6 (6,671 reads)
# A75 (9,704 reads)
# A77 (4,944 reads)
# A80 (27,797 reads)
# A86 (29,914 reads)


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


################
#Filter and Trim

#parameters explained here
#truncLen: the amount of bases you leave from F and R reads
# your truncLen (trimming) must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.
#maxN: Maximum number of N bases, dada doesn't allow Ns
#maxEE: reads with higher expected errors than X will be discarded, first number refers to forward. second number to reverse
#trunQ: truncate reads at first instance of quality less than X
# note modification maxEE=c(2,5) allows more errors on reverse reads. the bigger the number the longer it takes

#decided not to trim reads given the quite high quality of the sequences, and the small overlap of about 30 bp
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

#take a look at how the filtering went
head(out)
list(out)
#seems to take out approx 1/4-1/5 of sequences
colSums(out)
# reads.in reads.out 
# 4871494   3907637 

#################
#Learn Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
#plots look fine
#The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
#The red line shows the error rates expected under the nominal definition of the Q-score. 
#The estimated error rates (black line) should be good fit to the observed rates (points), and the error rates drop with increased quality as expected. 

#######################
#Sample Inference

#pooling can help identify rare sequence variants, but regardless >99% of reads are commonly assigned. 
#pooling will at least double the computing time
#pooling may be considered when the sequencing quality is poor.
#therefore, pooling not really necessary for this data set
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]


#minOverlap = 20 is ideal, but tried but didn't work out well
#using minOverlap of 16 which seems to be the highest that still works well, the standard for the dada2 pipeline is 12
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 16)

#check out a few of the samples
head(mergers[[1]])
head(mergers[[5]])
#looks good

############################
# Construct Sequence Table

seqtab <- makeSequenceTable(mergers)
#check out data frame: first number corresponds to number of samples second number corresponds to number of ASVs
dim(seqtab)
# samples: 80, ASVs: 9237

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))
#looks good, but few really short, should be 430 bp amplicon. Shorter reads potentially due to non-specific priming. 
#Change according to histogram for amplicons 390-440 bp.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 390:440]
dim(seqtab2)
#samples: 80, ASVs: 8738
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)))
#looks good!

##############################
# Remove Chimeras

#note change reference to seqtab2 per filtering in last step
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# samples:80, ASVs: 5504
# losing this many ASVs to bimeras is normal
sum(seqtab.nochim)/sum(seqtab2)
# remaining seq are 0.9871, therefore the bimeras were less than 2% of reads
sum(seqtab.nochim)
# total reads in dataset: 3,808,507

###################################
# Track reads through the pipeline

#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
#We should be all good here as we have checked our progress closely in each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#save as CSV
write.csv(track,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/output_files/tracked_reads.csv")

####################################
# Assign taxonomy

# the package uses naive Bayesian classifier method for this step
#Silva release 138.1, from March 2021, this is the latest release

#assign taxonomy of dataset WITH NEGATIVE CONTROLS
taxa <- assignTaxonomy(seqtab.nochim, "/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes//tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
#add species-level classification 
taxa <- addSpecies(taxa, "/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/tax/silva_species_assignment_v138.1.fa.gz")

#inspect the taxonomic assignments
taxa.print <- taxa 
# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#looks correct, top genera include Wolbachia, Fructilactobacillus, Neokomagataea, Bacillus

####################################
#Save all the files
#save all files
save(seqtab.nochim, seqtab2, mergers, dadaFs, dadaRs, errR, errR, filtFs, filtRs, sample.names, taxa, file ="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/antyogurt1_dada2.RData")
#load("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/antyogurt1_dada2.RData")

############################################
#SCRuB decontamination based on negative controls

#using the recently published package SCRuB
#SCRuB uses a simple metadata matrix with the sample names, 'is_control' which is TRUE or FALSE, 'sample_type', and the 'sample_well'
#we do not have the sample well because samples were processed as individual eppendorfs and sequenced at novogene (no plate data)
#SCRub then uses a data matrix that is a count matrix, based on the most granular phylogentetic level (ie ASV)
#therefore we make a simple metadata table here, and can use the seqtab.nochim before assigning taxonomy

#unlike the R package decontam, we do not need qPCR data! 
#based on the recent publication for SCRuB, it performs better than decontam, whether or not well data is included

#https://doi.org/10.1038/s41587-023-01696-w
#https://github.com/Shenhav-and-Korem-labs/SCRuB

#install and library necessary packages

#install.packages( c('glmnet', 'torch') )
library(glmnet)
library(torch)
#install.packages('devtools')
library(devtools)
#devtools::install_github("shenhav-and-korem-labs/SCRuB")
library(SCRuB)

#the 'data' here will be seqtab.nochim
dim(seqtab.nochim)
data<- seqtab.nochim
# sample/rows: 80, columns/ASVs: 5504

#read in metadata
SCRuB.metadata<- read_csv("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/SCRuB/Veronica_SCRuB/SCRuB_metadata.csv")

#check out metadata
str(SCRuB.metadata)
dim(SCRuB.metadata)
#make the column 1 the row names for the metadata
metadata1<-SCRuB.metadata %>% remove_rownames %>% column_to_rownames(var="...1")
#remove well locations, it isn't necessary to include as there is not data
metadata2<- subset(metadata1, select = c(is_control, sample_type) )

#note that the following samples have low read counts (below values are before pairing and filtering): 
# A24 (19,121 reads)
# A18 (40,757 reads)
# A6 (6,671 reads)
# A75 (9,704 reads)
# A77 (4,944 reads)
# A80 (27,797 reads)

#note the following samples are negative controls
#A73
#A74
#A75
#A76
#A77

#try out SCRuB
scr_out <- SCRuB(data, metadata2, control_order = NA)

#evaluate SCRuB outputs
#Let's take a look at the results. The estimated level of contamination was very low, as the fitted `p` parameters indicate our samples are close to `3%` contamination.
# 'p' is the estimated proportion of the sample that was not contamination
#look at proportion remaining for all samples
scr_out$p
#looks okay, most retain 0.95 for more, a few less but they were probelmatic samples, supernatant controls, etc
#look at box plot to confirm
scr_out$p %>% boxplot()

#make data frame with the proportion of reads remaining from all the samples (this excludes the 5 neg controls)
abund_remaining<-as.data.frame(scr_out$p)
abund_remaining1<-rownames_to_column(abund_remaining, var = "SampleID")

#make data frame of decontaminated samples in the same format as seqtab.nochim
decont_samp<- as.data.frame(scr_out$decontaminated_samples)
#great, 75 rows with 5504 ASVs

#lets get the estimated relative abundances of the contaminants, gamma
contam_abund<- as.data.frame(scr_out$inner_iterations$neg_control$gamma)
#plot it quick to see what the composition is
scr_out$inner_iterations$`neg_control`$gamma %>% plot()
#okay, its quite diverse, one ASV makes up .2 percent of the contaminants, others much less
#with quick BLAST online, can see that the main contaminant at 0.2 is Streptococcus, second most abundant at .12 is Vibrio. 

#save proportion of samples remaining after decontamination
write.csv(abund_remaining1,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/SCRuB/Veronica_SCRuB/Decontam_sample_abundances.csv")
#save ASV and read count dataframe for the decontaminated samples
write.csv(decont_samp,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/Decontam_sample_abundances.csv")
#save the estimated relative abundances of contaminants
write.csv(contam_abund,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/Contaminant_abundances.csv")
write.csv(contam_abund,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/SCRuB/Veronica_SCRuB/Contaminant_abundances.csv")

#have a quick look at the correlation between the bacterial load of samples form qPCR and the proportion remaining
decont_plot<- read_csv("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/decontam_sample_metadata.csv")

#you may predict that the lower the bacterial load, the higher the decontam abundance 
ggplot(decont_plot, aes(x = Ng_Bact_DNA, y = Decontam_abundance)) +
  geom_point()+
  geom_smooth(method=lm)

#great variation in ng DNA, try to log transform it to get a better idea
ggplot(decont_plot, aes(x = Ng_Bact_DNA, y = Decontam_abundance)) +
  geom_point()+
  scale_x_log10()+
  geom_smooth(method=lm)
#OK looks good, generally this is what comes out of Decontam R package pipeline

ps_contaminants<-phyloseq(otu_table(contam_abund, taxa_are_rows=TRUE), 
                          tax_table(taxa))
ps_contaminants1 <- ps_contaminants %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order) & Kingdom!="Eukaryota" | is.na(Kingdom ))
ps_contaminants2<- subset_taxa(ps_contaminants1, Phylum != "unassigned")

Contaminants<-prune_taxa(taxa_sums(ps_contaminants2) > 0.001 , ps_contaminants2)
plot_bar(tax_glom(Contaminants, taxrank = "Genus"), fill="Genus")
#okay, Streptococcus takes up most, makes sense, followed by Vigrio, Proteiniphilum. 
#only Enterococus is from the mock, no clear others from the positive controls. 
#Double-check if Enterococcus is removed from the mock when doing the analysis of its composition, then decide how to proceed

otu_table(Contaminants)
#create new object with taxonomy table of genera
Contaminants_Genus_level<-tax_glom(Contaminants, taxrank="Genus")
otu_table(Contaminants_Genus_level)
Contaminants_Genus_level_data<-otu_table(Contaminants_Genus_level)
head(Contaminants_Genus_level_data)
#make data frame where you take proportional data for each genus
Contam_ASV<-data.frame(rowSums(Contaminants_Genus_level_data))
colnames(Contam_ASV)<-"Contam_abundance"
#create new object with taxonomy table from proportional genera and % abundance just created
Contam_ab<-cbind(as.data.frame(Contaminants_Genus_level@tax_table@.Data), Contam_ASV)

##############################################
#Make Phyloseq Objects

#clear workspace
rm(list=ls())
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(DECIPHER)
library(microViz)
library(xlsx)
library(vegan)

#save phyloseq component data
#save(seqtab.nochim, taxa,decont_samp, file ="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_object_components.RData")
load("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_object_components.RData")

full_otu_tb<-otu_table(seqtab.nochim, taxa_are_rows = F)
decont_otu_tb<-otu_table(decont_samp, taxa_are_rows = F)

tax_tb<-tax_table(taxa)
#if you have phylogenetic data you can add this
#phy_t<-phy_tree(tree$tree)

#load metadata
#decontaminated data without negative controls
decont_metadata<- read.table(file="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/decont_seq_metadata.txt",header=T,sep="\t")
#full dataset without running SCRuB decontamination
full_metadata<-read.table(file="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/full_seq_metadata.txt",header=T,sep="\t")
#make the row names of the metadata correspond to the sample names
rownames(full_metadata)<-full_metadata$Sample_ID
rownames(decont_metadata)<-decont_metadata$Sample_ID
#colnames(decont_metadata)<- c('Sample.name','Treatment', 'Replicate', 'Nest', 'Plate', 'Metadata_missing')


#check out variables in metadata
head(decont_metadata)
head(full_metadata)

#check dimensions of everything
dim(decont_metadata)
# 75, 13
dim(full_metadata)
# 80, 13
dim(decont_otu_tb)
# 75, 5504
dim(full_otu_tb)
# 80, 5504
#great the dimensions of everything look like they should match

#convert the metadata dataframes to phyloseq sample data 
full_sample_data<-sample_data(full_metadata)
decont_sample_data<-sample_data(decont_metadata)

#make phyloseq object with full sequencing data
ps_full_v1<-phyloseq(sample_data = full_sample_data,
                     otu_table = full_otu_tb,
                     tax_table = tax_tb)
#make phyloseq object with decontaminated sequencing data
ps_decontam_v1<-phyloseq(sample_data = decont_sample_data,
                         otu_table = decont_otu_tb,
                         tax_table = tax_tb)


##############################################
#Remove Chloroplasts, Mitochondria, Eukarya from both phyloseq objects

ps_full_v2 <- ps_full_v1 %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order) & Kingdom!="Eukaryota" | is.na(Kingdom ))
ps_full_v3<- subset_taxa(ps_full_v2, Phylum != "unassigned")

#see how many ASVs classified to Genus before removing and after
tax_NA<-sum(is.na(tax_table(ps_full_v1)[,"Genus"]))
tax_all<-length(tax_table(ps_full_v1)[,"Genus"])
100-tax_NA/tax_all*100

tax_NA1<-sum(is.na(tax_table(ps_full_v2)[,"Genus"]))
tax_all1<-length(tax_table(ps_full_v2)[,"Genus"])
100-tax_NA1/tax_all1*100

tax_NA2<-sum(is.na(tax_table(ps_full_v3)[,"Genus"]))
tax_all2<-length(tax_table(ps_full_v3)[,"Genus"])
100-tax_NA2/tax_all2*100

# genus classification: 62.68% ASVS before
# genus classification: 63.87% ASVs after
# therefore not that many mitochondria and chloroplasts sequenced

#see how many ASVs classified to Family
tax_NA<-sum(is.na(tax_table(ps_full_v3)[,"Family"]))
tax_all<-length(tax_table(ps_full_v3)[,"Family"])
100-tax_NA/tax_all*100
# family classification: 82.23% ASVS
# therefore classification gets much better at the family level

#do the same for the decontaminated object
ps_decontam_v2 <- ps_decontam_v1 %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order) & Kingdom!="Eukaryota" | is.na(Kingdom ))
ps_decontam_v3<- subset_taxa(ps_decontam_v2, Phylum != "unassigned")

tax_NA<-sum(is.na(tax_table(ps_decontam_v3)[,"Genus"]))
tax_all<-length(tax_table(ps_decontam_v3)[,"Genus"])
100-tax_NA/tax_all*100
# genus classification: 63.87% ASVs after

tax_NA<-sum(is.na(tax_table(ps_decontam_v3)[,"Family"]))
tax_all<-length(tax_table(ps_decontam_v3)[,"Family"])
100-tax_NA/tax_all*100
# family classification: 82.23% ASVS
#Great, everything is the same for the decontam samples!

###############################################
# Evaluate Mock sample

#the Zymobiomics microbial community was used to determine extraction and sequencing bias
#the community contained 8 bacteria at the following expected relative abundances based on cells and 16S copy number:

#Pseudomonas aeruginosa: 4.2 
#Escherichia coli: 10.1
#Salmonella enterica: 10.4
#Lactobacillus fermentum: 18.4
#Enterococcus faecalis: 9.9
#Staphylococcus aureus: 15.5
#Listeria monocytogenes: 14.1
#Bacillus subtilis: 17.4

Mock<-prune_samples(sample_names(ps_full_v3) == "A78", ps_full_v3)
Mock<-prune_taxa(taxa_sums(Mock) > 0 , Mock)

plot_bar(tax_glom(Mock, taxrank = "Genus"), fill="Genus")
#cool, all genera are there except names slightly different

otu_table(Mock)
#we find 17 taxa, 5 with less than 8 reads, will be removed in 'decontam'
#not sure about the other 5... lets move on for now

#make the abundance proportional, out of 1.00
Mock.prop<-transform_sample_counts(Mock, function(x) x/sum(x))
Mock.prop.Gen<-tax_glom(Mock.prop, taxrank="Genus")

#plot to check, looks similar to that but proportional abundance
plot_bar(Mock.prop.Gen, fill="Genus")
otu_table(Mock.prop.Gen)
#note that when we tax_glom at the Genus level, we end up with 8 genera, which should be present in the mock

#take the OTU table to make a dataframe
Mock_data<-otu_table(Mock.prop.Gen)
head(Mock_data)

#make data frame where you take proportional data for each genus
Mock_ASV<-data.frame(colSums(Mock_data)*100/length(sample_names(Mock.prop.Gen)))
colnames(Mock_ASV)<-"Detected_abundance"
#now have % abundance

#create new object with taxonomy talbe from proportional genera and % abundance just created
Mock_ab<-cbind(as.data.frame(Mock.prop.Gen@tax_table@.Data), Mock_ASV)

#add the expected abundance of them, make sure the order of the values correspond to the dataframe
Mock_ab$Expected_abundance<-c(9.9, 17.4, 10.1, 10.4, 4.2, 14.1, 18.4, 15.5)
#double-check that they both equal 100
sum(Mock_ab$Detected_abundance)
sum(Mock_ab$Expected_abundance)
#good

#relative abundance deviation in average is <15% according to ZymoBIOMICS
#add column call Percent error
Mock_ab$Percent_error<-(Mock_ab$Detected_abundance-Mock_ab$Expected_abundance)*100/Mock_ab$Expected_abundance
#replace ASV row names with number row names
rownames(Mock_ab)<-c(1:length(rownames(Mock_ab)))
# error potentially introduced by extraction or primers
#save the plot, dataframe, etc for supplementary materials 

#select the key columns of the dataframe
Final_mock_table<- subset(Mock_ab, select= c(Genus, Detected_abundance, Expected_abundance, Percent_error))
#save as CSV
write.csv(Final_mock_table,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/Final_mock_table.csv")
#modified data slightly in excel to then plot
#read in csv again
mock_plot<- read_csv("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/Ant_Yogurt_Metabarcodes/dada2_data/Final_mock_forplot.csv")

ggplot(mock_plot, aes(x=Mock_Abundance, y=Detected_abundance, fill=Genus)) +
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_brewer(palette="Accent")+
  labs(y= "Abundance", x = "Mock Community Positive Control")+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#save the plot as a pdf

#Check how this looks for the decontaminated data
#Mock sample A78 only retains 92% abundance of reads, so good to check 
Mock1<-prune_samples(sample_names(ps_decontam_v3) == "A78", ps_decontam_v3)
Mock1<-prune_taxa(taxa_sums(Mock1) > 0 , Mock1)

plot_bar(tax_glom(Mock1, taxrank = "Genus"), fill="Genus")
#the decontam in SCRuB removed Entercoccus faecilis, which is also a human associated bacterium and could be a true contaminant in other samples. 
#therefore it is a good idea to check out the basic composition and things of both the contaminated and decontaminated datasets!
#especially the overall composition of positive controls

#save current phyloseq objects so ready to work tomorrow
save(ps_full_v3, ps_decontam_v3, file ="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_objects_8.1.24.RData")


###########################################
#Rename ASVs to ASV#

#clear workspace
rm(list=ls())

#load packages
rm(list=ls())
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(DECIPHER)
library(microViz)
library(xlsx)
library(vegan)

load("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_objects_8.1.24.RData")

#first remove any 0s ASVs from the decontaminated dataset
ps_decontam_v4 <- prune_taxa(taxa_sums(ps_decontam_v3) > 0, ps_decontam_v3)
#check the number of ASVs now
length(taxa_names(ps_decontam_v4))
ps_decontam_v4
#cool, 75 samples, 3,829 ASVs

#name ASVs instead of full DNA seqs
dna_decontam <- Biostrings::DNAStringSet(taxa_names(ps_decontam_v4))
names(dna_decontam) <- taxa_names(ps_decontam_v4)
ps_decontam_v4 <- merge_phyloseq(ps_decontam_v4, dna_decontam)
taxa_names(ps_decontam_v4) <- paste0("ASV", seq(ntaxa(ps_decontam_v4)))
ps_decontam_v4
#save the refseq object
refseq(ps_decontam_v4) %>%Biostrings::writeXStringSet("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/ps_decontam_v4_seqs.fa", 
                                                      append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
#save the phyloseq object
save(ps_decontam_v4, file ="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/decontam_phyloseq_object_12.1.24.RData")

#########################################################
#Considering pruning reads

load ("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/decontam_phyloseq_object_12.1.24.RData")

#pruning is sometimes performed to reduce really low-abundance ASVS
#this reduces the complexity of the data and makes it easier to build trees
#this is common for building trees, or to reduce potential contaminants
#see how it looks pruning ASVs with less than 20 or 50 total reads 


#see how these change when removing ASVs with less than 20 or 50 reads across the dataset

ps_decontam_v4_50plus <- prune_taxa(taxa_sums(ps_decontam_v4) > 50, ps_decontam_v4)
ps_decontam_v4_20plus <- prune_taxa(taxa_sums(ps_decontam_v4) > 20, ps_decontam_v4)
ps_decontam_v4_15plus <- prune_taxa(taxa_sums(ps_decontam_v4) > 15, ps_decontam_v4)

# 3829 ASVs before pruning
length(taxa_names(ps_decontam_v4_50plus))
# 441 ASVs
length(taxa_names(ps_decontam_v4_20plus))
# 804 ASVs
length(taxa_names(ps_decontam_v4_15plus))

#have a look at the composition of samples to see how it impacts things
ps_decontam_v4.prop <- transform_sample_counts(ps_decontam_v4, function(otu) otu/sum(otu))
decontam_data_plot<-prune_taxa((taxa_sums(ps_decontam_v4.prop) > 0.05), ps_decontam_v4.prop)
plot_bar(tax_glom(decontam_data_plot, taxrank = "Genus"), fill="Genus", title="Decontam Genus") 
plot_bar(tax_glom(decontam_data_plot, taxrank = "Family"), fill="Family", title="Decontam Family")

ps_decontam_v4.prop.plus20 <- transform_sample_counts(ps_decontam_v4_20plus, function(otu) otu/sum(otu))
decontam_data_plot_plus20<-prune_taxa((taxa_sums(ps_decontam_v4.prop.plus20) > 0.05), ps_decontam_v4.prop.plus20)
plot_bar(tax_glom(decontam_data_plot_plus20, taxrank = "Genus"), fill="Genus", title="Decontam Genus >20") 
plot_bar(tax_glom(decontam_data_plot_plus20, taxrank = "Family"), fill="Family", title="Decontam Family >20")

#things don't look very different at all pruning 20 ASVs. Lets continue with that dataset
ps_decontam_v5 <- prune_taxa(taxa_sums(ps_decontam_v4) > 20, ps_decontam_v4)
View(sample_data(ps_decontam_v5))

#lets subset the dataset without the controls
#remove conventional controls
ps_v5_1 <- prune_samples(!sample_data(ps_decontam_v5)$Ant_or_Yog_or_Ctl == "yogurt_control", ps_decontam_v5)
#remove supernatants
View(sample_data(ps_v5_1))
ps_v5_2 <- prune_samples(!sample_data(ps_v5_1)$Ant_or_Yog_or_Ctl == "supernatant", ps_v5_1)
View(sample_data(ps_v5_2))
#remove mock
ps_v5_3 <- prune_samples(!sample_data(ps_v5_2)$Ant_or_Yog_or_Ctl == "mock", ps_v5_2)
#remove no inoculum controls
ps_v5_big <- prune_samples(!sample_data(ps_v5_3)$Treatment == "ctrl", ps_v5_3)
View(sample_data(ps_v5_big))

#save the phyloseq object
save(ps_v5_big, file ="/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_experimental_samples.RData")

##########
#done
#########