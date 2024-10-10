#####################################
# author: Veronica M. Sinotte
# title: Ant yogurt 16S data analysis
# date: 2024-Oct-9
# output: R file
# R version: "R version 4.3.1 (2023-06-16)"
# manuscript: Making yogurt with the ant holobiont reveals bacteria, acids, and enzymes for food fermentation
# note: co-author VRV also contributed to the analysis, and some final figures were made by her in GraphPad Prisim, these cases are noted
######################################

##################################
# Micrombiome composition analysis
##################################

###################################
#Make figures to show seasonality of the Formica rufa and Formica polyctena microbiome composition

#using the microviz package, see more on visualizations below: 
#https://david-barnett.github.io/microViz/articles/shao19-analyses.html 

rm(list=ls())

library(microViz)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(microbiome)
library(knitr)

#load
load ("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/phyloseq_experimental_samples.RData")

#have a look at the data frame for the phyloseq object
View(sample_data(ps_v5_big))

#check the phyloseq object to see if you can use microviz with it
phyloseq_validate(ps_v5_big)
#it says that there are NAs detected, and it can be fixed with tax_fix()
#take a quick look at the tax table
tax_table(ps_v5_big)[40:54, 4:7]

#fix the errors with tax_fix
ps_big_fixed <- ps_v5_big %>%
  tax_fix(
    min_length = 4,
    unknowns = c("unidentified", "NA"),
    sep = " ", anon_unique = FALSE,
    suffix_rank = "classified"
  )

#check that the problems are resolved
tax_table(ps_big_fixed)[40:54, 4:7]
phyloseq_validate(ps_big_fixed)

#NOTE: the following visualizations are done with microviz package
#looks good, lets try to plot some things
#subset ant samples
Ants <- subset_samples(ps_big_fixed, Ant_or_Yog_or_Ctl == "ant")

#plots at the genus level, for the top 8 genera
plot_ants_genus_formatted<- Ants %>%
  ps_filter(Treatment == "live") %>%
  comp_barplot(tax_level= "Genus", facet_by = "Sample_type", n_taxa=8)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_ants_genus_formatted

#plots at the genus level, for the top 10 genera
plot_ants_genus_formatted<- Ants %>%
  ps_filter(Treatment == "live") %>%
  comp_barplot(tax_level= "Genus", facet_by = "Sample_type", n_taxa=10)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_ants_genus_formatted

#plots at the family level, for the top 8 families
plot_ants_family_formatted<- Ants %>%
  ps_filter(Treatment == "live") %>%
  comp_barplot(tax_level= "Family", facet_by = "Sample_type", n_taxa=8)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_ants_family_formatted

#plots at the family level, for the top 10 families
plot_ants_family_formatted<- Ants %>%
  ps_filter(Treatment == "live") %>%
  comp_barplot(tax_level= "Family", facet_by = "Sample_type", n_taxa=10)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_ants_family_formatted

####################################################
#Visualise the microbiome composition of the ant and the respective yogurts

#remove Formica rufa samples because they were not used for ant yogurt
ps_v5_no_F <- prune_samples(!sample_data(ps_big_fixed)$Sample_type == "li_ant_Fr_s", ps_big_fixed)
Ant_Yogurts <- prune_samples(!sample_data(ps_v5_no_F)$Sample_type == "li_ant_Fr_a", ps_v5_no_F)

#have a look at the data frame
View(sample_data(Ant_Yogurts))

# We will made the plots for the top 10 taxa at the genus and family level

#Genus-level plot
plot_antyogurt_genus_formatted<- Ant_Yogurts %>%
  comp_barplot(tax_level= "Genus", facet_by = "Sample_type", n_taxa=10)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.5, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_antyogurt_genus_formatted

#Family-level plot
plot_antyogurt_genus_formatted<- Ant_Yogurts %>%
  comp_barplot(tax_level= "Family", facet_by = "Sample_type", n_taxa=10)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.5, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_antyogurt_family_formatted

########################################################
# Alpha Diversity of the yogurts, figures and statistics

#first make a table of the alpha diversity indicies
#start with basic phyloseq functions
ant_yogurt_alpha <-estimate_richness(ps_big_fixed, measures = c("Chao1", "Shannon", "Simpson" ))
View(ant_yogurt_alpha)

#Have an initial look at all the samples
plot_richness(Ant_Yogurts, x="Sample_type", measures=c("Shannon", "Simpson"), color="Treatment")+
  geom_point(size=3)+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1),  geom="crossbar", alpha=0.7)

#adjust some colors/formats
plot_richness(Ant_Yogurts, x="Sample_type", measures=c("Shannon"), color="Treatment")+
  geom_point(size=3)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),  geom="crossbar", alpha=0.7)+
  labs(x= "", y="ng bacterial DNA")+
  theme(axis.text.x = element_text (size=10, colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=10, colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill ="white"),
        panel.grid.minor = element_blank())

#Okay, now lets focus on yogurts and subset those samples
Yogurts <- subset_samples(ps_big_fixed, Ant_or_Yog_or_Ctl == "yogurt")

#first recalculate alpha diversity indices and make a table
#start with basic phyloseq functions
yogurt_alpha <-estimate_richness(Yogurts, measures = c("Chao1", "Shannon", "Simpson" ))
View(yogurt_alpha)

#put alpha diversity in the metadata so we can run some stats later
alpha_metadata<-sample_data(Yogurts)<-cbind(sample_data(Yogurts), yogurt_alpha)

#now give a basic yogurt plot a go
plot_richness(Yogurts, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")+
  geom_point(size=3)+
  stat_summary(fun.data="mean_se", fun.args = list(mult=1),  geom="crossbar", alpha=0.7)
#great

#format it for the manuscript
plot_richness(Yogurts, x="Treatment", measures=c("Shannon"), color="Treatment")+
  geom_point(size=2)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),  geom="crossbar", alpha=0.7)+
  labs (y= "Alpha diveristy (Shannon)")+
  scale_colour_manual(values=c("#CC6699", "#009900","#FF9900"))+
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=10, colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill ="gray93"),
        panel.grid.minor = element_blank())
#gorgeous

#lets run some quick statistics on that

#libaray the stats packages you need
library(nlme)
library(emmeans)
library(lme4)
library(lmerTest)

#source the code to test the assumptions of the models we will run
source("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/diagnostic_fcns.r", chdir = FALSE)

lm_alpha_yogurts <- lm(data=alpha_metadata, Shannon~Treatment)
diagnostics.plot(lm_alpha_yogurts)
#looks pretty good
summary(lm_alpha_yogurts)
TukeyHSD(aov(lm_alpha_yogurts))
#great

###########################################################
#Beta Diversity figures and statistics for ants and yogurts 

#make the data proportional (microviz did that before for us, but we need to do it manually here)
ps.prop.prune <- transform_sample_counts(Ant_Yogurts, function(otu) otu/sum(otu))
#ordinate data with Bray Curtis method
ord.nmds.bray.2 <- ordinate(ps.prop.prune, method="NMDS", distance="bray")

#Make NMDS plot of the Bray-Curtis ordination for the ants and yogurts
plot_ordination(ps.prop.prune, ord.nmds.bray.2 , color="Treatment")+
  scale_colour_manual(values=c("#660033", "#CC6699", "#009900","#CC3300", "#FF9900"))+
  stat_ellipse(type = "norm", lwd=0.5)+
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=10, colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill ="gray93"),
        panel.grid.minor = element_blank())

#Now lets do the same, but just for the yogurts
#subset the yogurt samples
Yogurts <- subset_samples(ps_big_fixed, Ant_or_Yog_or_Ctl == "yogurt")
#make the data proportional
Yogurts.prop <- transform_sample_counts(Yogurts, function(otu) otu/sum(otu))
#ordinate with the Bray-Curtis method
ord.nmds.bray.yogurts <- ordinate(Yogurts.prop, method="NMDS", distance="bray")

#Make the plot for the manuscript
plot_ordination(Yogurts.prop, ord.nmds.bray.yogurts , color="Treatment")+
  scale_colour_manual(values=c("#CC6699", "#009900", "#FF9933"))+
  stat_ellipse(type = "norm", lwd=0.5)+
  theme(axis.text.x = element_text (size=10, colour = "black"),
        panel.border = element_rect(colour = "black", linewidth  = 1.0, linetype = "solid", fill = NA),
        axis.title.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=10, colour = "black"),
        legend.text = element_text (size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill ="gray93"),
        panel.grid.minor = element_blank())

#lets do the stats on it with an Adonis PERMANOVA test

#makure sure we have the right package
library(vegan)

#run the Adonis PERMANOVA
set.seed(123)
options(scipen=999)
adonis_yogurt<-adonis2(vegdist(otu_table(Yogurts), method = "bray") ~ 
                         sample_data(Yogurts)$Treatment,
                       permutations = 1e4)
#See output of model
adonis_yogurt
#treatment has a significant effect overall

#do a pairwise comparison with 'pairAdonis' package
# Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4
# https://github.com/pmartinezarbizu/pairwiseAdonis
library(pairwiseAdonis)

#run the parwise comparisons
pairwise_adonis_yogurt<-pairwise.adonis(vegdist(otu_table(Yogurts), method = "bray"), sample_data(Yogurts)$Treatment)
pairwise_adonis_yogurt
#even after p-adjustment, all comparisions are significant

################################################################
# Quick Look at Abundances of Select Bacterial Taxa per Sample

#we will generate a table for this data and then integrate it later when we look at the bacterial load from the qPCR experiments
#we do this for lactic acid bacteria, acetic acid bacteria, and bacillus
#additionally we do this for Streptococcus, because we find proteases from it in the proteomics, and want to determine the abundance in some samples

#using the big data set
ps_big_fixed

#agglomerate by families or genera
ps_big_families <- tax_glom (ps_big_fixed, taxrank = "Family")
ps_big_genera <- tax_glom (ps_big_fixed, taxrank = "Genus")

#calculate relative abundances
ps_families_prop <-transform_sample_counts(ps_big_families, function(otu) otu/sum(otu))
ps_genera_prop <-transform_sample_counts(ps_big_genera, function(otu) otu/sum(otu))

#subset the families "Bacillaceae", "Acetobacteraceae", "Lactobacillaceae"
ps_families_subset = subset_taxa(ps_families_prop, 
                                 Family=="Bacillaceae" |Family =="Acetobacteraceae"|Family== "Lactobacillaceae"|Family== "Streptococcaceae")
#View tax table to make sure it works
View(tax_table(ps_families_subset))

#subset the genera "Bacillus", "Fructilactobacillus", "Neokomagataea", "Streptococcus"
ps_genera_subset = subset_taxa(ps_genera_prop, 
                               Genus=="Bacillus" |Genus =="Neokomagataea"|Genus== "Fructilactobacillus"|Genus== "Streptococcus")
#View tax table to make sure it works
View(tax_table(ps_genera_subset))

#Make dataframe with abundance of the selected families and genera for each sample
Families_Summary <- psmelt(ps_families_subset)
View(Families_Summary)
Genera_Summary <-psmelt(ps_genera_subset)
View(Genera_Summary)

write.csv(Families_Summary,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/microbial load/Family_Abundances2.csv")
write.csv(Genera_Summary,"/Users/csp839/Documents/R scripts ant yogurt/R_scripts/microbial load/Genera_Abundances2.csv")

#these will be manually integrated (in excel) into the data for the next step, looking at bacterial loads from qpcr data

##########
#done
#########