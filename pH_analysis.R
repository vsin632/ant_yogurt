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
####Analysis of pH data
######################################
#clear workspace
rm(list=ls())

#load packages
library("ggplot2")
library("tidyverse")

#load data
data=read.table("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/pH/pH_metadata.txt",header=T,sep="\t")

#subset samples with only pH data
pH.data <- data[!is.na(data$pH_Quantification),]

pH.summary <- pH.data %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(pH_Quantification, na.rm = TRUE),
    pH_Quantification = mean(pH_Quantification))
pH.summary

#relevel the data
pH.relevel<- pH.data %>%
  mutate(Sample_type = fct_relevel(Sample_type, 
                                   "li_ant_yog_s", "li_ant_yog_a", "de_ant_yog_s", 
                                   "de_ant_yog_a", "Formic acid", "Conventional", 
                                   "No inoculum", "fr_ant_yog_s"))
#make plot with frozen sample
plot1 <- ggplot(pH.relevel, aes(x = Sample_type, y = pH_Quantification)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, fill=c("#FF9966","#FF9966","#990066","#990066","#0099CC", "#CCCCCC","#666666","#009900"))+
  geom_errorbar( aes(x=Sample_type, ymin=pH_Quantification-sd, ymax=pH_Quantification+sd), data=pH.summary, width = 0.3)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.15)+
  scale_y_continuous(limits = c(0,7.5), expand = c(0, 0)) +
  labs(title="Acidity", face="bold", x= "", y="pH")+
  theme_classic()
plot1
#note that VRV made the final version of this figure in GraphPad Prism

##########
#done
#########