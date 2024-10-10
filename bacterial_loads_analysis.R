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
#### Bacterial Loads Analysis
######################################

#this examines the data from qPCR experiments and microbiome composition

#clear workspace
rm(list=ls())

#read in data
qpcr.data=read.table("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/microbial load/load_metadata.txt",header=T,sep="\t")

#load packages
library("ggplot2")
library("tidyverse")

#check dataframe out
str(qpcr.data)
#number of variables and observations, 14 variables, 83 observations (correct)
head(qpcr.data)
#looks good

#create data frame for yogurt samples only
yogurt.data<-qpcr.data[(qpcr.data$Ant_or_Yog_or_Ctl=="yogurt"),]

#remove Cultura sample (A86) because the first wash-and-pellet did not work
yogurt.data.1<-yogurt.data[!yogurt.data$Sample_ID=="A86",]

#create plot of ng DNA/mL in yogurt samples, first look
ggplot(yogurt.data.1, aes(x = Sample_type, y = Ng_Bact_DNA_mL)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

ggplot(yogurt.data.1, aes(x = Treatment2, y = Ng_Bact_DNA_mL)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#remember that the conventional 
#conventional samples much higher than the rest, try to log-transform adding new column
yogurt.data.1$log.NgDNA=log10(yogurt.data.1$Ng_Bact_DNA_mL+1) 

#have another look at he basic plots
ggplot(yogurt.data.1, aes(x = Sample_type, y = log.NgDNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

ggplot(yogurt.data.1, aes(x = Treatment2, y = log.NgDNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#ok, not super great, but need to report it, lets make the plot
#reorder to make the fit VRVs others graphs
yogurt.relevel<- yogurt.data.1 %>%
  mutate(Sample_type = fct_relevel(Sample_type, 
                                   "li_ant_yog_s", "li_ant_yog_a", "de_ant_yog_s", 
                                   "de_ant_yog_a", "fr_ant_yog_s", "conventional", 
                                   "no_inoc"))
yogurt.relevel.2<- yogurt.data.1 %>%
  mutate(Treatment2 = fct_relevel(Treatment2, 
                                  "Alive", "Dehydrated", "Frozen", "Conventional", "No inoculum"))

#another look at plot
ggplot(yogurt.relevel, aes(x = Sample_type, y = log.NgDNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

ggplot(yogurt.relevel.2, aes(x = Treatment2, y = log.NgDNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#great, lets do a summary for error bars, then make the same theme plot as for pH
load.summary <- yogurt.relevel %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(log.NgDNA, na.rm = TRUE),
    log.NgDNA = mean(log.NgDNA))
load.summary

load.summary.2 <- yogurt.relevel.2 %>%
  group_by(Treatment2) %>%
  summarise(
    sd = sd(log.NgDNA, na.rm = TRUE),
    log.NgDNA = mean(log.NgDNA))
load.summary.2

#plot formated to theme

#plots with the yogurts seperated by seasons
plot1 <- ggplot(yogurt.relevel, aes(x = Sample_type, y = log.NgDNA)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, fill=c("#FF9966","#FF9966","#990066","#990066","#009900", "#CCCCCC","#666666"))+
  geom_errorbar( aes(x=Sample_type, ymin=log.NgDNA-sd, ymax=log.NgDNA+sd), data=load.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black",  binwidth = 0.03)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng bacterial DNA/mL (log10)")+
  scale_x_discrete(labels=c("li_ant_yog_s" = "Spring", "li_ant_yog_a" = "Autumn",
                            "de_ant_yog_s" = "Spring", "de_ant_yog_a" = "Autumn","fr_ant_yog_s" = "Spring",
                            "conventional" = "Conventional","no_inoc" = "No inoculum"))+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot1

#plots with the no seasonal seperation, only yogurt treatment groups
plot1.1 <- ggplot(yogurt.relevel.2, aes(x = Treatment2, y = log.NgDNA)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, fill=c("#FF9966","#990066","#009900", "#CCCCCC","#666666"))+
  geom_errorbar( aes(x=Treatment2, ymin=log.NgDNA-sd, ymax=log.NgDNA+sd), data=load.summary.2, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black",  binwidth = 0.03)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "Yogurts", y="ng bacterial DNA/mL (log10)")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot1.1

#make a plot without the conventional
#first remove conventional and no inoculum control from the dataset
yogurt.relevel.2<-yogurt.relevel[!yogurt.relevel$Sample_type=="conventional",]
ant.yogs<-yogurt.relevel.2[!yogurt.relevel.2$Sample_type=="no_inoc",]

#basic plot
ggplot(ant.yogs, aes(x = Sample_type, y = log.NgDNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#great, quickly calc summary stats
yog.load.summary <- ant.yogs %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(log.NgDNA, na.rm = TRUE),
    log.NgDNA = mean(log.NgDNA))
yog.load.summary

#cool, go for the formatted plot
plot2 <- ggplot(ant.yogs, aes(x = Sample_type, y = log.NgDNA)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, fill=c("#FF9966","#FF9966","#990066","#990066","#009900"))+
  geom_errorbar( aes(x=Sample_type, ymin=log.NgDNA-sd, ymax=log.NgDNA+sd), data=yog.load.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.005)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng bacterial DNA/mL (log10)")+
  scale_x_discrete(labels=c("li_ant_yog_s" = "Spring", "li_ant_yog_a" = "Autumn",
                            "de_ant_yog_s" = "Spring", "de_ant_yog_a" = "Autumn","fr_ant_yog_s" = "Spring"))+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot2

#we know from the metabarcoding data that these samples are dominated by Bacillaceae, Lactobacillaceae, Acetobacteraceae, 
#also Anaplasmatacae but that can't actually grow in the yogurt so we won't check it
#make a scatter plot with the x Bacillaceae abundance, y Ng Bact DNA
ant.yogs.Bacill <- ant.yogs[!is.na(ant.yogs$Bacillaceae),]
plot3<-ggplot(ant.yogs.Bacill, aes(x = Bacillaceae, y = log.NgDNA, fill=Treatment)) +
  geom_point(aes(color=Treatment))+
  geom_smooth(method=lm, aes(colour=Treatment))+
  scale_fill_manual(values=c("#990066", "#009900","#FF9966"))+
  scale_colour_manual(values=c("#990066", "#009900","#FF9966"))+
  labs(x= "Bacilliaceae relative abundance", y="ng bacterial DNA/mL (log10)")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x  = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot3

#libaray the stats packages you need
library(nlme)
library(emmeans)
library(lme4)
library(lmerTest)

#source the code to visualise the assumptions of the model (histogram of resid, qq-plot of resid, resid against fitted values)
source("/Users/csp839/Documents/R scripts ant yogurt/R_scripts/diagnostic_fcns.r", chdir = FALSE)

#first see what the overall model would look like
lm_bacill_yogs<- lm(data=ant.yogs.Bacill, log.NgDNA~Bacillaceae)
diagnostics.plot(lm_bacill_yogs)

#looks possible, lets see if it works when we subset

#make groups for each of the Treatment groups
dehyd.data<-ant.yogs.Bacill[(ant.yogs.Bacill$Treatment=="dehyd_ant"),]
froz.data<-ant.yogs.Bacill[(ant.yogs.Bacill$Treatment=="froz_ant"),]
live.data<-ant.yogs.Bacill[(ant.yogs.Bacill$Treatment=="live_ant"),]

#linear models for Bacillaceae
lm_bacill_dehyd<- lm(data=dehyd.data, log.NgDNA~Bacillaceae)
diagnostics.plot(lm_bacill_dehyd)
summary(lm_bacill_dehyd)
#looks okay

lm_bacill_froz<- lm(data=froz.data, log.NgDNA~Bacillaceae)
diagnostics.plot(lm_bacill_froz)
summary(lm_bacill_froz)
#looks okay

lm_bacill_live<- lm(data=live.data, log.NgDNA~Bacillaceae)
diagnostics.plot(lm_bacill_live)
summary(lm_bacill_live)
#looks okay

#linear models for Lactobacillaceae
lm_lacto_dehyd<- lm(data=dehyd.data, log.NgDNA~Lactobacillaceae)
diagnostics.plot(lm_lacto_dehyd)
summary(lm_lacto_dehyd)
#okay, hisotgram of resid not ideal

lm_lacto_froz<- lm(data=froz.data, log.NgDNA~Lactobacillaceae)
diagnostics.plot(lm_lacto_froz)
summary(lm_lacto_froz)
#looks okay

lm_lacto_live<- lm(data=live.data, log.NgDNA~Lactobacillaceae)
diagnostics.plot(lm_lacto_live)
summary(lm_lacto_live)
#looks okay

#cool the diagnostic plots look good for all of the linear models. Add these to the manuscript.

#try out for Lactobacillaceae
plot4<-ggplot(ant.yogs.Bacill, aes(x = Lactobacillaceae, y = log.NgDNA, fill=Treatment)) +
  geom_point(aes(color=Treatment))+
  geom_smooth(method=lm, aes(colour=Treatment))+
  scale_fill_manual(values=c("#990066", "#009900","#FF9966"))+
  scale_colour_manual(values=c("#990066", "#009900","#FF9966"))+
  labs(x= "Lactobacillaceae relative abundance", y="ng bacterial DNA/mL (log10)")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x  = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot4
#aligns with earlier compositional analysis that Lactobacillaceae is only really in live ant yogurts, sometimes in very low abundance in dehydrated

#Run for Acetobacteraceae
plot4.1<-ggplot(ant.yogs.Bacill, aes(x = Acetobacteraceae, y = log.NgDNA, fill=Treatment)) +
  geom_point(aes(color=Treatment))+
  geom_smooth(method=lm, aes(colour=Treatment))+
  scale_fill_manual(values=c("#990066", "#009900","#FF9966"))+
  scale_colour_manual(values=c("#990066", "#009900","#FF9966"))+
  labs(x= "Acetobacteraceae relative abundance", y="ng bacterial DNA/mL (log10)")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x  = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot4.1
#No clear patter here for Acetobacteraceae

#Just have a look at the alive ant yogurts, for Lactobacillaceae and Acetobacteraceae
#These dominate the community aside from Anaplasmataceae, which is obligately intracellular
alive.yogs <- ant.yogs[(ant.yogs$Treatment=="live_ant"),]

#make plots for Lactobacillaceae
plot5<-ggplot(alive.yogs, aes(x = Lactobacillaceae, y = Ng_Bact_DNA_mL, fill="#FF9966")) +
  geom_point(aes(shape=Season), colour="#FF9966", size=3)+
  geom_smooth(method=lm, colour="#FF9966")+
  labs(x= "Lactobacillaceae relative abundance", y="ng bacterial DNA/mL", title= "Alive ant yogurts")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        axis.title.x  = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot5
#looks great

#lets do it for Acetobacteraceae
plot6<-ggplot(alive.yogs, aes(x = Acetobacteraceae, y = Ng_Bact_DNA_mL, fill="#993300")) +
  geom_point(aes(shape=Season), colour="#993300", size=3)+
  geom_smooth(method=lm, colour="#993300")+
  labs(x= "Acetobacteraceae relative abundance", y="ng bacterial DNA/mL", title= "Alive ant yogurts")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        axis.title.x  = element_text(size=10, face="bold", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot6

###################################
#Does bacteria grow in the yogurts!?
#now lets think about how much bacteria is in the yogs from the ants to begin with
#In the dataframe, I have created a 'Hypothetical load of starter or yog'. 
#This is calculated by correcting for the total volume of the yogurt (30 mL) or number of ants (52), based on weight

#remove any samples that do not have data for that variable
starters.yogs <- qpcr.data[!is.na(qpcr.data$Hypothet_Ng_Bac_Starter_or_Yog),]
#and lets remove Formica rufa, because we only made yogs with Formica polyctena
#autumn F. rufa
starters.yogs.2<- starters.yogs[!starters.yogs$Sample_type=="li_ant_Fr_a",]
#spring F. rufa
starters.yogs.3<- starters.yogs.2[!starters.yogs.2$Sample_type=="li_ant_Fr_s",]
#no inoculum controls
starters.yogs.4<- starters.yogs.3[!starters.yogs.3$Sample_type=="no_inoc",]

#have a look at all bacteria, but also Bacillaceae, Lactobacillaceae, and Acetobacteraceae because they were in isolate data
starters.yogs.4$Hypothetical.Lacto = starters.yogs.4$Lactobacillaceae * starters.yogs.4$Hypothet_Ng_Bac_Starter_or_Yog
starters.yogs.4$Hypoth.Aceto = starters.yogs.4$Acetobacteraceae * starters.yogs.4$Hypothet_Ng_Bac_Starter_or_Yog
starters.yogs.4$Hypoth.Bacill = starters.yogs.4$Bacillaceae * starters.yogs.4$Hypothet_Ng_Bac_Starter_or_Yog
starters.yogs.4$Hypoth.Lacto.Aceto = starters.yogs.4$Hypothetical.Lacto + starters.yogs.4$Hypoth.Aceto

#need to remove 2 samples we don't have seq data for (double-check if true when redone with new abundances)
starters.yogs.5 <- starters.yogs.4[!is.na(starters.yogs.4$Bacillaceae),]
View (starters.yogs.5)

#do comparision across all live, dehy, froz, and seasons
#Relevel to make the plot formatted correct

big.comparison<- starters.yogs.5 %>%
  mutate(Sample_type = fct_relevel(Sample_type, 
                                   "li_ant_Fp_s","li_ant_yog_s", "li_ant_Fp_a","li_ant_yog_a", 
                                   "de_ant_Fp_s", "de_ant_yog_s", "de_ant_Fp_a", "de_ant_yog_a",
                                   "fr_ant_Fp_s", "fr_ant_yog_s"))  
#Try out a basic plot
ggplot(big.comparison, aes(x = Sample_type, y = Hypothet_Ng_Bac_Starter_or_Yog)) +
  geom_dotplot(binaxis = "y", stackdir = "center")


#run quick stats
big.comparison.summary <- big.comparison %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(Hypothet_Ng_Bac_Starter_or_Yog, na.rm = TRUE),
    Hypothet_Ng_Bac_Starter_or_Yog = mean(Hypothet_Ng_Bac_Starter_or_Yog))
big.comparison.summary

#####################
#run the formatted plots
#for total hypothetical 

plot7 <- ggplot(big.comparison, aes(x = Sample_type, y = Hypothet_Ng_Bac_Starter_or_Yog)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, size=1, color=c("#FF9966","#FF9966","#FF9966","#FF9966","#990066","#990066","#990066","#990066","#009900","#009900" ), fill=c("#FFFFFF","#FF9966","#FFFFFF","#FF9966","#FFFFFF","#990066","#FFFFFF","#990066", "#FFFFFF","#009900" ))+
  geom_errorbar( aes(x=Sample_type, ymin=Hypothet_Ng_Bac_Starter_or_Yog-sd, ymax=Hypothet_Ng_Bac_Starter_or_Yog+sd), data=big.comparison.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 1.0)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt", title = "Complete")+
  scale_x_discrete(labels=c("li_ant_Fp_s"="Spring ants","li_ant_yog_s"="Spring yogurts", "li_ant_Fp_a"="Autumn ants","li_ant_yog_a"="Autumn yogurts", 
                            "de_ant_Fp_s" = "Spring ants", "de_ant_yog_s"= "Spring yogurts", "de_ant_Fp_a"="Autumn ants", "de_ant_yog_a"="Autumn yogurts",
                            "fr_ant_Fp_s" = "Spring ants", "fr_ant_yog_s" = "Spring yogurts"))+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot7

#with facet wrap, you can crop them later
plot7.1 <- ggplot(big.comparison, aes(x = Sample_type, y = Hypothet_Ng_Bac_Starter_or_Yog)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, size=1, color=c("#FF9966","#FF9966","#FF9966","#FF9966","#990066","#990066","#990066","#990066","#009900","#009900" ), fill=c("#FFFFFF","#FF9966","#FFFFFF","#FF9966","#FFFFFF","#990066","#FFFFFF","#990066", "#FFFFFF","#009900" ))+
  geom_errorbar( aes(x=Sample_type, ymin=Hypothet_Ng_Bac_Starter_or_Yog-sd, ymax=Hypothet_Ng_Bac_Starter_or_Yog+sd), data=big.comparison.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 1.0)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt", title = "Complete")+
  facet_wrap(~Treatment2, scales="free_x")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot7.1

#we know that bacillaceae is only dominant in dehydrated and frozen, so lets plot that in starters and yogurts
#first we remove the alive ant samples
bacillaceae.comp<- big.comparison[!big.comparison$Treatment2=="Alive",]

#quick look
ggplot(bacillaceae.comp, aes(x = Sample_type, y = Hypoth.Bacill)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#summary stats for plot
bacillaceae.comp.summary <- bacillaceae.comp %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(Hypoth.Bacill, na.rm = TRUE),
    Hypoth.Bacill = mean(Hypoth.Bacill))
bacillaceae.comp.summary

#make a formatted plot for the Lactobacillaceae

plot8 <- ggplot(bacillaceae.comp, aes(x = Sample_type, y = Hypoth.Bacill)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, color=c("#990066","#990066","#990066","#990066","#009900","#009900"), fill=c("#FFFFFF","#990066","#FFFFFF","#990066", "#FFFFFF","#009900"), size=1)+
  geom_errorbar( aes(x=Sample_type, ymin=Hypoth.Bacill-sd, ymax=Hypoth.Bacill+sd), data=bacillaceae.comp.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = .80)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt")+
  facet_wrap(~Treatment2, scales="free_x")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot8

#lets run some stats on this 
bacillaceae.comp$log.Hypoth.Bacill=log10(bacillaceae.comp$Hypoth.Bacill+1) 

#run linear model for dehydrated alone
bacillaceae.comp.dehydrated<- bacillaceae.comp[bacillaceae.comp$Treatment2=="Dehydrated",]

lm_bacill_antyog_de1<- lm(data=bacillaceae.comp.dehydrated,log.Hypoth.Bacill ~ Treatment*Season)
diagnostics.plot(lm_bacill_antyog_de1)
#diagnostic plot looks much better
summary(lm_bacill_antyog_de1)

lm_bacill_antyog_de2<- lm(data=bacillaceae.comp.dehydrated,log.Hypoth.Bacill ~ Treatment+Season)
diagnostics.plot(lm_bacill_antyog_de2)
#diagnostic plot looks much better
summary(lm_bacill_antyog_de2)

anova(lm_bacill_antyog_de1, lm_bacill_antyog_de2, test="F")
#the interaction is not significant

#drop Treatment
lm_bacill_antyog_de3<- lm(data=bacillaceae.comp.dehydrated,log.Hypoth.Bacill ~ Season)
diagnostics.plot(lm_bacill_antyog_de3)
#diagnostic plot looks much better
summary(lm_bacill_antyog_de3)

#test effect of treatment
anova(lm_bacill_antyog_de2, lm_bacill_antyog_de3, test="F")
#not significant

#drop Season
lm_bacill_antyog_de4<- lm(data=bacillaceae.comp.dehydrated,log.Hypoth.Bacill ~ Treatment)
diagnostics.plot(lm_bacill_antyog_de4)
#diagnostic plot looks much better
summary(lm_bacill_antyog_de4)

#test effect of Season
anova(lm_bacill_antyog_de2, lm_bacill_antyog_de4, test="F")
#not significant

#subset the two seasons to run the pairwise comparison
bacillaceae.comp.dehydrated.spring<- bacillaceae.comp.dehydrated[bacillaceae.comp.dehydrated$Season=="spring",]
bacillaceae.comp.dehydrated.autumn<- bacillaceae.comp.dehydrated[bacillaceae.comp.dehydrated$Season=="autumn",]

lm_bacill_de_spring<- lm(data=bacillaceae.comp.dehydrated.spring, log.Hypoth.Bacill ~ Treatment)
diagnostics.plot(lm_bacill_de_spring)
summary(lm_bacill_de_spring)

#run two sample t-test for unequal variances
t.test(log.Hypoth.Bacill ~ Treatment, data=bacillaceae.comp.dehydrated.spring, var.equal=TRUE)


lm_bacill_de_autumn<- lm(data=bacillaceae.comp.dehydrated.autumn, log.Hypoth.Bacill ~ Treatment)
diagnostics.plot(lm_bacill_de_autumn)
summary(lm_bacill_de_autumn)

#Now look at the frozen
bacillaceae.comp.froz<- bacillaceae.comp[bacillaceae.comp$Treatment2=="Frozen",]

#run some linear models if they work...
#drop Treatment
lm_bacill_antyog_froz1<- lm(data=bacillaceae.comp.froz,log.Hypoth.Bacill ~ Treatment)
diagnostics.plot(lm_bacill_antyog_froz1)
#diagnostic plot looks much better
summary(lm_bacill_antyog_froz1)

#run two sample t-test for unequal variances
t.test(log.Hypoth.Bacill ~ Treatment, data=bacillaceae.comp.froz, var.equal=TRUE)

bacillaceae.comp.summary.1 <- bacillaceae.comp %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(log.Hypoth.Bacill, na.rm = TRUE),
    log.Hypoth.Bacill = mean(log.Hypoth.Bacill))
bacillaceae.comp.summary.1

#make a formatted plot for the Lactobacillaceae

plot8.1 <- ggplot(bacillaceae.comp, aes(x = Sample_type, y = log.Hypoth.Bacill)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, color=c("#990066","#990066","#990066","#990066","#009900","#009900"), fill=c("#FFFFFF","#990066","#FFFFFF","#990066", "#FFFFFF","#009900"), size=1)+
  geom_errorbar( aes(x=Sample_type, ymin=log.Hypoth.Bacill-sd, ymax=log.Hypoth.Bacill+sd), data=bacillaceae.comp.summary.1, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = .05)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt")+
  facet_wrap(~Treatment2, scales="free_x")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot8.1

#now to Lactobacillaceae and Acetobacteraceae
#make new data frame
lactacet.comp<- big.comparison[big.comparison$Treatment2=="Alive",]

#summary stats
lactacet.comp.summary <- lactacet.comp %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(Hypoth.Lacto.Aceto, na.rm = TRUE),
    Hypoth.Lacto.Aceto = mean(Hypoth.Lacto.Aceto))
lactacet.comp.summary


plot9 <- ggplot(lactacet.comp, aes(x = Sample_type, y = Hypoth.Lacto.Aceto)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, color=c("#FF9966","#FF9966","#FF9966","#FF9966"), fill=c("#FFFFFF","#FF9966","#FFFFFF","#FF9966"), size=1)+
  geom_errorbar( aes(x=Sample_type, ymin=Hypoth.Lacto.Aceto-sd, ymax=Hypoth.Lacto.Aceto+sd), data=lactacet.comp.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.85)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot9

#And just Lactobacillaceae
lact.comp.summary <- lactacet.comp %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(Hypothetical.Lacto, na.rm = TRUE),
    Hypothetical.Lacto = mean(Hypothetical.Lacto))
lact.comp.summary


plot10 <- ggplot(lactacet.comp, aes(x = Sample_type, y = Hypothetical.Lacto)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, color=c("#FF9966","#FF9966","#FF9966","#FF9966"), fill=c("#FFFFFF","#FF9966","#FFFFFF","#FF9966"), size=1)+
  geom_errorbar( aes(x=Sample_type, ymin=Hypothetical.Lacto-sd, ymax=Hypothetical.Lacto+sd), data=lact.comp.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.25)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot10
#and if we have a brief look at the ant species 

#lets run some quick stats for this plot

#run a linear model for Season and Ant/Yogurt
lm_lact_antyog<- lm(data=lactacet.comp, Hypothetical.Lacto ~ Treatment*Season)
diagnostics.plot(lm_lact_antyog)
summary(lm_lact_antyog)

#the diagnostic plots could look better, lets try a log transformation
lactacet.comp$log.Hypothetical.Lacto=log10(lactacet.comp$Hypothetical.Lacto+1) 

#try the plot again
lm_lact_antyog2<- lm(data=lactacet.comp,log.Hypothetical.Lacto ~ Treatment*Season)
diagnostics.plot(lm_lact_antyog2)
#diagnostic plot looks much better
summary(lm_lact_antyog2)


lm_lact_antyog3<- lm(data=lactacet.comp, log.Hypothetical.Lacto ~ Treatment+Season)
diagnostics.plot(lm_lact_antyog3)
anova(lm_lact_antyog2, lm_lact_antyog3, test="F")
#the interaction is the better model

#look pairwise within seasons
#subset the data
lactacet.comp.spring<- lactacet.comp[lactacet.comp$Season=="spring",]
lactacet.comp.autumn<- lactacet.comp[lactacet.comp$Season=="autumn",]

lm_lact_spring<- lm(data=lactacet.comp.spring, log.Hypothetical.Lacto ~ Treatment)
diagnostics.plot(lm_lact_spring)
summary(lm_lact_spring)

#run two sample t-test for unequal variances
t.test(log.Hypothetical.Lacto ~ Treatment, data=lactacet.comp.spring, var.equal=TRUE)

lm_lact_autumn<- lm(data=lactacet.comp.autumn, log.Hypothetical.Lacto ~ Treatment)
diagnostics.plot(lm_lact_autumn)
summary(lm_lact_autumn)

#run two sample t-test for unequal variances
t.test(log.Hypothetical.Lacto ~ Treatment, data=lactacet.comp.autumn, var.equal=TRUE)

#rerun the plot so it has the log-transformed axis
lact.comp.summary.1 <- lactacet.comp %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(log.Hypothetical.Lacto, na.rm = TRUE),
    log.Hypothetical.Lacto = mean(log.Hypothetical.Lacto))
lact.comp.summary.1

plot10.1 <- ggplot(lactacet.comp, aes(x = Sample_type, y = log.Hypothetical.Lacto)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, color=c("#FF9966","#FF9966","#FF9966","#FF9966"), fill=c("#FFFFFF","#FF9966","#FFFFFF","#FF9966"), size=1)+
  geom_errorbar( aes(x=Sample_type, ymin=log.Hypothetical.Lacto-sd, ymax=log.Hypothetical.Lacto+sd), data=lact.comp.summary.1, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.03)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "", y="ng Bact DNA in starter or yogurt")+
  theme(axis.text.x = element_text(size=10, face="bold", colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot10.1


#selects ants from the data
ant.data<-qpcr.data[(qpcr.data$Ant_or_Yog_or_Ctl=="ant"),]
#take out dehydrated, because we dont have them for Formica rufa
ant.data1<- ant.data[!ant.data$Treatment=="dehyd",]
#remove the frozen ant
ant.data2<- ant.data1[!ant.data1$Sample_type=="fr_ant_Fp_s",]

#do a quick relevel to make the plot
ant.data3<- ant.data2 %>%
  mutate(Sample_type = fct_relevel(Sample_type, 
                                   "li_ant_Fp_s","li_ant_Fp_a", "li_ant_Fr_s","li_ant_Fr_a"))  
#quick plot
ggplot(ant.data3, aes(x = Sample_type, y = Ng_Bact_DNA)) +
  geom_dotplot(binaxis = "y", stackdir = "center")

#big plot for ant bacterial loads

#quick stats
ant.summary <- ant.data3 %>%
  group_by(Sample_type) %>%
  summarise(
    sd = sd(Ng_Bact_DNA, na.rm = TRUE),
    Ng_Bact_DNA = mean(Ng_Bact_DNA))
ant.summary

#formatted ant plot
plot10 <- ggplot(ant.data3, aes(x = Sample_type, y = Ng_Bact_DNA)) +
  geom_bar(position='dodge', stat='summary', fun='mean', width=0.7, fill=c("#CCCCCC","#CCCCCC","#666666","#666666"))+
  geom_errorbar( aes(x=Sample_type, ymin=Ng_Bact_DNA-sd, ymax=Ng_Bact_DNA+sd), data=ant.summary, width = 0.4)+
  geom_dotplot(binaxis = "y", stackdir = "center", fill="black", binwidth = 0.08)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x= "Ants", y="ng bacterial DNA")+
  scale_x_discrete(labels=c("li_ant_Fp_s"="Spring","li_ant_Fp_a"="Autumn", "li_ant_Fr_s"="Spring","li_ant_Fr_a"="Autumn"))+ 
  theme(axis.text.x = element_text(size=10, colour = "black", angle=45, hjust=0.95),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.line = element_line(colour = "black", linewidth  = 0.7, linetype = "solid"), 
        axis.ticks = element_line(colour = "black",linewidth  = 0.7, linetype = "solid"),
        axis.title.y = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"),
        plot.title = element_text(size=14, face= "bold", colour= "black" , hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot10

##########
#done
##########