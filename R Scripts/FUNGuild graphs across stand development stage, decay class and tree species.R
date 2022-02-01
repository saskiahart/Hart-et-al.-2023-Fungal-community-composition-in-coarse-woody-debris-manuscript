rm(list=ls())

library(readr)
library(vegan)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyverse)
library(RColorBrewer)


# Read in tables from FUNGuild (they were big and downloaded in 2 files)
FG1 <- read_tsv(file="funGUILD.1.guilds.txt")

FG2<-read_tsv(file="funGUILD.2.guilds (1).txt")

# Making one big df
fundat<-bind_rows(FG1,FG2)

# Melting, specifying we don't want OTU_ID or taxonomic columns melted
funmelt <- melt(fundat, 
                id.vars=c("OTU_ID", colnames(fundat)[131:ncol(fundat)]), 
                variable.name = "SampleNo", 
                value.name ="abundance")

# read in mapping file
mapping <- read.csv("sampleMapregionspecies.csv")
names(mapping)[1] <- "SampleNo" 

#merging fungal data and mapping data together
funmelt <- merge(funmelt, mapping, by="SampleNo")


#creating df grouped by trait with mean & standard error & richness for plotting 
functionalsummary<- funmelt %>% 
  group_by(Trait) %>% 
  dplyr::summarise(mean.abundance = mean(abundance), Standard.Error = sd(abundance)/sqrt(n()), Richness = length(unique(Taxon)))

test<-funmelt %>% group_by(Trait, OTU_ID) %>% tally(abundance)


#Stacked Bar Plot w Decay Class
functionalsummarydc12<- funmelt %>% 
  group_by(Trait,DecayClass) %>% 
  dplyr::summarise(mean.abundance = mean(abundance), Standard.Error = sd(abundance)/sqrt(n()), Richness = length(unique(Taxon))) %>%
  # get rid of uninformative ones
  filter(!Trait %in% c("-", "Null", "NULL"))


abundanceplottestdc12 <- ggplot(functionalsummarydc12, aes(y=mean.abundance, x=DecayClass, fill = factor(Trait, levels=c("Brown Rot-White Rot","White Rot","Brown Rot","Soft Rot","Blue-Staining","Poisonous","Hypogeous"))))+
  geom_bar(position="fill", width=0.6, stat="identity", colour="black")+
  scale_fill_brewer(palette="Set1")+
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))

abundanceplottestdc12 + labs(y="Proportion (%)", fill="Trait") 

#Stacked Bar Plot w Tree Species
functionalsummarytree<- funmelt %>% 
  group_by(Trait,TreeSpecies) %>% 
  dplyr::summarise(mean.abundance = mean(abundance), Standard.Error = sd(abundance)/sqrt(n()), Richness = length(unique(Taxon))) %>%
  # get rid of uninformative ones
  filter(!Trait %in% c("-", "Null", "NULL"))


abundanceplottesttrees <- ggplot(functionalsummarytree, aes(y=mean.abundance, x=TreeSpecies, fill = factor(Trait, levels=c("Brown Rot-White Rot","White Rot","Brown Rot","Soft Rot","Blue-Staining","Poisonous","Hypogeous"))))+
  geom_bar(position="fill", width=0.6, stat="identity", colour="black")+
  scale_fill_brewer(palette="Set1")+
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))

abundanceplottesttrees + labs(y="Proportion (%)", fill="Trait") 

#Stacked Bar Plot w Stand Development Stage (aka Region)
functionalsummaryregion12<- funmelt %>% 
  group_by(Trait,Region) %>% 
  dplyr::summarise(mean.abundance = mean(abundance), Standard.Error = sd(abundance)/sqrt(n()), Richness = length(unique(Taxon))) %>%
  # get rid of uninformative ones
  filter(!Trait %in% c("-", "Null", "NULL"))


abundanceplottestregion12 <- ggplot(functionalsummaryregion12, aes(y=mean.abundance, x=Region, fill = factor(Trait, levels=c("Brown Rot-White Rot","White Rot","Brown Rot","Soft Rot","Blue-Staining","Poisonous","Hypogeous"))))+
  geom_bar(position="fill", width=0.6, stat="identity", colour="black")+
  scale_fill_brewer(palette="Set1")+
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))

abundanceplottestregion12 + labs(y="Proportion (%)", fill="Trait") 


RegionTally<-funmelt %>% group_by(Region) %>% tally(abudance)

#Arranging plots into one figure (Figure 5 in manuscript)

paneltotal <- ggarrange( abundanceplottestregion12 + rremove("xlab") +rremove("legend.title"), abundanceplottesttrees + rremove("xlab") +rremove("legend.title") +rremove("ylab"), abundanceplottestdc12 + rremove("xlab") +rremove("legend.title")+rremove("ylab"),
                         common.legend = TRUE, legend ="right",
                         ncol = 3, nrow = 1) +
  theme(plot.margin = margin(t = 2, unit ="cm"))
paneltotal

