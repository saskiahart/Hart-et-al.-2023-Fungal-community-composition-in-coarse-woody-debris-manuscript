rm(list=ls())

library(permute)
library(lattice)
library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table) # rowname to first col
library(devtools)
library(pairwiseAdonis)
library(goeveg) # scree
library(colorRamps)
library(dplyr)
library(scales)
library(RColorBrewer)
library(wesanderson)
library(dplyr)

##############
# This Script was co-written by Teresita M. Porter (https://github.com/terrimporter)
#############

# Read in taxonomy file
A<-read.table(file="Saskia_fITS74R_fITS94R_taxonomy.csv", head=TRUE, sep=",")

# Split GlobalESV to get amplicon as separate column
A2 <- data.frame(do.call('rbind', strsplit(as.character(A$GlobalESV),'_', fixed = TRUE)))

# Add amplicon column back to A2
A3 <-cbind(A,A2$X1)
names(A3)[length(names(A3))] <- "amplicon"

# Split SampleName to get SampleNumber as separate column
A4 <- data.frame(do.call('rbind', strsplit(as.character(A$SampleName), '_', fixed = TRUE)))

# Add SampleNumber column back to A3
A5 <- cbind(A3, A4$X3)
names(A5)[length(names(A5))] <- "SampleNumber"

# Read in sample mapping file
B <- read.table(file="sampleMapregionspecies.csv", head=TRUE, sep=",")
names(B)[1] <- "SampleNumber" 

# Left join (maps info from mapping file to Sample Number)
C2 <- merge(x = A5, y = B, by = "SampleNumber", all.x = TRUE)


#Filtering out unnecessary columns
C2 <- C2[ -c(3,5,6,7,8,9,10,11,13,16,19,22,25,28)]

#Filtering by 0.5 bp cutoff 
C2<-C2 %>% filter(pBP >= 0.5, cBP>=0.5, oBP>=0.5, fBP>=0.5, gBP>=0.5, sBP>=0.5)

# Get all unique ESV counts for region, decay class, tree species, etc
length(unique(C2$GlobalESV[C2$DecayClass=="DecayClass"]))

length(unique(C2$Species[C2$Region=="Northeast"]))
length(unique(C2$Species[C2$Region=="Northwest"]))

length(unique(C2$GlobalESV[C2$TreeSpecies=="Jack pine"]))
length(unique(C2$GlobalESV[C2$TreeSpecies=="Black spruce"]))
length(unique(C2$GlobalESV[C2$TreeSpecies=="Trembling aspen"]))

length(unique(C2$GlobalESV[C2$DecayClass=="1"]))
length(unique(C2$GlobalESV[C2$DecayClass=="2"]))
length(unique(C2$GlobalESV[C2$DecayClass=="3"]))
length(unique(C2$GlobalESV[C2$DecayClass=="4"]))
length(unique(C2$GlobalESV[C2$DecayClass=="5"]))

length(unique(C2$GlobalESV))

length(unique(C2$ESVsize))

#####

#Pivoting df for ggplot by Decay Class and Class
ClassDC <- dcast(C2, DecayClass+Class ~ . , value.var="GlobalESV", fun.aggregate=length)
names(ClassDC)<-c("DecayClass","Class", "ESVs")

#Total tally of ESV's in each Decay Class
SiteTallyDC<-ClassDC %>% group_by(DecayClass) %>% tally(ESVs) 

#Making object to remove incerates from list
remove<-c("Incertae_sedis_10","Incertae_sedis_14", "Incertae_sedis_4")

#Dataframe without uncertaes
ClassDCFilter<- ClassDC %>% select(DecayClass,Class,ESVs) %>% filter(!Class %in% remove)

#Melting for ggplot for big graph
phylumTable.longDCC<-melt(ClassDCFilter, id=c("DecayClass", "Class"))

#Plot of confidently ID'd ESVs to Class 
pDCC<-ggplot(data=phylumTable.longDCC, aes(x=phylumTable.longDCC$DecayClass, y=value, fill = factor(Class, levels=c("Agaricomycetes", "Leotiomycetes", "Eurotiomycetes", "Sordariomycetes", "Dothideomycetes", "Fungi_unidentified_1", "Ascomycota_unidentified", "Lecanoromycetes", "Microbotryomycetes", "Tremellomycetes", "Basidiomycota_unidentified", "Saccharomycetes", "Pezizomycetes", "Dacrymycetes", "Orbiliomycetes", "Exobasidiomycetes", "Chytridiomycetes", "Chytridiomycota_unidentified", "Glomeromycota_unidentified", "Cystobasidiomycetes", "Ustilaginomycetes", "Tritirachiomycetes", "Atractiellomycetes", "Agaricostilbomycetes"))))+
  geom_bar(position="fill", width=0.6, stat="identity", colour="black") +
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("#00CED1","#00FF00","#0000FF","#FFFF00","#00FFFF","#FF00FF","#800000","#008000","#800080","#008080",
                             "#8A2BE2","#9ACD32","#4B0082","#BC8F8F", "#FF6347","#87CEEB","#DAA520","#DDA0DD","#D2691E","#2F4F4F",
                             "#DC143C","#7CFC00","#FF1493","#808000"))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))
plot(pDCC)+labs(fill="Class")

#####

#Pivoting df for ggplot by Stand Development Stage (aka Region) and Class
ClassR <- dcast(C2, Region+Class ~ . , value.var="GlobalESV", fun.aggregate=length)
names(ClassR)<-c("Region","Class", "ESVs")

#Total tally of ESV's in each Region
RegionTallyClass<-ClassR %>% group_by(Region) %>% tally(ESVs) 

#Making object to remove incerates from list
remove<-c("Incertae_sedis_10","Incertae_sedis_14", "Incertae_sedis_4")

#Dataframe without uncertaes
ClassRFilter<- ClassR %>% filter(!Class %in% remove)

#Melting for ggplot
classTable.longR<-melt(ClassRFilter, id=c("Region", "Class"))

#Plot of confidently ID'd ESVs to Class
pCR<-ggplot(data=classTable.longR, aes(x=classTable.longR$Region, y=value, fill = factor(Class, levels=c("Agaricomycetes", "Leotiomycetes", "Eurotiomycetes", "Sordariomycetes", "Dothideomycetes", "Fungi_unidentified_1", "Ascomycota_unidentified", "Lecanoromycetes", "Microbotryomycetes", "Tremellomycetes", "Basidiomycota_unidentified", "Saccharomycetes", "Pezizomycetes", "Dacrymycetes", "Orbiliomycetes", "Exobasidiomycetes", "Chytridiomycetes", "Chytridiomycota_unidentified", "Glomeromycota_unidentified", "Cystobasidiomycetes", "Ustilaginomycetes", "Tritirachiomycetes", "Atractiellomycetes", "Agaricostilbomycetes"))))+
  geom_bar(position="fill", width=0.6, stat="identity", colour="black") +
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("#00CED1","#00FF00","#0000FF","#FFFF00","#00FFFF","#FF00FF","#800000","#008000","#800080","#008080",
                             "#8A2BE2","#9ACD32","#4B0082","#BC8F8F", "#FF6347","#87CEEB","#DAA520","#DDA0DD","#D2691E","#2F4F4F",
                             "#DC143C","#7CFC00","#FF1493","#808000"))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))
plot(pCR)+labs(fill="Class")

#####

#Pivoting df for ggplot by Tree Species and class
ClassS <- dcast(C2, TreeSpecies+Class ~ . , value.var="GlobalESV", fun.aggregate=length)
names(ClassS)<-c("TreeSpecies","Class", "ESVs")

#Total tally of ESV's in each Tree Species
SpeciesTallyClass<-ClassS %>% group_by(TreeSpecies) %>% tally(ESVs) 

#Making object to remove incerates from list
remove<-c("Incertae_sedis_10","Incertae_sedis_14", "Incertae_sedis_4")

#Dataframe without uncertaes
ClassSFilter<- ClassS  %>% filter(!Class %in% remove)

#Melting for ggplot
phylumTable.longSC<-melt(ClassSFilter, id=c("TreeSpecies", "Class"))

#Plot of confidently ID'd ESVs to Class 
pSC<-ggplot(data=phylumTable.longSC, aes(x=phylumTable.longSC$TreeSpecies, y=value, fill = factor(Class, levels=c("Agaricomycetes", "Leotiomycetes", "Eurotiomycetes", "Sordariomycetes", "Dothideomycetes", "Fungi_unidentified_1", "Ascomycota_unidentified", "Lecanoromycetes", "Microbotryomycetes", "Tremellomycetes", "Basidiomycota_unidentified", "Saccharomycetes", "Pezizomycetes", "Dacrymycetes", "Orbiliomycetes", "Exobasidiomycetes", "Chytridiomycetes", "Chytridiomycota_unidentified", "Glomeromycota_unidentified", "Cystobasidiomycetes", "Ustilaginomycetes", "Tritirachiomycetes", "Atractiellomycetes", "Agaricostilbomycetes"))))+
  ggtitle("b) Tree Species") +
  geom_bar(position="fill", width=0.6, stat="identity", colour="black") +
  labs(y="Proportion (%)", x="ESVs") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("#00CED1","#00FF00","#0000FF","#FFFF00","#00FFFF","#FF00FF","#800000","#008000","#800080","#008080",
                             "#8A2BE2","#9ACD32","#4B0082","#BC8F8F", "#FF6347","#87CEEB","#DAA520","#DDA0DD","#D2691E","#2F4F4F",
                             "#DC143C","#7CFC00","#FF1493","#808000"))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14))

plot(pSC)+labs(fill="Class")

#####
#Merging ggplot graphs together to make one big graph 
library(ggpubr)
paneltotal <- ggarrange( pCR + rremove("xlab") +rremove("legend.title"), pSC + rremove("xlab") +rremove("legend.title") +rremove("ylab"), pDCC + rremove("xlab") +rremove("legend.title")+rremove("ylab"),
                         common.legend = TRUE, legend ="bottom",
                         ncol = 3, nrow = 1) +
  theme(plot.margin = margin(t = 2, unit ="cm"))

paneltotal + guides(fill=guide_legend(title="Class"))
