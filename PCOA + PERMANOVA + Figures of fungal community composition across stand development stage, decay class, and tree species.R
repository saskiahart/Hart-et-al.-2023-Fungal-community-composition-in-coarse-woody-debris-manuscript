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
library(ape)

citation('goeveg')
citation('vegan')
citation('pairwiseAdonis')
citation('MuMIn')
citation('dunn.test')
citation('ape')

###################################################################
# This Script was written by Teresita M. Porter (https://github.com/terrimporter)
###################################################################
# Read in metadata for PCoA with environmental overlay 
###################################################################

ENV<-read.table(file="BigSampleMapWood2.csv", sep=",",head=TRUE)

### Summarize down to just one value for each variable for each unique combination of things... 

ENV <- melt(ENV, id.vars= c("Region", "Site", "DecayClass", "TreeSpecies","WoodType"))

ENV<- dcast(ENV, Region + Site + DecayClass + TreeSpecies + WoodType ~ variable, value.var="value", fun.aggregate = mean)

# Move grdiname to rownames
rownames(ENV)<-ENV$grdiname
rownames(ENV)<-paste(ENV$Region, ENV$Site, ENV$DecayClass, ENV$WoodType, ENV$TreeSpecies, sep="_")

# Rename to nicer names
names(ENV)[names(ENV) == "TC.."] <- "TotalCarbon"
names(ENV)[names(ENV) == "TKN.."] <- "TotalNitrogen"
names(ENV)[names(ENV) == "P.mg.kg"] <- "Phosphorous"
names(ENV)[names(ENV) == "K.mg.kg"] <- "Potassium"
names(ENV)[names(ENV) == "Ca.mg.kg"] <- "Calcium"
names(ENV)[names(ENV) == "Mg.mg.kg"] <- "Magnesium"
names(ENV)[names(ENV) == "Mn.mg.kg"] <- "Manganese"
names(ENV)[names(ENV) == "CN"] <- "C:N"
names(ENV)[names(ENV) == "NP"] <- "N:P"
names(ENV)[names(ENV) == "mgCO2.g.d"] <- "mgCO2/g/d"
names(ENV)[names(ENV)=="Moisture.Content...."]<-"MoistureContent"
names(ENV)[names(ENV)=="Density..g.cm3."]<-"Density"

# Make sure Decay Class and Sites are factors (as they are labelled "1, 2, 3, etc.")
ENV$DecayClass<-as.factor(ENV$DecayClass)
ENV$Site<-as.factor(ENV$Site)

# remove FreshJunction_1_hardwood, ESCAPE2_2_hardwood, ESCAPE2_2_softwood
remove <- c("Northeast_2_1_hardwood_Po", "Northwest_3_2_hardwood_Po","Northwest_2_1_softwood_Pj")

# remove outliers (samples with less than 5th percentile ESVs)
ENV2 <- ENV[!rownames(ENV) %in% remove,]

###################################################################
# Processing taxonomy file + creating matrix for NMDS
###################################################################

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
C <- merge(x = A5, y = B, by = "SampleNumber", all.x = TRUE)

# Pivot to make matrix for vegan, sites in rows, ESVs in columns
pivot <- dcast(C, Region + SiteID + DecayClass + WoodType +TreeSpecies  ~ GlobalESV, value.var="ESVsize", fun.aggregate = sum)

# Move SiteID + DecayClass + WoodType to row names, then remove
rownames(pivot)<-paste(pivot$Region, pivot$SiteID, pivot$DecayClass, pivot$WoodType, pivot$TreeSpecies, sep="_")
pivot <- pivot[,-c(1:5)]

# remove any columns with only zeros
pivot_notnull<-pivot[,colSums(pivot) !=0]

# remove any rows with only zeros
pivot_notnull2 <- pivot_notnull[rowSums(pivot_notnull) !=0,]

#calculate 15th percentile for rrarefy function
pivot_15percentile<-quantile(rowSums(pivot_notnull2), prob=0.15)

# Set random seed for rarefaction
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df <- rrarefy(pivot_notnull2, sample=pivot_15percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Check row totals, throw out samples with less than x ESVs
Z <- sort(rowSums(df))

# calculate 5th percentile of ESVs
Z_5percentile<-quantile(rowSums(df), prob=0.05)
# 55.3 ESVs

# identifying outliers (samples with less than 5th percentile ESVs)
remove <- c("Northeast_2_1_hardwood_Trembling aspen", "Northwest_3_2_hardwood_Trembling aspen", "Northwest_3_2_softwood_Jack pine", "Northwest_2_1_hardwood_Trembling aspen")

# remove outliers (samples with less than 5th percentile ESVs)
df.2 <- df[!rownames(df) %in% remove,]

# Scree plots to determine number of dimensions to use for PCoA
pdf("ScreeJan2020.pdf")
dimcheckMDS(df)
dev.off()

# Creating distance matrix for PCoA
dist1 <- vegdist(df.2, method="bray")

# PCoA
PCOA <- pcoa(dist1)
barplot(PCOA$values$Relative_eig[1:10])
biplot.pcoa(PCOA)
PCOAaxes<-PCOA$vectors[,c(1,2)]
par(mfrow=c(1,2))
biplot.pcoa(PCOA)
plot(PCOA)
#write.table(PCOAaxes, file = "PCOA2.csv", sep = ",", quote = FALSE, row.names = T,col.names = T)


# Make sure ENV table has only samples that are in microbial data frame (df.2)
ENV2 <- ENV[row.names(ENV) %in% row.names(df.2), ]

# check orders match! They DO!
order(rownames(ENV2), rownames(df.2))

# Run PCoA with environmental data included 
envdata<-envfit(PCOAdf, ENV2[, !colnames(ENV2) %in% c("SiteID")])


# Extract data from envdata for plotting on PCoA plots 
SegmentData <- data.frame(envdata$vectors$arrows)
SegmentData$label <- rownames(SegmentData)
SegmentData1<-SegmentData[-c(3,5,8,9,12,13),]

SegmentData <- data.frame(envdata$vectors$pvals)

# Extract data from envdata for plotting on PCoA plots 
SegmentData <- data.frame(envdata$vectors$arrows)
SegmentData$label <- rownames(SegmentData)
SegmentData1<-SegmentData[-c(3,5,6,8,9,12,13),]

# Back to taxonomy data:
# Create grouping matrix for samples by grabbing row names from above matrix
sample_df<-data.frame(row.names(df.2))

# Rename the column
names(sample_df)<-"SampleName"

# Copy column to row names
row.names(sample_df)<-sample_df$SampleName

# Split first column into their own fields
sample_df[,2:6]<-do.call('rbind', strsplit(as.character(sample_df$SampleName),'_',fixed=TRUE))

# Remove first column
sample_df<-sample_df[,-1]

# Rename columns
names(sample_df)<-c("Region", "SiteID","DecayClass","WoodType","TreeSpecies")

# Grab sites & scores from NMDS output
PCOAdf<- data.frame(PCOA$vectors[,c(1,2)])
PCOAaxes<-PCOA$vectors[,c(1,2)]


# Put it all in one df for ggplot
mergedB <- merge(PCOAdf, sample_df, by="row.names")
colnames(mergedB)[colnames(mergedB)=="Row.names"] <- "SampleName"


# Create factors
mergedB$Region<-factor(mergedB$Region, levels = c("Northwest","Northeast"))
mergedB$SiteID<-factor(mergedB$SiteID, levels=c("1","2","3"))
mergedB$DecayClass<-factor(mergedB$DecayClass, levels=c("1","2","3","4","5"))
mergedB$WoodType<-factor(mergedB$WoodType, levels=c("softwood","hardwood"))
mergedB$TreeSpecies<-factor(mergedB$TreeSpecies, levels=c("Black spruce","Jack pine","Trembling aspen"))

#write.table(mergedB, file = "mergedB.csv", sep = ",", quote = FALSE, row.names = T,col.names = T)

# Ready to plot PCoA plots using ggplot2

# Compile coord for convex hulls for Decay Class plots
chulls.2 <- ddply(mergedB, .(DecayClass), function(mergedB) mergedB[chull(mergedB$Axis.1, mergedB$Axis.2), ])

# Create scatter plot for Decay Class
p2<-ggplot(data=mergedB, aes(x=Axis.1, y=Axis.2)) + 
  ggtitle("b) Decay Class") +
  geom_polygon(data=chulls.2, aes(x=Axis.1, y=Axis.2, fill=DecayClass), alpha=0.5) +
  geom_point(data=mergedB, aes(x=Axis.1, y=Axis.2, color=DecayClass, shape=DecayClass)) +
  scale_fill_manual(values=c("#003D30","#005745","#009175","#00AF8E","#86FFDE")) +
  scale_color_manual(values=c("#003D30","#005745","#009175","#00AF8E","#86FFDE")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom")

plot(p2)

###########

# Compile coord for convex hulls for Stand Development Stage (aka Region) plots
chulls.5 <- ddply(mergedB, .(Region), function(mergedB) mergedB[chull(mergedB$Axis.1, mergedB$Axis.2), ])
 
# Create scatter plot for Stand Development Stage
p5<-ggplot(data=mergedB, aes(x=Axis.1, y=Axis.2)) + 
  ggtitle("c) Stand Development Stage") +
  geom_point(data=mergedB, aes(x=Axis.1, y=Axis.2, color=Region, shape=Region)) +
  geom_polygon(data=chulls.5, aes(x=Axis.1, y=Axis.2, colour=Region), fill=NA, alpha=0.5) +
  scale_fill_manual(values=c("#003D30","#00CBA7")) +
  scale_color_manual(values=c("#003D30","#00CBA7")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom")

plot(p5)



##########
###
# Compile coord for convex hulls for Tree Species plots
chulls.6 <- ddply(mergedB, .(TreeSpecies), function(mergedB) mergedB[chull(mergedB$Axis.1, mergedB$Axis.2), ])

# Create scatter plot for Tree Species
p6<-ggplot(data=mergedB, aes(x=Axis.1, y=Axis.2)) + 
  ggtitle("a) Tree Species") +
  geom_polygon(data=chulls.6, aes(x=Axis.1, y=Axis.2, fill=TreeSpecies),alpha=0.5) +
  geom_point(data=mergedB, aes(x=Axis.1, y=Axis.2, color=TreeSpecies, shape=TreeSpecies)) +
  scale_fill_manual(values=c("#003D30","#00735C","#00CBA7")) +
  scale_color_manual(values=c("#003D30","#00735C","#00CBA7")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom")

plot(p6)

######

# Create scatterplot for envdata
envplot<-ggplot(data=mergedB, aes(x=NMDS1, y=NMDS2)) + 
  ggtitle("d) Best predictors of fungal community composition") +
  geom_point(data=mergedB, aes(x=NMDS1, y=NMDS2), colour="#003D30",pch=0) +
  geom_segment(data=SegmentData, aes(x=0, y=0,xend=NMDS1, yend=NMDS2), 
               color="black",
               arrow=arrow(length = unit(0.1,"cm")))+
  geom_text(data=SegmentData, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=label), 
            color="white", size=1)+
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom")

plot(envplot)


# Time for Adonis test to do PERMANOVA 

# Only keep the samples that are in df
ENV <- data.frame(rownames(df.2))
names(ENV) <- "SampleName"

# separate into own columns
ENV[,2:6]<-data.frame(do.call('rbind', strsplit(as.character(ENV$SampleName),'_',fixed=TRUE)))
names(ENV)[2:6] <- c("Region", "SiteID", "DecayClass", "WoodType","TreeSpecies")

# # Remove first column
ENV<-ENV[,-1]

# Create factors
ENV$Region<-factor(ENV$Region, levels=c("Northeast","Northwest"))
ENV$SiteID<-factor(ENV$SiteID, levels=c("1","2","3"))
ENV$DecayClass<-factor(ENV$DecayClass, levels=c("1","2","3","4","5"))
ENV$WoodType<-factor(ENV$WoodType, levels=c("softwood","hardwood"))
ENV$TreeSpecies<-factor(ENV$TreeSpecies, levels=c("Black spruce","Jack pine","Trembling aspen"))

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray-Curtis dissimilarity 
sor <- vegdist(df.2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#Region
bd_region<-betadisper(sor, factor(ENV$Region))
# SiteID
bd_site<-betadisper(sor, factor(ENV$SiteID))
# DecayClass
bd_decayclass<-betadisper(sor, factor(ENV$DecayClass))
# WoodType
bd_woodtype<-betadisper(sor, factor(ENV$WoodType))
# WoodType
bd_treespecies<-betadisper(sor, factor(ENV$TreeSpecies))

# check for heterogeneity of beta dispersions within groups
anova(bd_region) #n/s
anova(bd_site) # n/s
anova(bd_decayclass) # 0.001975 ** missing one DC1 sample so not balanced design
anova(bd_woodtype) # n/s
anova(bd_treespecies)

pdf("BetaDispersionJan2020.pdf")
par(mfrow=c(2,2))
boxplot(bd_site, las=2)
mtext("Sites", side=3, adj=0, line=1.2)
boxplot(bd_decayclass)
mtext("Decay Class", side=3, adj=0, line=1.2)
boxplot(bd_woodtype)
mtext("Wood Type", side=3, adj=0, line=1.2)
boxplot(bd_treespecies)
mtext("TreeSpecies", side=3, adj=0, line=1.2)
dev.off()

# Shephards curve and goodness of fit calcs
# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplotJan2020.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

# Use ADONIS to test significance of groupings 

#Use strata so randomizations occur within each Stand Development Stage (aka Region)
adonis(sor~SiteID+TreeSpecies+DecayClass+TreeSpecies*DecayClass, data=ENV, permutations=999, strata=ENV$Region)

#Additionally testing the effect of sites nested within Stand Development Stage (aka Region)
adonis(sor~Region/SiteID, data=ENV, permutations=999)
 

# Merging scatterplots into one big plot (Figure 3 in manuscript)
lay<-rbind(c(1,2), 
           c(3,4))

g<-grid.arrange(p6, p2, p5, envplot, layout_matrix=lay)

g<-grid.arrange(p5, p6, p2, layout_matrix=lay)
# Save plot
ggsave("NMDStestwenv.pdf", g, width = 10, height = 9, units = c("in"))

ggsave("PCOAbig.pdf", g, width = 10, height = 9, units = c("in"))
