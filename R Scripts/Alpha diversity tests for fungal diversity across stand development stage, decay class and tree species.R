rm(list=ls())

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)
library(data.table) # rowname to first col
library(devtools)
library(pairwiseAdonis)
library(goeveg) # scree
library(permute)
library(lattice)
library(dunn.test)

###################################################################
# This Script was co-written by Teresita M. Porter (https://github.com/terrimporter)
###################################################################
# Processing taxonomy file 
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
pivot <- dcast(C, Region + SiteID + DecayClass + WoodType +TreeSpecies ~ GlobalESV, value.var="ESVsize", fun.aggregate = sum)

# Move SiteID + DecayClass + WoodType to row names, then remove
rownames(pivot)<-paste(pivot$Region, pivot$SiteID, pivot$DecayClass, pivot$WoodType, pivot$TreeSpecies, sep="_")
pivot <- pivot[,-c(1:5)]

# Remove any columns with only zeros
pivot_notnull<-pivot[,colSums(pivot) !=0]

# Remove any rows with only zeros
pivot_notnull2 <- pivot_notnull[rowSums(pivot_notnull) !=0,]

# Calculate 15th percentile for rrarefy function
pivot_15percentile<-quantile(rowSums(pivot_notnull2), prob=0.15)

# Set random seed for rarefaction
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
df <- rrarefy(pivot_notnull2, sample=pivot_15percentile)

# Convert to presence-absence matrix
df[df>0] <-1

# Check row totals, throw out samples with less than x ESVs
Z <- sort(rowSums(df))

# Calculate 5th percentile of ESVs
Z_5percentile<-quantile(rowSums(df), prob=0.05)
# 55.3 ESVs

# Remove FreshJunction_1_hardwood, ESCAPE2_2_hardwood, ESCAPE2_2_softwood
remove <- c("Northeast_2_1_hardwood_Trembling aspen", "Northwest_3_2_hardwood_Trembling aspen", "Northwest_3_2_softwood_Jack pine")

# Remove outliers (samples with less than 5th percentile ESVs)
df.2 <- df[!rownames(df) %in% remove,]

# Run diveristy metrics on df.2
H<-diversity(df.2,"shannon")
simp<-diversity(df.2,"simpson")
esvnumbers<-specnumber(df.2)
J <- H/log(specnumber(df.2))

# Create data frame with ESV richness 
esv.vibes<-data.frame(result=esvnumbers, stat="specnumber")

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

# Make one data frame with esv.vibes to do ESV richness 
mergedrichness <- merge(esv.vibes, sample_df, by="row.names")
colnames(mergedrichness)[colnames(mergedrichness)=="Row.names"] <- "SampleName"

# Create factors
mergedrichness$Region<-factor(mergedrichness$Region, levels = c("Northwest","Northeast"))
mergedrichness$SiteID<-factor(mergedrichness$SiteID, levels=c("1","2","3"))
mergedrichness$DecayClass<-factor(mergedrichness$DecayClass, levels=c("1","2","3","4","5"))
mergedrichness$WoodType<-factor(mergedrichness$WoodType, levels=c("softwood","hardwood"))
mergedrichness$TreeSpecies<-factor(mergedrichness$TreeSpecies, levels=c("Black spruce","Jack pine","Trembling aspen"))

# ESV Richness 
ESVrichnessDC<-mergedrichness %>% group_by(DecayClass) %>% summarise(mean = mean(result), n = n(),sd=sd(result),variance=(sd(result)/mean(result)*100),.groups = 'drop')
ESVrichnessTreeSpecies<-mergedrichness %>% group_by(TreeSpecies) %>% summarise(mean = mean(result), n = n(),sd=sd(result),variance=(sd(result)/mean(result)*100),.groups = 'drop')
ESVrichnessRegion<-mergedrichness %>% group_by(Region) %>% summarise(mean = mean(result), n = n(),d=sd(result),variance=(sd(result)/mean(result)*100),.groups = 'drop')

ESVtallyDC<-mergedrichness %>% group_by(DecayClass) %>% tally(result)
ESVtallyTreeSpecies<-mergedrichness %>% group_by(TreeSpecies) %>% tally(result)
ESVtallyRegion<-mergedrichness %>% group_by(Region) %>% tally(result)


# Make one data frame with J to do ESV evenness 
mergedevenness<-merge(J,sample_df,by="row.names")
colnames(mergedevenness)[colnames(mergedevenness)=="Row.names"] <- "SampleName"

# Create factors
mergedevenness$Region<-factor(mergedevenness$Region, levels = c("Northwest","Northeast"))
mergedevenness$SiteID<-factor(mergedevenness$SiteID, levels=c("1","2","3"))
mergedevenness$DecayClass<-factor(mergedevenness$DecayClass, levels=c("1","2","3","4","5"))
mergedevenness$WoodType<-factor(mergedevenness$WoodType, levels=c("softwood","hardwood"))
mergedevenness$TreeSpecies<-factor(mergedevenness$TreeSpecies, levels=c("Black spruce","Jack pine","Trembling aspen"))

#Evenness
ESVevennessDC<-mergedevenness %>% group_by(DecayClass) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')
ESVevennessTreeSpecies<-mergedevenness %>% group_by(TreeSpecies) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')
ESVevennessTreeSpecies<-mergedevenness %>% group_by(Region) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')

# Make one data frame with H for Shannon's diversity 
mergedshannon<-merge(H,sample_df,by="row.names")
colnames(mergedshannon)[colnames(mergedshannon)=="Row.names"] <- "SampleName"

# Create factors
mergedshannon$Region<-factor(mergedshannon$Region, levels = c("Northwest","Northeast"))
mergedshannon$SiteID<-factor(mergedshannon$SiteID, levels=c("1","2","3"))
mergedshannon$DecayClass<-factor(mergedshannon$DecayClass, levels=c("1","2","3","4","5"))
mergedshannon$WoodType<-factor(mergedshannon$WoodType, levels=c("softwood","hardwood"))
mergedshannon$TreeSpecies<-factor(mergedshannon$TreeSpecies, levels=c("Black spruce","Jack pine","Trembling aspen"))

#Shannon's diversity
ESVshannonDC<-mergedshannon %>% group_by(DecayClass) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')
ESVshannonTreeSpecies<-mergedshannon %>% group_by(TreeSpecies) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')
ESVshannonRegion<-mergedshannon %>% group_by(Region) %>% summarise(mean = mean(x), n = n(),sd=sd(x),variance=(sd(x)/mean(x)*100),.groups = 'drop')

#ANOVAS + post-hoc testing
richness<-lm(result~Region+SiteID+TreeSpecies+DecayClass+TreeSpecies*DecayClass, data=mergedrichness)
anova(richness)
richness<-aov(result~DecayClass, data=mergedrichness)
TukeyHSD(richness, "DecayClass", ordered=TRUE)

evenness<-lm(x~Region+SiteID+TreeSpecies+DecayClass+TreeSpecies*DecayClass, data=mergedevenness)
anova(evenness)

shannons<-lm(x~Region+SiteID+TreeSpecies+DecayClass+TreeSpecies*DecayClass, data=mergedshannon)
anova(shannons)
shannonst<-aov(x~DecayClass, data=mergedshannon)
TukeyHSD(shannonst, "DecayClass", ordered=TRUE)







