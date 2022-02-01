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
library(mvpart)

###################################################################

# Processing of taxonomy data was written by Teresita M. Porter (https://github.com/terrimporter)
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


# Left join (maps info from mapping file to Sample Number)
C <- merge(x = A5, y = B, by = "SampleNumber", all.x = TRUE)

# Pivot to make matrix for vegan, sites in rows, ESVs in columns
pivot <- dcast(C, SiteID + DecayClass + WoodType + TreeSpecies + Region ~ GlobalESV, value.var="ESVsize", fun.aggregate = sum)

# Move SiteID + DecayClass + WoodType to row names, then remove
rownames(pivot)<-paste(pivot$SiteID, pivot$DecayClass, pivot$WoodType, pivot$TreeSpecies, pivot$Region, sep="_")
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

#calculate 5th percentile of ESVs
Z_5percentile<-quantile(rowSums(df), prob=0.05)
# 55.3 ESVs

# identify outliers
remove <- c("2_1_hardwood_Trembling aspen_Northeast", "3_2_hardwood_Trembling aspen_Northwest", "3_2_softwood_Jack pine_Northwest", "2_1_hardwood_Trembling aspen_Northwest")

# remove outliers (samples with less than 5th percentile ESVs)
df2 <- df[!rownames(df) %in% remove,]

# read in mapping file
condensedmap<-read.table(file="condensemap3.csv", head=TRUE, sep=",")

# merge mapping and data files
MRTdf<- cbind(condensedmap, df2)
rownames(MRTdf) <- MRTdf$SampleName
MRTdf <- MRTdf[ -c(1) ]

MRTdf$DecayClass<-as.factor(MRTdf$DecayClass)


#### Making the Multivariate regression table used functions (Tokas and PrettyMRT) developed by Timothy Work and Jenna Jacobs. For details on these functions please contact directly.
# making MRT
CWDMrt<-mvpart(form = data.matrix(MRTdf[,6:11183]) ~Region + DecayClass + TreeSpecies, 
                MRTdf, xv="pick")

# using Tokas function to create variance table (Table 4 in manuscript). Tokas was developed by Jenna Jacobs and Timothy Working. 
CWDMrt_Tab<-tokas(CWDMrt)

#export MRT variance table as a csv file
#write.table(CWDMrt_Tab, file = "MRTtest.csv", sep = ",", quote = FALSE, row.names = T,col.names = T)


# create a label vector from the MRT to make a graph (labels are for end points)
rownames(CWDMrt$frame)
labv_CWDMrt=c("root","Jack pine,Black spruce","DC 1,2,3","Northeast","DC 1,2","DC 3","Northwest","DC 1,2","DC 2","DC 1","DC 3","DC 4,5","Black spruce","DC 5","DC 4","Jack pine","DC 5","Northeast","Northwest","DC 4","Trembling aspen","DC 1,2","DC 3,4,5","Northeast","DC 3,4","DC 3","DC 4","DC 5","Northwest","DC 4,5","DC 4","DC 5","DC 3")


# using PrettyMRT function to create a MRT with custom vector labels. PrettyMRT function was developed by Timothy Work. 
PrettyMRT(LisaMrt, labv_CWDMrt, c(5,1,1,1),0.2)

