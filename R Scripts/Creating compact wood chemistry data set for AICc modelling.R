rm(list=ls())

library(permute)
library(lattice)
library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(data.table)
library(devtools)
library(pairwiseAdonis)
library(goeveg) 

####Creating compact (n=59) data file for AICc modelling####

ENV<-read.table(file="BigSampleMapWood2.csv", sep=",",head=TRUE)

### Summarize down to just one value for each variable for each unique combination of things... 

ENV <- melt(ENV, id.vars= c("Region", "Site", "DecayClass", "WoodType","TreeSpecies"))

ENV<- dcast(ENV, Region + Site + DecayClass  + WoodType + TreeSpecies ~ variable, value.var="value", fun.aggregate = mean)

# Move grdiname to rownames
#rownames(ENV)<-ENV$grdiname
rownames(ENV)<-paste(ENV$Region, ENV$Site, ENV$DecayClass,ENV$WoodType, ENV$TreeSpecies, sep="_")

#remove outliers
remove <- c("Northeast_2_1_hardwood_Trembling aspen", "Northwest_3_2_hardwood_Trembling aspen", "Northwest_3_2_softwood_Jack pine", "Northwest_2_1_hardwood_Trembling aspen")
# remove outliers (samples with less than 5th percentile ESVs)
ENV <- ENV[!rownames(ENV) %in% remove,]


#Renaming weird column names
names(ENV)[names(ENV) == "SiteID"] <- "SiteID"
names(ENV)[names(ENV) == "SampleNo"] <- "SampleNumber"
names(ENV)[names(ENV) == "DecayClass"] <- "DecayClass"
names(ENV)[names(ENV) == "mg.CO2.g.d"] <- "mgCO2.g.d"
names(ENV)[names(ENV) == "TC.."] <- "TotalCarbon"
names(ENV)[names(ENV) == "TKN.."] <- "TotalNitrogen"
names(ENV)[names(ENV) == "P.mg.kg"] <- "Phosphorous"
names(ENV)[names(ENV) == "K.mg.kg"] <- "Potassium"
names(ENV)[names(ENV) == "Ca.mg.kg"] <- "Calcium"
names(ENV)[names(ENV) == "Mg.mg.kg"] <- "Magnesium"
names(ENV)[names(ENV) == "Mn.mg.kg"] <- "Manganese"
names(ENV)[names(ENV) == "C.N"] <- "C.N"
names(ENV)[names(ENV) == "N.P"] <- "N.P"
names(ENV)[names(ENV)=="Moisture.Content...."]<-"MoistureContent"

ENV$DecayClass<-as.factor(ENV$DecayClass)
ENV$Site<-as.factor(ENV$Site)

#write.csv(ENV,"compactenvapril2021.csv",row.names = T)

###################################################################