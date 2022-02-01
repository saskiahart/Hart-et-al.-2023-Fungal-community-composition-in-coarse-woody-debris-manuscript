rm(list=ls())
library(pastecs)
library(lattice)
library(ggplot2)
library(nortest)
library(tidyverse)
library(MuMIn)
library(emmeans)
library(nlme)
library(car)
library(dplyr)
library(dunn.test)
library(permute)
library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(plyr)
library(data.table) 
library(devtools)


####Testing significance of coarse woody debris physical, biological and chemical characteristics####

#Read in wood chemistry data
CWDChem <- read.csv(file="BigSampleMapWood.csv", header = TRUE, sep = ",")

#Removing N/A from density column 
CWDChem<-na.omit(CWDChem)

#Renaming column names
names(CWDChem)[names(CWDChem) == "SiteID"] <- "SiteID"
names(CWDChem)[names(CWDChem) == "SampleNo"] <- "SampleNumber"
names(CWDChem)[names(CWDChem) == "DecayClass"] <- "DecayClass"
names(CWDChem)[names(CWDChem) == "mg.CO2.g.d"] <- "mgCO2.g.d"
names(CWDChem)[names(CWDChem) == "TC.."] <- "TotalCarbon"
names(CWDChem)[names(CWDChem) == "TKN.."] <- "TotalNitrogen"
names(CWDChem)[names(CWDChem) == "P.mg.kg"] <- "Phosphorous"
names(CWDChem)[names(CWDChem) == "K.mg.kg"] <- "Potassium"
names(CWDChem)[names(CWDChem) == "Ca.mg.kg"] <- "Calcium"
names(CWDChem)[names(CWDChem) == "Mg.mg.kg"] <- "Magnesium"
names(CWDChem)[names(CWDChem) == "Mn.mg.kg"] <- "Manganese"
names(CWDChem)[names(CWDChem) == "C.N"] <- "C.N"
names(CWDChem)[names(CWDChem) == "N.P"] <- "N.P"
names(CWDChem)[names(CWDChem)=="Moisture.Content...."]<-"MoistureContent"
names(CWDChem)[names(CWDChem)=="Density..g.cm3."]<-"Density"
names(CWDChem)[names(CWDChem) == "Cmass"] <- "MassCarbon"
names(CWDChem)[names(CWDChem) == "Nmass"] <- "MassNitrogen"
names(CWDChem)[names(CWDChem) == "Pmass"] <- "MassPhosphorous"
names(CWDChem)[names(CWDChem) == "Kmass"] <- "MassPotassium"
names(CWDChem)[names(CWDChem) == "Camass"] <- "MassCalcium"
names(CWDChem)[names(CWDChem) == "Mgmass"] <- "MassMagnesium"
names(CWDChem)[names(CWDChem) == "Mnmass"] <- "MassManganese"

#Making sure decay class (1-5) and site ID (1-3) are considered factors
CWDChem$DecayClass<-as.factor(CWDChem$DecayClass)
CWDChem$Site<-as.factor(CWDChem$Site)

#Merging upper and lower CWD diameters
CWDChem$Diameter<-(CWDChem$LowerDiameter+CWDChem$UpperDiameter)/2

######################################

#Diameter
histogram(CWDChem$Diamter)
diametertest<-lm(Diameter~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(diametertest$residuals)
summary(diametertest)
anova(diametertest)

tukeyspeciesdc<-aov(Diameter~TreeSpecies, data=CWDChem)

TukeyHSD(tukeyspeciesdc, "TreeSpecies", ordered=TRUE)

tukeyspeciesdc<-aov(Diameter~DecayClass, data=CWDChem)

TukeyHSD(tukeyspeciesdc, "DecayClass", ordered=TRUE)


#Carbon Mineralization
histogram(CWDChem$mgCO2.g.d)
co2test<-lm(mgCO2.g.d~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(co2test$residuals)
summary(co2test)
anova(co2test)

tukeysco2c<-aov(mgCO2.g.d~TreeSpecies, data=CWDChem)

TukeyHSD(tukeysco2c, "TreeSpecies", ordered=TRUE)

tukeysco2region<-aov(mgCO2.g.d~Region, data=CWDChem)

TukeyHSD(tukeysco2region, "Region", ordered=TRUE)


#Moisture Content
histogram(CWDChem$MoistureContent)
moisturetest<-lm(MoistureContent~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(moisturetest$residuals)
summary(moisturetest)
anova(moisturetest)

tukeymoisturedc<-aov(MoistureContent~DecayClass, data=CWDChem)

TukeyHSD(tukeymoisturedc, "DecayClass", ordered=TRUE)

tukeymoistureeco<-aov(MoistureContent~Region, data=CWDChem)

TukeyHSD(tukeymoistureeco, "Region", ordered=TRUE)

tukeymoisturespecies<-aov(MoistureContent~TreeSpecies, data=CWDChem)

TukeyHSD(tukeymoisturespecies, "TreeSpecies", ordered=TRUE)


#Density
histogram(CWDChem$Density)
densitytest<-lm(Density~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(densitytest$residuals)
summary(densitytest)
anova(densitytest)

tukeydensitydc<-aov(Density~DecayClass, data=CWDChem)

TukeyHSD(tukeydensitydc, "DecayClass", ordered=TRUE)

tukeydensitysite<-aov(Density~Site, data=CWDChem)

TukeyHSD(tukeydensitysite, "Site", ordered=TRUE)

tukeydensityspecies<-aov(Density~TreeSpecies, data=CWDChem)

TukeyHSD(tukeydensityspecies, "TreeSpecies", ordered=TRUE)


#Total Carbon
histogram(CWDChem$TotalCarbon)
carbontest<-lm(TotalCarbon~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(carbontest$residuals)
summary(carbontest)
anova(carbontest)

tukeycarbonspecies<-aov(TotalCarbon~TreeSpecies, data=CWDChem)

TukeyHSD(tukeycarbonspecies, "TreeSpecies", ordered=TRUE)

#Total Nitrogen
histogram(CWDChem$TotalNitrogen)
nitrogentest<-lm(TotalNitrogen~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(nitrogentest$residuals)
summary(nitrogentest)
anova(nitrogentest)

tukeynitrogendc<-aov(TotalNitrogen~DecayClass, data=CWDChem)

TukeyHSD(tukeynitrogendc, "DecayClass", ordered=TRUE)

tukeynitrogenspecies<-aov(TotalNitrogen~TreeSpecies, data=CWDChem)

TukeyHSD(tukeynitrogenspecies, "TreeSpecies", ordered=TRUE)


#Phosphorous
phosphoroustest<-lm(Phosphorous~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(phosphoroustest$residuals)
summary(phosphoroustest)
anova(phosphoroustest)

tukeyphosdc<-aov(Phosphorous~DecayClass, data=CWDChem)

TukeyHSD(tukeyphosdc, "DecayClass", ordered=TRUE)


tukeyphossite<-aov(Phosphorous~Site, data=CWDChem)

TukeyHSD(tukeyphossite, "Site", ordered=TRUE)


tukeyphosphorousspecies<-aov(Phosphorous~TreeSpecies, data=CWDChem)

TukeyHSD(tukeyphosphorousspecies, "TreeSpecies", ordered=TRUE)


#Potassium
histogram(CWDChem$Potassium)
potassiumtest<-lm(Potassium~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(potassiumtest$residuals)
summary(potassiumtest)
anova(potassiumtest)

tukeypotdc<-aov(Potassium~DecayClass, data=CWDChem)

TukeyHSD(tukeypotdc, "DecayClass", ordered=TRUE)

tukeypotassiumspecies<-aov(Potassium~TreeSpecies, data=CWDChem)

TukeyHSD(tukeypotassiumspecies, "TreeSpecies", ordered=TRUE)


#Calcium
histogram(CWDChem$Calcium)
calciumtest<-lm(Calcium~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(calciumtest$residuals)
summary(calciumtest)
anova(calciumtest)

tukeycalciumdc<-aov(Calcium~DecayClass, data=CWDChem)
summary(tukeycalciumdc)

TukeyHSD(tukeycalciumdc, "DecayClass", ordered=TRUE)

tukeycalcsite<-aov(Calcium~Site, data=CWDChem)

TukeyHSD(tukeycalcsite, "Site", ordered=TRUE)

tukeycalciumspecies<-aov(Calcium~TreeSpecies, data=CWDChem)

TukeyHSD(tukeycalciumspecies, "TreeSpecies", ordered=TRUE)


#Magnesium
histogram(CWDChem$Magnesium)
magnesiumtest<-lm(Magnesium~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(magnesiumtest$residuals)
summary(magnesiumtest)
anova(magnesiumtest)

tukeymagnesiumspecies<-aov(Magnesium~TreeSpecies, data=CWDChem)

TukeyHSD(tukeymagnesiumspecies, "TreeSpecies", ordered=TRUE)


#Manganese
histogram(CWDChem$Manganese)
manganesetest<-lm(Manganese~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(manganesetest$residuals)
summary(manganesetest)
anova(manganesetest)

tukeymangdc<-aov(Manganese~DecayClass, data=CWDChem)

TukeyHSD(tukeymangdc, "DecayClass", ordered=TRUE)

tukeymanganesespecies<-aov(Manganese~TreeSpecies, data=CWDChem)

TukeyHSD(tukeymanganesespecies, "TreeSpecies", ordered=TRUE)

#C/N ratio
histogram(CWDChem$C.N)
CNtest<-lm(C.N ~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(CNtest$residuals)
summary(CNtest)
anova(CNtest)

tukeyCNdc<-aov(C.N~DecayClass, data=CWDChem)

TukeyHSD(tukeyCNdc, "DecayClass", ordered=TRUE)

tukeyCNspecies<-aov(C.N~TreeSpecies, data=CWDChem)

TukeyHSD(tukeyCNspecies, "TreeSpecies", ordered=TRUE)

#N/P ratio
histogram(CWDChem$N.P)
NPtest<-lm(N.P ~ Region + Site + TreeSpecies + DecayClass + TreeSpecies*DecayClass, data=CWDChem)
histogram(NPtest$residuals)
summary(NPtest)
anova(NPtest)

tukeyNPdc<-aov(N.P~DecayClass, data=CWDChem)

TukeyHSD(tukeyNPdc, "DecayClass", ordered=TRUE)

tukeyNPsite<-aov(N.P~Site, data=CWDChem)

TukeyHSD(tukeyNPsite, "Site", ordered=TRUE)

