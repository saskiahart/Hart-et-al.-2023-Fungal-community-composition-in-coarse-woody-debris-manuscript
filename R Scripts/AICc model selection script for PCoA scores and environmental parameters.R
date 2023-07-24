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
library(Hmisc)
library(MASS)

citation('Hmisc')
citation('MASS')

#Merging CSV file of PCoA scores environmental parameters 
#PCoA scores obtained from script "PCOA + PERMANOVA + Figures of fungal community composition across stand development stage, decay class, and tree species" 
#and environmental parameters obtained from "Creating compact wood chemistry data set for AICc modelling", both available from https://github.com/saskiahart
pcoachem <- read.csv(file="compactenvapril2021.csv", header = TRUE, sep = ",")
scores<-read.csv(file="PCOA2.csv", head=TRUE, sep=",")
big<- merge(x = nmdschem, y = scores, by = "SampleName", all.y=TRUE)

#Rename columns
colnames(big)
names(big)[names(big) == "Density..g.cm3."] <- "Density"
names(big)[names(big) == "Moisture.Content...."] <- "MoistureContent"
names(big)[names(big) == "TC.."] <- "TotalCarbon"
names(big)[names(big) == "TC.."] <- "TotalCarbon"
names(big)[names(big) == "TKN.."] <- "TotalNitrogen"
names(big)[names(big) == "P.mg.kg"] <- "Phosphorous"
names(big)[names(big) == "K.mg.kg"] <- "Potassium"
names(big)[names(big) == "Ca.mg.kg"] <- "Calcium"
names(big)[names(big) == "Mg.mg.kg"] <- "Magnesium"
names(big)[names(big) == "Mn.mg.kg"] <- "Manganese"
names(big)[names(big) == "C.N"] <- "CNratio"
names(big)[names(big) == "N.P"] <- "NPratio"
names(big)[names(big) == "mg.CO2.g.d"] <- "mgCO2.g.d"


big$DecayClass<-as.factor(big$DecayClass)
big$Site<-as.factor(big$Site)


#Assigning unique levels within each factor
big$DecayClass=factor(big$DecayClass, 
                      levels=unique(big$DecayClass))
big$WoodType=factor(big$WoodType, 
                    levels=unique(big$WoodType))
big$SiteID=factor(big$Site, 
                  levels=unique(big$Site))

big$Region=factor(big$Region, 
                  levels=unique(big$Region))

big$TreeSpecies=factor(big$TreeSpecies, 
                  levels=unique(big$TreeSpecies))



# Correltion matrix with all environmental data
cor1 <- rcorr(as.matrix(big[,c(7:21)]), type="pearson")
cor1p <- as.matrix(cor1$P)
cor1r <- as.matrix(cor1$r)

cor2 <- rcorr(as.matrix(big[,c(7:28)]), type="spearman")
cor2r<- as.matrix(cor2$r)
write.csv(cor1r, file="pcoacorrelation.csv")

#AIC models with correlated (>0.7 Pearson's correlation coefficient) factors removed
AICmod_nmds1<-lm(Axis1 ~MoistureContent+TotalCarbon+Phosphorous+Potassium+Manganese+CNratio+NPratio+Diameter, data=big)
options(na.action = "na.fail")
aic_nmds1<-dredge(AICmod_nmds1,rank="AICc",extra=c("R^2"))
aic_nmds1

aic_nmds1_tab<-as.data.frame(aic_nmds1[aic_nmds1$delta<2,])
nmds1_avg<-(model.avg(aic_nmds1, subset=delta<2))
summary(nmds1_avg)
summary(get.models(aic_nmds1,1)[[1]])


AICmod_nmds2<-lm(Axis2 ~MoistureContent+TotalCarbon+Phosphorous+Potassium+Magnesium+Manganese+CNratio+NPratio, data=big)
options(na.action = "na.fail")
aic_nmds2<-dredge(AICmod_nmds2,rank="AICc",extra=c("R^2"))
aic_nmds2

aic_nmds2_tab<-as.data.frame(aic_nmds2[aic_nmds2$delta<2,])
nmds2_avg<-(model.avg(aic_nmds2, subset=delta<2))
summary(nmds2_avg)
summary(get.models(aic_nmds2,1)[[1]])

