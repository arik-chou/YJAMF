#Reset R's Brain
rm(list=ls())

library(devtools)
#install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
#install.packages("phyloseq")
library(phyloseq)
library(vegan)
library(tidyverse)
library(plyr)
library(nortest)
library(ggpubr)
library(data.table)
library(dplyr)
library(MASS)
#########################################################################
#Bray curtis dissimilarity
#########################################################################

#make the bray curtis dissimilarity matrix
otu_AMF <- read.csv("YJAMF_Merged_47VTonly_norarefaction.csv", row.names = 1)
otu_AMF_t <-t(otu_AMF)


amf_avg_bray_e754_47VT_ <- avgdist(otu_AMF[rowSums(otu_AMF)>753,], 754, iterations=1000, meanfun=median, transf= sqrt, dmethod="bray")

class(amf_avg_bray_e754_47VT_)
amf_avg_bray_e754_47VT_.m=as.matrix(amf_avg_bray_e754_47VT_)
#convert to matrix
#save as csv
write.csv(amf_avg_bray_e754_47VT_.m, "YJAMF_Merged_avg_bray_matrixe754.csv")
amf_avg_bray_e754_47VT_.m = read.csv("YJAMF_Merged_avg_bray_matrixe754.csv", row.names = 1)
#
AMF_data = read.csv("YJAMFe754_AllVT_Richness_and_Metadata.csv", row.names = 1)
#Check the inner_join
AMF_data <- as.data.frame(AMF_data, row.names = AMF_data$id)
AMF_data$Timepoint = as.factor(AMF_data$Timepoint)
AMF_data$Timepoint <-revalue(AMF_data$Timepoint, c("T1"= "Summer", "T2"= "Fall", "T3" = "Winter", "T4" = "Spring"))
AMF_data$Timepoint = factor(AMF_data$Timepoint, levels = c("Summer", "Fall", "Winter", "Spring"))
AMF_data$TimeCol <-revalue(AMF_data$Timepoint, c("Summer" ="gold1" ,"Fall"= "firebrick1", "Winter"="aquamarine1", "Spring" = "green4"))
AMF_data$Type = as.factor(AMF_data$Type)
AMF_data$Treatment = as.factor(AMF_data$Treatment)
AMF_data$Treatment <- relevel(AMF_data$Treatment, ref = "Unburned")

#########################################################################
#Match bray matrix with map
#########################################################################

#discard samples in map but not in matrix

fungimetadata <- read.csv("YJ_Merged_Metadata-nodash.csv", row.names=1)
fungimetadata$id <-row.names(fungimetadata)
map2<-fungimetadata[fungimetadata$id %in% row.names(amf_avg_bray_e754_47VT_.m),]

dim(fungimetadata)
dim(map2)

#discard samples matrix but not in map
amf_avg_bray_e754_47VT_.m2<-amf_avg_bray_e754_47VT_.m[row.names(amf_avg_bray_e754_47VT_.m) %in% map2$id,row.names(amf_avg_bray_e754_47VT_.m) %in% map2$id]

dim(fungimetadata)
dim(map2)
dim(amf_avg_bray_e754_47VT_.m)
dim(amf_avg_bray_e754_47VT_.m2)

#####match up maps to dissimilarity matrices

#match up map to tables with all samples
mapnew<- match(row.names(amf_avg_bray_e754_47VT_.m2), map2$id)
map3<- map2[mapnew,]
mapnew
#check that they match
map3$id == row.names(amf_avg_bray_e754_47VT_.m2)

#check headers
head(amf_avg_bray_e754_47VT_.m2)
head(map3)

table(map3$Type)
#
# Root Soil
# 72   72
table(map3$Timepoint)
#
# T1 T2 T3 T4
# 32 39 35 38

#########################################################################
#Adonis for AMF
#########################################################################
###########Do adonis to determine effect of soil v root, timepoint
map3$Treatment = as.factor(map3$Treatment)
map3$Treatment <- relevel(map3$Treatment, ref = "Unburned")
map3$TreeNum = as.factor(map3$TreeNum)

# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# library(pairwiseAdonis)

adonis_amf_type <- adonis2(amf_avg_bray_e754_47VT_.m2 ~ Type, data=map3, permutations=999)
adonis_amf_type
# Type R2 = 0.097, p<0.001

adonis_amftime <- adonis2(amf_avg_bray_e754_47VT_.m2 ~ Timepoint, data=map3, permutations=999)
adonis_amftime
# Timepoint R2 = 0.039, p<0.008
adonis_amftimetype <- adonis2(amf_avg_bray_e754_47VT_.m2 ~ Timepoint*Type, data=map3, permutations=999)
adonis_amftimetype

adonis_amf_treat <- adonis2(amf_avg_bray_e754_47VT_.m2 ~ Treatment, data=map3, permutations=999)
adonis_amf_treat



# Beta dispersion
soilandroot = map3$Type
SampleType_Betadisper <- betadisper(amf_avg_bray_e754_47VT_, soilandroot)
SampleType_Betadisper

TukeyHSD(SampleType_Betadisper)

anova(SampleType_Betadisper)

plot(SampleType_Betadisper)
boxplot(SampleType_Betadisper)

# Time beta dispersion
seasons = map3$Timepoint
SampleTime_Betadisper <- betadisper(amf_avg_bray_e754_47VT_, seasons)
SampleTime_Betadisper

anova(SampleTime_Betadisper)

TukeyHSD(SampleTime_Betadisper)

plot(SampleTime_Betadisper)
boxplot(SampleTime_Betadisper)

#########################################################################
#NMDS for AMF
#########################################################################
#give colors to soil and root
map3$TypeCol <- ifelse(map3$Type=="Root","red","blue")


#do the NMDS metaMDS on the distance matrix
# k=2 results in too high stress
plotnmds_amf <- metaMDS(amf_avg_bray_e754_47VT_.m2, k=3, trymax=100)
plotnmds_amf
# k=3 stress = 0.181
#make the stressplot
stressplot(plotnmds_amf,amf_avg_bray_e754_47VT_.m2)

#make dataframe of NMDS scores
scores_amf <- as.data.frame(scores(plotnmds_amf ))
dim(scores_amf)

names(map3)
# Takes the values of column 0 and adds them to a new column called "id" which will allow for the inner_join()
scores_amf$id <- row.names(scores_amf)

#do a inner join to join by :id" to remove unneeded samples
betadiversitydata<- inner_join(scores_amf, map3, by="id")

mean_se_data <- betadiversitydata %>%
  group_by(Timepoint) %>%
  summarise(
    mean_NMDS1 = mean(NMDS1),
    mean_NMDS2 = mean(NMDS2),
    se_NMDS1 = sd(NMDS1) / sqrt(n()),
    se_NMDS2 = sd(NMDS2) / sqrt(n())
  )


ggplot(mean_se_data, aes(x = Timepoint)) +
  geom_bar(aes(y = mean_NMDS1), stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar(aes(y = mean_NMDS1, ymin = mean_NMDS1 - se_NMDS1, ymax = mean_NMDS1 + se_NMDS1),
                width = 0.4, color = "black") +
  geom_bar(aes(y = mean_NMDS2), stat = "identity", fill = "lightgreen", alpha = 0.7) +
  geom_errorbar(aes(y = mean_NMDS2, ymin = mean_NMDS2 - se_NMDS2, ymax = mean_NMDS2 + se_NMDS2),
                width = 0.4, color = "black") +
  labs(title = "NMDS Metrics by Timepoint",
       x = "Timepoint",
       y = "Mean NMDS Value") +
  theme_minimal()

dim(betadiversitydata)
betadiversitydata
betadiversitydata$Timepoint <-revalue(betadiversitydata$Timepoint, c("T1"= "Summer", "T2"= "Fall", "T3" = "Winter","T4" = "Spring"))
str(betadiversitydata)
betadiversitydata$Timepoint <-fct_relevel(betadiversitydata$Timepoint, "Summer", "Fall", "Winter", "Spring")

library(grDevices)
betadiversitydata$Type = as.factor(betadiversitydata$Type)
betadiversitydata$Treatment = as.factor(betadiversitydata$Treatment)
betadiversitydata$Timepoint = as.factor(betadiversitydata$Timepoint)

#######
#NMDS plot
NMDS1_treenum<- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2, col=TreeNum, group = TreeNum)) +
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Fungi Bray-Curtis Dissimilarity")+
  # scale_colour_manual(values = )+    #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=14,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #remove major and minor grid panels#+
  stat_ellipse();NMDS1_treenum  #include ellipse around the group


pdf("YJAMF_Merged_597_Timepoint_NMDS_treenum.pdf", width=8,height=8)
NMDS1_treenum
dev.off()


NMDS2_timepoint <- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2, shape = Type, group=Timepoint, col=Timepoint)) +
  geom_point(size=2) +
  theme_bw() +
  labs(title = "AMF VT Bray-Curtis Dissimilarity")+
  scale_color_manual(values=c("darkorange3","firebrick1", "cadetblue3","green4"))+
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=14,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #remove major and minor grid panels#+
  stat_ellipse();NMDS2_timepoint  #include ellipse around the group

#Combined root and soil dissimilarity across seasons
NMDS2_timepoint

NMDS2_timepoint + facet_wrap(~Type)


pdf("Figures/YJAMF_Merged_e754_47VT__TimepointshapedbyType_NMDS.pdf", width=10,height=8)
NMDS2_timepoint + facet_wrap(~Type)
dev.off()

pdf("Figures/YJAMF_Merged_e754_47VT__Timepoint_NMDS.pdf", width=7,height=7)
NMDS2_timepoint
dev.off()




NMDS3_treat <- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2,shape = Type, group=Treatment, col=Treatment)) +
  geom_point(size=2) +
  theme_bw() +
  labs(title = "AMF in Joshua trees - Treatment\nBray-Curtis Dissimilarity")+
  scale_colour_manual(values = c("red","blue"))+    #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=14,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #remove major and minor grid panels#+
  stat_ellipse();NMDS3_treat  #include ellipse around the group

NMDS3_treat + facet_wrap(~Timepoint)

pdf("YJAMF_Merged597_Treatment_NMDS.pdf", width=8,height=8)
NMDS3_treat
dev.off()


NMDS4_type <- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2, col=Type, shape = Type)) +
  geom_point(size=2) +
  theme_bw() +
  labs(title = "AMF VT Bray-Curtis Dissimilarity")+
  scale_colour_manual(values =c("green3","purple2") )+    #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=14,angle=0, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #remove major and minor grid panels#+
  stat_ellipse()  #include ellipse around the group

NMDS4_type + facet_wrap(~Timepoint)
pdf("YJAMF_Merged_e754_47VT_NMDSType.pdf", width=8,height=8)
NMDS4_type
dev.off()
pdf("YJAMF_Merged_e754_47VT_NMDSTypexTime.pdf", width=10,height=8)
NMDS4_type+ facet_wrap(~Timepoint)
dev.off()



## REDO with each timepoint
amf_avg_bray_e754.m2
t1 = "T1"
t2 = "T2"
t3 = "T3"
t4 = "T4"

summer_row <- grep(t1, colnames(amf_avg_bray_e754.m2))
summer_col <- grep(t1, row.names(amf_avg_bray_e754.m2))
summer_matrix <- amf_avg_bray_e754.m2[summer_row, summer_col]

fall_row <- grep(t2, colnames(amf_avg_bray_e754.m2))
fall_col <- grep(t2, row.names(amf_avg_bray_e754.m2))
fall_matrix <- amf_avg_bray_e754.m2[fall_row, fall_col]

winter_row <- grep(t3, colnames(amf_avg_bray_e754.m2))
winter_col <- grep(t3, row.names(amf_avg_bray_e754.m2))
winter_matrix <- amf_avg_bray_e754.m2[winter_row, winter_col]

spring_row <- grep(t4, colnames(amf_avg_bray_e754.m2))
spring_col <- grep(t4, row.names(amf_avg_bray_e754.m2))
spring_matrix <- amf_avg_bray_e754.m2[spring_row, spring_col]

#convert to matrix
summer_matrix <- as.matrix(summer_matrix)
#save as csv
#write.csv(amf_avg_bray_352.m_r, "amf_avg_bray_352.m570_r.csv")


# fungimetadata <- read.csv("YJAMF_metadata_Treenum.csv", row.names=1)
fungimetadata$id <-row.names(fungimetadata)
map2summer<-fungimetadata[fungimetadata$id %in% row.names(summer_matrix),]


dim(fungimetadata)
dim(map2summer)
dim(summer_matrix)

#discard samples matrix but not in map
summer_matrix2<-summer_matrix[row.names(summer_matrix) %in% map2summer$id,row.names(summer_matrix) %in% map2summer$id]

dim(fungimetadata)
dim(map2r)
dim(summer_matrix2)

#####match up maps to dissimilarity matrices

#match up map to tables with all samples: EMF_jac_noDSE, EMF_bray_noDSE, EMF_jac, EMF_bray
mapnew_summer<- match(row.names(summer_matrix2), map2summer$id)
map3summer<- map2summer[mapnew_summer,]

#check that they match
map3summer$id == row.names(summer_matrix2)
# ok good

#check headers
head(summer_matrix2)
head(map3summer)

table(map3summer$Type)
# Root Soil
# 15   15

#########################################################################
#Adonis for summer root v soil
adonis_summer <- adonis2(summer_matrix2 ~ Type, data=map3summer, permutations=999)
adonis_summer
# Df SumOfSqs      R2      F Pr(>F)
# Type      1   0.5743 0.07578 2.2957  0.036 *


########### Fall ###############

#convert to matrix
fall_matrix <- as.matrix(fall_matrix)
#save as csv
#write.csv(amf_avg_bray_352.m_r, "amf_avg_bray_352.m570_r.csv")


# fungimetadata <- read.csv("YJAMF_metadata_Treenum.csv", row.names=1)
fungimetadata$id <-row.names(fungimetadata)
map2fall<-fungimetadata[fungimetadata$id %in% row.names(fall_matrix),]


dim(fungimetadata)
dim(map2fall)
dim(fall_matrix)

#discard samples matrix but not in map
fall_matrix2<-fall_matrix[row.names(fall_matrix) %in% map2fall$id,row.names(fall_matrix) %in% map2fall$id]

dim(fungimetadata)
dim(map2r)
dim(fall_matrix2)

#####match up maps to dissimilarity matrices

#match up map to tables with all samples: EMF_jac_noDSE, EMF_bray_noDSE, EMF_jac, EMF_bray
mapnew_fall<- match(row.names(fall_matrix2), map2fall$id)
map3fall<- map2fall[mapnew_fall,]

#check that they match
map3fall$id == row.names(fall_matrix2)
# ok good

#check headers
head(fall_matrix2)
head(map3fall)

table(map3fall$Type)
# Root Soil
# 20   19

#########################################################################
#Adonis for fall root v soil
adonis_fall <- adonis2(fall_matrix2 ~ Type, data=map3fall, permutations=999)
adonis_fall
# Df SumOfSqs      R2      F Pr(>F)
# Type      1   0.9520 0.13066 5.5612  0.001 ***


######### Winter ###################################


#convert to matrix
winter_matrix <- as.matrix(winter_matrix)
#save as csv
#write.csv(amf_avg_bray_352.m_r, "amf_avg_bray_352.m570_r.csv")


# fungimetadata <- read.csv("YJAMF_metadata_Treenum.csv", row.names=1)
fungimetadata$id <-row.names(fungimetadata)
map2winter<-fungimetadata[fungimetadata$id %in% row.names(winter_matrix),]


dim(fungimetadata)
dim(map2winter)
dim(winter_matrix)

#discard samples matrix but not in map
winter_matrix2<-winter_matrix[row.names(winter_matrix) %in% map2winter$id,row.names(winter_matrix) %in% map2winter$id]

dim(fungimetadata)
dim(map2r)
dim(winter_matrix2)

#####match up maps to dissimilarity matrices

#match up map to tables with all samples: EMF_jac_noDSE, EMF_bray_noDSE, EMF_jac, EMF_bray
mapnew_winter<- match(row.names(winter_matrix2), map2winter$id)
map3winter<- map2winter[mapnew_winter,]

#check that they match
map3winter$id == row.names(winter_matrix2)
# ok good

#check headers
head(winter_matrix2)
head(map3winter)

table(map3winter$Type)
# Root Soil
# 18   17

#########################################################################
#Adonis for winter root v soil
adonis_winter <- adonis2(winter_matrix2 ~ Type, data=map3winter, permutations=999)
adonis_winter
# Df SumOfSqs      R2      F Pr(>F)
# Type      1   1.5062 0.19007 7.7445  0.001 ***




#################### SPRING ############################


#convert to matrix
spring_matrix <- as.matrix(spring_matrix)
#save as csv
#write.csv(amf_avg_bray_352.m_r, "amf_avg_bray_352.m570_r.csv")


# fungimetadata <- read.csv("YJAMF_metadata_Treenum.csv", row.names=1)
fungimetadata$id <-row.names(fungimetadata)
map2spring<-fungimetadata[fungimetadata$id %in% row.names(spring_matrix),]


dim(fungimetadata)
dim(map2spring)
dim(spring_matrix)

#discard samples matrix but not in map
spring_matrix2<-spring_matrix[row.names(spring_matrix) %in% map2spring$id,row.names(spring_matrix) %in% map2spring$id]

dim(fungimetadata)
dim(spring_matrix2)

#####match up maps to dissimilarity matrices

#match up map to tables with all samples: EMF_jac_noDSE, EMF_bray_noDSE, EMF_jac, EMF_bray
mapnew_spring<- match(row.names(spring_matrix2), map2spring$id)
map3spring<- map2spring[mapnew_spring,]

#check that they match
map3spring$id == row.names(spring_matrix2)
# ok good

#check headers
head(spring_matrix2)
head(map3spring)

table(map3spring$Type)
# Root Soil
# 19   19

#########################################################################
#Adonis for spring root v soil
adonis_spring <- adonis2(spring_matrix2 ~ Type, data=map3spring, permutations=999)
adonis_spring
# Df SumOfSqs      R2      F Pr(>F)
# Type      1   0.7202 0.07841 3.0631  0.007 **




adonis_summer
adonis_fall
adonis_winter
adonis_spring

spring_soilandroot = map3spring$Type
SampleType_springBetadisper <- betadisper(spring_matrix2, spring_soilandroot)
SampleType_springBetadisper

TukeyHSD(SampleType_Betadisper)

anova(SampleType_Betadisper)

plot(SampleType_Betadisper)
boxplot(SampleType_Betadisper)
