#Reset R's Brain
rm(list=ls())

# Load packages to test different models
library(plyr)
library(devtools)
library(tidyverse)
library(EcolUtils)
library(SPECIES)
library(vegan)
library(phyloseq)
library(BiodiversityR)
library(scales)
library(ggpubr)
library(nlme)
library(multcomp)
library(stringr)
library(dplyr)
library(car)
library(cowplot)
library(MASS)
library(stats)
library(lmerTest)
library(emmeans)
library(stargazer)

########alpha diversity analysis
asv_AMF_nonc = read.csv("YJAMF_andDFAMF_Merged_23097_collapsed_VT.csv", row.names = 1)



#Remove taxonomy
asvdf_AMF_nonc <- as.data.frame(asv_AMF_nonc)
row.names(asvdf_AMF_nonc) <- asvdf_AMF_nonc$taxonomy
#Remove non-VT species, and taxa column

# Remove originals of repeat samples
# ST1.1, ST1.7, ST2.8
# Also remove RT4.7, which got cross-contaminated with RT4.8 during DNA extraction
repeats=which(colnames(asvdf_AMF_nonc) %in%c("ST1.1", "ST1.7", "ST2.8", "RT4.7"))
repeats
# taxonomy = asvdf_AMF_nonc$Species
YJAMF_nonc_table = as.data.frame(asvdf_AMF_nonc[,-repeats])

AMF_VT_table <- YJAMF_nonc_table[-c(60),c(2:174)]
row.names(AMF_VT_table) = YJAMF_nonc_table[-c(60),]$Species
colnames(AMF_VT_table)

# Remove samples from the dome fire soil project
AMF_VT_table = AMF_VT_table[,-c(90:112)]




#transform table first
asv_AMF_t <- t(AMF_VT_table)
dim(asv_AMF_t)
VTnames = YJAMF_nonc_table[-c(60),]$Species
colnames(asv_AMF_t) = VTnames

# Remove VT from dome fire project
t_zeroes <- which(colSums(asv_AMF_t)==0)
t_zeroes
# s__VTX00085 s__VTX00135 s__VTX00222 s__VTX00244 s__VTX00279 s__VTX00352 s__VTX00377 s__VTX00412
#           9          18          29          31          34          40          47          50

asv_AMF_t_nozeroes <- asv_AMF_t[ ,-t_zeroes]
dim(asv_AMF_t_nozeroes)

getrowsums <- function(transdataframe){
  EMF_summary <- rowSums(transdataframe)
  EMF_summarysorted <- EMF_summary [order(EMF_summary, decreasing = TRUE)]
  boxplot(EMF_summarysorted)
  print(EMF_summarysorted)
  quantile10 <- quantile(EMF_summarysorted, 0.10)
  print(quantile10)
  quantile15 <- quantile(EMF_summarysorted, 0.15)
  print(quantile15)
}
# Order samples by row sum and plot median
getrowsums(asv_AMF_t_nozeroes)
rowSums(asv_AMF_t_nozeroes)
# Exclude the NA column since those aren't VT assigned ASVs
asv_AMF_t_nozeroes[,-c(48)]


colnames(asv_AMF_t_nozeroes)
write.csv(asv_AMF_t_nozeroes[,-c(48)],"YJAMF_Merged_47VTonly_norarefaction.csv")
asv_AMF_t_nozeroes = read.csv("YJAMF_Merged_47VTonly_norarefaction.csv", row.names = 1)
rootcurve = rarecurve(asv_AMF_t_nozeroes[c(1:89),])
soilcurve = rarecurve(asv_AMF_t_nozeroes[c(90:150),])

# Plot rarefaction curve
rarecurve(asv_AMF_t_nozeroes, label = NULL,cex.title = 0.5, cex = 0.5,
          xlab = "Sequence Depth", ylab = "Virtual Taxa",
          main = "AMF Rarefaction Curve")
abline(v = 754, col = "black", lwd = 2, lty = 2)

?specaccum
row.names(asv_AMF_t_nozeroes)
# Plot species accumulation curves
rootaccum = specaccum(asv_AMF_t_nozeroes[c(1:89),], method = "exact" )
soilaccum = specaccum(asv_AMF_t_nozeroes[c(90:150),], method = "exact" )
# plot(rootaccum, col = "green4",ylim=(c(0,70)),xlim=(c(0,100)))
# plot(soilaccum, col = "red4",ylim=(c(0,70)),xlim=(c(0,100)))

plot(rootaccum, col = "green3", ylim = c(0, 70), xlim = c(0, 90), lwd = 2, main = "VT Accumulation in Samples", xlab = "Number of Samples", ylab = "Virtual Taxa")
lines(soilaccum, col = "purple", lwd = 2)
abline(v = 10, col = "black", lwd = 2, lty = 2)
legend("bottomright", legend = c("Root", "Soil"), col = c("green3", "purple"), lty = 1, lwd = 2)

plot(rootcurve, col = "green3", ylim = c(0, 70), xlim = c(0, 90), lwd = 2, main = "VT Accumulation in Samples", xlab = "Number of Samples", ylab = "Virtual Taxa")
lines(soilcurve, col = "purple", lwd = 2)
abline(v = 10, col = "black", lwd = 2, lty = 2)
legend("bottomright", legend = c("Root", "Soil"), col = c("green3", "purple"), lty = 1, lwd = 2)

######## Make a table that's root and soil combined ############
# Filter VT table so that root and soil samples can be combined
# Also remove samples that are only present in 1 sample type
vt_combo = as.data.frame(asv_AMF_t_nozeroes[,-c(48)])
grouping <- sub(".", "", rownames(vt_combo))
grouping_counts <- table(grouping)
valid_groupings <- names(grouping_counts[grouping_counts > 1])
grouping_filtered <- grouping[grouping %in% valid_groupings]
filtered_row_numbers <- which(grouping %in% valid_groupings)
asv_AMF_filtered <- vt_combo[filtered_row_numbers, ]
filtered_data <- cbind(grouping_filtered, asv_AMF_filtered)
original_row_names <- rownames(filtered_data)

# Combine by the row count
combined_table <- aggregate(. ~ grouping_filtered, data = filtered_data, sum)
rownames(combined_table) <- original_row_names[match(combined_table$grouping_filtered, grouping_filtered)]
head(combined_table)

row.names(combined_table) = combined_table$grouping_filtered
ComboVT_Table = combined_table[,-1]
getrowsums(combined_table[,-1])

comboRich = as.data.frame(t(estimateR(combined_table[,-1])))
comboRich$grouping_filtered = row.names(comboRich)
# Extract the timepoint
comboRich$Timepoint <- sub("\\.(\\d+)", "", comboRich$grouping_filtered)
filtered_data$Timepoint <- timepoint_mapping[filtered_data$Timepoint]

# Look at richness
psych::describe(subset(comboRich, Timepoint == "T1")$S.obs)
psych::describe(subset(comboRich, Timepoint == "T2")$S.obs)
psych::describe(subset(comboRich, Timepoint == "T3")$S.obs)
psych::describe(subset(comboRich, Timepoint == "T4")$S.obs)

comboRich


# Test different rarefaction levels
asv_AMF.e4155_mean <- rrarefy.perm(combined_table[,-1], sample =4155, n = 1000, round.out = T)

dim(asv_AMF.e4155_mean) #68 47 Vt
colSums(asv_AMF.e418_mean)
asv_AMF.e4155_mean_no0 <- asv_AMF.e4155_mean[,colSums(asv_AMF.e4155_mean)>0]
dim(asv_AMF.e4155_mean) # 68  47 VT
OTU_AMF_richness <- estimateR(asv_AMF.e4155_mean_no0)
#run function on transformed table to switch rows and columns
OTU_AMF_richness<- t(estimateR(asv_AMF.e4155_mean_no0))

#make a function to get estimates and make TreeNums
estimates_plots_function <- function(datatrans,name){
  #TreeNum S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
  mtext("A",side=3,adj=0)
  #TreeNum S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  mtext("B",side=3,adj=0)
  dev.off()
}


#run function
estimates_plots_function(asv_AMF.e4155_mean_no0,"otu_enon")


#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
otu.H <- diversity(asv_AMF.e4155_mean_no0) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(asv_AMF.e4155_mean_no0, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness_Combo <- data.frame(otu.H,otu.N1,otu.N2)
Combo_richness <- cbind(otu.richness_Combo,OTU_AMF_richness)



### Get Row Sums ordered seqs #####
mean(rowSums(asv_AMF_t_nozeroes))
median(rowSums(asv_AMF_t_nozeroes))

getrowsums(asv_AMF_t_nozeroes[,-48])
# RT4.5    RT4.2    RT4.3    RT4.1   RT4.14    RT2.2   RT4.11   RT4.13   RT4.12   RT4.10    RT4.6
# 39799    37088    29461    27002    21245    20326    20127    19675    19122    18685    16606
# RT4.9   RT4.15    ST2.2   RT1.14   RT3.16    RT1.5    RT2.7    RT2.6   RT1.11    RT2.1    RT2.5
# 15986    15689    15666    15470    15181    14831    14565    14143    13652    13509    13432
# RT4.16    ST4.5   RT1.17   RT4.18    RT4.4   ST4.14    ST2.6   RT4.17    RT2.8   RT1.18   RT2.11
# 13429    12708    12603    11721    11187    10936    10804    10344    10174    10121     9969
# RT2.10   RT1.19   RT1.12   ST1.14   RT4.20   ST2.15   RT2.13   ST2.12   RT3.20    RT2.3    RT4.8
# 9963     9864     9727     9693     9292     9218     9098     9093     9077     8948     8867
# ST4.17   RT2.15   RT3.11   RT3.19   ST2.14   ST4.11   ST3.12   RT3.18   ST2.10   ST1.15   RT1.10
# 8747     8744     8123     8117     7849     7844     7826     7815     7755     7740     7333
# RT2.9   ST1.10    ST2.1   ST4.15    RT3.6    ST3.3   ST3.16    RT3.5    RT1.6   ST2.19   RT3.13
# 7323     7205     7031     7031     6705     6631     6518     6455     6436     6434     6374
# ST3.13    ST2.4   RT2.17   RT2.12   RT1.20    RT3.7   RT3.12    ST3.2   RT1.13    ST3.1    ST1.4
# 6230     6125     6110     6001     5990     5932     5845     5831     5584     5546     5529
# ST4.13    RT3.8    ST3.7    RT2.4   RT3.14   RT3.17    RT1.7    RT1.8   RT2.19   ST2.20   RT3.10
# 5357     5349     5211     5174     4984     4965     4942     4942     4647     4602     4591
# RT2.16   ST3.15 ST2.8.L2   RT4.19    ST3.9   ST4.20   ST2.13   ST1.12   ST1.11   ST3.11    RT1.3
# 4545     4305     4273     4163     4142     4085     4058     3967     3933     3882     3878
# ST4.9   ST4.10    ST4.7   RT1.15   ST4.12   ST4.16    ST2.7    ST1.3   ST2.11   ST2.18   ST3.14
# 3755     3616     3531     3483     3443     3327     3276     3169     3130     3118     3086
# ST2.3   ST1.19 ST1.1.L2   RT2.14   ST3.19   RT2.20    ST2.5   ST3.18    ST3.4    RT3.2   ST3.17
# 2974     2955     2900     2800     2785     2770     2724     2565     2527     2410     2395
# ST4.4   RT2.18    ST4.3   ST2.16    ST4.1    RT3.3    ST4.6    RT3.1    RT3.4    ST2.9   ST1.13
# 2328     2238     2234     2089     2053     1935     1791     1650     1628     1628     1613
# ST1.2    ST4.2   ST4.18    ST4.8    ST1.9    ST3.8   ST3.10 ST1.7.L2    ST1.5    ST1.8   ST1.20
# 1555     1489     1488     1453     1304     1243     1190     1166     1083     1056      889
# ST1.17    ST3.5   ST4.19   ST1.18    ST3.6   ST3.20   ST1.16
# 754      571      552      408      373      342      202
# 10%
# 1484.5
# 15%
# 1841.4

# Write out Dome Fire samples ####
colnames(YJAMF_nonc_table)
DF_VT_table = YJAMF_nonc_table[,c(91:113)]
write.csv(DF_VT_table,"DFSoil_AMF_Table.csv")
asv_AMF_t_nozeroes= asv_AMF_t_nozeroes[,-48]

###### Low  Rarefaction #######
# Try running at a low rarefaction level (ie ST1.16 202 ) to see if many OTUs are lost
# But limit it to a high seq# sample so you can compare them via Mantel test after
# matrix sets need to match
asv_AMF_t_nozeroes
# They are very correlated so just proceed with least rarefaction amount
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>1,]

#### 418 ######
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>417,]
asv_AMF.e418_mean <- rrarefy.perm(asv_AMF.888, sample =418, n = 100, round.out = T)

dim(asv_AMF.e418_mean) #143 47 Vt
colSums(asv_AMF.e418_mean)
asv_AMF.e418_mean_no0 <- asv_AMF.e418_mean[,colSums(asv_AMF.e418_mean)>0]
dim(asv_AMF.e418_mean_no0) # 146  45
# So excluding a few low samples also removes low abund VTs like VTX00067

#### 889 ######
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>888,]

asv_AMF.e889_mean <- rrarefy.perm(asv_AMF.888, sample =889, n = 100, round.out = T)
asv_AMF.e889_mean_no0 <- asv_AMF.e889_mean[,colSums(asv_AMF.e889_mean)>0]
dim(asv_AMF.e889_mean) # 143  47
dim(asv_AMF.e889_mean_no0) # 143  47

#### 1488 ######
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>1487,]
asv_AMF.e1488_mean <- rrarefy.perm(asv_AMF.888, sample =1488, n = 100, round.out = T)
asv_AMF.e1488_mean_no0 <- asv_AMF.e1488_mean[,colSums(asv_AMF.e1488_mean)>0]
dim(asv_AMF.e1488_mean) # 143  47
dim(asv_AMF.e1488_mean_no0) # 135  47


# 754 ####
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>753,]
asv_AMF.e754_mean <- rrarefy.perm(asv_AMF.888, sample =754, n = 100, round.out = T)
asv_AMF.e754_mean_no0 <- asv_AMF.e754_mean[,colSums(asv_AMF.e754_mean)>0]
dim(asv_AMF.e754_mean_no0)
# 144 samples  47 VT
# 754 still OK

# 571 #####
asv_AMF.888<- asv_AMF_t_nozeroes[rowSums(asv_AMF_t_nozeroes)>570,]
asv_AMF.e571_mean <- rrarefy.perm(asv_AMF.888, sample =571, n = 100, round.out = T)
asv_AMF.e571_mean_no0 <- asv_AMF.e571_mean[,colSums(asv_AMF.e571_mean)>0]
dim(asv_AMF.e571_mean_no0) # 145  46 VT
# So down to 571, we do actually lose a VT. Redo analysis at 754

norarDIV = as.data.frame(t(estimateR(asv_AMF_t_nozeroes)))
mean(norarDIV$S.obs)
median(norarDIV$S.obs)

no418DIV = as.data.frame(t(estimateR(asv_AMF.e418_mean_no0)))
mean(no418DIV$S.obs)
median(no418DIV$S.obs)

no754DIV = as.data.frame(t(estimateR(asv_AMF.e754_mean_no0)))
mean(no754DIV$S.obs)
median(no754DIV$S.obs)


no889DIV = as.data.frame(t(estimateR(asv_AMF.e889_mean_no0)))
mean(no889DIV$S.obs)
median(norarDIV$S.obs)

no1488DIV = as.data.frame(t(estimateR(asv_AMF.e1488_mean_no0)))
mean(no1488DIV$S.obs)
median(norarDIV$S.obs)

write.csv(asv_AMF.e754_mean_no0, "YJAMF_Merged_e754_nozeroes_47_VTonly.csv")
asv_AMF.e754_mean_no0 = read.csv("YJAMF_Merged_e754_nozeroes_47_VTonly.csv", row.names = 1)

######## Richness calc ##########
#now calculate richness with BioDiversityR package
#run function on transformed table to switch rows and columns
OTU_AMF_richness<- t(estimateR(asv_AMF.e754_mean_no0))

#make a function to get estimates
estimates_plots_function <- function(datatrans,name){
  #TreeNum S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
  mtext("A",side=3,adj=0)
  #TreeNum S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  mtext("B",side=3,adj=0)
  dev.off()
}
#run function
estimates_plots_function(asv_AMF.e754_mean_no0,"otu_enon")

#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
otu.H <- diversity(asv_AMF.e754_mean_no0) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(asv_AMF.e754_mean_no0, "inv") ## inv. Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness_fungi <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
asv_AMF.e754_richness <- cbind(otu.richness_fungi,OTU_AMF_richness)
#
# Read-write checkpoint ####
write.csv(asv_AMF.e754_richness, "VT_AMF-e754_richness_spponly.csv")
#
asv_AMF.e754_richness <- read.csv("VT_AMF-e754_richness_spponly.csv", row.names = 1)

# Check distribution of species richness
histogram(asv_AMF.e754_richness$S.obs)
shapiro.test(asv_AMF.e754_richness$S.obs)
#Not quite normal, despite its appearance

shapiro.test(log(asv_AMF.e754_richness$S.obs+1))
histogram(log(asv_AMF.e754_richness$S.obs+1))
# Log transform doesn't help. Use non-parametric stats going forward

############## Make some figures #######################################

#read in dataframe with METADATA
fungimetadata <- read.csv("YJ_Merged_Metadata-nodash.csv", row.names=1)

#give columns an id column for the sample ID to be able to join them by
asv_AMF.e754_richness$id <- row.names(asv_AMF.e754_richness)
fungimetadata$id <- row.names(fungimetadata)

#do a left join to join by id
fungidata <- inner_join(asv_AMF.e754_richness, fungimetadata, by="id")
# DF_data = subset(fungidata, Project == "DF")
# write.csv(DF_data,"DF_AMF_Richness_and_Metadata.csv")
AMF_data = subset(fungidata, Project == "YJAMF")
write.csv(AMF_data,"YJAMFe754_AllVT_Richness_and_Metadata.csv")
AMF_data = read.csv("YJAMFe754_AllVT_Richness_and_Metadata.csv", row.names = 1)
#Check the inner_join
AMF_data <- as.data.frame(AMF_data, row.names = AMF_data$id)
AMF_data$Timepoint = as.factor(AMF_data$Timepoint)
AMF_data$Timepoint <-revalue(AMF_data$Timepoint, c("T1"= "Summer", "T2"= "Fall", "T3" = "Winter", "T4" = "Spring"))
AMF_data$Timepoint = factor(AMF_data$Timepoint, levels = c("Summer", "Fall", "Winter", "Spring"))

# AMF_data$TimeCol <-revalue(AMF_data$Timepoint, c("Summer" ="gold1" ,"Fall"= "firebrick1", "Winter"="aquamarine1", "Spring" = "green4"))
AMF_data$Type = as.factor(AMF_data$Type)
AMF_data$Treatment = as.factor(AMF_data$Treatment)
AMF_data$Treatment <- relevel(AMF_data$Treatment, ref = "Unburned")
table(AMF_data$Type)
# Do pairs() on some samples to get an overview of comparisons
# Climatic data not as useful with 4 timepoints; its just 4 numbers instead of 4 seasons

# pairs(AMF_data[,-c(5:24)])
# write as csv
str(AMF_data)
#reorganize levels to be in order you want, regardless of alphabetical order
#look up colors in r: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

venn_colors =c("darkorange3","firebrick1","cadetblue3", "green4")

str(AMF_data)
Fig1_Time=ggplot(AMF_data, aes(x=Timepoint, y=S.obs, col = Timepoint ))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1)+ #  data points to show but make them transparent
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  theme_bw() + #make black and white
  ggtitle("AMF Richness by Season")+
  ylim(c(0,25))+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=45,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_Time    #make x axis text larger and angled
# Quick stat
Fig1_Time
Fig1_Time + facet_wrap(~Type)
pdf("Figures/YJAMF47_e754_Seasons.pdf", w=10, h=7)
Fig1_Time
dev.off()


pdf("Figures/YJAMF_e754_spponly_Merged_VT_Seasons_Typefacet.pdf", w=12, h=7)
Fig1_Time + facet_grid(~Type)
dev.off()



Fig1_timetype =ggplot(AMF_data, aes(x=Type, y=S.obs, col = Timepoint))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  theme_bw() + #make black and white
  ylim(c(0,22))+
  ggtitle("AMF Richness by Sample Type")+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=0,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_timetype

Fig1_type =ggplot(AMF_data, aes(x=Type, y=S.obs, col = Type))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("green2","purple"))+
  theme_bw() + #make black and white
  ylim(c(0,20))+
  # ggtitle("AMF Richness \nby Sample Type")+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    legend.title = element_blank(), #remove legend titles
    legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=0,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_type
Fig1_type
pdf("Figures/YJAMFe754_Merged_SoilandRoot.pdf", w=4,h=7)
Fig1_type
dev.off()

pdf("Figures/YJAMFe754_Merged_SoilandRoot_TimeSidexSide.pdf", w=10,h=7)
Fig1_timetype
dev.off()





Fig1_treatment =ggplot(AMF_data, aes(x=Treatment, y=S.obs, col = Treatment, group = Treatment ))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("blue", "red"))+
  theme_bw() + #make black and white
  # ggtitle("AMF Richness by Treatment")+
  ylim(c(0,20))+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    legend.title = element_blank(), #remove legend titles
    legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=0,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_treatment

Fig1_treatment + facet_grid(~Type)
Fig1_treatment
pdf("Figures/YJAMFe754_Merged_Burn_Thin.pdf", w=4,h=7)
Fig1_treatment
dev.off()

pdf("Figures/YJAMF_Merged47VT_BurnTreatment_TypeFacet.pdf", w=4,h=7)
Fig1_treatment + facet_grid(~Type)
dev.off()


AMF_data$Library = as.factor(AMF_data$Library)
Fig1_Lib  =ggplot(AMF_data, aes(x=Library, y=S.obs,group = Library ))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  # scale_color_manual(values=c("gold1","firebrick1","aquamarine1", "green4"))+
  theme_bw() + #make black and white
  ggtitle("AMF Richness by Season")+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=45,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_Lib +facet_wrap(~Type)

pdf("Figures/YJAMFe754Merged47spponly_LibraryTest.pdf",w=14,h=8)
Fig1_Lib + facet_grid(~Timepoint)
dev.off()

Fig1_Shannon =ggplot(AMF_data, aes(x=Timepoint, y=otu.N1,color = Timepoint ))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  theme_bw() + #make black and white
  ggtitle("AMF Richness by Season")+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=45,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_Shannon
Fig1_Shannon + facet_wrap(~Type)

Fig1_Simpson =ggplot(AMF_data, aes(x=Timepoint, y=otu.N2,color = Timepoint ))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  theme_bw() + #make black and white
  ggtitle("AMF Richness by Season")+
  ylab("AMF Simpson Diversity") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=45,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_Simpson

Fig1_SimpsonType =ggplot(AMF_data, aes(x=Type, y=otu.N2, color = Timepoint))+ #TreeNum qpcr data, x axis is the letter, y axis in the copy number, group by leter for mean and SE, shape by kingdom
  geom_boxplot(size=1, alpha = 0.4)+ #get all data points to show but make them transparent
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  theme_bw() + #make black and white
  ylim(c(0,15))+
  ggtitle("AMF Richness by Season")+
  ylab("AMF Simpson Diversity") +  #change yaxis label
  theme(
    # legend.title = element_blank(), #remove legend titles
    # legend.position = "none",
    text = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(), #remove TreeNum title
    axis.text.y = element_text(size=18, angle=45,hjust=0.5),
    axis.text.x = element_text(size=18));Fig1_SimpsonType

#### Version that's mean + se instead of boxes


Fig_SumRich_Seasons <-ggplot(AMF_data, aes(x=Timepoint, y=S.obs, color = Timepoint))+
  geom_point(size=1, alpha=0.2)+ #get all data points to show but make them transparent
  stat_summary(fun=mean,geom="point", size=4)+ #plot the mean
  stat_summary(fun.data = mean_se,geom = "errorbar", size=0.5, #plot the standard error
               alpha=0.7,position = position_dodge(0.01))+
  theme_bw() + #make black and white
  scale_color_manual(values=c("darkorange3","firebrick1","cadetblue3", "green4"))+
  ggtitle("Mean Richness of Soil and Root AMF VTs")+
  ylim(c(0,20))+
  ylab("AMF VT Richness") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(), #remove legend titles
        text = element_text(size=14),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=14, angle=0,hjust=0.5));Fig_SumRich_Seasons     #make x axis text larger and angled

#### Stats ####
library(rstatix)

std.error <- function(x) sd(x)/sqrt(length(x))
table(AMF_data$Type)
table(AMF_data$S.obs)
AMF_soildata <- subset(AMF_data, Type == "Soil")
AMF_rootdata <- subset(AMF_data, Type == "Root")


mean(AMF_data$S.obs)
std.error(AMF_data$S.obs)

mean(AMF_rootdata$S.obs)
mean(AMF_soildata$S.obs)

std.error(AMF_rootdata$S.obs)
std.error(AMF_soildata$S.obs)

table(AMF_data$Treatment)
AMF_burndata <- subset(AMF_data, Treatment == "Burned")
AMF_ubrndata <- subset(AMF_data, Treatment == "Unburned")

mean(AMF_burndata$S.obs)
mean(AMF_ubrndata$S.obs)

std.error(AMF_burndata$S.obs)

std.error(AMF_ubrndata$S.obs)


AMF_summdata <- subset(AMF_data, Timepoint == "Summer")
AMF_Falldata <- subset(AMF_data, Timepoint == "Fall")
AMF_Winterdata <- subset(AMF_data, Timepoint == "Winter")
AMF_Springdata <- subset(AMF_data, Timepoint == "Spring")

mean(AMF_summdata$S.obs)
std.error(AMF_summdata$S.obs)

mean(AMF_Falldata$S.obs)
std.error(AMF_Falldata$S.obs)

mean(AMF_Winterdata$S.obs)
std.error(AMF_Winterdata$S.obs)

mean(AMF_Springdata$S.obs)
std.error(AMF_Springdata$S.obs)

