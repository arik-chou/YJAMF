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

# #Import ASV Table ##########
asv_AMF <- read.csv("YJ_Merged_ASVs_wtaxa.csv", row.names =1)

# Many ASVs are genus level, not VT-matched
# Add additional BLAST matches based on e-50 and 90% seq similarity + less than 10 nt difference
asv_Genbank_matches = read.csv("YJAMF_nonhit_ASVs_Chaudhary_Parameters.csv")
asv_Genbank_match9 = subset(asv_Genbank_matches, mismatch < 10)

asv_VT_BLAST = as.data.frame(cbind(asv_Genbank_match9$id, asv_Genbank_match9$Species))
colnames(asv_VT_BLAST) = c("id", "Species")

# List of ASVs which have 2 matches
#### Check these ASVs that had equal matches to 2 hits after BLASTn ####


purples = c("f9b0844730df223873129465f909defa",
            "f8ce0b27adf731867dbf1dc5409a9db4",
            "ef2c793c217a64cbd1d497181c978af3",
            "e937d0b10167c48b92b37d034d70ce7e",
            "e4ca6dca62477c2fb5c3d8291c9fa257",
            "daa7365b7bd5ce9b0d1e774b167cc560",
            "da964523be7ae4579a549e9f0a441711",
            "d6a952f449eb0aa234a93d2a5e2c8dbe",
            "d4c8067e86bec553ae52bf59cd5f1ae8",
            "ce0138a57583487f958d4a51f8004668",
            "cb6e899dab1c35d2610b06bbbbb704f9",
            "bfb6167597944586b98628fb80fa358d",
            "b237bc7bb15958209ed79e05a0ab039d",
            "b0d4f7de198424e081d942d2dde5fddd",
            "a9b4f9728ee7b91bcdf9b426c1e41e19",
            "a831eb51c33481307fc53a27fb180b81",
            "a667a7efbb3eae9eca25439cc95858de",
            "9f720f3cf5dd673b09d4c700a5aea576",
            "9e289141bb3ababb090a08c86b756b7f",
            "97c00bc7387cdc517c6699be631fa8df",
            "939bcfe9712939526f9227f2d7fe1eae",
            "92169f440a2cc0185186526c9b0bdf6a",
            "8791f3388b5ce6ea5bb30a1d081eb4b0",
            "7c9d459347eefc31dcad3c53b4a61e74",
            "873aa6b675269f2a991df8716ebaf0a8",
            "7edd49b8b0ba4904438a92be161ae332",
            "76a7fe44a8b17f4654bc0d7fa3eac2af",
            "74daae9a03708ba4b22fd0826f598b3c",
            "64208b836125c05ef4a39baebe4c5360",
            "571c4497dc8d31d0c9c72db2846a2baf",
            "5311b5e4ef2e9c31ddf27ce09d8d48bc",
            "43e1537999215ca4d47a6230f8602083",
            "40e61734658870eabfab0addaac54530",
            "33d9cbbf89a8769d0081cf0056a00c25",
            "323ecec5d4dc7ff588a9e2a348f487f4",
            "2f05d0004c3629c4e10647fa13aa8d1a",
            "2df88257116ddbe003e533817755bc48",
            "29d29125a500c2400718bbd2cf35e775",
            "2684ad867875fdd63ad0c13da521cad2",
            "25b56d3d3055fb80bfd83af54612f3ef",
            "197683ade253a416dcf09e8ead56fbf7",
            "25e1097c58eacdb7013fd8a49580f774",
            "01be3edc9cfbe1065817b48e16fe0145")

pruplecount = asv_AMF[purples,]
purple_notaxa = pruplecount[,c(1:190)]
pruplecount$SUM = rowSums(purple_notaxa)
sum(pruplecount$SUM) # 143873

# A few high abundance ASVs, overall 143k reads from 43 ASVs
# Try to match them
# Check the ones that were actually over 10 mismatches in the BLAST

reds =
  c("07710e333b9be95efc5a92e7a4660758",
    "18b4c198c1e823c979bc6a46919fdbd4",
    "24e6099a63e006cd57b013e5434c95f2",
    "30b7e16d86a3c7c0db196d6f7821d3ad",
    "336c7fbad947b61cab0a95f5c4779d24",
    "4b16022e08a614601f3f4c70b2611d21",
    "566b85060d505bfd7138061c937e05ff",
    "654f1ac3a2bd69f6a66fd80f4d170e7a",
    "6a76b1b4761a51e1cb2b6a93e4aa8c01",
    "6ee1d7deb113f77d91600e4e292f9561",
    "703a59031df9da288358bd9eae0d9fe6",
    "82ee1b29a93cc70219bb20eabb0e41c6",
    "83fbc99482b9f7143e7280f04235d1a1",
    "8bed179a422373373d74e4794c6e02e1",
    "94a428a67cb7182c50361dda546cf7a0",
    "95d2d7513e5ee3cd3dc2ac710a8bc3c3",
    "9a6c77e3ac515acd4cdf09315fc912ae",
    "9f12f16e17f6c07a02479da53371f3ca",
    "a1a4865318b0bc12e292b1ccf7d06a11",
    "a37f39638013fa0ce11c3300ff7e72b8",
    "a90268724b8e16c5187ca533a7e230e4",
    "aaa7c6ae8dfd3a18303dd8a6b3560895",
    "b1bc0e2e6e625fffc99581a04952c400",
    "b6bcac85e7155d4426f2a259c12745b3",
    "bd89e8a087968f8a8f91d9caad3a8d28",
    "c93d7a11d9398ced3a40465d963a0bc6",
    "d0692a66cc4ef1b0a22d1254b2ce5046",
    "d034c5a615774f1b5a982137673f6923",
    "d2a6ad57986236488fcfcd7844a936db",
    "d2bd998b0726cd1b64c063b705b1aa39",
    "d96283dd0a05f79c44904e5caa8e8282",
    "dad6c47d68edc060f885e092c57a4b74",
    "e2c9e4d8b1127a2406c7f84145d9ad98",
    "e2dac4e80e12c9896047fe6e2eec0f78",
    "e698224311a81af8ba1141a7daea89cf",
    "e8ca0956c3f825fd1831caf6b717ceb9",
    "fb4dc690a792df0dee660bec5487c6d8")


redcount = asv_AMF[reds,]
sum(redcount[,-c(191)])
redcount$SUM = rowSums(redcount[,-c(191)])
# 13,140 reads represented by ASVs that matched but had >9 mismatches


lastcolumn <- asv_AMF$taxonomy
id <- row.names(asv_AMF)

#Data frame with the OTU IDs and the taxa associated with it
dd <- data.frame(id,lastcolumn)


otu_taxa1 <- ldply(str_split(string = dd$lastcolumn, pattern=";"), rbind)

# Divide a column using ";"and convert list to data frame
names(otu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_taxa2 <- as.data.frame(lapply(otu_taxa1, gsub, pattern=" ", replacement=""))
dd <- cbind(dd[,1:2 ],otu_taxa2)
rownames(dd) <- id

length(unique(dd$Species)) #43
# Combine blasted, add species where VT wasn't present ###
merged_df <- left_join(dd, asv_VT_BLAST, by = "id")
merged_df$Species <- ifelse(is.na(merged_df$Species.x), merged_df$Species.y, merged_df$Species.x)
dd_plusgenbank = merged_df[, -c(9,10 )]
length(unique(dd_plusgenbank$Species)) #56
# 56 taxa but also includes samples from Dome Fire soil - non-tree samples

# write.csv(dd, "YJAMFandDF_Taxa.csv")
asv_AMF$id = row.names(asv_AMF)
#Aggregate ASVs into virtual taxa
asv_BLAST_AMF = left_join(asv_AMF,dd_plusgenbank, by = "id" )

asv_BLAST_AMF[, -c(191:199)] %>% group_by(Species) %>% summarise_each(funs(sum)) -> asv_collapsed

collapsenames = asv_collapsed$Species





noDNA <- which(colnames(asv_collapsed) %in%c("RT1.NC","RT2.NC","RT2.PCRN","RT3.NC","ST1.NC","ST1.PCRN","ST2.NC",
                                             "RT3PCRN","Undetermined","L2Undetermined", "PCR.NEG1", "PCR.NEG2",
                                             "PCR.POS","RT4_NEG", "SOILT4.NEG", "T3BATCH.NEG"
))
rowSums(asv_collapsed[,-c(1)])
noDNA
#make otu table of no DNA controls
OTU_noDNA <-asv_collapsed[ ,noDNA]
OTU_noDNA
colSums(OTU_noDNA)
sum(OTU_noDNA)

### Collapsing into VT and deleting everything else method #####

#From 430 ASVs to 250 IDs via BLASTn with -10e50 cutoff to 44 Genbank IDs / Virtual Taxa
#Remove controls now
noDNA <- which(colnames(asv_collapsed) %in%c("RT1.NC","RT2.NC","RT2.PCRN","RT3.NC","ST1.NC","ST1.PCRN","ST2.NC",
                                             "RT3PCRN","Undetermined","L2Undetermined", "PCR.NEG1", "PCR.NEG2",
                                             "PCR.POS","RT4_NEG", "SOILT4.NEG", "T3BATCH.NEG"
))
# Version with only neg controls
negDNA <- which(colnames(asv_collapsed) %in%c("RT1.NC","RT2.NC","RT2.PCRN","RT3.NC","ST1.NC","ST1.PCRN","ST2.NC",
                                              "RT3PCRN", "PCR.NEG1", "PCR.NEG2","RT4_NEG", "SOILT4.NEG", "T3BATCH.NEG"
))



rowSums(asv_collapsed[,-c(1)])
sum(noDNA) # 1026
#make otu table of no DNA controls
OTU_noDNA <-asv_collapsed[ ,noDNA]
OTU_noDNA
row.names(OTU_noDNA) = collapsenames

OTU_negDNA <-asv_collapsed[ ,negDNA]
OTU_negDNA
row.names(OTU_negDNA) = collapsenames

colSums(OTU_negDNA)
rowSums(OTU_negDNA)

# > colSums(OTU_noDNA)
# RT1.NC      RT2.NC    RT2.PCRN      RT3.NC     RT3PCRN    PCR.NEG1    PCR.NEG2     RT4_NEG
# 44           2           2           2          37           2           3           2
# SOILT4.NEG T3BATCH.NEG
# 16           2

#pretty minimal spill-over

#Remove NCs and undetermined
asv_AMF_nonc <-asv_collapsed[ ,-noDNA]
colnames(asv_AMF_nonc)
row.names(asv_AMF_nonc)

# Check ASV count in ASV Table minus dome fire soil samples
ASV_jt=asv_BLAST_AMF[,-c(105:128)]
noDNA <- which(colnames(ASV_jt) %in%c("RT1.NC","RT2.NC","RT2.PCRN","RT3.NC","ST1.NC","ST1.PCRN","ST2.NC",
                                      "RT3PCRN","Undetermined","L2Undetermined", "PCR.NEG1", "PCR.NEG2",
                                      "PCR.POS","RT4_NEG", "SOILT4.NEG", "T3BATCH.NEG"
))

ASV_Treesonly=ASV_jt[,-noDNA]
names(ASV_Treesonly)
ASV_Treesonly


ASV_Treesonly_nozero <- ASV_Treesonly[rowSums(ASV_Treesonly[,c(1:154)]) != 0, ]
repeats=which(colnames(ASV_Treesonly_nozero) %in%c("ST1.1", "ST1.7", "ST2.8", "RT4.7"))
repeats

ASV_Treesonly_nodupes = as.data.frame(ASV_Treesonly_nozero[,-repeats])
colnames(ASV_Treesonly_nodupes)
ASV_Treesonly_count <- ASV_Treesonly_nodupes[rowSums(ASV_Treesonly_nodupes[,c(1:150)]) != 0, ]

dim(ASV_Treesonly) # 1706 ASV starting point with Dome Fire, undetermined, controls etc
dim(ASV_Treesonly_nozero) # 1282 ASV including repeat samples
dim(ASV_Treesonly_count) # 1270 ASV after removal of extra samples

# Check how many are NA on species
countofNA =which(is.na(ASV_Treesonly_count$Species))
JTSeq_NAsubset = ASV_Treesonly_count[countofNA,]
JTSeq_NAsubset$RowSum = rowSums(JTSeq_NAsubset[,c(1:150)])
JTSeq_NAsubset$RowSum

sum(JTSeq_NAsubset$RowSum) # Total non-assigned ASV reads: 174,239
# Many matched to non-AMF things after BLASTing so excluded from analysis

# Impromptu richness count of ASVs

tASV_Treesonly_count =as.data.frame(t(ASV_Treesonly_count[,c(1:150)]))
logical_matrix <- tASV_Treesonly_count > 1

# Sum the TRUE values for each row to get the richness
richness <- rowSums(logical_matrix)

# Add the richness as a new column to the data frame
tASV_Treesonly_count$Richness <- richness
mean(tASV_Treesonly_count$Richness)
median(tASV_Treesonly_count$Richness)
# Add the richness as a new column to the data frame
dim(tASV_Treesonly_count)

tASV_Treesonly_count$RowSum = rowSums(tASV_Treesonly_count[,c(1:1270)])
mean(tASV_Treesonly_count$RowSum) # 8427.533
median(tASV_Treesonly_count$RowSum) # 6817

sum(tASV_Treesonly_count$RowSum) # 1264133 -Total reads contributing to joshua tree project

###### End of sums, remove non-VTs ######
write.csv(asv_AMF_nonc, "YJAMF_andDFAMF_Merged_23097_collapsed_VT.csv")
asv_AMF_nonc = read.csv("YJAMF_andDFAMF_Merged_23097_collapsed_VT.csv", row.names = 1)
