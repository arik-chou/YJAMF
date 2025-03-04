#Reset R's Brain
rm(list=ls())

#Packages to load!
library(devtools)
library(tidyverse)
library("VennDiagram")
library(ggpubr)
library(plyr)

# Rarefied virtual taxa
YJAMF_Merged_47VTonly_norarefaction.csv
asv_AMF.e754_nozeroes = read.csv("YJAMF_Merged_e754_nozeroes_47_VTonly.csv", row.names = 1)
fungimetadata <- read.csv("YJ_Merged_Metadata-nodash.csv", row.names=1)
fungimetadata$id <- row.names(fungimetadata)
asv_AMF.e754_nozeroes <-as.data.frame(asv_AMF.e754_nozeroes)
asv_AMF.e754_nozeroes$id <-row.names(asv_AMF.e754_nozeroes)
meta_AMF <-inner_join(asv_AMF.e754_nozeroes, fungimetadata, by="id")
row.names(meta_AMF) <- meta_AMF$id
meta_AMF = subset(meta_AMF, Project == "YJAMF")

meta_AMF$Timepoint <-revalue(meta_AMF$Timepoint, c("T1"= "Summer", "T2"= "Fall", "T3" = "Winter", "T4" = "Spring"))
table(meta_AMF$Timepoint)
# Fall Spring Summer Winter
# 39     38     31     35
table(meta_AMF$Type)
# Root Soil
# 72   71
Meta_Root <- subset(meta_AMF, Type == "Root")
Meta_Soil <- subset(meta_AMF, Type == "Soil")

table(Meta_Root$Timepoint)
# Fall Spring Summer Winter
# 20     19     15     18


table(Meta_Soil$Timepoint)
# Fall Spring Summer Winter
# 19     19     16     17


# Check Venn Diagrams with matching sub-sampled N
# Subsample to match the N ####
meta_Rsub <- Meta_Root %>%
  group_by(Timepoint) %>%
  sample_n(size = 15)

meta_Ssub <- Meta_Soil %>%
  group_by(Timepoint) %>%
  sample_n(size = 15)
# Replace when comparing
################# Soil and Root Time Venn #######

# subset Timepoints for Root #####
Summer_spp <- subset(meta_AMF, Timepoint == "Summer")

#Remove extra metadata
Summer_OTU <- Summer_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozero<- Summer_OTU[, colSums(Summer_OTU != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTU <- t(Summer_OTUs_nozero)
Summ_OTU_list <- row.names(Summ_tOTU)


# subset
Fall_spp <- subset(meta_AMF, Timepoint == "Fall")

#Remove extra metadata
Fall_OTU <- Fall_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozero<- Fall_OTU[, colSums(Fall_OTU != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTU <- t(Fall_OTUs_nozero)
Fall_OTU_list<- row.names(Fall_tOTU)

# subset
Winter_spp <- subset(meta_AMF, Timepoint == "Winter")

#Remove extra metadata
Winter_OTU <- Winter_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Winter_OTUs_nozero<- Winter_OTU[, colSums(Winter_OTU != 0) > 0]

#Transform back and get list of OTUs
Winter_tOTU <- t(Winter_OTUs_nozero)
Winter_OTU_list <- row.names(Winter_tOTU)


# subset
Spring_spp <- subset(meta_AMF, Timepoint == "Spring")

#Remove extra metadata
Spring_OTU <- Spring_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Spring_OTUs_nozero<- Spring_OTU[, colSums(Spring_OTU != 0) > 0]

#Transform back and get list of OTUs
Spring_tOTU <- t(Spring_OTUs_nozero)
Spring_OTU_list <- row.names(Spring_tOTU)

## All Seasons, Root and Soil ##

venn_colors <- c("gold","aquamarine","firebrick","green")
### Venn Diagrams just output straight to the folder you're in
venn.diagram(
  x = list(Summ_OTU_list,Winter_OTU_list,Fall_OTU_list,Spring_OTU_list),
  category.names = c("Summer","Winter","Fall","Spring"),
  filename = 'YJAMF_Genbank_AllSamples_asv_AMF.e754_47nozeroes_venn.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1080 ,
  width = 1080 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  main = "AMF Richness Overlap",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "serif"
)

################################################
################################################
# Version with subsampling
# subset Timepoints for Root #####
Summer_sppr <- subset(meta_Rsub, Timepoint == "Summer")

#Remove extra metadata
Summer_OTUr <- Summer_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozeror<- Summer_OTUr[, colSums(Summer_OTUr != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTUr <- t(Summer_OTUs_nozeror)
Summ_OTU_listr <- row.names(Summ_tOTUr)


# subset
Fall_sppr <- subset(meta_Rsub, Timepoint == "Fall")

#Remove extra metadata
Fall_OTUr <- Fall_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozeror<- Fall_OTUr[, colSums(Fall_OTUr != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTUr <- t(Fall_OTUs_nozeror)
Fall_OTU_listr <- row.names(Fall_tOTUr)

# subset
Winter_sppr <- subset(meta_Rsub, Timepoint == "Winter")

#Remove extra metadata
Winter_OTUr <- Winter_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Winter_OTUs_nozeror<- Winter_OTUr[, colSums(Winter_OTUr != 0) > 0]

#Transform back and get list of OTUs
Winter_tOTUr <- t(Winter_OTUs_nozeror)
Winter_OTU_listr <- row.names(Winter_tOTUr)


# subset
Spring_sppr <- subset(meta_Rsub, Timepoint == "Spring")

#Remove extra metadata
Spring_OTUr <- Spring_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Spring_OTUs_nozeror<- Spring_OTUr[, colSums(Spring_OTUr != 0) > 0]

#Transform back and get list of OTUs
Spring_tOTUr <- t(Spring_OTUs_nozeror)
Spring_OTU_listr <- row.names(Spring_tOTUr)

## All Seasons, Root subsampled ##

venn_colors <- c("gold3","aquamarine4","firebrick4","green4")
### Venn Diagrams just output straight to the folder you're in
venn.diagram(
  x = list(Summ_OTU_listr,Winter_OTU_listr,Fall_OTU_listr,Spring_OTU_listr),
  category.names = c("Summer","Winter","Fall","Spring"),
  filename = 'YJAMF_Genbank_Root_venn_subset_n_of_15.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1080 ,
  width = 1080 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  main = "Root AMF Richness Overlap",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# subset Timepoints for Soil with subsampling #####
Summer_spp <- subset(meta_Ssub, Timepoint == "Summer")

#Remove extra metadata
Summer_OTU <- Summer_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozero<- Summer_OTU[, colSums(Summer_OTU != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTU <- t(Summer_OTUs_nozero)
Summ_OTU_list <- row.names(Summ_tOTU)


# subset
Fall_spp <- subset(meta_Ssub, Timepoint == "Fall")

#Remove extra metadata
Fall_OTU <- Fall_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozero<- Fall_OTU[, colSums(Fall_OTU != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTU <- t(Fall_OTUs_nozero)
Fall_OTU_list <- row.names(Fall_tOTU)

# subset
Winter_spp <- subset(meta_Ssub, Timepoint == "Winter")

#Remove extra metadata
Winter_OTU <- Winter_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Winter_OTUs_nozero<- Winter_OTU[, colSums(Winter_OTU != 0) > 0]

#Transform back and get list of OTUs
Winter_tOTU <- t(Winter_OTUs_nozero)
Winter_OTU_list <- row.names(Winter_tOTU)


# subset
Spring_spp <- subset(meta_Ssub, Timepoint == "Spring")

#Remove extra metadata
Spring_OTU <- Spring_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Spring_OTUs_nozero<- Spring_OTU[, colSums(Spring_OTU != 0) > 0]

#Transform back and get list of OTUs
Spring_tOTU <- t(Spring_OTUs_nozero)
Spring_OTU_list <- row.names(Spring_tOTU)


## All Seasons, Soil ##

venn_colors <- c("gold1","aquamarine1","firebrick1","green2")
### Venn Diagrams just output straight to the folder you're in
venn.diagram(
  x = list(Summ_OTU_list,Winter_OTU_list,Fall_OTU_list,Spring_OTU_list),
  category.names = c("Summer","Winter","Fall","Spring"),
  filename = 'YJAMF_Genbank_Soil_venn_subset_n_of_15.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1080 ,
  width = 1080 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  main = "Soil AMF Richness Overlap",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)









#List unique taxa if you want to dig deeper ####
#ASVs listed are in the first list (ie in x, not in y)
setdiff(Summ_OTU_list, Fall_OTU_list)
setdiff(Summ_OTU_list, Winter_OTU_list)
setdiff(Fall_OTU_list, Winter_OTU_list)
setdiff(Fall_OTU_list, Spring_OTU_list)
setdiff(Summ_OTU_list, Spring_OTU_list)
setdiff(Winter_OTU_list, Spring_OTU_list)

#Prints OTU list
sink("Summer_OTU_gen.txt")
print(Summ_OTU_list)
sink()

#List unique OTUs in Fall

sink("Fall_OTU_gen.txt")
print(Fall_OTU_list)
sink()

sink("Winter_OTU_gen.txt")
print(Winter_OTU_list)
sink()

sink("Spring_OTU_gen.txt")
print(Spring_OTU_list)
sink()

## Function to just look at the diagrams when it saves
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

##### Root V Soil Total #######
root_meta_AMF = subset(meta_AMF, Type == "Root")

soil_meta_AMF <-subset(meta_AMF, Type == "Soil")
#Root
#Remove extra metadata
root_OTU <- root_meta_AMF[, c(1:47)]

#Remove OTUs that aren't present in any samples
root_OTU_nozero<- root_OTU[, colSums(root_OTU != 0) > 0]

#Transform back and get list of OTUs
root_OTU_nozerot <- t(root_OTU_nozero)
Root_List <- row.names(root_OTU_nozerot)


#soil
#Remove extra metadata
soil_OTU <- soil_meta_AMF[, c(1:47)]

#Remove OTUs that aren't present in any samples
soil_OTU_nozero<- soil_OTU[, colSums(soil_OTU != 0) > 0]

#Transform back and get list of OTUs
soil_OTU_nozerot <- t(soil_OTU_nozero)
soil_List <- row.names(soil_OTU_nozerot)


table(meta_AMF$Type)


venn.diagram(
  x = list(soil_List,Root_List),
  category.names = c("Soil\nn=71","Root\nn=72"),
  filename = 'YJAMF__Merged_47VT_asv_AMF.e754_nozeroes_RootvSoil.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1000 ,
  width = 1000 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("purple2","green3"),


  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)
#List unique OTUs in PC if you want to dig deeper
setdiff(Root_List, soil_List)

setdiff(soil_List,Root_List)




### Comparisons of all samples, not subsampled for equal # #######


# subset Timepoints for Soil without subsampling #####
Summer_spp <- subset(Meta_Soil, Timepoint == "Summer")

#Remove extra metadata
Summer_OTU <- Summer_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozero<- Summer_OTU[, colSums(Summer_OTU != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTU <- t(Summer_OTUs_nozero)
Summ_OTU_list <- row.names(Summ_tOTU)


# subset
Fall_spp <- subset(Meta_Soil, Timepoint == "Fall")

#Remove extra metadata
Fall_OTU <- Fall_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozero<- Fall_OTU[, colSums(Fall_OTU != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTU <- t(Fall_OTUs_nozero)
Fall_OTU_list <- row.names(Fall_tOTU)

# subset
Winter_spp <- subset(Meta_Soil, Timepoint == "Winter")

#Remove extra metadata
Winter_OTU <- Winter_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Winter_OTUs_nozero<- Winter_OTU[, colSums(Winter_OTU != 0) > 0]

#Transform back and get list of OTUs
Winter_tOTU <- t(Winter_OTUs_nozero)
Winter_OTU_list <- row.names(Winter_tOTU)


# subset
Spring_spp <- subset(Meta_Soil, Timepoint == "Spring")

#Remove extra metadata
Spring_OTU <- Spring_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Spring_OTUs_nozero<- Spring_OTU[, colSums(Spring_OTU != 0) > 0]

#Transform back and get list of OTUs
Spring_tOTU <- t(Spring_OTUs_nozero)
Spring_OTU_list <- row.names(Spring_tOTU)


## All Seasons, Soil ##

venn_colors <- c("gold1","aquamarine1","firebrick1","green2")
### Venn Diagrams just output straight to the folder you're in
venn.diagram(
  x = list(Summ_OTU_list,Winter_OTU_list,Fall_OTU_list,Spring_OTU_list),
  category.names = c("Summer","Winter","Fall","Spring"),
  filename = 'YJAMF_47VT_Genbank_Soil_venn_no_subset.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1080 ,
  width = 1080 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  main = "Soil AMF Richness Overlap",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)











# subset Timepoints for Root #####
Summer_sppr <- subset(Meta_Root, Timepoint == "Summer")

#Remove extra metadata
Summer_OTUr <- Summer_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozeror<- Summer_OTUr[, colSums(Summer_OTUr != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTUr <- t(Summer_OTUs_nozeror)
Summ_OTU_listr <- row.names(Summ_tOTUr)


# subset
Fall_sppr <- subset(Meta_Root, Timepoint == "Fall")

#Remove extra metadata
Fall_OTUr <- Fall_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozeror<- Fall_OTUr[, colSums(Fall_OTUr != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTUr <- t(Fall_OTUs_nozeror)
Fall_OTU_listr <- row.names(Fall_tOTUr)

# subset
Winter_sppr <- subset(Meta_Root, Timepoint == "Winter")

#Remove extra metadata
Winter_OTUr <- Winter_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Winter_OTUs_nozeror<- Winter_OTUr[, colSums(Winter_OTUr != 0) > 0]

#Transform back and get list of OTUs
Winter_tOTUr <- t(Winter_OTUs_nozeror)
Winter_OTU_listr <- row.names(Winter_tOTUr)


# subset
Spring_sppr <- subset(Meta_Root, Timepoint == "Spring")

#Remove extra metadata
Spring_OTUr <- Spring_sppr[, c(1:47)]

#Remove OTUs that aren't present in any samples
Spring_OTUs_nozeror<- Spring_OTUr[, colSums(Spring_OTUr != 0) > 0]

#Transform back and get list of OTUs
Spring_tOTUr <- t(Spring_OTUs_nozeror)
Spring_OTU_listr <- row.names(Spring_tOTUr)

## All Seasons, Root ##

venn_colors <- c("gold1","aquamarine1","firebrick1","green2")
### Venn Diagrams just output straight to the folder you're in
venn.diagram(
  x = list(Summ_OTU_listr,Winter_OTU_listr,Fall_OTU_listr,Spring_OTU_listr),
  category.names = c("Summer","Winter","Fall","Spring"),
  filename = 'YJAMF_Genbank_Root_venn_47VT_Nosubsampling.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 1080 ,
  width = 1080 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  main = "Root AMF Richness Overlap",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)


####### DF vs YJAMF #########

#
#
# # DF_AMF_table = read.csv("DFSoil_AMF_Table.csv", row.names = 1)
# row.names(DF_AMF_table)
# DF_AMF_tbl = DF_AMF_table[-c(1,2,4,5,13,14,20,48),]
# DF_AMF_t = as.data.frame(t(DF_AMF_tbl))
# DF_AMF_t$id = row.names(DF_AMF_t)
# DF_meta = read.csv("metadata/DomeFire-SamplingInfo-Metadata.csv", row.names = 1)
# DF_meta$id = row.names(DF_meta)
# DF_meta$id = gsub("DF","S", DF_meta$id)
#
# DF_VT_meta = inner_join(DF_AMF_t, DF_meta, by = "id")
#
# # Create VT lists
#
#
# Burn_spp <- subset(DF_VT_meta, Treatment == "Burned")
# #Remove extra metadata
# Burn_VT <- Burn_spp[, c(1:23)]
#
# #Remove OTUs that aren't present in any samples
# Burn_VTs_nozero<- Burn_VT[, colSums(Burn_VT != 0) > 0]
#
# #Transform back and get list of OTUs
# BURN_VTt <- t(Burn_VTs_nozero)
# BURN_vt_list <- row.names(BURN_VTt)
#
#
# Unburn_spp <- subset(DF_VT_meta, Treatment == "Unburned")
# #Remove extra metadata
# Unburn_VT <- Unburn_spp[, c(1:23)]
#
# #Remove OTUs that aren't present in any samples
# Unburn_VTs_nozero<- Unburn_VT[, colSums(Unburn_VT != 0) > 0]
#
# #Transform back and get list of OTUs
# Unburn_VTt <- t(Unburn_VTs_nozero)
# Unburn_vt_list <- row.names(Unburn_VTt)
#
#
#
# venn.diagram(
#   x = list(Unburn_vt_list,BURN_vt_list),
#   category.names = c("Unburned","Burned"),
#   filename = 'DFAMF_BvU_Alltime_venn_diagramm.png',
#   output=TRUE,
#
#   # Output features
#   imagetype="png" ,
#   height = 720 ,
#   width = 720 ,
#   resolution = 300,
#   compression = "lzw",
#
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = c("red","blue"),
#
#   cex = 2,
#   fontface = "bold",
#   fontfamily = "sans",
#
#   # Set names
#   cat.cex = 0.6,
#   cat.dist = 0.001,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.fontfamily = "sans"
# )
#
#
# # General DF list
#
# DF_VT <- DF_VT_meta[, c(1:23)]
#
# #Remove OTUs that aren't present in any samples
# DF_VT_nozero<- DF_VT[, colSums(DF_VT != 0) > 0]
#
# #Transform back and get list of OTUs
# DFN_VTt <- t(DF_VT_nozero)
# Dome_vt_list <- row.names(DFN_VTt)
#
# meta_AMF
#
# YJ_VT <- meta_AMF[, c(1:47)]
#
# #Remove OTUs that aren't present in any samples
# YJ_VT_nozero<- YJ_VT[, colSums(YJ_VT != 0) > 0]
#
# #Transform back and get list of OTUs
# YJ_VTt <- t(YJ_VT_nozero)
# YJ_vt_list <- row.names(YJ_VTt)
#
# ### UGH Why are the Joshua Tree VTs formatted different. Christ.
#
# x <- "some text in a string"
# Dome_VTonly_list=str_sub(Dome_vt_list,-8,-1)
# YJ_VTonly_list=str_sub(YJ_vt_list,-8,-1)
#
#
#
# venn.diagram(
#   x = list(Dome_VTonly_list,YJ_VTonly_list),
#   category.names = c("Dome Fire","Joshua Trees"),
#   filename = 'DFvYJAMF_Alltime_venn_diagramm.png',
#   output=TRUE,
#
#   # Output features
#   imagetype="png" ,
#   height = 720 ,
#   width = 720 ,
#   resolution = 300,
#   compression = "lzw",
#
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = c("brown","green4"),
#
#   cex = 2,
#   fontface = "bold",
#   fontfamily = "sans",
#
#   # Set names
#   cat.cex = 0.6,
#   cat.dist = 0.001,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.fontfamily = "sans"
# )
# setdiff(YJ_VTonly_list, Unburn_vt_list)
# setdiff(YJ_VTonly_list, BURN_vt_list)
# setdiff(YJ_VTonly_list, Dome_VTonly_list)
# Unburn_vt_list
#

#### Check Root and Soil vs DF Burned and Unburned Plots #######


# subset
Root_spp <- subset(meta_AMF, Type == "Root")

#Remove extra metadata
Root_OTU <- Root_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Root_OTU_nozero<- Root_OTU[, colSums(Root_OTU != 0) > 0]

#Transform back and get list of OTUs
Root_tOTU <- t(Root_OTU_nozero)
Root_tOTU_list <- row.names(Root_tOTU)


Burn_VTonly_list=str_sub(BURN_vt_list,-8,-1)
Unburn_VTonly_list=str_sub(Unburn_vt_list,-8,-1)
YJRoot_VTonly_list=str_sub(Root_tOTU_list,-8,-1)



venn.diagram(
  x = list(Burn_VTonly_list,YJRoot_VTonly_list),
  category.names = c("Dome Fire Burned Soils","Joshua Tree Roots"),
  filename = 'DFvYJAMF_burnVroot_venn_diagramm.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 720 ,
  width = 720 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("brown","green4"),

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.dist = 0.001,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(Unburn_VTonly_list,YJRoot_VTonly_list),
  category.names = c("Dome Fire Unburned Soils","Joshua Tree Roots"),
  filename = 'DFvYJAMF_unburnVroot_venn_diagramm.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 720 ,
  width = 720 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("brown","green4"),

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.dist = 0.001,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)



# subset Treatment #####
Burned_YJAMF <- subset(meta_AMF, Treatment == "Burned")

#Remove extra metadata
Summer_OTU <- Burned_YJAMF[, c(1:47)]

#Remove OTUs that aren't present in any samples
Summer_OTUs_nozero<- Summer_OTU[, colSums(Summer_OTU != 0) > 0]

#Transform back and get list of OTUs
Summ_tOTU <- t(Summer_OTUs_nozero)
BurnedYJAMF_list <- row.names(Summ_tOTU)


# subset treatment for unburned
Unburned_spp <- subset(meta_AMF, Treatment == "Unburned")

#Remove extra metadata
Fall_OTU <- Unburned_spp[, c(1:47)]

#Remove OTUs that aren't present in any samples
Fall_OTUs_nozero<- Fall_OTU[, colSums(Fall_OTU != 0) > 0]

#Transform back and get list of OTUs
Fall_tOTU <- t(Fall_OTUs_nozero)
Unburned_spp_list <- row.names(Fall_tOTU)

table(meta_AMF$Treatment)
venn.diagram(
  x = list(Unburned_spp_list,BurnedYJAMF_list),
  category.names = c("Unburned Trees\nn=116","Burned Trees\nn=26"),
  filename = 'YJAMF_Treatment_Overlap_venn_diagramm.png',
  output=TRUE,

  # Output features
  imagetype="png" ,
  height = 720 ,
  width = 720 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("blue","red"),

  cex = 2,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.5,
  cat.dist = 0.001,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)


setdiff(Unburned_spp_list,BurnedYJAMF_list)
#Exclusive to Unburned Trees (vastly outnumbered burned trees)
# [1] "k__Fungi..p__Mucoromycota..c__Archaeosporomycetes..o__Archaeosporales..f__Archaeosporaceae..g__Archaeospora..s__VTX00376"
# [2] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Diversisporales..f__Diversisporaceae..g__Diversispora..s__VTX00354"
# [3] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Diversisporales..f__Diversisporaceae..g__Diversispora..s__VTX00355"
# [4] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Claroideoglomeraceae..g__Claroideoglomus..s__VTX00402"
# [5] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00104"
# [6] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00150"
# [7] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00165"
# [8] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00174"
# [9] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00188"
# [10] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00248"
# [11] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00364"
# [12] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00414

setdiff(BurnedYJAMF_list,Unburned_spp_list)
# [1] "k__Fungi..p__Mucoromycota..c__Glomeromycetes..o__Glomerales..f__Glomeraceae..g__Glomus..s__VTX00093"
meta_AMF[,c(14)]
# VTX00093 only present as 12 reads anyway



############################ Harrower Comparison #####################
# used image extractor to pull VTs from Harrower and Gilbert 2021
# They assessed Y. brevifolia AMF from 100km and across an elevation gradient
# Used NS31-AML2 which is about as close to WANDA-AML2 as you can get
# Ambispora sp.:VTX00242
# Ambispora callosa:VTX00283
# Paraglomus sp.:VTX00001
# Paraglomus sp.:VTX00281
# Paraglomus sp.:VTX00336
# Claroideoglomus sp.:VTX00193
# Claroideoglomus sp.:VTX00056
# Claroideoglomus sp.:VTX00278
# Diversispora sp.:VTX00062
# Diversispora sp.:VTX00061
# Diversispora sp.:VTX00306
# Acaulospora sp.:VTX00231
# Acaulospora sp.:VTX00028
# Scutellospora sp.VTX00049
# Scutellospora sp.:VTX00052
# Glomus sp.:VTX00103
# Glomus MO-G31:VTX00191
# Glomus MO-G7:VTX00199
# Glomus sp.:VTX00196
# Glomus sp.:VTX00149
# Glomus sp.:VTX00129
# Glomus sp.:VTX00140
# Glomus MO-G20:VTX00143
# Glomus sp.:VTX00137
# Glomus sp.:VTX00053
# Glomus sp.:VTX00093
# Glomus sp.:VTX00069
# Glomus sp.:VTX00089
# Glomus sp.:VTX00074
# Glomus:VTX00404
# Glomus Yamato08:VTX00100
# Glomus MO-G54:VTX00092
# Glomus sp.:VTX00105
# Glomus sp.:VTX00087
# Glomus sp.:VTX00114
# Glomus sp.:VTX00113
# Glomus Glo-A:VTX00295

# Yucca brevifolia virtual taxa codes
YBHGAMF = c("VTX00242",
            "VTX00283",
            "VTX00001",
            "VTX00281",
            "VTX00336",
            "VTX00193",
            "VTX00056",
            "VTX00278",
            "VTX00062",
            "VTX00061",
            "VTX00306",
            "VTX00231",
            "VTX00028",
            "VTX00052",
            "VTX00103",
            "VTX00191",
            "VTX00199",
            "VTX00196",
            "VTX00149",
            "VTX00129",
            "VTX00140",
            "VTX00049",
            "VTX00143",
            "VTX00137",
            "VTX00053",
            "VTX00093",
            "VTX00069",
            "VTX00089",
            "VTX00074",
            "VTX00404",
            "VTX00100",
            "VTX00092",
            "VTX00105",
            "VTX00087",
            "VTX00114",
            "VTX00113",
            "VTX00295"
)
# get virtual taxa from Yucca jaegeriana
#Import ASV Table ##########
colnames(asv_AMF.e754_nozeroes)
lastcolumn <- colnames(asv_AMF.e754_nozeroes[,c(1:47)])

#Data frame with the OTU IDs and the taxa associated with it
dd <- data.frame(c(1:47),lastcolumn)


otu_taxa1 <- ldply(str_split(string = dd$lastcolumn, pattern="__"), rbind)

# Divide a column using ";"and convert list to data frame
names(otu_taxa1) <- c("Genus", "Species")
# VT names from another script; print here
SpeciesNames = otu_taxa1$Species
SpeciesNames
# [1] "VTX00376" "VTX00054" "VTX00058" "VTX00061" "VTX00354" "VTX00355" "VTX00056" "VTX00193"
# [9] "VTX00357" "VTX00402" "VTX00063" "VTX00069" "VTX00092" "VTX00093" "VTX00098" "VTX00099"
# [17] "VTX00104" "VTX00114" "VTX00130" "VTX00150" "VTX00155" "VTX00165" "VTX00174" "VTX00188"
# [25] "VTX00195" "VTX00214" "VTX00234" "VTX00248" "VTX00294" "VTX00326" "VTX00364" "VTX00409"
# [33] "VTX00414" "VTX00442" "VTX00448" "VTX00444"
str(SpeciesNames)

duplicates <- intersect(YBHGAMF, SpeciesNames)

# Print the duplicates
cat("The duplicates are: ", duplicates, "\n")
# The duplicates are:
# VTX00193 VTX00056 VTX00061 VTX00306 VTX00093 VTX00069 VTX00092 VTX00105 VTX00114




######## Do tree core community per timepoint; start with Fall


# Don't use rarefied data for this
asv_AMF_norar = read.csv("YJAMF_Merged_47VTonly_norarefaction.csv", row.names = 1)
fungimetadata <- read.csv("YJ_Merged_Metadata-nodash.csv", row.names=1)
fungimetadata$id <- row.names(fungimetadata)
asv_AMF_norar <-as.data.frame(asv_AMF_norar)
asv_AMF_norar$id <-row.names(asv_AMF_norar)
meta_AMF <-inner_join(asv_AMF_norar, fungimetadata, by="id")
row.names(meta_AMF) <- meta_AMF$id
meta_AMF = subset(meta_AMF, Project == "YJAMF")

meta_AMF$Timepoint <-revalue(meta_AMF$Timepoint, c("T1"= "Summer", "T2"= "Fall", "T3" = "Winter", "T4" = "Spring"))
table(meta_AMF$Timepoint)


# Subset to Fall samples
Fall_Fungi= subset(meta_AMF, Timepoint == "Fall")
str(Fall_Fungi)
Fallsub1=Fall_Fungi[,-c(48:49)]
Fallsub2=Fallsub1[,-c(49:73)]

FallRS_TreeNums <- aggregate(. ~ TreeNum, data = Fallsub2, FUN = sum)

Core_Fall <- FallRS_TreeNums %>%
  select(where(~ all(. != 0)))


# Subset to Winter samples
Winter_Fungi= subset(meta_AMF, Timepoint == "Winter")
str(Winter_Fungi)
Wintersub1=Winter_Fungi[,-c(48:49)]
Wintersub2=Wintersub1[,-c(49:73)]

WinterRS_TreeNums <- aggregate(. ~ TreeNum, data = Wintersub2, FUN = sum)

Core_Winter <- WinterRS_TreeNums %>%
  dplyr::select(where(~ all(. != 0)))
# none


# Subset to Spring samples
Spring_Fungi= subset(meta_AMF, Timepoint == "Spring")
str(Spring_Fungi)
Springsub1=Spring_Fungi[,-c(48:49)]
Springsub2=Springsub1[,-c(49:73)]

SpringRS_TreeNums <- aggregate(. ~ TreeNum, data = Springsub2, FUN = sum)

Core_Spring <- SpringRS_TreeNums %>%
  select(where(~ all(. != 0)))

# Subset to Summer samples
Summer_Fungi= subset(meta_AMF, Timepoint == "Summer")
str(Summer_Fungi)
Summersub1=Summer_Fungi[,-c(48:49)]
Summersub2=Summersub1[,-c(49:73)]

SummerRS_TreeNums <- aggregate(. ~ TreeNum, data = Summersub2, FUN = sum)

Core_Summer <- SummerRS_TreeNums %>%
  select(where(~ all(. != 0)))


#### Turns out Fall is the only time Glomus 294 is in every tree
#### Same results whether or not you use the rarefied data
