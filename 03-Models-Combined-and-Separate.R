
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

psych::describe(AMF_data$S.obs)
psych::describe(subset(AMF_soildata, Timepoint == "Summer")$S.obs)
psych::describe(subset(AMF_soildata, Timepoint == "Fall")$S.obs)

psych::describe(AMF_rootdata$S.obs)
psych::describe(subset(AMF_rootdata, Timepoint == "Summer")$S.obs)
psych::describe(subset(AMF_rootdata, Timepoint == "Fall")$S.obs)
psych::describe(subset(AMF_rootdata, Timepoint == "Winter")$S.obs)
psych::describe(subset(AMF_rootdata, Timepoint == "Spring")$S.obs)


psych::describe(subset(AMF_data, Timepoint == "Summer")$S.obs)
psych::describe(subset(AMF_data, Timepoint == "Fall")$S.obs)
psych::describe(subset(AMF_data, Timepoint == "Winter")$S.obs)
psych::describe(subset(AMF_data, Timepoint == "Spring")$S.obs)

#### Manual Model Selection for Richness comparisons ####

# Log transform doesn't normalize data
shapiro.test(AMF_data$S.obs)
histogram(AMF_data$S.obs)

car::leveneTest(S.obs ~ interaction(Timepoint,Type), AMF_data)
# Though Levenes test not significant

str(AMF_data)
# Shapiro-Wilk test significant so check distribution
norm = fitdistr(AMF_data$S.obs, "normal")
lognorm = fitdistr(AMF_data$S.obs, "lognormal")
normmodel = qqp(AMF_data$S.obs, "norm")
lognormmodel = qqp(AMF_data$S.obs, "lnorm")
nbinom <- fitdistr(AMF_data$S.obs, "Negative Binomial")
negbimodel=qqp(AMF_data$S.obs, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
poisson <- fitdistr(AMF_data$S.obs, "Poisson")
poismodel = qqp(AMF_data$S.obs, "pois", poisson$estimate, lambda = poisson$estimate)
gamma <- fitdistr(AMF_data$S.obs, "gamma")
gammamodel = qqp(AMF_data$S.obs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#pick the distribution for which the largest number of observations falls between the dashed lines.

qqp(AMF_data$S.obs, "norm")
qqp(AMF_data$S.obs, "lnorm")
qqp(AMF_data$S.obs, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
qqp(AMF_data$S.obs, "pois", poisson$estimate, lambda = poisson$estimate)
qqp(AMF_data$S.obs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

AIC(norm, lognorm, poisson, nbinom, gamma)

# Poisson distribution best fit
str(AMF_data)
# Check random effects (Different trees, 2 Illumina Libraries)
gm1 = glmer(S.obs ~  (1 | TreeNum) +  (1|Library),  data = AMF_data, family = poisson(link = "log"))
gm2 = glmer(S.obs ~  (1 | TreeNum),  data = AMF_data, family = poisson(link = "log"))
gm3 = glmer(S.obs ~  (1 | Library),  data = AMF_data, family = poisson(link = "log"))

AIC(gm1, gm2,gm3)
#     df      AIC
# gm1  3 724.0033
# gm2  2 724.0223
# gm3  2 737.6252
# Very similar, just do tree number

str(AMF_data)
# Account for TreeNum in poisson-distributed model
Model1 =glmer(S.obs ~ Treatment*Type*Timepoint + PrecipThreeMonths.in + (1 | TreeNum) ,
              data = AMF_data, family = poisson(link = "log"))

Model2 =glmer(S.obs ~ Treatment*Timepoint +(1 | TreeNum) ,  data = AMF_data,
              family = poisson(link = "log"))

Modeltwo =glmer(S.obs ~ Treatment*Type+ PrecipTwoMonths +(1 | TreeNum) ,
                data = AMF_data, family = poisson(link = "log"))

Modelthree =glmer(S.obs ~ Treatment*Type+ PrecipThreeMonths.in +(1 | TreeNum) ,
                  data = AMF_data, family = poisson(link = "log"))

Modelrain3 =glmer(S.obs ~ Type + PrecipThreeMonths.in +(1 | TreeNum) ,
                  data = AMF_data, family = poisson(link = "log"));summary(Modelrain3)

Modelrain4 =glmer(S.obs ~ Type+ PrecipFourMonths +(1 | TreeNum) ,
                  data = AMF_data, family = poisson(link = "log"));summary(Modelrain4)

AIC(Model1,Model2,Modeltwo,Modelrain3,Modelrain4)
summary(Modelthree)
anova(Modelthree)

# Burn effect (Treatment) is not significant so test without it
model_time = glmer(S.obs ~ Timepoint  +(Timepoint | TreeNum),
                   data = AMF_data, family = poisson(link = "log"))

model_type = glmer(S.obs ~ Type +(1 | TreeNum) ,
                   data = AMF_data, family = poisson(link = "log"))

model_typetime = glmer(S.obs ~ Type*Timepoint + (1 | TreeNum),
                       data = AMF_data, family = poisson(link = "log"))

model_typeandtime = glmer(S.obs ~ Type + Timepoint + (1 | TreeNum),
                          data = AMF_data, family = poisson(link = "log"))

model_typerain3 = glmer(S.obs ~ Type + PrecipThreeMonths.in + (1 | TreeNum),
                        data = AMF_data, family = poisson(link = "log"))

model_typetimerain3 = glmer(S.obs ~ Type*Timepoint + PrecipThreeMonths.in + (1 | TreeNum),
                            data = AMF_data, family = poisson(link = "log"))


AIC(model_typerain3,model_typetimerain3,Model1,Model2,Modeltwo,Modelthree,
    model_typeandtime,model_typetime,model_type,model_time)

summary(model_typerain3)

summary(model_typetimerain3)

summary(model_typeandtime)
summary(model_typetime)
AIC(model_typetime,model_typeandtime)

emmtime <- emmeans(model_typeandtime, ~  Timepoint)
comparisontime <- pairs(emmtime,by = NULL, adjust = "bonf")
comparisontime2 <- pairs(emmtime)
summary(comparisontime)
summary(comparisontime2)
model_means_cld <- multcomp::cld(object = emmtime,
                                 # adjust = "tukey",
                                 Letters = letters,
                                 alpha = 0.05)
model_means_cld

summary(model_typetime)
summary(model_typetime4)
summary(model_typeandtime)
summary(model_typeandtime4)
summary(model_type)
summary(model_time)

anova(Model1,Model2,Model3)

stargazer(model_typetime, model_typeandtime, type ="html", out = "AMF_Rich_Modcompare.htm")

summary(model_typeandtime)
summary(Modelfour)

anova(Model1,Modeltwo,Modelfour)

# Compare to the null
anova(model_typeandtime,model_typetime,gm2)

# Get all pairwise comparisons between Timepoint and Type
summary(model_typeandtime)
emm <- emmeans(model_typeandtime, ~ Type)
comparisons <- pairs(emm)
summary(comparisons)
# contrast    estimate     SE  df z.ratio p.value
# Root - Soil   -0.153 0.0536 Inf  -2.861  0.0042

emmtime <- emmeans(model_typeandtime, ~  Timepoint)
comparisontime <- pairs(emmtime)
summary(comparisontime)
# contrast        estimate     SE  df z.ratio p.value
# Summer - Fall    -0.1938 0.0788 Inf  -2.459  0.0665
# Summer - Winter  -0.1794 0.0812 Inf  -2.208  0.1211
# Summer - Spring  -0.1564 0.0800 Inf  -1.956  0.2048
# Fall - Winter     0.0144 0.0728 Inf   0.198  0.9973
# Fall - Spring     0.0373 0.0713 Inf   0.524  0.9534
# Winter - Spring   0.0229 0.0739 Inf   0.310  0.9896
emmtypetime <- emmeans(model_typeandtime, ~ Type +  Timepoint)
comparisontime <- pairs(emmtypetime)
summary(comparisontime)

sink(file = "YJAMF_MergedVT49_FinalModele754_Stats.txt", append=TRUE)
summary(model_typeandtime)
summary(comparisons)
sink()

model_means_cld <- cld(object = emm,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# show output
model_means_cld
sink(file = "YJAMF_MergedVT49_FinalModele754_CLD_Letters.txt", append=TRUE)
model_means_cld
sink()

model_time_cld <- cld(object = emmtime,
                      adjust = "sidak",
                      Letters = letters,
                      alpha = 0.05)
model_time_cld


### Subset to do Root only #####

#Test different probability distributions for data

AMF_rootdata = subset(AMF_data, Type == "Root")
shapiro.test(AMF_rootdata$S.obs)
leveneTest(S.obs ~ Timepoint, AMF_rootdata)

norm = fitdistr(AMF_rootdata$S.obs, "normal")
lognorm = fitdistr(AMF_rootdata$S.obs, "lognormal")

normmodel = qqp(AMF_rootdata$S.obs, "norm")
lognormmodel = qqp(AMF_rootdata$S.obs, "lnorm")
nbinom <- fitdistr(AMF_rootdata$S.obs, "Negative Binomial")
negbimodel=qqp(AMF_rootdata$S.obs, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
poisson <- fitdistr(AMF_rootdata$S.obs, "Poisson")
poismodel = qqp(AMF_rootdata$S.obs, "pois", poisson$estimate, lambda = poisson$estimate)
gamma <- fitdistr(AMF_rootdata$S.obs, "gamma")
gammamodel = qqp(AMF_rootdata$S.obs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

qqp(AMF_rootdata$S.obs, "norm")
qqp(AMF_rootdata$S.obs, "lnorm")
qqp(AMF_rootdata$S.obs, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
qqp(AMF_rootdata$S.obs, "pois", poisson$estimate, lambda = poisson$estimate)
qqp(AMF_rootdata$S.obs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

AIC(norm, lognorm, poisson, nbinom, gamma)
str(AMF_rootdata)

gm1 =  lmer(S.obs ~  (1 | TreeNum) +  (1|Library),  data = AMF_rootdata )
gm2 =  lmer(S.obs ~  (1 | TreeNum),  data = AMF_rootdata )
gm3 =  lmer(S.obs ~  (1|Library),  data = AMF_rootdata )

AIC(gm1, gm2,gm3)
summary(gm2)

# Account for TreeNum
Model1 = lmer(S.obs ~ Treatment*Timepoint +(1 | TreeNum) ,  data = AMF_rootdata )
model_time =  lmer(S.obs ~ Timepoint  +(1 | TreeNum) ,  data = AMF_rootdata )
model_burn =  lmer(S.obs ~ Treatment +(1 | TreeNum) ,  data = AMF_rootdata )
model_rain3root =  lmer(S.obs ~ PrecipThreeMonths.in  +(1 | TreeNum) ,  data = AMF_rootdata );summary(model_rain3root)
model_rain4root =  lmer(S.obs ~ PrecipFourMonths  +(1 | TreeNum) ,  data = AMF_rootdata );summary(model_rain4root)
model_burntime =  lmer(S.obs ~ Treatment*Timepoint +(1 | TreeNum) ,  data = AMF_rootdata )

AIC(Model1,model_time,model_burn)

# model_type lower AIC but best one including time is the time*type interaction
summary(Model1)
summary(model_rain)
summary(model_time)
summary(model_burn)

# Compare to the null
anova(model_time,gm2)
anova(Model1,gm2)

emm <- emmeans(Model1, ~ Timepoint*Treatment)
comparisons <- pairs(emm)
summary(comparisons)

model_means_cld <- cld(object = emm,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# show output
model_means_cld

# So the Season*Burn shows that everything is the same in roots
## Compare to just soil seasonal model

emmtime <- emmeans(model_time, ~ Timepoint)
comparisonstime <- pairs(emmtime)
summary(comparisonstime)

model_means_cldTime <- cld(object = emmtime,
                           adjust = "sidak",
                           Letters = letters,
                           alpha = 0.05)
model_means_cldTime
### In root, no effect
# Timepoint emmean    SE   df lower.CL upper.CL .group
# Summer      8.57 0.727 60.2     6.71     10.4  a
# Winter      8.93 0.673 55.2     7.20     10.7  a
# Fall        9.45 0.644 52.0     7.79     11.1  a
# Spring      9.74 0.658 53.6     8.04     11.4  a


#### Report root data ########
# Treat time as correlated and seasonal sampling equidistant from one another
rootmodel_repeated_measure <- lme(S.obs ~ Treatment * Timepoint,
                                  random = ~ 1 | TreeNum,
                                  correlation = corAR1(form = ~ as.numeric(Timepoint) | TreeNum),
                                  data = AMF_rootdata)
anova(rootmodel_repeated_measure)


### Soil only #####


#Test different probability distributions for data

AMF_soildata= subset(AMF_data, Type =="Soil")
shapiro.test(AMF_soildata$S.obs)
leveneTest(S.obs ~ Timepoint, AMF_soildata)

qqp(AMF_soildata$S.obs, "norm")
qqp(AMF_soildata$S.obs, "lnorm")
qqp(AMF_soildata$S.obs, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
qqp(AMF_soildata$S.obs, "pois", poisson$estimate, lambda = poisson$estimate)
qqp(AMF_soildata$S.obs, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

AIC(norm, lognorm, poisson, nbinom, gamma)
str(AMF_soildata)

gm1 = lmer(S.obs ~  (1 | TreeNum) +  (1|Library),  data = AMF_soildata)
gm2 = lmer(S.obs ~  (1 | TreeNum),  data = AMF_soildata)
gm3 = lmer(S.obs ~  (1|Library),  data = AMF_soildata)

AIC(gm1, gm2,gm3)

Model1 =lmer(S.obs ~ Treatment*Timepoint +(1 | TreeNum) ,  data = AMF_soildata )
model_time = lmer(S.obs ~ Timepoint  +(1 | TreeNum) ,  data = AMF_soildata )
model_burn = lmer(S.obs ~ Treatment +(1 | TreeNum) ,  data = AMF_soildata )
model_rain = lmer(S.obs ~ PrecipThreeMonths.in  +(1 | TreeNum) ,  data = AMF_soildata );summary(model_rain)
model <- lmer(S.obs ~ Treatment * Timepoint + (Timepoint | TreeNum), data = AMF_rootdata)

summary(Model1)
summary(model_rain)
summary(model_time)
summary(model_burn)


AIC(Model1,model_time,model_rain,model_burn)
# df      AIC
# Model1     10 349.7708
# model_time  6 358.6562
# model_rain  4 371.5778
# model_burn  4 369.5477

# Model 1 with a Time and Treatment interaction is the lowest AIC
# but the model with only burn shows no effect
# So there is no overall burn effect on soil AMF richness but an interaction with time

summary(Model1)

# Compare to the null
anova(Model1,gm2)
anova(model_time,gm2)

# Better than null
library(emmeans)
library(multcomp)
emm <- emmeans(Model1, ~ Timepoint*Treatment)
comparisons <- pairs(emm)
summary(comparisons)

model_means_cld <- cld(object = emm,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

# show output
model_means_cld

# So the Season*Burn shows that everything is close,
# except Summer unburned is lower and Fall unburned is higher

## Compare to just soil seasonal model

emmtime <- emmeans(model_time, ~ Timepoint)
comparisonstime <- pairs(emmtime)
summary(comparisonstime)

model_means_cldTime <- cld(object = emmtime,
                           adjust = "tukey",
                           Letters = letters,
                           alpha = 0.05)
model_means_cldTime
# Timepoint emmean    SE   df lower.CL upper.CL .group

# Summer      8.71 0.789 53.4     6.68     10.7  a
# Spring     10.46 0.737 48.0     8.55     12.4  ab
# Fall       11.55 0.738 48.0     9.64     13.5   b
# Winter     11.68 0.770 51.6     9.70     13.7   b

### Looking at seasons alone,
# Summer is low
# Spring is ab to summer and rest of seasons
# But fall and winter are different from summer
# Report lme model

soilmodel_repeated_measure <- lme(S.obs ~  Timepoint,
                                  random = ~ 1 | TreeNum,
                                  correlation = corAR1(form = ~ as.numeric(Timepoint) | TreeNum),
                                  data = AMF_soildata)
asoil = anova(soilmodel_repeated_measure)
rootmodel_repeated_measure <- lme(S.obs ~  Timepoint,
                                  random = ~ 1 | TreeNum,
                                  correlation = corAR1(form = ~ as.numeric(Timepoint) | TreeNum),
                                  data = AMF_rootdata)
summary(rootmodel_repeated_measure)
aroot = anova(rootmodel_repeated_measure)

stargazer(soilmodel_repeated_measure,
          rootmodel_repeated_measure,
          dep.var.labels=c("Soil VT Richness","Root VT Richness"),
          type = "html"
)

emmtime <- emmeans(soilmodel_repeated_measure, ~ Timepoint)
comparisonstime <- pairs(emmtime)
summary(comparisonstime)

model_means_cldTime <- cld(object = emmtime,
                           adjust = "tukey",
                           Letters = letters,
                           alpha = 0.05)
model_means_cldTime
# Soil
# Timepoint emmean    SE df lower.CL upper.CL .group
# Summer      8.68 0.789 19     6.51     10.9  a
# Spring     10.47 0.739 19     8.44     12.5  ab
# Fall       11.55 0.739 19     9.51     13.6   b
# Winter     11.69 0.769 19     9.58     13.8   b

emmtime <- emmeans(rootmodel_repeated_measure, ~ Timepoint)
comparisonstime <- pairs(emmtime)
summary(comparisonstime)

model_means_cldTime <- cld(object = emmtime,
                           adjust = "tukey",
                           Letters = letters,
                           alpha = 0.05)
model_means_cldTime
# Root
# Timepoint emmean    SE df lower.CL upper.CL .group
# Summer      8.94 0.842 18     6.61     11.3  a
# Winter      9.40 0.824 18     7.12     11.7  a
# Spring      9.48 0.818 18     7.21     11.7  a
# Fall        9.75 0.814 18     7.50     12.0  a

