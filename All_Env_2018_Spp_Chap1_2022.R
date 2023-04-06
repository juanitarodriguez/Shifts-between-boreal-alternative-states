## Script information
### Script of the article entitled: Drivers of contrasting boreal understory vegetation in coniferous and broadleaf deciduous alternative states
### Article submitted to Ecological Monographs
### Authors of the article: Juanita C. Rodríguez-Rodríguez, Nicole J. Fenton, Steven W. Kembel, Evick Mestre, Mélanie Jean, and Yves Bergeron
### Script developed on 2021 by: Juanita C. Rodríguez-Rodríguez (e-mail: juanitacarolina.rodriguezrodriguez@uqat.ca), for my PhD in Environmental Sciences. A part of the script was developed during the statistical cours of Pierre Legendre.
### Environmental data analysis for Chap1 data (I just have soil samples from 4 treatments but the other information for all treatments)

library(ggplot2)
library(stats)
library(nlme)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(lme4)
library(dplyr)
library(rcompanion) #For cldList 
library(lsmeans)
library(emmeans)
library(lmerTest) #To include p values in lmer!

### Soil data
setwd("Z:/UQAT/PhD PROJECT/5. Chapter 1 - OsV-UsV/0. Field work data/8. R Analysis/Env_variables")

## Load all data
#DBEnv <- read.csv2("All_Env_2018_Spp.csv")
DBEnv <- read.csv2("All_Env_2018_Spp(DatInput).csv")
str(DBEnv)
#Example for more changes: #B_EPN_Herbs <- mutate(EPN_Herbs, Site = as.factor(Site), Block = as.factor(Block))

## ANOVA tests for each environmental variable

## Filter table to get only soil physicochemical properties from treatments C,1F,Ti,To
### Separate by Treatments
Treats <- split(DBEnv, DBEnv$Treatment)
Tr_C <- Treats$'C'
Tr_1F <- Treats$'1F'
Tr_Ti <- Treats$'Ti'
Tr_To <- Treats$'To'

Treats4 <- rbind(Tr_C,Tr_1F,Tr_Ti,Tr_To)

#Order levels of the factor  so I do  comparisons of C versus all
DBEnv$Treatment <- factor(DBEnv$Treatment, levels=c("C","Li","Nu","1F","2F","To","Ti"))
DBEnv$Site <- as.factor(DBEnv$Site)
DBEnv$Block <- as.factor(DBEnv$Block)

### Separate by Canopy (Vordach)
Vor <- split(DBEnv, DBEnv$Canopy)
Vor_BS <- Vor$'BS'
Vor_TA <- Vor$'TA'

Farben <- c("azure3","gold","mediumpurple1","tan1","lightcoral","chartreuse","paleturquoise2","azure4","goldenrod","darkorchid3","orangered","red3","lightseagreen","forestgreen")
Farben_Treats4 <- c("azure3","lightcoral","paleturquoise2","chartreuse")

### Analysis per variable
## STEPS:
# STEP1: Get mean and SD for each variable per Canopy for the different treatments
# STEP2: Do lineal model and anova to get differences between Canopy, Treatment and the interaction of both
# STEP3: Make contrasts emmeans to assign letters and if Treatments or the interaction were also different (not only Canopy), then check the contrasts for treatments to assign capital letters
# STEP4: Figure

#### Hum ####
## STEP1 (meand and SD) Hum
# BS Hum
aggregate(Vor_BS$Hum, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Hum, list(Vor_BS$Treatment), FUN=sd)
# TA Hum
aggregate(Vor_TA$Hum, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Hum, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lmer model)) Hum
Hum_mod <- lmer(Hum ~ Canopy*Treatment + (1|Site/Block), data = DBEnv)
anova(Hum_mod)

### STEP3 (contrasts emmeans of Contro vs. Treatments for each Canopy) Hum
Hum_emm <- emmeans(Hum_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
#Letters
Hum_emm_let <- multcomp::cld(object = Hum_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
Hum_emm_let

### STEP4 (Figure)
Hum_Fig <- DBEnv %>%
  ggplot(aes(Treatment, Hum, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Temp ####
## STEP1 (meand and SD) Temp
# BS Temp
aggregate(Vor_BS$Temp, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Temp, list(Vor_BS$Treatment), FUN=sd)
# TA Temp
aggregate(Vor_TA$Temp, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Temp, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lmer model)) Temp
Temp_mod <- lmer(Temp ~ Canopy*Treatment + (1|Site/Block), data = DBEnv)
anova(Temp_mod)

### STEP3 (contrasts emmeans of Contro vs. Treatments for each Canopy) Temp
Temp_emm <- emmeans(Temp_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
#Letters
Temp_emm_let <- multcomp::cld(object = Temp_emm$emmeans,
                              
                              Letters = letters,no.readonly=TRUE)
Temp_emm_let

### STEP4 (Figure)
Temp_Fig <- DBEnv %>%
  ggplot(aes(Treatment, Temp, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))


#### Light ####
## STEP1 (meand and SD) Light
# BS Light
aggregate(Vor_BS$Light, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Light, list(Vor_BS$Treatment), FUN=sd)
# TA Light
aggregate(Vor_TA$Light, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Light, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lmer model)) Light
Light_mod <- lmer(Light ~ Canopy*Treatment + (1|Site/Block), data = DBEnv)
anova(Light_mod)

### STEP3 (contrasts emmeans of Contro vs. Treatments for each Canopy) Light
Light_emm <- emmeans(Light_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
#Letters
Light_emm_let <- multcomp::cld(object = Light_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
Light_emm_let

### STEP4 (Figure)
Light_Fig <- DBEnv %>%
  ggplot(aes(Treatment, Light, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

Light_Fig2 <- DBEnv %>%
  ggplot(aes(Canopy, Light, group=Treatment, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Treatment)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white")) # PDF 4x6 landscape (To extract light contrast)

#### OverDen ####
## STEP1 (meand and SD) OverDen
# BS OverDen
aggregate(Vor_BS$OverDen, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$OverDen, list(Vor_BS$Treatment), FUN=sd)
# TA OverDen
aggregate(Vor_TA$OverDen, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$OverDen, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lmer model)) OverDen
OverDen_mod <- lmer(OverDen ~ Canopy*Treatment + (1|Site/Block), data = DBEnv)
anova(OverDen_mod)

### STEP3 (contrasts emmeans of Contro vs. Treatments for each Canopy) OverDen
OverDen_emm <- emmeans(OverDen_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
#Letters
OverDen_emm_let <- multcomp::cld(object = OverDen_emm$emmeans,
                                 
                                 Letters = letters,no.readonly=TRUE)
OverDen_emm_let

### STEP4 (Figure)
OverDen_Fig <- DBEnv %>%
  ggplot(aes(Treatment, OverDen, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

OverDen_Fig2 <- DBEnv %>%
  ggplot(aes(Canopy, OverDen, group=Treatment, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Treatment)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white")) # PDF 4x6 landscape (To extract light contrast)


## All Figs together
ggarrange(C_per_Fig,N_per_Fig,K_per_Fig,Mg_per_Fig,Ca_per_Fig,Na_per_Fig,H_per_Fig,P_mg_Fig,Al_mg_Fig,Mn_mg_Fig,Fe_mg_Fig,S_mg_Fig,CEC_Fig,C_N_Fig,N_P_Fig,P_Al_Fig,P_Ca_Fig,pH_H2O_Fig,pH_buf_Fig,Temp_Fig,Hum_Fig,Light_Fig,OverDen_Fig, ncol = 4, nrow = 6, common.legend = TRUE, labels="auto")
## All Figs with sig. differences together PDF 7x10 landscape
All_Env_Fig <- ggarrange(Temp_Fig,Hum_Fig,Light_Fig,OverDen_Fig, ncol = 2, nrow = 2, common.legend = TRUE, labels="auto")
All_Env_Fig

#### Soil physico-chemical properties ####
### Analysis per variable
## STEPS:
# STEP1: Get mean and SD for each variable per Canopy for the different treatments
# STEP2: Do lineal model and anova to get differences between Canopy, Treatment and the interaction of both
# STEP3: Make contrasts emmeans to assign letters and if Treatments or the interaction were also different (not only Canopy), then check the contrasts for treatments to assign capital letters
# STEP4: Figure

#### C_per (%) ####
## STEP1 (meand and SD) C_per
# BS C_per (%)
aggregate(Vor_BS$C_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$C_per, list(Vor_BS$Treatment), FUN=sd)
# TA C_per (%)
aggregate(Vor_TA$C_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$C_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) C_per
C_per_mod <- lmer(C_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(C_per_mod)

### STEP3 (contrasts emmeans) C_per
C_per_emm <- emmeans(C_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
C_per_emm
emmeans(C_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
C_per_emm_let <- multcomp::cld(object = C_per_emm$emmeans,
                                 
                                 Letters = letters,no.readonly=TRUE)
C_per_emm_let

### STEP4 (Figure)
C_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, C_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))


#### N_per (%) ####
## STEP1 (meand and SD) N_per
# BS N_per (%)
aggregate(Vor_BS$N_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$N_per, list(Vor_BS$Treatment), FUN=sd)
# TA N_per (%)
aggregate(Vor_TA$N_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$N_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) N_per
N_per_mod <- lmer(N_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(N_per_mod)

### STEP3 (contrasts emmeans) N_per
N_per_emm <- emmeans(N_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
N_per_emm
emmeans(N_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
N_per_emm_let <- multcomp::cld(object = N_per_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
N_per_emm_let

### STEP4 (Figure)
N_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, N_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### K_per (%) ####
## STEP1 (meand and SD) K_per
# BS K_per (%)
aggregate(Vor_BS$K_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$K_per, list(Vor_BS$Treatment), FUN=sd)
# TA K_per (%)
aggregate(Vor_TA$K_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$K_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) K_per
K_per_mod <- lmer(K_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(K_per_mod)

### STEP3 (contrasts emmeans) K_per
K_per_emm <- emmeans(K_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
K_per_emm
emmeans(K_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
K_per_emm_let <- multcomp::cld(object = K_per_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
K_per_emm_let

### STEP4 (Figure)
K_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, K_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Mg_per (%) ####
## STEP1 (meand and SD) Mg_per
# BS Mg_per (%)
aggregate(Vor_BS$Mg_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Mg_per, list(Vor_BS$Treatment), FUN=sd)
# TA Mg_per (%)
aggregate(Vor_TA$Mg_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Mg_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Mg_per
Mg_per_mod <- lmer(Mg_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Mg_per_mod)

### STEP3 (contrasts emmeans) Mg_per
Mg_per_emm <- emmeans(Mg_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Mg_per_emm
emmeans(Mg_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Mg_per_emm_let <- multcomp::cld(object = Mg_per_emm$emmeans,
                                
                                Letters = letters,no.readonly=TRUE)
Mg_per_emm_let

### STEP4 (Figure)
Mg_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Mg_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Ca_per (%) ####
## STEP1 (meand and SD) Ca_per
# BS Ca_per (%)
aggregate(Vor_BS$Ca_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Ca_per, list(Vor_BS$Treatment), FUN=sd)
# TA Ca_per (%)
aggregate(Vor_TA$Ca_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Ca_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Ca_per
Ca_per_mod <- lmer(Ca_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Ca_per_mod)

### STEP3 (contrasts emmeans) Ca_per
Ca_per_emm <- emmeans(Ca_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Ca_per_emm
emmeans(Ca_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Ca_per_emm_let <- multcomp::cld(object = Ca_per_emm$emmeans,
                                
                                Letters = letters,no.readonly=TRUE)
Ca_per_emm_let

### STEP4 (Figure)
Ca_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Ca_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Na_per (%) ####
## STEP1 (meand and SD) Na_per
# BS Na_per (%)
aggregate(Vor_BS$Na_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Na_per, list(Vor_BS$Treatment), FUN=sd)
# TA Na_per (%)
aggregate(Vor_TA$Na_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Na_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Na_per
Na_per_mod <- lmer(Na_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Na_per_mod)

### STEP3 (contrasts emmeans) Na_per
Na_per_emm <- emmeans(Na_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Na_per_emm
emmeans(Na_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Na_per_emm_let <- multcomp::cld(object = Na_per_emm$emmeans,
                                
                                Letters = letters,no.readonly=TRUE)
Na_per_emm_let

### STEP4 (Figure)
Na_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Na_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### H_per (%) ####
## STEP1 (meand and SD) H_per
# BS H_per (%)
aggregate(Vor_BS$H_per, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$H_per, list(Vor_BS$Treatment), FUN=sd)
# TA H_per (%)
aggregate(Vor_TA$H_per, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$H_per, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) H_per
H_per_mod <- lmer(H_per ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(H_per_mod)

### STEP3 (contrasts emmeans) H_per
H_per_emm <- emmeans(H_per_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
H_per_emm
emmeans(H_per_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
H_per_emm_let <- multcomp::cld(object = H_per_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
H_per_emm_let

### STEP4 (Figure)
H_per_Fig <- Treats4 %>%
  ggplot(aes(Treatment, H_per, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### P_mg (%) ####
## STEP1 (meand and SD) P_mg
# BS P_mg (%)
aggregate(Vor_BS$P_mg, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$P_mg, list(Vor_BS$Treatment), FUN=sd)
# TA P_mg (%)
aggregate(Vor_TA$P_mg, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$P_mg, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) P_mg
P_mg_mod <- lmer(P_mg ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(P_mg_mod)

### STEP3 (contrasts emmeans) P_mg
P_mg_emm <- emmeans(P_mg_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
P_mg_emm
emmeans(P_mg_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
P_mg_emm_let <- multcomp::cld(object = P_mg_emm$emmeans,
                              
                              Letters = letters,no.readonly=TRUE)
P_mg_emm_let

### STEP4 (Figure)
P_mg_Fig <- Treats4 %>%
  ggplot(aes(Treatment, P_mg, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Al_mg (%) ####
## STEP1 (meand and SD) Al_mg
# BS Al_mg (%)
aggregate(Vor_BS$Al_mg, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Al_mg, list(Vor_BS$Treatment), FUN=sd)
# TA Al_mg (%)
aggregate(Vor_TA$Al_mg, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Al_mg, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Al_mg
Al_mg_mod <- lmer(Al_mg ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Al_mg_mod)

### STEP3 (contrasts emmeans) Al_mg
Al_mg_emm <- emmeans(Al_mg_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Al_mg_emm
emmeans(Al_mg_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Al_mg_emm_let <- multcomp::cld(object = Al_mg_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
Al_mg_emm_let

### STEP4 (Figure)
Al_mg_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Al_mg, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Mn_mg (%) ####
## STEP1 (meand and SD) Mn_mg
# BS Mn_mg (%)
aggregate(Vor_BS$Mn_mg, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Mn_mg, list(Vor_BS$Treatment), FUN=sd)
# TA Mn_mg (%)
aggregate(Vor_TA$Mn_mg, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Mn_mg, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Mn_mg
Mn_mg_mod <- lmer(Mn_mg ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Mn_mg_mod)

### STEP3 (contrasts emmeans) Mn_mg
Mn_mg_emm <- emmeans(Mn_mg_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Mn_mg_emm
emmeans(Mn_mg_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Mn_mg_emm_let <- multcomp::cld(object = Mn_mg_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
Mn_mg_emm_let

### STEP4 (Figure)
Mn_mg_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Mn_mg, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Fe_mg (%) ####
## STEP1 (meand and SD) Fe_mg
# BS Fe_mg (%)
aggregate(Vor_BS$Fe_mg, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$Fe_mg, list(Vor_BS$Treatment), FUN=sd)
# TA Fe_mg (%)
aggregate(Vor_TA$Fe_mg, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$Fe_mg, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) Fe_mg
Fe_mg_mod <- lmer(Fe_mg ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(Fe_mg_mod)

### STEP3 (contrasts emmeans) Fe_mg
Fe_mg_emm <- emmeans(Fe_mg_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
Fe_mg_emm
emmeans(Fe_mg_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
Fe_mg_emm_let <- multcomp::cld(object = Fe_mg_emm$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
Fe_mg_emm_let

### STEP4 (Figure)
Fe_mg_Fig <- Treats4 %>%
  ggplot(aes(Treatment, Fe_mg, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### S_mg (%) ####
## STEP1 (meand and SD) S_mg
# BS S_mg (%)
aggregate(Vor_BS$S_mg, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$S_mg, list(Vor_BS$Treatment), FUN=sd)
# TA S_mg (%)
aggregate(Vor_TA$S_mg, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$S_mg, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) S_mg
S_mg_mod <- lmer(S_mg ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(S_mg_mod)

### STEP3 (contrasts emmeans) S_mg
S_mg_emm <- emmeans(S_mg_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
S_mg_emm
emmeans(S_mg_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
S_mg_emm_let <- multcomp::cld(object = S_mg_emm$emmeans,
                              
                              Letters = letters,no.readonly=TRUE)
S_mg_emm_let

### STEP4 (Figure)
S_mg_Fig <- Treats4 %>%
  ggplot(aes(Treatment, S_mg, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### CEC (%) ####
## STEP1 (meand and SD) CEC
# BS CEC (%)
aggregate(Vor_BS$CEC, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$CEC, list(Vor_BS$Treatment), FUN=sd)
# TA CEC (%)
aggregate(Vor_TA$CEC, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$CEC, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) CEC
CEC_mod <- lmer(CEC ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(CEC_mod)

### STEP3 (contrasts emmeans) CEC
CEC_emm <- emmeans(CEC_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
CEC_emm
emmeans(CEC_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
CEC_emm_let <- multcomp::cld(object = CEC_emm$emmeans,
                             
                             Letters = letters,no.readonly=TRUE)
CEC_emm_let

### STEP4 (Figure)
CEC_Fig <- Treats4 %>%
  ggplot(aes(Treatment, CEC, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### C_N (%) ####
## STEP1 (meand and SD) C_N
# BS C_N (%)
aggregate(Vor_BS$C_N, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$C_N, list(Vor_BS$Treatment), FUN=sd)
# TA C_N (%)
aggregate(Vor_TA$C_N, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$C_N, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) C_N
C_N_mod <- lmer(C_N ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(C_N_mod)

### STEP3 (contrasts emmeans) C_N
C_N_emm <- emmeans(C_N_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
C_N_emm
emmeans(C_N_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
C_N_emm_let <- multcomp::cld(object = C_N_emm$emmeans,
                             
                             Letters = letters,no.readonly=TRUE)
C_N_emm_let

### STEP4 (Figure)
C_N_Fig <- Treats4 %>%
  ggplot(aes(Treatment, C_N, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### N_P (%) ####
## STEP1 (meand and SD) N_P
# BS N_P (%)
aggregate(Vor_BS$N_P, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$N_P, list(Vor_BS$Treatment), FUN=sd)
# TA N_P (%)
aggregate(Vor_TA$N_P, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$N_P, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) N_P
N_P_mod <- lmer(N_P ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(N_P_mod)

### STEP3 (contrasts emmeans) N_P
N_P_emm <- emmeans(N_P_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
N_P_emm
emmeans(N_P_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#N_P_emm2 <- emmeans(N_P_mod, ~ Canopy:Treatment, data = Treats4)
#pairs(N_P_emm2, simple = list("Canopy", "Treatment")) #Gives different results for control (BS-TA) P = 0.0501! and in Ti (BS-TA) P = 0.0261

#Letters
N_P_emm_let <- multcomp::cld(object = N_P_emm$emmeans,
                             
                             Letters = letters,no.readonly=TRUE)
N_P_emm_let

### STEP4 (Figure)
N_P_Fig <- Treats4 %>%
  ggplot(aes(Treatment, N_P, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### P_Al (%) ####
## STEP1 (meand and SD) P_Al
# BS P_Al (%)
aggregate(Vor_BS$P_Al, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$P_Al, list(Vor_BS$Treatment), FUN=sd)
# TA P_Al (%)
aggregate(Vor_TA$P_Al, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$P_Al, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) P_Al
P_Al_mod <- lmer(P_Al ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(P_Al_mod)

### STEP3 (contrasts emmeans) P_Al
P_Al_emm <- emmeans(P_Al_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
P_Al_emm
emmeans(P_Al_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
P_Al_emm_let <- multcomp::cld(object = P_Al_emm$emmeans,
                              
                              Letters = letters,no.readonly=TRUE)
P_Al_emm_let

### STEP4 (Figure)
P_Al_Fig <- Treats4 %>%
  ggplot(aes(Treatment, P_Al, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### P_Ca (%) ####
## STEP1 (meand and SD) P_Ca
# BS P_Ca (%)
aggregate(Vor_BS$P_Ca, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$P_Ca, list(Vor_BS$Treatment), FUN=sd)
# TA P_Ca (%)
aggregate(Vor_TA$P_Ca, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$P_Ca, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) P_Ca
P_Ca_mod <- lmer(P_Ca ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(P_Ca_mod)

### STEP3 (contrasts emmeans) P_Ca
P_Ca_emm <- emmeans(P_Ca_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
P_Ca_emm
emmeans(P_Ca_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
P_Ca_emm_let <- multcomp::cld(object = P_Ca_emm$emmeans,
                              
                              Letters = letters,no.readonly=TRUE)
P_Ca_emm_let

### STEP4 (Figure)
P_Ca_Fig <- Treats4 %>%
  ggplot(aes(Treatment, P_Ca, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### pH_H2O (%) ####
## STEP1 (meand and SD) pH_H2O
# BS pH_H2O (%)
aggregate(Vor_BS$pH_H2O, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$pH_H2O, list(Vor_BS$Treatment), FUN=sd)
# TA pH_H2O (%)
aggregate(Vor_TA$pH_H2O, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$pH_H2O, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) pH_H2O
pH_H2O_mod <- lmer(pH_H2O ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(pH_H2O_mod)

### STEP3 (contrasts emmeans) pH_H2O
pH_H2O_emm <- emmeans(pH_H2O_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
pH_H2O_emm
emmeans(pH_H2O_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
pH_H2O_emm_let <- multcomp::cld(object = pH_H2O_emm$emmeans,
                                
                                Letters = letters,no.readonly=TRUE)
pH_H2O_emm_let

### STEP4 (Figure)
pH_H2O_Fig <- Treats4 %>%
  ggplot(aes(Treatment, pH_H2O, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### pH_buf (%) ####
## STEP1 (meand and SD) pH_buf
# BS pH_buf (%)
aggregate(Vor_BS$pH_buf, list(Vor_BS$Treatment), FUN=mean)
aggregate(Vor_BS$pH_buf, list(Vor_BS$Treatment), FUN=sd)
# TA pH_buf (%)
aggregate(Vor_TA$pH_buf, list(Vor_TA$Treatment), FUN=mean)
aggregate(Vor_TA$pH_buf, list(Vor_TA$Treatment), FUN=sd)

### STEP2 (aov(lme model)) pH_buf
pH_buf_mod <- lmer(pH_buf ~ Canopy*Treatment + (1|Site/Block), data = Treats4)
anova(pH_buf_mod)

### STEP3 (contrasts emmeans) pH_buf
pH_buf_emm <- emmeans(pH_buf_mod, specs = trt.vs.ctrl ~ Treatment|Canopy, type = "response")
pH_buf_emm
emmeans(pH_buf_mod, specs = trt.vs.ctrl ~ Canopy*Treatment, type = "response") #To check BS-TA for C

#Letters
pH_buf_emm_let <- multcomp::cld(object = pH_buf_emm$emmeans,
                                
                                Letters = letters,no.readonly=TRUE)
pH_buf_emm_let

### STEP4 (Figure)
pH_buf_Fig <- Treats4 %>%
  ggplot(aes(Treatment, pH_buf, group=Canopy, color=Treatment)) +
  stat_summary(geom="pointrange", fun.data = mean_cl_boot, size=1) + 
  facet_grid(cols=vars(Canopy)) + 
  scale_colour_manual(values = Farben)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

## All Figs together
ggarrange(C_per_Fig,N_per_Fig,K_per_Fig,Mg_per_Fig,Ca_per_Fig,Na_per_Fig,H_per_Fig,P_mg_Fig,Al_mg_Fig,Mn_mg_Fig,Fe_mg_Fig,S_mg_Fig,CEC_Fig,C_N_Fig,N_P_Fig,P_Al_Fig,P_Ca_Fig,pH_H2O_Fig,pH_buf_Fig,Temp_Fig,Hum_Fig,Light_Fig,OverDen_Fig, ncol = 4, nrow = 6, common.legend = TRUE, labels="auto")
## All Figs with sig. differences together
#ggarrange(C_per_Fig,N_per_Fig,Ca_per_Fig,Na_per_Fig,H_per_Fig,P_mg_Fig,Al_mg_Fig,Mn_mg_Fig,Fe_mg_Fig,S_mg_Fig,C_N_Fig,N_P_Fig,P_Al_Fig,pH_H2O_Fig,pH_buf_Fig, ncol = 5, nrow = 3, common.legend = TRUE, labels="auto")


ggarrange(C_per_Fig,Na_per_Fig,P_mg_Fig,Al_mg_Fig,Fe_mg_Fig,S_mg_Fig,C_N_Fig,N_P_Fig,P_Al_Fig, ncol = 3, nrow = 3, common.legend = TRUE, labels="auto")
