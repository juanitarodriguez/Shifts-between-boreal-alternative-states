## Script information
### Script of the article entitled: Drivers of contrasting boreal understory vegetation in coniferous and broadleaf deciduous alternative states
### Article submitted to Ecological Monographs
### Authors of the article: Juanita C. Rodríguez-Rodríguez, Nicole J. Fenton, Steven W. Kembel, Evick Mestre, Mélanie Jean, and Yves Bergeron
### Script developed on 2021 by: Juanita C. Rodríguez-Rodríguez (e-mail: juanitacarolina.rodriguezrodriguez@uqat.ca), for my PhD in Environmental Sciences. A part of the script was developed during the statistical cours of Pierre Legendre.

## Load packages
library(dplyr) #For droplevels()
library(tidyr) #For spread()
library(stats) #For filter(), For PCA: prcomp()
library(vegan) #For PCA: decostand(), PCA: rda()
library(BiodiversityR) #For ordiequilibriumcircle() in PCA
library(ggplot2) #For Figures
library(nlme) # (for lme)
library(devtools)
install_github("sdray/adespatial")
library(adespatial) # For beta.div (TBI)
library(ade4) # For s.value (Beta div plot) and is.euclid()
library(factoextra) # For fviz_pca_biplot
library(cowplot) # For ggdraw
library(ggpubr) # For ggarrange
library(fastDummies)
library(ape) # For pcoa()
library(emmeans) #For lsmeans()
library(ggbeeswarm) #For geom_quasirandom()
library(lmerTest) #To include p values in lmer

## Load all data and prepare tables
DB_All <- read.csv2("DB_AllData2.csv")
DB_All2 <- DB_All

#### Prepare a species table with all variables
DB_All2$Species <- NULL # Leaving just species Abbreviations (Abbr)
DB_All2$Func_gr <- NULL # Not for now

Spp_All <- spread(DB_All2, Abbr, Coverage, fill = DB_All2$Coverage)
which(is.na(DB_All2$Coverage)) # I checked it because after spread, I get a Warning message: In if (!is.na(fill)) { : the condition has length > 1 and only the first element will be used
head(Spp_All) # Species table with all variables

#### Prepare tables from Spp_All
Env <- select(Spp_All, 1:5)
Spp <- select(Spp_All, 6:70)

Env$Year <- as.factor(Env$Year)
Env$Block <- as.factor(Env$Block)

#### Organize levels of the factor Treatment
Env$Treatment <- factor(Env$Treatment, levels = c("Li","Nu","1F","2F","To","Ti","C"))

#### Separate matrices by forest type (Canopy) to do alpha analysis (with all years included)
Canopy <- split(Spp_All, Spp_All$Canopy)
BS <- Canopy$'BS'
TA <- Canopy$'TA'

#	Alpha diversity of understory communities ####
## Figure 4
### Alpha diversity (Shannon-Weiner index) of understory vegetation in 2018 among treatments in each forest type (Black spruce - BS and Trembling aspen - TA), over time (2013-2018) (n=9). Treatments correspond to light (Li, in yellow), nutrients (Nu, in purple), single amount of leaves (1F, in orange) and double amount of leaves (2F, in red), transplants-out (To, in blue), Transplants-in (Ti, in green) and control conditions (C, in grey).

### Diversity plot with Shannon index - facet_grid with both Year and Canopy
####Ecosystem approach: "gold","mediumpurple1","tan1","lightcoral","paleturquoise2","chartreuse","azure3"
####Community approach: "goldenrod","darkorchid3","orangered","red3","lightseagreen","forestgreen","azure4"

div.plot2 <- ggplot(Spp, aes(x = Env$Treatment, y = diversity(Spp, index = "shannon"), fill=Env$Treatment)) +
  geom_boxplot()+
  stat_summary(fun = mean, geom="point", pch=21, col = "black", bg = "white", size = 2.8)+
  labs(x="Treatments", y = "Alpha diversity (Shannon index)")+
  facet_grid(cols=vars(Env$Year), vars(Env$Canopy)) +
  scale_fill_manual(values=alpha(c("gold","mediumpurple1","tan1","lightcoral","paleturquoise2","chartreuse","azure3")), name = "Treatments") + theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))
div.plot2

## ANOVA of alpha diversity
## Table S2. 
### Differences in α-diversity (Shannon-Weiner index) based on the ANOVA of the linear mixed effect model of understory species abundance to evaluate the factors of Year, Forest type, Treatment and all interactions, using Sites and Blocks as nested random factors.
set.seed(1989)
div.aov5 <- lmer(Div.spp ~ Year*Canopy*Treatment + (1|Site/Block), data = Env)
summary(div.aov5)
anova(div.aov5)
plot(div.aov5)

### Comparisons per Year and Canopy, between C and each Treatment
div.emmeans <- emmeans(div.aov5, specs = trt.vs.ctrl ~ Treatment|Year*Canopy, type = "response")
div.emmeans
#### Obtain the letters for those comparisons
div.mod_means <- multcomp::cld(object = div.emmeans$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
div.mod_means

## Post-hoc test
## Table S3. 
### Post-hoc emmeans contrasts of the lineral mixed effect model (Table S2) of the α-diversity (Shannon index) of plant understory species in each forest type (Black spruce – BS, Trembling aspen - TA), comparing each treatment to the control. Significant differences correspond to P < 0.0001 ***, <0.001 **, < 0.05 *, denoted in letters for Figure 4.
lsmeans(div.aov2, pairwise~Year*Canopy*Treatment, adjust="tukey")
t.test(diversity(Spp) ~ Env$Canopy)

### Re-Order levels of the factor  so I do  comparisons of C versus all
Env$Treatment <- factor(Env$Treatment, levels=c("C","Li","Nu","1F","2F","To","Ti"))

# Beta diversity ####
#### Make matrix of 2018 data
Spp_Year <- split(Spp_All, Spp_All$Year)
Spp_2018 <- Spp_Year$`2018`
Spp_18 <- select(Spp_2018, 6:70)
Env_18 <- select(Spp_2018, 1:7)

# Non-directional approach: Partitioning beta diversity ####
## Compute beta diversity with function beta.div() of adepsatial
#### B.div for each objective separately and each canopy

#### Separate by Canopy
Spp_18Can <- split(Spp_2018, Spp_2018$Canopy)
Spp_18BS <- Spp_18Can$`BS`
Spp_18TA <- Spp_18Can$`TA`

####Separate by Treatment - Objective approach
####BS
Spp_18BS_Tbio <- filter(Spp_18BS, Treatment %in% c("C","1F","2F","Li","Nu"))
Spp_18BS_Tss <- filter(Spp_18BS, Treatment %in% c("C","Ti","To"))
####TA
Spp_18TA_Tbio <- filter(Spp_18TA, Treatment %in% c("C","1F","2F","Li","Nu"))
Spp_18TA_Tss <- filter(Spp_18TA, Treatment %in% c("C","Ti","To"))

####Divide into Spp and Env tables
####BS, Tbio
BS_Tbio_spp <- select(Spp_18BS_Tbio, 6:70)
BS_Tbio_env <- select(Spp_18BS_Tbio, 1:5)
####BS, Tss
BS_Tss_spp <- select(Spp_18BS_Tss, 6:70)
BS_Tss_env <- select(Spp_18BS_Tss, 1:5)
####TA, Tbio
TA_Tbio_spp <- select(Spp_18TA_Tbio, 6:70)
TA_Tbio_env <- select(Spp_18TA_Tbio, 1:5)
####BS, Tss
TA_Tss_spp <- select(Spp_18TA_Tss, 6:70)
TA_Tss_env <- select(Spp_18TA_Tss, 1:5)

Tbio_cols <- c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid")
Tss_cols <- c("antiquewhite4","chartreuse3","lightseagreen")

## Beta-diversity for each approach
## Figure 5
### Partitioning of β-diversity into local contributions of single plots (LCBD = Local Contributions to Beta Diversity) for treatments in each forest type (Black spruce and Trembling aspen) of (a, b) the ecosystem approach corresponding to Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red), and for (c, d) the community approach corresponding to Transplants-out (To, in blue), Transplants-in (Ti, in green), as well as Control conditions (C, in grey) for both approaches. Treatments in black spruce forests in light colors and in trembling aspen forests in dark colors. Tables below each figure present the total sum-of-squares (SStotal), total β-diversity (BDtotal) and the contributions (decreasing order) of individual species (SCBD = Species Contributions to Beta Diversity), for each forest type in each approach. Data correspond to 2018, n = 9, from the combination of Site (A, B, C) and Blocks (1, 2, 3) in each forest type.

####BS, Tbio = BS_Tbio_spp, BS_Tbio_env
res_BS_Tbio <- beta.div(BS_Tbio_spp, method="hellinger", nperm = 999) # Hell transformation of data
summary(res_BS_Tbio)
res_BS_Tbio$beta
res_BS_Tbio$LCBD
which(res_BS_Tbio$p.LCBD <= 0.05) # Which are the significant LCBD indices?
res_BS_Tbio$SCBD
res_BS_Tbio$SCBD[res_BS_Tbio$SCBD >= mean(res_BS_Tbio$SCBD)]
res_BS_Tbio$p.adj

#### Figure 5a
LCBD.plot_BST_bio <- ggplot(BS_Tbio_spp, aes(x = BS_Tbio_env$Treatment, y = res_BS_Tbio$LCBD), fill=res_BS_Tbio$LCBD) +
  geom_boxplot(fill=Tbio_cols) +
  labs(x="Treatments", y = "LCBD") +
  facet_grid(cols=vars(BS_Tbio_env$Canopy)) + 
  ylim(0, 0.08) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

LCBD.plot_BST_bio_aov <- aov(res_BS_Tbio$LCBD ~ BS_Tbio_env$Treatment / BS_Tbio_env$Site*BS_Tbio_env$Block)
summary(LCBD.plot_BST_bio_aov)


#### Example:https://www.flutterbys.com.au/stats/tut/tut9.2a.html
####data.nest.agg = data.nest %>% group_by(A, Sites) %>% summarize(y = mean(y))
ejemp <- how(nperm = 9999, 
             blocks = BS_Tbio_env$Site,  
             plots = Plots(BS_Tbio_env$Block), 
             Within(type = "free"))

ejemp_ran <- adonis2(res_BS_Tbio$LCBD ~ BS_Tbio_env$Treatment,
                     permutations = ejemp,
                     method = "euc",
                     by = "term")
ejemp_ran

####TA, Tbio = TA_Tbio_spp, TA_Tbio_env
res_TA_Tbio <- beta.div(TA_Tbio_spp, method="hellinger", nperm = 999) # Hell transformation of data
summary(res_TA_Tbio)
res_TA_Tbio$beta
res_TA_Tbio$LCBD
which(res_TA_Tbio$p.LCBD <= 0.05)
res_TA_Tbio$SCBD
res_TA_Tbio$SCBD[res_TA_Tbio$SCBD >= mean(res_TA_Tbio$SCBD)]
res_TA_Tbio$p.adj

#### Figure 5b
LCBD.plot_TAT_bio <- ggplot(TA_Tbio_spp, aes(x = TA_Tbio_env$Treatment, y = res_TA_Tbio$LCBD), fill=res_TA_Tbio$LCBD) +
  geom_boxplot(fill=Tbio_cols) +
  labs(x="Treatments", y = "LCBD") +
  facet_grid(cols=vars(TA_Tbio_env$Canopy)) + 
  ylim(0, 0.08) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

LCBD.plot_BST_bio_aov <- aov(res_BS_Tbio$LCBD ~ BS_Tbio_env$Treatment / BS_Tbio_env$Site*BS_Tbio_env$Block)
summary(LCBD.plot_BST_bio_aov)

ejemp_TA <- how(nperm = 9999, 
                blocks = TA_Tbio_env$Site,  
                plots = Plots(TA_Tbio_env$Block), 
                Within(type = "free"))

ejemp_ran_TA <- adonis2(res_TA_Tbio$LCBD ~ TA_Tbio_env$Treatment,
                        permutations = ejemp,
                        method = "euc",
                        by = "term")
ejemp_ran_TA


####BS, Tss = BS_Tss_spp, BS_Tss_env
res_BS_Tss <- beta.div(BS_Tss_spp, method="hellinger", nperm = 999) # Hell transformation of data
summary(res_BS_Tss)
res_BS_Tss$beta
res_BS_Tss$LCBD
which(res_BS_Tss$p.LCBD <= 0.05)
res_BS_Tss$SCBD
res_BS_Tss$SCBD[res_BS_Tss$SCBD >= mean(res_BS_Tss$SCBD)]
res_BS_Tss$p.adj

#### Figure 5c
LCBD.plot_BST_ss <- ggplot(BS_Tss_spp, aes(x = BS_Tss_env$Treatment, y = res_BS_Tss$LCBD), fill=res_BS_Tss$LCBD) +
  geom_boxplot(fill=Tss_cols)+
  labs(x="Treatments", y = "LCBD") +
  facet_grid(cols=vars(BS_Tss_env$Canopy)) + 
  ylim(0, 0.08) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#TA, Tss = TA_Tss_spp, TA_Tss_env
res_TA_Tss <- beta.div(TA_Tss_spp, method="hellinger", nperm = 999) # Hell transformation of data
summary(res_TA_Tss)
res_TA_Tss$beta
res_TA_Tss$LCBD
which(res_TA_Tss$p.LCBD <= 0.05)
res_TA_Tss$SCBD
res_TA_Tss$SCBD[res_TA_Tss$SCBD >= mean(res_TA_Tss$SCBD)]
res_TA_Tss$p.adj

#### Figure 5d
LCBD.plot_TAT_ss <- ggplot(TA_Tss_spp, aes(x = TA_Tss_env$Treatment, y = res_TA_Tss$LCBD), fill=res_TA_Tss$LCBD) +
  geom_boxplot(fill=Tss_cols)+
  labs(x="Treatments", y = "LCBD") +
  facet_grid(cols=vars(TA_Tss_env$Canopy)) + 
  ylim(0, 0.08) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Plot all 3 Figures together. Save PDF as 5x7 landscape
#### Figure 5a-d
#### Plot all four Figures together. Save PDF as 5x7 landscape
(LCBD_Allplots <- ggarrange(LCBD.plot_BST_bio, LCBD.plot_TAT_bio, LCBD.plot_BST_ss, LCBD.plot_TAT_ss + rremove("x.text"), 
                            labels = c("a)", "b)", "c)", "d)"),
                            ncol = 2, nrow = 2))

(LCBD_Allplots_Tbio <- ggarrange(LCBD.plot_BST_bio, LCBD.plot_TAT_bio + rremove("x.text"), 
                                 labels = c("a)", "b)"),
                                 ncol = 2, nrow = 2))

(LCBD_Allplots_Tss <- ggarrange(LCBD.plot_BST_ss, LCBD.plot_TAT_ss + rremove("x.text"), 
                                labels = c("c)", "d)"),
                                ncol = 2, nrow = 2))


# Species turnover: Temporal analysis of multivariate data (TBI) ####
####Note: As adespatial package did not worked, I had to upload the whole script of TBI from Legendre script.

### TBI analysis with all results together for each forest type, between years 2013 and 2018
#### Preparing data tables
DBTbio <- DB_All

#### Filter by canopy dominance (BS or TA)
#### BS Matrix
BS.mat <- DBTbio %>%
  filter(Canopy %in% c("BS")) %>% 
  droplevels()
#### TA Matrix
TA.mat <- DBTbio %>%
  filter(Canopy %in% c("TA")) %>% 
  droplevels()

#### Make two matrices, one with 2013 data and the second with 2018 data. I'll get two matrices for each canopy composition.
#### Matrices of BS, filtering by year.
#### 2013 Matrix of BS
BS.2013 <- BS.mat %>%
  filter(Year %in% c("2013")) %>% 
  droplevels()
#### 2018 Matrix of BS
BS.2018 <- BS.mat %>%
  filter(Year %in% c("2018")) %>% 
  droplevels()
#### Matrices of TA, filtering by year.
#### 2013 Matrix of TA
TA.2013 <- TA.mat %>%
  filter(Year %in% c("2013")) %>% 
  droplevels()
#### 2018 Matrix of TA
TA.2018 <- TA.mat %>%
  filter(Year %in% c("2018")) %>% 
  droplevels()
#### Now, I already have my four matrices of 2013 and 2018 for BS and 2013 and 2018 for TA. 

#### Join columns of Site, Block and Treatment
#### For BS, 2013:
Mat1.BS <- unite(BS.2013, TP,c("Site","Block","Treatment"),  sep = ".", remove = TRUE)
Mat1.BS$Year <- NULL
Mat1.BS$Canopy <- NULL
Mat1.BS$Species <- NULL
Mat1.BS$Func_gr <- NULL
#### For BS, 2018:
Mat2.BS <- unite(BS.2018, TP,c("Site", "Block","Treatment"),  sep = ".", remove = TRUE)
Mat2.BS$Year <- NULL
Mat2.BS$Canopy <- NULL
Mat2.BS$Species <- NULL
Mat2.BS$Func_gr <- NULL
#### For TA, 2013:
Mat1.TA <- unite(TA.2013, TP,c("Site", "Block","Treatment"),  sep = ".", remove = TRUE)
Mat1.TA$Year <- NULL
Mat1.TA$Canopy <- NULL
Mat1.TA$Species <- NULL
Mat1.TA$Func_gr <- NULL
#### For TA, 2018:
Mat2.TA <- unite(TA.2018, TP,c("Site", "Block","Treatment"),  sep = ".", remove = TRUE)
Mat2.TA$Year <- NULL
Mat2.TA$Canopy <- NULL
Mat2.TA$Species <- NULL
Mat2.TA$Func_gr <- NULL
#### Now, I have all matrices (Mat1.BS, Mat2.BS, Mat1.TA, Mat2.TA) that I need for the analysis (filtered by Year, Canopy, with Abbr and with one single column with TP). However, I need a "species table", where the first column is my "sites" <- and the other columns are each Func_gr, filled by the Coverage in each cell.

####BS
Mat1.BS.spp <- Mat1.BS %>%
  pivot_wider(names_from = Abbr, values_from = Coverage)
Mat2.BS.spp <- Mat2.BS %>%
  pivot_wider(names_from = Abbr, values_from = Coverage)
####TA
Mat1.TA.spp <- Mat1.TA %>%
  pivot_wider(names_from = Abbr, values_from = Coverage)
Mat2.TA.spp <- Mat2.TA %>%
  pivot_wider(names_from = Abbr, values_from = Coverage)
#### Now I have spp tables ((Mat1.BS, Mat2.BS, Mat1.TA, Mat2.TA).spp) as I need.

#### I take a new table with the Treatment-plot names, because I need to eliminate it to do the analysis (must be numeric). I do the same for all the tables but they are the same, so I take just one in TP at the end.
TP.BS.Mat1 <- as.matrix(Mat1.BS.spp$TP)
TP.BS.Mat2 <- as.matrix(Mat2.BS.spp$TP)
TP.TA.Mat1 <- as.matrix(Mat1.TA.spp$TP)
TP.TA.Mat2 <- as.matrix(Mat2.TA.spp$TP)
(TP <- as.matrix(Mat1.BS.spp$TP))

(TP2 <- Mat1.BS.spp$TP)

####  Then, I need to add an ID to each row and then eliminate the TP columns of each DB that I'll use in the analysis.
ID <- data.frame("ID" = c(1:63)) #I create an ID with 9 numbers

#### I add the ID to all tables.
BS1 <- Mat1.BS.spp
BS1$TP <- NULL

BS2 <- Mat2.BS.spp
BS2$TP <- NULL

TA1 <- Mat1.TA.spp
TA1$TP <- NULL

TA2 <- Mat2.TA.spp
TA2$TP <- NULL

#### Now, I got finally all the tables for the analysis: BS1.spe, BS2.spe, TA1, TA2

#### TBI Analysis per forest type
BS1 <- as.data.frame(BS1)
BS2 <- as.data.frame(BS2)
TA1 <- as.data.frame(TA1)
TA2 <- as.data.frame(TA2)

# TBI - BS
TBI.BS <- TBI(BS1,BS2,method="%diff", pa.tr=FALSE, nperm=9999, BCD=TRUE, replace=FALSE, test.BC=TRUE, test.t.perm=TRUE, save.BC=TRUE, seed.=NULL, clock=FALSE)
TBI.BS # Non Sig - change: 0.3238 (P>0.05)
summary(TBI.BS)
BS.BCD.mat <- TBI.BS$BCD.mat

# TBI - TA
TBI.TA <- TBI(TA1,TA2,method="%diff", pa.tr=FALSE, nperm=9999, BCD=TRUE, replace=FALSE, test.BC=TRUE, test.t.perm=TRUE, save.BC=TRUE, seed.=NULL, clock=FALSE)
TBI.TA # Sig + change: 1e-04 * (P<0.05)
summary(TBI.TA)
TA.BCD.mat <- TBI.TA$BCD.mat

## Figure 6. 
### Analysis of Temporal Beta-diversity Index (TBI) of understory communities in black spruce (left panel) and trembling aspen (right panel) stands. Figures correspond to B-C plots of understory vegetation comparing all treatments from 2013 to 2018, where 63 plots per forest type were plotted using losses (B/den statistics) and gains (C/den statistics) computed from the abundance data of understory species. These sites correspond to the combination of sampling design of Sites (A,B,C), Blocks (1, 2, 3) and seven Treatments: Control conditions (C, in grey), to treatments of the ecosystem approach: Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red), and of the community approach: Transplants-out (To, in blue), Transplants-in (Ti, in green). Treatments in black spruce forests in light colors and in trembling aspen forests in dark colors. Notice also that axes are different to allow clearer data visualisation. The green line (with a slope of 1) is drawn to through the origin and it represent the theoretical positions of a net-zero species gain/loss scenario. The position of green line respective to red line indicates an average of species losses (green line above red line, as in left panel) or species gains (green line below red line, as in right panel) across the sites. Distinctive symbols are used for the sites dominated by gains (squares) and by losses (circles). Symbols sizes represent the values of the D = (B+C) statistics: larger points found in the upper-right corner of each plot represent a higher temporal β-diversity than the smaller points in the lower-left corner.
#### B-C Plots
#####plot.TBI # Appeler la fonction avant! From: https://github.com/sdray/adespatial/blob/master/R/plot.TBI.R
plot.TBI(TBI.BS)
plot.TBI(TBI.TA)

#### Save both Figs together as 12x6 portrait PDF
par(mfrow = c(1, 2))

#### BS
plot.TBI(TBI.BS, type = "BC", 
         s.names = TP2, pch.loss = 21, pch.gain = 22, 
         cex.names = 0.5, col.rim = "black", col.bg = "dodgerblue3", 
         cex.symb = 3, diam = TRUE, 
         main = "Black spruce forests", cex.main = 1, cex.lab = 1, 
         xlim = NULL, ylim = NULL, silent = TRUE)
#### TA
plot.TBI(TBI.TA, type = "BC", 
         s.names = TP2, pch.loss = 21, pch.gain = 22, 
         cex.names = 0.5, col.rim = "black", col.bg = "hotpink", 
         cex.symb = 3, diam = TRUE, 
         main = "Trembling aspen forests", cex.main = 1, cex.lab = 1, 
         xlim = NULL, ylim = NULL, silent = TRUE)

par(mfrow = c(1, 1))

# Ordinations for each objective (ecosystem and community approaches) ####

## Figure 7.
### Principal Component Analysis (PCA) for the Objective 1 (a - Ecosystem approach) and the Objective 2 (b - Community approach) illustrating community change between 2013 and 2018. Points correspond to the centroids per treatment (means of 9 replicate plots from nested design of 3 Sites x 3 Blocks) of the Hellinger-transformed abundance of understory species from the beginning (2013, circles) and the end of the experiment (2018, triangles), united by lines indicating change over time. Solid lines correspond to black spruce stands (BS) and dotted lines correspond to trembling aspen stands (TA). Colors correspond to treatments of the ecosystem approach: Single-litter (1F, in orange), Double-litter (2F, in red), Control (C, in grey), Light (Li, in yellow) and Nutrients (Nu, in purple) and of treatments of the community approach: Control (C, in grey), Transplants-in (Ti, in green) and Transplants-out (To, in blue). Treatments in black spruce forests in light colors and in trembling aspen forests in dark colors. Species ordination is in the middle of the panel, with only the main species to allow a clear visualization. Notice that the scales are different.

## Table S4.
###PERMANOVA based on the Hellinger-transformed data (Bray-Curtis distance) of understory vegetation abundance for the ecosystem and community approaches to test the variables of Year, Forest type, Treatment and all possible interactions, using Site and Blocks as nested random factors.

### Prepare a species table with all variables
DB_All$Species <- NULL
DB_All$Func_gr <- NULL

Spp_All <- spread(DB_All, Abbr, Coverage, fill = DB_All$Coverage)
Spp_All$Y.C.T <- paste0(Spp_All$Year,"_", Spp_All$Canopy, "_", Spp_All$Treatment)
Spp_All <- Spp_All[, c(1:5, 71, 6:70)] # To order columns

#### Ecosytem approach (Tbio) - Analysis with both Canopies
#### Select treatments of Tbio (Exclude treatments Ti and To)
Spp.Tbio <- Spp_All %>%
  filter(Treatment %in% c("C", "Li", "Nu", "1F", "2F")) %>% 
  droplevels()

## Separate years
Year <- split(Spp.Tbio, Spp.Tbio$Year)
Y2013 <- Year$`2013`
Y2018 <- Year$`2018`

BSTA.Year <- rbind(Y2013, Y2018)
BSTA.sxs <- BSTA.Year[,6:71] #From column Y.C.T to the end of spp columns
BSTA.env <- BSTA.Year[,1:5] #All sampling info

# Combining the data of 2013 and 2018 for each stand
BSTA <- rbind(Y2013, Y2018)
BSTA.PCA <- BSTA[,6:71] #From column Y.C.T to the end of spp columns

#Make a matrix with Site and Block columns to add random factors in the following PERMANOVAs
BSTA_SB <- rbind(Y2013, Y2018)
BSTA_SB <- BSTA_SB[,2:3]

#### Centroid for Cntrl
Data.BSTA<-droplevels(subset(BSTA.PCA))
Data0.BSTA<- Data.BSTA[, colSums(Data.BSTA != 0) > 0] # To eliminate empty (0) columns

# Data transformation and PCA
PCA.crd.dec <- decostand(Data0.BSTA[,-1], "hellinger")
PCA.BSTA <- prcomp(PCA.crd.dec, center = TRUE, scale. = FALSE)
PCA.crd.dec <- cbind(Data0.BSTA$Y.C.T, PCA.crd.dec) #Adding column with names
colnames(PCA.crd.dec)
names(PCA.crd.dec)[names(PCA.crd.dec) == "Data0.BSTA$Y.C.T"] <- "Y.C.T"
PCAD.BSTA <- as.data.frame(PCA.BSTA$x)

PCAD.BSTA$Year <- with(PCA.crd.dec, substr(Y.C.T, 1,4))
PCAD.BSTA$Canopy <- with(PCA.crd.dec, substr(Y.C.T, 6,7))
PCAD.BSTA$Treatment <- with(PCA.crd.dec, substr(Y.C.T, 9,nchar(as.character(Y.C.T))))

# Mean for each treatment and year
pca.centroids.BSTA <- aggregate(PCAD.BSTA[,1:2], list(Year = PCAD.BSTA$Year, Canopy=PCAD.BSTA$Canopy, Treatment=PCAD.BSTA$Treatment), mean)
pca.centroids.BSTA$LAB<-with(pca.centroids.BSTA, ifelse(Year==2018, "18", "13"))

#### Figure 7 (part), Save as 5x5 Landscape
Main.BSTA<-ggplot(data=pca.centroids.BSTA, aes(PC1, PC2, col=Treatment))+
  theme_bw()+
  geom_hline(yintercept = 0, lty=2, col="gray")+
  geom_vline(xintercept = 0, lty=2, col="gray")+
  geom_point(aes(shape = Year, color=Treatment),size = 2)+
  scale_color_manual(values=c('darkorange1', 'red3', 'gray53', 'gold', 'darkmagenta'))+
  geom_line(aes(PC1, PC2,  group=interaction(Canopy, Treatment), linetype=Canopy)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom",
        plot.title = element_text(size=12,hjust = 0.5,face="bold"),
        axis.title=element_text(size=11),text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5, size=8),
        axis.text.y = element_text(angle=0,hjust=0.5,size=8))+ 
  labs(x= "PC1 (37.8%)", y ="PC2 (9.7%)", title="Ecosystem approach")

summary(PCA.BSTA)$importance[2,1] # PCA1: 0.37872
summary(PCA.BSTA)$importance[2,2] # PCA2: 0.09751
summary(PCA.BSTA)$importance[2,3] # PCA2: 0.07542
summary(PCA.BSTA)$importance[2,4] # PCA2: 0.05592

# Screen plot and broken stick model to test the significance of PCA axes
screeplot(PCA.BSTA, bstick = TRUE, npcs = length(PCA.BSTA$center))
bstick(PCA.BSTA) # I get four axes

### Subplot of species
scores(PCA.BSTA$center) #To get the scores of species and select the most important ones

#### Figure 7 (part-Species)
#Ancient spp selected: "AUR", "ARN", "CLB", "CON", "EQP", "GAH", "LIB", "LYA", "MAC", "PLS", "POA", "PTC", "LEG", "RIG", "RUP", "SPS", "VAM", "VIE", "VIS"
#Selecteed species with the highest scores (0.39-0.04): PLS,CON,GAH,PTC,RUP,LYA,MAC,VIS,LIB,CLB,VAM,VIE,ARN
Sub.BSTA<-fviz_pca_biplot(PCA.BSTA, label = "var",labelsize = 3, habillage="none",title="",
                          xlab="", ylab="", addEllipses=FALSE, geom.ind = "",col.var = "black", select.var=list(name=c('PLS','CON','GAH','PTC','RUP','LYA','MAC','VIS','LIB','CLB','VAM','VIE','ARN')), repel =TRUE, lable="ind", fill.ind = "white",  legend.title = "SP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13), legend.title=element_text(size=15, face="bold"))+
  theme_void()

#### Figure 7 (both joined), Save as 5x5 Landscape
plot.with.inset.BSTA <-ggdraw() + draw_plot(Main.BSTA) +
  draw_plot(Sub.BSTA, x = 0.191, y = 0.22, width = .60, height = .60)

## Table S4.(Tbio) - PERMANOVA
#### PERMANOVA Including both blocks and Sites in permutations
perm_Tbio <- how(nperm = 9999, 
                 blocks = BSTA_SB$Site,  
                 plots = Plots(BSTA_SB$Block), 
                 Within(type = "free"))

set.seed(1989)
ado.BSTA_random2 <- adonis2(PCA.crd.dec[, -1] ~ PCAD.BSTA$Year * PCAD.BSTA$Canopy * PCAD.BSTA$Treatment,
                            permutations = perm_Tbio,
                            method = "euc",
                            by = "term"
)
ado.BSTA_random2

#### Community approach (Tss) - Analysis with both Canopies
#### Select treatments of Tss (Exclude treatments "C", "Li", "Nu", "1F", "2F")
Spp.Tss <- Spp_All %>%
  filter(Treatment %in% c("C", "Ti", "To")) %>% 
  droplevels()

#### Separate years
Year.Tss <- split(Spp.Tss, Spp.Tss$Year)
Y2013.Tss <- Year.Tss$`2013`
Y2018.Tss <- Year.Tss$`2018`

### Combining the data of 2013 and 2018 for each stand
BSTA.Tss <- rbind(Y2013.Tss, Y2018.Tss)
BSTA.PCA.Tss <- BSTA.Tss[,6:71] #From column Y.C.T to the end of spp columns

####Make a matrix with Site and Block columns to add random factors in the following PERMANOVAs
BSTA.Tss_SB <- rbind(Y2013.Tss, Y2018.Tss)
BSTA.Tss_SB <- BSTA.Tss_SB[,2:3]

#### Centroid for Cntrl
Data.BSTA.Tss<-droplevels(subset(BSTA.PCA.Tss))
Data0.BSTA.Tss<- Data.BSTA.Tss[, colSums(Data.BSTA.Tss != 0) > 0] # To eliminate empty (0) columns

#### Data transformation and PCA
PCA.crd.dec.Tss <- decostand(Data0.BSTA.Tss[,-1], "hellinger")
PCA.BSTA.Tss <- prcomp(PCA.crd.dec.Tss, center = TRUE, scale. = FALSE)
PCA.crd.dec.Tss <- cbind(Data0.BSTA.Tss$Y.C.T, PCA.crd.dec.Tss) #Adding colmn with names
colnames(PCA.crd.dec.Tss)
names(PCA.crd.dec.Tss)[names(PCA.crd.dec.Tss) == "Data0.BSTA.Tss$Y.C.T"] <- "Y.C.T"
str(PCA.crd.dec.Tss)
PCAD.BSTA.Tss <- as.data.frame(PCA.BSTA.Tss$x)

PCAD.BSTA.Tss$Year <- with(PCA.crd.dec.Tss, substr(Y.C.T, 1,4))
PCAD.BSTA.Tss$Canopy <- with(PCA.crd.dec.Tss, substr(Y.C.T, 6,7))
PCAD.BSTA.Tss$Treatment <- with(PCA.crd.dec.Tss, substr(Y.C.T, 9,nchar(as.character(Y.C.T))))

#### Mean for each treatment and year
pca.centroids.BSTA.Tss <- aggregate(PCAD.BSTA.Tss[,1:2], list(Year = PCAD.BSTA.Tss$Year, Canopy=PCAD.BSTA.Tss$Canopy, Treatment=PCAD.BSTA.Tss$Treatment), mean)
pca.centroids.BSTA.Tss$LAB<-with(pca.centroids.BSTA.Tss, ifelse(Year==2018, "18", "13"))

#### Figure 7 (part), Save as 5x5 Landscape
Main.BSTA.Tss<-ggplot(data=pca.centroids.BSTA.Tss, aes(PC1, PC2, col=Treatment))+
  theme_bw()+
  geom_hline(yintercept = 0, lty=2, col="gray")+
  geom_vline(xintercept = 0, lty=2, col="gray")+
  geom_point(aes(shape = Year, color=Treatment),size = 2)+
  scale_color_manual(values=c('gray53', 'green', 'blue'))+
  geom_line(aes(PC1, PC2,  group=interaction(Canopy, Treatment), linetype=Canopy)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="bottom",
        plot.title = element_text(size=12,hjust = 0.5,face="bold"),
        axis.title=element_text(size=11),text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5, size=8),
        axis.text.y = element_text(angle=0,hjust=0.5,size=8))+ 
  labs(x= "PC1 (38%)", y ="PC2 (9%)", title="Species abundance under BS and TA stands")

summary(PCA.BSTA.Tss)$importance[2,1] # PCA1: 0.38404
summary(PCA.BSTA.Tss)$importance[2,2] # PCA2: 0.09073

#### Subplot of species
scores(PCA.BSTA.Tss$center) #To get the scores of species and select the most important ones

#### Figure 7 (part-Species)
#Selected species with the highest scores (0.36-0.04): PLS,CON,GAH,PTC,RUP,LIB,MAC,VIS,LYA,CLB,VAM,DIP,PES,POA,RIG,POC,GAA
Sub.BSTA.Tss<-fviz_pca_biplot(PCA.BSTA.Tss, label = "var",labelsize = 3, col.var = "black",habillage="none",title="", xlab="", ylab="", addEllipses=FALSE, geom.ind = "", select.var=list(name=c('PLS','CON','GAH','PTC','RUP','LIB','MAC','VIS','LYA','CLB','VAM','DIP','PES','POA','RIG','POC','GAA')), repel =TRUE, lable="ind", fill.ind = "white",  legend.title = "SP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13), legend.title=element_text(size=15, face="bold"))+
  theme_void()

#### Figure 7 (both joined), Save as 5x5 Landscape
plot.with.inset.BSTA.Tss <-ggdraw() + draw_plot(Main.BSTA.Tss) +
  draw_plot(Sub.BSTA.Tss, x = 0.212, y = 0.326, width = .60, height = .60)

## Table S4.(Tss) - PERMANOVA
#### PERMANOVA Including both blocks and Sites in permutations! (08/08/2022). If Year is included in the equation, I don't need to include  it as random factor (As discussed with Steve and Valentina)
perm <- how(nperm = 9999, 
            blocks = BSTA.Tss_SB$Site,  
            plots = Plots(BSTA.Tss_SB$Block), 
            Within(type = "free"))

set.seed(1989)
ado.BSTA.Tss_random2 <- adonis2(PCA.crd.dec.Tss[, -1] ~ PCAD.BSTA.Tss$Year * PCAD.BSTA.Tss$Canopy * PCAD.BSTA.Tss$Treatment,
                                permutations = perm,
                                method = "euc",
                                by = "term"
)
ado.BSTA.Tss_random2


## Ordinations for each objective and each forest type ####

## Figure 8. 
### Principal Coordinate Analysis (PCoA) based on the Bray-Curtis distance with Cailliez correction of the understory species abundance in 2013 (crossed squares) and 2018 (filled squares), with their corresponding percentage of variances in each axis. Treatments of the ecosystem approach in black spruce (a) and trembling aspen forests (b) correspond to Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange), Double-litter (2F, in red) and Control (C, in grey), with the corresponding centroids (standard deviation). Treatments of the community approach in black spruce (c) and trembling aspen forests (d) correspond to Control (C, in grey), Transplants-in (Ti, in green) and Transplants-out (To, in blue). Treatments in black spruce forests in light colors and in trembling aspen forests in dark colors. Arrows correspond to the most important understory species defining the community composition (Envfit, P < 0.05). Notice that scales are different.

## Table 3. 
### Differences in composition of understory plant communities between treatments and forest types.

## Table 4.
###Homogeneity of multivariance dispersion of understory plant abundance among treatments and years for each forest type.

## Tbio (Ecosystem approach)
#### PCoA separated by Canopy with 2013 and 2018 data, Treatments Tbio (Ecosystem approach)####
View(BSTA.Year)
BSTA.Year0<-droplevels(subset(BSTA.Year))
BSTA_All<- BSTA.Year0[, colSums(BSTA.Year0 != 0) > 0] # To eliminate empty (0) columns (1 was eliminated)
BSTA_All$Y.C.T <- NULL
View(BSTA_All) #Table with Year 2013 and 2018, with all other info, to separate between canopies

###Separate by Canopy
BSTA_Canopy<- split(BSTA_All, BSTA_All$Canopy)
BS_BegEnd <- BSTA_Canopy$'BS'
TA_BegEnd <- BSTA_Canopy$'TA'

#### Separate by Table of Vars and Spp
####BS
BS_BegEnd_Var <- select(BS_BegEnd, 1:5)
BS_BegEnd_Spp <- select(BS_BegEnd, 6:68)

####TA
TA_BegEnd_Var <- select(TA_BegEnd, 1:5)
TA_BegEnd_Spp <- select(TA_BegEnd, 6:68)

####Preparing data for PCoA
Farben <- c("gold","mediumorchid","salmon1","firebrick2","antiquewhite4","lightseagreen","chartreuse3")

## Tbio-TA
BS.bray <- vegdist(BS_BegEnd_Spp, method = "bray") # Distance matrix with method = bray (it's good for spp abundance)
is.euclid(BS.bray) #FALSE (Bad because then it has negative eigenvalues)

#### PCoA with cmdscale (by adding 'add = TRUE', I get the caillez correction)
BS.bray.pcoa <- cmdscale(BS.bray, k=2, eig=TRUE, add = TRUE) #Multidimensional scaling (also known as PCoA)

## Figure 8a 
#### PCoA plot (with bray-curtis distance) - Save PDF as portrait 6x5 (name InR_PCoA_BS_Tbio)
ordi.pcoa_BS <- ordiplot(BS.bray.pcoa, display="sites", type="points", cex = 0.2)

points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2013' & BS_BegEnd_Var$Treatment == '1F',], pch=7, col="salmon1", bg="salmon1", cex = 1.3)
points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2018' & BS_BegEnd_Var$Treatment == '1F',], pch=15, col="salmon1", bg="salmon1", cex = 1.3)

points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2013' & BS_BegEnd_Var$Treatment == '2F',], pch=7, col="firebrick2", bg="firebrick2", cex = 1.3)
points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2018' & BS_BegEnd_Var$Treatment == '2F',], pch=15, col="firebrick2", bg="firebrick2", cex = 1.3)

points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2013' & BS_BegEnd_Var$Treatment == 'C',], pch=7, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)
points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2018' & BS_BegEnd_Var$Treatment == 'C',], pch=15, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)

points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2013' & BS_BegEnd_Var$Treatment == 'Li',], pch=7, col="gold", bg="gold", cex = 1.3)
points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2018' & BS_BegEnd_Var$Treatment == 'Li',], pch=15, col="gold", bg="gold", cex = 1.3)

points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2013' & BS_BegEnd_Var$Treatment == 'Nu',], pch=7, col="mediumorchid", bg="mediumorchid", cex = 1.3)
points(ordi.pcoa_BS$sites[BS_BegEnd_Var$Year == '2018' & BS_BegEnd_Var$Treatment == 'Nu',], pch=15, col="mediumorchid", bg="mediumorchid", cex = 1.3)

ordiellipse(ordi.pcoa_BS, BS_BegEnd_Var$Treatment, label=TRUE, col = Farben, cex=1.2, font=2, lty=2, lwd=1)

plot(envfit(ordi.pcoa_BS, BS_BegEnd_Spp), p.max=0.05, col="gray25", cex=1, font=2)

#Legend 1
legend(0.58, 0.4, 
       legend = c("2013","2018"), #names to display
       col = c("gray25", "gray25"), 
       pch = c(7,15), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = "gray25", 
       horiz = F , 
       inset = c(1.5, 1.5)) #% (from 0 to 1) to draw the legend away from x and y axis. You can also give the X and Y coordinate of the legend: legend(3, 5, ...)
#Legend 2
legend(0.58, 0.4, 
       legend = c("Li","Nu","1F","2F","C","Ti","To"), #names to display
       col = Farben, 
       pch = c(0,0,0,0,0,0,0), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = Farben, 
       horiz = F)

#### Eigenvalues and their proportion to the sum
BS.bray.pcoa$eig[1:3] # 14.650921 10.682960  7.390236
#### Calculate the percent of variance explained by first two axes
100*BS.bray.pcoa$eig[1:3]/sum(BS.bray.pcoa$eig) # 15.020097 10.952151  7.576457

## Table 3 - PERMANOVA Tbio-BS
#### Including both blocks and Sites in permutations
perm_BS_Tbio <- how(nperm = 9999, 
                    blocks = BS_BegEnd_Var$Site,  
                    plots = Plots(BS_BegEnd_Var$Block), 
                    Within(type = "free"))

set.seed(1989)
BS.bray.pcoa_perm_bray2 <- adonis2(BS.bray ~ Treatment*Year, data=BS_BegEnd_Var, 
                                   permutations = perm_BS_Tbio,
                                   method = "bray")
BS.bray.pcoa_perm_bray2

### Table 4. Tbio-BS
#### Posthoc tests - PERMUTEST
BS_BegEnd_Var2 <- BS_BegEnd_Var
BS_BegEnd_Var2$YT <- paste(BS_BegEnd_Var2$Year,BS_BegEnd_Var2$Treatment)

dispersion_BS_Tbio <- betadisper(BS.bray, group=BS_BegEnd_Var2$YT,type = c("median","centroid"), bias.adjust = TRUE)
dispersion_BS_Tbio # !!! To get the contrast (average distance to median)
perm_BSTbio <- permutest(dispersion_BS_Tbio, pairwise = TRUE, permutations = 9999) # !!! TO get the P value in contrasts

anova(dispersion_BS_Tbio)

mod.HSD_BS_Tbio <- TukeyHSD(dispersion_BS_Tbio)
mod.HSD_BS_Tbio
plot(mod.HSD_BS_Tbio)
pairwise.adonis2(BS.bray ~ YT, data = BS_BegEnd_Var2)

## Tbio-TA
TA.bray <- vegdist(TA_BegEnd_Spp, method = "bray") # Distance matrix with method = bray (it's good for spp abundance)
is.euclid(TA.bray) #FALSE

## PCoA with cmdscale (by adding 'add = TRUE', I get the caillez correction)
TA.bray.pcoa <- cmdscale(TA.bray, k=2, eig=TRUE, add = TRUE) #Multidimensional scaling (also known as PCoA)

## Figure 8b
#### PCoA plot (with bray-curtis distance) - Save PDF as portrait 6x5 (name InR_PCoA_TA_Tbio)
ordi.pcoa_TA <- ordiplot(TA.bray.pcoa, display="sites", type="points", cex = 0.2)

points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2013' & TA_BegEnd_Var$Treatment == '1F',], pch=7, col="salmon1", bg="salmon1", cex = 1.3)
points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2018' & TA_BegEnd_Var$Treatment == '1F',], pch=15, col="salmon1", bg="salmon1", cex = 1.3)

points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2013' & TA_BegEnd_Var$Treatment == '2F',], pch=7, col="firebrick2", bg="firebrick2", cex = 1.3)
points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2018' & TA_BegEnd_Var$Treatment == '2F',], pch=15, col="firebrick2", bg="firebrick2", cex = 1.3)

points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2013' & TA_BegEnd_Var$Treatment == 'C',], pch=7, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)
points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2018' & TA_BegEnd_Var$Treatment == 'C',], pch=15, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)

points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2013' & TA_BegEnd_Var$Treatment == 'Li',], pch=7, col="gold", bg="gold", cex = 1.3)
points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2018' & TA_BegEnd_Var$Treatment == 'Li',], pch=15, col="gold", bg="gold", cex = 1.3)

points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2013' & TA_BegEnd_Var$Treatment == 'Nu',], pch=7, col="mediumorchid", bg="mediumorchid", cex = 1.3)
points(ordi.pcoa_TA$sites[TA_BegEnd_Var$Year == '2018' & TA_BegEnd_Var$Treatment == 'Nu',], pch=15, col="mediumorchid", bg="mediumorchid", cex = 1.3)

ordiellipse(ordi.pcoa_TA, TA_BegEnd_Var$Year, label=TRUE, col = "gray25", cex=1.2, font=2, lty=1, lwd=1) #the default for the ellipse is the standard deviation (kind="sd")
ordiellipse(ordi.pcoa_TA, TA_BegEnd_Var$Treatment, label=TRUE, col = Farben, cex=1.2, font=2, lty=2, lwd=1)

plot(envfit(ordi.pcoa_TA, TA_BegEnd_Spp), p.max=0.05, col="gray25", cex=1, font=2)

#Legend 1
legend(0.58, 0.4, 
       legend = c("2013","2018"), #names to display
       col = c("gray25", "gray25"), 
       pch = c(7,15), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = "gray25", 
       horiz = F , 
       inset = c(1.5, 1.5)) #% (from 0 to 1) to draw the legend away from x and y axis. You can also give the X and Y coordinate of the legend: legend(3, 5, ...)
#Legend 2
legend(1.5, 0.4, 
       legend = c("Li","Nu","1F","2F","C","Ti","To"), #names to display
       col = Farben, 
       pch = c(0,0,0,0,0,0,0), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = Farben, 
       horiz = F)

#### Eigenvalues and their proportion to the sum
TA.bray.pcoa$eig[1:3] # 10.419424  7.126796  4.062613
#### Calculate the percent of variance explained by first two axes
100*TA.bray.pcoa$eig[1:3]/sum(TA.bray.pcoa$eig) # 14.647095 10.018487  5.711014

## Table 3 - PERMANOVA Tbio-TA
#### Including both blocks and Sites in permutations
perm_TA_Tbio <- how(nperm = 9999, 
                    blocks = TA_BegEnd_Var$Site,  
                    plots = Plots(TA_BegEnd_Var$Block), 
                    Within(type = "free"))

set.seed(1989)
TA.bray.pcoa_perm_bray2 <- adonis2(TA.bray ~ Treatment*Year, data=TA_BegEnd_Var, 
                                   permutations = perm_TA_Tbio,
                                   method = "bray")
TA.bray.pcoa_perm_bray2

### Table 4. Tbio-TA
#### Posthoc tests - PERMUTEST
TA_BegEnd_Var2 <- TA_BegEnd_Var
TA_BegEnd_Var2$YT <- paste(TA_BegEnd_Var2$Year,TA_BegEnd_Var2$Treatment)

dispersion_TA_Tbio <- betadisper(TA.bray, group=TA_BegEnd_Var2$YT,type = c("median","centroid"), bias.adjust = TRUE)
dispersion_TA_Tbio # !!! To get the contrast (average distance to median)
perm_TATbio <- permutest(dispersion_TA_Tbio, pairwise = TRUE, permutations = 9999) # !!! TO get the P value in contrasts

anova(dispersion_TA_Tbio)

mod.HSD_TA_Tbio <- TukeyHSD(dispersion_TA_Tbio)
mod.HSD_TA_Tbio
plot(mod.HSD_TA_Tbio)
pairwise.adonis2(TA.bray ~ YT, data = TA_BegEnd_Var2)

## Tss (Community approach)
#### PCoA separated by Canopy with 2013 and 2018 data, Treatments Tss (Community approach)####
View(BSTA.Tss)
BSTA.Tss0<-droplevels(subset(BSTA.Tss))
BSTA_All_Tss<- BSTA.Tss0[, colSums(BSTA.Tss0 != 0) > 0] # To eliminate empty (0) columns (12 were eliminated)
BSTA_All_Tss$Y.C.T <- NULL
View(BSTA_All_Tss) #Table with Year 2013 and 2018, with all other info, to separate between canopies

####Separate by Canopy
BSTA_Tss_Canopy<- split(BSTA_All_Tss, BSTA_All_Tss$Canopy)
BS_BegEnd_Tss <- BSTA_Tss_Canopy$'BS'
TA_BegEnd_Tss <- BSTA_Tss_Canopy$'TA'

#### Separate by Table of Vars and Spp
####BS
BS_BegEnd_Tss_Var <- select(BS_BegEnd_Tss, 1:5)
BS_BegEnd_Tss_Spp <- select(BS_BegEnd_Tss, 6:59)

####TA
TA_BegEnd_Tss_Var <- select(TA_BegEnd_Tss, 1:5)
TA_BegEnd_Tss_Spp <- select(TA_BegEnd_Tss, 6:59)

####Preparing data for PCoA
Farben_Tss <- c("lightseagreen","chartreuse3","antiquewhite4")

## Tss-BS
BS.bray_Tss <- vegdist(BS_BegEnd_Tss_Spp, method = "bray") # Distance matrix with method = bray (it's good for spp abundance)
is.euclid(BS.bray_Tss) #FALSE

#### PCoA with cmdscale (by adding 'add = TRUE', I get the caillez correction)
BS.bray.pcoa_Tss <- cmdscale(BS.bray_Tss, k=2, eig=TRUE, add = TRUE) #Multidimansional scaling (also known as PCoA)

## Figure 8c
#### PCoA plot (with bray-curtis distance) - Save PDF as portrait 6x5 (name InR_PCoA_BS_Tss)
ordi.pcoa_BS_Tss <- ordiplot(BS.bray.pcoa_Tss, display="sites", type="points", cex = 0.2)

points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2013' & BS_BegEnd_Tss_Var$Treatment == 'C',], pch=7, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)
points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2018' & BS_BegEnd_Tss_Var$Treatment == 'C',], pch=15, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)

points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2013' & BS_BegEnd_Tss_Var$Treatment == 'Ti',], pch=7, col="chartreuse3", bg="chartreuse3", cex = 1.3)
points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2018' & BS_BegEnd_Tss_Var$Treatment == 'Ti',], pch=15, col="chartreuse3", bg="chartreuse3", cex = 1.3)

points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2013' & BS_BegEnd_Tss_Var$Treatment == 'To',], pch=7, col="lightseagreen", bg="lightseagreen", cex = 1.3)
points(ordi.pcoa_BS_Tss$sites[BS_BegEnd_Tss_Var$Year == '2018' & BS_BegEnd_Tss_Var$Treatment == 'To',], pch=15, col="lightseagreen", bg="lightseagreen", cex = 1.3)

ordiellipse(ordi.pcoa_BS_Tss, BS_BegEnd_Tss_Var$Treatment, label=TRUE, col = Farben_Tss, cex=1.2, font=2, lty=2, lwd=1)

plot(envfit(ordi.pcoa_BS_Tss, BS_BegEnd_Tss_Spp), p.max=0.05, col="gray25", cex=1, font=2)

#Legend 1
legend(0.58, 0.4, 
       legend = c("2013","2018"), #names to display
       col = c("gray25", "gray25"), 
       pch = c(7,15), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = "gray25", 
       horiz = F , 
       inset = c(1.5, 1.5)) #% (from 0 to 1) to draw the legend away from x and y axis. You can also give the X and Y coordinate of the legend: legend(3, 5, ...)
#Legend 2
legend(0.58, 0.4, 
       legend = c("C","Ti","To"), #names to display
       col = Farben_Tss, 
       pch = c(0,0,0), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = Farben_Tss, 
       horiz = F)

#### Eigenvalues and their proportion to the sum 
BS.bray.pcoa_Tss$eig[1:3] # 12.605008  4.847055  2.556525
#### Calculate the percent of variance explained by first two axes
100*BS.bray.pcoa_Tss$eig[1:3]/sum(BS.bray.pcoa_Tss$eig) # 27.376608 10.527239  5.552475

## Table 3 - PERMANOVA Tss-BS
#### Including both blocks and Sites in permutations
perm_BS_Tss <- how(nperm = 9999, 
                   blocks = BS_BegEnd_Tss_Var$Site,  
                   plots = Plots(BS_BegEnd_Tss_Var$Block), 
                   Within(type = "free"))

set.seed(1989)
BS.bray_Tss_perm_bray2 <- adonis2(BS.bray_Tss ~ Treatment*Year, data=BS_BegEnd_Tss_Var, 
                                  permutations = perm_BS_Tss,
                                  method = "bray")
BS.bray_Tss_perm_bray2

### Table 4. Tbio-TA
#### Posthoc tests - PERMUTEST
BS_BegEnd_Tss_Var2 <- BS_BegEnd_Tss_Var
BS_BegEnd_Tss_Var2$YT <- paste(BS_BegEnd_Tss_Var2$Year,BS_BegEnd_Tss_Var2$Treatment)

dispersion_BS_Tss <- betadisper(BS.bray_Tss, group=BS_BegEnd_Tss_Var2$YT,type = c("median","centroid"), bias.adjust = TRUE)
dispersion_BS_Tss # !!! To get the contrast (average distance to median)
perm_BSTss <- permutest(dispersion_BS_Tss, pairwise = TRUE, permutations = 9999) # !!! TO get the P value in contrasts

anova(dispersion_BS_Tss)

mod.HSD_BS_Tss <- TukeyHSD(dispersion_BS_Tss)
mod.HSD_BS_Tss
plot(mod.HSD_BS_Tss)
pairwise.adonis2(BS.bray_Tss ~ YT, data = BS_BegEnd_Var2)

## Tss-TA
TA.bray_Tss <- vegdist(TA_BegEnd_Tss_Spp, method = "bray") # Distance matrix with method = bray (it's good for spp abundance)
is.euclid(TA.bray_Tss) #FALSE

#### PCoA with cmdscale (by adding 'add = TRUE', I get the caillez correction)
TA.bray.pcoa_Tss <- cmdscale(TA.bray_Tss, k=2, eig=TRUE, add = TRUE) #Multidimensional scaling (also known as PCoA)

## Figure 8d
#### PCoA plot (with bray-curtis distance) - Save PDF as portrait 6x5 (name InR_PCoA_TA_Tss)
ordi.pcoa_TA_Tss <- ordiplot(TA.bray.pcoa_Tss, display="sites", type="points", cex = 0.2, col="cornsilk2")

points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2013' & TA_BegEnd_Tss_Var$Treatment == 'C',], pch=7, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)
points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2018' & TA_BegEnd_Tss_Var$Treatment == 'C',], pch=15, col="antiquewhite4", bg="antiquewhite4", cex = 1.3)

points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2013' & TA_BegEnd_Tss_Var$Treatment == 'Ti',], pch=7, col="chartreuse3", bg="chartreuse3", cex = 1.3)
points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2018' & TA_BegEnd_Tss_Var$Treatment == 'Ti',], pch=15, col="chartreuse3", bg="chartreuse3", cex = 1.3)

points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2013' & TA_BegEnd_Tss_Var$Treatment == 'To',], pch=7, col="lightseagreen", bg="lightseagreen", cex = 1.3)
points(ordi.pcoa_TA_Tss$sites[TA_BegEnd_Tss_Var$Year == '2018' & TA_BegEnd_Tss_Var$Treatment == 'To',], pch=15, col="lightseagreen", bg="lightseagreen", cex = 1.3)

ordiellipse(ordi.pcoa_TA_Tss, TA_BegEnd_Tss_Var$Year, label=TRUE, col = "gray25", cex=1.2, font=2, lty=1, lwd=1) #the default for the ellipse is the standard deviation (kind="sd")
ordiellipse(ordi.pcoa_TA_Tss, TA_BegEnd_Tss_Var$Treatment, label=TRUE, col = Farben_Tss, cex=1.2, font=2, lty=2, lwd=1)

plot(envfit(ordi.pcoa_TA_Tss, TA_BegEnd_Tss_Spp), p.max=0.05, col="gray25", cex=1, font=2)

#Legend 1
legend(0.58, 0.4, 
       legend = c("2013","2018"), #names to display
       col = c("gray25", "gray25"), 
       pch = c(7,15), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = "gray25", 
       horiz = F , 
       inset = c(1.5, 1.5)) #% (from 0 to 1) to draw the legend away from x and y axis. You can also give the X and Y coordinate of the legend: legend(3, 5, ...)
#Legend 2
legend(-0.5, 0.4, 
       legend = c("C","Ti","To"), #names to display
       col = Farben_Tss, 
       pch = c(0,0,0), #symbol type
       bty = "n", #type of box around the legend
       pt.cex = 1.3, #symbol size
       cex = 1, #text size
       text.col = Farben_Tss, 
       horiz = F)

#### Eigenvalues and their proportion to the sum 
TA.bray.pcoa_Tss$eig[1:3] # 10.131209  4.704789  3.861210
#### Calculate the percent of variance explained by first two axes
100*TA.bray.pcoa_Tss$eig[1:3]/sum(TA.bray.pcoa_Tss$eig) # 22.128320 10.276077  8.433552

## Table 3 - PERMANOVA Tss-TA
#### Including both blocks and Sites in permutations
perm_TA_Tss <- how(nperm = 9999, 
                   blocks = TA_BegEnd_Tss_Var$Site,  
                   plots = Plots(TA_BegEnd_Tss_Var$Block), 
                   Within(type = "free"))

set.seed(1989)
TA.bray_Tss_perm_bray2 <- adonis2(TA.bray_Tss ~ Treatment*Year, data=TA_BegEnd_Tss_Var, 
                                  permutations = perm_TA_Tss,
                                  method = "bray")
TA.bray_Tss_perm_bray2

### Table 4. Tbio-TA
#### Posthoc tests - PERMUTEST
TA_BegEnd_Tss_Var2 <- TA_BegEnd_Tss_Var
TA_BegEnd_Tss_Var2$YT <- paste(TA_BegEnd_Tss_Var2$Year,TA_BegEnd_Tss_Var2$Treatment)

dispersion_TA_Tss <- betadisper(TA.bray_Tss, group=TA_BegEnd_Tss_Var2$YT,type = c("median","centroid"), bias.adjust = TRUE)
dispersion_TA_Tss # !!! To get the contrast (average distance to median)
perm_TATss <- permutest(dispersion_TA_Tss, pairwise = TRUE, permutations = 9999) # !!! TO get the P value in contrasts

anova(dispersion_TA_Tss)

mod.HSD_TA_Tss <- TukeyHSD(dispersion_TA_Tss)
mod.HSD_TA_Tss
plot(mod.HSD_TA_Tss)
pairwise.adonis2(TA.bray_Tss ~ YT, data = TA_BegEnd_Var2)

## Abundance of species/Functional groups in 2018 ####
#### Supp.Table 1 and Boxplots plot (Figure 2) for functional groups in 2018
DB_All3 <- DB_All
DB_All3_Year <- split(DB_All3, DB_All3$Year)
DB_All3_2018 <- DB_All3_Year$`2018`

## Table S1.
### Average abundance of understory plants of the cumulative effect of treatments in 2018 in black spruce and trembling aspen forests, for the Control (C, in grey) and the different treatments of the ecosystem approach (Objective 1): Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red), and of the community approach (Objective 2): Transplants-out (To, in blue), Transplants-in (Ti, in green). Species, with their corresponding abbreviations, were classified in functional groups: Bryophytes, Sphagnaceae, Ericaceae, Pteridophyta, Graminoids, Herbs, and Shrubs, and Trees. Data correspond to the average and standard deviation (in grey numbers) of blocks and sites for each species among treatments and forest types. Values for each functional group (in colors) correspond to the average abundance per functional group for each treatment.
####Do a sum of Coverage per Func_gr
Func18 <- DB_All3_2018
Func18$Year <- NULL

####make a joint column
Func18_united <- unite(Func18, Site, Canopy, Treatment, Species, Abbr, Func_gr, col = "SCTSAF", sep = "_")

####Do a mean of Coverage per Block
Func18_mean <- aggregate(Func18_united$Coverage, list(Func18_united$SCTSAF), FUN=mean)

####Arrange columns again
Func18_mean2 <- separate(Func18_mean, Group.1, into = c("Site","Canopy","Treatment","Species","Abbr","Func_gr"), sep = "_")
Func18_mean3 <- unite(Func18_mean2, Canopy, Treatment, Species, Abbr, Func_gr, col = "CTSAF", sep = "_")

####Do a mean of Coverage per Site (from mean of blocks)
Func18_mean4 <- aggregate(Func18_mean3$x, list(Func18_mean3$CTSAF), FUN=mean)
Func18_SD <- aggregate(Func18_mean3$x, list(Func18_mean3$CTSAF), FUN=sd)

####Arrange columns again
Func18_mean5 <- separate(Func18_mean4, Group.1, into = c("Canopy","Treatment","Species","Abbr","Func_gr"), sep = "_")
Func18_mean6 <- unite(Func18_mean5, Canopy, Species, Abbr, Func_gr, col = "CSAF", sep = "_")

Func18_SD2 <- separate(Func18_SD, Group.1, into = c("Canopy","Treatment","Species","Abbr","Func_gr"), sep = "_")
Func18_SD3 <- unite(Func18_SD2, Canopy, Species, Abbr, Func_gr, col = "CSAF", sep = "_")

####Spread Treatmetns
Func18_mean7 <- spread(Func18_mean6, Treatment, x, fill = Func18_mean6$x)
Func18_SD4 <- spread(Func18_SD3, Treatment, x, fill = Func18_SD3$x)

####Merge to get a tbale with mean and SD for each treatment
Func18_merged <- merge(Func18_mean7,Func18_SD4, by="CSAF")
Func18_merged2 <- separate(Func18_merged, CSAF, into = c("Canopy","Species","Abbr","Func_gr"), sep = "_")

####Merge again from tables before spreading to get a table with mean and SD for each treatment (but witouth spread to arrange in Excel)
Func18_merged3 <- merge(Func18_mean4,Func18_SD, by="Group.1")
Func18_merged4 <- separate(Func18_merged3, Group.1, into = c("Canopy","Treatment","Species","Abbr","Func_gr"), sep = "_")
####write.csv(Func18_merged4, "Func18_merged4.csv")

#### Prepare data for Boxplots (Figure 2) for functional groups in 2018 only C treatment
DB3_2018_Treat <- split(DB_All3_2018, DB_All3_2018$Treatment)
DB3_2018_C <- DB3_2018_Treat$'C'

#####Calculate cover by C plot for each of each Func-gr
DB3_2018_C2 <- DB3_2018_C
DB3_2018_C2$Species <- NULL
DB3_2018_C2$Abbr <- NULL
DB3_2018_C2$SBC <- paste0(DB3_2018_C2$Site,"_", DB3_2018_C2$Block, "_", DB3_2018_C2$Canopy) #Create a single column for variables
DB3_2018_C2$Year <- NULL
DB3_2018_C2$Treatment <- NULL

####Do a sum of Coverage per Func_gr
Func18_C <- DB3_2018_C2 %>%
  group_by(Site, Block, Canopy, Func_gr) %>%
  summarise(Coverage = sum(Coverage))

#### make a column joining Site and Block
Func18_C2 <- unite(Func18_C, Site, Canopy, Func_gr, col = "Site_Can_Func", sep = " ")
####Do a mean of Coverage per Block
Func18_C_mean <- aggregate(Func18_C2$Coverage, list(Func18_C2$Site_Can_Func), FUN=mean)
####Arrange columns again
Func18_C_mean2 <- separate(Func18_C_mean, Group.1, into = c("Site", "Canopy", "Func_gr"), sep = " ")
Func18_C_mean3 <- unite(Func18_C_mean2, Canopy, Func_gr, col = "Can_Func", sep = " ")
####Do a mean and SD of Coverage (x) per Site
Func18_C_mean4 <- aggregate(Func18_C_mean3$x, list(Func18_C_mean3$Can_Func), FUN=mean)
Func18_C_SD <- aggregate(Func18_C_mean3$x, list(Func18_C_mean3$Can_Func), FUN=sd)
Func18_C_mean5 <- merge(Func18_C_mean4,Func18_C_SD, by="Group.1")
####Arrange columns again
colnames(Func18_C_mean5)[2] <- "Mean"
colnames(Func18_C_mean5)[3] <- "SD"
Func18_C_mean6 <- separate(Func18_C_mean5, Group.1, into = c("Canopy", "Func_gr"), sep = " ")
####With Func18_C_mean6 I got the table with the values of Mean and SD per Func_gr per Canopy, but for the figure, I must use data before calculating the mean, so (Func18_C_mean3)

####Arrange Func18_C_mean3 (with Mean per Blocks) for the boxplot
Func18_C_mean7 <- separate(Func18_C_mean3, Can_Func, into = c("Canopy", "Func_gr"), sep = " ")
colnames(Func18_C_mean7)[4] <- "Mean"
Func18_C_mean8 <- Func18_C_mean7
Func18_C_mean8$Site <- NULL

## Figure 2. 
### Experimental design for overstory – understory relationships. a) The study area is located in the boreal forest of eastern Canada; b) above the parallel 49⁰N in western Quebec. c) The experiment was conducted in three sites A, B and C (Sites B and C are 0.5 km apart and Site A is 5.3 km away). d,e,f) Each site has adjacent stands (separated from 34 to 115 m) dominated by black spruce (BS, as triangles) and trembling aspen (TA, as circles), in which we placed 3 blocks (B1, B2, B3) separated from 16 to 49 m. g). In each block, one plot per treatment was placed to obtain a nested block experimental design with a total of 126 plots (3 sites x 2 forest types x 3 blocks x 7 treatments = 126 plots).. Vegetation abundance in all treatments was measured every year from 2013 to 2018 by using a 50 x 50 cm grid placed in the middle of each plot. Treatments in different colors correspond to Control conditions (C, in grey) and to treatments of the ecosystem approach: Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red); and of the community approach: Transplants-out (To, in blue), Transplants-in (Ti, in green). Treatments in black spruce forests in light colors and in trembling aspen forests in dark colors.
#### Boxplot Coverage of Control plots in 2018 (average per Site (from averages per Blocks)) Save as PDF 8x6 portrait
cov.plot <- ggplot(Func18_C_mean8, aes(x = Func_gr, y =Mean, fill=Func_gr)) +
  geom_boxplot()+
  facet_grid(cols=vars(Canopy)) +
  scale_fill_brewer(palette="Dark2") +
  theme_light(base_size = 12)

cov.plot.lm <- lm(Mean ~ Func_gr + Canopy + Func_gr:Canopy, data = Func18_C_mean8)
summary(cov.plot.lm)

cov.emmeans <- emmeans(cov.plot.lm, specs = pairwise ~ Func_gr|Canopy, type = "response")
cov.emmeans

plot(cov.emmeans, comparisons = TRUE)

cov.mod_means <- multcomp::cld(object = cov.emmeans$emmeans,
                               
                               Letters = letters,no.readonly=TRUE)
cov.mod_means

## Figure S3 and Figure S4 - Smoothlines ####
#### Variation in abundance of functional groups (Bryophytes, Ericaceae, Pteridophyta, Herbs and Shrubs) over time (from 2013 to 2018) for black spruce (left panels) and trembling aspen stands (right panels). 

#### Load main Table with all information and with the five Func-gr
DB_All_FG <- DB_All

#### Reduce table to leave just Func_gr and sum coverage (eliminating Species and Abbr)
DB_All_FG2 <- DB_All_FG %>%
  group_by(Year, Site, Block, Canopy, Treatment, Func_gr) %>%
  summarise(Coverage = sum(Coverage))
head(DB_All_FG2)
####Rename levels of Canopy
DB_All_FG2$Canopy <- recode_factor(DB_All_FG2$Canopy, BS  = "Black spruce", TA = "Trembling aspen")

## Tbio ###
## Figure S2.
### Variation in abundance of functional groups over time (from 2013 to 2018) for black spruce (left panels) and trembling aspen stands (right panels). Functional groups: (a) Bryophytes, (b) Sphagnaceae, (c) Ericaceae, (d) Peridophyta, (e) Herbs, (f) Shrubs, and (g) Trees. Treatments correspond to ecosystem approach: Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red) and Control conditions (C, in grey). Smooth lines are based on the linear model of understory species abundances, each point in different colors (treatments) from data each year of the study.
#### Filter by treatments Tbio to leave only (C, Li, Nu, 1F and 2F)
DB_Tbio_FG <- DB_All_FG2 %>%
  filter(Treatment %in% c("C", "Li", "Nu", "1F", "2F")) %>% 
  droplevels()
unique(DB_Tbio_FG$Treatment)

### Tbio, Func_gr: Herbs
DB_Tbio_Her <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Herbs")) %>% 
  droplevels()
unique(DB_Tbio_Her$Func_gr)

####Herbs (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Her <- ggplot(DB_Tbio_Her, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Herbs")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Her

### Tbio, Func_gr: Bryophytes
DB_Tbio_Bry <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Bryophytes")) %>% 
  droplevels()
unique(DB_Tbio_Bry$Func_gr)

####Bryophytes (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Bry <- ggplot(DB_Tbio_Bry, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Bryophytes")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Bry

### Tbio, Func_gr: Ericaceae
DB_Tbio_Eri <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Ericaceae")) %>% 
  droplevels()
unique(DB_Tbio_Eri$Func_gr)

####Ericaceae (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Eri <- ggplot(DB_Tbio_Eri, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Ericaceae")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Eri

### Tbio, Func_gr: Pteridophyta
DB_Tbio_Pte <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Pteridophyta")) %>% 
  droplevels()
unique(DB_Tbio_Pte$Func_gr)

####Pteridophyta (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Pte <- ggplot(DB_Tbio_Pte, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Pteridophyta")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Pte

### Tbio, Func_gr: Shrubs
DB_Tbio_Shr <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Shrubs")) %>% 
  droplevels()
unique(DB_Tbio_Shr$Func_gr)

####Shrubs (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Sch <- ggplot(DB_Tbio_Shr, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Shrubs")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Sch

### Tbio, Func_gr: Sphagnaceae
DB_Tbio_Sph <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Sphagnaceae")) %>% 
  droplevels()
unique(DB_Tbio_Sph$Func_gr)

####Sphagnaceae (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Sph <- ggplot(DB_Tbio_Sph, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Sphagnaceae")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Sph

### Tbio, Func_gr: Trees
DB_Tbio_Tre <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Trees")) %>% 
  droplevels()
unique(DB_Tbio_Tre$Func_gr)

####Trees (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Tre <- ggplot(DB_Tbio_Tre, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Trees")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Tre

### Tbio, Func_gr: Graminoids
DB_Tbio_Gra <- DB_Tbio_FG %>% 
  filter(Func_gr %in% c("Graminoids")) %>% 
  droplevels()
unique(DB_Tbio_Gra$Func_gr)

####Graminoids (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tbio_Gra <- ggplot(DB_Tbio_Gra, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Graminoids")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("salmon1","firebrick2","antiquewhite4","gold","mediumorchid"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tbio_Gra

## Tss ###
## Figure S3.
### Variation in abundance of functional groups over time (from 2013 to 2018) for black spruce (left panels) and trembling aspen stands (right panels). Functional groups: (a) Bryophytes, (b) Sphagnaceae, (c) Ericaceae, (d) Peridophyta, (e) Herbs, (f) Shrubs, and (g) Trees. Treatments correspond to ecosystem approach: Light (Li, in yellow), Nutrients (Nu, in purple), Single-litter (1F, in orange) and Double-litter (2F, in red) and Control conditions (C, in grey). Smooth lines are based on the linear model of understory species abundances, each point in different colors (treatments) from data each year of the study.
#### Filter by treatments Tss to leave only (C, Ti, To)
DB_Tss_FG <- DB_All_FG2 %>%
  filter(Treatment %in% c("C", "Ti", "To")) %>% 
  droplevels()
unique(DB_Tss_FG$Treatment)

### Tss, Func_gr: Herbs
DB_Tss_Her <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Herbs")) %>% 
  droplevels()
unique(DB_Tss_Her$Func_gr)

####Herbs (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Her <- ggplot(DB_Tss_Her, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Herbs")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Her

### Tss, Func_gr: Bryophytes
DB_Tss_Bry <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Bryophytes")) %>% 
  droplevels()
unique(DB_Tss_Bry$Func_gr)

####Bryophytes (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Bry <- ggplot(DB_Tss_Bry, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Bryophytes")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Bry

### Tss, Func_gr: Ericaceae
DB_Tss_Eri <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Ericaceae")) %>% 
  droplevels()
unique(DB_Tss_Eri$Func_gr)

####Ericaceae (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Eri <- ggplot(DB_Tss_Eri, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Ericaceae")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Eri

### Tss, Func_gr: Pteridophyta
DB_Tss_Pte <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Pteridophyta")) %>% 
  droplevels()
unique(DB_Tss_Pte$Func_gr)

####Pteridophyta (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Pte <- ggplot(DB_Tss_Pte, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Pteridophyta")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Pte

### Tss, Func_gr: Shrubs
DB_Tss_Shr <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Shrubs")) %>% 
  droplevels()
unique(DB_Tss_Shr$Func_gr)

####Shrubs (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Sch <- ggplot(DB_Tss_Shr, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Shrubs")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Sch

### Tss, Func_gr: Sphagnaceae
DB_Tss_Sph <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Sphagnaceae")) %>% 
  droplevels()
unique(DB_Tss_Sph$Func_gr)

####Sphagnaceae (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Sph <- ggplot(DB_Tss_Sph, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Sphagnaceae")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Sph

### Tss, Func_gr: Trees
DB_Tss_Tre <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Trees")) %>% 
  droplevels()
unique(DB_Tss_Tre$Func_gr)

####Trees (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Tre <- ggplot(DB_Tss_Tre, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Trees")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Tre

### Tss, Func_gr: Graminoids
DB_Tss_Gra <- DB_Tss_FG %>% 
  filter(Func_gr %in% c("Graminoids")) %>% 
  droplevels()
unique(DB_Tss_Gra$Func_gr)

####Graminoids (ggplot with smoothed lines - Abundance per Func_gr over time) (PDF 4x10 landscape)
Plot_Tss_Gra <- ggplot(DB_Tss_Gra, aes(x=Year, y=Coverage, group=Treatment, color=Treatment, pch=Treatment)) +
  facet_grid(~Canopy) +
  geom_point (size = 1.5, alpha = 0.5)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=1.5)+
  theme_bw()+
  labs(x="Year", y = "Abundance of Graminoids")+
  scale_shape_manual(values=c(17, 16, 17, 16, 17))+
  scale_color_manual(values=c("antiquewhite4","lightseagreen","chartreuse3"))+
  theme(axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14,face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
Plot_Tss_Gra

#### -- END
