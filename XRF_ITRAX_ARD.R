
# STILL TO DO :-  

# 1 - add PCA code to ARD and YAN from Lago Pato analysis 
# 2 - write code to add age-depth model columns from BACON output
# 3 - recreate imported csv structure for matching ITRAX & XRF datasets done in excel into R / Tidyverse
# 4 - reproduce Figure 3 plots from original dataset

# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Load libraries & colours ----------------------------------------------------------

##load libraries
library(tidyverse)
library(tidypaleo)
library(readr)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(cowplot) # for plotting
library(vegan)
library(rioja)
library(ellipse)  # for PCA and cluster
library(factoextra) # for PCA and cluster
library(reshape2)
library(GGally)
library(ggsci)
library(ggdendro)
library(dendextend)
library(dynamicTreeCut)
library(colorspace)
library(cluster)
library(magrittr) # for piping %>%
library(mgcv)
library(gtable)
library(repr)
library(bestNormalize)
library(sjmisc)
library(chemometrics)
library(compositions)
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)
#RColorBrewer
display.brewer.all()
display.brewer.pal(11,"BrBG")
display.brewer.all(colorblindFriendly = TRUE)
# Show BrBG colour palette with 11 colours & get colour codes to copy  --------
nb.cols <- 11
PiYG1 <- colorRampPalette(brewer.pal(11, "PiYG"))(nb.cols)
PiYG1
nb.cols <- 9
Greys1 <- colorRampPalette(brewer.pal(9, "Greys"))(nb.cols)
Greys1
nb.cols <- 9
Greens1 <- colorRampPalette(brewer.pal(9, "Greens"))(nb.cols)
Greens1
nb.cols <- 9
Blues1 <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)
Blues1

# PARTS -------------------------------------------------------------------

# PART 0 - Mass% to Elemental conversion & clr for core and reference data
# PART 1 - Ardley Lake - ARD-ITRAX-SH2
# PART 2 - Yanou Lake - YAN-ITRAX-SH20
# PART 3 - Calibration - Matched ARD & YAN ITRAX and XRF 1 cm dataset
# PART 4 - ITRAX & XRF comparison vs Age plots
# PART 5 - Linear regression models 1-5 (Ca, P, Cu, Zn, Sr)

# -------------------------------------------------------------------------

# Data and code are located here:
# https://github.com/stever60/Bertrand_2022 

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Bertrand_2022")
#check working directory
getwd() 

# -------------------------------------------------------------------------

#  PART 0: Convert subsample XRF mass oxide to elemental data  STILL TO FINISH -------------------------------------------------------

# Location of ARD, Yan and SSI subsample XRF data in Dropbox - add figure doi when done 
# https://github.com/stever60/Bertrand_2022/tree/main/Data

# Import subsample XRF data as oxide mass% and convert to elemental molecular weight
# for use later on in correlations etc.



# Oxide to elemental conversion factors
# Datafile: https://github.com/stever60/Bertrand_2022/tree/main/Data
ecf <- read_csv("Data/Oxide_elemental_conv_factors.csv")
ecf

# Ardley Lake - mass oxide to elemental conversion
# Datafile: https://github.com/stever60/Bertrand_2022/tree/main/Data
ARD_mass <- read_csv("Data/ARD_mass.csv")
ARD_mass

ARD_elemental <- ARD_mass %>% 
  mutate(across(c("SiO2"), .fns = ~.*ecf$Si)) %>% rename(Si = SiO2) %>% 
  mutate(across(c("Al2O3"), .fns = ~.*ecf$Al)) %>% rename(Al = Al2O3) %>%
  mutate(across(c("K2O"), .fns = ~.*ecf$K)) %>% rename(K = K2O) %>%
  mutate(across(c("CaO"), .fns = ~.*ecf$Ca)) %>% rename(Ca = CaO) %>%
  mutate(across(c("TiO2"), .fns = ~.*ecf$Ti)) %>% rename(Ti = TiO2) %>%
  mutate(across(c("MnO"), .fns = ~.*ecf$Mn)) %>% rename(Mn = MnO) %>%
  mutate(across(c("FeO"), .fns = ~.*ecf$Fe)) %>% rename(Fe = FeO) %>%
  mutate(across(c("Na2O"), .fns = ~.*ecf$Na)) %>% rename(Na = Na2O) %>%
  mutate(across(c("MgO"), .fns = ~.*ecf$Mg)) %>% rename(Mg = MgO) %>% 
  mutate(across(c("P2O5"), .fns = ~.*ecf$P)) %>% rename(P = P2O5)
  #mutate(across(c("Cr2O3"), .fns = ~.*ecf$Cr)) %>% rename(Cr = Cr2O3) %>%
  #mutate(across(c("NiO"), .fns = ~.*ecf$Ni)) %>% rename(Ni = NiO) %>%
  #mutate(across(c("SO2"), .fns = ~.*ecf$S)) %>% rename(S = SO2)
ARD_elemental 
write.csv(ARD_elemental,"Output/ARD_elemental.csv", row.names = FALSE)

# Yanou Lake - mass oxide to elemental conversion
# Datafile: https://github.com/stever60/Bertrand_2022/tree/main/Data
YAN_mass <- read_csv("Data/YAN_mass.csv")
YAN_mass

YAN_elemental <- YAN_mass %>% 
  mutate(across(c("SiO2"), .fns = ~.*ecf$Si)) %>% rename(Si = SiO2) %>%
  mutate(across(c("Al2O3"), .fns = ~.*ecf$Al)) %>% rename(Al = Al2O3) %>%
  mutate(across(c("K2O"), .fns = ~.*ecf$K)) %>% rename(K = K2O) %>%
  mutate(across(c("CaO"), .fns = ~.*ecf$Ca)) %>% rename(Ca = CaO) %>%
  mutate(across(c("TiO2"), .fns = ~.*ecf$Ti)) %>% rename(Ti = TiO2) %>%
  mutate(across(c("MnO"), .fns = ~.*ecf$Mn)) %>% rename(Mn = MnO) %>%
  mutate(across(c("FeO"), .fns = ~.*ecf$Fe)) %>% rename(Fe = FeO) %>%
  mutate(across(c("Na2O"), .fns = ~.*ecf$Na)) %>% rename(Na = Na2O) %>%
  mutate(across(c("MgO"), .fns = ~.*ecf$Mg)) %>% rename(Mg = MgO) %>%
  mutate(across(c("P2O5"), .fns = ~.*ecf$P)) %>% rename(P = P2O5) %>%
#mutate(across(c("Cr2O3"), .fns = ~.*ecf$Cr)) %>% rename(Cr = Cr2O3) %>%
#mutate(across(c("NiO"), .fns = ~.*ecf$Ni)) %>% rename(Ni = NiO) %>%
  mutate(across(c("SO2"), .fns = ~.*ecf$S)) %>% rename(S = SO2)
YAN_elemental 
write.csv(YAN_elemental,"Output/YAN_elemental.csv", row.names = FALSE)

# Fildes and Ardley reference data - mass oxide to elemental conversion
SSI_ref_mass <- read_csv("Data/SSI_ref.csv")

SSI_ref_elemental <- SSI_ref_mass %>% 
  rename(FeO = FeOT) %>% 
  mutate(across(c("SiO2"), .fns = ~.*ecf$Si)) %>% rename(Si = SiO2) %>% 
  mutate(across(c("Al2O3"), .fns = ~.*ecf$Al)) %>% rename(Al = Al2O3) %>%
  mutate(across(c("K2O"), .fns = ~.*ecf$K)) %>% rename(K = K2O) %>%
  mutate(across(c("CaO"), .fns = ~.*ecf$Ca)) %>% rename(Ca = CaO) %>%
  mutate(across(c("TiO2"), .fns = ~.*ecf$Ti)) %>% rename(Ti = TiO2) %>%
  mutate(across(c("MnO"), .fns = ~.*ecf$Mn)) %>% rename(Mn = MnO) %>%
  mutate(across(c("FeO"), .fns = ~.*ecf$Fe)) %>% rename(Fe = FeO) %>%
  mutate(across(c("Na2O"), .fns = ~.*ecf$Na)) %>% rename(Na = Na2O) %>%
  mutate(across(c("MgO"), .fns = ~.*ecf$Mg)) %>% rename(Mg = MgO) %>% 
  mutate(across(c("P2O5"), .fns = ~.*ecf$P)) %>% rename(P = P2O5) 
  #mutate(across(c("Cr2O3"), .fns = ~.*ecf$Cr)) %>% rename(Cr = Cr2O3) %>%
  #mutate(across(c("NiO"), .fns = ~.*ecf$Ni)) %>% rename(Ni = NiO) %>%
  #mutate(across(c("SO2"), .fns = ~.*ecf$S)) %>% rename(S = SO2)
SSI_ref_elemental
tail(SSI_ref_elemental)
write.csv(SSI_ref_elemental,"Output/SSI_ref_elemental.csv", row.names = FALSE)

# Centered log ratio (clr) conversion for all three datasets --------------------------------------------
library(compositions)

# SSI ref elemental data ------------------------------------------

# remove columns with any NAs - will remove Si, Ti, LOI, TC, TN
SSI_ref_clr <- SSI_ref_elemental %>% 
  select_if(~!any(is.na(.))) %>% 
  select(Al:P) %>% 
  clr()
SSI_ref_clr <- as_tibble(SSI_ref_clr)
SSI_ref_clr
write.csv(SSI_ref_clr,"Output/SSI_ref_clr.csv", row.names = FALSE)

# How to remove columns with only NAs - will remove TC and TN 
#SSI_ref_elemental1 <- SSI_ref_elemental[, colSums(is.na(SSI_ref_elemental)) != nrow(SSI_ref_elemental)]
#tail(SSI_ref_elemental1)

# How to set threshold for removing columns with NAs in them based on number or percentage of columns containing NA columns
# e.g., remove whole columns where eg >50% or 10% of data in columns are NA – 0.1 will remove all NA columns 
#SSI_ref_elemental1 <- SSI_ref_elemental[, colSums(is.na(SSI_ref_elemental)) < nrow(SSI_ref_elemental) * 0.1]
#tail(SSI_ref_elemental1)

# How to remove all rows with NAs -  this leaves nothing as all rows have NAs in them
#SSI_ref_elemental1 <- SSI_ref_elemental %>% 
#  filter(!if_any(everything(), is.na))
#SSI_ref_elemental1

# Merge clr dataframe with Sample_code and location dataframe columns
SSI_ref_elemental1 <- SSI_ref_elemental %>% 
  select(Sample_code: Guano_Sum, Location: Reference, Location, Location1, Location2, Reference)
SSI_ref_elemental1

SSI_ref_elemental_clr <- bind_cols(SSI_ref_elemental1, SSI_ref_clr) %>% 
  relocate(Al:P, .before = Location)
SSI_ref_elemental_clr

write.csv(SSI_ref_elemental_clr,"Output/SSI_ref_elemental_clr.csv", row.names = FALSE)

# ARD elemental data ----------------------------------------------------------------

# remove columns with any NAs - will remove Si, Ti, LOI, TC, TN
ARD_clr <- ARD_elemental %>% 
  select(Al:P) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
ARD_clr <- as_tibble(ARD_clr)
ARD_clr 
write.csv(ARD_clr,"Output/ARD_clr.csv", row.names = FALSE)

# Merge clr dataframe with Sample_code and location dataframe columns
ARD_elemental1 <- ARD_elemental %>% 
  select(Record: Guano_Sum, Location: Reference)
ARD_elemental1

ARD_elemental_clr <- bind_cols(ARD_elemental1, ARD_clr) %>% 
  relocate(Al:P, .before = Location) 
ARD_elemental_clr
write.csv(ARD_elemental_clr,"Output/ARD_elemental_clr.csv", row.names = FALSE)

# YAN elemental data ----------------------------------------------------------------

# remove columns with any NAs - will remove Si, Ti, LOI, TC, TN
YAN_clr <- YAN_elemental %>% 
  select(Al:P) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
YAN_clr <- as_tibble(YAN_clr)
YAN_clr 
write.csv(YAN_clr,"Output/YAN_clr.csv", row.names = FALSE)

# Merge clr dataframe with Sample_code and location dataframe columns
YAN_elemental1 <- YAN_elemental %>% 
  select(Record: Guano_Sum, Location: Reference)
YAN_elemental1

YAN_elemental_clr <- bind_cols(YAN_elemental1, YAN_clr) %>% 
  relocate(Al:P, .before = Location)
YAN_elemental_clr
write.csv(YAN_elemental_clr,"Output/YAN_elemental_clr.csv", row.names = FALSE)

# Merging the clr datasets for plotting ---------------------------------------
ARD_elemental_clr2 <- ARD_elemental_clr %>% 
  mutate(Depth_cm = Strat_Depth_cm*1) %>% 
  mutate(Record1 = Record) %>%
  relocate(Depth_cm, .before = Strat_Depth_cm) %>% 
  relocate(Record1, .before = Section) %>% 
  unite("Sample_code", Record, Strat_Depth_cm, sep = "_") %>% 
  select(-c(SH20_95CI_min_age, SH20_95CI_max_age, SH20_median_age))
ARD_elemental_clr2

YAN_elemental_clr2 <- YAN_elemental_clr %>% 
  mutate(Depth_cm = Strat_Depth_cm*1) %>% 
  mutate(Record1 = Record) %>%
  relocate(Depth_cm, .before = Strat_Depth_cm) %>% 
  relocate(Record1, .before = Section) %>% 
  unite("Sample_code", Record, Strat_Depth_cm, sep = "_") %>% 
  select(-c(SH20_95CI_min_age, SH20_95CI_max_age, SH20_median_age))
YAN_elemental_clr2

SSI_ref_elemental_clr2 <- SSI_ref_elemental_clr %>% 
  rename(Record1 = Record)
SSI_ref_elemental_clr2

ARD_YAN_SSI_clr <- bind_rows(ARD_elemental_clr2,YAN_elemental_clr2,SSI_ref_elemental_clr2) %>% 
  rename(Record = Record1)
head(ARD_YAN_SSI_clr)
tail(ARD_YAN_SSI_clr)
write.csv(ARD_YAN_SSI_clr,"Output/ARD_YAN_SSI_elemental_clr.csv", row.names = FALSE)


# Bi-plots and correlation ----------------------------------------------------------------
library(ggpubr)

# Show BrBG colour palette with 5 colours & get colour codes to copy  --------
display.brewer.pal(6,"BrBG")
nb.cols <- 6
BrBG1 <- colorRampPalette(brewer.pal(6, "BrBG"))(nb.cols)
BrBG1
BrBG2 <- c("#A6611A", "#DFC27D", "#018571")
BrBG3 <- c("#8C510A", "#D8B365", "#F6E8C3", "#C7EAE5", "#5AB4AC", "#01665E")

# formula1 <- y ~ poly(x, 1, raw = TRUE)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

# correlation all data 
p1 <- ggplot(ARD_YAN_SSI_clr, aes(x=Al, y=Ca)) +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE) +
  stat_cor(label.y = c(2, 1.5, 1, 0.5, 0, -0.5), label.x = c(2, 2, 2, 2, 2)) + 
  #stat_regline_equation(label.y = c(68), label.x = c(20)) +
  scale_shape_manual(values = c(21, 22, 23, 21, 22, 23)) +
  scale_fill_manual(values = BrBG3) +
  scale_color_manual(values = BrBG3) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL)) +
  labs(y=expression('clr(P)'), x=expression('clr(Al)'))

p1

ggsave("Figures/clr/Fig 1_clr_XRFelemental_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# correlation by record
p2 <- ggplot(ARD_YAN_SSI_clr, aes(x=Al, y=Ca, color=Record)) +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, aes(group=Record)) +
  stat_cor(label.y = c(2, 1.5, 1, 0.5, 0, -0.5), label.x = c(2, 2, 2, 2, 2)) + 
  #stat_regline_equation(label.y = c(68), label.x = c(20)) +
  scale_shape_manual(values = c(21, 22, 23, 21, 22, 23)) +
  scale_fill_manual(values = BrBG3) +
  scale_color_manual(values = BrBG3) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL)) +
  labs(y=expression('clr(P)'), x=expression('clr(Al)'))

p2

ggsave("Figures/clr/Fig 1_clr_XRFelemental_byrecord.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# PART 1 - Ardley Lake - ARD-ITRAX-SH20  ------------------------------------------------------------------

# ARD-ITRAX-SH20 - Import cps data ---------------------------------------------
library(dplyr)

# Location of ARD-ITRAX-SH20 composite dataset
# https://github.com/stever60/Bertrand_2022/tree/main/Data 

# Composite ARD - ITRAX cps data & Sh20 BACON age-depth model (added in Excel)
ARD_x.raw <- read_csv("Data/ARD_ITRAX_COMP_SH20.csv")
ARD_x.raw

# Filter data to keep only valid data and remove kcps <-2sd and MSE >2sd ----------------------------
ARD_kcps.mean <- mean(ARD_x.raw$kcps)
ARD_kcps.sd <- 2*sd(ARD_x.raw$kcps)
ARD_kcps.thres <- ARD_kcps.mean - ARD_kcps.sd 
ARD_kcps.thres

ARD_MSE.mean <- mean(ARD_x.raw$MSE)
ARD_MSE.sd <- 2*sd(ARD_x.raw$MSE)
ARD_MSE.thres <- ARD_MSE.mean + ARD_MSE.sd 
ARD_MSE.thres

# filter by validity and then remove analyses above MSE threshold and below kcps - 2s & 
ARD_COMP <- filter(ARD_x.raw, validity == "1") %>% 
  filter(MSE < ARD_MSE.thres) %>%
  filter(kcps > ARD_kcps.thres)

# Remove unwanted columns - might need to change these - and write to file
ARD_COMP_filter <- select(ARD_COMP, 
                          -c(filename, position_mm, 
                             core_depth_cm, min, max,
                             median, mean_2s_range:Foffset,
                             D1, S1, S2, S3)) %>% 
  rename(SH20_age = mean)

write.csv(ARD_COMP_filter,"Output/ARD/1_cps_filter.csv", row.names = FALSE)

# Define elements lists -------------------------------------------------
library(PeriodicTable)

# list of all possible elements
elements <- c(symb(1:117), "Mo_inc", "Mo_coh", "TS_sum", "cps_sum")

# list of all elements from Q-spec matching
ARD_elements <- select(ARD_COMP_filter, c(Mg:Mo_coh)) %>% 
  names()
ARD_elements

# Ar (tube gas), Ta & W (splutter)
machine_elements <- select(ARD_COMP_filter, c(Ar, Ta, W)) %>% 
  names()

# REE 
REE <- select(ARD_COMP_filter, c(La:Ho)) %>% 
  names()
REE

# Machine generated elements Ar (tube gas), Ta & W (splutter), REE removed
# ARD_elements1 <- select(ARD_COMP_filter, 
#                        -c(Core:MSE, all_of(machine_elements), 
#                           Nb:Cs, REE, Ir, Pt)) %>% names()
# ARD_elements1

# Calculate as % of normalising factors TS (Total Scatter), cps_sum & inc/coh -------------------------------------------------------

# Add to columns - find a way to replace [9:67] with column headings as "Mg":"Mo_coh"
ARD_COMP_filter

# Type in row numbers of elements and scatter to calculate  across 
ARD_rowsums <- 9:67

ARD_COMP_filter1 <- ARD_COMP_filter %>% 
  replace(is.na(.), 0) %>%
  mutate(TS_sum = Mo_inc + Mo_coh) %>% 
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  mutate(coh_inc = Mo_coh / Mo_inc) %>%
  mutate(cps_sum = rowSums(.[ARD_rowsums]))
  #mutate(cps_sum = row_sums(ARD_elements))

ARD_COMP_filter1

# Standardise and centre cps data  - using column names --------------------------------------

# list of element column names for plotting
ARD_col <- select(ARD_COMP_filter1, c(Mg:inc_coh)) %>% names()
ARD_col

# list of all element and other parameter column names 
ARD_col1 <- select(ARD_COMP_filter1, c(kcps:cps_sum)) %>% names()
ARD_col1

ARD_COMP_filter1.Z <- ARD_COMP_filter1
ARD_COMP_filter1.Z[, ARD_col] <- scale(ARD_COMP_filter1[, ARD_col], center = TRUE, scale = TRUE)
ARD_COMP_filter1.Z

# Convert cps and cps Z-scores to long format for plotting -------------------------------------------------

ARD_COMP_filter1_long <- select(ARD_COMP_filter1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
ARD_COMP_filter1_long 

ARD_COMP_filter1_long.Z <- select(ARD_COMP_filter1.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
ARD_COMP_filter1_long.Z 


# Generate cps summary stats  --------------------------------------------
library(psych)

ARD_summary <- ARD_COMP_filter1 %>%
  select(all_of(ARD_col1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter cps dataset and stats to file --------------------------------------------------
write.csv(ARD_COMP_filter1,"Output/ARD/2_cps_filter1.csv", row.names = FALSE)
write.csv(ARD_summary,"Output/ARD/2.1_cps_filter1_stats.csv", row.names = FALSE)
write.csv(ARD_COMP_filter1.Z,"Output/ARD/3_cps_filter1_Z.csv", row.names = FALSE)

# Filter cps element column list based on mean and max cps values -------------------

ARD_mean50 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 50)
}
ARD_mean200 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 200)
}
ARD_max50 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 50)
}

ARD_m50 <- select(ARD_COMP_filter1, 
                  Core, depth_cm, SH20_age, kcps, MSE, 
                  Mg:Mo_coh & where(ARD_mean50)) %>% print()
ARD_m200 <- select(ARD_COMP_filter1, 
                   Core, depth_cm, SH20_age, kcps, MSE, 
                   Mg:Mo_coh & where(ARD_mean200)) %>% print()
ARD_mx50 <- select(ARD_COMP_filter1, 
                   Core, depth_cm, SH20_age, kcps, MSE, 
                   Mg:Mo_coh & where(ARD_max50), TS_sum:cps_sum) %>% print()

# Create lists of column names to take forward
ARD_col_m50 <- select(ARD_COMP_filter1, Mg:Mo_coh 
                      & -REE & -machine_elements 
                      & where(ARD_mean50)) %>% names() %>% print()
ARD_col_m200 <- select(ARD_COMP_filter1, Mg:Mo_coh 
                       & -REE & -machine_elements 
                       & where(ARD_mean200)) %>% names() %>% print()
ARD_col_m50_mx50 <- select(ARD_COMP_filter1, Mg:Mo_coh 
                        & -REE & -machine_elements 
                        & where(ARD_max50) 
                        & where(ARD_mean50)) %>% names() %>% print()

# Correlation matrices -------------------------------------------------------

library(GGally)
library(dplyr)
# Correlation plot with max200 elements - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=2))

ggcorr(ARD_COMP_filter1[, ARD_col1], method = c("everything", "pearson"), 
       size = 2, label = FALSE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/ARD/Fig 0_Corr_matrix_cps_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Manual list of elements with significant correlations - determined visually from correlation matrix
plot_elements1_ARD <-c("Si", "P", "S", "K", "Ca", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Cu", "Zn",
                       "Br", "Rb", "Sr", "Zr", "Mo", "Ba", "Pb", "Bi","Mo_inc", "Mo_coh", "inc_coh", "coh_inc")

# Stats generated list of elements to take forward based on m50
plot_elements2_ARD <- c("Si", "P", ARD_col_m50, "inc_coh", "coh_inc")

ggcorr(ARD_COMP_filter1[,  plot_elements2_ARD], method = c("everything", "pearson"), 
       size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/ARD/Fig 0_Corr_matrix_cps_m50.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
theme_set(theme_bw(base_size=8))
ggpairs(ARD_COMP_filter1, columns = plot_elements2_ARD, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Figures/ARD/Fig 0__Corr-den_matrix_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# Summary cps plots vs depth --------------------------------------------------------

library(tidypaleo)

# Figure 1 - cps elements1
theme_set(theme_bw(8))
ARD_Fig1 <- ARD_COMP_filter1_long  %>%
  filter(param %in% plot_elements1_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements1_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Signifcantly correlated elements cps summary")
ARD_Fig1
ggsave("Figures/ARD/Fig 1_cps_elements1.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 2 - cps elements2 based on mean>50 cps stats
theme_set(theme_bw(8))
ARD_Fig2 <- ARD_COMP_filter1_long  %>%
  filter(param %in% plot_elements2_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements2_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps based on mean>50 cps stats")
ARD_Fig2
ggsave("Figures/ARD/Fig 2_cps_elements2.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 3 - cps elements2 Z-scores
theme_set(theme_bw(8))
ARD_Fig3 <- ARD_COMP_filter1_long.Z  %>%
  filter(param %in% plot_elements2_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements2_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps Z-scores based on mean>50 cps stats")
ARD_Fig3
ggsave("Figures/ARD/Fig 3_cps_elements2_Z.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Create a cps centered log ratio (clr) elemental datafile and plot  -------

# load compositions package
library(compositions)
ARD_COMP_filter1_text <- select(ARD_COMP_filter1, Core:MSE, inc_coh, coh_inc)
ARD_COMP_filter1_text

# look at element list to use
plot_elements2_ARD

# remove scatter ratios to leave only subset of elements and scatter parameters i.e., closed sum measurements 
# can change elements used in clr here!
ARD_clr_elements2 <- purrr::discard(plot_elements2_ARD,.p = ~stringr::str_detect(.x,"inc_coh"))
ARD_clr_elements2 <- purrr::discard(ARD_clr_elements2 ,.p = ~stringr::str_detect(.x,"coh_inc"))
ARD_clr_elements2

# apply clr to cps data and convert all 0 to NA for plotting clarity
ARD_COMP_filter1_clr1 <- select(ARD_COMP_filter1, all_of(ARD_clr_elements2)) %>% 
  select(all_of(ARD_clr_elements2)) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
ARD_COMP_filter1_clr1 <- as_tibble(ARD_COMP_filter1_clr1) %>% 
  na_if(0)
ARD_COMP_filter1_clr1

ARD_COMP_filter1_clr <- bind_cols(ARD_COMP_filter1_text, ARD_COMP_filter1_clr1) %>% 
  relocate(inc_coh, .after = Mo_coh) %>% 
  relocate(coh_inc, .after = inc_coh)
ARD_COMP_filter1_clr

write.csv(ARD_COMP_filter1_clr,"Output/ARD/3.1_cps_filter1_clr.csv", row.names = FALSE)

# convert clr datafile to long format for plotting
ARD_COMP_filter1_clr_long <- select(ARD_COMP_filter1_clr,  Core, depth_cm, SH20_age, 
                                    kcps, MSE, inc_coh, coh_inc, all_of(plot_elements2_ARD)) %>%
  pivot_longer(all_of(plot_elements2_ARD), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
head(ARD_COMP_filter1_clr_long)
tail(ARD_COMP_filter1_clr_long)

# Figure 3a - cps clr for clr_elements
theme_set(theme_bw(8))
ARD_Fig3A <- ARD_COMP_filter1_clr_long  %>%
  filter(param %in% plot_elements2_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements2_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps clr based on mean>50 cps stats")
ARD_Fig3A
ggsave("Figures/ARD/Fig 3A_cps_clr_elements.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Normalization (elements and scatter parameters) of cps data  ---------------------------------------------------------
  
# Inc normalised
ARD_inc_norm <- ARD_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Mo_inc))
# coh/inc normalised - Boyle 2015
ARD_coh_inc_norm <- ARD_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ (Mo_coh/Mo_inc)))
# Total scatter (TS) normalised - Kylander et al 2011/2012
ARD_TS_norm <- ARD_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ TS_sum))
# cps_sum normalised
ARD_cps_sum_norm <- ARD_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ cps_sum))
# Ti normalised
ARD_Ti_norm <- ARD_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Ti))
# Natural log Ti normalized & replace -Inf with NA --------
ARD_Ln_Ti_norm <- ARD_Ti_norm  %>% 
  mutate(across(c("Mg":"Mo_coh"),log)) 
is.na(ARD_Ln_Ti_norm)<-sapply(ARD_Ln_Ti_norm, is.infinite)
# Replace NA with 0 - if needed
# ARD_Ln_Ti_norm[is.na(ARD_Ln_Ti_norm)]<-0

#Check = 1
ARD_inc_norm$Mo_inc
ARD_Ti_norm$Ti
#Check = 0
ARD_Ln_Ti_norm$Ti

# Standardize and center Ti-normalized data  --------------------------------------
ARD_Ln_Ti_norm.Z <- ARD_Ln_Ti_norm
ARD_Ln_Ti_norm.Z[, ARD_col1] <- scale(ARD_Ln_Ti_norm[, ARD_col1], center = TRUE, scale = TRUE)
ARD_Ln_Ti_norm.Z
# can replace Ti normalized with other normalisation parameters above

# write ARD_Ln_Ti_norm.Z to file for Part 3 & 4
write.csv(ARD_Ln_Ti_norm.Z,"Output/ARD/ARD_Ln_Ti_norm.Z.csv", row.names = FALSE)

# Elements and scatter as %cps sum - Bertrand et al 2021 -------- produces the same result as %TSN sum

ARD_cps_sum_norm_pc <- ARD_cps_sum_norm %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(cps_pc_sum = rowSums(.[ARD_rowsums]))

#Check everything = 100
head(ARD_cps_sum_norm_pc$cps_pc_sum)

# define new element list to include %cps sum (check=100 in output)
ARD_col3 <- select(ARD_cps_sum_norm_pc, c(Mg:cps_pc_sum)) %>% 
  names()
ARD_col3

# Create %cps_sum summary stats table
ARD_cps_sum_pc_summary <- ARD_cps_sum_norm_pc %>%
  select(all_of(ARD_col3)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Elements and scatter as percentage of TSN sum - Roberts et al 2017 --------- produces the same result as %cps sum
ARD_TSN_sum <- ARD_TS_norm %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_sum = rowSums(.[ARD_rowsums]))
ARD_TSN_sum

ARD_TSN_pc <- ARD_TSN_sum %>% mutate(across(c("Mg":"Mo_coh"),
                                         .fns = ~./TSN_sum)) %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum = rowSums(.[ARD_rowsums]))
ARD_TSN_pc

#Check everything = 100
head(ARD_TSN_pc$TSN_pc_sum)

# define new element list to include %TSN sum (check=100 in output)
ARD_col4 <- select(ARD_TSN_pc, c(Mg:TSN_pc_sum)) %>% 
  names()
ARD_col4

# Create TSN stats table
ARD_TSN_pc_summary <- ARD_TSN_pc %>%
  select(all_of(ARD_col4)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Filter elements based on %cps_sum (%TSN) max > 0.5 and mean >0.1 - *** replace %TSN with %cps sum here *** ----------------
ARD_max0.5 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 0.5)
}
ARD_mean0.1 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 0.1)
}

ARD_mx0.5 <- select(ARD_TSN_pc, Mg:Mo_coh & where(ARD_max0.5)) %>% print()
ARD_m0.1 <- select(ARD_TSN_pc, Mg:Mo_coh & where(ARD_mean0.1)) %>% print()

# define new element list of filtered elements and their %TSN sum (check=100 in output)
ARD_TSN_mx0.5 <- select(ARD_TSN_pc, Mg:Mo_coh 
                      & -REE & -machine_elements 
                      & where(ARD_max0.5)) %>% names()

ARD_TSN_m0.1 <- select(ARD_TSN_pc, Mg:Mo_coh 
                             & -REE & -machine_elements 
                             & where(ARD_mean0.1)) %>% names()
ARD_TSN_mx0.5
ARD_TSN_m0.1

# New element list based on  ARD_col_m0.1 + Si, P , S
ARD_TSN_m0.1_list <- c("Si", "P", "S", ARD_TSN_m0.1)
ARD_TSN_m0.1_list

# Add filtered element TSN sum to TSN_pc dataset & rename
ARD_TSN_pc1 <- ARD_TSN_pc %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum1 = rowSums(.[ARD_TSN_m0.1_list]))
ARD_TSN_pc1
# check filtered element total <100 but >95%
head(ARD_TSN_pc1$TSN_pc_sum1)

# Create TSN stats table with filtered element %TSN Sum added to end 
ARD_col5 <- select(ARD_TSN_pc1, Mg:TSN_pc_sum1) %>% names()

ARD_TSN_pc1_summary <- ARD_TSN_pc1 %>%
  select(all_of(ARD_col5)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write normalised datasets, %cps & %TSN and stats summaries to file --------------------------------------------------
write.csv(ARD_Ti_norm,"Output/ARD/4_Ti_norm.csv", row.names = FALSE)
write.csv(ARD_inc_norm,"Output/ARD/5_inc_norm.csv", row.names = FALSE)
write.csv(ARD_TS_norm,"Output/ARD/6_TS_norm.csv", row.names = FALSE)
write.csv(ARD_cps_sum_norm,"Output/ARD/7_cps_sum_norm.csv", row.names = FALSE)
write.csv(ARD_TSN_pc,"Output/ARD/8_TSN_pc.csv", row.names = FALSE)
write.csv(ARD_TSN_pc_summary,"Output/ARD/8.1_TSN_pc_stats.csv", row.names = FALSE)
write.csv(ARD_TSN_pc1,"Output/ARD/8.2_TSN_pc1.csv", row.names = FALSE)
write.csv(ARD_TSN_pc1_summary,"Output/ARD/8.3_TSN_pc1_stats.csv", row.names = FALSE)
write.csv(ARD_cps_sum_norm_pc,"Output/ARD/9_cps_sum_pc.csv", row.names = FALSE)
write.csv(ARD_cps_sum_pc_summary,"Output/ARD/9.1_cps_sum_pc_stats.csv", row.names = FALSE)
write.csv(ARD_Ln_Ti_norm,"Output/ARD/10_Ln_Ti_norm.csv", row.names = FALSE)
write.csv(ARD_Ln_Ti_norm.Z,"Output/ARD/11_Ln_Ti_norm_Z.csv", row.names = FALSE)

# Convert to long format ------------------------------------------------
ARD_inc_norm_long <- select(ARD_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
ARD_inc_norm_long

ARD_coh_inc_norm_long <- select(ARD_coh_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
ARD_coh_inc_norm_long

ARD_Ln_Ti_norm_long <- select(ARD_Ln_Ti_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
ARD_Ln_Ti_norm_long

ARD_Ln_Ti_norm.Z_long <- select(ARD_Ln_Ti_norm.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col1)) %>%
  pivot_longer(all_of(ARD_col1), names_to = "param", values_to = "value")
ARD_Ln_Ti_norm.Z_long

ARD_cps_sum_norm_pc_long <- select(ARD_cps_sum_norm_pc,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col3)) %>%
  pivot_longer(all_of(ARD_col3), names_to = "param", values_to = "value")
ARD_cps_sum_norm_pc_long

ARD_TSN_pc1_long <- select(ARD_TSN_pc1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(ARD_col5)) %>%
  pivot_longer(all_of(ARD_col5), names_to = "param", values_to = "value")
ARD_TSN_pc1_long

# Summary normalised plots--------------------------------------

# # Figure 4 cps/inc. -----------------------------------------------------

# manually user defined
plot_elements3_ARD <- c("K", "Ca", "Ti", "Mn", "Fe", "Cu", "Zn", 
                        "Br", "Sr", "Zr", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")
# OR 

# stats defined based on mean %TSN >0.1% - %TSN output is the same as %cps sum output
plot_elements3_ARD <- c(ARD_TSN_m0.1_list, "inc_coh", "coh_inc")

# can replace %TSN_M0.1 element list with cps ARD_col_m50 element list
# plot_elements3_ARD <- c(ARD_col_m50, "inc_coh", "coh_inc")

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
# Roberts et al 2017 - CONISS with broken stick Hellingers dist defined zones - red dashed lines
theme_set(theme_bw(8))
ARD_Fig4a <- ARD_inc_norm_long  %>%
  filter(param %in% plot_elements3_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements3_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc.", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Incoherent (inc.) scatter normalised (filtered elements): CONISS XRF")

ARD_Fig4b <- ARD_coh_inc_norm_long  %>%
  filter(param %in% plot_elements3_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements3_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/(coh./inc.)", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Coherent/Incoherent (coh./inc.) scatter ratio normalised (filtered elements): CONISS XRF")
ggarrange(ARD_Fig4a, ARD_Fig4b, nrow = 2)
ggsave("Figures/ARD/Fig 4__cps_inc&coh_inc.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 5 %TSN----------------------------------------------------------------

# stats defined 
plot_elements4_ARD <- c(ARD_TSN_m0.1_list, "inc_coh", "coh_inc")

# OR

# as m0.1 defined list but with noisy elements eg Si, P, S removed
plot_elements4_ARD <- c("K", "Ca", "Ti", "Mn", "Fe", "Cu", "Zn", 
                        "Br", "Sr", "Zr","Mo_inc", "Mo_coh", "inc_coh", "coh_inc")

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
# Roberts et al 2017 - CONISS with broken stick Hellingers dist defined zones - red dashed lines
theme_set(theme_bw(8))
ARD_Fig5 <- ARD_TSN_pc1_long  %>%
  filter(param %in% plot_elements4_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements4_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%TSN", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%TSN (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss1_ARD <- ARD_TSN_pc1_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 6 - %TSN with CONISS
ARD_Fig6 <- ARD_Fig5 +
  layer_dendrogram(coniss1_ARD, aes(y = depth_cm), param = "CONISS") +
  layer_zone_boundaries(coniss1_ARD, aes(y = depth_cm, col = "blue", lty = 2, alpha = 0.7)) +
  ggtitle("%TSN (filtered elements): CONISS XRF(red) ITRAX (black)")

ggarrange(ARD_Fig5, ARD_Fig6, nrow = 2)
ggsave("Figures/ARD/Fig 5_6_TSN_pc_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

#Figure 7 - TSN with CONISS and Roberts et al (2017) zone boundaries
#ARD_Fig6 +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%TSN (filtered elements): Coniss zone comparison")
#ggsave("Figures/ARD/Fig 7_TSN_pc_CONISS_comparison.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 8  %cps_sum
#define colour scheme for 6 groups (5 guano and 1 non-guano)
guano_zone_colours <- c("#BDBDBD", "#00441B", "#00441B", "#006D2C", "#238B45", "#74C476")
theme_set(theme_bw(8))
ARD_Fig8 <- ARD_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements3_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements3_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm, aes(colour = guano_zone_colours))) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%cps_sum", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%cps sum (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss2_ARD <- ARD_cps_sum_norm_pc_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 9 %cps_sum with CONISS
ARD_Fig9 <- ARD_Fig8 +
  layer_dendrogram(coniss2_ARD, aes(y = depth_cm), param = "CONISS") +
  layer_zone_boundaries(coniss2_ARD, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
  ggtitle("%cps sum (filtered elements): CONISS XRF(red) ITRAX (black)")
ggarrange(ARD_Fig8, ARD_Fig9, nrow = 2)
ggsave("Figures/ARD/Fig 8_9_cps_sum_pc_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 10 - compare to Roberts et al. 2017 CONISS zoning 
#ARD_Fig9 +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%cps sum (filtered elements): CONISS comparison")
#ggsave("Figures/ARD/Fig 10_cps_sum_pc_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")


# create a new element list without Ti
plot_elements5_ARD <- purrr::discard(plot_elements4_ARD,.p = ~stringr::str_detect(.x,"Ti"))
plot_elements5_ARD

# Figure 11 Ti normalised
theme_set(theme_bw(8))
ARD_Fig11 <- ARD_Ln_Ti_norm_long  %>%
  filter(param %in% plot_elements5_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements5_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Ln(cps/Ti)", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised (filtered elements): CONISS XRF")
ARD_Fig11
ggsave("Figures/ARD/Fig 11_Ti_norm.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 12 - Ti normalised as Z-scores 
theme_set(theme_bw(8))
ARD_Fig12 <- ARD_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements5_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements5_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "LN(cps/Ti) Z-scores", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")

# Add CONISS Zones using in built tidypalaeo CONISS 
coniss3_ARD <- ARD_Ln_Ti_norm.Z_long %>%
  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

# Figure 13
ARD_Fig13 <- ARD_Fig12 + labs(x = "Ln(cps/Ti Z-scores)", y = "Depth (cm)") +
  layer_dendrogram(coniss3_ARD, aes(y = depth_cm), param = "CONISS ITRAX") +
  layer_zone_boundaries(coniss3_ARD, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF(red) ITRAX (black)")
ggarrange(ARD_Fig12, ARD_Fig13, nrow = 2)
ggsave("Figures/ARD/Fig 12_13_Ti_Z_CONISS.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 14
#ARD_Fig12 + labs(x = "cps/Ti Z-scores", y = "Depth (cm)") +
#  layer_dendrogram(coniss3_ARD, aes(y = depth_cm), param = "CONISS") +
#  layer_zone_boundaries(coniss3_ARD, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
#  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS comparison")
#ggsave("Figures/ARD/Fig 14_Ti_Z_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")


# Centered log ratio (clr) elemental datafile from original cps data and normalised element selections -------

# load compositions package
library(compositions)

# Figure 14 - cps clr for clr_elements - based on mean > 50 cps - same as Fig3A
theme_set(theme_bw(8))
ARD_Fig14 <- ARD_COMP_filter1_clr_long  %>%
  filter(param %in% plot_elements2_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements2_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Filtered elements cps clr based on mean>50 cps stats")
ARD_Fig14


# Change element list to m0.1 (cps_sum mean >0.1%) list or ARD_elements4 filtered list
# remove ratios to leave only elements and scatter - can change elements used in clr here!
ARD_clr_elements4 <- purrr::discard(plot_elements4_ARD,.p = ~stringr::str_detect(.x,"inc_coh"))
ARD_clr_elements4 <- purrr::discard(ARD_clr_elements4 ,.p = ~stringr::str_detect(.x,"coh_inc"))
ARD_clr_elements4

# apply clr to cps data and convert all 0 to NA for plotting clarity
ARD_COMP_filter1_clr4 <- select(ARD_COMP_filter1, all_of(ARD_clr_elements4)) %>% 
  select(all_of(ARD_clr_elements4)) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
ARD_COMP_filter1_clr4 <- as_tibble(ARD_COMP_filter1_clr4) %>% 
  na_if(0)
ARD_COMP_filter1_clr4

ARD_COMP_filter1_clr4 <- bind_cols(ARD_COMP_filter1_text, ARD_COMP_filter1_clr4) %>% 
  relocate(inc_coh, .after = Mo_coh) %>% 
  relocate(coh_inc, .after = inc_coh)
ARD_COMP_filter1_clr4

write.csv(ARD_COMP_filter1_clr4,"Output/ARD/14A_cps_filter1_clr_elements4.csv", row.names = FALSE)

# convert clr datafile to long format for plotting
ARD_COMP_filter1_clr4_long <- select(ARD_COMP_filter1_clr4,  Core, depth_cm, SH20_age, kcps, MSE, all_of(plot_elements4_ARD)) %>%
  pivot_longer(all_of(plot_elements4_ARD), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
head(ARD_COMP_filter1_clr4_long)
tail(ARD_COMP_filter1_clr4_long)

# Figure 14a - cps clr for clr_elements
theme_set(theme_bw(8))
ARD_Fig14A <- ARD_COMP_filter1_clr4_long  %>%
  filter(param %in% plot_elements4_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements4_ARD)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  geom_hline(yintercept = c(13, 35, 64.5, 155, 175, 190, 284, 310, 320, 326), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Filtered elements cps clr based on mean>0.1% %cps_sum elements")
ARD_Fig14A

ggarrange(ARD_Fig14, ARD_Fig14A, nrow = 2)
ggsave("Figures/ARD/Fig 14_cps_clr_elements2&4.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Plots vs age ------------------------------------------------------------

# Figure 15 %cps sum vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
ARD_Fig15 <- ARD_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements4_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements4_ARD)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "%cps sum") +
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%cps sum (filtered elements): CONISS XRF")
ggsave("Figures/ARD/Fig 15_cps_sum_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Figure 16 cps/Ti as Z-scores vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
ARD_Fig16 <- ARD_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements4_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements4_ARD)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps/Ti Z-score") +
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), colour = "red", lty = 2, alpha = 0.7)+
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")
ggsave("Figures/ARD/Fig 16_Ti_norm_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Figure 17 cps clr (clr_elements4) vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
ARD_Fig17 <- ARD_COMP_filter1_clr4_long  %>%
  filter(param %in% plot_elements4_ARD) %>%
  mutate(param = fct_relevel(param, plot_elements4_ARD)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps clr") +
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), colour = "red", lty = 2, alpha = 0.7)+
  ggtitle("Centered log ratio (clr) cps (filtered elements): CONISS XRF")
ggsave("Figures/ARD/Fig 17_clr_cps_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Comparison plots 

#plot Figure 15 and 16 side by side
ggarrange(ARD_Fig15, ARD_Fig16, nrow = 1)
ggsave("Figures/ARD/Fig 15&16_cps_sum_&_Ti_Z_CONISS_2017_age.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

#plot Figure 15 and 17 side by side
ggarrange(ARD_Fig15, ARD_Fig17, nrow = 1)
ggsave("Figures/ARD/Fig 15&17_clr_&_Ti_Z_CONISS_2017_age.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")



# PART 2 - Yanou Lake - YAN-ITRAX-SH20  ------------------------------------------------------------------

# YAN-ITRAX-SH20 - Import cps data ---------------------------------------------
library(dplyr)

# Location of YAN-ITRAX-SH20 composite dataset 
# https://github.com/stever60/Bertrand_2022/tree/main/Data

# Composite YAN - ITRAX cps data & Sh20 BACON age-depth model (added in Excel)
YAN_x.raw <- read_csv("Data/YAN_ITRAX_COMP_SH20.csv")
YAN_x.raw

# Filter data to keep only valid data and remove kcps <-2sd and MSE >2sd ----------------------------
YAN_kcps.mean <- mean(YAN_x.raw$kcps)
YAN_kcps.sd <- 2*sd(YAN_x.raw$kcps)
YAN_kcps.thres <- YAN_kcps.mean - YAN_kcps.sd 
YAN_kcps.thres

YAN_MSE.mean <- mean(YAN_x.raw$MSE)
YAN_MSE.sd <- 2*sd(YAN_x.raw$MSE)
YAN_MSE.thres <- YAN_MSE.mean + YAN_MSE.sd 
YAN_MSE.thres

# filter by validity and then remove analyses above MSE threshold and below kcps - 2s & 
YAN_COMP <- filter(YAN_x.raw, validity == "1") %>% 
  filter(MSE < YAN_MSE.thres) %>%
  filter(kcps > YAN_kcps.thres)

# Remove unwanted columns - might need to change these - and write to file
YAN_COMP_filter <- select(YAN_COMP, 
                          -c(filename, position_mm, 
                             core_depth_cm, min, max,
                             median, mean_2s_range:Foffset,
                             D1, S1, S2, S3)) %>% 
  rename(SH20_age = mean)

write.csv(YAN_COMP_filter,"Output/YAN/1_cps_filter.csv", row.names = FALSE)

# Define elements lists -------------------------------------------------

library(PeriodicTable)

#list of all possible elements
elements <- c(symb(1:117), "Mo_inc", "Mo_coh", "TS_sum", "cps_sum")

YAN_elements <- select(YAN_COMP_filter, c(Mg:Mo_coh)) %>% 
  names()
YAN_elements

# define elements associated with detector gas and splutter
machine_elements <- select(YAN_COMP_filter, c(Ar, Ta, W)) %>% 
  names()

# define REE - not used in YAN
#REE <- select(YAN_COMP_filter, c(La:Ho)) %>% 
#  names()
#REE

# Remove machine generated elements Ar (tube gas), Ta & W (splutter)
YAN_elements1 <- select(YAN_COMP_filter, -c(Core:MSE, all_of(machine_elements), Ir, Pt))  %>% 
  names()
YAN_elements1

# Calculate normalising factors TS (Total Scatter), cps_sum & inc/coh -------------------------------------------------------

# Add to columns - find a way to replace [9:67] with column headings as "Mg":"Mo_coh"
YAN_COMP_filter

# Type in row numbers of elements and scatter to calculate  across 
YAN_rowsums <- 9:49

YAN_COMP_filter1 <- YAN_COMP_filter %>% 
  replace(is.na(.), 0) %>%
  mutate(TS_sum = Mo_inc + Mo_coh) %>% 
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  mutate(coh_inc = Mo_coh / Mo_inc) %>%
  mutate(cps_sum = rowSums(.[YAN_rowsums]))
#mutate(cps_sum = row_sums(YAN_elements))

YAN_COMP_filter1

# Standardise and centre cps data  - using column names --------------------------------------

# list of element column names for plotting
YAN_col <- select(YAN_COMP_filter1, c(Mg:inc_coh)) %>% names()
YAN_col

# list of all element and other parameter column names 
YAN_col1 <- select(YAN_COMP_filter1, c(kcps:cps_sum)) %>% names()
YAN_col1

YAN_COMP_filter1.Z <- YAN_COMP_filter1
YAN_COMP_filter1.Z[, YAN_col] <- scale(YAN_COMP_filter1[, YAN_col], center = TRUE, scale = TRUE)
YAN_COMP_filter1.Z

# Convert to long format for plotting -------------------------------------------------

YAN_COMP_filter1_long <- select(YAN_COMP_filter1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
YAN_COMP_filter1_long 

YAN_COMP_filter1_long.Z <- select(YAN_COMP_filter1.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
YAN_COMP_filter1_long.Z 

# Generate cps summary stats & write to file  --------------------------------------------
library(psych)

# Create an extended version with a bunch of stats 
YAN_summary <- YAN_COMP_filter1 %>%
  select(all_of(YAN_col1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter cps dataset and stats to file --------------------------------------------------
write.csv(YAN_COMP_filter1,"Output/YAN/2_cps_filter1.csv", row.names = FALSE)
write.csv(YAN_summary,"Output/YAN/2.1_cps_filter1_stats.csv", row.names = FALSE)
write.csv(YAN_COMP_filter1.Z,"Output/YAN/3_cps_filter1_Z.csv", row.names = FALSE)

# Filter cps element column list based on mean and max cps values -------------------

YAN_mean50 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 50)
}
YAN_mean200 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 200)
}
YAN_max50 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 50)
}

YAN_m50 <- select(YAN_COMP_filter1, 
                  Core, depth_cm, SH20_age, kcps, MSE, 
                  Mg:Mo_coh & where(YAN_mean50)) %>% print()
YAN_m200 <- select(YAN_COMP_filter1, 
                   Core, depth_cm, SH20_age, kcps, MSE, 
                   Mg:Mo_coh & where(YAN_mean200)) %>% print()
YAN_mx50 <- select(YAN_COMP_filter1, 
                   Core, depth_cm, SH20_age, kcps, MSE, 
                   Mg:Mo_coh & where(YAN_max50), TS_sum:cps_sum) %>% print()

# Create lists of column names to take forward
YAN_col_m50 <- select(YAN_COMP_filter1, Mg:Mo_coh 
                      & -machine_elements 
                      & where(YAN_mean50)) %>% names() %>% print()
YAN_col_m200 <- select(YAN_COMP_filter1, Mg:Mo_coh 
                       &  -machine_elements 
                       & where(YAN_mean200)) %>% names() %>% print()
YAN_col_m50_mx50 <- select(YAN_COMP_filter1, Mg:Mo_coh 
                           & -machine_elements 
                           & where(YAN_max50) 
                           & where(YAN_mean50)) %>% names() %>% print()

# Correlation matrices -------------------------------------------------------

library(GGally)
library(dplyr)
# Correlation plot with max200 elements - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=2))

ggcorr(YAN_COMP_filter1[, YAN_col1], method = c("everything", "pearson"), 
       size = 2, label = FALSE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/YAN/Fig 0_Corr_matrix_cps_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Manual list of elements with significant correlations - made visually from correlation matrix
plot_elements1_YAN <- c("Si", "P", "S", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Cu", "Zn", 
                        "Br", "Rb", "Sr", "Zr", "Y", "Ba", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")

# Stats generated list of elements to take forward based on m50
plot_elements2_YAN <- c(YAN_col_m50, "inc_coh", "coh_inc")

ggcorr(YAN_COMP_filter1[,  plot_elements2_YAN], method = c("everything", "pearson"), 
       size = 4, label = TRUE, label_size = 3, label_alpha = TRUE, label_round=2) 
ggsave("Figures/YAN/Fig 0_Corr_matrix_cps_m50.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
theme_set(theme_bw(base_size=8))
ggpairs(YAN_COMP_filter1, columns = plot_elements2_YAN, upper = list(continuous = wrap("cor", size = 2)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Figures/YAN/Fig 0__Corr-den_matrix_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary cps plots vs depth --------------------------------------------------------

library(tidypaleo)

# Figure 1 - cps elements1
theme_set(theme_bw(8))
YAN_Fig1 <- YAN_COMP_filter1_long  %>%
  filter(param %in% plot_elements1_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements1_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Signifcantly correlated elements cps summary")
YAN_Fig1
ggsave("Figures/YAN/Fig 1_cps_elements1.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 2 - cps elements2
theme_set(theme_bw(8))
YAN_Fig2 <- YAN_COMP_filter1_long  %>%
  filter(param %in% plot_elements2_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements2_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm), group = ) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps based on mean>50 cps stats")
YAN_Fig2
ggsave("Figures/YAN/Fig 2_cps_elements2.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 3 - cps elements2 Z-scores
theme_set(theme_bw(8))
YAN_Fig3 <- YAN_COMP_filter1_long.Z  %>%
  filter(param %in% plot_elements2_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements2_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps Z-scores based on mean>50 cps stats")
YAN_Fig3
ggsave("Figures/YAN/Fig 3_cps_elements2_Z.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Create a cps centered log ratio (clr) elemental datafile and plot  -------

# load compositions package
library(compositions)
YAN_COMP_filter1_text <- select(YAN_COMP_filter1, Core:MSE, inc_coh, coh_inc)
YAN_COMP_filter1_text

# look at element list to use
plot_elements2_YAN

# remove scatter ratios to leave only subset of elements and scatter parameters i.e., closed sum measurements 
# can change elements used in clr here!
YAN_clr_elements2 <- purrr::discard(plot_elements2_YAN,.p = ~stringr::str_detect(.x,"inc_coh"))
YAN_clr_elements2 <- purrr::discard(YAN_clr_elements2 ,.p = ~stringr::str_detect(.x,"coh_inc"))
YAN_clr_elements2

# apply clr to cps data and convert all 0 to NA for plotting clarity
YAN_COMP_filter1_clr1 <- select(YAN_COMP_filter1, all_of(YAN_clr_elements2)) %>% 
  select(all_of(YAN_clr_elements2)) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
YAN_COMP_filter1_clr1 <- as_tibble(YAN_COMP_filter1_clr1) %>% 
  na_if(0)
YAN_COMP_filter1_clr1

YAN_COMP_filter1_clr <- bind_cols(YAN_COMP_filter1_text, YAN_COMP_filter1_clr1) %>% 
  relocate(inc_coh, .after = Mo_coh) %>% 
  relocate(coh_inc, .after = inc_coh)
YAN_COMP_filter1_clr

write.csv(YAN_COMP_filter1_clr,"Output/YAN/3.1_cps_filter1_clr.csv", row.names = FALSE)

# convert clr datafile to long format for plotting
YAN_COMP_filter1_clr_long <- select(YAN_COMP_filter1_clr,  Core, depth_cm, SH20_age, 
                                    kcps, MSE, inc_coh, coh_inc, all_of(plot_elements2_YAN)) %>%
  pivot_longer(all_of(plot_elements2_YAN), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
head(YAN_COMP_filter1_clr_long)
tail(YAN_COMP_filter1_clr_long)

# Figure 3a - cps clr for clr_elements
theme_set(theme_bw(8))
YAN_Fig3A <- YAN_COMP_filter1_clr_long  %>%
  filter(param %in% plot_elements2_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements2_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps clr based on mean>50 cps stats")
YAN_Fig3A
ggsave("Figures/YAN/Fig 3A_cps_clr_elements.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Normalization (elements and scatter parameters) ---------------------------------------------------------

# Inc normalised
YAN_inc_norm <- YAN_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Mo_inc))
# coh/inc normalised - Boyle 2015
YAN_coh_inc_norm <- YAN_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ (Mo_coh/Mo_inc)))
# Total scatter (TS) normalised - Kylander et al 2011/12
YAN_TS_norm <- YAN_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ TS_sum))
# cps_sum normalised
YAN_cps_sum_norm <- YAN_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ cps_sum))
# Ti normalised
YAN_Ti_norm <- YAN_COMP_filter1 %>% 
  mutate(across(c("Mg":"Mo_coh"), .fns = ~./ Ti))
# Natural log Ti normalized replacing -Inf with NA --------
YAN_Ln_Ti_norm <- YAN_Ti_norm  %>% mutate(across(c("Mg":"Mo_coh"),
                                                 log)) 
is.na(YAN_Ln_Ti_norm)<-sapply(YAN_Ln_Ti_norm, is.infinite)
# Replace NA with 0 - if needed
# YAN_Ln_Ti_norm[is.na(YAN_Ln_Ti_norm)]<-0

# check = 1
YAN_Ti_norm$Ti
# check = 0
YAN_Ln_Ti_norm$Ti

# Standardize and center Ti-normalized data - can rpelace Ti normlaised with others if needed --------------------------------------
YAN_Ln_Ti_norm.Z <- YAN_Ln_Ti_norm
YAN_Ln_Ti_norm.Z[, YAN_col1] <- scale(YAN_Ln_Ti_norm[, YAN_col1], center = TRUE, scale = TRUE)
YAN_Ln_Ti_norm.Z

# write ARD_Ln_Ti_norm.Z to file for Part 3 & 4
write.csv(YAN_Ln_Ti_norm.Z,"Output/YAN/YAN_Ln_Ti_norm.Z.csv", row.names = FALSE)

# Elements and scatter as %cps sum - Bertrand et al 2021 -------- same result as %TSN
YAN_cps_sum_norm_pc <- YAN_cps_sum_norm %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(cps_pc_sum = rowSums(.[YAN_rowsums]))

#Check everything = 100
head(YAN_cps_sum_norm_pc$cps_pc_sum)

YAN_col3 <- select(YAN_cps_sum_norm_pc, c(Mg:cps_pc_sum)) %>% 
  names()
YAN_col3

# Create %cps_sum summary stats table
YAN_cps_sum_pc_summary <- YAN_cps_sum_norm_pc %>%
  select(all_of(YAN_col3)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Elements and scatter as percentage of TSN sum - Roberts et al 2017-------- same result as %cps sum
YAN_TSN_sum <- YAN_TS_norm %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_sum = rowSums(.[YAN_rowsums]))
YAN_TSN_sum

YAN_TSN_pc <- YAN_TSN_sum %>% mutate(across(c("Mg":"Mo_coh"),
                                            .fns = ~./TSN_sum)) %>% 
  mutate(across(c("Mg":"Mo_coh"),.fns = ~.*100)) %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum = rowSums(.[YAN_rowsums]))
YAN_TSN_pc

#Check everything = 100
head(YAN_TSN_pc$TSN_pc_sum)

YAN_col4 <- select(YAN_TSN_pc, c(Mg:TSN_pc_sum)) %>% 
  names()
YAN_col4

# Create TSN stats table
YAN_TSN_pc_summary <- YAN_TSN_pc %>%
  select(all_of(YAN_col4)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Filter elements based on %TSN max > 0.5 and mean >0.1 - *** replace %TSN with %cps sum *** ----------------
YAN_max0.5 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 0.5)
}
YAN_mean0.1 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 0.1)
}

YAN_mx0.5 <- select(YAN_TSN_pc, Mg:Mo_coh & where(YAN_max0.5)) %>% print()
YAN_m0.1 <- select(YAN_TSN_pc, Mg:Mo_coh & where(YAN_mean0.1)) %>% print()

# Create names lists
YAN_TSN_mx0.5 <- select(YAN_TSN_pc, Mg:Mo_coh 
                        & -machine_elements 
                        & where(YAN_max0.5)) %>% names()

YAN_TSN_m0.1 <- select(YAN_TSN_pc, Mg:Mo_coh 
                       & -machine_elements 
                       & where(YAN_mean0.1)) %>% names()
YAN_TSN_mx0.5
YAN_TSN_m0.1

# New element list based on  YAN_col_m0.1 + Si, P, S 
YAN_TSN_m0.1_list <- c("Si", "P", "S", YAN_TSN_m0.1)
YAN_TSN_m0.1_list

# Add filtered element TSN sum to TSN_pc dataset
YAN_TSN_pc1 <- YAN_TSN_pc %>% 
  replace(is.na(.), 0) %>%
  mutate(TSN_pc_sum1 = rowSums(.[YAN_TSN_m0.1_list]))
YAN_TSN_pc1
# check new total <100 but >95%
head(YAN_TSN_pc1$TSN_pc_sum1)

# Create TSN stats table with filtered element TSN_Sum added to end 
YAN_col5 <- select(YAN_TSN_pc1, Mg:TSN_pc_sum1) %>% names()

YAN_TSN_pc1_summary <- YAN_TSN_pc1 %>%
  select(all_of(YAN_col5)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write normalised datasets to file --------------------------------------------------
write.csv(YAN_Ti_norm,"Output/YAN/4_Ti_norm.csv", row.names = FALSE)
write.csv(YAN_inc_norm,"Output/YAN/5_inc_norm.csv", row.names = FALSE)
write.csv(YAN_coh_inc_norm,"Output/YAN/5.1_coh_inc_norm.csv", row.names = FALSE)
write.csv(YAN_TS_norm,"Output/YAN/6_TS_norm.csv", row.names = FALSE)
write.csv(YAN_cps_sum_norm,"Output/YAN/7_cps_sum_norm.csv", row.names = FALSE)
write.csv(YAN_TSN_pc,"Output/YAN/8_TSN_pc.csv", row.names = FALSE)
write.csv(YAN_TSN_pc_summary,"Output/YAN/8.1_TSN_pc_stats.csv", row.names = FALSE)
write.csv(YAN_TSN_pc1,"Output/YAN/8.2_TSN_pc1.csv", row.names = FALSE)
write.csv(YAN_TSN_pc1_summary,"Output/YAN/8.3_TSN_pc1_stats.csv", row.names = FALSE)
write.csv(YAN_cps_sum_norm_pc,"Output/YAN/9_cps_sum_norm_pc.csv", row.names = FALSE)
write.csv(YAN_cps_sum_pc_summary,"Output/YAN/9.1_cps_sum_pc_stats.csv", row.names = FALSE)
write.csv(YAN_Ln_Ti_norm,"Output/YAN/10_Ln_Ti_norm.csv", row.names = FALSE)
write.csv(YAN_Ln_Ti_norm.Z,"Output/YAN/11_Ln_Ti_norm_Z.csv", row.names = FALSE)

# Convert to long format ------------------------------------------------
YAN_inc_norm_long <- select(YAN_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
YAN_inc_norm_long

YAN_coh_inc_norm_long <- select(YAN_coh_inc_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
YAN_inc_norm_long

YAN_Ln_Ti_norm_long <- select(YAN_Ln_Ti_norm,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
YAN_Ln_Ti_norm_long

YAN_Ln_Ti_norm.Z_long <- select(YAN_Ln_Ti_norm.Z,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col1)) %>%
  pivot_longer(all_of(YAN_col1), names_to = "param", values_to = "value")
YAN_Ln_Ti_norm.Z_long

YAN_cps_sum_norm_pc_long <- select(YAN_cps_sum_norm_pc,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col3)) %>%
  pivot_longer(all_of(YAN_col3), names_to = "param", values_to = "value")
YAN_cps_sum_norm_pc_long

YAN_TSN_pc1_long <- select(YAN_TSN_pc1,  Core, depth_cm, SH20_age, kcps, MSE, all_of(YAN_col5)) %>%
  pivot_longer(all_of(YAN_col5), names_to = "param", values_to = "value")
YAN_TSN_pc1_long

# Summary normalised plots--------------------------------------

# # Figure 4 cps/inc. -----------------------------------------------------

# manually defined
plot_elements3_YAN <- c("K", "Ca", "Ti", "Mn", "Fe", "Cu", "Zn", 
                        "Br", "Sr", "Zr", "Mo_inc", "Mo_coh", "inc_coh", "coh_inc")
# OR 

# stats defined 
plot_elements3_YAN <- c(YAN_TSN_m0.1_list, "inc_coh", "coh_inc")

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
# Roberts et al 2017 - CONISS with broken stick Hellingers dist defined zones - red dashed lines
theme_set(theme_bw(8))
YAN_Fig4a <- YAN_inc_norm_long  %>%
  filter(param %in% plot_elements3_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements3_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc.", y = "Depth (cm)") +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Incoherent (inc.) scatter normalised (filtered elements): CONISS XRF")

theme_set(theme_bw(8))
YAN_Fig4b <- YAN_coh_inc_norm_long  %>%
  filter(param %in% plot_elements3_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements3_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc.", y = "Depth (cm)") +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Incoherent (inc.) scatter normalised (filtered elements): CONISS XRF")
ggarrange(YAN_Fig4a, YAN_Fig4b, nrow = 2)
ggsave("Figures/YAN/Fig 4_cps_inc.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 5 %TSN----------------------------------------------------------------

plot_elements4_YAN <- c("K", "Ca", "Ti", "Mn", "Fe", "Cu", "Zn", 
                        "Br", "Sr", "Zr","Mo_inc", "Mo_coh", "inc_coh", "coh_inc")
# OR

# stats defined 
plot_elements4_YAN <- c(YAN_TSN_m0.1_list, "inc_coh", "coh_inc")
#plot_elements4_YAN <- c(YAN_m50, "inc_coh", "coh_inc")

# Plot multiple variables & add CONISS zone boundaries defined from subsample data 
# Roberts et al 2017 - CONISS with broken stick Hellingers dist defined zones - red dashed lines
theme_set(theme_bw(8))
YAN_Fig5 <- YAN_TSN_pc1_long  %>%
  filter(param %in% plot_elements4_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements4_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%TSN", y = "Depth (cm)")
YAN_Fig5 +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%TSN (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 5_TSN_pc.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Add CONISS Zones using in built tidypalaeo CONISS - THIS TAKES TOO LONG - NEED A BETTER COMPUTER!
#coniss1_YAN <- YAN_TSN_pc1_long %>%
#  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
#  nested_chclust_coniss()

# Figure 6 - %TSN with CONISS
#YAN_Fig6 <- YAN_Fig5 +
#  layer_dendrogram(coniss1_YAN, aes(y = depth_cm), param = "CONISS") +
#  layer_zone_boundaries(coniss1_YAN, aes(y = depth_cm, col = "blue", lty = 2, alpha = 0.7)) +
#  ggtitle("%TSN (filtered elements): CONISS ITRAX")
#YAN_Fig6
#ggsave("Figures/YAN/Fig 6_TSN_pc_CONISS.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

#Figure 7 - TSN with CONISS and Roberts et al (2017) zone boundaries
#YAN_Fig6 +
#  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%TSN (filtered elements): Coniss zone comparison")
#ggsave("Figures/YAN/Fig 7_TSN_pc_CONISS_comparison.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 8  %cps_sum
#define colour scheme for 6 groups (5 guano and 1 non-guano)
guano_zone_colours <- c("#BDBDBD", "#00441B", "#00441B", "#006D2C", "#238B45", "#74C476")
theme_set(theme_bw(8))
YAN_Fig8 <- YAN_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements3_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements3_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm, aes(colour = guano_zone_colours))) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%cps_sum", y = "Depth (cm)")
YAN_Fig8 +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%cps sum (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 8_cps_sum_pc.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Add CONISS Zones using in built tidypalaeo CONISS 
#coniss2_YAN <- YAN_cps_sum_norm_pc_long %>%
#  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
#  nested_chclust_coniss()

# Figure 9 %cps_sum with CONISS
#YAN_Fig9 <- YAN_Fig8 +
#  layer_dendrogram(coniss2_YAN, aes(y = depth_cm), param = "CONISS") +
#  layer_zone_boundaries(coniss2_YAN, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
#  ggtitle("%cps sum (filtered elements): CONISS ITRAX")
#YAN_Fig9
#ggsave("Figures/YAN/Fig 9_cps_sum_pc_CONISS.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 10 - compare to Roberts et al. 2017 CONISS zoning 
#YAN_Fig9 +
#  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("%cps sum (filtered elements): CONISS comparison")
#ggsave("Figures/YAN/Fig 10_cps_sum_pc_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# create a new elemenet list without Ti
plot_elements5_YAN <- purrr::discard(plot_elements4_YAN,.p = ~stringr::str_detect(.x,"Ti"))
plot_elements5_YAN

# Figure 11 Ti normalised
theme_set(theme_bw(8))
YAN_Fig11 <- YAN_Ln_Ti_norm_long  %>%
  filter(param %in% plot_elements5_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements5_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_lineh(size = 0.5) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Ln(cps/Ti)", y = "Depth (cm)") +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 11_Ti_norm.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 12 - Ti normalised as Z-scores 
theme_set(theme_bw(8))
YAN_Fig12 <- YAN_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements5_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements5_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_color_manual(values = guano_zone_colours) +
  scale_fill_manual(values = guano_zone_colours) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Ln(cps/Ti) Z-scores", y = "Depth (cm)") +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 12_Ti_norm_Z.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

ggarrange(YAN_Fig11, YAN_Fig12, nrow=2)
ggsave("Figures/YAN/Fig 11_12_Ti_norm & Ti_norm_Z.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Add CONISS Zones using in built tidypalaeo CONISS 
#coniss3_YAN <- YAN_Ln_Ti_norm.Z_long %>%
#  nested_data(qualifiers = c(SH20_age, depth_cm), key = param, value = value, trans = scale) %>%
#  nested_chclust_coniss()

# Figure 13
#YAN_Fig12 + labs(x = "Ln(cps/Ti Z-scores)", y = "Depth (cm)") +
#  layer_dendrogram(coniss3_YAN, aes(y = depth_cm), param = "CONISS ITRAX") +
#  layer_zone_boundaries(coniss3_YAN, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
#  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS ITRAX")
#ggsave("Figures/YAN/Fig 13_Ti_Z_CONISS.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")

# Figure 14
#YAN_Fig12 + labs(x = "cps/Ti Z-scores", y = "Depth (cm)") +
#  layer_dendrogram(coniss3_YAN, aes(y = depth_cm), param = "CONISS") +
#  layer_zone_boundaries(coniss3_YAN, aes(y = depth_cm, col = "red", lty = 2, alpha = 0.7)) +
#  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
#  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS comparison")
#ggsave("Figures/YAN/Fig 14_Ti_Z_CONISS_comp.pdf",
#       height = c(15), width = c(30), dpi = 600, units = "cm")


# Centered log ratio (clr) elemental datafile from original cps data and normalised element selections -------

# load compositions package
library(compositions)

# Figure 14 - cps clr for clr_elements - based on mean > 50 cps - same as Fig3A
theme_set(theme_bw(8))
YAN_Fig14 <- YAN_COMP_filter1_clr_long  %>%
  filter(param %in% plot_elements2_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements2_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps clr based on mean>50 cps stats")
YAN_Fig14

# Change element list to m0.1 (cps_sum mean >0.1%) list or YAN_elements4 filtered list
# remove ratios to leave only elements and scatter - can change elements used in clr here!
YAN_clr_elements4 <- purrr::discard(plot_elements4_YAN,.p = ~stringr::str_detect(.x,"inc_coh"))
YAN_clr_elements4 <- purrr::discard(YAN_clr_elements4 ,.p = ~stringr::str_detect(.x,"coh_inc"))
YAN_clr_elements4

# apply clr to cps data and convert all 0 to NA for plotting clarity
YAN_COMP_filter1_clr4 <- select(YAN_COMP_filter1, all_of(YAN_clr_elements4)) %>% 
  select(all_of(YAN_clr_elements4)) %>%
  filter(!if_any(everything(), is.na)) %>% 
  clr()
YAN_COMP_filter1_clr4 <- as_tibble(YAN_COMP_filter1_clr4) %>% 
  na_if(0)
YAN_COMP_filter1_clr4

YAN_COMP_filter1_clr4 <- bind_cols(YAN_COMP_filter1_text, YAN_COMP_filter1_clr4) %>% 
  relocate(inc_coh, .after = Mo_coh) %>% 
  relocate(coh_inc, .after = inc_coh)
YAN_COMP_filter1_clr4

write.csv(YAN_COMP_filter1_clr4,"Output/YAN/14A_cps_filter1_clr_elements4.csv", row.names = FALSE)

# convert clr datafile to long format for plotting
YAN_COMP_filter1_clr4_long <- select(YAN_COMP_filter1_clr4,  Core, depth_cm, SH20_age, kcps, MSE, all_of(plot_elements4_YAN)) %>%
  pivot_longer(all_of(plot_elements4_YAN), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
head(YAN_COMP_filter1_clr4_long)
tail(YAN_COMP_filter1_clr4_long)

# Figure 14a - cps clr for clr_elements
theme_set(theme_bw(8))
YAN_Fig14A <- YAN_COMP_filter1_clr4_long  %>%
  filter(param %in% plot_elements4_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements4_YAN)) %>%
  ggplot(aes(x = value, y = depth_cm)) +
  geom_point(size = 0.01) +
  geom_lineh(size = 0.5) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  geom_hline(yintercept = c(19, 33, 192, 214, 249, 261, 279), colour = "red", lty = 2, alpha = 0.7) +
  labs(x = "cps", y = "Depth (cm)") +
  ggtitle("Filtered elements cps clr based on mean>0.1% %cps_sum elements")
YAN_Fig14A

ggarrange(YAN_Fig14, YAN_Fig14A, nrow = 2)
ggsave("Figures/YAN/Fig 14_cps_clr_elements2&4.pdf",
       height = c(15), width = c(30), dpi = 600, units = "cm")



# Plots vs age ------------------------------------------------------------

# Figure 15 %cps sum vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
YAN_Fig15 <- YAN_cps_sum_norm_pc_long  %>%
  filter(param %in% plot_elements4_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements4_YAN)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "%cps sum") +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("%cps sum (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 15_cps_sum_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Figure 16 cps/Ti as Z-scores vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
YAN_Fig16 <- YAN_Ln_Ti_norm.Z_long  %>%
  filter(param %in% plot_elements4_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements4_YAN)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps/Ti Z-score") +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Log Ti-normalised Z-scores (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 16_Ti_norm_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Figure 17 cps clr (clr_elements4) vs SH20_age  ----------------------------------
theme_set(theme_bw(7))
YAN_Fig17 <- YAN_COMP_filter1_clr4_long  %>%
  filter(param %in% plot_elements4_YAN) %>%
  mutate(param = fct_relevel(param, plot_elements4_YAN)) %>%
  ggplot(aes(x = SH20_age, y = value)) +
  geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps clr") +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), colour = "red", lty = 2, alpha = 0.7) +
  ggtitle("Centered log ratio (clr) cps (filtered elements): CONISS XRF")
ggsave("Figures/YAN/Fig 17_clr_cps_CONISS_2017.pdf",
       height = c(30), width = c(15), dpi = 600, units = "cm")

# Comparison plots 

#plot Figure 15 and 16 side by side
ggarrange(YAN_Fig15, YAN_Fig16, nrow = 1)
ggsave("Figures/YAN/Fig 15&16_cps_sum_&_Ti_Z_CONISS_2017_age.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

#plot Figure 15 and 17 side by side
ggarrange(YAN_Fig15, YAN_Fig17, nrow = 1)
ggsave("Figures/YAN/Fig 15&17_clr_&_Ti_Z_CONISS_2017_age.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")



# PART 3 - Calibration - Using Pre-matched ARD & YAN ITRAX and XRF 1 cm datasets in Roberts et al. (2017)  ------------------------------------------------------------------

# Currently uses 1cm ITRAX smoothed dataset matched to XRF & ICPMS ratio data from original composite ARD & YAN excel files
# Correlation and calibration/regression code below utlimately uses log-ratio Z-scores

# Still to do: 
# 1) match dataset from original cps and subsample data using itrax.r
# 2) clr on cps (ITRAX) and XRF, ICPMS %, ppm subsample data
# 3) elemental correlation, glm - 


# Set up and Import data  ----------------------------------------------------------

# Data are located here:
# https://github.com/stever60/Bertrand_2022/tree/main/Data

# Dataframe db contains all oxide & elemental subsample and XRF-CS matched data for P, Ca, Cu, Zn, 
# for ARD and YAN (terrestrial and marine)

# data as Ti ratios
db <- read_csv("Data/ARD_YAN_ratio.csv") %>% 
  filter(!if_any(everything(), is.na)) #looking for row without anything in them - wont find any in this dataset
db

# Set up lists of variables ---------------------------------------------

# ratio data
XRF_oxide <- c("P2O5_TiO2", "CaO_TiO2", "Cu_TiO2", "Zn_TiO2", "Sr_TiO2")

XRF_elemental <- c("P_Ti", "Ca_Ti", "Cu_Ti", "Zn_Ti", "Sr_Ti")

XRF_CS <- c("P_Ti_CS","Ca_Ti_CS", "Cu_Ti_CS", "Zn_Ti_CS", "Sr_Ti_CS")

all_oxide <- c("P2O5_TiO2", "CaO_TiO2", "Cu_TiO2", "Zn_TiO2", "Sr_TiO2",
               "P_Ti_CS","Ca_Ti_CS", "Cu_Ti_CS", "Zn_Ti_CS", "Sr_Ti_CS")

all_elemental <- c("P_Ti", "Ca_Ti", "Cu_Ti", "Zn_Ti", "Sr_Ti",
                   "P_Ti_CS","Ca_Ti_CS", "Cu_Ti_CS", "Zn_Ti_CS", "Sr_Ti_CS")

all_variables <- c("P2O5_TiO2", "CaO_TiO2", "Cu_TiO2", "Zn_TiO2", "Sr_TiO2",
              "P_Ti", "Ca_Ti", "Cu_Ti", "Zn_Ti", "Sr_Ti",
              "P_Ti_CS","Ca_Ti_CS", "Cu_Ti_CS", "Zn_Ti_CS", "Sr_Ti_CS")

# create log ratio, centered and centred, standardised (Z-scores) log ratios for all data --------------------------
# to do - clr version of below from raw cps (ITRAX) and , ppm (XRF, ICPMS) elemenatal data
# clr would need to be done before depth matching in itrax.R first

db_ln <- db
db_ln[, all_variables]<- log(db_ln[all_variables]) %>% 
  as_tibble()
db_ln

db_cln <- db_ln
db_cln[, all_variables] <- scale(db_cln[all_variables], center = TRUE) %>% 
  as_tibble()
db_cln

db_cln.Z <- db_ln
db_cln.Z[, all_variables] <- scale(db_cln.Z[all_variables], center = TRUE, scale = TRUE) %>% 
  as_tibble()
db_cln.Z

# Note1: centred oxide and elemental data should be the same when centered
# Note2: centered and standardised (scaled) produces the same output as centering alone for log ratio data

#  Write data without blanks to Output folder --------------------------------------------------
write.csv(db,"Output/Calibration/ARD_YAN_ratio.csv", row.names = FALSE)
write.csv(db_ln,"Output/Calibration/ARD_YAN_ln.csv", row.names = FALSE)
write.csv(db_cln,"Output/Calibration/ARD_YAN_cln.csv", row.names = FALSE)
write.csv(db_cln.Z,"Output/Calibration/ARD_YAN_cln_Z.csv", row.names = FALSE)

# Generate stats & write to file  --------------------------------------------
library(psych)

dbln_summary <- db_ln %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

dbln_guano <- filter(db_ln, Guano_Sum == 'Guano')
dbln_summary_guano <- dbln_guano %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

dbln_non_guano <- filter(db_ln, Guano_Sum == 'Non Guano')
dbln_summary_non_guano <- dbln_non_guano %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter cps dataset and stats to file --------------------------------------------------
write.csv(dbln_summary,"Output/Calibration/ARD-YAN_stats_all_ln.csv", row.names = FALSE)
write.csv(dbln_summary_guano,"Output/Calibration/ARD-YAN_all_stats_guano_ln.csv", row.names = FALSE)
write.csv(dbln_summary_non_guano,"Output/Calibration/ARD-YAN_all_stats_non_guano_ln.csv", row.names = FALSE)

# Reorder the data to plot_order column and convert to long format --------------------------------

# YAN and ARD non guano plot first then ARD GZ1-5 on top
db1 <- arrange(db_cln.Z, Plot_order)
head(db1)
tail(db1)

# Convert to long format - i.e., one value per row per variable - for later
db1_long <- select(db1,  Record, strat_depth, SH20_age, Group, Guano_Zone, Guano_Sum, Plot_order, all_of(all_variables)) %>%
  pivot_longer(all_of(all_variables), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
db1_long


# Correlation and covariance matrices -------------------------------------

# Examine correlation matrix for PCA df
cor <- cor(db1)
round(cor, 2)
# p-values
library("Hmisc")
cor1 <- rcorr(as.matrix(db1))
cor1

# Examine co-variance in the PCA df
cov <- cov(db1)
round(cov,2)
# p-values
cov1 <- rcov(as.matrix(db1))
cov1


# Correlation summary for all data - terrestrial and marine ---------------
library(GGally)

# Correlation plot - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=8))
ggcorr(db1[, all_variables], method = c("everything", "pearson"), 
       size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/Calibration/Fig 1A_Corr_matrix_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
theme_set(theme_bw(base_size=8))
ggpairs(db1, columns = all_variables, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Figures/Calibration/Fig 1B_Corr-den_matrix_all_data.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Correlation density matrix - by Guano Zone or Unit
theme_set(theme_bw(base_size=8))
ggpairs(db1, columns = all_variables, upper = list(continuous = wrap("cor", size = 1)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.5)),
        ggplot2::aes(colour = Guano_Zone, title="Correlation plot: Guano Zone", alpha = 0.5))
ggsave("Figures/Calibration/Fig 1C_Corr-den-unit_matrix_all_data.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# get full stats for correlation tests if needed
res_P <- cor.test(db1$P_Ti_CS, db1$P_Ti, method = "pearson")
print(res_P)
res_Ca <- cor.test(db1$Ca_Ti_CS, db1$Ca_Ti, method = "pearson")
print(res_Ca)
res_Cu <- cor.test(db1$Cu_Ti_CS, db1$Cu_Ti, method = "pearson")
print(res_Cu)
res_Zn <- cor.test(db1$Zn_Ti_CS, db1$Zn_Ti, method = "pearson")
print(res_Zn)
res_Sr <- cor.test(db1$Sr_Ti_CS, db1$Sr_Ti, method = "pearson")
print(res_Sr)

# Save console outputs above to file
sink("Output/Calibration/ARD-YAN_all_corr_stats.txt")
print(res_P)
print(res_Ca)
print(res_Cu)
print(res_Zn)
print(res_Sr)
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
res_P

#axis labels 
P_CS <-  "P/Ti [XRF-CS]"
P_WD <- "P/Ti [WD-XRF]"
Ca_CS <-  "Ca/Ti [XRF-CS]"
Ca_WD <- "Ca/Ti [WD-XRF]"
Cu_CS <-  "Cu/Ti [XRF-CS]"
Cu_WD <- "Cu/Ti [WD-XRF]"
Zn_CS <-  "Zn/Ti [XRF-CS]"
Zn_WD <- "Zn/Ti [WD-XRF]"
Sr_CS <-  "Sr/Ti [XRF-CS]"
Sr_WD <- "Sr/Ti [WD-XRF]"


#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
# original zone colors
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# guano zones green
guano_colours1 <- c("#BDBDBD", "#74C476", "#74C476", "#74C476", "#74C476", "#74C476", "#969696")
# guano zones graduational green
guano_colours2 <- c("#BDBDBD", "#00441B", "#006D2C", "#238B45", "#74C476", "#74C476", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

library(ggpubr)
theme_set(theme_bw(base_size=10) + theme(
  plot.title = element_text(color="black", size=10, face="bold.italic"),
  axis.title.x = element_text(color="black", size=10),
  axis.title.y = element_text(color="black", size=10)
))
db1_corr_P <- ggscatter(db1, x = "P_Ti_CS", y = "P_Ti", shape ="Guano_Sum", colour = "Guano_Sum", 
                        palette = c("#74C476", "#BDBDBD"), size = 2,
                        add = "reg.line", 
                        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        conf.int = TRUE,
                        cor.coef = TRUE,
                        cor.coeff.args = list(method = "pearson", 
                                              label.x.npc = "left", label.y.npc = "top", size = 3),
                        ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                        xlab = P_CS, ylab = P_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

db1_corr_Ca <- ggscatter(db1, x = "Ca_Ti_CS", y = "Ca_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                        palette = c("#74C476", "#BDBDBD"), size = 2,
                        add = "reg.line", 
                        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        conf.int = TRUE,
                        cor.coef = TRUE, 
                        cor.coeff.args = list(method = "pearson", 
                                              label.x.npc = "left", label.y.npc = "top", size = 3),
                        ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                        xlab = Ca_CS, ylab = Ca_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

db1_corr_Cu <- ggscatter(db1, x = "Cu_Ti_CS", y = "Cu_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Cu_CS, ylab = Cu_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

db1_corr_Zn <- ggscatter(db1, x = "Zn_Ti_CS", y = "Zn_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Zn_CS, ylab = Zn_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

db1_corr_Sr <- ggscatter(db1, x = "Sr_Ti_CS", y = "Sr_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Sr_CS, ylab = Sr_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))
# plot as 2x3 grid
corrplot <- ggarrange(db1_corr_P, db1_corr_Ca, db1_corr_Cu, db1_corr_Zn, db1_corr_Sr,
          ncol = 2, nrow = 3, labels = c("a", "b", "c", "d", "e"))
corrplot
annotate_figure(corrplot, top = text_grob("Correlation: ARD & YAN Freshwater & Marine [log-n Z-scores; 95% ellispses]", 
                                           color = "black", face = "bold", size = 12))
ggsave("Figures/Calibration/Fig 1D_Correlation_plots_all_data.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")


# REMOVE YANOU MARINE DATA - i.e, only use freshwater data only from this point onwards for analysis and plotting   -----------------
library(sjmisc)
df <- select(db_ln, Record:Plot_order) %>% 
  filter(!(Record == 'YAN' & strat_depth>270))
df
# is.na(df)<-sapply(df, is.infinite)
# df contains only terrestrial matched data - for regression analysis
# replace(is.na(.), 0) %>% #convert NA to 0 for correlation plotting - not with log data
# mutate_if(is.logical,as.numeric) #to convert overcomes parsing problem
# To convert NA to 0 in base R: EPMA_df[is.na(EPMA_df)] <- 0

# Convert to long format 
df_long <- select(df,  Record, strat_depth, SH20_age, Group, Guano_Zone, Plot_order, all_of(all_variables)) %>%
  pivot_longer(all_of(all_variables), names_to = "param", values_to = "value")
#relocate(param, .before = Type)
df_long

# Generate stats for terrestrial data & write to file  --------------------------------------------
library(psych)

df_summary <- df %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

df_guano <- filter(df, Guano_Sum == 'Guano')
df_summary_guano <- df_guano %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

df_non_guano <- filter(df, Guano_Sum == 'Non Guano')
df_summary_non_guano <- df_non_guano %>%
  select(all_of(all_variables)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#  Write filter cps dataset and stats to file --------------------------------------------------
write.csv(df_summary,"Output/Calibration/ARD-YAN_Terr_stats_all_ln.csv", row.names = FALSE)
write.csv(df_summary_guano,"Output/Calibration/ARD-YAN_Terr_stats_guano_ln.csv", row.names = FALSE)
write.csv(df_summary_non_guano,"Output/Calibration/ARD-YAN_Terr_stats_non_guano_ln.csv", row.names = FALSE)

# Transform data ----------------------------------------------------------
# Plot original ln data
plot(df[, all_variables], pch=19, cex = 0.05)

# Standardise and centre variables (cln dataset is already centred)
# Z-scores using scale() function - values are mean of 0 and +/-1 of 1 std dev
# copy the original filename to a new filename, center and apply a z-score transform, look at it,plot it
df.cln <- df
df.cln[, all_variables] <- scale(df.cln[all_variables], center = TRUE)
df.cln

df.Z <- df
df.Z[, all_variables] <- scale(df.Z[all_variables], center = TRUE, scale = TRUE)
df.Z

# Replot transformed (standardised and centered) data
plot(df.Z[, all_variables], pch=19, cex = 0.05)

#  Write transformed data to file --------------------------------------------------
write.csv(df.cln,"Output/Calibration/ARD-YAN_Terr_cln.csv", row.names = FALSE)
write.csv(df.Z,"Output/Calibration/ARD-YAN_Terr_Z-scores.csv", row.names = FALSE)

# Correlation matrices -------------------------------------------------------
library(GGally)

# Correlation plot - use this to see where positive/significant correlations as an overview
theme_set(theme_bw(base_size=8))
ggcorr(df.Z[, all_variables], method = c("everything", "pearson"), size = 4, label = TRUE, label_alpha = TRUE, label_round=2) 
ggsave("Figures/Calibration/Fig 1E_Corr_matrix_terrestrial.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plot
theme_set(theme_bw(base_size=8))
ggpairs(df.Z, columns = all_variables, upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot")
ggsave("Figures/Calibration/Fig 1F_Corr-den_matrix.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Correlation density matrix - by Guano Zone or Unit
theme_set(theme_bw(base_size=8))
ggpairs(df.Z, columns = all_variables, upper = list(continuous = wrap("cor", size = 1)),
                       lower = list(continuous = wrap("points", alpha = 0.5, size=0.5)),
                       ggplot2::aes(colour = Guano_Zone, title="Correlation plot: Guano Zone", alpha = 0.5))
ggsave("Figures/Calibration/Fig 1G_Corr-den-unit_matrix.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# get full stats for correlation tests if needed
res_P_terr <- cor.test(df.Z$P_Ti_CS, df.Z$P_Ti, method = "pearson")
print(res_P_terr)
res_Ca_terr <- cor.test(df.Z$Ca_Ti_CS, df.Z$Ca_Ti, method = "pearson")
print(res_Ca_terr)
res_Cu_terr <- cor.test(df.Z$Cu_Ti_CS, df.Z$Cu_Ti, method = "pearson")
print(res_Cu_terr)
res_Zn_terr <- cor.test(df.Z$Zn_Ti_CS, df.Z$Zn_Ti, method = "pearson")
print(res_Zn_terr)
res_Sr_terr <- cor.test(df.Z$Sr_Ti_CS, df.Z$Sr_Ti, method = "pearson")
print(res_Sr_terr)

# Save console outputs above to file
sink("Output/Calibration/ARD-YAN_terr_corr_stats.txt")
print(res_P_terr)
print(res_Ca_terr)
print(res_Cu_terr)
print(res_Zn_terr)
print(res_Sr_terr)
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
res_P

library(ggpubr)
df.Z_corr_P <- ggscatter(df.Z, x = "P_Ti_CS", y = "P_Ti", shape ="Guano_Sum", colour = "Guano_Sum", 
                        palette = c("#74C476", "#BDBDBD"), size = 2,
                        add = "reg.line", 
                        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        conf.int = TRUE,
                        cor.coef = TRUE,
                        cor.coeff.args = list(method = "pearson", 
                                              label.x.npc = "left", label.y.npc = "top", size = 3),
                        ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                        xlab = P_CS, ylab = P_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

df.Z_corr_Ca <- ggscatter(df.Z, x = "Ca_Ti_CS", y = "Ca_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Ca_CS, ylab = Ca_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

df.Z_corr_Cu <- ggscatter(df.Z, x = "Cu_Ti_CS", y = "Cu_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Cu_CS, ylab = Cu_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

df.Z_corr_Zn <- ggscatter(df.Z, x = "Zn_Ti_CS", y = "Zn_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Zn_CS, ylab = Zn_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

df.Z_corr_Sr <- ggscatter(df.Z, x = "Sr_Ti_CS", y = "Sr_Ti", colour = "Guano_Sum", shape ="Guano_Sum",
                         palette = c("#74C476", "#BDBDBD"), size = 2,
                         add = "reg.line", 
                         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE,
                         cor.coef = TRUE, 
                         cor.coeff.args = list(method = "pearson", 
                                               label.x.npc = "left", label.y.npc = "top", size = 3),
                         ellipse = TRUE, ellipse.level = .95, mean.point = FALSE,
                         xlab = Sr_CS, ylab = Sr_WD) + 
  theme(legend.title = element_blank(),legend.text = element_text(size = 10, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL))

# plot as 2x3 grid
corrplot1 <- ggarrange(df.Z_corr_P, df.Z_corr_Ca, df.Z_corr_Cu, df.Z_corr_Zn, df.Z_corr_Sr,
          ncol = 2, nrow = 3, labels = c("a", "b", "c", "d", "e"))
corrplot1 
annotate_figure(corrplot1, top = text_grob("Correlation & 95% ellipses: ARD & YAN Freshwater [log-n Z-scores; 95% ellispses]", 
                                      color = "black", face = "bold", size = 12))
ggsave("Figures/Calibration/Fig 1H_Correlation_plots_terr.pdf", 
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Further normality and stats tests for each ratio are at the start of PART 5: Regression Models


# PART 4 ITRAX & XRF comparison vs Age plots ----------------------------------------------------------------

library(ggpubr)
library(dplyr)

# Import ARD & YAN ITRAX data (ARD_Ln_Ti_norm.Z & YAN_Ln_Ti_norm.Z) from Parts 1 and 2 without running them all again
ARD_Ln_Ti_norm.Z <- read.csv("Output/ARD/ARD_Ln_Ti_norm.Z.csv") %>% 
  as_tibble()
ARD_Ln_Ti_norm.Z

YAN_Ln_Ti_norm.Z <- read.csv("Output/YAN/YAN_Ln_Ti_norm.Z.csv") %>% 
  as_tibble()
YAN_Ln_Ti_norm.Z

# ARD
# filter to remove YAN data from ITRAX-XRF matched df - df is log transformed and reordered dataset
df_ARD <- filter(df, Record=='ARD')

# Standardise and centre variables -  Z-scores using scale() function - values are mean of 0 and +/-1 of 1 std dev
# copy the original filename to a new filename, apply a z-score transform, look at it,plot it
df_ARD.Z <- df_ARD
df_ARD.Z [, all_variables] <- scale(df_ARD.Z[all_variables], center = TRUE, scale = TRUE)
df_ARD.Z 

tail(df_ARD.Z)

#  YAN

# filter to remove marine data >7.5ka or >9 ka
YAN_Ln_Ti_norm.Z_nonmarine  <- YAN_Ln_Ti_norm.Z %>% filter(SH20_age<9000)
tail(YAN_Ln_Ti_norm.Z_nonmarine)

# filter to remove ARD data & YAN data >7.4 ka - dbln is log transfromed and includes marine data for YAN
# df is fwater data only for YAN
db_YAN <- filter(db_ln, Record=='YAN' & SH20_age<7400)
all_variables <- all_variables
db_YAN

# Standardise and centre variables -  Z-scores using scale() function - values are mean of 0 and +/-1 of 1 std dev
# copy the original filename to a new filename, apply a z-score transform, look at it,plot it
db_YAN.Z <- db_YAN
db_YAN.Z [, all_variables] <- scale(db_YAN.Z[all_variables], center = TRUE, scale = TRUE)
db_YAN.Z 

Greys1
Greens1

# VARIABLE 1 - Ca - uXRF and XRF data versus age  ---------------------------------

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
ARD_age <- ggplot(data = ARD_Ln_Ti_norm.Z, aes(SH20_age, Ca)) + 
  geom_line(colour = "#74C476", alpha = 1, size = 0.5) +
  geom_point(data = df_ARD.Z, aes(SH20_age, Ca_Ti), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = df_ARD.Z, aes(SH20_age, Ca_Ti), colour = "#006D2C", size = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Ca", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
ARD_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, Ca)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, Ca_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, Ca_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN): Ca", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p1.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p1.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp1 <- list(p1.1, p1.2)
plot_grid(plotlist = pp1, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_Ca.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")

# VARIABLE 2 - P - uXRF and XRF data versus age  ---------------------------------

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
ARD_age <- ggplot(data = ARD_Ln_Ti_norm.Z, aes(SH20_age, P)) + 
  geom_line(colour = "#74C476", alpha = 1, size = 0.5) +
  geom_point(data = df_ARD.Z, aes(SH20_age, P_Ti), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = df_ARD.Z, aes(SH20_age, P_Ti), colour = "#006D2C", size = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): P", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
ARD_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, P)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, P_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, P_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-6,4,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN: P", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p2.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p2.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp2 <- list(p2.1, p2.2)
plot_grid(plotlist = pp2, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_P.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")

# VARIABLE 3 - Cu - uXRF and XRF data versus age  ---------------------------------

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
ARD_age <- ggplot(data = ARD_Ln_Ti_norm.Z, aes(SH20_age, Cu)) + 
  geom_line(colour = "#74C476", alpha = 1, size = 0.5) +
  geom_point(data = df_ARD.Z, aes(SH20_age, Cu_Ti), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = df_ARD.Z, aes(SH20_age, Cu_Ti), colour = "#006D2C", size = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Cu", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
ARD_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, P)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, Cu_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, Cu_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-6,4,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN: Cu", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p3.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p3.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp3 <- list(p3.1, p3.2)
plot_grid(plotlist = pp3, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_Cu.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")

# VARIABLE 4 - Zn - uXRF and XRF data versus age  ---------------------------------

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
ARD_age <- ggplot(data = ARD_Ln_Ti_norm.Z, aes(SH20_age, Zn)) + 
  geom_line(colour = "#74C476", alpha = 1, size = 0.5) +
  geom_point(data = df_ARD.Z, aes(SH20_age, Zn_Ti), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = df_ARD.Z, aes(SH20_age, Zn_Ti), colour = "#006D2C", size = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Zn", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
ARD_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, Zn)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, Zn_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, Zn_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-6,4,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN: Zn", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p4.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p4.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp4 <- list(p4.1, p4.2)
plot_grid(plotlist = pp4, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_Zn.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")


# Model 5 - Sr - uXRF and XRF data versus age - TO DO  ---------------------------------

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
ARD_age <- ggplot(data = ARD_Ln_Ti_norm.Z, aes(SH20_age, Zn)) + 
  geom_line(colour = "#74C476", alpha = 1, size = 0.5) +
  geom_point(data = df_ARD.Z, aes(SH20_age, Sr_Ti), colour = "#006D2C", fill = "blue", size = 0.5) + 
  geom_line(data = df_ARD.Z, aes(SH20_age, Sr_Ti), colour = "#006D2C", size = 0.5) +
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-4,6,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Ardley Lake (ARD): Zn", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black"))+
  geom_vline(xintercept = c(1257, 2552, 2933, 3800, 4163, 4418, 5298, 5874, 6538, 6936), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
ARD_age

YAN_age <- ggplot(data = YAN_Ln_Ti_norm.Z_nonmarine, aes(SH20_age, Zn)) + 
  geom_line(colour = "#969696", alpha = 1, size = 0.5) +
  geom_point(data = db_YAN.Z, aes(SH20_age, Sr_Ti), colour = "#252525", fill = "#252525", size = 0.5) + 
  geom_line(data = db_YAN.Z, aes(SH20_age, Sr_Ti), colour = "#252525", size = 0.5) + 
  expand_limits(x = -100, y = -3) +
  scale_x_continuous(breaks=seq(0,9000,1000), minor_breaks = seq(NULL), expand = c(0.05,0)) +
  scale_y_continuous(breaks=seq(-6,4,2), minor_breaks = seq(NULL), expand = c(0,0)) +
  labs(title = "Yanou Lake (YAN: Zn", x = "Age (cal a BP)", y = "Z-score") + 
  theme(axis.ticks.length=unit(0.15, "cm"), axis.text = element_text(colour = "black")) +
  geom_vline(xintercept = c(3369.89, 4836.78, 5520.18, 6058.28, 6876.51, 7024.81, 7550.8), 
             size = 0.5, colour = "red", lty = "dotted", alpha = 0.5)
YAN_age

# align axes exactly for both graphs
p4.1 <- ARD_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
p4.2 <- YAN_age + coord_cartesian(xlim = c(0,9000), ylim = c(-4,6))
pp4 <- list(p4.1, p4.2)
plot_grid(plotlist = pp4, ncol = 1, nrow = 2, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_Sr.pdf",
       height = c(6), width = c(10), dpi = 600, units = "cm")

# Final age GRID - Models 1-4 -------------------------------------------
ggarrange(p1.1,  p2.1, p3.1, p4.1, p1.2, p2.2, p3.2, p4.2, ncol = 2, nrow = 4, align = "v")
ggsave("Figures/Fig 1C_ITRAX-XRF_Age_final_CONISS_ALL.pdf",
       height = c(12), width = c(20), dpi = 600, units = "cm")


# PART 5: Linear regression using ggplot ----------------------------------------------

# MODEL 1 -----------------------------------------------------------------

# set up x and y variables and axis, title labels 
x1.reg <- df.Z$Ca_Ti_CS
y1.reg <- df.Z$Ca_Ti

x1 <- "Ca_Ti_CS"
y1 <- "Ca_Ti"

xtitle1 <- bquote('Ln'~'('*Ca/Ti*') XRF-CS [Z-score]')
ytitle1 <- bquote('Ln'~'('*Ca/Ti*') WD-XRF [Z-score]')

xlab1 <- xlab(bquote('Ln'~'('*Ca/Ti*') XRF-CS [Z-score]'))
xlab1_y <- ylab(bquote('Ln'~'('*Ca/Ti*') XRF-CS [Z-score]'))
xlab1_y1 <- ylab(bquote('Ln'~'('*Ca/Ti*') XRF-CS'))

ylab1 <- ylab(bquote('Ln'~'('*Ca/Ti*') WD-XRF [Z-score]'))
ylab1_y1 <- ylab(bquote('Ln'~'('*Ca/Ti*') WD-XRF'))
ylab1_x <- xlab(bquote('Ln'~'('*Ca/Ti*') WD-XRF [Z-score]'))
ylab1_x1 <- xlab(bquote('Ln'~'('*Ca/Ti*') WD-XRF'))

# Build linear regression model 1
model_1 <- lm(y1.reg~x1.reg, data = df.Z)
# y~x for y = mx + c linear reg model
# uXRF data used to predict subsample data values
# y and x can both be Ln Z scores (log of centred ratio data) or y (subsample) can be % or ppm Z scores and x (uXRF) Ln Z scores 

# Produce sumary statistics for model 1 in console
model_1
summary(model_1)
confint(model_1)

# Save console outputs above to file
sink("Output/Calibration/model1_summary.txt")
print(summary(model_1))
sink("Output/Calibration/model1_confint.txt")
print(confint(model_1))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
model_1
summary(model_1)
confint(model_1)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x1.reg = c(1, 2, 3))
predict(model_1, newdata = new.y)
predict(model_1, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int1 <- predict(model_1, interval = "prediction")
data_1 <- cbind(df.Z, pred.int1)
data_1
as_tibble(data_1)

# Plotting set up ---------------------------------------------------------

#clear plot window
dev.off()

#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
# original zone colors
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# guano zones green
guano_colours1 <- c("#BDBDBD", "#74C476", "#74C476", "#74C476", "#74C476", "#74C476", "#969696")
# guano zones graduational green
guano_colours2 <- c("#BDBDBD", "#00441B", "#006D2C", "#238B45", "#74C476", "#74C476", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

# Density and qq plots and residuals plots to assess normality --------------------------------------------

library(ggpubr)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12),
  axis.ticks.length=unit(0.15, "cm")))

den.x1 <- ggdensity(data_1, x = x1) +
  xlab1 +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

den.y1 <- ggdensity(data_1, x = y1) +
  ylab1_x +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

qq.x1 <- ggqqplot(data_1, x = x1,
                  colour = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(xtitle1) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))

qq.y1 <- ggqqplot(data_1, x = y1,
                  colour = "Guano_Zone",
                  fill = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(ytitle1) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.ticks.length=unit(0.15, "cm"))

# Residuals vs fitted values plot
modf <- fortify(model_1)
res <- ggplot(modf, aes(x = .fitted, y = .resid)) + geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle('Residuals plot') +
  theme(axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm"))

# Violin + box plots for x and y - to add as F below in 3x3 grid --------

# draw violin plot for log Ca/Ti ITRAX data - use df for Ln(Ca/Ti) data-1 to plot Z-score data
# Guano_sum defines the groups - only 3 groups - Guano, non-guano and Yanou
# set Guano_sum column to factor rather than numerical plot
df$Guano_Sum <- as.factor(df$Guano_Sum)

v1 <- ggplot(df, aes(x=Guano_Sum, y=Ca_Ti_CS, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v1

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter1 <- v1 + geom_jitter(shape=16, position=position_jitter(0.05))

# Final violin + box plot plot
vbp1 <- vjitter1 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab1_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ylim (-1, 4) +
  ggtitle("Violin Plot")
vbp1

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter <- v1 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp1_nolegend <- vjitter1 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab1_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ylim (-1, 4) +
  ggtitle("Violin Plot")
vbp1_nolegend

# draw violin plot for XRF subsample data
v1.2 <- ggplot(df, aes(x=Guano_Sum, y=Ca_Ti, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v1.2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter1.2 <- v1.2 + geom_jitter(shape=16, position=position_jitter(0.05))

# final violin + box plot plot
vbp1.2 <- vjitter1.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab1_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ylim (-1, 4) +
  ggtitle("Violin Plot") 
vbp1.2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter1 <- v1.2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp1.2_nolegend <- vjitter1 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab1_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ylim (-1, 4) +
  ggtitle("Violin Plot") 
vbp1.2_nolegend

# plot as 3 x3 grid - one example of violin plot 
ggarrange(den.x1, qq.x1, den.y1, qq.y1, res, vbp1_nolegend, 
          ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("Figures/Calibration/Model1_Ca/Fig 1_Regression_summary1.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# plot as 2 x2 grid - qq and violin plots
ggarrange( qq.x1, qq.y1, vbp1_nolegend, vbp2_nolegend, 
           ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("Figures/Calibration/Model1_Ca/Fig 2_Regression_summary2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Shapiro-Wilk test ----------

# Use this to test normality of one variable (univariate) at a time 
shapiro.test(x1.reg)
shapiro.test(y1.reg)
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed ...
# Note that 'perfect' normality is hard to achieve in most downcore geochemical datasets!

# Save console outputs above to file
sink("Output/Calibration/model1_shapiro_x.txt")
print(shapiro.test(x1.reg))
sink("Output/Calibration/model1_shapiro_y.txt")
print(shapiro.test(y1.reg))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
shapiro.test(x1.reg)
shapiro.test(y1.reg)

# Durbin-Watson Test  ---------------------------------------------------

# Use this to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
library(lmtest)
dwtest(model_1)
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
# Save console outputs above to file
sink("Output/Calibration/model1_dw_test.txt")
print(dwtest(model_1))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
dwtest(model_1)

# Plot Regression & confidence intervals, prediction intervals ---------------------------

# Regression Plot Model 1 ----------------------------------------------------------------
library(ggpubr)
library(dplyr)

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
p1 <- ggplot(data_1,aes(x1.reg, y1.reg)) + 
  geom_point(data = data_1, aes(x1.reg, y1.reg, colour = Guano_Zone, fill = Guano_Zone, shape = Guano_Zone), size = 2, alpha = 1) +
  scale_shape_manual(values = plot_shapes) +
  scale_fill_manual(values = guano_colours1) +
  scale_color_manual(values = guano_colours1) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 7, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "right", #c(1,0 ),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=7, colour = "black"), axis.title=element_text(size=7, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2, byrow = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "#9ECAE1") +
  xlab1 +
  ylab1
p1

# Add prediction intervals
p1.1 <- p1 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle("Linear regression (solid blue), 95% CI (blue shaded) & prediction (dashed blue)") + 
  theme(axis.ticks.length=unit(0.15, "cm"))
p1.1

# Final MODEL 1 regression plot as 1x2 grid
ggarrange(p1.1,  ncol = 1, nrow = 1, labels = c("A"))
ggsave("Figures/Calibration/Model1_Ca/Fig 3_Regression_final.pdf",
       height = c(10), width = c(10), dpi = 600, units = "cm")


# MODEL 2 for P -----------------------------------------------------------------

# set up x and y variables and axis, title labels 
x2.reg <- df.Z$P_Ti_CS
y2.reg <- df.Z$P_Ti

x2 <- "P_Ti_CS"
y2 <- "P_Ti"

xtitle2 <- bquote('Ln'~'('*P/Ti*') XRF-CS [Z-score]')
ytitle2 <- bquote('Ln'~'('*P/Ti*') WD-XRF [Z-score]')

xlab2 <- xlab(bquote('Ln'~'('*P/Ti*') XRF-CS [Z-score]'))
xlab2_y <- ylab(bquote('Ln'~'('*P/Ti*') XRF-CS [Z-score]'))
xlab2_y1 <- ylab(bquote('Ln'~'('*P/Ti*') XRF-CS'))

ylab2 <- ylab(bquote('Ln'~'('*P/Ti*') WD-XRF [Z-score]'))
ylab2_y1 <- ylab(bquote('Ln'~'('*P/Ti*') WD-XRF'))
ylab2_x <- xlab(bquote('Ln'~'('*P/Ti*') WD-XRF [Z-score]'))
ylab2_x1 <- xlab(bquote('Ln'~'('*P/Ti**') WD-XRF'))

# Build linear regression model 1
model_2 <- lm(y2.reg~x2.reg, data = df.Z)
# y~x for y = mx + c linear reg model
# uXRF data used to predict subsample data values
# y and x can both be Ln Z scores (clr - centred log ratio - method) or y (subsample) can be % or ppm Z scores and x (uXRF) Ln Z scores 

# Produce sumary statistics for model 1 in console
model_2
summary(model_2)
confint(model_2)

# Save console outputs above to file
sink("Output/Calibration/model2_summary.txt")
print(summary(model_2))
sink("Output/Calibration/model2_confint.txt")
print(confint(model_2))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
model_2
summary(model_2)
confint(model_2)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x2.reg = c(1, 2, 3))
predict(model_2, newdata = new.y)
predict(model_2, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int2 <- predict(model_2, interval = "prediction")
data_2 <- cbind(df.Z, pred.int2)
data_2
as_tibble(data_2)

# Plotting set up ---------------------------------------------------------

#clear plot window
dev.off()

#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
# original zone colors
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# guano zones green
guano_colours1 <- c("#BDBDBD", "#74C476", "#74C476", "#74C476", "#74C476", "#74C476", "#969696")
# guano zones graduational green
guano_colours2 <- c("#BDBDBD", "#00441B", "#006D2C", "#238B45", "#74C476", "#74C476", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

# Density and qq plots and residuals plots to assess normality --------------------------------------------

library(ggpubr)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12),
  axis.ticks.length=unit(0.15, "cm")))

den.x2 <- ggdensity(data_2, x = x2) +
  xlab2 +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

den.y2 <- ggdensity(data_2, x = y2) +
  ylab2_x +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

qq.x2 <- ggqqplot(data_2, x = x2,
                  colour = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(xtitle2) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))

qq.y2 <- ggqqplot(data_2, x = y2,
                  colour = "Guano_Zone",
                  fill = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(ytitle2) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.ticks.length=unit(0.15, "cm"))

# Residuals vs fitted values plot
modf <- fortify(model_2)
res <- ggplot(modf, aes(x = .fitted, y = .resid)) + geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle('Residuals plot') +
  theme(axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm"))

# Violin + box plots for x and y - to add as F below in 3x3 grid --------

# set Guano_sum column to factor rather than numerical plot
df$Guano_Sum <- as.factor(df$Guano_Sum)

v2 <- ggplot(df, aes(x=Guano_Sum, y=P_Ti_CS, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v2 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter2 <- v2 + geom_jitter(shape=16, position=position_jitter(0.05))

# Final violin + box plot plot
vbp2 <- vjitter2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab2_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter2 <- v2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp2_nolegend <- vjitter2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab2_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp2_nolegend

# draw violin plot for XRF subsample data
v2.2 <- ggplot(df, aes(x=Guano_Sum, y=P_Ti, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v2.2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter2.2 <- v2.2 + geom_jitter(shape=16, position=position_jitter(0.05))

# final violin + box plot plot
vbp2.2 <- vjitter2.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab2_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp2.2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter2.2 <- v2.2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp2.2_nolegend <- vjitter2.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab2_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp2.2_nolegend

# plot as 3 x3 grid - one example of violin plot 
ggarrange(den.x2, qq.x2, den.y2, qq.y2, res, vbp2.2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("Figures/Calibration/Model2_P/Fig 1_Regression_summary1.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# plot as 2 x2 grid - qq and violin plots
ggarrange( qq.x2, qq.y2, vbp2_nolegend, vbp2.2_nolegend, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("Figures/Calibration/Model2_p/Fig 2_Regression_summary2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Shapiro-Wilk test ----------

# Use this to test normality of one variable (univariate) at a time 
shapiro.test(x2.reg)
shapiro.test(y2.reg)
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed ...
# Note that 'perfect' normality is hard to achieve in most downcore geochemical datasets!

# Save console outputs above to file
sink("Output/Calibration/model2_shapiro_x.txt")
print(shapiro.test(x2.reg))
sink("Output/Calibration/model2_shapiro_y.txt")
print(shapiro.test(y2.reg))
sink(file = NULL)
# reset back to console output
sink(file = NULL)
shapiro.test(x2.reg)
shapiro.test(y2.reg)

# Durbin-Watson Test  ---------------------------------------------------

# Use this to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
library(lmtest)
dwtest(model_2)
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
# Save console outputs above to file
sink("Output/Calibration/model2_dw_test.txt")
print(dwtest(model_2))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
dwtest(model_2)

# Plot Regression & confidence intervals, prediction intervals ---------------------------

# Regression Plot Model 2 ----------------------------------------------------------------
library(ggpubr)
library(dplyr)

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
p2 <- ggplot(data_2,aes(x2.reg, y2.reg)) + 
  geom_point(data = data_2, aes(x2.reg, y2.reg, colour = Guano_Zone, fill = Guano_Zone, shape = Guano_Zone), size = 2, alpha = 1) +
  scale_shape_manual(values = plot_shapes) +
  scale_fill_manual(values = guano_colours1) +
  scale_color_manual(values = guano_colours1) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 7, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "right", #c(1,0 ),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=7, colour = "black"), axis.title=element_text(size=7, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2, byrow = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "#9ECAE1") +
  xlab2 +
  ylab2
p2

# Add prediction intervals
p2.1 <- p2 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle("Linear regression (solid blue), 95% CI (blue shaded) & prediction (dashed blue)") + 
  theme(axis.ticks.length=unit(0.15, "cm"))
p2.1

# Final regression plot as 1x2 grid
ggarrange(p2.1,  ncol = 1, nrow = 1, labels = c("A"))
ggsave("Figures/Calibration/Model2_P/Fig 3_Regression_final.pdf",
       height = c(10), width = c(10), dpi = 600, units = "cm")


# MODEL 3 for Cu -----------------------------------------------------------------

# set up x and y variables and axis, title labels 
x3.reg <- df.Z$Cu_Ti_CS
y3.reg <- df.Z$Cu_Ti

x3 <- "Cu_Ti_CS"
y3 <- "Cu_Ti"

xtitle3 <- bquote('Ln'~'('*Cu/Ti*') XRF-CS [Z-score]')
ytitle3 <- bquote('Ln'~'('*Cu/Ti*') WD-XRF [Z-score]')

xlab3 <- xlab(bquote('Ln'~'('*Cu/Ti*') XRF-CS [Z-score]'))
xlab3_y <- ylab(bquote('Ln'~'('*Cu/Ti*') XRF-CS [Z-score]'))
xlab3_y1 <- ylab(bquote('Ln'~'('*Cu/Ti*') XRF-CS'))

ylab3 <- ylab(bquote('Ln'~'('*Cu/Ti*') WD-XRF [Z-score]'))
ylab3_y1 <- ylab(bquote('Ln'~'('*Cu/Ti*') WD-XRF'))
ylab3_x <- xlab(bquote('Ln'~'('*Cu/Ti*') WD-XRF [Z-score]'))
ylab3_x1 <- xlab(bquote('Ln'~'('*Cu/Ti**') WD-XRF'))

# Build linear regression model 1
model_3 <- lm(y3.reg~x3.reg, data = df.Z)
# y~x for y = mx + c linear reg model
# uXRF data used to predict subsample data values
# y and x can both be Ln Z scores (clr - centred log ratio - method) or y (subsample) can be % or ppm Z scores and x (uXRF) Ln Z scores 

# Produce sumary statistics for model 1 in console
model_3
summary(model_3)
confint(model_3)

# Save console outputs above to file
sink("Output/Calibration/model3_summary.txt")
print(summary(model_3))
sink("Output/Calibration/model3_confint.txt")
print(confint(model_3))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
model_3
summary(model_3)
confint(model_3)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x3.reg = c(1, 2, 3))
predict(model_3, newdata = new.y)
predict(model_3, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int3 <- predict(model_3, interval = "prediction")
data_3 <- cbind(df.Z, pred.int3)
data_3
as_tibble(data_3)

# Plotting set up ---------------------------------------------------------

#clear plot window
dev.off()

#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
# original zone colors
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# guano zones green
guano_colours1 <- c("#BDBDBD", "#74C476", "#74C476", "#74C476", "#74C476", "#74C476", "#969696")
# guano zones graduational green
guano_colours2 <- c("#BDBDBD", "#00441B", "#006D2C", "#238B45", "#74C476", "#74C476", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

# Density and qq plots and residuals plots to assess normality --------------------------------------------

library(ggpubr)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12),
  axis.ticks.length=unit(0.15, "cm")))

den.x3 <- ggdensity(data_3, x = x3) +
  xlab3 +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

den.y3 <- ggdensity(data_3, x = y3) +
  ylab3_x +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

qq.x3 <- ggqqplot(data_3, x = x3,
                  colour = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(xtitle3) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))

qq.y3 <- ggqqplot(data_3, x = y3,
                  colour = "Guano_Zone",
                  fill = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(ytitle3) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.ticks.length=unit(0.15, "cm"))

# Residuals vs fitted values plot
modf <- fortify(model_3)
res <- ggplot(modf, aes(x = .fitted, y = .resid)) + geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle('Residuals plot') +
  theme(axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm"))

# Violin + box plots for x and y - to add as F below in 3x3 grid --------

# set Guano_sum column to factor rather than numerical plot
df$Guano_Sum <- as.factor(df$Guano_Sum)

v3 <- ggplot(df, aes(x=Guano_Sum, y=Cu_Ti_CS, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v3

# violin plot with dot plot - using jitter to offset the data
# vdot <- v3 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter3 <- v3 + geom_jitter(shape=16, position=position_jitter(0.05))

# Final violin + box plot plot
vbp3 <- vjitter3 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab3_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp3

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter3 <- v3 + geom_jitter(shape=16, position=position_jitter(0.03))
vbp3_nolegend <- vjitter3 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab3_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp3_nolegend

# draw violin plot for XRF subsample data
v3.2 <- ggplot(df, aes(x=Guano_Sum, y=Cu_Ti, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v3.2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter3.2 <- v3.2 + geom_jitter(shape=16, position=position_jitter(0.05))

# final violin + box plot plot
vbp3.2 <- vjitter3.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab3_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp3.2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter3.2 <- v3.2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp3.2_nolegend <- vjitter3.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab3_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp3.2_nolegend

# plot as 3 x3 grid - one example of violin plot 
ggarrange(den.x3, qq.x3, den.y3, qq.y3, res, vbp3.2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("Figures/Calibration/Model3_Cu/Fig 1_Regression_summary1.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# plot as 2 x2 grid - qq and violin plots
ggarrange(qq.x3, qq.y3, vbp3_nolegend, vbp3.2_nolegend, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("Figures/Calibration/Model3_Cu/Fig 2_Regression_summary2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Shapiro-Wilk test ----------

# Use this to test normality of one variable (univariate) at a time 
shapiro.test(x3.reg)
shapiro.test(y3.reg)
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed ...
# Note that 'perfect' normality is hard to achieve in most downcore geochemical datasets!

# Save console outputs above to file
sink("Output/Calibration/model3_shapiro_x.txt")
print(shapiro.test(x3.reg))
sink("Output/Calibration/model3_shapiro_y.txt")
print(shapiro.test(y3.reg))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
shapiro.test(x3.reg)
shapiro.test(y3.reg)

# Durbin-Watson Test  ---------------------------------------------------

# Use this to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
library(lmtest)
dwtest(model_3)
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
# Save console outputs above to file
sink("Output/Calibration/model2_dw_test.txt")
print(dwtest(model_3))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
dwtest(model_3)

# Plot Regression & confidence intervals, prediction intervals ---------------------------

# Regression Plot Model 3 ----------------------------------------------------------------
library(ggpubr)
library(dplyr)

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
p3 <- ggplot(data_3,aes(x3.reg, y3.reg)) + 
  geom_point(data = data_3, aes(x3.reg, y3.reg, colour = Guano_Zone, fill = Guano_Zone, shape = Guano_Zone), size = 2, alpha = 1) +
  scale_shape_manual(values = plot_shapes) +
  scale_fill_manual(values = guano_colours1) +
  scale_color_manual(values = guano_colours1) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 7, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "right", #c(1,0 ),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=7, colour = "black"), axis.title=element_text(size=7, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2, byrow = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "#9ECAE1") +
  xlab3 +
  ylab3
p3

# Add prediction intervals
p3.1 <- p3 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle("Linear regression (solid blue), 95% CI (blue shaded) & prediction (dashed blue)") + 
  theme(axis.ticks.length=unit(0.15, "cm"))
p3.1

# Final regression plot as 1x2 grid
ggarrange(p3.1,  ncol = 1, nrow = 1, labels = c("A"))
ggsave("Figures/Calibration/Model3_Cu/Fig 3_Regression_final.pdf",
       height = c(10), width = c(10), dpi = 600, units = "cm")


# MODEL 4 for Zn -----------------------------------------------------------------

# set up x and y variables and axis, title labels 
x4.reg <- df.Z$Zn_Ti_CS
y4.reg <- df.Z$Zn_Ti

x4 <- "Zn_Ti_CS"
y4 <- "Zn_Ti"

xtitle4 <- bquote('Ln'~'('*Zn/Ti*') XRF-CS [Z-score]')
ytitle4 <- bquote('Ln'~'('*Zn/Ti*') WD-XRF [Z-score]')

xlab4 <- xlab(bquote('Ln'~'('*Zn/Ti*') XRF-CS [Z-score]'))
xlab4_y <- ylab(bquote('Ln'~'('*Zn/Ti*') XRF-CS [Z-score]'))
xlab4_y1 <- ylab(bquote('Ln'~'('*Zn/Ti*') XRF-CS'))

ylab4 <- ylab(bquote('Ln'~'('*Zn/Ti*') WD-XRF [Z-score]'))
ylab4_y1 <- ylab(bquote('Ln'~'('*Zn/Ti*') WD-XRF'))
ylab4_x <- xlab(bquote('Ln'~'('*Zn/Ti*') WD-XRF [Z-score]'))
ylab4_x1 <- xlab(bquote('Ln'~'('*Zn/Ti**') WD-XRF'))

# Build linear regression model 4
model_4 <- lm(y4.reg~x4.reg, data = df.Z)
# y~x for y = mx + c linear reg model
# uXRF data used to predict subsample data values
# y and x can both be Ln Z scores (clr - centred log ratio - method) or y (subsample) can be % or ppm Z scores and x (uXRF) Ln Z scores 

# Produce sumary statistics for model 1 in console
model_4
summary(model_4)
confint(model_4)

# Save console outputs above to file
sink("Output/Calibration/model4_summary.txt")
print(summary(model_4))
sink("Output/Calibration/model4_confint.txt")
print(confint(model_4))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
model_4
summary(model_4)
confint(model_4)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x4.reg = c(1, 2, 3))
predict(model_4, newdata = new.y)
predict(model_4, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int4 <- predict(model_4, interval = "prediction")
data_4 <- cbind(df.Z, pred.int4)
data_4
as_tibble(data_4)

# Plotting set up ---------------------------------------------------------

#clear plot window
dev.off()

#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

# Density and qq plots and residuals plots to assess normality --------------------------------------------

library(ggpubr)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12),
  axis.ticks.length=unit(0.15, "cm")))

den.x4 <- ggdensity(data_4, x = x4) +
  xlab4 +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

den.y4 <- ggdensity(data_4, x = y4) +
  ylab4_x +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

qq.x4 <- ggqqplot(data_4, x = x4,
                  colour = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(xtitle4) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))

qq.y4 <- ggqqplot(data_4, x = y4,
                  colour = "Guano_Zone",
                  fill = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(ytitle4) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.ticks.length=unit(0.15, "cm"))

# Residuals vs fitted values plot
modf <- fortify(model_4)
res <- ggplot(modf, aes(x = .fitted, y = .resid)) + geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle('Residuals plot') +
  theme(axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm"))

# Violin + box plots for x and y - to add as F below in 4x3 grid --------

# set Guano_sum column to factor rather than numerical plot
df$Guano_Sum <- as.factor(df$Guano_Sum)

v4 <- ggplot(df, aes(x=Guano_Sum, y=Zn_Ti_CS, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v4

# violin plot with dot plot - using jitter to offset the data
# vdot <- v4 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter4 <- v4 + geom_jitter(shape=16, position=position_jitter(0.05))

# Final violin + box plot plot
vbp4 <- vjitter4 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab4_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp4

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter4 <- v4 + geom_jitter(shape=16, position=position_jitter(0.04))
vbp4_nolegend <- vjitter4 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab4_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp4_nolegend

# draw violin plot for XRF subsample data
v4.2 <- ggplot(df, aes(x=Guano_Sum, y=Zn_Ti, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v4.2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter4.2 <- v4.2 + geom_jitter(shape=16, position=position_jitter(0.05))

# final violin + box plot plot
vbp4.2 <- vjitter4.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab4_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp4.2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter4.2 <- v4.2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp4.2_nolegend <- vjitter4.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab4_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp4.2_nolegend

# plot as 3 x3 grid - one example of violin plot 
ggarrange(den.x4, qq.x4, den.y4, qq.y4, res, vbp4.2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("Figures/Calibration/Model4_Zn/Fig 1_Regression_summary1.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# plot as 2 x2 grid - qq and violin plots
ggarrange(qq.x4, qq.y4, vbp4_nolegend, vbp4.2_nolegend, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("Figures/Calibration/Model4_Zn/Fig 2_Regression_summary2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Shapiro-Wilk test ----------

# Use this to test normality of one variable (univariate) at a time 
shapiro.test(x4.reg)
shapiro.test(y4.reg)
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed ...
# Note that 'perfect' normality is hard to achieve in most downcore geochemical datasets!

# Save console outputs above to file
sink("Output/Calibration/model4_shapiro_x.txt")
print(shapiro.test(x4.reg))
sink("Output/Calibration/model4_shapiro_y.txt")
print(shapiro.test(y4.reg))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
shapiro.test(x4.reg)
shapiro.test(y4.reg)

# Durbin-Watson Test  ---------------------------------------------------

# Use this to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
library(lmtest)
dwtest(model_4)
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
# Save console outputs above to file
sink("Output/Calibration/model2_dw_test.txt")
print(dwtest(model_4))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
dwtest(model_4)

# Plot Regression & confidence intervals, prediction intervals ---------------------------

# Regression Plot Model 4 ----------------------------------------------------------------
library(ggpubr)
library(dplyr)

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
p4 <- ggplot(data_4,aes(x4.reg, y4.reg)) + 
  geom_point(data = data_4, aes(x4.reg, y4.reg, colour = Guano_Zone, fill = Guano_Zone, shape = Guano_Zone), size = 2, alpha = 1) +
  scale_shape_manual(values = plot_shapes) +
  scale_fill_manual(values = guano_colours1) +
  scale_color_manual(values = guano_colours1) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 7, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "right", #c(1,0 ),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=7, colour = "black"), axis.title=element_text(size=7, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2, byrow = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "#9ECAE1") +
  xlab4 +
  ylab4
p4

# Add prediction intervals
p4.1 <- p4 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle("Linear regression (solid blue), 95% CI (blue shaded) & prediction (dashed blue)") + 
  theme(axis.ticks.length=unit(0.15, "cm"))
p4.1

# Final regression plot as 1x2 grid
ggarrange(p4.1,  ncol = 1, nrow = 1, labels = c("A"))
ggsave("Figures/Calibration/Model4_Zn/Fig 3_Regression_final.pdf",
       height = c(10), width = c(10), dpi = 600, units = "cm")


# MODEL 5 for Sr -----------------------------------------------------------------

# set up x and y variables and axis, title labels 
x5.reg <- df.Z$Sr_Ti_CS
y5.reg <- df.Z$Sr_Ti

x5 <- "Sr_Ti_CS"
y5 <- "Sr_Ti"

xtitle5 <- bquote('Ln'~'('*Sr/Ti*') XRF-CS [Z-score]')
ytitle5 <- bquote('Ln'~'('*Sr/Ti*') WD-XRF [Z-score]')

xlab5 <- xlab(bquote('Ln'~'('*Sr/Ti*') XRF-CS [Z-score]'))
xlab5_y <- ylab(bquote('Ln'~'('*Sr/Ti*') XRF-CS [Z-score]'))
xlab5_y1 <- ylab(bquote('Ln'~'('*Sr/Ti*') XRF-CS'))

ylab5 <- ylab(bquote('Ln'~'('*Sr/Ti*') WD-XRF [Z-score]'))
ylab5_y1 <- ylab(bquote('Ln'~'('*Sr/Ti*') WD-XRF'))
ylab5_x <- xlab(bquote('Ln'~'('*Sr/Ti*') WD-XRF [Z-score]'))
ylab5_x1 <- xlab(bquote('Ln'~'('*Sr/Ti**') WD-XRF'))

# Build linear regression model 1
model_5 <- lm(y5.reg~x5.reg, data = df.Z)
# y~x for y = mx + c linear reg model
# uXRF data used to predict subsample data values
# y and x can both be Ln Z scores (clr - centred log ratio - method) or y (subsample) can be % or ppm Z scores and x (uXRF) Ln Z scores 

# Produce sumary statistics for model 1 in console
model_5
summary(model_5)
confint(model_5)

# Save console outputs above to file
sink("Output/Calibration/model5_summary.txt")
print(summary(model_5))
sink("Output/Calibration/model5_confint.txt")
print(confint(model_5))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
model_5
summary(model_5)
confint(model_5)

# Create predicted values and upper/lower CI to check model is working
new.y <- data.frame(x5.reg = c(1, 2, 3))
predict(model_5, newdata = new.y)
predict(model_5, newdata = new.y, interval = "confidence")

# Add prediction intervals to model data frame
pred.int5 <- predict(model_5, interval = "prediction")
data_5 <- cbind(df.Z, pred.int5)
data_5
as_tibble(data_5)

# Plotting set up ---------------------------------------------------------

#clear plot window
dev.off()

#plot set up
#define colour scheme for 7 groups (5 guano and 2 non-guano) - ORIGINAL ORDER -
guano_colours <- c("#00441B", "#00441B", "#006D2C", "#238B45", "#74C476", "#BDBDBD", "#969696")
# define plot shapes for up to 12 groups 
plot_shapes <- c(21, 21, 21, 21, 21, 21, 25)
# 21:25 are outline = colour and fill = fill shapes

# Density and qq plots and residuals plots to assess normality --------------------------------------------

library(ggpubr)
theme_set(theme_bw(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold.italic"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12),
  axis.ticks.length=unit(0.15, "cm")))

den.x5 <- ggdensity(data_5, x = x5) +
  xlab5 +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

den.y5 <- ggdensity(data_5, x = y5) +
  ylab5_x +
  ylab(bquote('Density')) +
  ylim (0, 0.6)

qq.x5 <- ggqqplot(data_5, x = x5,
                  colour = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(xtitle5) +
  theme(plot.title = element_text(color="black", size=12, face="bold"))

qq.y5 <- ggqqplot(data_5, x = y5,
                  colour = "Guano_Zone",
                  fill = "Guano_Zone",
                  palette = guano_colours) +
  ggtitle(ytitle5) +
  theme(plot.title = element_text(color="black", size=12, face="bold"),
        axis.ticks.length=unit(0.15, "cm"))

# Residuals vs fitted values plot
modf <- fortify(model_5)
res <- ggplot(modf, aes(x = .fitted, y = .resid)) + geom_point() +
  xlab('Residuals') +
  ylab('Fitted Values Residuals') + 
  ggtitle('Residuals plot') +
  theme(axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm"))

# Violin + box plots for x and y - to add as F below in 3x3 grid --------

# set Guano_sum column to factor rather than numerical plot
df$Guano_Sum <- as.factor(df$Guano_Sum)

v5 <- ggplot(df, aes(x=Guano_Sum, y=Sr_Ti_CS, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v5

# violin plot with dot plot - using jitter to offset the data
# vdot <- v5 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter5 <- v5 + geom_jitter(shape=16, position=position_jitter(0.05))

# Final violin + box plot plot
vbp5 <- vjitter5 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab5_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp5

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter5 <- v5 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp5_nolegend <- vjitter5 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  xlab5_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot")
vbp5_nolegend

# draw violin plot for XRF subsample data
v5.2 <- ggplot(df, aes(x=Guano_Sum, y=Sr_Ti, fill=Guano_Sum)) + 
  geom_violin(trim=FALSE) + 
  scale_fill_manual(values=c("#969696", "#BDBDBD", "#74C476"))
v5.2

# violin plot with dot plot - using jitter to offset the data
# vdot <- v1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
vjitter5.2 <- v5.2 + geom_jitter(shape=16, position=position_jitter(0.05))

# final violin + box plot plot
vbp5.2 <- vjitter5.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab5_y1 +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(1, 0), legend.position = c(1,0),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=12, colour = "black"),
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp5.2

# final violin + box plot plot without legend for 3x3 matrix plot
vjitter5.2 <- v5.2 + geom_jitter(shape=16, position=position_jitter(0.02))
vbp5.2_nolegend <- vjitter5.2 + geom_boxplot(width=0.1, fill = "white", alpha = 0.8) +
  ylab5_y1 +
  theme(legend.position = "none",
        axis.title=element_text(size=12, colour = "black"), 
        axis.text=element_text(size=12, colour = "black"),
        axis.ticks.length=unit(0.15, "cm")) +
  ggtitle("Violin Plot") 
vbp5.2_nolegend

# plot as 3 x3 grid - one example of violin plot 
ggarrange(den.x5, qq.x5, den.y5, qq.y5, res, vbp5.2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
ggsave("Figures/Calibration/Model5_Sr/Fig 1_Regression_summary1.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# plot as 2 x2 grid - qq and violin plots
ggarrange( qq.x5, qq.y5, vbp5_nolegend, vbp5.2_nolegend, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave("Figures/Calibration/Model5_Sr/Fig 2_Regression_summary2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

# Shapiro-Wilk test ----------

# Use this to test normality of one variable (univariate) at a time 
shapiro.test(x5.reg)
shapiro.test(y5.reg)
# If the p-value > 0.05 it implies that the distribution of the data 
# are not significantly different from normal distribution - i.e., can assume normality
# < 0.05 means data are not normally distributed ...
# Note that 'perfect' normality is hard to achieve in most downcore geochemical datasets!

# Save console outputs above to file
sink("Output/Calibration/model5_shapiro_x.txt")
print(shapiro.test(x5.reg))
sink("Output/Calibration/model5_shapiro_y.txt")
print(shapiro.test(y5.reg))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
shapiro.test(x5.reg)
shapiro.test(y5.reg)

# Durbin-Watson Test  ---------------------------------------------------

# Use this to test for correlation among residuals
# H0 (null hypothesis) = no correlation among the residuals
# HA (alternative hypothesis) = residuals are autocorrelated
library(lmtest)
dwtest(model_5)
# Value of 2.0 indicates there is no autocorrelation detected in the sample. 
# 0-2 = positive autocorrelation; 2-4 = negative autocorrelation.
# Save console outputs above to file
sink("Output/Calibration/model5_dw_test.txt")
print(dwtest(model_5))
sink(file = NULL)
# reset back to consoole output
sink(file = NULL)
dwtest(model_5)


# Model 5 ----------------------------------------------------------------
library(ggpubr)
library(dplyr)

theme_set(theme_bw(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold.italic"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7)
))
p5 <- ggplot(data_5,aes(x5.reg, y5.reg)) + 
  geom_point(data = data_5, aes(x5.reg, y5.reg, colour = Guano_Zone, fill = Guano_Zone, shape = Guano_Zone), size = 2, alpha = 1) +
  scale_shape_manual(values = plot_shapes) +
  scale_fill_manual(values = guano_colours1) +
  scale_color_manual(values = guano_colours1) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 7, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "right", #c(1,0 ),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text=element_text(size=7, colour = "black"), axis.title=element_text(size=7, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2, byrow = TRUE)) +
  geom_smooth(method = "lm", se = TRUE, col = "blue", fill = "#9ECAE1") +
  xlab5 +
  ylab5
p5

# Add prediction intervals
p5.1 <- p5 + geom_line(aes(y = lwr), color = "blue", linetype = "dashed") +
  geom_line(aes(y = upr), color = "blue", linetype = "dashed")  +
  ggtitle("Linear regression (solid blue), 95% CI (blue shaded) & prediction (dashed blue)") + 
  theme(axis.ticks.length=unit(0.15, "cm"))
p5.1

# Final regression plot as 1x2 grid
ggarrange(p5.1,  ncol = 1, nrow = 1, labels = c("A"))
ggsave("Figures/Calibration/Model5_Sr/Fig 3_Regression_final.pdf",
       height = c(10), width = c(10), dpi = 600, units = "cm")







# THE END! --------------------------------------------------------------------



