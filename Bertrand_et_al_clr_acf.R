# load packages needed
library(compositions)
library(tidyverse)
library(PeriodicTable)

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

#set working directory (example)
setwd("D:/Dropbox/BAS/Data/R/Papers/Bertrand_2022/Data") # windows
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Bertrand_2022/Data") # mac
setwd("/Users/sjro/Dropbox/BAS/Data/R/Papers/Bertrand_2022/Data") # Macbook M2
getwd()

# assign elements to from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb
rm(periodicTable)

# Import data
Chile_xrf <- read_csv("Bertrand_Chile.csv")
Chile_xrf

# Apply clr 
Chile_xrf_clr <- Chile_xrf %>% 
  select(any_of(elementsList)) %>% # remove non element and scatter columns 
  filter(!if_any(everything(), is.na)) %>% # remove rows with all NAs
  clr()
Chile_xrf_clr <- as_tibble(Chile_xrf_clr) %>% 
  mutate(checksum = rowSums(across(Si:Pb))) %>% 
  mutate_if(is.numeric, round, digits = 3)
Chile_xrf_clr

# Save output
write.csv(Chile_xrf_clr,"Bertrand_Chile_clr.csv", row.names = FALSE)


# Use acf to remove 'noisy' elements from XRF-CS dataset

# Use autocorrelation function (acf) and plots to explore noise in a time-series
library(forecast)
library(ggrepel)
library(directlabels)

# Filter elements based on acf lag thresholds
# set lag threshold to 20 for whole dataset
# 0.2 or 0.1 as  minimum threshold
# 0.5 as maximum threshold - i.e, lag time to half correlation coefficient 

# define filter and lag thresholds
acf_thres_min <- 0.1
acf_thres_max <- 0.5
lag_thres <- 20

# MINIMUM ACF threshold element filtering ACF > 0.1 - ALL SITES DATA ----------------------------------------

Fig1a <- apply(Chile_xrf %>% select(any_of(elementsList)), 
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", linewidth=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig1a)
ggsave("Figures/Fig1a_ACF_all_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# apply acf based filtering to Chile_xrf elmement list - leaving acf elements >0.1 min threshold
apply(Chile_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_min) %>% 
  pull(elements) %>% 
  ordered() -> acfElements_min 
acfElements_min

acfElementsList_min <- select(Chile_xrf, any_of(acfElements_min)) %>% 
  names()
acfElementsList_min

# Replot with min acf filtered elements only
Fig1b <- apply(Chile_xrf %>% select(any_of(acfElements_min)), 
                 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 50) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3.5b)
ggsave("Figures/Fig1b_ACF_elements.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Summary ACF element plot vs Depth_cm
Fig2 <- Chile_xrf %>% 
  #filter(qc == TRUE) %>% 
  select(any_of(acfElements_min), Depth_cm, label) %>%
  pivot_longer(any_of(acfElements_min), names_to = "param", values_to = "element") %>%
  filter(param %in% acfElements_min) %>%
  mutate(param = fct_relevel(param, acfElementsList_min)) %>%
  ggplot(aes(x = element, y = Depth_cm)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 5)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth_cm [cm]") +
  ggtitle("Intensity (cps), ACF filtered elements ACF min >0.1 at lag=20")
print(Fig3.6)
ggsave("Figures/Fig2_All_Sites_ACFmin.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig2, # First row
          ggarrange(Fig1a, Fig1b, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Figures/Fig2_ACF_min.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# MAXIMUM ACF threshold element filtering - all sites ----------------------------------------

# ACF > 0.5 - ALL SITES DATA
# apply acf based filtering to Chile_xrf elmement list - leaving acf elements >0.1 min threshold
apply(Chile_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_max) %>% 
  pull(elements) %>% 
  ordered() -> acfElements_max 
acfElements_max

acfElementsList_max <- select(Chile_xrf, any_of(acfElements_max)) %>% 
  names()
acfElementsList_max

# Replot with max acf filtered elements only
Fig3 <- apply(Chile_xrf %>% select(any_of(acfElements_max)), 
                2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 50) +
  ylim(0.3, 1) +
  theme(legend.position="NULL") +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.5)) # adds element labels to end
print(Fig3)
ggsave("Figures/Fig3_ACF_max_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary acf element plot vs Depth_cm
library(tidypaleo)
theme_set(theme_bw(base_size=8))
Fig4 <- Chile_xrf %>% 
  #filter(qc == TRUE) %>% 
  select(any_of(acfElements_max), Depth_cm, label) %>%
  pivot_longer(any_of(acfElements_max), names_to = "param", values_to = "element") %>%
  filter(param %in% acfElements_max) %>%
  mutate(param = fct_relevel(param, acfElementsList_max)) %>%
  ggplot(aes(x = element, y = Depth_cm)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 5)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth_cm [cm]") +
  ggtitle("Intensity (cps), ACF filtered elements ACF max >0.5 at lag=20")
print(Fig3.9)
ggsave("Figures/Fig3.9_ACFmax_key_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots per row
ggarrange(Fig4, # First row
          ggarrange(Fig1a, Fig3, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Figures/Fig4_ACF_max.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")