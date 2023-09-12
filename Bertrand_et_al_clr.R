# clear previous console
remove (list = ls())
# clear plot window (if needed)
dev.off()

# load packages needed
library(tidyverse) # R working environment - better/easier than base R
library(PeriodicTable) # list of all elements in periodic table
library(compositions) # package for analysing geochemcial compositional data
# for details see: https://cran.r-project.org/web/packages/compositions/index.html 
# Reference: van den Boogaart (2023) version 2.0.6: http://www.stat.boogaart.de/compositions/ 

# Set working directory (example for Steve’s folder structure - replace with your own folder structure)
# setwd("D:/Dropbox/BAS/Data/R/Papers/Bertrand_2022/Data") # windows
setwd("/Users/sjro/Dropbox/BAS/Data/R/Papers/Bertrand_2022/Data") # mac

# Use PeriodicTable package to make 'elementsList' of elements to use from imported dataset (then unload the package) 
data(periodicTable)
elementsList <- periodicTable$symb
rm(periodicTable)

# Import data
Chile_xrf <- read_csv("Bertrand_Chile1.csv")
Chile_xrf

# Apply clr() using Compositions package
Chile_xrf_clr0 <- Chile_xrf %>% 
  select(any_of(elementsList)) %>% # removes any non-element and scatter columns 
  mutate_at(vars(any_of(elementsList)), ## Replace zeros in each data column with half minimum value of that column to allow linear modelling to work
            ~ (. == 0) * min(.[. != 0])/2 + .) %>% # Recommended procedure from Bertrand et al. (submitted) - retains dataframe structure
  #mutate_if(is.numeric, ~na_if(., 0)) %>%  # alternate method - replace any remaining zeroes with NA – comment out mutate code line above if using this
  #filter(!if_any(everything(), is.na)) %>% # alternate method - remove rows with NAs – comment out mutate code line above if using this
  clr()

# Check sum = zero
Chile_xrf_clr0 <- as_tibble(Chile_xrf_clr0) %>% 
  mutate(checksum = rowSums(across(any_of(elementsList)))) %>% 
  mutate_if(is.numeric, round, digits = 3) # round to 3 d.p. 

# Add depth column from original dataframe to clr dataframe
Chile_xrf_depth <- Chile_xrf %>%
  #mutate_if(is.numeric, ~na_if(., 0)) %>%  # alternate method - replace any zeroes with NA to match alternate method dataframe structure
  #filter(!if_any(everything(), is.na)) %>% # alternate method - remove rows with NAs to match alternate method dataframe structure
  select(Depth_cm)
Chile_xrf_clr <-bind_cols(Chile_xrf_depth, Chile_xrf_clr0)
Chile_xrf_clr

# Save output to file in current working directory
write.csv(Chile_xrf_clr,"Bertrand_Chile1_clr.csv", row.names = FALSE)

