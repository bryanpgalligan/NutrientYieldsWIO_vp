## Load packages

packages <- c(
  "dplyr", #for data wrangling
  "corrplot", #visualize correlation matrices
  "forcats", #for reordering categorical variables
  "ggeffects", #generate model predictions
  "ggmap", #mapping tools
  "ggplot2", #plotting
  "ggpubr", #combining multiple ggplots
  "ggsci", #fun palettes for ggplot
  "glmmTMB", #linear mixed models
  "gtable", #extract legend from ggplot and save as separate ggplot object
  "randomForest", #random forest models
  "readr", #loading *.csv files
  "readxl", #loading *.xlsx files
  "rfishbase", #programmatic access to FishBase
  "scales", #format ratios as percent on plot axis
  "sf", #maps in ggplot2
  "stringdist", #fuzzy text search for character strings
  "stringr", #data wrangling for character strings
  "taxize", #taxonomic data from around the web - package is down as of May 30, 2023
  "tidyr", #for tidy data
  "worrms" #World Register of Marine Species
  ) 
lapply(packages, require, character.only = TRUE)

